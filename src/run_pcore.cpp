#include "run_pcore.hpp"
#include <genopack/archive.hpp>
#include <genopack/aamer.hpp>
#include <genopack/pcore.hpp>
#include <genopack/gcov.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/section_checksum.hpp>
#include <genopack/toc.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <fstream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

namespace genopack {

namespace {

std::vector<std::filesystem::path> collect_gpk_paths(const std::filesystem::path& p) {
    std::vector<std::filesystem::path> v;
    if (p.extension() == ".gpk") { v.push_back(p); return v; }
    for (const auto& e : std::filesystem::directory_iterator(p))
        if (e.path().extension() == ".gpk") v.push_back(e.path());
    std::sort(v.begin(), v.end());
    return v;
}

std::string genus_of(std::string_view tax) {
    constexpr std::string_view needle = ";g__";
    auto pos = tax.find(needle);
    std::string_view g;
    if (pos != std::string_view::npos) {
        auto start = pos + 1;
        auto end   = tax.find(';', start);
        g = tax.substr(start, end == std::string_view::npos ? end : end - start);
    } else if (tax.starts_with("g__")) {
        auto end = tax.find(';', 3);
        g = tax.substr(0, end);
    }
    if (g.empty() || g == "g__") return {};
    return std::string(g);
}

std::string family_of(std::string_view tax) {
    constexpr std::string_view needle = ";f__";
    auto pos = tax.find(needle);
    std::string_view f;
    if (pos != std::string_view::npos) {
        auto start = pos + 1;
        auto end   = tax.find(';', start);
        f = tax.substr(start, end == std::string_view::npos ? end : end - start);
    } else if (tax.starts_with("f__")) {
        auto end = tax.find(';', 3);
        f = tax.substr(0, end);
    }
    if (f.empty() || f == "f__") return {};
    return std::string(f);
}

void genome_aamers(const std::string& fasta, int k, int min_seg_aa, uint64_t frac_max,
                   std::vector<uint64_t>& out) {
    out.clear();
    const char* p = fasta.data();
    const char* end = p + fasta.size();
    std::string seq;
    while (p < end) {
        while (p < end && *p != '>') ++p;
        if (p >= end) break;
        while (p < end && *p != '\n') ++p;
        if (p < end) ++p;
        seq.clear();
        while (p < end && *p != '>') {
            if (*p != '\n' && *p != '\r') seq.push_back(*p);
            ++p;
        }
        if (!seq.empty()) extract_aamers_dna_into(seq, k, min_seg_aa, frac_max, out);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
}

void append_pcore(const std::filesystem::path& gpk, PcoreWriter& pw, uint64_t frac_max) {
    int fd = ::open(gpk.c_str(), O_RDWR);
    if (fd < 0) throw std::runtime_error("pcore: cannot open for locking: " + gpk.string());
    struct LG { int fd; explicit LG(int f):fd(f){::flock(fd,LOCK_EX);} ~LG(){::flock(fd,LOCK_UN);::close(fd);} } lg(fd);

    MmapFileReader mmap;
    mmap.open(gpk);
    auto toc = TocReader::read(mmap);
    const auto* tail = mmap.ptr_at<TailLocator>(mmap.size() - sizeof(TailLocator));
    uint64_t prev_toc_offset = tail->toc_offset;
    uint64_t generation      = toc.header.generation + 1;
    mmap.close();

    AppendWriter writer;
    writer.open_append(gpk);
    SectionDesc sd = pw.finalize(writer, toc.next_section_id(), frac_max);

    TocWriter nt;
    for (const auto& s : toc.sections) if (s.type != SEC_PCORE) nt.add_section(s);
    nt.add_section(sd);

    writer.flush();
    stamp_section_checksums(gpk.c_str(), nt.sections(), /*only_if_zero=*/true);
    nt.finalize(writer, generation, toc.header.live_genome_count, toc.header.total_genome_count,
                prev_toc_offset, toc.header.catalog_root_section_id,
                toc.header.accession_root_section_id, toc.header.tombstone_root_section_id);
}

} // namespace

int cmd_pcore(const std::filesystem::path& pack_path, int threads,
              const std::filesystem::path& members_list,
              float bootstrap_theta, uint32_t bootstrap_min_seg_aa,
              uint64_t bootstrap_frac_max_hash) {
    auto gpks = collect_gpk_paths(pack_path);
    if (gpks.empty()) throw std::runtime_error("pcore: no .gpk found at " + pack_path.string());
    if (threads < 1) threads = 1;

    std::unordered_set<std::string> ref_members;
    if (!members_list.empty()) {
        std::ifstream f(members_list);
        if (!f) throw std::runtime_error("pcore: cannot open members list " + members_list.string());
        std::string line;
        while (std::getline(f, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (!line.empty()) ref_members.insert(line);
        }
        spdlog::info("pcore: restricting to {} reference members", ref_members.size());
    }

    for (const auto& gp : gpks) {
        ArchiveReader ar;
        ar.open(gp);
        // Build config comes from an existing SEC_CORE header when present. SEC_CORE has
        // no in-repo builder (see docs/aamer.md), so when absent we bootstrap from the
        // CLI-supplied canonical params instead of skipping the archive.
        int k; int min_seg; uint64_t frac_max; float theta;
        if (ar.has_core()) {
            const CoreReader* cr = ar.core_reader();
            k = cr->k(); min_seg = cr->min_seg_aa();
            frac_max = cr->frac_max_hash();
            theta = cr->theta();               // capture before ar.close() (cr dangles after)
        } else {
            k = AAMER_K; min_seg = static_cast<int>(bootstrap_min_seg_aa);
            frac_max = bootstrap_frac_max_hash; theta = bootstrap_theta;
            spdlog::info("pcore: {} has no SEC_CORE; bootstrapping build config "
                         "(k={} min_seg_aa={} theta={:.2f} frac_max_hash={:#018x})",
                         gp.filename().string(), k, min_seg, theta, frac_max);
        }

        std::unordered_map<std::string, std::string> acc_genus, acc_family;
        ar.scan_taxonomy([&](std::string_view acc, std::string_view tax) {
            std::string g = genus_of(tax);
            if (!g.empty()) acc_genus.emplace(std::string(acc), std::move(g));
            std::string f = family_of(tax);
            if (!f.empty()) acc_family.emplace(std::string(acc), std::move(f));
        });

        std::vector<std::string> accessions, idx_genus, idx_family;
        ar.scan_genome_accessions([&](std::string_view acc, GenomeId gid) {
            if (ar.is_deleted(gid)) return;
            if (!ref_members.empty() && !ref_members.count(std::string(acc))) return;
            auto it = acc_genus.find(std::string(acc));
            if (it == acc_genus.end()) return;
            accessions.emplace_back(acc);
            idx_genus.push_back(it->second);
            auto fit = acc_family.find(std::string(acc));
            idx_family.push_back(fit != acc_family.end() ? fit->second : std::string());
        });
        if (accessions.empty()) { spdlog::warn("pcore: {} no genus-labeled members; skipping", gp.filename().string()); ar.close(); continue; }

        // Accumulate per-genus AND per-family member aamer sets (both tiers, one section).
        std::unordered_map<std::string, std::vector<std::vector<uint64_t>>> genus_members, family_members;
        std::mutex mtx;
        ar.visit_shard_batches_parallel(accessions, threads,
            [&](ArchiveReader::ShardBatch& batch) {
                thread_local std::vector<uint64_t> qm;
                for (auto& [idx, g] : batch) {
                    genome_aamers(g.fasta, k, min_seg, frac_max, qm);
                    if (qm.empty()) continue;
                    const std::string& genus  = idx_genus[idx];
                    const std::string& family = idx_family[idx];
                    std::lock_guard<std::mutex> lk(mtx);
                    if (!family.empty()) family_members[family].push_back(qm);  // copy
                    genus_members[genus].push_back(std::move(qm));               // move
                }
            });
        ar.close();

        PcoreWriter pw(static_cast<uint32_t>(k), static_cast<uint32_t>(min_seg), theta);
        size_t n_genus = 0, n_family = 0;
        for (auto& [genus, members] : genus_members)  { pw.add_from_members(GcovWriter::hash_genus(genus), members);  ++n_genus; }
        for (auto& [family, members] : family_members){ pw.add_from_members(GcovWriter::hash_genus(family), members); ++n_family; }
        if (pw.n_entries() == 0) { spdlog::warn("pcore: {} produced no entries", gp.filename().string()); continue; }
        append_pcore(gp, pw, frac_max);
        spdlog::info("pcore: {} — {} genus + {} family references (dense, count-annotated, theta={:.2f}), model_hash={:#018x}",
                     gp.filename().string(), n_genus, n_family, theta, pw.model_hash());
    }
    return 0;
}

} // namespace genopack
