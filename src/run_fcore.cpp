#include "run_fcore.hpp"
#include <genopack/archive.hpp>
#include <genopack/aamer.hpp>
#include <genopack/core_section.hpp>
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

std::vector<std::filesystem::path> collect_gpk_paths(const std::filesystem::path& pack_path) {
    std::vector<std::filesystem::path> paths;
    if (pack_path.extension() == ".gpk") { paths.push_back(pack_path); return paths; }
    for (const auto& e : std::filesystem::directory_iterator(pack_path))
        if (e.path().extension() == ".gpk") paths.push_back(e.path());
    std::sort(paths.begin(), paths.end());
    return paths;
}

// Parse the f__ family token from a GTDB-style taxonomy string (matches pass_a).
std::string family_of(std::string_view tax) {
    constexpr std::string_view fneedle = ";f__";
    auto fpos = tax.find(fneedle);
    std::string_view fam;
    if (fpos != std::string_view::npos) {
        auto fstart = fpos + 1;
        auto fend   = tax.find(';', fstart);
        fam = tax.substr(fstart, fend == std::string_view::npos ? fend : fend - fstart);
    } else if (tax.starts_with("f__")) {
        auto fend = tax.find(';', 3);
        fam = tax.substr(0, fend);
    }
    if (fam.empty() || fam == "f__") return {};
    return std::string(fam);
}

// Per-genome sorted-unique aamer set, extracted per-contig (no cross-contig frames).
void genome_aamers(const std::string& fasta, int k, int min_seg_aa, uint64_t frac_max,
                   std::vector<uint64_t>& out) {
    out.clear();
    const char* p = fasta.data();
    const char* end = p + fasta.size();
    std::string seq;
    while (p < end) {
        while (p < end && *p != '>') ++p;             // to header
        if (p >= end) break;
        while (p < end && *p != '\n') ++p;            // skip header line
        if (p < end) ++p;
        seq.clear();
        while (p < end && *p != '>') {                // collect contig sequence
            if (*p != '\n' && *p != '\r') seq.push_back(*p);
            ++p;
        }
        if (!seq.empty())
            extract_aamers_dna_into(seq, k, min_seg_aa, frac_max, out);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
}

// Append a fresh SEC_FCORE (dropping any prior) to a single .gpk.
void append_fcore(const std::filesystem::path& gpk, CoreWriter& fcw, uint64_t frac_max) {
    int lock_fd = ::open(gpk.c_str(), O_RDWR);
    if (lock_fd < 0) throw std::runtime_error("fcore: cannot open for locking: " + gpk.string());
    struct LockGuard { int fd; explicit LockGuard(int f):fd(f){::flock(fd,LOCK_EX);}
                       ~LockGuard(){::flock(fd,LOCK_UN);::close(fd);} } lock(lock_fd);

    MmapFileReader mmap;
    mmap.open(gpk);
    auto toc = TocReader::read(mmap);
    const auto* tail = mmap.ptr_at<TailLocator>(mmap.size() - sizeof(TailLocator));
    uint64_t prev_toc_offset = tail->toc_offset;
    uint64_t generation      = toc.header.generation + 1;
    mmap.close();

    AppendWriter writer;
    writer.open_append(gpk);
    uint64_t section_id = toc.next_section_id();
    SectionDesc fcore_sd = fcw.finalize(writer, section_id, frac_max);

    TocWriter new_toc;
    for (const auto& sd : toc.sections)
        if (sd.type != SEC_FCORE) new_toc.add_section(sd);
    new_toc.add_section(fcore_sd);

    writer.flush();
    stamp_section_checksums(gpk.c_str(), new_toc.sections(), /*only_if_zero=*/true);

    new_toc.finalize(writer, generation,
                     toc.header.live_genome_count, toc.header.total_genome_count,
                     prev_toc_offset, toc.header.catalog_root_section_id,
                     toc.header.accession_root_section_id, toc.header.tombstone_root_section_id);
}

} // namespace

int cmd_fcore(const std::filesystem::path& pack_path, int threads,
              float theta_override, unsigned min_members,
              const std::filesystem::path& members_list) {
    auto gpks = collect_gpk_paths(pack_path);
    if (gpks.empty()) throw std::runtime_error("fcore: no .gpk found at " + pack_path.string());
    if (threads < 1) threads = 1;

    // Optional reference-member allowlist (one accession per line).
    std::unordered_set<std::string> ref_members;
    if (!members_list.empty()) {
        std::ifstream f(members_list);
        if (!f) throw std::runtime_error("fcore: cannot open members list " + members_list.string());
        std::string line;
        while (std::getline(f, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (!line.empty()) ref_members.insert(line);
        }
        spdlog::info("fcore: restricting cores to {} reference members from {}",
                     ref_members.size(), members_list.filename().string());
    }

    for (const auto& gp : gpks) {
        ArchiveReader ar;
        ar.open(gp);
        if (!ar.has_core()) {
            spdlog::warn("fcore: {} has no CORE section; skipping (build genus cores first)",
                         gp.filename().string());
            ar.close();
            continue;
        }
        const CoreReader* cr = ar.core_reader();
        const int      k        = cr->k();
        const int      min_seg  = cr->min_seg_aa();
        const uint64_t frac_max = cr->frac_max_hash();
        const float    theta    = theta_override > 0.0f ? theta_override : cr->theta();

        // acc -> family (only genomes with a resolvable family).
        std::unordered_map<std::string, std::string> acc_family;
        ar.scan_taxonomy([&](std::string_view acc, std::string_view tax) {
            std::string fam = family_of(tax);
            if (!fam.empty()) acc_family.emplace(std::string(acc), std::move(fam));
        });

        // Build the live accession list + parallel idx -> (family, gid) maps.
        std::vector<std::string> accessions;
        std::vector<std::string> idx_family;
        std::vector<uint64_t>    idx_gid;
        ar.scan_genome_accessions([&](std::string_view acc, GenomeId gid) {
            if (ar.is_deleted(gid)) return;
            if (!ref_members.empty() && !ref_members.count(std::string(acc))) return;
            auto it = acc_family.find(std::string(acc));
            if (it == acc_family.end()) return;
            accessions.emplace_back(acc);
            idx_family.push_back(it->second);
            idx_gid.push_back(static_cast<uint64_t>(gid));
        });
        if (accessions.empty()) {
            spdlog::warn("fcore: {} has no family-labeled live genomes; skipping",
                         gp.filename().string());
            ar.close();
            continue;
        }

        // Stream shards (memory-bounded per shard), extract per-genome aamers, group by family.
        std::unordered_map<std::string, std::vector<std::vector<uint64_t>>> fam_members;
        std::unordered_map<std::string, std::vector<uint64_t>>              fam_ids;
        std::mutex mtx;
        ar.visit_shard_batches_parallel(accessions, threads,
            [&](ArchiveReader::ShardBatch& batch) {
                thread_local std::vector<uint64_t> qm;
                for (auto& [idx, g] : batch) {
                    genome_aamers(g.fasta, k, min_seg, frac_max, qm);
                    if (qm.empty()) continue;
                    const std::string& fam = idx_family[idx];
                    const uint64_t     gid = idx_gid[idx];
                    std::lock_guard<std::mutex> lk(mtx);
                    fam_members[fam].push_back(std::move(qm));   // qm refilled next genome
                    fam_ids[fam].push_back(gid);
                }
            });

        CoreWriter fcw(static_cast<uint32_t>(k), static_cast<uint32_t>(min_seg), theta, SEC_FCORE);
        size_t n_fam = 0, n_skipped = 0;
        for (auto& [fam, members] : fam_members) {
            if (members.size() < min_members) { ++n_skipped; continue; }
            fcw.add_from_members(GcovWriter::hash_genus(fam), members, fam_ids[fam]);
            ++n_fam;
        }
        ar.close();

        if (fcw.n_genera() == 0) {
            spdlog::warn("fcore: {} produced no family cores ({} families < min-members={}); not writing",
                         gp.filename().string(), n_skipped, min_members);
            continue;
        }
        append_fcore(gp, fcw, frac_max);
        spdlog::info("fcore: {} — {} family cores written ({} families skipped < {} members), "
                     "theta={:.2f}, model_hash={:#018x}",
                     gp.filename().string(), n_fam, n_skipped, min_members, theta, fcw.core_model_hash());
    }
    return 0;
}

} // namespace genopack
