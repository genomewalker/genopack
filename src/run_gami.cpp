#include "run_gami.hpp"
#include <genopack/archive.hpp>
#include <genopack/gami.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/section_checksum.hpp>
#include <genopack/toc.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <filesystem>
#include <stdexcept>
#include <vector>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

namespace genopack {

namespace {

std::vector<std::filesystem::path> collect_gpks(const std::filesystem::path& p) {
    std::vector<std::filesystem::path> v;
    if (p.extension() == ".gpk") { v.push_back(p); return v; }
    // Caller may have passed the output stem (e.g. "pack" → actual file is "pack.gpk")
    auto with_ext = std::filesystem::path(p.string() + ".gpk");
    if (std::filesystem::is_regular_file(with_ext)) { v.push_back(with_ext); return v; }
    for (const auto& e : std::filesystem::directory_iterator(p))
        if (e.path().extension() == ".gpk") v.push_back(e.path());
    std::sort(v.begin(), v.end());
    return v;
}

void append_gami(const std::filesystem::path& gpk, const GlobalMultiplicityIndex& gmi,
                 uint64_t frac_max) {
    int fd = ::open(gpk.c_str(), O_RDWR);
    if (fd < 0) throw std::runtime_error("gami: cannot open for locking: " + gpk.string());
    struct LG {
        int fd; explicit LG(int f) : fd(f) { ::flock(fd, LOCK_EX); }
        ~LG() { ::flock(fd, LOCK_UN); ::close(fd); }
    } lg(fd);

    MmapFileReader mmap;
    mmap.open(gpk);
    auto toc = TocReader::read(mmap);
    const auto* tail = mmap.ptr_at<TailLocator>(mmap.size() - sizeof(TailLocator));
    uint64_t prev_toc = tail->toc_offset;
    uint64_t gen      = toc.header.generation + 1;
    mmap.close();

    AppendWriter writer;
    writer.open_append(gpk);

    GamiV2Writer gw;
    SectionDesc sd = gw.finalize(writer, toc.next_section_id(), gmi, frac_max);

    TocWriter nt;
    for (const auto& s : toc.sections) if (s.type != SEC_GAMI) nt.add_section(s);
    nt.add_section(sd);

    writer.flush();
    stamp_section_checksums(gpk.c_str(), nt.sections(), /*only_if_zero=*/true);
    nt.finalize(writer, gen,
                toc.header.live_genome_count, toc.header.total_genome_count,
                prev_toc, toc.header.catalog_root_section_id,
                toc.header.accession_root_section_id, toc.header.tombstone_root_section_id);
}

} // namespace

int cmd_gami_build(const std::filesystem::path& pack_path, int threads) {
    auto gpks = collect_gpks(pack_path);
    if (gpks.empty())
        throw std::runtime_error("gami: no .gpk found at " + pack_path.string());
    if (threads < 1) threads = 1;

    for (const auto& gp : gpks) {
        ArchiveReader ar;
        ar.open(gp);
        if (!ar.has_pcore() && !ar.has_core()) {
            spdlog::warn("gami: {} has no PCORE/CORE; skipping", gp.filename().string());
            ar.close(); continue;
        }

        const bool use_pcore = ar.has_pcore();
        const PcoreReader* pr = use_pcore ? ar.pcore_reader() : nullptr;
        const CoreReader*  cr = (!use_pcore && ar.has_core()) ? ar.core_reader() : nullptr;
        const uint64_t frac_max = use_pcore ? pr->frac_max_hash() : cr->frac_max_hash();
        const uint32_t n_entries = use_pcore ? pr->n_entries() : cr->n_genera();

        spdlog::info("gami: {} — building GMI from {} entries ({}) …",
                     gp.filename().string(), n_entries, use_pcore ? "PCORE" : "CORE");

        GlobalMultiplicityIndex gmi;
        // Pre-sum total aamers across all entries so reserve() starts with the right capacity.
        size_t total_aamers = 0;
        if (use_pcore) {
            for (uint32_t i = 0; i < n_entries; ++i) {
                PcoreView v = pr->lookup(pr->key_hash_at(i));
                if (v.valid()) total_aamers += v.total();
            }
        } else {
            for (uint32_t i = 0; i < n_entries; ++i) {
                CoreView v = cr->lookup(cr->genus_hash_at(i));
                if (v.valid()) total_aamers += v.n_aamers;
            }
        }
        gmi.reserve(total_aamers);
        const int par = std::min(threads, 8);
        std::vector<std::vector<uint64_t>> buf(par);
        std::vector<std::vector<float>>    prevf_buf(par);

        for (uint32_t bs = 0; bs < n_entries; bs += static_cast<uint32_t>(par)) {
            const uint32_t be  = std::min(bs + static_cast<uint32_t>(par), n_entries);
            const int      bsz = static_cast<int>(be - bs);
            #pragma omp parallel for schedule(static) num_threads(bsz)
            for (int bi = 0; bi < bsz; ++bi) {
                buf[bi].clear(); prevf_buf[bi].clear();
                if (use_pcore) {
                    const uint64_t gh = pr->key_hash_at(bs + static_cast<uint32_t>(bi));
                    PcoreView v = pr->lookup(gh);
                    if (v.valid()) v.materialize(buf[bi], prevf_buf[bi]);
                } else {
                    const uint64_t gh = cr->genus_hash_at(bs + static_cast<uint32_t>(bi));
                    CoreView v = cr->lookup(gh);
                    if (v.valid()) buf[bi].assign(v.aamers, v.aamers + v.n_aamers);
                }
            }
            for (int bi = 0; bi < bsz; ++bi) {
                for (uint64_t h : buf[bi]) gmi.increment(h);
                buf[bi].clear();
                gmi.maybe_resize();
            }
        }
        ar.close();

        spdlog::info("gami: {} — GMI: {} M unique aamers, {:.1f} GB; writing SEC_GAMI v2 …",
                     gp.filename().string(), gmi.count() / 1'000'000, gmi.bytes() / 1e9);
        append_gami(gp, gmi, frac_max);
        spdlog::info("gami: {} — done", gp.filename().string());
    }
    return 0;
}

} // namespace genopack
