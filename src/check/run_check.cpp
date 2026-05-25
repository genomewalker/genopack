#include "run_check.hpp"
#include "pass_a.hpp"
#include "pass_b.hpp"
#include "pack_iface.hpp"
#include <genopack/archive.hpp>
#include <genopack/qual.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/toc.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

namespace genopack::check {

namespace {

std::vector<std::string> read_accession_list(const std::filesystem::path& p) {
    std::vector<std::string> result;
    std::ifstream f(p);
    if (!f) throw std::runtime_error("cannot open accession list: " + p.string());
    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty()) result.push_back(std::move(line));
    }
    return result;
}

std::vector<std::filesystem::path> collect_gpk_paths(const std::filesystem::path& pack_path) {
    std::vector<std::filesystem::path> paths;
    if (pack_path.extension() == ".gpk") {
        paths.push_back(pack_path);
        return paths;
    }
    for (const auto& entry : std::filesystem::directory_iterator(pack_path)) {
        if (entry.path().extension() == ".gpk")
            paths.push_back(entry.path());
    }
    std::sort(paths.begin(), paths.end());
    return paths;
}

void write_qual_to_archive(const std::filesystem::path& gpk_path,
                           const std::unordered_map<std::string, GenomeQuality>& quality,
                           const std::unordered_map<std::string, GenomeId>& acc_to_id)
{
    int lock_fd = ::open(gpk_path.c_str(), O_RDWR);
    if (lock_fd < 0)
        throw std::runtime_error("check: cannot open for locking: " + gpk_path.string());
    struct LockGuard {
        int fd;
        explicit LockGuard(int fd) : fd(fd) { ::flock(fd, LOCK_EX); }
        ~LockGuard() { ::flock(fd, LOCK_UN); ::close(fd); }
    } lock(lock_fd);

    MmapFileReader mmap;
    mmap.open(gpk_path);
    auto toc = TocReader::read(mmap);
    const auto* tail = mmap.ptr_at<TailLocator>(mmap.size() - sizeof(TailLocator));
    uint64_t prev_toc_offset = tail->toc_offset;
    uint64_t generation      = toc.header.generation + 1;
    mmap.close();

    QualWriter qw;
    for (const auto& [acc, q] : quality) {
        auto it = acc_to_id.find(acc);
        if (it == acc_to_id.end()) continue;

        QualRecord r{};
        r.genome_id                     = it->second;
        r.completeness_cluster_relative = q.completeness_cluster_relative;
        r.completeness_fragmentation    = q.completeness_fragmentation;
        r.completeness_post_decontam    = q.completeness_post_decontam;
        r.contamination_leakage         = q.contamination_leakage;
        r.contamination_tnf_excess      = q.contamination_tnf_excess;
        r.chromosome_skew_closure       = q.chromosome_skew_closure;
        r.leakage_residual              = q.leakage_residual;
        r.support_tier                  = static_cast<uint8_t>(q.support_tier);
        r.interval_width                = q.interval_width;
        r.self_coherence                = q.self_coherence;
        r.contamination_mixture         = q.contamination_mixture;
        r.mixture_sources               = static_cast<int16_t>(q.mixture_sources);
        r.n_mix_windows                 = q.n_mix_windows;
        r.fiedler_u16                   = static_cast<uint16_t>(
            std::min(1.0f, std::isnan(q.fiedler_value) ? 0.0f : q.fiedler_value) * 65535.0f);
        r.contig_outlier_u8             = static_cast<uint8_t>(
            std::min(1.0f, std::isnan(q.contamination_contig_outlier) ? 0.0f : q.contamination_contig_outlier) * 255.0f);
        r.spe_outlier_u8                = static_cast<uint8_t>(
            std::min(1.0f, std::isnan(q.contamination_spe) ? 0.0f : q.contamination_spe) * 255.0f);

        r.qual_flags                    = q.qual_flags;
        r.sketch_breadth                = NAN; // not computable at check time (no consensus sketch)
        // Set policy flags from continuous values (thresholds tunable without rebuild).
        if (!std::isnan(q.completeness_fragmentation) && q.completeness_fragmentation < 0.65f)
            r.qual_flags |= QualRecord::QUAL_FLAG_FRAG_GATED;
        if (!std::isnan(q.self_coherence) && q.self_coherence < 0.80f)
            r.qual_flags |= QualRecord::QUAL_FLAG_COHERENCE_GATED;
        qw.add(r);
    }

    AppendWriter writer;
    writer.open_append(gpk_path);

    uint64_t section_id = toc.next_section_id();
    SectionDesc qual_sd = qw.finalize(writer, section_id);

    TocWriter new_toc;
    for (const auto& sd : toc.sections) {
        if (sd.type != SEC_QUAL) new_toc.add_section(sd);
    }
    new_toc.add_section(qual_sd);

    new_toc.finalize(writer,
                     generation,
                     toc.header.live_genome_count,
                     toc.header.total_genome_count,
                     prev_toc_offset,
                     toc.header.catalog_root_section_id,
                     toc.header.accession_root_section_id,
                     toc.header.tombstone_root_section_id);

    spdlog::info("check: wrote {} QUAL records to {}", qw.size(), gpk_path.string());
}

} // namespace

int cmd_check(const std::filesystem::path& pack_path,
              const std::filesystem::path& genomes_file,
              int threads,
              int min_genus_size,
              float leakage_threshold,
              const std::filesystem::path& output,
              bool recompute)
{
    auto gpk_paths = collect_gpk_paths(pack_path);
    if (gpk_paths.empty())
        throw std::runtime_error("check: no .gpk files found at " + pack_path.string());

    std::unordered_map<std::string, GenomeQuality> all_quality;

    for (const auto& gp : gpk_paths) {
        spdlog::info("check: processing {}", gp.string());

        ArchiveReader ar;
        ar.open(gp);
        SingleArchiveCheckReader pack(ar);

        // Collect accessions for this archive
        std::vector<std::string> accessions;
        if (!genomes_file.empty()) {
            auto all_acc = read_accession_list(genomes_file);
            std::unordered_set<std::string> filter(all_acc.begin(), all_acc.end());
            ar.scan_genome_accessions([&](std::string_view acc, GenomeId) {
                if (filter.count(std::string(acc)))
                    accessions.emplace_back(acc);
            });
        } else {
            ar.scan_genome_accessions([&](std::string_view acc, GenomeId) {
                accessions.emplace_back(acc);
            });
        }

        if (accessions.empty()) {
            spdlog::debug("check: no accessions in {}, skipping", gp.string());
            continue;
        }
        spdlog::info("check: {} accessions in {}", accessions.size(), gp.filename().string());

        // Load existing QUAL section as a cache to skip the SKCH scan.
        std::unordered_map<uint64_t, QualRecord> qual_cache;
        const std::unordered_map<uint64_t, QualRecord>* qual_cache_ptr = nullptr;
        if (!recompute && ar.has_qual()) {
            ar.scan_qual([&](const QualRecord& r) {
                qual_cache.emplace(r.genome_id, r);
            });
            spdlog::info("check: loaded {} cached QUAL records (use --recompute to force rescan)",
                         qual_cache.size());
            qual_cache_ptr = &qual_cache;
        }

        auto pass_a = run_pass_a(pack, accessions, threads, min_genus_size,
                                 leakage_threshold, qual_cache_ptr);
        auto quality = pass_a.quality;

        PassBConfig pb_cfg;
        run_pass_b(pack, pass_a, quality, pb_cfg, threads, qual_cache_ptr);

        // Build acc → GenomeId for QUAL write
        std::unordered_map<std::string, GenomeId> acc_to_id;
        acc_to_id.reserve(accessions.size());
        ar.scan_genome_accessions([&](std::string_view acc, GenomeId gid) {
            acc_to_id.emplace(std::string(acc), gid);
        });

        ar.close();

        write_qual_to_archive(gp, quality, acc_to_id);

        for (auto& [acc, q] : quality)
            all_quality.emplace(acc, std::move(q));
    }

    if (output.empty() || all_quality.empty()) {
        spdlog::info("check: done");
        return 0;
    }

    // TSV output (aggregated across all parts)
    std::ofstream tsv(output);
    if (!tsv) throw std::runtime_error("check: cannot open output: " + output.string());
    tsv << "accession\tcompleteness_cluster_relative\tcompleteness_fragmentation"
           "\tcompleteness_post_decontam\tcontamination_leakage\tcontamination_tnf_excess"
           "\tchromosome_skew_closure\tleakage_residual\tself_coherence"
           "\tchargaff_parity\tspectral_gap\tscale_kink\tcontamination_mixture\tmixture_sources"
           "\tn_mix_windows\tfiedler_value\tcontamination_contig_outlier\tcontamination_spe\tqual_flags\tsupport_tier\tinterval_width\n";

    auto tier_str = [](SupportTier t) -> const char* {
        switch (t) {
        case SupportTier::GenusSaturated: return "GenusSaturated";
        case SupportTier::GenusSparse:    return "GenusSparse";
        case SupportTier::Singleton:      return "Singleton";
        }
        return "Unknown";
    };
    auto fmt = [](float v) -> std::string {
        return std::isnan(v) ? "NA" : std::to_string(v);
    };

    for (const auto& [acc, q] : all_quality) {
        tsv << acc << '\t'
            << fmt(q.completeness_cluster_relative) << '\t'
            << fmt(q.completeness_fragmentation) << '\t'
            << fmt(q.completeness_post_decontam) << '\t'
            << fmt(q.contamination_leakage) << '\t'
            << fmt(q.contamination_tnf_excess) << '\t'
            << fmt(q.chromosome_skew_closure) << '\t'
            << fmt(q.leakage_residual) << '\t'
            << fmt(q.self_coherence) << '\t'
            << fmt(q.chargaff_parity) << '\t'
            << fmt(q.spectral_gap) << '\t'
            << fmt(q.scale_kink) << '\t'
            << fmt(q.contamination_mixture) << '\t'
            << q.mixture_sources << '\t'
            << q.n_mix_windows << '\t'
            << fmt(q.fiedler_value) << '\t'
            << fmt(q.contamination_contig_outlier) << '\t'
            << fmt(q.contamination_spe) << '\t'
            << static_cast<int>(q.qual_flags) << '\t'
            << tier_str(q.support_tier) << '\t'
            << q.interval_width << '\n';
    }
    spdlog::info("check: TSV written to {}", output.string());
    return 0;
}

} // namespace genopack::check
