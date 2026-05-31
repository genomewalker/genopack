#include "run_check.hpp"
#include "pass_a.hpp"
#include "pass_b.hpp"
#include "pass_marker.hpp"
#include "pack_iface.hpp"
#include <genopack/archive.hpp>
#include <genopack/qual.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/toc.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cmath>
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
        r.sibling_outlier_u8            = static_cast<uint8_t>(
            std::min(1.0f, std::isnan(q.contamination_sibling_outlier) ? 0.0f : q.contamination_sibling_outlier) * 255.0f);
        r.rho_outlier_u8                = static_cast<uint8_t>(
            std::min(1.0f, std::isnan(q.contamination_rho_outlier) ? 0.0f : q.contamination_rho_outlier) * 255.0f);
        r.cross_genus_u8                = static_cast<uint8_t>(
            std::min(1.0f, std::isnan(q.contamination_cross_genus) ? 0.0f : q.contamination_cross_genus) * 255.0f);
        r.fmh_minority_u8               = static_cast<uint8_t>(
            std::min(1.0f, std::isnan(q.fmh_minority_fraction) ? 0.0f : q.fmh_minority_fraction) * 255.0f);
        // marker_completeness: 0=not_scored, 1-255=(value-1)/254 → range [0,1] with sentinel 0
        if (!std::isnan(q.marker_completeness))
            r.marker_completeness_u8 = static_cast<uint8_t>(
                std::clamp(q.marker_completeness, 0.0f, 1.0f) * 254.0f + 1.0f);
        else
            r.marker_completeness_u8 = 0;

        r.chargaff_parity               = q.chargaff_parity;
        r.spectral_gap                  = q.spectral_gap;
        r.scale_kink                    = q.scale_kink;

        r.qual_flags                    = q.qual_flags;
        // Mark GCOV scoring as complete. Set whenever pass-B ran for this genome
        // (contamination signals may be NAN when no GCOV model exists — still correct).
        if (!std::isnan(q.self_coherence))
            r.qual_flags |= QualRecord::QUAL_FLAG_GCOV_SCORED;
        if (!std::isnan(q.marker_redundancy)) {
            const float clamped = std::min(0.09999f, q.marker_redundancy);
            r.marker_redundancy_u16 = static_cast<uint16_t>(clamped / 0.1f * 65534.0f);
            r.qual_flags |= QualRecord::QUAL_FLAG_MARKER_SCORED;
        } else {
            r.marker_redundancy_u16 = 0xFFFFu;
        }
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
              bool recompute,
              const std::filesystem::path& markers_path)
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
        if (!markers_path.empty()) pb_cfg.markers_path = markers_path.string();
        run_pass_b(pack, pass_a, quality, pb_cfg, threads, qual_cache_ptr);

        // Unconditional marker pass: score every genome with a valid genus assignment,
        // not just Pass-B-flagged ones. Fills marker_completeness for clean/complete
        // genomes that pass_b skips due to low TNF excess. Makes the SCG completeness
        // metric available for all genomes (comparable to CheckM2's unconditional scoring).
        if (!markers_path.empty()) {
            PassMarkerConfig pm_cfg;
            pm_cfg.markers_path   = markers_path.string();
            pm_cfg.marker_min_hits = pb_cfg.marker_min_hits;
            run_pass_marker(pack, pass_a, quality, pm_cfg, threads);
        }

        // Re-score marker_redundancy_z using in-distribution per-genus calibration.
        // MSA-based calibration (in .mrk file) underestimates by ~16× because gap-spanning
        // k-mers from MSA columns differ from 6-frame translation k-mers used at check time.
        // Genera with ≥10 scored observations override the MSA z-score with a robust estimate.
        for (const auto& [genus, accs] : pass_a.genus_members) {
            std::vector<float> vals;
            for (const auto& acc : accs) {
                auto it = quality.find(acc);
                if (it != quality.end() && !std::isnan(it->second.marker_redundancy))
                    vals.push_back(it->second.marker_redundancy);
            }
            if (vals.size() < 10) continue;

            const size_t n = vals.size();
            std::nth_element(vals.begin(), vals.begin() + n / 2, vals.end());
            const float med = vals[n / 2];

            std::vector<float> devs(n);
            for (size_t i = 0; i < n; ++i) devs[i] = std::abs(vals[i] - med);
            std::nth_element(devs.begin(), devs.begin() + n / 2, devs.end());
            const float mad_eff = std::max({devs[n / 2], 0.01f * med, 1e-4f});

            for (const auto& acc : accs) {
                auto it = quality.find(acc);
                if (it != quality.end() && !std::isnan(it->second.marker_redundancy))
                    it->second.marker_redundancy_z =
                        (it->second.marker_redundancy - med) / (1.4826f * mad_eff);
            }
        }

        // Per-genus p90 calibration of marker_completeness.
        // Fixed denominator (120) underestimates completeness for lineages that
        // naturally lack some universal markers. Use p90 of observed n_present
        // within the genus (≥10 scored members) as the denominator; fall back to
        // the archive-wide p90 for small genera and singletons — this avoids the
        // self-calibration trap (n_present / n_present = 1.0) for singletons while
        // still correcting lineage bias better than a fixed 120.
        {
            std::vector<float> all_np;
            all_np.reserve(quality.size());
            for (const auto& [acc, q] : quality)
                if (q.marker_n_expected > 0 && q.marker_n_present > 0)
                    all_np.push_back(static_cast<float>(q.marker_n_present));

            float global_p90 = 120.0f;
            if (!all_np.empty()) {
                const size_t idx = static_cast<size_t>(all_np.size() * 0.90f);
                std::nth_element(all_np.begin(), all_np.begin() + idx, all_np.end());
                global_p90 = std::max(1.0f, all_np[idx]);
            }

            int n_genus_cal = 0, n_fallback = 0;
            for (const auto& [genus, accs] : pass_a.genus_members) {
                std::vector<float> np;
                for (const auto& acc : accs) {
                    auto it = quality.find(acc);
                    if (it != quality.end() && it->second.marker_n_expected > 0
                            && it->second.marker_n_present > 0)
                        np.push_back(static_cast<float>(it->second.marker_n_present));
                }

                float denom = global_p90;
                if (np.size() >= 10) {
                    const size_t idx = static_cast<size_t>(np.size() * 0.90f);
                    std::nth_element(np.begin(), np.begin() + idx, np.end());
                    denom = std::max(1.0f, np[idx]);
                    ++n_genus_cal;
                } else {
                    ++n_fallback;
                }

                for (const auto& acc : accs) {
                    auto it = quality.find(acc);
                    if (it != quality.end() && it->second.marker_n_expected > 0)
                        it->second.marker_completeness =
                            std::min(1.0f, static_cast<float>(it->second.marker_n_present) / denom);
                }
            }
            spdlog::info("check: marker_completeness calibrated: {} genera (genus p90), "
                         "{} small/singleton genera (global p90={:.0f})",
                         n_genus_cal, n_fallback, global_p90);
        }

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
    tsv << "accession\tquality_tier\tcompleteness_effective\tcompleteness_cluster_relative\tcompleteness_fragmentation"
           "\tcompleteness_post_decontam\tcontamination_leakage\tcontamination_tnf_excess"
           "\tchromosome_skew_closure\tleakage_residual\tself_coherence"
           "\tchargaff_parity\tspectral_gap\tscale_kink\tcontamination_mixture\tmixture_sources"
           "\tn_mix_windows\tfiedler_value\tcontamination_contig_outlier\tcontamination_spe"
           "\tcontamination_sibling_outlier\tcontamination_rho_outlier\tcontamination_cross_genus\tcontamination_contig_split"
           "\tcontamination_self_outlier\tfiedler_oph_split"
           "\tfiedler_tnf_bimod\tfiedler_tnf_gap"
           "\tfmh_minority_fraction"
           "\tmarker_completeness\tmarker_redundancy\tmarker_redundancy_z\tmarker_joint_contamination\tmarker_n_present\tmarker_n_expected"
           "\tmarker_present_hex"
           "\tqual_flags\tsupport_tier\tinterval_width\n";

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
        const float comp_eff = (!std::isnan(q.marker_completeness) && q.marker_completeness > q.completeness_cluster_relative)
                                   ? q.marker_completeness
                                   : q.completeness_cluster_relative;
        const float tnf = std::isnan(q.contamination_tnf_excess) ? 0.0f : q.contamination_tnf_excess;
        const char* qtier = "LQ";
        if (!std::isnan(comp_eff) && comp_eff >= 0.90f && tnf < 0.05f) qtier = "HQ";
        else if (!std::isnan(comp_eff) && comp_eff >= 0.50f && tnf < 0.20f) qtier = "MQ";
        tsv << acc << '\t'
            << qtier << '\t'
            << fmt(comp_eff) << '\t'
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
            << fmt(q.contamination_sibling_outlier) << '\t'
            << fmt(q.contamination_rho_outlier) << '\t'
            << fmt(q.contamination_cross_genus) << '\t'
            << fmt(q.contamination_contig_split) << '\t'
            << fmt(q.contamination_self_outlier) << '\t'
            << fmt(q.fiedler_oph_split) << '\t'
            << fmt(q.fiedler_tnf_bimod) << '\t'
            << fmt(q.fiedler_tnf_gap) << '\t'
            << fmt(q.fmh_minority_fraction) << '\t'
            << fmt(q.marker_completeness) << '\t'
            << fmt(q.marker_redundancy) << '\t'
            << fmt(q.marker_redundancy_z) << '\t'
            << fmt(q.marker_joint_contamination) << '\t'
            << q.marker_n_present << '\t'
            << q.marker_n_expected << '\t';
        // Per-SCG presence bitmask as lowercase hex string (empty when not scored).
        if (!q.marker_present_bits.empty()) {
            static constexpr char hex[] = "0123456789abcdef";
            for (uint8_t b : q.marker_present_bits) {
                tsv << hex[b >> 4] << hex[b & 0xf];
            }
        }
        tsv << '\t'
            << static_cast<int>(q.qual_flags) << '\t'
            << tier_str(q.support_tier) << '\t'
            << q.interval_width << '\n';
    }
    spdlog::info("check: TSV written to {}", output.string());
    return 0;
}

} // namespace genopack::check
