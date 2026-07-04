#include "run_check.hpp"
#include "pass_a.hpp"
#include "pass_aamer.hpp"
#include "pass_b.hpp"
#include "pack_iface.hpp"
#include <genopack/archive.hpp>
#include <genopack/qual.hpp>
#include <genopack/qual_columns.hpp>
#include <genopack/qcontig.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/section_checksum.hpp>
#include <genopack/toc.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <future>
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

// Continuous genome quality score in [0,1] for threshold-free ranking, replacing the
// discrete HQ/MQ/LQ cutoffs. Multiplicative because completeness, contamination, and
// contig-coherence (fiedler) are near-independent failure modes — a complete-but-chimeric
// bin is still a poor representative, and additive mixing would let high completeness mask
// contamination. NaN completeness → 0.5 neutral prior; NaN contamination/fiedler → no penalty.
//   comp_eff already encodes the cr_unreliable insight, so sparse-reference MAGs are not crushed.
//   cont_factor: rational decay centred on 0.05 (=0.5 there, =0.2 at 0.10) — a smooth analog of
//     the old HQ/MQ contamination cliffs, so ranking stays stable across them.
//   coh_penalty: bounded ≤50% cut, engaged only for HIGH fiedler_split (>0.5, = chimera/bin-split);
//     capped rather than a hard gate because the fiedler threshold is uncalibrated (no test coverage).
//     NB: high fiedler_oph_split = poorly-connected contig graph = split; low = coherent (score_bin.cpp).
float compute_quality_score(float comp_eff, float cont_val, float fiedler_split) {
    const float C = std::isnan(comp_eff) ? 0.5f : comp_eff;
    const float X = std::isnan(cont_val) ? 0.0f : cont_val;
    const float s = std::isnan(fiedler_split) ? 0.0f : fiedler_split;
    const float xr = X / 0.05f;
    const float cont_factor = 1.0f / (1.0f + xr * xr);
    float ramp = (s - 0.5f) / 0.4f;
    ramp = ramp < 0.0f ? 0.0f : (ramp > 1.0f ? 1.0f : ramp);
    const float coh_penalty = 1.0f - 0.5f * ramp;
    return C * cont_factor * coh_penalty;
}

// Single contamination determinant: NA-safe max over every observable contamination axis.
// HQ used to hinge on the FMH minority axis alone, so a missing (NA) FMH call was scored as
// zero contamination and a genuinely dirty genome slipped through HQ (GCMeta_00156701: fmh
// NA, CheckM2 12.2%). Taking the max over all axes and ignoring NA means no single absent
// axis can lower the score. All axes are [0,1] minority/outlier fractions (mixture util.hpp:81,
// spe/rho = u8/255 pass_b.cpp:196-198, duplication/contig_outlier_adj fractions types.hpp:70,81)
// except tnf_excess (dist/ref-1 ratio, builder.cpp:1415) which is clamped into the same range.
// cross_genus is deliberately excluded — it is a hard chimera veto handled at the call site.
//   contamination_tnf_minor is the near-clade detector (util.cpp: BIC-free TNF-GMM minority
//   mass, multi-contig guarded) — catches same-species/genus mergers that stay k-mer-similar
//   so fmh/cross_genus read ~0 but a compositional minority still exists.
//   observed: FMH/mixture/tnf_minor/contig_outlier_adj/spe/rho/duplication carry a NAN sentinel
//     when unmeasured, so any non-NAN one proves contamination was actually observed. leakage and
//     tnf_excess default to 0.0f when unmeasured (indistinguishable from clean), so they raise
//     the max but never by themselves establish observability. When observed stays false, no
//     axis was available → contamination is unconfirmed and HQ must be capped at MQ.
float contamination_aggregate(const GenomeQuality& q, bool& observed) {
    float m = 0.0f;
    observed = false;
    auto clamp01 = [](float v) { return v < 0.0f ? 0.0f : (v > 1.0f ? 1.0f : v); };
    auto acc = [&](float v) {                         // NA-capable axes: establish observability
        if (std::isnan(v)) return;
        observed = true;
        const float c = clamp01(v);
        if (c > m) m = c;
    };
    auto acc_max = [&](float v) {                     // default-0 axes: raise max only
        if (std::isnan(v)) return;
        const float c = clamp01(v);
        if (c > m) m = c;
    };
    acc(q.fmh_minority_fraction);
    acc(q.contamination_mixture);
    acc(q.contamination_tnf_minor);
    acc(q.contamination_contig_outlier_adj);
    acc(q.contamination_spe);
    acc(q.contamination_rho_outlier);
    acc(q.contamination_duplication);
    acc_max(q.contamination_leakage);
    acc_max(q.contamination_tnf_excess);
    return m;
}

// Cross-genus chimera veto: a bin whose contigs fit a foreign genus better than the assigned
// one is the most dangerous defect for a reference DB. Kept independent of the aggregate (a
// chimera can look clean on every other axis, e.g. ZAHE052 cross_genus=1.0). >=0.10 of scored
// bp cross-assigned = chimera, above the incidental-HGT range.
bool cross_genus_chimera(const GenomeQuality& q) {
    return !std::isnan(q.contamination_cross_genus) && q.contamination_cross_genus >= 0.10f;
}

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
                           const std::unordered_map<std::string, GenomeId>& acc_to_id,
                           uint64_t core_model_hash,
                           uint64_t fcore_model_hash)
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

    std::vector<QualRecord> recs;
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
        // Build-time SCC core-duplication (encode_dup maps NaN -> 0 sentinel = not scored).
        // Must be re-persisted here: a zero-initialized QualRecord would otherwise wipe the
        // build value on every QCOL write-back, flipping a contaminated genome back to HQ.
        r.contamination_duplication_u16 = QualRecord::encode_dup(q.contamination_duplication);
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
        {
            const float cf = std::min(1.0f, std::isnan(q.contamination_contig_outlier) ? 0.0f : q.contamination_contig_outlier);
            r.contig_outlier_u8   = static_cast<uint8_t> (cf * 255.0f);
            r.contig_outlier_u16  = static_cast<uint16_t>(cf * 65535.0f + 0.5f);
        }
        r.spe_outlier_u8                = static_cast<uint8_t>(
            std::min(1.0f, std::isnan(q.contamination_spe) ? 0.0f : q.contamination_spe) * 255.0f);
        r.sibling_outlier_u8            = static_cast<uint8_t>(
            std::min(1.0f, std::isnan(q.contamination_sibling_outlier) ? 0.0f : q.contamination_sibling_outlier) * 255.0f);
        r.rho_outlier_u8                = static_cast<uint8_t>(
            std::min(1.0f, std::isnan(q.contamination_rho_outlier) ? 0.0f : q.contamination_rho_outlier) * 255.0f);
        r.cross_genus_u8                = static_cast<uint8_t>(
            std::min(1.0f, std::isnan(q.contamination_cross_genus) ? 0.0f : q.contamination_cross_genus) * 255.0f);
        // sketch_fill: 0=not_scored, 1-255 = clamp(value,0,1.27)*200 (200=100%)
        r.sketch_fill_u8 = std::isnan(q.completeness_sketch_fill) ? 0u
            : static_cast<uint8_t>(std::clamp(q.completeness_sketch_fill, 0.0f, 1.27f) * 200.0f + 0.5f);
        {
            const float fmhf = std::min(1.0f, std::isnan(q.fmh_minority_fraction) ? 0.0f : q.fmh_minority_fraction);
            r.fmh_minority_u8  = static_cast<uint8_t> (fmhf * 255.0f);   // legacy: kept for old readers
            r.fmh_minority_u16 = static_cast<uint16_t>(fmhf * 65535.0f + 0.5f);
        }
        // marker_completeness: 0=not_scored, 1-255=(value-1)/254 → range [0,1] with sentinel 0
        if (!std::isnan(q.marker_completeness))
            r.marker_completeness_u8 = static_cast<uint8_t>(
                std::clamp(q.marker_completeness, 0.0f, 1.0f) * 254.0f + 1.0f);
        else
            r.marker_completeness_u8 = 0;

        r.chargaff_parity               = q.chargaff_parity;
        r.spectral_gap                  = q.spectral_gap;
        r.scale_kink                    = q.scale_kink;
        r.completeness_aamer_core        = q.completeness_aamer_core;
        r.completeness_aamer_family_core = q.completeness_aamer_family_core;

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
        // Encode the quality tier (mirrors TSV output logic).
        {
            const bool aamer_only = std::isnan(q.completeness_post_decontam) &&
                                    std::isnan(q.completeness_cluster_relative) &&
                                    !std::isnan(q.completeness_aamer_core);
            const float _pd = q.completeness_post_decontam;
            const float _cr = q.completeness_cluster_relative;
            const float _ac = q.completeness_aamer_core;
            // cr_unreliable: cluster_relative near zero means centroid sparsity in a diverse
            // genus, not low completeness. Safe to bypass geomean and use pd directly when:
            //   (a) classic: pd≥0.90 + ac≥0.90 + cr<0.50 (original guard), OR
            //   (b) degenerate cr<0.05 + all contamination axes quiet — a clean genome at
            //       cr=0.02 is centroid-geometry noise, not misassignment; but if ANY
            //       contamination signal fires, low cr IS a misbin fingerprint and must stand, OR
            //   (c) no marker-core data (ac NaN) + contamination quiet — MAGs from undersampled
            //       biomes (marine/soil/wastewater) have sparse reference clusters; a low cr there
            //       reflects a bad comparator, not incomplete sequence. Without ac to fall back on,
            //       geomean(pd,cr) mislabeled ~1.2M otherwise-clean MAGs as LQ in the r232_v5 audit.
            const bool cont_quiet_cr =
                (std::isnan(q.fmh_minority_fraction)        || q.fmh_minority_fraction        < 0.10f) &&
                (std::isnan(q.contamination_mixture)        || q.contamination_mixture        < 0.10f) &&
                (std::isnan(q.contamination_leakage)        || q.contamination_leakage        < 0.02f) &&
                (std::isnan(q.contamination_contig_outlier) || q.contamination_contig_outlier < 0.05f) &&
                (std::isnan(q.contamination_contig_split)   || q.contamination_contig_split   < 0.05f);
            const bool cr_unreliable = !std::isnan(_pd) && _pd >= 0.90f
                                    && !std::isnan(_cr) && _cr < 0.50f
                                    && ((!std::isnan(_ac) && _ac >= 0.90f)
                                        || (_cr < 0.05f && cont_quiet_cr)
                                        || (std::isnan(_ac) && cont_quiet_cr));
            const float ce = cr_unreliable ? _pd
                           : (!std::isnan(_pd) && !std::isnan(_cr) && std::fabs(_pd - _cr) > 0.30f)
                           ? std::sqrt(_pd * _cr)
                           : (!std::isnan(_pd) ? _pd : (!std::isnan(_cr) ? _cr : _ac));
            const float fmh = !std::isnan(q.fmh_minority_fraction) ? q.fmh_minority_fraction : NAN;
            const bool leakage_fires = (q.contamination_leakage >= 0.02f);
            const int nsig =
                (!std::isnan(fmh)                          && fmh                          >= 0.10f) +
                (!std::isnan(q.contamination_mixture)      && q.contamination_mixture      >= 0.10f) +
                (!std::isnan(q.contamination_contig_outlier) && q.contamination_contig_outlier >= 0.05f) +
                (!std::isnan(q.contamination_contig_split) && q.contamination_contig_split >= 0.05f) +
                leakage_fires;
            bool cont_observed = false;
            const float cv = contamination_aggregate(q, cont_observed);
            const char* t = "LQ";
            if (!std::isnan(ce) && ce >= 0.90f && cv < 0.05f) t = "HQ";
            else if (!std::isnan(ce) && ce >= 0.50f && cv < 0.10f) t = "MQ";
            // fiedler_oph_split has no validated threshold or test coverage (audited 2026-07-02);
            // require at least one corroborating contamination/split signal before trusting it
            // alone to override an otherwise-clean comp/cont verdict.
            if (!std::isnan(q.fiedler_oph_split) && q.fiedler_oph_split < 0.1f && nsig >= 1) t = "LQ";
            if (aamer_only && t[0] == 'H') t = "MQ";
            // Cross-genus chimera: hard veto, independent of the aggregate.
            if (cross_genus_chimera(q) && t[0] != 'L') t = "LQ";
            // HQ requires a positive, observed clean contamination axis. If nothing was
            // observable, contamination is unconfirmed → cap at MQ, never HQ.
            if (t[0] == 'H' && !cont_observed) t = "MQ";
            r.quality_tier_u8 = QualRecord::encode_qtier(t);
            r.set_quality_score(compute_quality_score(ce, cv, q.fiedler_oph_split));
        }
        recs.push_back(r);
    }

    AppendWriter writer;
    writer.open_append(gpk_path);

    // AAMER_GENUS_CORE intrinsic completeness → its own provenance-carrying column
    // (references core_model_hash; absent for genomes without a genus core).
    std::unordered_map<uint64_t, float> aamer_core_cov;
    for (const auto& [acc, q] : quality) {
        if (std::isnan(q.completeness_aamer_core)) continue;
        auto it = acc_to_id.find(acc);
        if (it == acc_to_id.end()) continue;
        aamer_core_cov.emplace(it->second, q.completeness_aamer_core);
    }
    // AAMER_FAMILY_CORE intrinsic completeness (genus-core fallback) → its own column.
    std::unordered_map<uint64_t, float> aamer_family_cov;
    for (const auto& [acc, q] : quality) {
        if (std::isnan(q.completeness_aamer_family_core)) continue;
        auto it = acc_to_id.find(acc);
        if (it == acc_to_id.end()) continue;
        aamer_family_cov.emplace(it->second, q.completeness_aamer_family_core);
    }
    QcolExtraColumns extra;
    extra.core_model_hash = core_model_hash;
    if (!aamer_core_cov.empty()) extra.aamer_genus_core = &aamer_core_cov;
    extra.family_core_model_hash = fcore_model_hash;
    if (!aamer_family_cov.empty()) extra.aamer_family_core = &aamer_family_cov;

    // SEC_QCONTIG — per-contig overlay: persist pass_b's per-contig flags so a
    // consumer can see which contigs drive each genome's contamination score.
    QcontigWriter qcw;
    std::unordered_set<uint64_t> qcw_seen;   // dedup aliased accessions → one genome_id
    for (const auto& [acc, q] : quality) {
        if (q.contig_flags.empty()) continue;
        auto it = acc_to_id.find(acc);
        if (it == acc_to_id.end()) continue;
        if (!qcw_seen.insert(static_cast<uint64_t>(it->second)).second) continue;
        qcw.add_genome(static_cast<uint64_t>(it->second), q.contig_flags);
    }

    uint64_t section_id = toc.next_section_id();
    SectionDesc qcol_sd = qcol_write(writer, section_id, std::move(recs), extra);
    SectionDesc qcontig_sd = qcw.finalize(writer, section_id + 1);

    TocWriter new_toc;
    for (const auto& sd : toc.sections) {
        if (sd.type != SEC_QUAL && sd.type != SEC_QCOL && sd.type != SEC_QCONTIG)
            new_toc.add_section(sd);
    }
    new_toc.add_section(qcol_sd);
    if (qcontig_sd.item_count > 0) new_toc.add_section(qcontig_sd);

    // Content-checksum the appended QUAL (and any still-unstamped section) so
    // verify can validate it.
    writer.flush();
    stamp_section_checksums(gpk_path.c_str(), new_toc.sections(), /*only_if_zero=*/true);

    new_toc.finalize(writer,
                     generation,
                     toc.header.live_genome_count,
                     toc.header.total_genome_count,
                     prev_toc_offset,
                     toc.header.catalog_root_section_id,
                     toc.header.accession_root_section_id,
                     toc.header.tombstone_root_section_id);

    spdlog::info("check: wrote {} QCOL records + {} QCONTIG rows to {}",
                 qcol_sd.item_count, qcontig_sd.item_count, gpk_path.string());
}

} // namespace

int cmd_check(const std::filesystem::path& pack_path,
              const std::filesystem::path& genomes_file,
              int threads,
              int min_genus_size,
              float leakage_threshold,
              const std::filesystem::path& output,
              bool recompute,
              const std::filesystem::path& markers_path,
              bool scan_all)
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
            ar.scan_genome_accessions([&](std::string_view acc, GenomeId gid) {
                if (!ar.is_deleted(gid) && filter.count(std::string(acc)))
                    accessions.emplace_back(acc);
            });
        } else {
            ar.scan_genome_accessions([&](std::string_view acc, GenomeId gid) {
                if (!ar.is_deleted(gid)) accessions.emplace_back(acc);
            });
        }

        if (accessions.empty()) {
            spdlog::debug("check: no accessions in {}, skipping", gp.string());
            continue;
        }
        spdlog::info("check: {} accessions in {}", accessions.size(), gp.filename().string());

        // Load existing QUAL section.
        // Used to: (a) skip SKCH scan when !recompute, (b) compute CCO baseline always.
        std::unordered_map<uint64_t, QualRecord> qual_cache;
        const std::unordered_map<uint64_t, QualRecord>* qual_cache_ptr = nullptr;
        if (ar.has_qual()) {
            ar.scan_qual([&](const QualRecord& r) {
                qual_cache.emplace(r.genome_id, r);
            });
            if (!recompute) {
                spdlog::info("check: loaded {} cached QUAL records (use --recompute to force rescan)",
                             qual_cache.size());
                qual_cache_ptr = &qual_cache;
            }
        }

        // Preload GMI concurrently with pass-A: SEC_GAMI is independent of flagged-genome results.
        genopack::GlobalMultiplicityIndex preloaded_gmi;
        std::future<void> gmi_future;
        if (pack.has_gami_v2())
            gmi_future = std::async(std::launch::async,
                                    [&]{ pack.load_gami_into(preloaded_gmi); });

        auto pass_a = run_pass_a(pack, accessions, threads, min_genus_size,
                                 leakage_threshold, /*min_gstx_members=*/50,
                                 qual_cache_ptr);
        auto quality = pass_a.quality;

        run_pass_aamer(pack, pass_a, quality, threads);

        PassBConfig pb_cfg;
        pb_cfg.markers_path = markers_path.string();   // empty = marker scoring disabled
        pb_cfg.scan_all     = scan_all;
        const auto* baseline_ptr = qual_cache.empty() ? nullptr : &qual_cache;
        run_pass_b(pack, pass_a, quality, pb_cfg, threads, qual_cache_ptr, baseline_ptr,
                   pack.has_gami_v2() ? &preloaded_gmi : nullptr,
                   gmi_future.valid()  ? &gmi_future   : nullptr);

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

        // core_dup is a build-time-only signal: the per-genus SCC set and ceiling live in
        // the .csp panel, not the pack, so check cannot recompute it. Preserve the build
        // value from the QUAL cache in ALL modes — pass_a skips it under --recompute (cache
        // ptr nulled) and a fresh scan never fills it, which otherwise drops the axis to NA
        // and lets a same-species-contaminated genome climb back to HQ under --scan-all.
        for (auto& [acc, q] : quality) {
            if (!std::isnan(q.contamination_duplication)) continue;
            auto idit = acc_to_id.find(acc);
            if (idit == acc_to_id.end()) continue;
            auto cit = qual_cache.find(idit->second);
            if (cit != qual_cache.end())
                q.contamination_duplication =
                    QualRecord::decode_dup(cit->second.contamination_duplication_u16);
        }

        const uint64_t core_model_hash =
            ar.core_reader() ? ar.core_reader()->core_model_hash() : 0;
        const uint64_t fcore_model_hash =
            ar.fcore_reader() ? ar.fcore_reader()->core_model_hash() : 0;
        ar.close();

        // Contamination-axis provenance: when the FMH minority axis is absent
        // (no FMHR section / not computed -> NaN), the quality tier falls back to
        // the sketch-leakage axis. Flag it per genome so the axis used is
        // auditable in both the on-disk QUAL record and the TSV.
        for (auto& [acc, q] : quality)
            if (std::isnan(q.fmh_minority_fraction))
                q.qual_flags |= QualRecord::QUAL_FLAG_FMH_AXIS_ABSENT;

        write_qual_to_archive(gp, quality, acc_to_id, core_model_hash, fcore_model_hash);

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
    tsv << "accession\tquality_tier\tquality_score\tcompleteness_effective\tcompleteness_cluster_relative\tcompleteness_sketch_fill\tcompleteness_fragmentation"
           "\tcompleteness_post_decontam"
           "\tcompleteness_aamer_core"
           "\tcompleteness_aamer_family_core"
           "\tfmh_contamination"
           "\tcontamination_leakage\tcontamination_tnf_excess"
           "\tchromosome_skew_closure\tleakage_residual\tself_coherence"
           "\tchargaff_parity\tspectral_gap\tscale_kink\tcontamination_mixture\tmixture_sources"
           "\tn_mix_windows\tcontamination_contig_outlier\tcontamination_spe"
           "\tcontamination_rho_outlier\tcontamination_cross_genus\tcontamination_contig_split"
           "\tcontamination_duplication\tmarker_present_bits"
           "\tqual_flags\tsupport_tier\tinterval_width"
           "\tcontamination_contig_outlier_adj\tassembly_tier\n";

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

    // TODO(D1): Replace fixed thresholds with a trained ensemble scorer.
    // Wire comp_eff, cont_val, fiedler_oph_split, and marker signals into
    // EnsembleScorer::predict() and replace the qtier chain below when ready.
    struct EnsembleScorer {
        float w_comp = 0.0f, w_cont = 0.0f, w_fiedler = 0.0f; // zero → no-op
        const char* predict(float /*comp*/, float /*cont*/, float /*fiedler*/) const {
            return nullptr; // nullptr → fall through to rule-based assignment
        }
    };

    for (const auto& [acc, q] : all_quality) {
        const bool aamer_only = std::isnan(q.completeness_post_decontam) &&
                                std::isnan(q.completeness_cluster_relative) &&
                                !std::isnan(q.completeness_aamer_core);
        const float _pd2 = q.completeness_post_decontam;
        const float _cr2 = q.completeness_cluster_relative;
        const float _ac2 = q.completeness_aamer_core;
        const bool cont_quiet_cr2 =
            (std::isnan(q.fmh_minority_fraction)        || q.fmh_minority_fraction        < 0.10f) &&
            (std::isnan(q.contamination_mixture)        || q.contamination_mixture        < 0.10f) &&
            (std::isnan(q.contamination_leakage)        || q.contamination_leakage        < 0.02f) &&
            (std::isnan(q.contamination_contig_outlier) || q.contamination_contig_outlier < 0.05f) &&
            (std::isnan(q.contamination_contig_split)   || q.contamination_contig_split   < 0.05f);
        const bool cr2_unreliable = !std::isnan(_pd2) && _pd2 >= 0.90f
                                 && !std::isnan(_cr2) && _cr2 < 0.50f
                                 && ((!std::isnan(_ac2) && _ac2 >= 0.90f)
                                     || (_cr2 < 0.05f && cont_quiet_cr2)
                                     || (std::isnan(_ac2) && cont_quiet_cr2));
        const float comp_eff = cr2_unreliable ? _pd2
                             : (!std::isnan(_pd2) && !std::isnan(_cr2) && std::fabs(_pd2 - _cr2) > 0.30f)
                             ? std::sqrt(_pd2 * _cr2)
                             : (!std::isnan(_pd2) ? _pd2 : (!std::isnan(_cr2) ? _cr2 : _ac2));
        const float fmh_cont = !std::isnan(q.fmh_minority_fraction) ? q.fmh_minority_fraction : NAN;
        // Count independent contamination signals that fire. FMH alone is insufficient because
        // it captures HGT/mobile elements as "foreign" k-mers (leakage disagrees for genuine HGT).
        // LQ on contamination grounds requires >=2 orthogonal signals to fire simultaneously,
        // AND leakage must co-fire before FMH is used as primary proxy (prevents FMH+CCO-only FPs).
        const bool leak_fires = (q.contamination_leakage >= 0.02f);
        const int cont_signals =
            (!std::isnan(fmh_cont)                         && fmh_cont                         >= 0.10f) +
            (!std::isnan(q.contamination_mixture)          && q.contamination_mixture          >= 0.10f) +
            (!std::isnan(q.contamination_contig_outlier)   && q.contamination_contig_outlier   >= 0.05f) +
            (!std::isnan(q.contamination_contig_split)     && q.contamination_contig_split     >= 0.05f) +
            leak_fires;
        // TODO(D1): const EnsembleScorer ens; if (auto t = ens.predict(...)) qtier = t;
        // NA-safe contamination aggregate (mirrors QCOL): max over all observable axes.
        bool cont_observed = false;
        const float cont_val = contamination_aggregate(q, cont_observed);
        const char* qtier = "LQ";
        if (!std::isnan(comp_eff) && comp_eff >= 0.90f && cont_val < 0.05f) qtier = "HQ";
        else if (!std::isnan(comp_eff) && comp_eff >= 0.50f && cont_val < 0.10f) qtier = "MQ";
        // fiedler_oph_split has no validated threshold or test coverage (audited 2026-07-02);
        // require at least one corroborating contamination/split signal before trusting it
        // alone to override an otherwise-clean comp/cont verdict.
        if (!std::isnan(q.fiedler_oph_split) && q.fiedler_oph_split < 0.1f && cont_signals >= 1) qtier = "LQ";
        if (aamer_only && qtier[0] == 'H') qtier = "MQ";
        // Cross-genus chimera: hard veto, independent of the aggregate.
        if (cross_genus_chimera(q) && qtier[0] != 'L') qtier = "LQ";
        // HQ requires a positive, observed clean contamination axis; if none was observable,
        // contamination is unconfirmed → cap at MQ, never HQ.
        if (qtier[0] == 'H' && !cont_observed) qtier = "MQ";
        const float quality_score = compute_quality_score(comp_eff, cont_val, q.fiedler_oph_split);
        tsv << acc << '\t'
            << qtier << '\t'
            << fmt(quality_score) << '\t'
            << fmt(comp_eff) << '\t'
            << fmt(q.completeness_cluster_relative) << '\t'
            << fmt(q.completeness_sketch_fill) << '\t'
            << fmt(q.completeness_fragmentation) << '\t'
            << fmt(q.completeness_post_decontam) << '\t'
            << fmt(q.completeness_aamer_core) << '\t'
            << fmt(q.completeness_aamer_family_core) << '\t'
            << fmt(fmh_cont) << '\t'
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
            << fmt(q.contamination_contig_outlier) << '\t'
            << fmt(q.contamination_spe) << '\t'
            << fmt(q.contamination_rho_outlier) << '\t'
            << fmt(q.contamination_cross_genus) << '\t'
            << fmt(q.contamination_contig_split) << '\t'
            << fmt(q.contamination_duplication) << '\t';
        // Per-SCG presence bitmask as lowercase hex string (empty when not scored).
        if (!q.marker_present_bits.empty()) {
            static constexpr char hex[] = "0123456789abcdef";
            for (uint8_t b : q.marker_present_bits) {
                tsv << hex[b >> 4] << hex[b & 0xf];
            }
        }
        // Assembly quality tier using completeness_effective + adjusted CCO.
        const float cont_axis = !std::isnan(q.contamination_contig_outlier_adj)
                                    ? q.contamination_contig_outlier_adj
                                    : (std::isnan(q.contamination_contig_outlier) ? 0.0f
                                                                                  : q.contamination_contig_outlier);
        const char* assembly_tier = "LQ";
        if (!std::isnan(comp_eff) && comp_eff >= 0.90f && cont_axis < 0.05f) assembly_tier = "HQ";
        else if (!std::isnan(comp_eff) && comp_eff >= 0.50f && cont_axis < 0.10f) assembly_tier = "MQ";

        tsv << '\t'
            << static_cast<int>(q.qual_flags) << '\t'
            << tier_str(q.support_tier) << '\t'
            << q.interval_width << '\t'
            << fmt(q.contamination_contig_outlier_adj) << '\t'
            << assembly_tier << '\n';
    }
    spdlog::info("check: TSV written to {}", output.string());
    return 0;
}

} // namespace genopack::check
