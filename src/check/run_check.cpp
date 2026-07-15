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

// Continuous genome quality score in [0,1] for threshold-free ranking. COMPLETENESS-ONLY:
// contamination no longer folds into this score. The old multiplicative cont_factor was a
// back-door gate -- a 0.10 knee halves the score on a signal whose PPV vs CheckM2 never exceeds
// ~15%, silently reordering the catalog on noise. Contamination is now reported as three
// separate channels (contamination_channels) and consumed by derep as an explicit D->S->G
// tiebreak, so the ranking stays interpretable and the noisy channels can never override the
// calibrated one. NaN completeness → 0.5 neutral prior.
float compute_quality_score(float comp_eff) {
    return std::isnan(comp_eff) ? 0.5f : comp_eff;
}

// Intrinsic completeness_effective for tier/score, kept identical between the QCOL
// (quality_tier_u8) and TSV paths. INTRINSIC completeness only: fraction-tracking
// marker_completeness (present/expected SCG, genus-calibrated) is the primary signal, with
// aamer_core (CheckM2-aligned but presence-saturating) then post_decontam (bp-retention after
// the contig contamination scan) as fallbacks when marker_completeness is unavailable.
//   completeness_cluster_relative is NOT intrinsic completeness — it is the fraction of the
//   genus PANGENOME a genome covers. A finished isolate in a diverse genus recovers ~100% of
//   its own core but only a small slice of the genus accessory pangenome, so a low cluster_relative
//   there is genus diversity, not missing sequence. It must never by itself drive comp_eff down.
//   It is admitted only as a SOFT CORROBORATOR: when the intrinsic signal is ALSO low
//   (genuine agreement that the genome is partial) AND contamination is quiet, geomean(intrinsic, cr)
//   lets a real strain-partial genome settle below the intrinsic estimate. When intrinsic
//   completeness is high, cr is ignored (its low value is diversity, not incompleteness).
float completeness_effective(const GenomeQuality& q) {
    const float mc = q.marker_completeness;
    const float ac = q.completeness_aamer_core;
    const float pd = q.completeness_post_decontam;
    const float cr = q.completeness_cluster_relative;
    // Intrinsic estimate priority: fraction-tracking marker_completeness (present/expected SCG,
    // genus-calibrated; slope ≈0.49 vs CheckM2's 0.86 on the fragmentation panel) first. Fall back
    // to aamer_core when marker_completeness is NaN — it is only populated on the pass-B route, so
    // QCOL-cache genomes (marker_completeness_u8==0 → NaN) must still get a sane estimate. aamer_core
    // saturates (0.047 dynamic range: 1.00/0.97/0.95 at true 1.0/0.7/0.5) but never NaN when scored.
    // post_decontam (bp-retention) is the last resort.
    const float intrinsic = !std::isnan(mc) ? mc
                          : (!std::isnan(ac) ? ac : pd);
    if (std::isnan(intrinsic))
        return cr;  // no intrinsic signal at all → cluster_relative is the only completeness proxy
    // Soft corroboration: only when intrinsic ALSO reads genuinely partial (below the MQ line)
    // does cr pull it down. A genome with a supported intrinsic (≥0.50, MQ-complete) keeps its
    // estimate regardless of a low (diverse-genus) cr — a low cr there is diversity, not missing
    // sequence, and must never geomean-crush an MQ genome below the MQ threshold.
    if (intrinsic < 0.50f && !std::isnan(cr) && cr < intrinsic && (intrinsic - cr) > 0.30f)
        return std::sqrt(intrinsic * cr);
    return intrinsic;
}

// Duplication contamination in CheckM2 units. Phase-2: prefer the non-saturating SCC dup
// MASS (contamination_core_dup_mass = Σ(c-1)/Σc, cladesplit.cpp:610) over the saturating
// EXCESS. The mass tracks true contamination linearly on the spike panel: measured means
// 0.031/0.056/0.136/0.198 at true CheckM2 0.05/0.10/0.20/0.30. A least-squares line
// (slope 1.449, intercept 0.010) maps mass → CheckM2 fraction with residuals
// 0.055/0.091/0.207/0.297 — a calibrated graded term replacing the ad-hoc ÷8 hack.
//   When core_dup_mass is unavailable (old-stride pack / no SCC set → NaN), fall back to
// the legacy excess ÷8 so v3 packs and unscored genera keep a duplication signal.
constexpr float kDupMassSlope     = 1.448954f;  // spike-panel OLS: core_dup_mass → CheckM2 contam
constexpr float kDupMassIntercept = 0.009998f;
constexpr float kDupToContamScale = 0.125f;     // legacy ÷8 fallback for excess when mass is NaN
float duplication_contamination(float excess, float mass) {
    if (!std::isnan(mass)) {
        const float c = kDupMassSlope * mass + kDupMassIntercept;
        return c < 0.0f ? 0.0f : c;
    }
    return std::isnan(excess) ? NAN : excess * kDupToContamScale;
}

// Reference-free contamination as THREE near-independent channels, replacing the old NA-safe MAX
// over 7 axes. That MAX let the single noisiest of ~7 CORRELATED axes set the score: the four
// geometry axes (rho/spe/contig_outlier_adj/tnf_minor) co-fire at 20-50x independence -- they are
// ONE TNF/GCOV signal -- so a MAX counted that signal up to four times. The calibration (all
// 9.53M genomes, CheckM2 as truth) showed the resulting gate at 6.9% PPV / 5.7% recall, demoting
// 362,935 genomes to catch ~1,200 real contaminants. So:
//   D = duplication_contamination(...)  -- the ONLY CheckM2-calibrated channel (marker duplication)
//   S = fmh_minority_fraction           -- sketch minority (k-mer)
//   G = median(present geometry axes)   -- the correlated geometry signal, collapsed to one vote
// leakage is DROPPED (mathematically dead: max 0.005 << its 0.02 threshold; the alpha=1 forcing
// is not worth refitting) and tnf_excess is DROPPED (default-0, untrusted, only ever entered
// through the MAX). cross_genus stays out -- it is a ranker (PPV 2.8%), never a gate.
// Contamination is REPORTED, never a discard gate (see quality_tier): a genome is LQ on
// completeness alone, and these channels drive derep as a reliability-ordered tiebreak, so
// genuine heterogeneity (prophage/plasmid/HGT island -- real signal CheckM2 rightly calls clean)
// is preserved as evidence instead of deleted. `observed` is true iff any channel was measured.
struct ContamChannels {
    float D = NAN, S = NAN, G = NAN;
    float display = NAN;   // noisy-OR union of present channels; TSV/human display ONLY
    bool  observed = false;
    int   geom_fired = 0;  // geometry axes over their per-axis threshold; channel counts iff >=2
};

ContamChannels contamination_channels(const GenomeQuality& q) {
    auto clamp01 = [](float v) { return v < 0.0f ? 0.0f : (v > 1.0f ? 1.0f : v); };
    ContamChannels c;

    const float d = duplication_contamination(q.contamination_duplication, q.contamination_core_dup_mass);
    if (!std::isnan(d)) { c.D = clamp01(d); c.observed = true; }
    if (!std::isnan(q.fmh_minority_fraction)) { c.S = clamp01(q.fmh_minority_fraction); c.observed = true; }

    // Geometry: median of the present axes so the one signal votes once, not four times. Per-axis
    // fire thresholds (contig_outlier_adj 0.05, others 0.10) accumulate geom_fired; the channel is
    // only trusted to have fired when >=2 of the 4 agree -- D5 in the calibration (require two
    // channels to corroborate) was the WORST row, so a lone geometry axis must never carry weight.
    const std::pair<float, float> geom[] = {
        {q.contamination_rho_outlier,        0.10f},
        {q.contamination_spe,                0.10f},
        {q.contamination_contig_outlier_adj, 0.05f},
        {q.contamination_tnf_minor,          0.10f},
    };
    std::vector<float> present;
    present.reserve(4);
    for (const auto& [v, thr] : geom) {
        if (std::isnan(v)) continue;
        present.push_back(clamp01(v));
        if (v >= thr) ++c.geom_fired;
    }
    if (!present.empty()) {
        std::sort(present.begin(), present.end());
        const size_t n = present.size();
        c.G = (n & 1) ? present[n / 2] : 0.5f * (present[n / 2 - 1] + present[n / 2]);
        c.observed = true;
    }

    // Display scalar: probabilistic union (noisy-OR) over PRESENT channels, for TSV/human reading
    // only. Derep must use the D->S->G lexicographic priority, NEVER this scalar: noisy-OR inflates
    // several weak channels into false corroboration (three 0.2s -> 0.49), and RMS would dilute the
    // one calibrated channel against the silence of the others. Neither is a safe ordering key.
    float keep = 1.0f; bool any = false;
    for (float ch : {c.D, c.S, c.G})
        if (!std::isnan(ch)) { keep *= (1.0f - ch); any = true; }
    if (any) c.display = 1.0f - keep;
    return c;
}

// Diagnostic count of channels that fired (TSV only). Geometry counts as one channel and only
// when >=2 of its axes agree. Not a tier input -- tiering is completeness-only (quality_tier).
int contam_channels_fired(const ContamChannels& c) {
    int n = 0;
    if (!std::isnan(c.D) && c.D >= 0.05f) ++n;
    if (!std::isnan(c.S) && c.S >= 0.10f) ++n;
    if (c.geom_fired >= 2) ++n;
    return n;
}

// The tier decision, shared by the QCOL (quality_tier_u8) and TSV paths so they can never drift.
// COMPLETENESS-ONLY: contamination is no longer grounds to demote to LQ. The old gate demoted
// 362,935 genomes at 6.9% PPV; decoupling removes every one of those wrongful discards and leaves
// only the genuine-incompleteness floor (ce<0.50). Contamination survives one narrow role: the
// single CheckM2-CALIBRATED channel (duplication, D) may cap HQ->MQ, so the top tier still means
// "clean and complete". This cap NEVER forces LQ, so it cannot reintroduce wrongful discard -- it
// only distinguishes HQ from MQ, and only on the one axis validated against CheckM2.
const char* quality_tier(const GenomeQuality& q) {
    const float ce = completeness_effective(q);
    if (std::isnan(ce)) return "LQ";
    if (ce >= 0.90f) {
        // aamer_core is presence-saturating (0.047 dynamic range); it cannot support HQ alone.
        const bool aamer_only = std::isnan(q.completeness_post_decontam) &&
                                std::isnan(q.completeness_cluster_relative) &&
                                !std::isnan(q.completeness_aamer_core);
        if (aamer_only) return "MQ";
        const float D = duplication_contamination(q.contamination_duplication,
                                                  q.contamination_core_dup_mass);
        if (!std::isnan(D) && D >= 0.05f) return "MQ";
        return "HQ";
    }
    if (ce >= 0.50f) return "MQ";
    return "LQ";
}

// cross_genus is a RANKER, not a gate. GTDB-Tk ground truth on a 2,036-genome stratified
// sample (900 top-bin, 450 veto, 300 mid, 450 clean controls): it orders misplacement
// correctly -- 23x enrichment over the clean baseline, monotonic in every bin -- but its
// positive predictive value as a demotion rule is 2.8% at >=0.10 and only 20.4% even at a
// saturated 1.0. It was demoting ~97 correctly-labelled genomes for every 3 misplaced ones.
//
// The reason is in pass_b.cpp:864-895: a contig is flagged when ANY of ~52 candidate foreign
// Gaussians beats the host log-likelihood by more than cross_genus_lr_margin, which defaulted
// to 0 -- an uncorrected maximum over 52 competitors. Whether a foreign model out-scores the
// host is a GENOME-level property (the genome sits slightly off-centre inside its own
// legitimate genus), so the winner beats the host on essentially every contig at once and the
// statistic collapses to 0-or-1: 70% of the top bin sits at exactly 1.0000.
//
// So: no hard veto. Set a calibrated cross_genus_lr_margin to correct the multiple comparison,
// and let the >=2-trusted-axes gate decide, like every other axis.

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

// Re-injection path for the build-time-only dup axes. Accepts either a genopack quality TSV
// (contamination_duplication / core_dup_mass) or a `cladesplit score` TSV (core_dup /
// core_dup_mass) -- both carry the axis keyed by accession, so a pack whose QUAL lost it can be
// repaired without a full rebuild.
std::unordered_map<std::string, std::pair<float, float>>
read_dup_restore(const std::filesystem::path& tsv)
{
    std::ifstream in(tsv);
    if (!in) throw std::runtime_error("check: cannot open --dup-restore: " + tsv.string());

    std::string line;
    if (!std::getline(in, line))
        throw std::runtime_error("check: --dup-restore is empty: " + tsv.string());

    int c_acc = -1, c_dup = -1, c_mass = -1, col = 0;
    for (size_t p = 0; p <= line.size(); ) {
        const size_t e = line.find('\t', p);
        const std::string h = line.substr(p, (e == std::string::npos ? line.size() : e) - p);
        if      (h == "accession")                 c_acc  = col;
        else if (h == "contamination_duplication" ||
                 h == "core_dup")                  c_dup  = col;
        else if (h == "core_dup_mass" ||
                 h == "contamination_core_dup_mass") c_mass = col;
        ++col;
        if (e == std::string::npos) break;
        p = e + 1;
    }
    if (c_acc < 0 || (c_dup < 0 && c_mass < 0))
        throw std::runtime_error("check: --dup-restore needs an 'accession' column plus at least "
                                 "one of core_dup/contamination_duplication/core_dup_mass: " +
                                 tsv.string());

    auto parse = [](const std::string& s) -> float {
        if (s.empty() || s == "NA" || s == "nan") return std::numeric_limits<float>::quiet_NaN();
        try { return std::stof(s); }
        catch (...) { return std::numeric_limits<float>::quiet_NaN(); }
    };

    std::unordered_map<std::string, std::pair<float, float>> out;
    const int need = std::max({c_acc, c_dup, c_mass});
    std::vector<std::string> f;
    while (std::getline(in, line)) {
        f.clear();
        for (size_t p = 0; p <= line.size(); ) {
            const size_t e = line.find('\t', p);
            f.push_back(line.substr(p, (e == std::string::npos ? line.size() : e) - p));
            if (e == std::string::npos) break;
            p = e + 1;
        }
        if (static_cast<int>(f.size()) <= need) continue;
        const float d = c_dup  >= 0 ? parse(f[c_dup])  : std::numeric_limits<float>::quiet_NaN();
        const float m = c_mass >= 0 ? parse(f[c_mass]) : std::numeric_limits<float>::quiet_NaN();
        if (std::isnan(d) && std::isnan(m)) continue;
        out.emplace(f[c_acc], std::make_pair(d, m));
    }
    spdlog::info("check: --dup-restore loaded {} genomes with a dup axis from {}",
                 out.size(), tsv.string());
    return out;
}

// `prior` is the QUAL section this write supersedes (nullptr if the pack had none). A QUAL
// section is a REPLACEMENT, not a delta: the reader takes the newest one only. So a run that
// scores a subset (-g) and writes only what it scored silently deletes every genome it did not
// touch -- a 2k-genome spot-check over a 953k pack destroyed 951k records exactly this way, and
// because core_dup is build-time-only (see the restore block below), a later --recompute could
// not rebuild them and wrote NaN over the primary contamination estimator. Two invariants now
// make that unrepresentable: (1) records in `prior` but not in `quality` are carried forward, so
// the section can only grow; (2) a write that would lower core_dup_mass coverage on genomes the
// cache already had is refused outright rather than silently persisted.
void write_qual_to_archive(const std::filesystem::path& gpk_path,
                           const std::unordered_map<std::string, GenomeQuality>& quality,
                           const std::unordered_map<std::string, GenomeId>& acc_to_id,
                           uint64_t core_model_hash,
                           const std::unordered_map<uint64_t, QualRecord>* prior)
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
        r.completeness_post_decontam    = q.completeness_post_decontam;
        r.contamination_leakage         = q.contamination_leakage;
        r.contamination_tnf_excess      = q.contamination_tnf_excess;
        // Build-time SCC core-duplication (encode_dup maps NaN -> 0 sentinel = not scored).
        // Must be re-persisted here: a zero-initialized QualRecord would otherwise wipe the
        // build value on every QCOL write-back, flipping a contaminated genome back to HQ.
        r.contamination_duplication_u16 = QualRecord::encode_dup(q.contamination_duplication);
        r.support_tier                  = static_cast<uint8_t>(q.support_tier);
        {
            const float cf = std::min(1.0f, std::isnan(q.contamination_contig_outlier) ? 0.0f : q.contamination_contig_outlier);
            r.contig_outlier_u8   = static_cast<uint8_t> (cf * 255.0f);
            r.contig_outlier_u16  = static_cast<uint16_t>(cf * 65535.0f + 0.5f);
        }
        r.spe_outlier_u8                = static_cast<uint8_t>(
            std::min(1.0f, std::isnan(q.contamination_spe) ? 0.0f : q.contamination_spe) * 255.0f);
        r.rho_outlier_u8                = static_cast<uint8_t>(
            std::min(1.0f, std::isnan(q.contamination_rho_outlier) ? 0.0f : q.contamination_rho_outlier) * 255.0f);
        r.cross_genus_u8                = static_cast<uint8_t>(
            std::min(1.0f, std::isnan(q.contamination_cross_genus) ? 0.0f : q.contamination_cross_genus) * 255.0f);
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
        r.completeness_aamer_core        = q.completeness_aamer_core;
        // Phase-2 estimators re-persisted so a QCOL write-back keeps the build/check value
        // (raw float; NAN = not scored). Without this a zero-init record would wipe them.
        r.contamination_core_dup_mass    = q.contamination_core_dup_mass;
        r.accessory_ratio                = q.accessory_ratio;
        // qual_flags carry GCOV_SCORED (set in pass-B) / FMH_AXIS_ABSENT provenance.
        r.qual_flags                    = q.qual_flags;
        // Tier (completeness-only) and score (completeness-only), shared with the TSV path.
        r.quality_tier_u8 = QualRecord::encode_qtier(quality_tier(q));
        r.set_quality_score(compute_quality_score(completeness_effective(q)));
        recs.push_back(r);
    }

    if (prior) {
        std::unordered_map<uint64_t, size_t> written;
        written.reserve(recs.size());
        for (size_t i = 0; i < recs.size(); ++i) written.emplace(recs[i].genome_id, i);

        size_t carried = 0, dup_lost = 0;
        for (const auto& [gid, pr] : *prior) {
            auto it = written.find(gid);
            if (it == written.end()) { recs.push_back(pr); ++carried; continue; }
            if (!std::isnan(pr.contamination_core_dup_mass) &&
                std::isnan(recs[it->second].contamination_core_dup_mass)) ++dup_lost;
        }
        if (carried)
            spdlog::info("check: carried {} cached QUAL records not rescored in this run", carried);
        if (dup_lost)
            throw std::runtime_error(
                "check: refusing to write QUAL -- would drop build-time core_dup_mass on " +
                std::to_string(dup_lost) + " genomes that already had it. core_dup is "
                "build-time-only (SCC set lives in the .csp panel, not the pack) and cannot be "
                "recomputed here. Re-run with --markers, or re-inject the axis with "
                "--dup-restore <quality-or-cladesplit-score.tsv>.");
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
    QcolExtraColumns extra;
    extra.core_model_hash = core_model_hash;
    if (!aamer_core_cov.empty()) extra.aamer_genus_core = &aamer_core_cov;

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
              bool scan_all,
              float cross_genus_margin,
              const std::filesystem::path& dup_restore_path)
{
    auto gpk_paths = collect_gpk_paths(pack_path);
    if (gpk_paths.empty())
        throw std::runtime_error("check: no .gpk files found at " + pack_path.string());

    std::unordered_map<std::string, std::pair<float, float>> dup_restore;
    if (!dup_restore_path.empty()) dup_restore = read_dup_restore(dup_restore_path);

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
        pb_cfg.cross_genus_lr_margin = cross_genus_margin;
        const auto* baseline_ptr = qual_cache.empty() ? nullptr : &qual_cache;
        run_pass_b(pack, pass_a, quality, pb_cfg, threads, qual_cache_ptr, baseline_ptr,
                   pack.has_gami_v2() ? &preloaded_gmi : nullptr,
                   gmi_future.valid()  ? &gmi_future   : nullptr);

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
        // Both dup axes are build-time-only: the legacy u16 contamination_duplication and
        // the non-saturating contamination_core_dup_mass (Σ(c-1)/Σc, cladesplit). Without
        // restoring the latter, --recompute would rewrite QUAL with core_dup_mass=NaN and
        // silently destroy the primary graded contamination estimator.
        size_t dup_from_cache = 0, dup_from_tsv = 0, dup_missing = 0;
        for (auto& [acc, q] : quality) {
            bool need_dup      = std::isnan(q.contamination_duplication);
            bool need_dup_mass = std::isnan(q.contamination_core_dup_mass);
            if (!need_dup && !need_dup_mass) continue;

            auto idit = acc_to_id.find(acc);
            if (idit != acc_to_id.end()) {
                auto cit = qual_cache.find(idit->second);
                if (cit != qual_cache.end()) {
                    if (need_dup) {
                        const float d = QualRecord::decode_dup(cit->second.contamination_duplication_u16);
                        if (!std::isnan(d)) { q.contamination_duplication = d; need_dup = false; }
                    }
                    if (need_dup_mass && !std::isnan(cit->second.contamination_core_dup_mass)) {
                        q.contamination_core_dup_mass = cit->second.contamination_core_dup_mass;
                        need_dup_mass = false;
                    }
                    if (!need_dup || !need_dup_mass) ++dup_from_cache;
                }
            }

            // Cache lost the axis (a prior subset/no-marker run overwrote QUAL): fall back to the
            // external table. This is the only way back short of a full rebuild.
            if ((need_dup || need_dup_mass) && !dup_restore.empty()) {
                auto rit = dup_restore.find(acc);
                if (rit != dup_restore.end()) {
                    if (need_dup && !std::isnan(rit->second.first)) {
                        q.contamination_duplication = rit->second.first;
                        need_dup = false;
                    }
                    if (need_dup_mass && !std::isnan(rit->second.second)) {
                        q.contamination_core_dup_mass = rit->second.second;
                        need_dup_mass = false;
                    }
                    ++dup_from_tsv;
                }
            }
            if (need_dup && need_dup_mass) ++dup_missing;
        }
        spdlog::info("check: dup axis restored from cache={} tsv={} still-absent={}",
                     dup_from_cache, dup_from_tsv, dup_missing);

        const uint64_t core_model_hash =
            ar.core_reader() ? ar.core_reader()->core_model_hash() : 0;
        ar.close();

        // Contamination-axis provenance: when the FMH minority axis is absent
        // (no FMHR section / not computed -> NaN), the quality tier falls back to
        // the sketch-leakage axis. Flag it per genome so the axis used is
        // auditable in both the on-disk QUAL record and the TSV.
        for (auto& [acc, q] : quality)
            if (std::isnan(q.fmh_minority_fraction))
                q.qual_flags |= QualRecord::QUAL_FLAG_FMH_AXIS_ABSENT;

        write_qual_to_archive(gp, quality, acc_to_id, core_model_hash,
                              qual_cache.empty() ? nullptr : &qual_cache);

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
    tsv << "accession\tquality_tier\tquality_score\tcompleteness_effective\tcompleteness_cluster_relative"
           "\tcompleteness_post_decontam"
           "\tcompleteness_aamer_core"
           "\tfmh_contamination"
           "\tcontamination_leakage\tcontamination_tnf_excess"
           "\tcontamination_contig_outlier\tcontamination_spe"
           "\tcontamination_rho_outlier\tcontamination_cross_genus\tcontamination_contig_split"
           "\tcontamination_duplication\tmarker_present_bits"
           "\tqual_flags\tsupport_tier"
           "\tcontamination_contig_outlier_adj\tassembly_tier"
           // pangenome_fraction: honest alias of completeness_cluster_relative — the fraction of
           // the genus pangenome (union aamer content across genus members) this genome covers,
           // i.e. accessory-genome breadth / pangenome representativeness. NOT genome completeness
           // and deliberately NOT an input to quality_score/quality_tier. Interpretable only when
           // support_tier == GenusSaturated; in a Sparse/Singleton genus the denominator is tiny
           // so the fraction saturates near 1 trivially.
           "\tpangenome_fraction"
           // Phase-1 estimators (in-memory only, no QUAL/format bump): non-saturating SCC
           // duplication mass (build/score-time; NA in pure check where the .csp panel isn't
           // loaded) and relative-conspecific OPH-containment ratio/z (no-GSTX GenusSaturated only).
           "\tcore_dup_mass\taccessory_ratio\taccessory_z"
           // Phase-1 completeness ESTIMATOR (in-memory + TSV only, no QUAL/format bump).
           // completeness_marker = present/expected single-copy markers — the CheckM2-style
           // fraction that tracks fraction-of-genome (declines with fragmentation), unlike the
           // presence-saturating completeness_aamer_core. marker_n_present/expected are the raw
           // counts behind it. NA when no --markers panel or no genus marker calibration.
           "\tcompleteness_marker\tmarker_n_present\tmarker_n_expected"
           // Contamination is reported, not gated. The three near-independent channels feed
           // derep's D->S->G lexicographic tiebreak: contam_D (duplication, the only
           // CheckM2-calibrated channel), contam_S (fmh sketch), contam_G (median of the present
           // correlated geometry axes -- one vote, not four). contam_score is the noisy-OR union
           // for DISPLAY ONLY; channels_fired counts how many of the three fired (geometry only
           // when >=2 of its axes agree). tnf_minor is emitted raw so contam_G is reconstructible.
           "\tcontamination_tnf_minor\tcontam_D\tcontam_S\tcontam_G\tcontam_score\tchannels_fired\n";

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
        const float comp_eff = completeness_effective(q);
        const float fmh_cont = !std::isnan(q.fmh_minority_fraction) ? q.fmh_minority_fraction : NAN;
        // Contamination is reported (three channels), not gated: tier is completeness-only and
        // identical to the QCOL path via the shared helper.
        const ContamChannels cc = contamination_channels(q);
        const char* qtier = quality_tier(q);
        const float quality_score = compute_quality_score(comp_eff);
        tsv << acc << '\t'
            << qtier << '\t'
            << fmt(quality_score) << '\t'
            << fmt(comp_eff) << '\t'
            << fmt(q.completeness_cluster_relative) << '\t'
            << fmt(q.completeness_post_decontam) << '\t'
            << fmt(q.completeness_aamer_core) << '\t'
            << fmt(fmh_cont) << '\t'
            << fmt(q.contamination_leakage) << '\t'
            << fmt(q.contamination_tnf_excess) << '\t'
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
            << fmt(q.contamination_contig_outlier_adj) << '\t'
            << assembly_tier << '\t'
            << fmt(q.completeness_cluster_relative) << '\t'
            << fmt(q.contamination_core_dup_mass) << '\t'
            << fmt(q.accessory_ratio) << '\t'
            << fmt(q.accessory_z) << '\t'
            << fmt(q.marker_completeness) << '\t'
            << (q.marker_n_expected >= 0 ? std::to_string(q.marker_n_present) : std::string("NA")) << '\t'
            << (q.marker_n_expected >= 0 ? std::to_string(q.marker_n_expected) : std::string("NA")) << '\t'
            << fmt(q.contamination_tnf_minor) << '\t'
            << fmt(cc.D) << '\t'
            << fmt(cc.S) << '\t'
            << fmt(cc.G) << '\t'
            << fmt(cc.display) << '\t'
            << contam_channels_fired(cc) << '\n';
    }
    spdlog::info("check: TSV written to {}", output.string());
    return 0;
}

} // namespace genopack::check
