#pragma once
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>

namespace genopack::check {

enum class SupportTier : uint8_t {
    GenusSaturated = 0,
    GenusSparse    = 1,
    Singleton      = 2,
};

struct ConsistencyVector {
    float kmer_anchor_agreement = NAN;
    float containment_k_ratio   = NAN;
    float skew_closure          = NAN;

    float norm() const {
        float d = 0.0f;
        if (!std::isnan(kmer_anchor_agreement)) d += (1.0f - kmer_anchor_agreement) * (1.0f - kmer_anchor_agreement);
        if (!std::isnan(containment_k_ratio))   d += (1.0f - containment_k_ratio)   * (1.0f - containment_k_ratio);
        if (!std::isnan(skew_closure))           d += (1.0f - skew_closure)           * (1.0f - skew_closure);
        return std::sqrt(d);
    }
};

struct ContigFlag {
    uint32_t contig_offset = 0;
    uint32_t contig_length = 0;
    float    tnf_score     = NAN;  // centroid-distance ratio (coarse; set for contigs >=5kb)
    float    leakage_score = NAN;
    float    gcov_t2_pct   = NAN;  // GCOV Hotelling-T² percentile (within-genus Mahalanobis)
    float    gcov_spe_pct  = NAN;  // GCOV SPE percentile (residual outside genus subspace; catches cross-genus)
    // Protein-aamer foreign-containment channel (works on SMALL contigs where TNF fails).
    float    prot_foreign_frac = NAN;  // foreign_specific / (host_specific + foreign_specific)
    float    prot_loglr        = NAN;  // log((foreign_specific+1)/(host_specific+1))  [v1 presence-only]
    uint16_t prot_classifiable    = 0; // host_specific + foreign_specific aamers
    uint16_t prot_host_specific   = 0;
    uint16_t prot_foreign_specific= 0;
    uint16_t prot_family_hits     = 0;
    uint32_t prot_best_genus      = 0; // 22-bit tag of the dominant foreign genus
    uint8_t  prot_flags           = 0; // bit0 ABSTAIN_LOW_N, bit1 MOBILE_NATIVE, bit2 NOVEL_GENUS_FALLBACK
};
static constexpr uint8_t PROT_ABSTAIN_LOW_N      = 0x1;
static constexpr uint8_t PROT_MOBILE_NATIVE      = 0x2;
static constexpr uint8_t PROT_NOVEL_GENUS_FALLBACK = 0x4;

struct GenomeQuality {
    float completeness_cluster_relative = NAN;
    // Phase-1 relative-conspecific-containment estimators (OPH-containment, NOT pangenome-union).
    // Computed only in the no-GSTX path when the genus is GenusSaturated with ≥10 members; NaN otherwise (abstain).
    float accessory_ratio               = NAN;  // c0_query / median(c0_all)
    float accessory_z                   = NAN;  // (c0_query - median) / (1.4826 * MAD)
    float completeness_post_decontam    = NAN;
    float completeness_aamer_core        = NAN;  // genus prevalence-core coverage (intrinsic; CheckM2-aligned)
    float completeness_aamer_family_core = NAN;  // family prevalence-core coverage (genus-core fallback; SEC_FCORE)
    float contamination_leakage    = 0.0f;
    float contamination_tnf_excess = 0.0f;
    float contamination_tnf_minor = NAN;  // near-clade TNF-GMM minority mass (BIC-free, multi-contig guarded)
    float    contamination_contig_outlier = NAN;
    float    contamination_contig_outlier_adj = NAN; // CCO minus per-genus clean-genome baseline
    float    contamination_spe              = NAN;
    float    contamination_rho_outlier      = NAN;  // ρ* Mahalanobis outlier fraction
    float    contamination_contig_split     = NAN;  // per-contig containment split minority fraction
    float    fmh_minority_fraction         = NAN;  // FMH k=21,c=125 minority bp / scored_bp
    float    contamination_cross_genus     = NAN;  // fraction bp where any foreign genus fits better than assigned
    float    contamination_duplication     = NAN;  // .csp diagnostic-aamer duplication fraction (build-time; CheckM2-aligned)
    float    contamination_core_dup_mass   = NAN;  // Phase-1: non-saturating SCC dup mass Σ(c-1)/Σc (build/score-time only)

    // Marker-panel completeness (protein k-mer based; NAN if no .mrk file supplied)
    float    marker_completeness         = NAN;  // present / expected markers
    int      marker_n_present            = -1;   // raw present count
    int      marker_n_expected           = -1;   // expected marker count for this lineage
    // Per-SCG presence bitmask: bit i set if marker i had ≥1 contig vote.
    // Max 173 markers (120 bac + 53 arc); stored as packed bytes (22 bytes covers 176 bits).
    // Empty when marker scoring not run.
    std::vector<uint8_t> marker_present_bits;

    uint8_t  qual_flags           = 0;
    float    sketch_breadth       = NAN; // NAN in check path (requires build-time consensus sketch)
    std::vector<ContigFlag> contig_flags;
    SupportTier support_tier   = SupportTier::GenusSaturated;
};

} // namespace genopack::check
