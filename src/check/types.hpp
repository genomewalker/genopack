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
    float completeness_sketch_fill      = NAN;
    float completeness_fragmentation    = NAN;
    float completeness_post_decontam    = NAN;
    float completeness_aamer_core        = NAN;  // genus prevalence-core coverage (intrinsic; CheckM2-aligned)
    float completeness_aamer_family_core = NAN;  // family prevalence-core coverage (genus-core fallback; SEC_FCORE)
    float contamination_leakage    = 0.0f;
    float contamination_tnf_excess = 0.0f;
    float chromosome_skew_closure = NAN;
    float leakage_residual        = NAN;
    float self_coherence          = NAN;
    float chargaff_parity         = NAN;
    float spectral_gap            = NAN;
    float scale_kink              = NAN;
    float contamination_mixture   = NAN;
    float contamination_tnf_minor = NAN;  // near-clade TNF-GMM minority mass (BIC-free, multi-contig guarded)
    int   mixture_sources         = 1;
    uint16_t n_mix_windows        = 0;
    float    fiedler_value              = NAN;
    float    contamination_contig_outlier = NAN;
    float    contamination_contig_outlier_adj = NAN; // CCO minus per-genus clean-genome baseline
    float    contamination_spe              = NAN;
    float    contamination_sibling_outlier  = NAN;
    float    contamination_rho_outlier      = NAN;  // ρ* Mahalanobis outlier fraction
    float    contamination_contig_split     = NAN;  // per-contig containment split minority fraction
    float    contamination_self_outlier    = NAN;  // self-containment z-score outlier fraction
    float    fiedler_oph_split             = NAN;  // sketch Fiedler 1−λ₂/2 of pairwise OPH Jaccard
    float    fiedler_tnf_bimod             = NAN;  // TNF-kernel v₂ normalized max-gap bimodality
    float    fiedler_tnf_gap               = NAN;  // TNF-kernel λ₃−λ₂ eigengap
    float    fmh_minority_fraction         = NAN;  // FMH k=21,c=125 minority bp / scored_bp
    float    contamination_cross_genus     = NAN;  // fraction bp where any foreign genus fits better than assigned
    float    contamination_duplication     = NAN;  // .csp diagnostic-aamer duplication fraction (build-time; CheckM2-aligned)

    // Marker-panel completeness (protein k-mer based; NAN if no .mrk file supplied)
    float    marker_completeness         = NAN;  // present / expected markers
    float    marker_redundancy           = NAN;  // markers with ≥2 contig votes / expected
    float    marker_redundancy_z         = NAN;  // (observed - genus_median) / (1.4826 * MAD_eff)
    float    marker_joint_contamination  = NAN;  // markers with both native+cross_genus votes / expected
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
    float       interval_width = 1.0f;
};

} // namespace genopack::check
