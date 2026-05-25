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
    float    tnf_score     = NAN;
    float    leakage_score = NAN;
};

struct GenomeQuality {
    float completeness_cluster_relative = NAN;
    float completeness_fragmentation    = NAN;
    float completeness_post_decontam    = NAN;
    float contamination_leakage    = 0.0f;
    float contamination_tnf_excess = 0.0f;
    float chromosome_skew_closure = NAN;
    float leakage_residual        = NAN;
    float self_coherence          = NAN;
    float chargaff_parity         = NAN;
    float spectral_gap            = NAN;
    float scale_kink              = NAN;
    float contamination_mixture   = NAN;
    int   mixture_sources         = 1;
    uint16_t n_mix_windows        = 0;
    float    fiedler_value              = NAN;
    float    contamination_contig_outlier = NAN;
    float    contamination_spe           = NAN;

    uint8_t  qual_flags           = 0;
    float    sketch_breadth       = NAN; // NAN in check path (requires build-time consensus sketch)
    std::vector<ContigFlag> contig_flags;
    SupportTier support_tier   = SupportTier::GenusSaturated;
    float       interval_width = 1.0f;
};

} // namespace genopack::check
