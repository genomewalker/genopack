#pragma once
#include "types.hpp"
#include "pass_a.hpp"
#include "pack_iface.hpp"
#include <genopack/qual.hpp>
#include <cstdint>
#include <string>
#include <unordered_map>

namespace genopack::check {

struct PassBConfig {
    float contamination_flag_threshold = 0.02f;
    float tnf_flag_threshold           = 0.02f;  // lowered from 0.10 — fires at ~2% phylum spike
    float completeness_flag_threshold  = 0.70f;
    float contig_tnf_threshold         = 3.0f;
    float contig_leakage_threshold     = 0.20f;
    float gcov_outlier_percentile      = 0.99f; // GCOV-based outlier threshold (percentile of calibration dist)
    float cross_genus_lr_margin        = 0.0f;  // log-LR margin for cross_genus/sibling; 0 = flag if any foreign LL > host LL

    // Marker panel scoring (protein k-mer completeness via .mrk sidecar)
    std::string markers_path;       // path to markers.mrk; empty = disabled
    int         marker_min_hits = 5; // min pool-hit count to call a marker present (pre-calibration)
};

void run_pass_b(ICheckReader& pack,
                const PassAResult& pass_a,
                std::unordered_map<std::string, GenomeQuality>& quality,
                const PassBConfig& cfg,
                int threads,
                const std::unordered_map<uint64_t, QualRecord>* qual_cache = nullptr);

} // namespace genopack::check
