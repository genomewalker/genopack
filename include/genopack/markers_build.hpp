#pragma once
#include "markers.hpp"
#include <filesystem>

namespace genopack {

struct MarkersBuildConfig {
    std::filesystem::path gtdbtk_db;   // path to GTDB-Tk reference database root
    std::filesystem::path output;      // output .mrk file
    int threads = 1;
    float default_threshold = 0.30f;   // min containment fraction → marker present
    uint16_t frac_scale = 1;           // FracMinHash divisor: keep hashes ≤ UINT64_MAX/frac_scale
};

// Build a markers.mrk sidecar from GTDB-Tk r232 MSA files.
// Reads:  <gtdbtk_db>/msa/gtdb_r232_bac120.faa
//         <gtdbtk_db>/msa/gtdb_r232_ar53.faa
//         <gtdbtk_db>/taxonomy/bac120_taxonomy_r232_reps.tsv
//         <gtdbtk_db>/taxonomy/ar53_taxonomy_r232_reps.tsv
// Writes: cfg.output (.mrk file)
void build_markers_panel(const MarkersBuildConfig& cfg);

} // namespace genopack
