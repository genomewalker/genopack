#pragma once
#include "markers.hpp"
#include <filesystem>

namespace genopack {

struct MarkersBuildConfig {
    std::filesystem::path gtdbtk_db;   // path to GTDB-Tk reference database root
    std::filesystem::path output;      // output .mrk file
    int threads = 1;
    float default_threshold = 0.30f;   // min containment fraction → marker present
    float expected_min_frac = 0.95f;   // marker counted as expected only if detectable (≥1 IC syncmer)
                                       // in ≥this fraction of genus reference genomes — mirrors
                                       // CheckM2's ~97% single-copy universality criterion
    uint16_t frac_scale = 1;           // FracMinHash divisor (ignored when dayhoff6=true)
    bool dayhoff6 = false;             // build Dayhoff-6 k=12 syncmer profile pool
    float ic_threshold = 0.25f;        // min per-column IC (0–1) to include a k-mer position
};

// Build a markers.mrk sidecar from GTDB-Tk r232 MSA files.
// Reads:  <gtdbtk_db>/msa/gtdb_r232_bac120.faa
//         <gtdbtk_db>/msa/gtdb_r232_ar53.faa
//         <gtdbtk_db>/taxonomy/bac120_taxonomy_r232_reps.tsv
//         <gtdbtk_db>/taxonomy/ar53_taxonomy_r232_reps.tsv
// Writes: cfg.output (.mrk file)
void build_markers_panel(const MarkersBuildConfig& cfg);

} // namespace genopack
