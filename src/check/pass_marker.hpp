#pragma once
#include "types.hpp"
#include "pass_a.hpp"
#include "pack_iface.hpp"
#include <string>
#include <unordered_map>

namespace genopack::check {

struct PassMarkerConfig {
    std::string markers_path;
    int         marker_min_hits = 5;
};

// Unconditional marker completeness pass: scores every genome with a valid genus
// assignment regardless of TNF excess or Pass B flagging.  Fills
// GenomeQuality::marker_completeness, marker_redundancy, marker_n_present,
// marker_n_expected, and marker_present_bits for any genome where those fields
// are still NaN/empty after Pass B.
//
// This is the CheckM2-equivalent completeness estimate — same bac120/ar53
// marker set, Dayhoff-6 k=12 syncmer detection, no GCOV dependency.
void run_pass_marker(ICheckReader&                                  pack,
                     const PassAResult&                             pass_a,
                     std::unordered_map<std::string, GenomeQuality>& quality,
                     const PassMarkerConfig&                        cfg,
                     int                                            threads);

} // namespace genopack::check
