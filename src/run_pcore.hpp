#pragma once
#include <filesystem>

namespace genopack {

// `genopack pcore` — build SEC_PCORE: the unified prevalence-annotated per-genus
// aamer reference (every aamer in ≥1 member + u8 prevalence percent). Dense enough
// to detect small foreign contigs (the conserved-only CORE is too sparse for
// 1-2kb). Reuses CORE's k/min_seg_aa/frac_max so query extraction matches.
// members_list (optional): restrict the members to a reference accession list.
// bootstrap_* are used only when an input archive has no SEC_CORE section to
// source the build config from (SEC_CORE has no in-repo builder; see docs/aamer.md).
// Canonical values: min_seg_aa = k = 8, theta = 0.90, frac_max_hash = UINT64_MAX
// (dense, no FracMinHash subsampling — PCORE must be dense to catch 1–2 kb contigs).
int cmd_pcore(const std::filesystem::path& pack_path, int threads,
              const std::filesystem::path& members_list,
              float bootstrap_theta = 0.90f,
              uint32_t bootstrap_min_seg_aa = 8,
              uint64_t bootstrap_frac_max_hash = UINT64_MAX);

} // namespace genopack
