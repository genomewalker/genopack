#pragma once
// ── genopack tier — IDF tier table build tools ──────────────────────────────
// Phase 2: merge .ptier side-channel files into a global .gtier table.
// Phase 3: stamp .gtier tiers into a PCORE v1 section producing PCORE v2.

#include <filesystem>
#include <string>
#include <vector>

namespace genopack {

// genopack tier merge -i p0.ptier p1.ptier ... -o global.gtier
// Merges N per-part .ptier files into one global .gtier IDF table.
// All input files must share the same fmh_seed, frac_max_hash, and k.
int cmd_tier_merge(const std::vector<std::string>& ptier_files,
                   const std::string& out_path);

// genopack tier stamp -i part.gpk --table global.gtier [-o out.gpk]
// Reads V1 PCORE from part.gpk, looks up per-aamer tier in global.gtier,
// and rewrites the PCORE section as V2 (tier bytes appended after prevq_c).
// If -o is omitted, writes to <input>.tier.gpk.
int cmd_tier_stamp(const std::string& gpk_path,
                   const std::string& gtier_path,
                   const std::string& out_path = {});

} // namespace genopack
