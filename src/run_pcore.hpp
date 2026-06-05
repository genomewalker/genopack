#pragma once
#include <filesystem>

namespace genopack {

// `genopack pcore` — build SEC_PCORE: the unified prevalence-annotated per-genus
// aamer reference (every aamer in ≥1 member + u8 prevalence percent). Dense enough
// to detect small foreign contigs (the conserved-only CORE is too sparse for
// 1-2kb). Reuses CORE's k/min_seg_aa/frac_max so query extraction matches.
// members_list (optional): restrict the members to a reference accession list.
int cmd_pcore(const std::filesystem::path& pack_path, int threads,
              const std::filesystem::path& members_list);

} // namespace genopack
