#pragma once
#include <filesystem>
#include <vector>

namespace genopack {

// Merge N .gpk archives into a single output archive.
// Shard bytes are copied verbatim (no recompression).
// CATL and ACCX are rebuilt from all inputs.
// Does not deduplicate — run 'genopack dedup' after merging if needed.
//
// remap_genome_ids: when true (default), genome IDs and shard IDs are
//   renumbered globally. When false, assumes all inputs already have
//   globally unique, non-overlapping genome_id ranges (e.g. from a
//   parallel build that used per-part starting_genome_id).
void merge_archives(const std::vector<std::filesystem::path>& inputs,
                    const std::filesystem::path& output,
                    bool remap_genome_ids = true,
                    bool merge_cidx = true);

} // namespace genopack
