#pragma once
#include <filesystem>
#include <vector>

namespace genopack {

// Merge N .gpk archives into a single output archive.
// Shard bytes are copied verbatim (no recompression).
// Genome IDs and shard IDs are renumbered globally to avoid conflicts.
// CATL and ACCX are rebuilt from all inputs.
// Does not deduplicate — run 'genopack dedup' after merging if needed.
void merge_archives(const std::vector<std::filesystem::path>& inputs,
                    const std::filesystem::path& output);

} // namespace genopack
