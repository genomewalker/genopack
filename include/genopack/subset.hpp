#pragma once
#include "shard.hpp"
#include <filesystem>
#include <string>
#include <unordered_set>

namespace genopack {

struct SubsetConfig {
    std::unordered_set<std::string> accessions;   // target accession set (populated by caller)
    ShardWriterConfig               shard_cfg;
    size_t                          threads            = 1;
    uint64_t                        max_bucket_bytes   = 32ULL << 30;
    bool                            verbose            = false;
};

// Read accessions from a plain-text file (one per line) or a TSV file whose
// header contains a column named "accession" or "name" (case-insensitive).
// Auto-detection: if the first line contains a tab it is treated as TSV.
std::unordered_set<std::string> load_accession_set(const std::filesystem::path& path);

// Copy only the genomes whose accession appears in cfg.accessions from
// input_gpk into a new archive at output_gpk.
void subset_archive(const std::filesystem::path& input_gpk,
                    const std::filesystem::path& output_gpk,
                    const SubsetConfig& cfg);

} // namespace genopack
