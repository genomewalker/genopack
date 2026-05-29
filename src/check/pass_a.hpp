#pragma once
#include "types.hpp"
#include "pack_iface.hpp"
#include <genopack/qual.hpp>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace genopack::check {

struct PassAResult {
    std::vector<std::string>                       accessions;
    std::unordered_map<std::string, GenomeQuality> quality;
    std::unordered_map<std::string, std::vector<std::string>> genus_members;
    std::unordered_map<std::string, std::vector<std::string>> family_members;
    std::unordered_map<std::string, std::vector<std::string>> family_to_genera; // family → genera in that family
};

// qual_cache: gid → QualRecord from an existing QUAL section.
// When non-null, Phase 3 (SKCH scan) is skipped and sketch-derived scores are
// loaded from the cache instead.  TNF + taxonomy phases still run (needed for
// pass-B genus centroids).
PassAResult run_pass_a(ICheckReader& pack,
                       const std::vector<std::string>& accessions,
                       int threads,
                       int min_genus_size,
                       float leakage_threshold,
                       const std::unordered_map<uint64_t, QualRecord>* qual_cache = nullptr);

} // namespace genopack::check
