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
// ncbi_taxdump_dir: when non-empty, the merged TXDB is reconciled against the
//   NCBI taxdump there (nodes.dmp+names.dmp) — taxa whose scientific name is
//   unambiguous in NCBI get native NCBI taxids (NCBI precedence across sources),
//   so e.g. Bacteria/Archaea are taxid 2/2157 in both GTDB- and NCBI-derived packs.
void merge_archives(const std::vector<std::filesystem::path>& inputs,
                    const std::filesystem::path& output,
                    bool remap_genome_ids = true,
                    bool merge_cidx = true,
                    const std::string& spill_dir = {},
                    const std::string& ncbi_taxdump_dir = {});

} // namespace genopack
