#pragma once
#include <filesystem>

namespace genopack {

// Ingest external quality measurements (CheckM2 / anvi'o) into an archive's SEC_XQAL
// columnar section, keyed by genome_id. Additive: merges with any existing XQAL.
// At least one of checkm2_tsv / anvio_tsv must be non-empty.
int cmd_ingest(const std::filesystem::path& pack_path,
               const std::filesystem::path& checkm2_tsv = {},
               const std::filesystem::path& anvio_tsv = {});

} // namespace genopack
