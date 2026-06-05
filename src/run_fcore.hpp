#pragma once
#include <filesystem>

namespace genopack {

// `genopack fcore` — build SEC_FCORE per-family prevalence cores: the genus-core
// fallback for novel/sparse genera. Reuses the existing CORE params (k, min_seg_aa,
// theta, frac_max) so family-core coverage is directly comparable to genus-core
// coverage at check time. A second FASTA pass over the archive's shards extracts
// per-genome aamer sets grouped by family. theta_override > 0 overrides the
// prevalence threshold; families with fewer than min_members genomes are skipped
// (a single-member "core" is just that genome and carries no prevalence signal).
// members_list: optional file of reference accessions (one per line). When given,
// only those genomes contribute to the family cores — the rest of the archive's
// live genomes (e.g. query/test genomes co-resident in the same .gpk) are excluded
// from core construction. Empty = every live family-labeled genome contributes.
int cmd_fcore(const std::filesystem::path& pack_path, int threads,
              float theta_override, unsigned min_members,
              const std::filesystem::path& members_list);

} // namespace genopack
