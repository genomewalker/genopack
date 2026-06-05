#pragma once
#include <filesystem>
#include <string>
#include <vector>

namespace genopack {

// `genopack report` — resolve every quality axis through a named profile and emit
// a unified per-genome table with a provenance column per axis (which tool/method
// produced each reported value). Fusion is never silent: the chosen column's
// identity is what the profile (built-in rule or stored policy) selected.
//   profile_name: built-in {intrinsic,external,calibrated,best} or a stored name.
//   output: TSV path; empty → stdout.
//   list:   ignore profile_name; enumerate built-in + stored profiles and the
//           available provenance columns per axis, then return.
int cmd_report(const std::filesystem::path& pack_path,
               const std::string& profile_name,
               const std::filesystem::path& output,
               bool list);

// `genopack profile add` — author and store a named SEC_PROF policy. Each selector
// is "axis=tool:method[@cal][/fallback_tool:fallback_method[@cal]]"; tool:method is
// resolved to the exact on-disk column identity present in the archive, so the
// stored profile pins content, not names. Written into each .gpk of the target.
int cmd_profile_add(const std::filesystem::path& pack_path,
                    const std::string& name,
                    const std::string& description,
                    const std::vector<std::string>& selectors);

// `genopack profile list` — show stored profiles (and their resolved provenance).
int cmd_profile_list(const std::filesystem::path& pack_path);

// `genopack qcontig` — dump the per-contig quality overlay (SEC_QCONTIG): one row
// per (genome, contig) with offset/length + GCOV T²/SPE percentiles + leakage and a
// `suspicious` flag, so a consumer can track which contigs drive a genome's
// contamination. accession_filter empty = all genomes. min_outlier > 0 restricts
// output to suspicious contigs — a contig is suspicious if its T²∪SPE percentile
// >= min_outlier (the cross-genus/foreign signal) or its leakage >= min_leakage.
int cmd_qcontig(const std::filesystem::path& pack_path,
                const std::string& accession_filter,
                const std::filesystem::path& output,
                float min_outlier, float min_leakage,
                float min_foreign, float min_lr);

} // namespace genopack
