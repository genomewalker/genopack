#pragma once
#include <genopack/types.hpp>
#include <array>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace genopack {

// Decompress gzipped (or plain) FASTA file into a string.
std::string decompress_gz(const std::filesystem::path& path);

// Same but reads from an already-open file descriptor (closes it before returning).
// Use when the caller pre-opened with fadvise(WILLNEED) for NFS prefetch.
std::string decompress_gz_fd(int fd, const std::filesystem::path& path);

// MinHash minimum over k-mers; locality-sensitive fingerprint.
uint64_t genome_minhash(const std::string& fasta, int k = 21);

// Per-genome stats computed from a FASTA string.
struct FastaStats {
    uint64_t genome_length   = 0;
    uint32_t n_contigs       = 0;
    uint16_t gc_pct_x100     = 0;   // 0-10000, e.g. 5234 = 52.34%
    uint64_t oph_fingerprint = 0;   // MinHash minimum (k=21)
    // Canonical k=4 tetranucleotide frequencies, L2-normalised.
    // 136 = number of distinct canonical 4-mers over {A,C,G,T}.
    std::array<float, 136> kmer4_profile = {};
};

FastaStats compute_fasta_stats(const std::string& fasta, int k = 21);

// Days elapsed since 2024-01-01 (the archive epoch).
uint32_t days_since_epoch();

// Parse a TSV file into BuildRecord list.
// Recognised columns: accession/acc, file_path/path/fasta_path/fasta/file,
//   completeness, contamination, genome_length, n_contigs.
// Unknown columns go to BuildRecord::extra_fields.
std::vector<BuildRecord> parse_tsv_records(const std::filesystem::path& tsv_path);

} // namespace genopack
