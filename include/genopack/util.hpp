#pragma once
#include <genopack/types.hpp>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <span>
#include <string>
#include <string_view>
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
    float    gc_skew_closure = NAN; // 1 - |final_cum| / max_abs; NAN if genome_length < 10000
    // Canonical k=4 tetranucleotide frequencies, L2-normalised.
    // 136 = number of distinct canonical 4-mers over {A,C,G,T}.
    std::array<float, 136> kmer4_profile = {};
};

FastaStats compute_fasta_stats(const std::string& fasta, int k = 21);

// Per-contig TNF contamination check. Returns fraction of total bp from contigs whose
// L2-normalised TNF distance to tnf_centroid[136] is below threshold (default 3.0).
// Contigs < 5000 bp are always counted clean (matching check pass-B behaviour).
// Returns 1.0f if tnf_centroid is nullptr.
float compute_completeness_post_decontam(
    const std::string& fasta,
    const float*       tnf_centroid,
    float              threshold = 3.0f);

// ── Multi-field contamination signal ──────────────────────────────────────────

// Fused single-pass computation of FASTA-level quality signals.
// Currently surfaces the near-clade TNF-GMM minority mass (span overload only).
struct AllSignals {
    float    contamination_tnf_minor = NAN; // near-clade TNF-GMM minority mass (BIC-free, multi-contig guarded); span overload only
};

// Pre-parsed contig variant: avoids raw FASTA re-scan and duplicate TNF computation.
// Pass one entry per contig; tnf136 is the L2-normalised TNF from compute_tnf (null if
// contig was too short to score). Contigs shorter than min_contig_len are skipped.
struct ContigAccum {
    std::string_view seq;      // raw sequence (needed for the mix-window walk)
    const float*     tnf136;  // pre-computed L2-norm TNF[136], or nullptr
    uint32_t         len;
};
AllSignals compute_all_signals(std::span<const ContigAccum> contigs,
                                int  min_contig_len  = 5000,
                                int  mix_window      = 16384,
                                int  min_mix_windows = 5,
                                bool skip_mixture    = false);

// Thread-aggregated sub-phase nanosecond totals from compute_all_signals (ContigAccum overload).
struct CasTimings { int64_t scan_ns, coh_ns, chargaff_ns, spectral_ns, mix_ns; };
CasTimings get_cas_timings();

// L2-normalised k=4 TNF profile per contig, for contigs >= min_bp.
// Uses the same canonical k=4 index as compute_fasta_stats().
struct ContigProfile {
    float    p[136];       // L2-normalised k=4 TNF profile
    float    freq_norm;    // ‖raw_frequency‖₂ before L2-norm (recover f_raw = p * freq_norm)
    uint32_t bp;
    float    rho[16];      // Karlin dinucleotide relative abundance: rho[i*4+j]=P(ij)/(P(i)*P(j))
};
std::vector<ContigProfile> compute_long_contig_profiles(const std::string& fasta,
                                                         uint32_t min_bp = 20000);

// Days elapsed since 2024-01-01 (the archive epoch).
uint32_t days_since_epoch();

// Parse a TSV file into BuildRecord list.
// Recognised columns: accession/acc, file_path/path/fasta_path/fasta/file,
//   completeness, contamination, genome_length, n_contigs.
// Unknown columns go to BuildRecord::extra_fields.
std::vector<BuildRecord> parse_tsv_records(const std::filesystem::path& tsv_path);

} // namespace genopack
