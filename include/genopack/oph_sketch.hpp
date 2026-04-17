#pragma once
#include <cstdint>
#include <vector>

namespace genopack {

// OPH (One-Permutation Hashing) sketch result.
// m bins, each holds the minimum hash among k-mers assigned to that bin.
// Empty bins are densified deterministically.
struct OPHSketchResult {
    std::vector<uint32_t> signature;        // m bins (uint32_t; caller truncates to uint16_t for storage)
    std::vector<uint64_t> real_bins_bitmask; // ceil(m/64) words; bit t=1 iff bin t has a real k-mer
    uint64_t genome_length = 0;
    uint32_t n_real_bins   = 0;             // non-empty bins before densification
    uint32_t n_contigs     = 0;             // number of FASTA headers
};

// Dual-seed result: two OPH signatures from the same k-mer stream, one per seed.
// real_bins_bitmask is seed-independent (which bins were hit by any real k-mer),
// so it is stored once and shared by both signatures.
struct OPHDualSketchResult {
    std::vector<uint32_t> signature1;       // seed1
    std::vector<uint32_t> signature2;       // seed2
    std::vector<uint64_t> real_bins_bitmask;
    uint64_t genome_length = 0;
    uint32_t n_real_bins   = 0;
    uint32_t n_contigs     = 0;
};

// Sketcher configuration — mirrors geodesic's MinHasher::Config.
struct OPHSketchConfig {
    int      kmer_size   = 16;      // k=16 gives J~0.29 at 95% ANI
    int      sketch_size = 10000;   // number of OPH bins
    int      syncmer_s   = 0;       // 0 = disabled; >0 = open syncmer prefilter
    uint64_t seed        = 42;
};

// Compute OPH sketch from a pre-decompressed FASTA buffer.
// Thread-safe: no shared mutable state.
OPHSketchResult sketch_oph_from_buffer(const char* data, size_t len,
                                       const OPHSketchConfig& cfg);

// Compute BOTH seed=seed1 and seed=seed2 sketches in a SINGLE k-mer scan.
// The shared canonical k-mer stream is hashed twice per k-mer (once per
// seed), and each result feeds its own signature. real_bins_bitmask is
// identical to what a single-seed scan would produce.
OPHDualSketchResult sketch_oph_dual_from_buffer(const char* data, size_t len,
                                                int kmer_size, int sketch_size,
                                                int syncmer_s,
                                                uint64_t seed1, uint64_t seed2);

} // namespace genopack
