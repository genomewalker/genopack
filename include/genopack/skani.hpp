#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace genopack {

struct SkaniSketch {
    std::string accession;
    std::vector<uint64_t> hashes;  // sorted canonical k-mer hashes, density ~1/c
    uint64_t genome_length = 0;
};

struct SkaniResult {
    double ani  = 0.0;
    double af   = 0.0;
    double c_ab = 0.0;
    double c_ba = 0.0;
};

// Build a sketch using skani's compression model: select k-mer if canonical_hash % c == 0.
// Sketch size scales with genome_length / c, identical to skani -c parameter.
// fasta may be a raw sequence (no FASTA header) or a full multi-record FASTA.
SkaniSketch build_sketch(std::string_view accession, std::string_view fasta,
                          int k = 21, int c = 125);

SkaniResult compute_ani(const SkaniSketch& a, const SkaniSketch& b, int k = 21);

// Containment of query in ref: |query ∩ ref| / |query|. Both arrays must be sorted.
double fmh_containment(const std::vector<uint64_t>& query,
                        const std::vector<uint64_t>& ref) noexcept;
double fmh_containment(const std::vector<uint64_t>& query,
                        const uint64_t* ref_data, uint32_t ref_n) noexcept;

// One ntHash scan with k, applying nc c-filters simultaneously.
// Returns nc sorted/deduped hash vectors.
std::vector<std::vector<uint64_t>>
fmh_multi_c(std::string_view seq, int k, const int* c_vals, int nc);

} // namespace genopack
