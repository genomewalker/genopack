#pragma once
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/types.hpp>
#include <array>
#include <cstdint>
#include <utility>
#include <vector>

namespace genopack {

// ── KmrxHeader: 32 bytes ─────────────────────────────────────────────────────
// Layout (from section start):
//   [KmrxHeader (32 bytes)]
//   [uint64_t genome_ids[n_genomes]]   ← index_offset from section start
//   [float profiles[n_genomes * 136]]  ← data_offset from section start
//
// genome_ids[] is sorted ascending. profiles[] is parallel to genome_ids[].
// Each float[136] is an L2-normalised canonical k=4 tetranucleotide frequency vector.
// Stored uncompressed (float data compresses poorly).

struct KmrxHeader {
    uint32_t magic;           // SEC_KMRX
    uint32_t n_genomes;
    uint32_t k;               // 4
    uint32_t n_features;      // 136
    uint64_t index_offset;    // from section start → uint64_t[n_genomes] sorted genome_ids
    uint64_t data_offset;     // from section start → float[n_genomes * 136] profiles
};
static_assert(sizeof(KmrxHeader) == 32);

// ── KmrxWriter ───────────────────────────────────────────────────────────────

class KmrxWriter {
public:
    void add(GenomeId genome_id, const std::array<float, 136>& profile);

    // Write section and return its descriptor.
    SectionDesc finalize(AppendWriter& writer, uint64_t section_id);

private:
    std::vector<std::pair<GenomeId, std::array<float, 136>>> entries_;
};

// ── KmrxReader ───────────────────────────────────────────────────────────────

class KmrxReader {
public:
    void open(const uint8_t* base, uint64_t offset, uint64_t size);

    // Returns pointer to float[136] L2-normalised profile, nullptr if not found.
    // Binary search on sorted genome_ids[].
    const float* profile_for(GenomeId genome_id) const;

    size_t n_genomes() const { return n_genomes_; }

private:
    const uint64_t* ids_      = nullptr;
    const float*    profiles_ = nullptr;
    uint32_t        n_genomes_ = 0;
};

// ── cosine_similarity ─────────────────────────────────────────────────────────
// Dot product of two L2-normalised vectors of length 136.
inline float cosine_similarity(const float* a, const float* b) {
    float dot = 0.0f;
    for (int i = 0; i < 136; ++i) dot += a[i] * b[i];
    return dot;
}

} // namespace genopack
