#pragma once
#include "types.hpp"
#include "toc.hpp"
#include "format.hpp"
#include <array>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace genopack {

// On-disk HNSW section header (64 bytes)
struct HnswSectionHeader {
    uint32_t magic;             // SEC_HNSW
    uint32_t version;           // 1
    uint32_t n_elements;        // number of vectors indexed
    uint32_t dim;               // vector dimension (136)
    uint32_t M;                 // hnswlib M param
    uint32_t ef_construction;   // hnswlib efConstruction param
    uint64_t index_offset;      // from section start -> hnswlib serialized blob
    uint64_t index_size;        // bytes
    uint64_t label_map_offset;  // from section start -> uint64_t[n_elements] genome_ids
    uint8_t  reserved[16];
};
static_assert(sizeof(HnswSectionHeader) == 64);

// Writes the HNSW index section.
class HnswSectionWriter {
public:
    void build(const std::vector<std::pair<GenomeId, std::array<float, 136>>>& profiles,
               uint32_t M = 16, uint32_t ef_construction = 200);
    SectionDesc finalize(AppendWriter& w, uint64_t section_id);

private:
    std::vector<uint8_t> index_blob_;
    std::vector<uint64_t> label_map_;   // label i -> genome_id
    uint32_t n_elements_ = 0;
    uint32_t M_ = 16, ef_ = 200;
};

// Loads and queries the HNSW index section.
class HnswSectionReader {
public:
    void open(const uint8_t* data, uint64_t offset, uint64_t size);
    bool is_open() const { return open_; }

    // Returns up to k (genome_id, cosine_similarity) pairs, sorted descending.
    std::vector<std::pair<GenomeId, float>>
    find_similar(const float* query_profile_136, size_t k) const;

private:
    bool open_ = false;
    struct Impl;
    std::shared_ptr<Impl> impl_;
};

} // namespace genopack
