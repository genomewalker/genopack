#pragma once
#include "types.hpp"
#include "format.hpp"
#include "mmap_file.hpp"
#include <cstdint>
#include <vector>

namespace genopack {

struct GidxHeader {        // 32 bytes
    uint32_t magic;        // SEC_GIDX
    uint32_t n_entries;
    uint32_t n_buckets;    // power of 2
    uint32_t _pad;
    uint64_t entries_offset; // from section start
    uint64_t buckets_offset; // from section start
};
static_assert(sizeof(GidxHeader) == 32);

struct GidxEntry {         // 32 bytes
    uint64_t genome_id;
    uint32_t shard_section_id;
    uint32_t dir_index;
    uint64_t catl_row_index;
    uint64_t _reserved;
};
static_assert(sizeof(GidxEntry) == 32);

class GidxWriter {
public:
    void add(GenomeId genome_id, uint32_t shard_section_id,
             uint32_t dir_index, uint64_t catl_row_index);
    SectionDesc finalize(AppendWriter& w, uint64_t section_id);

private:
    std::vector<GidxEntry> entries_;
};

class GidxReader {
public:
    void open(const uint8_t* data, uint64_t offset, uint64_t size);
    const GidxEntry* lookup(GenomeId genome_id) const;
    uint32_t n_entries() const { return header_ ? header_->n_entries : 0; }
    bool is_open() const { return data_ != nullptr; }

private:
    const uint8_t*    data_      = nullptr;
    const GidxHeader* header_    = nullptr;
    const GidxEntry*  entries_   = nullptr;
    const uint32_t*   buckets_   = nullptr;
    uint32_t          n_buckets_ = 0;
};

} // namespace genopack
