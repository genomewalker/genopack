#pragma once
#include <genopack/types.hpp>
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace genopack {

// ── AccxHeader: 32 bytes ──────────────────────────────────────────────────────
// 4+4+4+4+8+8 = 32
struct AccxHeader {
    uint32_t magic;           // SEC_ACCX = 0x58434341
    uint32_t n_entries;
    uint32_t n_buckets;       // power of 2
    uint32_t _pad;
    uint64_t strings_offset;  // from section payload start
    uint64_t buckets_offset;  // from section payload start
};
static_assert(sizeof(AccxHeader) == 32);

// ── AccessionIndexWriter ──────────────────────────────────────────────────────

class AccessionIndexWriter {
public:
    void add(std::string accession, GenomeId genome_id);

    // Sort, pack strings, build hash table, write to writer.
    // Returns a filled SectionDesc (file_offset etc. set).
    SectionDesc finalize(AppendWriter& writer, uint64_t section_id);

    // Same as finalize() but skips the sort step (entries already sorted by accession).
    SectionDesc finalize_presorted(AppendWriter& writer, uint64_t section_id);

private:
    struct Entry {
        std::string accession;
        GenomeId    genome_id;
    };
    std::vector<Entry> entries_;
};

// ── AccessionIndexReader ──────────────────────────────────────────────────────

class AccessionIndexReader {
public:
    // base: start of the mmap'd file buffer; offset/size: section payload position.
    void open(const uint8_t* base, uint64_t offset, uint64_t size);

    std::optional<GenomeId> find(std::string_view accession) const;

    // Iterate all entries: callback(accession_sv, genome_id).
    void scan(const std::function<void(std::string_view, GenomeId)>& cb) const;

    size_t size() const;

private:
    const AccxHeader* hdr_     = nullptr;
    const uint64_t*   ids_     = nullptr;  // n_entries
    const uint32_t*   offsets_ = nullptr;  // n_entries (string offsets from strings_area)
    const uint32_t*   buckets_ = nullptr;  // n_buckets
    const char*       strings_ = nullptr;  // packed null-terminated strings
};

} // namespace genopack
