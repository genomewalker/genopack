#pragma once
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace genopack {

// ── TaxnHeader: 32 bytes ──────────────────────────────────────────────────────
// Mirrors AccxHeader. Values are taxonomy strings instead of genome_ids.
// Layout: [TaxnHeader(32)] [acc_offsets:4*n] [tax_offsets:4*n]
//         [buckets:4*n_buckets] [acc_strings] [tax_strings]
struct TaxnHeader {
    uint32_t magic;              // SEC_TAXN
    uint32_t n_entries;
    uint32_t n_buckets;          // power of 2
    uint32_t _pad;
    uint64_t acc_strings_offset; // from section payload start
    uint64_t tax_strings_offset; // from section payload start
};
static_assert(sizeof(TaxnHeader) == 32);

// ── TaxonomyIndexWriter ───────────────────────────────────────────────────────

class TaxonomyIndexWriter {
public:
    struct Entry {
        std::string accession;
        std::string taxonomy;
    };

    void add(std::string accession, std::string taxonomy);

    SectionDesc finalize(AppendWriter& writer, uint64_t section_id);

    // Same as finalize() but skips the sort step (entries already sorted by accession).
    SectionDesc finalize_presorted(AppendWriter& writer, uint64_t section_id);

private:
    std::vector<Entry> entries_;
};

// ── TaxonomyIndexReader ───────────────────────────────────────────────────────

class TaxonomyIndexReader {
public:
    void open(const uint8_t* base, uint64_t offset, uint64_t size);

    std::optional<std::string_view> find(std::string_view accession) const;

    void scan(const std::function<void(std::string_view acc,
                                       std::string_view taxonomy)>& cb) const;

    size_t size() const;

private:
    const TaxnHeader* hdr_          = nullptr;
    const uint32_t*   acc_offsets_  = nullptr; // n_entries
    const uint32_t*   tax_offsets_  = nullptr; // n_entries
    const uint32_t*   buckets_      = nullptr; // n_buckets
    const char*       acc_strings_  = nullptr;
    const char*       tax_strings_  = nullptr;
};

} // namespace genopack
