#include <genopack/taxn.hpp>
#include <algorithm>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace genopack {

// FNV-1a 32-bit hash (same as accx.cpp)
static uint32_t fnv1a(std::string_view s) {
    uint32_t h = 2166136261u;
    for (unsigned char c : s)
        h = (h ^ c) * 16777619u;
    return h;
}

// ── MSD radix sort ───────────────────────────────────────────────────────────

// Byte-level MSD radix sort over a range of entries.
// sort_by_taxonomy=true  → sort key is entry.taxonomy
// sort_by_taxonomy=false → sort key is entry.accession
// Falls back to std::sort for small ranges (≤ 32 entries).
static void msd_radix_sort(std::vector<TaxonomyIndexWriter::Entry>& entries,
                           size_t lo, size_t hi,
                           int depth, bool sort_by_taxonomy)
{
    if (hi - lo <= 32) {
        if (sort_by_taxonomy) {
            std::sort(entries.begin() + static_cast<ptrdiff_t>(lo),
                      entries.begin() + static_cast<ptrdiff_t>(hi),
                      [](const TaxonomyIndexWriter::Entry& a,
                         const TaxonomyIndexWriter::Entry& b) {
                          return a.taxonomy < b.taxonomy;
                      });
        } else {
            std::sort(entries.begin() + static_cast<ptrdiff_t>(lo),
                      entries.begin() + static_cast<ptrdiff_t>(hi),
                      [](const TaxonomyIndexWriter::Entry& a,
                         const TaxonomyIndexWriter::Entry& b) {
                          return a.accession < b.accession;
                      });
        }
        return;
    }

    // Count occurrences of each byte at position depth; 256 + 1 (sentinel for exhausted).
    // Bucket 0 = string exhausted at this depth (shorter strings come first).
    // Buckets 1..256 = byte value 0..255.
    size_t cnt[257] = {};
    for (size_t i = lo; i < hi; ++i) {
        const std::string& key = sort_by_taxonomy ? entries[i].taxonomy : entries[i].accession;
        if (depth < static_cast<int>(key.size()))
            ++cnt[1 + static_cast<unsigned char>(key[depth])];
        else
            ++cnt[0];
    }

    // Compute prefix sums → bucket start positions.
    size_t start[258];
    start[0] = lo;
    for (int b = 0; b <= 256; ++b)
        start[b + 1] = start[b] + cnt[b];

    // Scatter into a temp buffer, then move back. O(N) space, unconditionally correct.
    std::vector<TaxonomyIndexWriter::Entry> tmp(hi - lo);
    size_t cur[257];
    for (int b = 0; b <= 256; ++b) cur[b] = start[b] - lo;  // relative to tmp base
    for (size_t i = lo; i < hi; ++i) {
        const std::string& key = sort_by_taxonomy ? entries[i].taxonomy : entries[i].accession;
        int b = (depth < static_cast<int>(key.size()))
                    ? 1 + static_cast<unsigned char>(key[depth])
                    : 0;
        tmp[cur[b]++] = std::move(entries[i]);
    }
    for (size_t i = lo; i < hi; ++i)
        entries[i] = std::move(tmp[i - lo]);

    // Recurse into each non-trivial bucket (skip bucket 0 — already terminal at depth).
    for (int b = 1; b <= 256; ++b) {
        size_t blo = start[b], bhi = start[b + 1];
        if (bhi - blo > 1)
            msd_radix_sort(entries, blo, bhi, depth + 1, sort_by_taxonomy);
    }
}

// ── TaxonomyIndexWriter ───────────────────────────────────────────────────────

void TaxonomyIndexWriter::add(std::string accession, std::string taxonomy) {
    entries_.push_back({std::move(accession), std::move(taxonomy)});
}

// Internal helper that packs and writes entries in their current order.
static SectionDesc write_taxn_entries(
    const std::vector<TaxonomyIndexWriter::Entry>& entries,
    AppendWriter& writer, uint64_t section_id)
{
    uint32_t n = static_cast<uint32_t>(entries.size());

    std::string acc_area;
    std::vector<uint32_t> acc_offsets;
    acc_offsets.reserve(n);
    for (const auto& e : entries) {
        acc_offsets.push_back(static_cast<uint32_t>(acc_area.size()));
        acc_area.append(e.accession);
        acc_area.push_back('\0');
    }

    std::string tax_area;
    std::vector<uint32_t> tax_offsets;
    tax_offsets.reserve(n);
    for (const auto& e : entries) {
        tax_offsets.push_back(static_cast<uint32_t>(tax_area.size()));
        tax_area.append(e.taxonomy);
        tax_area.push_back('\0');
    }

    uint32_t n_buckets = 1;
    while (n_buckets < 2 * n) n_buckets <<= 1;
    std::vector<uint32_t> buckets(n_buckets, UINT32_MAX);
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t h = fnv1a(entries[i].accession) & (n_buckets - 1);
        while (buckets[h] != UINT32_MAX)
            h = (h + 1) & (n_buckets - 1);
        buckets[h] = i;
    }

    uint64_t acc_off_start = sizeof(TaxnHeader);
    uint64_t tax_off_start = acc_off_start + 4ULL * n;
    uint64_t buckets_start = tax_off_start + 4ULL * n;
    uint64_t acc_str_start = buckets_start + 4ULL * n_buckets;
    uint64_t tax_str_start = acc_str_start + acc_area.size();

    TaxnHeader hdr{};
    hdr.magic              = SEC_TAXN;
    hdr.n_entries          = n;
    hdr.n_buckets          = n_buckets;
    hdr._pad               = 0;
    hdr.acc_strings_offset = acc_str_start;
    hdr.tax_strings_offset = tax_str_start;

    uint64_t section_start = writer.current_offset();

    writer.append(&hdr,               sizeof(hdr));
    writer.append(acc_offsets.data(), 4ULL * n);
    writer.append(tax_offsets.data(), 4ULL * n);
    writer.append(buckets.data(),     4ULL * n_buckets);
    writer.append(acc_area.data(),    acc_area.size());
    writer.append(tax_area.data(),    tax_area.size());

    uint64_t section_end  = writer.current_offset();
    uint64_t payload_size = section_end - section_start;

    SectionDesc desc{};
    desc.type              = SEC_TAXN;
    desc.version           = 1;
    desc.flags             = 0;
    desc.section_id        = section_id;
    desc.file_offset       = section_start;
    desc.compressed_size   = payload_size;
    desc.uncompressed_size = payload_size;
    desc.item_count        = n;
    desc.aux0              = 0;
    desc.aux1              = 0;
    return desc;
}

SectionDesc TaxonomyIndexWriter::finalize(AppendWriter& writer, uint64_t section_id) {
    // MSD radix sort by taxonomy VALUE — exploits massive shared prefixes (80-130 chars)
    // for O(S) behaviour vs O(S log N) for comparison sort.
    // Hash-table lookup uses accession key regardless of storage order, so this is safe.
    msd_radix_sort(entries_, 0, entries_.size(), 0, /*sort_by_taxonomy=*/true);
    return write_taxn_entries(entries_, writer, section_id);
}

SectionDesc TaxonomyIndexWriter::finalize_presorted(AppendWriter& writer, uint64_t section_id) {
    // Entries are already sorted by accession — skip the sort step entirely.
    return write_taxn_entries(entries_, writer, section_id);
}

// ── TaxonomyIndexReader ───────────────────────────────────────────────────────

void TaxonomyIndexReader::open(const uint8_t* base, uint64_t offset, uint64_t size) {
    if (size < sizeof(TaxnHeader))
        throw std::runtime_error("TaxonomyIndexReader: section too small");

    const uint8_t* payload = base + offset;
    hdr_ = reinterpret_cast<const TaxnHeader*>(payload);

    if (hdr_->magic != SEC_TAXN)
        throw std::runtime_error("TaxonomyIndexReader: bad magic");

    uint32_t n         = hdr_->n_entries;
    uint32_t n_buckets = hdr_->n_buckets;

    uint64_t acc_off_start = sizeof(TaxnHeader);
    uint64_t tax_off_start = acc_off_start + 4ULL * n;
    uint64_t needed_min    = tax_off_start + 4ULL * n + 4ULL * n_buckets;
    if (needed_min > size)
        throw std::runtime_error("TaxonomyIndexReader: section truncated");

    uint64_t buckets_start = tax_off_start + 4ULL * n;
    acc_offsets_ = reinterpret_cast<const uint32_t*>(payload + acc_off_start);
    tax_offsets_ = reinterpret_cast<const uint32_t*>(payload + tax_off_start);
    buckets_     = reinterpret_cast<const uint32_t*>(payload + buckets_start);
    acc_strings_ = reinterpret_cast<const char*>(payload + hdr_->acc_strings_offset);
    tax_strings_ = reinterpret_cast<const char*>(payload + hdr_->tax_strings_offset);
}

std::optional<std::string_view> TaxonomyIndexReader::find(std::string_view accession) const {
    if (!hdr_ || hdr_->n_entries == 0) return std::nullopt;

    uint32_t n_buckets = hdr_->n_buckets;
    uint32_t h = fnv1a(accession) & (n_buckets - 1);

    for (uint32_t probe = 0; probe < n_buckets; ++probe) {
        uint32_t idx = buckets_[(h + probe) & (n_buckets - 1)];
        if (idx == UINT32_MAX) return std::nullopt;
        const char* s = acc_strings_ + acc_offsets_[idx];
        if (accession == s)
            return std::string_view(tax_strings_ + tax_offsets_[idx]);
    }
    return std::nullopt;
}

void TaxonomyIndexReader::scan(
    const std::function<void(std::string_view, std::string_view)>& cb) const
{
    if (!hdr_) return;
    uint32_t n = hdr_->n_entries;
    for (uint32_t i = 0; i < n; ++i) {
        std::string_view acc(acc_strings_ + acc_offsets_[i]);
        std::string_view tax(tax_strings_ + tax_offsets_[i]);
        cb(acc, tax);
    }
}

size_t TaxonomyIndexReader::size() const {
    return hdr_ ? hdr_->n_entries : 0;
}

} // namespace genopack
