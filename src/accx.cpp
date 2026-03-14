#include <genopack/accx.hpp>
#include <algorithm>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace genopack {

// FNV-1a 32-bit hash
static uint32_t fnv1a(std::string_view s) {
    uint32_t h = 2166136261u;
    for (unsigned char c : s)
        h = (h ^ c) * 16777619u;
    return h;
}

// ── AccessionIndexWriter ──────────────────────────────────────────────────────

void AccessionIndexWriter::add(std::string accession, GenomeId genome_id) {
    entries_.push_back({std::move(accession), genome_id});
}

SectionDesc AccessionIndexWriter::finalize(AppendWriter& writer, uint64_t section_id) {
    // Sort by accession for binary-search capability and deterministic output.
    std::sort(entries_.begin(), entries_.end(),
              [](const Entry& a, const Entry& b) { return a.accession < b.accession; });

    uint32_t n = static_cast<uint32_t>(entries_.size());

    // Pack null-terminated strings and build offset table.
    std::string strings_area;
    std::vector<uint32_t> str_offsets;
    str_offsets.reserve(n);
    for (const auto& e : entries_) {
        str_offsets.push_back(static_cast<uint32_t>(strings_area.size()));
        strings_area.append(e.accession);
        strings_area.push_back('\0');
    }

    // Build power-of-2 bucket table (linear probing, load factor ≤ 0.5).
    uint32_t n_buckets = 1;
    while (n_buckets < 2 * n) n_buckets <<= 1;
    std::vector<uint32_t> buckets(n_buckets, UINT32_MAX);
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t h = fnv1a(entries_[i].accession) & (n_buckets - 1);
        while (buckets[h] != UINT32_MAX)
            h = (h + 1) & (n_buckets - 1);
        buckets[h] = i;
    }

    // Compute section payload layout.
    // Payload: [AccxHeader(32)] [genome_ids: 8*n] [str_offsets: 4*n] [buckets: 4*n_buckets] [strings]
    uint64_t ids_offset      = sizeof(AccxHeader);
    uint64_t strofft_offset  = ids_offset + 8ULL * n;
    uint64_t buckets_offset  = strofft_offset + 4ULL * n;
    uint64_t strings_offset  = buckets_offset + 4ULL * n_buckets;

    AccxHeader hdr{};
    hdr.magic          = SEC_ACCX;
    hdr.n_entries      = n;
    hdr.n_buckets      = n_buckets;
    hdr._pad           = 0;
    hdr.strings_offset = strings_offset;
    hdr.buckets_offset = buckets_offset;

    uint64_t section_start = writer.current_offset();

    writer.append(&hdr, sizeof(hdr));

    // genome_id array
    for (const auto& e : entries_)
        writer.append(&e.genome_id, sizeof(e.genome_id));

    // string_offset array
    writer.append(str_offsets.data(), 4ULL * n);

    // FNV bucket array
    writer.append(buckets.data(), 4ULL * n_buckets);

    // packed strings
    writer.append(strings_area.data(), strings_area.size());

    uint64_t section_end = writer.current_offset();
    uint64_t payload_size = section_end - section_start;

    SectionDesc desc{};
    desc.type              = SEC_ACCX;
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

// ── AccessionIndexReader ──────────────────────────────────────────────────────

void AccessionIndexReader::open(const uint8_t* base, uint64_t offset, uint64_t size) {
    if (size < sizeof(AccxHeader))
        throw std::runtime_error("AccessionIndexReader: section too small");

    const uint8_t* payload = base + offset;
    hdr_ = reinterpret_cast<const AccxHeader*>(payload);

    if (hdr_->magic != SEC_ACCX)
        throw std::runtime_error("AccessionIndexReader: bad magic");

    uint32_t n         = hdr_->n_entries;
    uint32_t n_buckets = hdr_->n_buckets;

    uint64_t ids_offset     = sizeof(AccxHeader);
    uint64_t strofft_offset = ids_offset + 8ULL * n;
    uint64_t needed_min     = strofft_offset + 4ULL * n + 4ULL * n_buckets;
    if (needed_min > size)
        throw std::runtime_error("AccessionIndexReader: section truncated");

    ids_     = reinterpret_cast<const uint64_t*>(payload + ids_offset);
    offsets_ = reinterpret_cast<const uint32_t*>(payload + strofft_offset);
    buckets_ = reinterpret_cast<const uint32_t*>(payload + hdr_->buckets_offset);
    strings_ = reinterpret_cast<const char*>(payload + hdr_->strings_offset);
}

std::optional<GenomeId> AccessionIndexReader::find(std::string_view accession) const {
    if (!hdr_ || hdr_->n_entries == 0) return std::nullopt;

    uint32_t n_buckets = hdr_->n_buckets;
    uint32_t h = fnv1a(accession) & (n_buckets - 1);

    for (uint32_t probe = 0; probe < n_buckets; ++probe) {
        uint32_t idx = buckets_[(h + probe) & (n_buckets - 1)];
        if (idx == UINT32_MAX) return std::nullopt;
        const char* s = strings_ + offsets_[idx];
        if (accession == s)
            return ids_[idx];
    }
    return std::nullopt;
}

void AccessionIndexReader::scan(
    const std::function<void(std::string_view, GenomeId)>& cb) const
{
    if (!hdr_) return;
    uint32_t n = hdr_->n_entries;
    for (uint32_t i = 0; i < n; ++i) {
        const char* s = strings_ + offsets_[i];
        cb(std::string_view(s), ids_[i]);
    }
}

size_t AccessionIndexReader::size() const {
    return hdr_ ? hdr_->n_entries : 0;
}

} // namespace genopack
