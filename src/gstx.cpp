#include <genopack/gstx.hpp>
#include <algorithm>
#include <cstring>
#include <numeric>
#include <stdexcept>

namespace genopack {

// FNV-1a 64-bit hash; remap 0 → 1 so 0 stays the "empty" sentinel.
uint64_t GstxWriter::hash_genus(std::string_view s) noexcept {
    uint64_t h = 0xcbf29ce484222325ULL;
    for (unsigned char c : s) {
        h ^= c;
        h *= 0x100000001b3ULL;
    }
    return h == 0 ? 1 : h;
}

uint64_t GstxReader::hash_genus(std::string_view s) noexcept {
    uint64_t h = 0xcbf29ce484222325ULL;
    for (unsigned char c : s) {
        h ^= c;
        h *= 0x100000001b3ULL;
    }
    return h == 0 ? 1 : h;
}

uint32_t GstxWriter::next_pow2(uint32_t v) noexcept {
    if (v == 0) return 1;
    --v;
    v |= v >> 1; v |= v >> 2; v |= v >> 4; v |= v >> 8; v |= v >> 16;
    return v + 1;
}

// ── GstxWriter ────────────────────────────────────────────────────────────────

void GstxWriter::add_genus(std::string_view genus,
                            uint32_t         n_members,
                            uint32_t         n_k,
                            const std::vector<std::vector<uint16_t>>& consensus,
                            const float*     p90,
                            const float*     tnf_mu,
                            const uint32_t*  kmer_sizes,
                            float            nrb_p90)
{
    GstxEntry e{};
    e.genus_hash  = hash_genus(genus);
    e.n_members   = n_members;
    e.n_k_stored  = static_cast<uint8_t>(std::min(n_k, GSTX_MAX_K));
    e.flags       = tnf_mu ? 0x1u : 0x0u;

    for (uint32_t ki = 0; ki < e.n_k_stored; ++ki) {
        e.p90_containment[ki] = p90[ki];
        if (ki < consensus.size() && consensus[ki].size() == GSTX_BINS)
            std::memcpy(e.consensus[ki], consensus[ki].data(),
                        GSTX_BINS * sizeof(uint16_t));
    }
    if (tnf_mu)
        std::memcpy(e.tnf_mu, tnf_mu, 136 * sizeof(float));
    e.nrb_p90 = nrb_p90;

    // Update header k-metadata (all genera in a single build have the same k-set)
    if (n_k_hdr_ == 0 && n_k > 0) {
        n_k_hdr_ = e.n_k_stored;
        for (uint32_t ki = 0; ki < e.n_k_stored && ki < 4; ++ki)
            kmer_sizes_hdr_[ki] = kmer_sizes[ki];
    }

    entry_order_.push_back(static_cast<uint32_t>(entries_.size()));
    entries_.push_back(e);
}

SectionDesc GstxWriter::finalize(AppendWriter& w, uint64_t section_id) {
    // Sort entries by arrival index so parallel builds produce the same byte
    // sequence as sequential builds (entry_order_ is 0,1,2,... when sequential).
    if (!entry_order_.empty()) {
        std::vector<size_t> idx(entries_.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(),
            [&](size_t a, size_t b){ return entry_order_[a] < entry_order_[b]; });
        std::vector<GstxEntry> sorted;
        sorted.reserve(entries_.size());
        for (size_t i : idx) sorted.push_back(entries_[i]);
        entries_ = std::move(sorted);
        entry_order_.clear();
    }
    const uint32_t n = static_cast<uint32_t>(entries_.size());
    const uint32_t min_buckets = (n == 0) ? 1
        : static_cast<uint32_t>(static_cast<double>(n) / 0.7 + 1.0);
    const uint32_t n_buckets = next_pow2(min_buckets);
    const uint32_t mask      = n_buckets - 1;

    static constexpr uint32_t EMPTY = UINT32_MAX;
    std::vector<uint32_t> buckets(n_buckets, EMPTY);
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t slot = static_cast<uint32_t>(entries_[i].genus_hash) & mask;
        while (buckets[slot] != EMPTY) slot = (slot + 1) & mask;
        buckets[slot] = i;
    }

    w.align(8);
    const uint64_t section_start = w.current_offset();

    const uint64_t entries_off = sizeof(GstxHeader);
    const uint64_t buckets_off = entries_off + static_cast<uint64_t>(n) * sizeof(GstxEntry);

    GstxHeader hdr{};
    hdr.magic          = SEC_GSTX;
    hdr.n_genera       = n;
    hdr.n_buckets      = n_buckets;
    hdr.n_k            = n_k_hdr_;
    hdr.sketch_size    = GSTX_BINS;
    for (int i = 0; i < 4; ++i) hdr.kmer_sizes[i] = kmer_sizes_hdr_[i];
    hdr.entries_offset = entries_off;
    hdr.buckets_offset = buckets_off;
    hdr.entry_stride   = sizeof(GstxEntry);

    w.append(&hdr, sizeof(hdr));
    if (n > 0)
        w.append(entries_.data(), static_cast<uint64_t>(n) * sizeof(GstxEntry));
    w.append(buckets.data(), static_cast<uint64_t>(n_buckets) * sizeof(uint32_t));

    const uint64_t section_end = w.current_offset();

    SectionDesc sd{};
    sd.type              = SEC_GSTX;
    sd.version           = 1;
    sd.flags             = 0;
    sd.section_id        = section_id;
    sd.file_offset       = section_start;
    sd.compressed_size   = section_end - section_start;
    sd.uncompressed_size = section_end - section_start;
    sd.item_count        = n;
    return sd;
}

// ── GstxReader ────────────────────────────────────────────────────────────────

void GstxReader::open(const uint8_t* data, uint64_t offset, uint64_t size) {
    if (size < sizeof(GstxHeader))
        throw std::runtime_error("GSTX section too small");

    data_   = data + offset;
    header_ = reinterpret_cast<const GstxHeader*>(data_);

    if (header_->magic != SEC_GSTX)
        throw std::runtime_error("GSTX: bad magic");
    if (header_->entry_stride != sizeof(GstxEntry))
        throw std::runtime_error("GSTX: entry_stride mismatch — rebuild required");
    if (header_->sketch_size != GSTX_BINS)
        throw std::runtime_error("GSTX: sketch_size mismatch");

    const uint64_t entries_end = header_->entries_offset
        + static_cast<uint64_t>(header_->n_genera) * sizeof(GstxEntry);
    const uint64_t buckets_end = header_->buckets_offset
        + static_cast<uint64_t>(header_->n_buckets) * sizeof(uint32_t);

    if (entries_end > size || buckets_end > size)
        throw std::runtime_error("GSTX section truncated");

    entries_   = reinterpret_cast<const GstxEntry*>(data_ + header_->entries_offset);
    buckets_   = reinterpret_cast<const uint32_t*>(data_ + header_->buckets_offset);
    n_buckets_ = header_->n_buckets;
}

const GstxEntry* GstxReader::lookup(std::string_view genus) const {
    if (!data_ || n_buckets_ == 0) return nullptr;
    static constexpr uint32_t EMPTY = UINT32_MAX;
    const uint64_t h    = hash_genus(genus);
    const uint32_t mask = n_buckets_ - 1;
    uint32_t slot = static_cast<uint32_t>(h) & mask;
    for (;;) {
        const uint32_t idx = buckets_[slot];
        if (idx == EMPTY) return nullptr;
        if (entries_[idx].genus_hash == h) return &entries_[idx];
        slot = (slot + 1) & mask;
    }
}

} // namespace genopack
