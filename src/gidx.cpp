#include <genopack/gidx.hpp>
#include <cstring>
#include <stdexcept>

namespace genopack {

// FNV-1a hash for uint64_t genome_id
static uint64_t hash_genome_id(uint64_t id) {
    uint64_t h = 0xcbf29ce484222325ULL;
    for (int i = 0; i < 8; ++i) {
        h ^= (id >> (i * 8)) & 0xFF;
        h *= 0x100000001b3ULL;
    }
    return h;
}

static uint32_t next_power_of_2(uint32_t v) {
    if (v == 0) return 1;
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    return v + 1;
}

// ── GidxWriter ──────────────────────────────────────────────────────────────

void GidxWriter::add(GenomeId genome_id, uint32_t shard_section_id,
                     uint32_t dir_index, uint64_t catl_row_index)
{
    GidxEntry e{};
    e.genome_id        = genome_id;
    e.shard_section_id = shard_section_id;
    e.dir_index        = dir_index;
    e.catl_row_index   = catl_row_index;
    e._reserved        = 0;
    entries_.push_back(e);
}

SectionDesc GidxWriter::finalize(AppendWriter& w, uint64_t section_id) {
    uint32_t n = static_cast<uint32_t>(entries_.size());

    // n_buckets >= n / 0.7, rounded up to power of 2
    uint32_t min_buckets = (n == 0) ? 1 : static_cast<uint32_t>(
        static_cast<double>(n) / 0.7 + 1.0);
    uint32_t n_buckets = next_power_of_2(min_buckets);

    // Build bucket array (open addressing, linear probing)
    static constexpr uint32_t EMPTY = UINT32_MAX;
    std::vector<uint32_t> buckets(n_buckets, EMPTY);

    uint32_t mask = n_buckets - 1;
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t slot = static_cast<uint32_t>(hash_genome_id(entries_[i].genome_id)) & mask;
        while (buckets[slot] != EMPTY)
            slot = (slot + 1) & mask;
        buckets[slot] = i;
    }

    // Write section
    w.align(8);
    uint64_t section_start = w.current_offset();

    GidxHeader hdr{};
    hdr.magic          = SEC_GIDX;
    hdr.n_entries      = n;
    hdr.n_buckets      = n_buckets;
    hdr._pad           = 0;
    hdr.entries_offset = sizeof(GidxHeader);
    hdr.buckets_offset = sizeof(GidxHeader) + static_cast<uint64_t>(n) * sizeof(GidxEntry);

    w.append(&hdr, sizeof(hdr));
    if (n > 0)
        w.append(entries_.data(), static_cast<uint64_t>(n) * sizeof(GidxEntry));
    w.append(buckets.data(), static_cast<uint64_t>(n_buckets) * sizeof(uint32_t));

    uint64_t section_end = w.current_offset();

    SectionDesc sd{};
    sd.type              = SEC_GIDX;
    sd.version           = 1;
    sd.flags             = 0;
    sd.section_id        = section_id;
    sd.file_offset       = section_start;
    sd.compressed_size   = section_end - section_start;
    sd.uncompressed_size = section_end - section_start;
    sd.item_count        = n;
    sd.aux0              = 0;
    sd.aux1              = 0;
    std::memset(sd.checksum, 0, sizeof(sd.checksum));
    return sd;
}

// ── GidxReader ──────────────────────────────────────────────────────────────

void GidxReader::open(const uint8_t* data, uint64_t offset, uint64_t size) {
    if (size < sizeof(GidxHeader))
        throw std::runtime_error("GIDX section too small");

    data_ = data + offset;
    header_ = reinterpret_cast<const GidxHeader*>(data_);

    if (header_->magic != SEC_GIDX)
        throw std::runtime_error("GIDX: bad magic");

    uint64_t entries_end = header_->entries_offset +
        static_cast<uint64_t>(header_->n_entries) * sizeof(GidxEntry);
    uint64_t buckets_end = header_->buckets_offset +
        static_cast<uint64_t>(header_->n_buckets) * sizeof(uint32_t);

    if (entries_end > size || buckets_end > size)
        throw std::runtime_error("GIDX section truncated");

    entries_   = reinterpret_cast<const GidxEntry*>(data_ + header_->entries_offset);
    buckets_   = reinterpret_cast<const uint32_t*>(data_ + header_->buckets_offset);
    n_buckets_ = header_->n_buckets;
}

const GidxEntry* GidxReader::lookup(GenomeId genome_id) const {
    if (!data_ || n_buckets_ == 0) return nullptr;

    static constexpr uint32_t EMPTY = UINT32_MAX;
    uint32_t mask = n_buckets_ - 1;
    uint32_t slot = static_cast<uint32_t>(hash_genome_id(genome_id)) & mask;

    for (;;) {
        uint32_t idx = buckets_[slot];
        if (idx == EMPTY) return nullptr;
        if (entries_[idx].genome_id == genome_id)
            return &entries_[idx];
        slot = (slot + 1) & mask;
    }
}

} // namespace genopack
