#include <genopack/tombstone.hpp>
#include <algorithm>
#include <stdexcept>

namespace genopack {

// ── TombstoneWriter ───────────────────────────────────────────────────────────

void TombstoneWriter::add(GenomeId id) {
    ids_.push_back(id);
}

SectionDesc TombstoneWriter::finalize(AppendWriter& writer, uint64_t section_id) {
    std::sort(ids_.begin(), ids_.end());
    ids_.erase(std::unique(ids_.begin(), ids_.end()), ids_.end());

    uint32_t n = static_cast<uint32_t>(ids_.size());

    TombHeader hdr{};
    hdr.magic     = SEC_TOMB;
    hdr.n_entries = n;
    hdr.reserved  = 0;

    uint64_t section_start = writer.current_offset();

    writer.append(&hdr, sizeof(hdr));
    if (n > 0)
        writer.append(ids_.data(), 8ULL * n);

    uint64_t section_end  = writer.current_offset();
    uint64_t payload_size = section_end - section_start;

    SectionDesc desc{};
    desc.type              = SEC_TOMB;
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

// ── TombstoneReader ───────────────────────────────────────────────────────────

void TombstoneReader::open(const uint8_t* base, uint64_t offset, uint64_t size) {
    if (size < sizeof(TombHeader))
        throw std::runtime_error("TombstoneReader: section too small");

    const uint8_t* payload = base + offset;
    const auto* hdr = reinterpret_cast<const TombHeader*>(payload);

    if (hdr->magic != SEC_TOMB)
        throw std::runtime_error("TombstoneReader: bad magic");

    uint32_t n = hdr->n_entries;
    if (sizeof(TombHeader) + 8ULL * n > size)
        throw std::runtime_error("TombstoneReader: section truncated");

    count_ = n;
    ids_   = (n > 0)
        ? reinterpret_cast<const GenomeId*>(payload + sizeof(TombHeader))
        : nullptr;
}

bool TombstoneReader::is_deleted(GenomeId id) const {
    if (!ids_ || count_ == 0) return false;
    return std::binary_search(ids_, ids_ + count_, id);
}

} // namespace genopack
