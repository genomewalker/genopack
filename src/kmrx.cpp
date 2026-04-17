#include <genopack/kmrx.hpp>
#include <algorithm>
#include <cstring>
#include <stdexcept>

namespace genopack {

// ── KmrxWriter ────────────────────────────────────────────────────────────────

void KmrxWriter::add(GenomeId genome_id, const std::array<float, 136>& profile) {
    entries_.emplace_back(genome_id, profile);
}

SectionDesc KmrxWriter::finalize(AppendWriter& writer, uint64_t section_id) {
    // Sort entries by genome_id so the reader can binary-search
    std::sort(entries_.begin(), entries_.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    uint32_t n = static_cast<uint32_t>(entries_.size());

    // Layout:
    //   [KmrxHeader (32)]
    //   [uint64_t genome_ids[n]]
    //   [float profiles[n * 136]]
    uint64_t index_offset = sizeof(KmrxHeader);
    uint64_t data_offset  = index_offset + 8ULL * n;

    KmrxHeader hdr{};
    hdr.magic        = SEC_KMRX;
    hdr.n_genomes    = n;
    hdr.k            = 4;
    hdr.n_features   = 136;
    hdr.index_offset = index_offset;
    hdr.data_offset  = data_offset;

    uint64_t section_start = writer.current_offset();

    writer.append(&hdr, sizeof(hdr));

    // Write genome_id index
    for (const auto& [gid, _] : entries_) {
        uint64_t id = gid;
        writer.append(&id, sizeof(id));
    }

    // Write profiles
    for (const auto& [_, prof] : entries_) {
        writer.append(prof.data(), 136 * sizeof(float));
    }

    uint64_t section_end  = writer.current_offset();
    uint64_t payload_size = section_end - section_start;

    SectionDesc desc{};
    desc.type              = SEC_KMRX;
    desc.version           = 1;
    desc.flags             = 0;
    desc.section_id        = section_id;
    desc.file_offset       = section_start;
    desc.compressed_size   = payload_size;
    desc.uncompressed_size = payload_size;
    desc.item_count        = n;
    desc.aux0              = 4;    // k
    desc.aux1              = 136;  // n_features
    std::memset(desc.checksum, 0, sizeof(desc.checksum));
    return desc;
}

// ── KmrxReader ────────────────────────────────────────────────────────────────

void KmrxReader::open(const uint8_t* base, uint64_t offset, uint64_t size) {
    if (size < sizeof(KmrxHeader))
        throw std::runtime_error("KmrxReader: section too small");

    const uint8_t* payload = base + offset;
    const auto* hdr = reinterpret_cast<const KmrxHeader*>(payload);

    if (hdr->magic != SEC_KMRX)
        throw std::runtime_error("KmrxReader: bad magic");
    if (hdr->n_features != 136)
        throw std::runtime_error("KmrxReader: unexpected n_features");

    n_genomes_ = hdr->n_genomes;
    ids_       = reinterpret_cast<const uint64_t*>(payload + hdr->index_offset);
    profiles_  = reinterpret_cast<const float*>(payload + hdr->data_offset);
}

const float* KmrxReader::profile_for(GenomeId genome_id) const {
    if (!ids_ || n_genomes_ == 0) return nullptr;

    uint32_t lo = 0, hi = n_genomes_;
    while (lo < hi) {
        uint32_t mid = lo + (hi - lo) / 2;
        if (ids_[mid] == genome_id)
            return profiles_ + static_cast<size_t>(mid) * 136;
        if (ids_[mid] < genome_id) lo = mid + 1;
        else                        hi = mid;
    }
    return nullptr;
}

} // namespace genopack
