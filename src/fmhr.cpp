#include <genopack/fmhr.hpp>
#include <cstring>

namespace genopack {

SectionDesc FmhrWriter::finalize(AppendWriter& w, uint64_t section_id) {
    const uint32_t n = static_cast<uint32_t>(entries_.size());
    const uint32_t min_buckets = (n == 0) ? 1
        : static_cast<uint32_t>(static_cast<double>(n) / 0.7 + 1.0);
    const uint32_t n_buckets = next_pow2(min_buckets);
    const uint32_t mask      = n_buckets - 1;

    std::vector<uint32_t> buckets(n_buckets, FMHR_EMPTY);
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t slot = static_cast<uint32_t>(entries_[i].genus_hash) & mask;
        while (buckets[slot] != FMHR_EMPTY) slot = (slot + 1) & mask;
        buckets[slot] = i;
    }

    w.align(8);
    const uint64_t section_start = w.current_offset();

    const uint64_t entries_off = sizeof(FmhrHeader);
    const uint64_t buckets_off = entries_off + static_cast<uint64_t>(n) * sizeof(FmhrEntry);
    // align hashes to 8 bytes
    const uint64_t hashes_off_unaligned = buckets_off + static_cast<uint64_t>(n_buckets) * 4;
    const uint64_t hashes_off = (hashes_off_unaligned + 7) & ~uint64_t{7};

    // Compute per-entry hash offsets (relative to section start)
    std::vector<FmhrEntry> disk_entries(n);
    uint64_t cur_hash_off = hashes_off;
    for (uint32_t i = 0; i < n; ++i) {
        disk_entries[i].genus_hash   = entries_[i].genus_hash;
        disk_entries[i].hashes_offset = cur_hash_off;  // relative to section start
        disk_entries[i].n_hashes     = static_cast<uint32_t>(entries_[i].hashes.size());
        disk_entries[i]._pad         = 0;
        cur_hash_off += static_cast<uint64_t>(entries_[i].hashes.size()) * sizeof(uint64_t);
    }

    FmhrHeader hdr{};
    hdr.magic          = SEC_FMHR;
    hdr.n_genera       = n;
    hdr.n_buckets      = n_buckets;
    hdr.k              = 21;
    hdr.c              = 125;
    hdr.entries_offset = entries_off;
    hdr.buckets_offset = buckets_off;
    hdr.hashes_offset  = hashes_off;

    w.append(&hdr, sizeof(hdr));
    if (n > 0)
        w.append(disk_entries.data(), static_cast<uint64_t>(n) * sizeof(FmhrEntry));
    w.append(buckets.data(), static_cast<uint64_t>(n_buckets) * sizeof(uint32_t));

    // Pad to 8-byte alignment before hash data
    if (const uint64_t pad = hashes_off - (buckets_off + n_buckets * 4); pad > 0) {
        uint64_t zero = 0;
        w.append(&zero, pad);
    }

    for (uint32_t i = 0; i < n; ++i) {
        if (!entries_[i].hashes.empty())
            w.append(entries_[i].hashes.data(),
                     static_cast<uint64_t>(entries_[i].hashes.size()) * sizeof(uint64_t));
    }

    const uint64_t section_end = w.current_offset();

    SectionDesc sd{};
    sd.type              = SEC_FMHR;
    sd.version           = 1;
    sd.flags             = 0;
    sd.section_id        = section_id;
    sd.file_offset       = section_start;
    sd.compressed_size   = section_end - section_start;
    sd.uncompressed_size = section_end - section_start;
    sd.item_count        = n;
    return sd;
}

} // namespace genopack
