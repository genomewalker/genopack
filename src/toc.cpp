#include <genopack/toc.hpp>
#include <genopack/checksum.hpp>
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <xxhash.h>

namespace genopack {

// ── Toc ───────────────────────────────────────────────────────────────────────

std::vector<const SectionDesc*> Toc::find_by_type(uint32_t type) const {
    std::vector<const SectionDesc*> result;
    for (const auto& s : sections)
        if (s.type == type) result.push_back(&s);
    return result;
}

const SectionDesc* Toc::find_by_id(uint64_t id) const {
    for (const auto& s : sections)
        if (s.section_id == id) return &s;
    return nullptr;
}

uint64_t Toc::next_section_id() const {
    if (sections.empty()) return 1;
    uint64_t max_id = 0;
    for (const auto& s : sections)
        if (s.section_id > max_id) max_id = s.section_id;
    return max_id + 1;
}

// ── TocWriter ─────────────────────────────────────────────────────────────────

void TocWriter::add_section(SectionDesc desc) {
    sections_.push_back(desc);
}

uint64_t TocWriter::finalize(AppendWriter& writer,
                              uint64_t generation,
                              uint64_t live_count,
                              uint64_t total_count,
                              uint64_t prev_toc_offset,
                              uint64_t catalog_root_id,
                              uint64_t accession_root_id,
                              uint64_t tombstone_root_id,
                              uint64_t gpk_uuid_lo,
                              uint64_t gpk_uuid_hi)
{
    uint64_t n = sections_.size();

    // ── Embedded MetaBundle ───────────────────────────────────────────────────
    // Written BEFORE the TOC so its offset is known when we write the TailLocator.
    // One pread of this region gives the complete section directory — the reader
    // never needs to mmap-fault near EOF of a multi-TB archive.
    writer.align(8);
    const uint64_t meta_bundle_offset = writer.current_offset();

    {
        const uint32_t nsec = static_cast<uint32_t>(sections_.size());
        const size_t total_bytes = sizeof(MetaHeader)
                                 + static_cast<size_t>(nsec) * sizeof(SectionDesc)
                                 + sizeof(MetaTrailer);
        std::vector<uint8_t> mbuf(total_bytes, 0);

        MetaHeader* mhdr = reinterpret_cast<MetaHeader*>(mbuf.data());
        mhdr->magic                     = GMET_MAGIC;
        mhdr->version                   = 1;
        mhdr->n_sections                = nsec;
        mhdr->generation                = generation;
        mhdr->gpk_uuid_lo               = gpk_uuid_lo;
        mhdr->gpk_uuid_hi               = gpk_uuid_hi;
        mhdr->live_genome_count         = live_count;
        mhdr->total_genome_count        = total_count;
        mhdr->catalog_root_section_id   = catalog_root_id;
        mhdr->accession_root_section_id = accession_root_id;
        mhdr->tombstone_root_section_id = tombstone_root_id;

        uint8_t* sec_base = mbuf.data() + sizeof(MetaHeader);
        for (uint32_t i = 0; i < nsec; ++i)
            std::memcpy(sec_base + i * sizeof(SectionDesc),
                        &sections_[i], sizeof(SectionDesc));

        MetaTrailer* trail = reinterpret_cast<MetaTrailer*>(
            mbuf.data() + sizeof(MetaHeader) + static_cast<size_t>(nsec) * sizeof(SectionDesc));
        trail->magic      = GMET_TRAIL_MAGIC;
        trail->n_sections = nsec;
        trail->generation = generation;

        compute_checksum(mbuf.data(), total_bytes, offsetof(MetaHeader, checksum));
        writer.append(mbuf.data(), static_cast<uint64_t>(total_bytes));
    }

    const uint64_t meta_bundle_size = writer.current_offset() - meta_bundle_offset;

    TocHeader hdr{};
    hdr.magic                    = TOCB_MAGIC;
    hdr.version                  = FORMAT_MAJOR;
    hdr.flags                    = 0;
    hdr.generation               = generation;
    hdr.prev_toc_offset          = prev_toc_offset;
    hdr.section_count            = n;
    hdr.live_genome_count        = live_count;
    hdr.total_genome_count       = total_count;
    hdr.catalog_root_section_id  = catalog_root_id;
    hdr.accession_root_section_id = accession_root_id;
    hdr.tombstone_root_section_id = tombstone_root_id;
    std::memset(hdr.checksum, 0, sizeof(hdr.checksum));
    std::memset(hdr.reserved,  0, sizeof(hdr.reserved));

    // Compute XXH128 over (header || all SectionDescs) with checksum field zeroed.
    // toc_buf is kept in scope so we can hash the same bytes for toc_checksum below.
    size_t toc_block_size = sizeof(TocHeader) + sections_.size() * sizeof(SectionDesc);
    std::vector<uint8_t> toc_buf(toc_block_size);
    std::memcpy(toc_buf.data(), &hdr, sizeof(hdr));
    for (size_t i = 0; i < sections_.size(); ++i)
        std::memcpy(toc_buf.data() + sizeof(TocHeader) + i * sizeof(SectionDesc),
                    &sections_[i], sizeof(SectionDesc));
    compute_checksum(toc_buf.data(), toc_block_size, offsetof(TocHeader, checksum));
    std::memcpy(&hdr, toc_buf.data(), sizeof(hdr));

    uint64_t toc_start = writer.current_offset();

    writer.append(&hdr, sizeof(hdr));
    for (const auto& s : sections_)
        writer.append(&s, sizeof(s));

    uint64_t toc_end  = writer.current_offset();
    uint64_t toc_size = toc_end - toc_start;

    // Flush all pending writes to the server before writing the TailLocator.
    // On NFS, write-back caching means ENOSPC may be reported on a later write
    // for data buffered earlier. fsync() forces all dirty pages to the server,
    // ensuring the file size and content are stable before we commit the footer.
    writer.flush();

    TailLocator tail{};
    tail.magic       = GPKT_MAGIC;
    tail.version     = FORMAT_MAJOR;
    tail.flags       = 0;
    tail.toc_offset  = toc_start;
    tail.toc_size    = toc_size;
    tail.generation  = generation;
    tail.set_meta_bundle(meta_bundle_offset, meta_bundle_size);
    // toc_checksum covers the same TOC bytes we just wrote (toc_buf has the
    // final content including the computed TocHeader.checksum).
    {
        XXH128_hash_t h = XXH3_128bits(toc_buf.data(), toc_block_size);
        XXH128_canonical_t canon;
        XXH128_canonicalFromHash(&canon, h);
        std::memcpy(tail.toc_checksum, canon.digest, 16);
    }

    writer.append(&tail, sizeof(tail));

    return toc_start;
}

// ── TocReader ─────────────────────────────────────────────────────────────────

Toc TocReader::read(const MmapFileReader& mmap) {
    uint64_t sz = mmap.size();
    if (sz < sizeof(TailLocator))
        throw std::runtime_error("TocReader: file too small to contain TailLocator");

    uint64_t tail_offset = sz - sizeof(TailLocator);
    const TailLocator* tail = mmap.ptr_at<TailLocator>(tail_offset);

    if (tail->magic != GPKT_MAGIC)
        throw std::runtime_error("TocReader: TailLocator magic mismatch");

    return read_at(mmap.data(), tail->toc_offset, tail->toc_size);
}

Toc TocReader::read_at(const uint8_t* base, uint64_t offset, uint64_t size) {
    if (size < sizeof(TocHeader))
        throw std::runtime_error("TocReader: TOC block too small");

    const TocHeader* hdr = reinterpret_cast<const TocHeader*>(base + offset);

    if (hdr->magic != TOCB_MAGIC)
        throw std::runtime_error("TocReader: TocHeader magic mismatch");

    uint64_t n = hdr->section_count;
    uint64_t required = sizeof(TocHeader) + n * sizeof(SectionDesc);
    if (required > size)
        throw std::runtime_error("TocReader: TOC block truncated");

    Toc toc;
    toc.header = *hdr;

    const SectionDesc* descs =
        reinterpret_cast<const SectionDesc*>(base + offset + sizeof(TocHeader));
    toc.sections.assign(descs, descs + n);

    return toc;
}

} // namespace genopack
