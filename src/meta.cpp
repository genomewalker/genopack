#include <genopack/meta.hpp>
#include <genopack/checksum.hpp>
#include <cstring>
#include <stdexcept>
#include <unistd.h>
#include <vector>

namespace genopack {

bool MetaBundleReader::open(int fd, uint64_t offset, uint64_t size) {
    open_ = false;
    sections_.clear();

    const size_t min_size = sizeof(MetaHeader) + sizeof(MetaTrailer);
    if (size < min_size) return false;

    std::vector<uint8_t> buf(static_cast<size_t>(size));
    {
        uint8_t* p   = buf.data();
        off_t    pos = static_cast<off_t>(offset);
        size_t   rem = buf.size();
        while (rem > 0) {
            ssize_t n = ::pread(fd, p, rem, pos);
            if (n <= 0) return false;
            p   += static_cast<size_t>(n);
            pos += static_cast<off_t>(n);
            rem -= static_cast<size_t>(n);
        }
    }

    const MetaHeader* hdr = reinterpret_cast<const MetaHeader*>(buf.data());
    if (hdr->magic != GMET_MAGIC || hdr->version != 1) return false;

    const uint32_t n = hdr->n_sections;
    const size_t expected = sizeof(MetaHeader)
                          + static_cast<size_t>(n) * sizeof(SectionDesc)
                          + sizeof(MetaTrailer);
    if (static_cast<size_t>(size) != expected) return false;

    const MetaTrailer* trail = reinterpret_cast<const MetaTrailer*>(
        buf.data() + sizeof(MetaHeader) + static_cast<size_t>(n) * sizeof(SectionDesc));
    if (trail->magic      != GMET_TRAIL_MAGIC
        || trail->n_sections != n
        || trail->generation != hdr->generation)
        return false;

    if (!verify_checksum(buf.data(), buf.size(), offsetof(MetaHeader, checksum)))
        return false;

    hdr_ = *hdr;
    sections_.resize(n);
    const uint8_t* base = buf.data() + sizeof(MetaHeader);
    for (uint32_t i = 0; i < n; ++i)
        std::memcpy(&sections_[i], base + i * sizeof(SectionDesc), sizeof(SectionDesc));

    open_ = true;
    return true;
}

Toc MetaBundleReader::to_toc() const {
    Toc toc;
    toc.header.magic                     = TOCB_MAGIC;
    toc.header.version                   = FORMAT_MAJOR;
    toc.header.generation                = hdr_.generation;
    toc.header.live_genome_count         = hdr_.live_genome_count;
    toc.header.total_genome_count        = hdr_.total_genome_count;
    toc.header.section_count             = static_cast<uint64_t>(hdr_.n_sections);
    toc.header.catalog_root_section_id   = hdr_.catalog_root_section_id;
    toc.header.accession_root_section_id = hdr_.accession_root_section_id;
    toc.header.tombstone_root_section_id = hdr_.tombstone_root_section_id;
    toc.sections                         = sections_;
    return toc;
}

} // namespace genopack
