#include <genopack/qcontig.hpp>
#include <algorithm>

namespace genopack {

SectionDesc QcontigWriter::finalize(AppendWriter& w, uint64_t section_id) {
    std::sort(recs_.begin(), recs_.end(),
              [](const QcontigRecord& a, const QcontigRecord& b) {
                  if (a.genome_id != b.genome_id) return a.genome_id < b.genome_id;
                  return a.contig_index < b.contig_index;
              });

    w.align(8);
    const uint64_t section_start = w.current_offset();

    QcontigHeader hdr{};
    hdr.magic          = SEC_QCONTIG;
    hdr.format_version = QCONTIG_FORMAT_VERSION;
    hdr.flags          = 0;
    hdr.n_rows         = recs_.size();
    hdr.records_off    = sizeof(QcontigHeader);
    w.append(&hdr, sizeof(hdr));
    if (!recs_.empty())
        w.append(recs_.data(), recs_.size() * sizeof(QcontigRecord));

    const uint64_t section_end = w.current_offset();
    SectionDesc sd{};
    sd.type              = SEC_QCONTIG;
    sd.version           = QCONTIG_FORMAT_VERSION;
    sd.flags             = 0;
    sd.section_id        = section_id;
    sd.file_offset       = section_start;
    sd.compressed_size   = section_end - section_start;
    sd.uncompressed_size = section_end - section_start;
    sd.item_count        = recs_.size();
    sd.header_size       = sizeof(QcontigHeader);
    return sd;
}

} // namespace genopack
