#pragma once
#include "format.hpp"
#include "toc.hpp"
#include <cstdint>
#include <vector>

namespace genopack {

// ── MetaBundleReader ──────────────────────────────────────────────────────────
// Reads the embedded MetaBundle (MetaHeader + SectionDesc[] + MetaTrailer) from
// the .gpk fd in a single pread. Used by ArchiveReader::open_gpk() so the
// complete section directory is loaded in one NFS round-trip, bypassing the
// mmap page-faults that occur near EOF of a multi-TB archive.
//
// MetaBundleWriter is part of TocWriter::finalize() — no standalone writer API.

class MetaBundleReader {
public:
    // Read bundle from fd at the offset/size stored in TailLocator.reserved.
    // Returns true if the data is present, checksum-valid, and magic-correct.
    bool open(int fd, uint64_t offset, uint64_t size);

    bool is_open() const { return open_; }

    const std::vector<SectionDesc>& sections()             const { return sections_; }
    uint64_t generation()                                  const { return hdr_.generation; }
    uint64_t live_genome_count()                           const { return hdr_.live_genome_count; }
    uint64_t total_genome_count()                          const { return hdr_.total_genome_count; }
    uint64_t catalog_root_section_id()                     const { return hdr_.catalog_root_section_id; }
    uint64_t accession_root_section_id()                   const { return hdr_.accession_root_section_id; }
    uint64_t tombstone_root_section_id()                   const { return hdr_.tombstone_root_section_id; }

    // Synthesise a Toc (callers that still expect one use this).
    Toc to_toc() const;

private:
    MetaHeader           hdr_{};
    std::vector<SectionDesc> sections_;
    bool                 open_ = false;
};

} // namespace genopack
