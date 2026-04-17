#pragma once
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <cstdint>
#include <vector>

namespace genopack {

// ── In-memory TOC representation ─────────────────────────────────────────────

struct Toc {
    TocHeader             header{};
    std::vector<SectionDesc> sections;

    std::vector<const SectionDesc*> find_by_type(uint32_t type) const;
    const SectionDesc*              find_by_id(uint64_t id) const;

    // Returns max section_id + 1, or 1 if sections is empty.
    uint64_t next_section_id() const;
};

// ── TocWriter ─────────────────────────────────────────────────────────────────

class TocWriter {
public:
    void add_section(SectionDesc desc);

    // Writes TocHeader + all SectionDescs to writer, then appends a TailLocator.
    // Returns the offset where the TocHeader was written.
    uint64_t finalize(AppendWriter& writer,
                      uint64_t generation,
                      uint64_t live_count,
                      uint64_t total_count,
                      uint64_t prev_toc_offset,
                      uint64_t catalog_root_id,
                      uint64_t accession_root_id,
                      uint64_t tombstone_root_id);

private:
    std::vector<SectionDesc> sections_;
};

// ── TocReader ─────────────────────────────────────────────────────────────────

class TocReader {
public:
    // Read from a fully memory-mapped file using the TailLocator at the end.
    static Toc read(const MmapFileReader& mmap);

    // Read from a known offset and size (offset points to TocHeader).
    static Toc read_at(const uint8_t* base, uint64_t offset, uint64_t size);
};

} // namespace genopack
