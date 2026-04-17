#pragma once
#include "types.hpp"
#include "format.hpp"
#include "mmap_file.hpp"
#include <filesystem>
#include <functional>
#include <memory>
#include <vector>

namespace genopack {

// ── CatalogWriter ─────────────────────────────────────────────────────────────
// Writes catalog.gpkc: memory-mappable columnar metadata store.
// Rows sorted by oph_fingerprint for locality-sensitive similarity lookup.
// Column layout: SoA (struct of arrays) for fast predicate filtering.

class CatalogWriter {
public:
    explicit CatalogWriter(const std::filesystem::path& path);
    ~CatalogWriter();

    void add(const GenomeMeta& meta);
    void finalize();  // sort by oph_fingerprint, write footer, close

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// ── CatalogReader ─────────────────────────────────────────────────────────────
// Memory-mapped read-only view of catalog.gpkc.
// Supports: lookup by genome_id, filter by oph_fingerprint range, predicate scan.

class CatalogReader {
public:
    CatalogReader();
    ~CatalogReader();
    CatalogReader(const CatalogReader&) = delete;
    CatalogReader& operator=(const CatalogReader&) = delete;

    void open(const std::filesystem::path& path);
    void close();

    size_t n_rows() const;
    size_t n_genomes() const;   // excludes deleted

    // Exact lookup by genome_id. Returns nullptr if not found.
    const GenomeMeta* find_genome(GenomeId id) const;

    // All genomes in a MinHash fingerprint range [lo, hi]
    std::vector<const GenomeMeta*> for_oph_range(uint64_t lo, uint64_t hi) const;

    // General predicate scan
    std::vector<const GenomeMeta*> filter(const ExtractQuery& q) const;
    void scan(std::function<bool(const GenomeMeta&)> predicate) const;

    // Statistics
    struct Stats {
        size_t   n_total;
        size_t   n_deleted;
        uint64_t total_genome_length;
        double   mean_completeness;
        double   mean_contamination;
        double   mean_gc;
    };
    Stats compute_stats() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// ── On-disk catalog file layout ──────────────────────────────────────────────
//
// [GPKC magic 4B] [version 2B] [n_rows 8B] [n_row_groups 4B]
// [catalog_data_offset 8B] [footer_offset 8B] [schema_hash 16B]
// --- column arrays (SoA, sorted by oph_fingerprint) ---
// Column 0: genome_id[]       uint64  x n_rows
// Column 1: _reserved0[]      uint32  x n_rows  (was: taxon_id; kept for binary compat)
// Column 2: shard_id[]        uint32  x n_rows
// Column 3: blob_offset[]     uint64  x n_rows
// Column 4: blob_len_cmp[]    uint32  x n_rows
// Column 5: blob_len_raw[]    uint32  x n_rows
// Column 6: genome_length[]   uint64  x n_rows
// Column 7: n_contigs[]       uint32  x n_rows
// Column 8: gc_pct_x100[]     uint16  x n_rows
// Column 9: completeness_x10[] uint16 x n_rows
// Column 10: contam_x10[]     uint16  x n_rows
// Column 11: oph_fingerprint[] uint64 x n_rows
// Column 12: date_added[]     uint32  x n_rows
// Column 13: flags[]          uint32  x n_rows
// --- row group stats (min/max per 4096 rows for pushdown) ---
// RowGroupStats x n_row_groups
// --- footer: column offsets + checksums ---
//
// All columns are 8-byte aligned.

struct RowGroupStats {
    uint32_t first_row;
    uint32_t last_row;
    uint64_t oph_min;
    uint64_t oph_max;
    uint64_t genome_length_min;
    uint64_t genome_length_max;
    uint16_t completeness_min;
    uint16_t completeness_max;
    uint32_t flags_any;
};
static_assert(sizeof(RowGroupStats) == 48, "RowGroupStats layout changed");

// ── V2 row-group stats (used by CatalogSectionWriter/Reader) -----------------

struct RowGroupStatsV2 {
    uint32_t first_row;
    uint32_t last_row;           // inclusive
    uint32_t live_count;         // non-deleted rows in this group
    uint32_t _pad;
    uint64_t oph_min;
    uint64_t oph_max;
    uint64_t genome_length_min;
    uint64_t genome_length_max;
    uint16_t completeness_min;
    uint16_t completeness_max;
    uint16_t contamination_min;
    uint16_t contamination_max;
    uint32_t flags_any;
    uint8_t  reserved[12];
};
static_assert(sizeof(RowGroupStatsV2) == 72, "RowGroupStatsV2 layout changed");

// ── CATL section header (32 bytes) -------------------------------------------

struct CatlHeader {
    uint32_t magic;          // SEC_CATL = 0x4C544143
    uint32_t n_rows;
    uint32_t n_groups;
    uint32_t row_group_size;
    uint64_t stats_offset;   // from section start to RowGroupStatsV2 array
    uint64_t rows_offset;    // from section start to GenomeMeta array
};
static_assert(sizeof(CatlHeader) == 32, "CatlHeader layout changed");

// ── CatalogSectionWriter (v2) -------------------------------------------------
// Accumulates GenomeMeta rows, then writes a CATL section into an AppendWriter.
// Caller must add rows in oph_fingerprint order (already sorted).

class CatalogSectionWriter {
public:
    static constexpr size_t CATL_DEFAULT_ROW_GROUP_SIZE = 32768;

    explicit CatalogSectionWriter(size_t row_group_size = CATL_DEFAULT_ROW_GROUP_SIZE);

    // Add a genome metadata row. Must be called in oph_fingerprint order.
    void add(const GenomeMeta& m);

    // Finalize and write CATL section to AppendWriter.
    // Returns SectionDesc describing where it was written.
    SectionDesc finalize(AppendWriter& writer, uint64_t section_id);

    size_t n_rows() const { return rows_.size(); }

private:
    std::vector<GenomeMeta>      rows_;
    std::vector<RowGroupStatsV2> groups_;
    size_t row_group_size_;

    void flush_group(size_t first, size_t last);
};

// ── CatalogSectionReader (v2) -------------------------------------------------
// Zero-copy read of a CATL section from a memory-mapped region.

class CatalogSectionReader {
public:
    // Open from a memory-mapped region.
    void open(const uint8_t* base, uint64_t section_offset, uint64_t section_size);

    // Scan all rows matching a query (zero-copy, returns pointers into mmap).
    std::vector<const GenomeMeta*> filter(const ExtractQuery& q) const;

    // Find by genome_id (linear scan).
    const GenomeMeta* find_genome(GenomeId id) const;

    // Direct zero-copy row access.
    const GenomeMeta* row_at(uint32_t index) const;

    // Full scan; return false from cb to stop.
    void scan(std::function<bool(const GenomeMeta&)> cb) const;

    size_t n_rows() const;
    size_t n_groups() const;

private:
    const CatlHeader*      header_ = nullptr;
    const RowGroupStatsV2* groups_ = nullptr;
    const GenomeMeta*      rows_   = nullptr;

    bool row_matches(const GenomeMeta& m, const ExtractQuery& q) const;
    bool group_may_match(const RowGroupStatsV2& g, const ExtractQuery& q) const;
};

// ── MergedCatalogReader -------------------------------------------------------
// Merges multiple CATL fragments (oldest first). Newer fragment wins on
// duplicate genome_id.

class MergedCatalogReader {
public:
    // Add a fragment in generation order (oldest first).
    void add_fragment(const uint8_t* base, uint64_t offset, uint64_t size);

    // Same API as CatalogSectionReader.
    std::vector<const GenomeMeta*> filter(const ExtractQuery& q) const;
    const GenomeMeta* find_genome(GenomeId id) const;
    void scan(std::function<bool(const GenomeMeta&)> cb) const;
    size_t n_rows() const;

private:
    std::vector<CatalogSectionReader> fragments_;
};

} // namespace genopack
