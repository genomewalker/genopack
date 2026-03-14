#pragma once
#include "types.hpp"
#include <filesystem>
#include <functional>
#include <memory>
#include <vector>

namespace genopack {

// ── CatalogWriter ─────────────────────────────────────────────────────────────
// Writes catalog.gpkc: memory-mappable columnar metadata store.
// Rows sorted by taxon_id for binary-search taxonomy lookup.
// Column layout: SoA (struct of arrays) for fast predicate filtering.

class CatalogWriter {
public:
    explicit CatalogWriter(const std::filesystem::path& path);
    ~CatalogWriter();

    void add(const GenomeMeta& meta);
    void finalize();  // sort by taxon_id, write footer, close

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// ── CatalogReader ─────────────────────────────────────────────────────────────
// Memory-mapped read-only view of catalog.gpkc.
// Supports: lookup by genome_id, filter by taxon_id range, predicate scan.

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

    // All genomes for a taxon_id (species or genus range)
    std::vector<const GenomeMeta*> for_taxon(TaxonId taxon_id) const;

    // Range of genomes for a taxon_id subtree [lo, hi)
    std::vector<const GenomeMeta*> for_taxon_range(TaxonId lo, TaxonId hi) const;

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
    Stats compute_stats(TaxonId taxon_id = INVALID_TAXON_ID) const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// ── On-disk catalog file layout ──────────────────────────────────────────────
//
// [GPKC magic 4B] [version 2B] [n_rows 8B] [n_row_groups 4B]
// [catalog_data_offset 8B] [footer_offset 8B] [schema_hash 16B]
// --- column arrays (SoA, sorted by taxon_id) ---
// Column 0: genome_id[]       uint64  × n_rows
// Column 1: taxon_id[]        uint32  × n_rows
// Column 2: shard_id[]        uint32  × n_rows
// Column 3: blob_offset[]     uint64  × n_rows
// Column 4: blob_len_cmp[]    uint32  × n_rows
// Column 5: blob_len_raw[]    uint32  × n_rows
// Column 6: genome_length[]   uint64  × n_rows
// Column 7: n_contigs[]       uint32  × n_rows
// Column 8: gc_pct_x100[]     uint16  × n_rows
// Column 9: completeness_x10[] uint16 × n_rows
// Column 10: contam_x10[]     uint16  × n_rows
// Column 11: oph_fingerprint[] uint64 × n_rows
// Column 12: date_added[]     uint32  × n_rows
// Column 13: flags[]          uint32  × n_rows
// --- row group stats (min/max per 4096 rows for pushdown) ---
// RowGroupStats × n_row_groups
// --- footer: column offsets + checksums ---
//
// All columns are 8-byte aligned. The genome_id column is also used as the
// lookup key for a binary-searchable auxiliary index within the file.

struct RowGroupStats {
    uint32_t first_row;
    uint32_t last_row;
    TaxonId  taxon_id_min;
    TaxonId  taxon_id_max;
    uint64_t genome_length_min;   // placed before uint16_t to avoid implicit padding
    uint64_t genome_length_max;
    uint16_t completeness_min;
    uint16_t completeness_max;
    uint32_t flags_any;           // OR of all flags in group (nonzero → has deleted rows)
};
static_assert(sizeof(RowGroupStats) == 40, "RowGroupStats layout changed");

} // namespace genopack
