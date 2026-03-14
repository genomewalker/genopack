#pragma once
#include "types.hpp"
#include "catalog.hpp"
#include "shard.hpp"
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace genopack {

// ── Archive directory layout ──────────────────────────────────────────────────
//
// dataset.gpk/
//   MANIFEST.bin      — version, shard list, generation, checksums
//   taxon.idx         — taxon name → taxon_id + shard_ids + catalog row range
//   catalog.gpkc      — columnar metadata (all genomes, SoA, sorted by taxon_id)
//   tombstones.bits   — deleted genome_id bitset (indexed by genome_id)
//   shards/
//     genus_0001_00001.gpks
//     genus_0001_00002.gpks   ← split when shard exceeds 512 MB
//     genus_0042_00001.gpks
//   sketches/           ← optional, separate from hot read path
//     sketch_00001.gpkh

struct ManifestHeader {
    uint32_t magic;           // GPKM_MAGIC
    uint16_t version;
    uint16_t _pad0;
    uint64_t generation;      // monotonically increasing, incremented on each write
    uint64_t created_at;      // unix timestamp
    uint32_t n_shards;
    uint32_t n_taxons;
    uint64_t n_genomes;       // total (including deleted)
    uint64_t n_genomes_live;  // excluding deleted
    uint64_t catalog_size;
    uint64_t taxon_idx_size;
    uint8_t  catalog_checksum[16];
    uint8_t  taxon_idx_checksum[16];
    uint8_t  _pad1[32];
};
static_assert(sizeof(ManifestHeader) == 128, "ManifestHeader layout changed");

// ── TaxonIndex ────────────────────────────────────────────────────────────────
// Binary-searchable index: taxon_name → {taxon_id, shard_ids[], catalog_rows[]}

struct TaxonIndexEntry {
    TaxonId  taxon_id;
    TaxonId  parent_id;
    uint8_t  rank;           // TaxonRank
    uint8_t  _pad[3];
    uint32_t n_shards;
    uint32_t shard_list_offset;  // offset into shard_list[] array
    uint32_t catalog_row_lo;     // first row in catalog for this taxon
    uint32_t catalog_row_hi;     // one past last row
    uint32_t n_genomes;
    uint32_t _pad2;
    uint32_t _pad3;
};
static_assert(sizeof(TaxonIndexEntry) == 40, "TaxonIndexEntry layout changed");

// ── ArchiveReader ─────────────────────────────────────────────────────────────

class ArchiveReader {
public:
    ArchiveReader();
    ~ArchiveReader();
    ArchiveReader(const ArchiveReader&) = delete;
    ArchiveReader& operator=(const ArchiveReader&) = delete;

    void open(const std::filesystem::path& archive_dir);
    void close();
    bool is_open() const;

    // Metadata without sequence decompression
    std::optional<GenomeMeta>  genome_meta(GenomeId id)         const;
    std::optional<TaxonInfo>   taxon_info(TaxonId id)           const;
    std::optional<TaxonInfo>   taxon_by_name(std::string_view name) const;
    std::optional<TaxonId>     resolve_taxon(std::string_view lineage_or_name) const;

    // Count genomes matching query (no decompression)
    size_t count(const ExtractQuery& q) const;

    // Metadata-only filter result
    std::vector<GenomeMeta>  filter_meta(const ExtractQuery& q) const;

    // Full extraction with FASTA decompression
    std::vector<ExtractedGenome> extract(const ExtractQuery& q) const;

    // Fetch one genome by id
    std::optional<ExtractedGenome> fetch_genome(GenomeId id) const;

    // Fetch all genomes for a taxon (species most common)
    std::vector<ExtractedGenome> fetch_taxon(TaxonId taxon_id) const;
    std::vector<ExtractedGenome> fetch_taxon(std::string_view name) const;

    // Summary statistics
    struct TaxonStats {
        TaxonId     taxon_id;
        std::string name;
        size_t      n_genomes;
        size_t      n_deleted;
        double      mean_completeness;
        double      mean_contamination;
        uint64_t    total_bp;
        double      mean_gc;
    };
    std::vector<TaxonStats> taxon_summary(TaxonRank rank) const;

    // Archive-wide stats
    struct ArchiveStats {
        uint64_t generation;
        size_t   n_shards;
        size_t   n_taxons;
        size_t   n_genomes_total;
        size_t   n_genomes_live;
        uint64_t total_raw_bp;
        uint64_t total_compressed_bytes;
        double   compression_ratio;
    };
    ArchiveStats archive_stats() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// ── ArchiveBuilder ────────────────────────────────────────────────────────────
// Builds a new archive from scratch. Not thread-safe (single writer).

struct ArchiveBuilderConfig {
    ShardWriterConfig shard_cfg;
    size_t io_threads = 4;      // parallel FASTA decompression
    bool   verbose    = false;
    size_t batch_size = 1000;   // genomes buffered before shard flush
};

class ArchiveBuilder {
public:
    using Config = ArchiveBuilderConfig;

    explicit ArchiveBuilder(const std::filesystem::path& archive_dir,
                            Config cfg = Config{});
    ~ArchiveBuilder();

    // Add all genomes from a TSV (accession, taxonomy, file_path [, completeness, contamination])
    void add_from_tsv(const std::filesystem::path& tsv_path);

    // Add individual record
    void add(const BuildRecord& rec);

    // Write all shards, catalog, taxon index, manifest
    void finalize();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// ── ArchiveAppender ───────────────────────────────────────────────────────────
// Append genomes to an existing archive (new shard generation).

class ArchiveAppender {
public:
    explicit ArchiveAppender(const std::filesystem::path& archive_dir);
    ~ArchiveAppender();

    void add_from_tsv(const std::filesystem::path& tsv_path);
    void add(const BuildRecord& rec);

    // Tombstone genomes (mark deleted; does not remove from shards)
    void remove(GenomeId id);
    void remove_by_accession(std::string_view accession);

    // Write new shard(s), update catalog, increment generation
    void commit();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace genopack
