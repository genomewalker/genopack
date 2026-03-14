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
//   catalog.gpkc      — columnar metadata (all genomes, SoA, sorted by oph_fingerprint)
//   meta.tsv          — sidecar: accession <-> genome_id mapping (optional)
//   tombstones.bits   — deleted genome_id bitset (indexed by genome_id)
//   shards/
//     shard_00001.gpks
//     shard_00002.gpks   ← split when shard exceeds 512 MB

struct ManifestHeader {
    uint32_t magic;           // GPKM_MAGIC
    uint16_t version;
    uint16_t _pad0;
    uint64_t generation;      // monotonically increasing, incremented on each write
    uint64_t created_at;      // unix timestamp
    uint32_t n_shards;
    uint32_t _reserved0;      // was: n_taxons
    uint64_t n_genomes;       // total (including deleted)
    uint64_t n_genomes_live;  // excluding deleted
    uint64_t catalog_size;
    uint64_t _reserved1;      // was: taxon_idx_size
    uint8_t  catalog_checksum[16];
    uint8_t  _reserved2[16]; // was: taxon_idx_checksum
    uint8_t  _pad1[32];
};
static_assert(sizeof(ManifestHeader) == 128, "ManifestHeader layout changed");

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
    std::optional<GenomeMeta> genome_meta(GenomeId id) const;

    // Count genomes matching query (no decompression)
    size_t count(const ExtractQuery& q) const;

    // Metadata-only filter result
    std::vector<GenomeMeta> filter_meta(const ExtractQuery& q) const;

    // Full extraction with FASTA decompression
    std::vector<ExtractedGenome> extract(const ExtractQuery& q) const;

    // Fetch one genome by id
    std::optional<ExtractedGenome> fetch_genome(GenomeId id) const;

    // Fetch one genome by accession (requires meta.tsv sidecar)
    std::optional<ExtractedGenome> fetch_by_accession(std::string_view accession) const;

    // Archive-wide stats
    struct ArchiveStats {
        uint64_t generation;
        size_t   n_shards;
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

    // Add all genomes from a TSV (accession, file_path [, completeness, contamination])
    void add_from_tsv(const std::filesystem::path& tsv_path);

    // Add individual record
    void add(const BuildRecord& rec);

    // Write all shards, catalog, manifest, meta.tsv sidecar
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
