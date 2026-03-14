#pragma once
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace genopack {

// ── Stable identifiers ───────────────────────────────────────────────────────

using GenomeId = uint64_t;   // globally unique, never reused
using TaxonId  = uint32_t;   // internal taxonomy node ID
using ShardId  = uint32_t;   // index into shard table in MANIFEST

static constexpr GenomeId INVALID_GENOME_ID = 0;
static constexpr TaxonId  INVALID_TAXON_ID  = 0;
static constexpr ShardId  INVALID_SHARD_ID  = UINT32_MAX;

// ── File format magic numbers & version ──────────────────────────────────────

static constexpr uint32_t GPKM_MAGIC  = 0x4D4B5047u; // "GPKM"
static constexpr uint32_t GPKC_MAGIC  = 0x434B5047u; // "GPKC"
static constexpr uint32_t GPKS_MAGIC  = 0x534B5047u; // "GPKS"
static constexpr uint32_t GPKI_MAGIC  = 0x494B5047u; // "GPKI" (taxon index)
static constexpr uint16_t FORMAT_VERSION = 1;

// ── Per-genome metadata ───────────────────────────────────────────────────────

struct GenomeMeta {
    GenomeId genome_id    = INVALID_GENOME_ID;
    TaxonId  taxon_id     = INVALID_TAXON_ID;
    ShardId  shard_id     = INVALID_SHARD_ID;

    // Sequence stats
    uint64_t genome_length    = 0;   // total bp
    uint32_t n_contigs        = 0;
    uint16_t gc_pct_x100      = 0;   // 0-10000, e.g. 5234 = 52.34%

    // Quality (from CheckM2; 0 = missing)
    uint16_t completeness_x10 = 0;   // 0-1000, e.g. 987 = 98.7%
    uint16_t contamination_x10 = 0;

    // Compact OPH fingerprint (XOR-folded from full 10k-bin sketch)
    uint64_t oph_fingerprint  = 0;

    // Location within shard
    uint64_t blob_offset  = 0;   // byte offset in .gpks blob area
    uint32_t blob_len_cmp = 0;   // compressed size
    uint32_t blob_len_raw = 0;   // uncompressed FASTA size

    // Bookkeeping
    uint32_t date_added = 0;   // days since 2024-01-01
    uint32_t flags      = 0;   // bit 0: deleted (tombstoned)

    static constexpr uint32_t FLAG_DELETED = 1u;

    bool is_deleted() const { return (flags & FLAG_DELETED) != 0; }
};
static_assert(sizeof(GenomeMeta) == 72, "GenomeMeta layout changed");

// ── Taxonomy ──────────────────────────────────────────────────────────────────

struct TaxonInfo {
    TaxonId  taxon_id   = INVALID_TAXON_ID;
    TaxonId  parent_id  = INVALID_TAXON_ID;
    uint8_t  rank       = 0;   // 0=domain … 6=species
    std::string name;          // e.g. "Escherichia coli"
    std::string full_lineage;  // GTDB string: "d__Bacteria;…;s__Escherichia coli"
};

// GTDB rank indices
enum class TaxonRank : uint8_t {
    Domain  = 0, // d__
    Life    = 1, // l__
    Kingdom = 2, // k__
    Phylum  = 3, // p__
    Class   = 4, // c__
    Order   = 5, // o__
    Family  = 6, // f__
    Genus   = 7, // g__
    Species = 8, // s__
};

// ── Shard descriptor (stored in MANIFEST) ────────────────────────────────────

struct ShardDescriptor {
    ShardId  shard_id         = INVALID_SHARD_ID;
    TaxonId  primary_taxon_id = INVALID_TAXON_ID;  // genus or higher
    uint32_t n_genomes        = 0;
    uint32_t n_deleted        = 0;
    uint64_t file_size        = 0;
    uint64_t blob_area_offset = 0;
    char     filename[64]     = {};  // relative path under shards/
    uint8_t  checksum[16]     = {};  // MD5 of shard file
};
static_assert(sizeof(ShardDescriptor) == 112, "ShardDescriptor layout changed");

// ── Extract query ─────────────────────────────────────────────────────────────

struct ExtractQuery {
    // Taxonomy filters (empty = no filter)
    std::optional<std::string> genus;
    std::optional<std::string> species;
    std::optional<std::string> family;
    std::optional<TaxonId>     taxon_id;

    // Quality filters (0 = no filter)
    float min_completeness   = 0.0f;   // %
    float max_contamination  = 100.0f; // %
    uint64_t min_genome_length = 0;
    uint64_t max_genome_length = UINT64_MAX;
    uint32_t max_contigs     = UINT32_MAX;

    // Output control
    bool include_deleted = false;
    size_t limit = SIZE_MAX;
};

// ── Per-archive result ────────────────────────────────────────────────────────

struct ExtractedGenome {
    GenomeMeta  meta;
    std::string taxonomy;       // full GTDB lineage string
    std::string fasta;          // decompressed FASTA content
};

// ── Build record (input to ArchiveBuilder) ────────────────────────────────────

struct BuildRecord {
    std::string             accession;
    std::string             taxonomy;   // full GTDB lineage
    std::filesystem::path   file_path;  // path to gzipped FASTA

    // Optional quality metadata
    float completeness  = 0.0f;
    float contamination = 0.0f;
    uint64_t genome_length = 0;
    uint32_t n_contigs     = 0;
};

} // namespace genopack
