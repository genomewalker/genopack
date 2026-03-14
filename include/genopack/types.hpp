#pragma once
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace genopack {

// ── Stable identifiers ───────────────────────────────────────────────────────

using GenomeId = uint64_t;
using ShardId  = uint32_t;

static constexpr GenomeId INVALID_GENOME_ID = 0;
static constexpr ShardId  INVALID_SHARD_ID  = UINT32_MAX;

// ── File format magic numbers & version ──────────────────────────────────────

static constexpr uint32_t GPKM_MAGIC    = 0x4D4B5047u; // "GPKM"
static constexpr uint32_t GPKC_MAGIC    = 0x434B5047u; // "GPKC"
static constexpr uint32_t GPKS_MAGIC    = 0x534B5047u; // "GPKS"
static constexpr uint16_t FORMAT_VERSION = 2;

// ── Per-genome metadata ───────────────────────────────────────────────────────
// Taxonomy-free. Shard grouping is by MinHash similarity, not taxon label.

struct GenomeMeta {
    GenomeId genome_id    = INVALID_GENOME_ID;
    uint32_t _reserved0   = 0;    // was: TaxonId taxon_id (kept for binary compat)
    ShardId  shard_id     = INVALID_SHARD_ID;

    uint64_t genome_length     = 0;   // total bp
    uint32_t n_contigs         = 0;
    uint16_t gc_pct_x100       = 0;   // 0–10000, e.g. 5234 = 52.34%
    uint16_t completeness_x10  = 0;   // 0–1000,  e.g. 987  = 98.7%
    uint16_t contamination_x10 = 0;

    // MinHash minimum (k=21 k-mer); locality-sensitive similarity sort key.
    // Genomes with nearby oph_fingerprint values share k-mer content.
    uint64_t oph_fingerprint   = 0;

    uint64_t blob_offset  = 0;
    uint32_t blob_len_cmp = 0;
    uint32_t blob_len_raw = 0;

    uint32_t date_added = 0;   // days since 2024-01-01
    uint32_t flags      = 0;   // bit 0: deleted
    // Note: 6 bytes of implicit padding exist before oph_fingerprint (compiler-inserted)
    // Total struct size: 72 bytes (alignment 8)

    static constexpr uint32_t FLAG_DELETED = 1u;
    bool is_deleted() const { return (flags & FLAG_DELETED) != 0; }
};
static_assert(sizeof(GenomeMeta) == 72, "GenomeMeta layout changed");

// ── Shard descriptor (stored in MANIFEST) ────────────────────────────────────

struct ShardDescriptor {
    ShardId  shard_id   = INVALID_SHARD_ID;
    uint32_t cluster_id = 0;      // similarity cluster index
    uint32_t n_genomes  = 0;
    uint32_t n_deleted  = 0;
    uint64_t file_size        = 0;
    uint64_t blob_area_offset = 0;
    char     filename[64]     = {};
    uint8_t  checksum[16]     = {};
};
static_assert(sizeof(ShardDescriptor) == 112, "ShardDescriptor layout changed");

// ── Extract query ─────────────────────────────────────────────────────────────

struct ExtractQuery {
    // MinHash similarity range filter
    uint64_t min_oph = 0;
    uint64_t max_oph = UINT64_MAX;

    // Quality filters
    float    min_completeness  = 0.0f;
    float    max_contamination = 100.0f;
    uint64_t min_genome_length = 0;
    uint64_t max_genome_length = UINT64_MAX;
    uint32_t max_contigs       = UINT32_MAX;

    // Accession-based lookup (resolved via meta.tsv sidecar)
    std::vector<std::string> accessions;

    bool   include_deleted = false;
    size_t limit           = SIZE_MAX;
};

// ── Extraction result ────────────────────────────────────────────────────────

struct ExtractedGenome {
    GenomeMeta  meta;
    std::string accession;
    std::string fasta;
};

// ── Build record (input to ArchiveBuilder) ────────────────────────────────────

struct BuildRecord {
    std::string           accession;
    std::filesystem::path file_path;

    float    completeness  = 0.0f;
    float    contamination = 0.0f;
    uint64_t genome_length = 0;
    uint32_t n_contigs     = 0;

    // Extra TSV columns written verbatim to meta.tsv sidecar
    std::vector<std::pair<std::string, std::string>> extra_fields;
};

} // namespace genopack
