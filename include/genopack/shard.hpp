#pragma once
#include "types.hpp"
#include <filesystem>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace genopack {

// ── On-disk shard format (.gpks) ─────────────────────────────────────────────
//
// Shards are immutable after build. Each covers one similarity cluster
// (~128–512 MB compressed). Genomes are ordered by oph_fingerprint within the
// shard for zstd dictionary reuse.
//
// [GPKS header 128B]
// [genome directory: GenomeDirEntry × n_genomes]
// [zstd shared dictionary (optional)]
// [blob area: compressed genome blobs, back-to-back]
//   each blob: independently zstd-decompressible FASTA chunk
//              2-bit normalized for ACGT, N-runs side channel
//              per-blob checkpoints every 64K bases for random access
// [checkpoint area: CheckpointEntry[] per genome with >64K bases]
// [footer: offsets, checksums]
//
// Genome ordering within shard:
//   Sort by oph_fingerprint (approximates sketch similarity).
//   This maximises zstd LDM reuse across nearby genomes.

struct ShardHeader {
    uint32_t magic;                    // GPKS_MAGIC
    uint16_t version;
    uint16_t _pad0;
    uint32_t shard_id;
    uint32_t n_genomes;
    uint32_t n_deleted;                // genomes with FLAG_DELETED
    uint32_t cluster_id;               // similarity cluster index
    uint32_t dict_size;                // 0 = no shared dictionary
    uint32_t _pad1;                    // explicit alignment to 8-byte boundary
    uint64_t genome_dir_offset;
    uint64_t dict_offset;
    uint64_t blob_area_offset;
    uint64_t checkpoint_area_offset;
    uint64_t footer_offset;
    uint8_t  file_checksum[16];        // MD5 of everything before footer
    uint8_t  _pad2[40];
};
static_assert(sizeof(ShardHeader) == 128, "ShardHeader layout changed");

struct GenomeDirEntry {
    GenomeId genome_id;
    uint32_t _reserved0;
    uint32_t flags;
    uint64_t blob_offset;    // byte offset within blob area
    uint32_t blob_len_cmp;   // compressed blob length (bytes)
    uint32_t blob_len_raw;   // uncompressed FASTA length (bytes)
    uint32_t n_checkpoints;  // checkpoints stored for this genome
    uint32_t checkpoint_idx; // index of first checkpoint in checkpoint area
    uint64_t oph_fingerprint;
    uint8_t  _pad[16];
};
static_assert(sizeof(GenomeDirEntry) == 64, "GenomeDirEntry layout changed");

// One checkpoint per 65536 bases. Allows random access within large genomes.
struct CheckpointEntry {
    uint64_t symbol_offset;   // base position (0-indexed)
    uint32_t block_offset;    // byte offset within the decompressed blob
    uint32_t _pad;
};
static_assert(sizeof(CheckpointEntry) == 16, "CheckpointEntry layout changed");

// ── ShardWriter ───────────────────────────────────────────────────────────────

struct ShardWriterConfig {
    int    zstd_level           = 6;
    bool   use_long_match       = true;    // zstd LDM for cross-genome repetition
    int    zstd_wlog            = 27;      // window log (1 << 27 = 128 MB)
    bool   train_dict           = true;    // build shared dictionary from first N genomes
    size_t dict_samples         = 100;     // genomes sampled for dictionary training
    size_t dict_size            = 112640;  // 110 KB dictionary
    size_t max_shard_size_bytes = 512ULL << 20; // 512 MB
};

class ShardWriter {
public:
    using Config = ShardWriterConfig;

    explicit ShardWriter(const std::filesystem::path& path, ShardId shard_id,
                         uint32_t cluster_id, Config cfg = Config{});
    ~ShardWriter();
    ShardWriter(const ShardWriter&) = delete;
    ShardWriter& operator=(const ShardWriter&) = delete;

    // Add a genome. FASTA content is the decompressed FASTA string.
    // Caller must add genomes in oph_fingerprint order for best compression.
    void add_genome(GenomeId id, uint64_t oph_fingerprint,
                    const char* fasta_data, size_t fasta_len,
                    uint32_t flags = 0);

    void finalize();  // write checkpoint area + footer + header

    size_t n_genomes() const;
    size_t compressed_size() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// ── ShardReader ───────────────────────────────────────────────────────────────

class ShardReader {
public:
    ShardReader();
    ~ShardReader();
    ShardReader(const ShardReader&) = delete;
    ShardReader& operator=(const ShardReader&) = delete;
    ShardReader(ShardReader&&) noexcept;
    ShardReader& operator=(ShardReader&&) noexcept;

    void open(const std::filesystem::path& path);
    void close();

    ShardId  shard_id()   const;
    uint32_t cluster_id() const;
    size_t   n_genomes()  const;

    // Fetch decompressed FASTA for one genome. O(1) seek + decompress.
    std::string fetch_genome(GenomeId id) const;

    // Iterate genome directory
    const GenomeDirEntry* dir_begin() const;
    const GenomeDirEntry* dir_end()   const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace genopack
