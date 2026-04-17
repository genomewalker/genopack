#pragma once
#include "mmap_file.hpp"
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
//   each blob: independently decompressible FASTA payload
// [reserved for future checkpoint/footer data]
//
// Genome ordering within shard:
//   Sort by oph_fingerprint (approximates sketch similarity).
//   This maximises zstd LDM reuse across nearby genomes.

// Shard section header. All offsets are relative to the shard section start.
// Layout: 32B scalar fields, 7x8B offsets/counters, checksum(16), reserved(24).
struct ShardHeader {
    uint32_t magic;                      // GPKS_MAGIC
    uint16_t version;                    // 2=legacy blobs, 3=chunk checkpoints
    uint16_t flags;
    uint32_t shard_id;
    uint32_t cluster_id;
    uint32_t n_genomes;
    uint32_t n_deleted;
    uint32_t codec;                      // 0=zstd, 1=zstd+dict, 2=ref-dict, 3=delta, 4=mem-delta
    uint32_t dict_size;                  // bytes, 0 if no dict
    uint64_t genome_dir_offset;          // relative to shard section start
    uint64_t dict_offset;                // relative to shard section start
    uint64_t blob_area_offset;           // relative to shard section start
    uint64_t checkpoint_area_offset;     // relative to shard section start, 0 if absent
    uint64_t checkpoint_count;           // total CheckpointEntry records in shard
    uint64_t shard_raw_bp;               // total uncompressed genome bytes
    uint64_t shard_compressed_bytes;     // total compressed bytes
    uint8_t  checksum[16];
    uint8_t  reserved[24];
};
static_assert(sizeof(ShardHeader) == 128, "ShardHeader layout changed");

// Directory entry. Layout: 8+8+8+4+4+4+4+4+4+16=64 bytes.
struct GenomeDirEntry {
    uint64_t genome_id;
    uint64_t oph_fingerprint;
    uint64_t blob_offset;      // relative to shard blob_area_offset
    uint32_t blob_len_cmp;
    uint32_t blob_len_raw;
    uint32_t checkpoint_idx;
    uint32_t n_checkpoints;
    uint32_t flags;
    uint32_t meta_row_id;      // row index in logical catalog
    uint8_t  reserved[16];
};
static_assert(sizeof(GenomeDirEntry) == 64, "GenomeDirEntry layout changed");

// Reserved for future intra-genome checkpoints.
struct CheckpointEntry {
    uint64_t symbol_offset;   // base position (0-indexed)
    uint32_t block_offset;    // byte offset within the decompressed blob
    uint32_t _pad;            // packed slice metadata: low16=seq_byte_offset, high16=line_bases
};
static_assert(sizeof(CheckpointEntry) == 16, "CheckpointEntry layout changed");

struct ShardWriterConfig {
    int    zstd_level           = 6;
    bool   use_long_match       = false;   // LDM only helps in streaming mode; per-blob has no cross-genome context
    int    zstd_wlog            = 27;      // window log used when use_long_match=true
    bool   train_dict           = true;    // build shared dictionary from first N genomes
    size_t dict_samples         = 100;     // genomes sampled for dictionary training (train_dict only)
    size_t dict_size            = 112640;  // 110 KB dictionary (train_dict only)
    bool   use_reference_dict   = false;   // use first genome as reference content dict
    size_t ref_dict_max_size    = 4 << 20; // max reference dict size: 4 MB
    bool   use_delta            = false;   // compress non-reference blobs against the first genome
                                           // using zstd refPrefix; requires reference to be decompressed
                                           // before fetching any delta blob (cached in ShardReader)
    bool   auto_codec           = false;   // choose between plain chunked zstd and delta per shard
    bool   use_mem_delta        = false;   // MEM-delta: seed with k=31 k-mers, store MEM list +
                                           // zstd-compressed verbatim residue; codec=4 in ShardHeader
    size_t mem_delta_ref_panel  = 5;       // number of reference genomes forming the panel
                                           // (stored as plain zstd; all others delta against
                                           //  their concatenated sequences); for codec=4,
                                           //  ShardHeader.dict_size is repurposed as n_ref_genomes
    size_t checkpoint_bases     = 65536;   // target sequence bases per independently-compressed chunk
    size_t max_shard_size_bytes = 512ULL << 20; // 512 MB
    size_t compress_threads     = 0;       // 0 = hardware_concurrency
    bool   use_2bit_pack        = false;   // pack nucleotide sequence to 2 bits/base before zstd;
                                           // skipped per-blob if ambiguous base fraction > 5%
};

// Fully serialized shard section bytes returned by ShardWriter::freeze().
// The bytes can be appended directly to an AppendWriter in the calling thread
// while the freeze ran in a background thread.
struct FrozenShard {
    std::vector<uint8_t> bytes;     // header + dir + dict + compressed blobs
    uint32_t             shard_id;
    uint32_t             n_genomes;
    uint64_t             raw_bytes; // total uncompressed genome bytes
};

// ── ShardWriter ─────────────────────────────────────────────────────────────
// Two-pass writer: buffers raw FASTA, trains a shared zstd dictionary at
// finalize() time, then compresses and writes into an AppendWriter section.

class ShardWriter {
public:
    using Config = ShardWriterConfig;

    explicit ShardWriter(uint32_t shard_id, uint32_t cluster_id, Config cfg = {});
    ~ShardWriter();
    ShardWriter(const ShardWriter&) = delete;
    ShardWriter& operator=(const ShardWriter&) = delete;

    // Buffer raw FASTA for this genome. oph_fingerprint is the MinHash value.
    void add_genome(GenomeId id, uint64_t oph_fingerprint,
                    const char* fasta_data, size_t fasta_len,
                    uint32_t flags = 0);

    // Serialize: train dict (if enabled), compress all genomes, return fully
    // serialized shard bytes. Safe to call from a background thread.
    FrozenShard freeze();

    // Finalize: calls freeze() then appends to writer. Returns section start offset.
    uint64_t finalize(AppendWriter& writer);

    size_t n_genomes()   const;
    size_t n_bytes_raw() const;  // sum of all raw FASTA sizes

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// ── ShardReader ─────────────────────────────────────────────────────────────
// Reads a shard section from a memory-mapped .gpk file or a standalone .gpks
// file. All accesses are zero-copy into the mapped region.

class ShardReader {
public:
    ShardReader();
    ~ShardReader();
    ShardReader(ShardReader&&) noexcept;
    ShardReader& operator=(ShardReader&&) noexcept;
    ShardReader(const ShardReader&) = delete;
    ShardReader& operator=(const ShardReader&) = delete;

    // Open from a memory-mapped .gpk file at the given shard section.
    // base: start of the entire file mapping.
    // section_offset: byte offset of this shard section in the file.
    // section_size: byte length of this shard section.
    void open(const uint8_t* base, uint64_t section_offset, uint64_t section_size);

    // Open a standalone .gpks file (loads into owned memory).
    void open_file(const std::filesystem::path& path);

    bool     is_open()   const;
    uint32_t shard_id()  const;
    uint32_t n_genomes() const;

    // Fetch decompressed FASTA for a genome. Throws if not found.
    std::string fetch_genome(GenomeId id) const;

    // Fetch by directory index in O(1).
    std::string fetch_genome_at(uint32_t dir_index) const;

    // Fetch a subsequence on pure-sequence coordinates [start, start + length).
    std::string fetch_sequence_slice(GenomeId id, uint64_t start, uint64_t length) const;
    std::string fetch_sequence_slice_at(uint32_t dir_index, uint64_t start, uint64_t length) const;

    // Returns raw concatenated sequence (no FASTA headers/newlines) and contig boundary
    // positions. contig_ends[i] is the exclusive end (in seq coords) of contig i.
    // The contig_ends pointer is valid until the next call to this method or open().
    // For v4 shards: seq_buf gets the pure sequence; contig_ends points into the mmap.
    // Fallback (v2/v3 or codec=4): seq_buf stripped from FASTA; contig_ends is nullptr.
    void fetch_sequence_at_into(uint32_t dir_index,
                                std::string& seq_buf,
                                const uint32_t*& contig_ends,
                                uint32_t& n_contigs) const;

    // Release mmap page cache for the entire shard section (MADV_DONTNEED).
    // Call after all genomes from this shard are fetched in a batch to
    // avoid evicting pages mid-batch. No-op if not mmap-backed.
    void release_pages() const noexcept;

    // Iterate directory entries.
    const GenomeDirEntry* dir_begin() const;
    const GenomeDirEntry* dir_end()   const;
    const GenomeDirEntry* dir_entry(uint32_t dir_index) const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace genopack
