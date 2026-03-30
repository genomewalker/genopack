# API Reference

All public symbols live in the `genopack` namespace. Include `<genopack/archive.hpp>` for the high-level reader/builder/appender. Low-level types are in `<genopack/types.hpp>`.

---

## Types

### `GenomeMeta`

Per-genome metadata record. 72 bytes, 8-byte aligned.

```cpp
#include <genopack/types.hpp>
```

| Field | Type | Description |
|-------|------|-------------|
| `genome_id` | `uint64_t` | Stable integer identifier (1-based) |
| `shard_id` | `uint32_t` | Shard section this genome lives in |
| `genome_length` | `uint64_t` | Total assembly length in bp |
| `n_contigs` | `uint32_t` | Number of contigs |
| `gc_pct_x100` | `uint16_t` | GC% × 100 (e.g. 5234 = 52.34%) |
| `completeness_x10` | `uint16_t` | CheckM completeness × 10 (0–1000) |
| `contamination_x10` | `uint16_t` | CheckM contamination × 10 |
| `oph_fingerprint` | `uint64_t` | MinHash minimum (k=21); similarity sort key |
| `blob_offset` | `uint64_t` | Byte offset of compressed blob within shard blob area |
| `blob_len_cmp` | `uint32_t` | Compressed blob size |
| `blob_len_raw` | `uint32_t` | Decompressed FASTA size |
| `date_added` | `uint32_t` | Days since 2024-01-01 |
| `flags` | `uint32_t` | Bitfield; bit 0 = deleted |

**Methods:**

```cpp
bool is_deleted() const;   // (flags & FLAG_DELETED) != 0
```

---

### `ExtractQuery`

Filter for `ArchiveReader::extract`, `filter_meta`, `count`, and `ScanEngine::scan_filtered`.

```cpp
struct ExtractQuery {
    uint64_t min_oph = 0;
    uint64_t max_oph = UINT64_MAX;
    float    min_completeness  = 0.0f;
    float    max_contamination = 100.0f;
    uint64_t min_genome_length = 0;
    uint64_t max_genome_length = UINT64_MAX;
    uint32_t max_contigs       = UINT32_MAX;
    std::vector<std::string> accessions;  // if non-empty, restrict to these
    bool     include_deleted   = false;
    size_t   limit             = SIZE_MAX;
};
```

---

### `ExtractedGenome`

Result of a fetch or extract operation.

```cpp
struct ExtractedGenome {
    GenomeMeta  meta;
    std::string accession;
    std::string fasta;    // decompressed FASTA (header + sequence)
};
```

---

### `BuildRecord`

Input record for `ArchiveBuilder::add`.

```cpp
struct BuildRecord {
    std::string           accession;
    std::filesystem::path file_path;
    float    completeness  = 0.0f;
    float    contamination = 0.0f;
    uint64_t genome_length = 0;
    uint32_t n_contigs     = 0;
    std::vector<std::pair<std::string, std::string>> extra_fields;
};
```

---

## ArchiveReader

High-level read-only access to a `.gpk` archive. All methods are thread-safe after `open()` (the underlying mmap is read-only).

```cpp
#include <genopack/archive.hpp>
```

### Construction

```cpp
genopack::ArchiveReader reader;
reader.open("mydb.gpk");  // memory-maps the file; cheap
```

### Fetch by accession

```cpp
// Single fetch - decompresses the genome's shard section
std::optional<ExtractedGenome> fetch_by_accession(std::string_view accession) const;

// Metadata only (no FASTA decompression) - O(1) via ACCX + CATL
std::optional<GenomeMeta> genome_meta_by_accession(std::string_view accession) const;
```

### Batch fetch

```cpp
// Batch fetch: groups by shard, reads each shard once.
// results[i] = nullopt if accession not found.
std::vector<std::optional<ExtractedGenome>>
batch_fetch_by_accessions(const std::vector<std::string>& accessions) const;
```

### Streaming shard visitors

These are the primary interface for large-scale genome processing. Each shard is read once and its pages released.

```cpp
// Streaming visitor: cb(original_index, genome) in shard order.
// Peak memory: O(max_genomes_per_shard).
void visit_by_shard(
    const std::vector<std::string>& accessions,
    const std::function<void(size_t idx, ExtractedGenome)>& cb) const;

// Batch visitor: delivers all genomes from one shard at a time.
// Allows OMP parallelism within each shard batch.
using ShardBatch = std::vector<std::pair<size_t, ExtractedGenome>>;
void visit_shard_batches(
    const std::vector<std::string>& accessions,
    const std::function<void(ShardBatch&)>& cb) const;
```

Example - embedding with OMP within each batch:

```cpp
reader.visit_shard_batches(accessions, [&](auto& batch) {
    #pragma omp parallel for num_threads(16)
    for (int i = 0; i < (int)batch.size(); ++i) {
        auto& [idx, genome] = batch[i];
        results[idx] = embed(genome.fasta);
    }
});
```

### Taxonomy

```cpp
// Lineage string for one accession (from TAXN section)
std::optional<std::string> taxonomy_for_accession(std::string_view accession) const;

// Iterate all (accession, taxonomy) pairs
void scan_taxonomy(const std::function<void(std::string_view acc,
                                            std::string_view tax)>& cb) const;

// Full taxonomy tree from TXDB section
std::optional<TaxonomyTree> taxonomy_tree() const;
```

### Contig lookup

```cpp
// Single contig accession → genome_id (UINT32_MAX if not found)
uint32_t find_contig_genome_id(std::string_view contig_acc) const;

// Batch lookup - parallelised over n_threads
// out_genome_ids[i] = genome_id for accs[i], or UINT32_MAX if not found
void batch_find_contig_genome_ids(const std::string_view* accs,
                                  uint32_t* out_genome_ids,
                                  size_t n,
                                  size_t n_threads = 0) const;
```

### KMRX similarity

```cpp
// Pointer to float[136] L2-normalised k=4 profile, nullptr if not stored
const float* kmer_profile(GenomeId id) const;
const float* kmer_profile_by_accession(std::string_view accession) const;

// Linear scan: up to k most similar genomes by cosine similarity
std::vector<std::pair<GenomeId, float>>
find_similar(GenomeId id, size_t k = 10) const;

std::vector<std::pair<GenomeId, float>>
find_similar_by_accession(std::string_view accession, size_t k = 10) const;
```

### Genome iteration

```cpp
// Reverse lookup: genome_id → accession string
std::string accession_for_genome_id(GenomeId id) const;

// Iterate all (accession, genome_id) pairs
void scan_genome_accessions(
    const std::function<void(std::string_view, GenomeId)>& cb) const;

// Low-level shard scan (rarely needed)
void scan_shards(const std::function<void(const uint8_t* data,
                                          uint64_t offset,
                                          uint64_t compressed_size,
                                          uint32_t shard_id)>& cb) const;
```

### Statistics

```cpp
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
```

---

## ArchiveBuilder

Builds a new archive from scratch. Not thread-safe.

```cpp
#include <genopack/archive.hpp>
```

### Config

```cpp
struct ArchiveBuilderConfig {
    ShardWriterConfig shard_cfg;
    size_t io_threads           = 16;    // parallel FASTA decompression
    bool   verbose              = false;
    size_t batch_size           = 1000;  // genomes buffered before shard flush
    GenomeId starting_genome_id = 1;
    bool   build_hnsw           = true;  // disable with --no-hnsw for > 1M genomes
    bool   build_cidx           = true;  // disable with --no-cidx for large archives
    bool   kmer_nn_sort         = false; // sort within shards by kmer NN chain
    bool   taxonomy_group       = false; // bucket by taxonomy before shard formation
    std::string taxonomy_rank   = "g";   // "g" = genus, "f" = family
};
```

### Usage

```cpp
genopack::ArchiveBuilder builder("mydb.gpk");
builder.add_from_tsv("genomes.tsv");   // or call add() per record
builder.finalize();
```

```cpp
// Add individual records
BuildRecord rec;
rec.accession = "GCA_000008085.1";
rec.file_path = "/data/GCA_000008085.1.fna.gz";
rec.completeness = 98.7f;
rec.contamination = 0.3f;
builder.add(rec);
```

---

## ArchiveAppender

Appends genomes to an existing archive without full rebuild.

```cpp
genopack::ArchiveAppender app("mydb.gpk");
app.add_from_tsv("new_genomes.tsv");
app.commit();  // writes new shard generation, updates catalog
```

Soft-delete genomes:

```cpp
app.remove_by_accession("GCA_000001405.1");
app.commit();
```

---

## ScanEngine

High-throughput streaming scan optimised for NFS and NVMe. Decouples I/O threads (sequential reads) from worker threads (decompression + processing).

```cpp
#include <genopack/scan_engine.hpp>
```

### Config

```cpp
struct ScanEngine::Config {
    size_t io_threads     = 1;    // dedicated I/O threads
    size_t worker_threads = 15;   // decompression + processing threads
    size_t slab_size_mb   = 64;   // per-slab buffer (should match shard size)
    size_t prefetch_depth = 3;    // shards prefetched ahead of current
    bool   nfs_mode       = true; // sequential I/O; false = NVMe parallel reads
};
```

### Usage

```cpp
genopack::ScanEngine engine({.io_threads = 2, .worker_threads = 22});

engine.scan_all(reader, [](GenomeId id, std::string_view fasta,
                            const genopack::GenomeMeta& meta) {
    // called concurrently from worker threads
    process(id, fasta);
});

// Filtered scan - applies CATL + taxonomy pre-filter before decompression
genopack::ExtractQuery q;
q.min_completeness = 95.0f;
q.max_contamination = 5.0f;
engine.scan_filtered(reader, q, callback);
```

---

## repack_archive

Re-shards an existing archive by taxonomy for fast per-taxon NFS access.

```cpp
#include <genopack/repack.hpp>
```

```cpp
struct RepackConfig {
    char              taxonomy_rank    = 'g';         // 'g'=genus, 'f'=family
    ShardWriterConfig shard_cfg;
    size_t            threads          = 1;            // OMP decompression threads
    uint64_t          max_bucket_bytes = 32ULL << 30; // eviction cap (32 GB)
    bool              verbose          = false;
};

void repack_archive(const std::filesystem::path& input_gpk,
                    const std::filesystem::path& output_gpk,
                    const RepackConfig& cfg = {});
```

The three-phase algorithm:

| Phase | I/O | Description |
|-------|-----|-------------|
| 1 - Directory scan | ~300 MB | Reads only `GenomeDirEntry` headers; builds genome→taxonomy index |
| 2 - Sort | 0 | Sorts all records by `(taxonomy, oph_fingerprint)` in memory |
| 3 - Decompress + route | Full archive | Single sequential pass; smart eviction of largest writer at memory cap |

---

## ShardReader

Zero-copy, mmap-backed shard reader. Used internally by `ArchiveReader` but exposed for low-level access.

```cpp
#include <genopack/shard.hpp>
```

```cpp
ShardReader shard;
shard.open(base_ptr, section_offset, section_size);

// Fetch by directory index - O(1), reads only the one genome's blob
std::string fasta = shard.fetch_genome_at(dir_index);

// Iterate directory entries (no decompression)
for (auto* de = shard.dir_begin(); de != shard.dir_end(); ++de) {
    // de->genome_id, de->oph_fingerprint, de->blob_len_raw, ...
}

// Release mmap page cache after batch fetch
shard.release_pages();
```

### `GenomeDirEntry`

64-byte on-disk directory entry within a shard section.

| Field | Type | Description |
|-------|------|-------------|
| `genome_id` | `uint64_t` | Stable genome identifier |
| `oph_fingerprint` | `uint64_t` | MinHash minimum (k=21) |
| `blob_offset` | `uint64_t` | Offset relative to shard blob area start |
| `blob_len_cmp` | `uint32_t` | Compressed blob bytes |
| `blob_len_raw` | `uint32_t` | Decompressed FASTA bytes |
| `checkpoint_idx` | `uint32_t` | First checkpoint entry index |
| `n_checkpoints` | `uint32_t` | Number of checkpoint entries |
| `flags` | `uint32_t` | Codec flags |
| `meta_row_id` | `uint32_t` | Row index in the logical catalog |

---

## KmrxReader

Direct access to the KMRX k=4 tetranucleotide profile section.

```cpp
#include <genopack/kmrx.hpp>

KmrxReader kr;
kr.open(base_ptr, section_offset, section_size);

// Returns pointer to float[136] L2-normalised profile, nullptr if not found
const float* p = kr.profile_for(genome_id);

// Cosine similarity between two profiles
float sim = genopack::cosine_similarity(p1, p2);
```

Profiles are stored uncompressed as `float[n_genomes × 136]` with sorted `genome_ids[]` for O(log n) binary search.
