# Binary Format

A `.gpk` file is a **seekable single-file container** inspired by Parquet. Sections can be appended without rewriting existing data, and the TOC at the end is updated with each generation. A 64-byte `TailLocator` at EOF points to the current TOC, allowing the reader to open the archive with two seeks.

---

## File layout

```
┌──────────────────────────────┐
│  FileHeader  (128 B)         │  magic, version, UUID, timestamp
├──────────────────────────────┤
│  SHRD × N                    │  zstd-compressed genome shards
│  SHRD × N  (generation 2)    │  appended; old shards untouched
├──────────────────────────────┤
│  CATL                        │  columnar genome metadata (SoA)
│  GIDX                        │  genome_id → (section_id, dir_index, catl_row)
│  ACCX                        │  FNV-1a hash table: accession → genome_id
│  CIDX                        │  sorted FNV-1a-64(contig) → genome_id pairs
│  TAXN                        │  FNV-1a hash table: accession → lineage string
│  TXDB                        │  full taxonomy tree
│  KMRX                        │  float[n × 136] k=4 profiles
│  HNSW                        │  hnswlib ANN index blob
│  TOMB                        │  tombstone (soft-delete) records
├──────────────────────────────┤
│  TOC                         │  section descriptor table (zstd-compressed)
│  TailLocator  (64 B)         │  fixed-size footer: points to TOC offset
└──────────────────────────────┘
```

---

## Section types

| Magic | Name | Description |
|-------|------|-------------|
| `SHRD` | Shard | Compressed genome blobs + directory |
| `CATL` | Catalog | Columnar genome metadata (SoA, sorted by oph) |
| `GIDX` | Genome index | genome_id → (section_id, dir_index, catl_row) |
| `ACCX` | Accession index | FNV-1a hash table: accession string → genome_id |
| `CIDX` | Contig index | Sorted (FNV-1a-64(contig_acc), genome_id) array |
| `TAXN` | Taxonomy strings | FNV-1a hash table: accession → lineage string |
| `TXDB` | Taxonomy tree | Parsed taxid/parent/rank/name nodes + acc→taxid |
| `KMRX` | K-mer profiles | float[n × 136] L2-normalised k=4 tetranucleotide frequencies |
| `HNSW` | HNSW index | hnswlib serialised blob + label map (genome_id per vector) |
| `TOMB` | Tombstone | Soft-deleted genome_id records |

---

## Shard layout (`SHRD`)

Each shard section starts with a 128-byte `ShardHeader`, followed by the genome directory, an optional dictionary, and the blob area.

```
┌─────────────────────────────────────────────┐
│  ShardHeader  (128 B)                        │
│    magic, version, shard_id, n_genomes       │
│    codec, dict_size                          │
│    genome_dir_offset                         │
│    dict_offset                               │
│    blob_area_offset                          │
│    checkpoint_area_offset                    │
├─────────────────────────────────────────────┤
│  GenomeDirEntry[n_genomes]  (64 B each)      │
│    genome_id, oph_fingerprint                │
│    blob_offset, blob_len_cmp, blob_len_raw   │
│    checkpoint_idx, n_checkpoints             │
├─────────────────────────────────────────────┤
│  zstd dictionary  (dict_size B, optional)    │
├─────────────────────────────────────────────┤
│  Blob area                                   │
│    blob[0]  - independently decompressible   │
│    blob[1]                                   │
│    ...                                       │
├─────────────────────────────────────────────┤
│  CheckpointEntry[]  (16 B each, optional)    │
│    symbol_offset, block_offset               │
└─────────────────────────────────────────────┘
```

Genomes are sorted by `oph_fingerprint` within each shard. Nearby OPH values indicate similar k-mer content, which maximises zstd LDM reuse and shared dictionary effectiveness.

### Codec field

| Value | Codec | Description |
|-------|-------|-------------|
| 0 | `PLAIN` | Each blob is independent zstd |
| 1 | `ZSTD_DICT` | Shared dictionary trained on first N genomes |
| 2 | `REF_DICT` | First genome used as reference content dictionary |
| 3 | `DELTA` | Non-reference blobs zstd-compressed with `refPrefix` from genome 0 |
| 4 | `MEM_DELTA` | Seed with k=31 k-mers; store MEM list + zstd verbatim residue |

---

## Catalog section (`CATL`)

The catalog stores `GenomeMeta` rows in a columnar struct-of-arrays layout, sorted by `oph_fingerprint`. Row-group statistics (min/max OPH, completeness, genome_length) enable predicate pushdown - scans can skip entire row groups without accessing individual rows.

```
┌───────────────────────────────┐
│  CatlHeader  (32 B)           │
│    magic, n_rows, n_groups    │
│    stats_offset, rows_offset  │
├───────────────────────────────┤
│  RowGroupStatsV2[n_groups]    │  72 B each; min/max per 32768 rows
├───────────────────────────────┤
│  GenomeMeta[n_rows]           │  72 B each; in oph_fingerprint order
└───────────────────────────────┘
```

Multiple CATL fragments (one per generation) are merged by `MergedCatalogReader` at read time; newer fragments take precedence on duplicate `genome_id`.

---

## Genome index (`GIDX`)

Maps `genome_id` to its physical location for O(1) fetch:

```
genome_id  →  (section_id, dir_index, catl_row_index)
```

`section_id` is looked up in the TOC to find the shard section's file offset. `dir_index` is the position in `GenomeDirEntry[]`. `catl_row_index` is the row in the merged catalog.

---

## TOC and TailLocator

The Table of Contents is a zstd-compressed array of `SectionDesc` records. Each record describes one section:

```cpp
struct SectionDesc {
    uint32_t type;              // section magic (e.g. SEC_SHRD)
    uint16_t version;
    uint16_t flags;
    uint64_t section_id;        // unique, monotonically increasing
    uint64_t file_offset;
    uint64_t compressed_size;
    uint64_t uncompressed_size;
    uint64_t item_count;        // genomes in a shard, rows in a catalog, etc.
    uint64_t aux0;              // type-specific (shard_id for SHRD)
    uint64_t aux1;
    uint8_t  checksum[16];
};
```

The 64-byte `TailLocator` at EOF contains the TOC file offset and a file UUID, allowing the reader to open the archive with:

1. `lseek(-64, SEEK_END)` → read `TailLocator`
2. `lseek(toc_offset, SEEK_SET)` → read and decompress TOC
3. Parse section descriptors → mmap the entire file

---

## Versioning

| Field | Description |
|-------|-------------|
| `FileHeader.version_major` | Breaking format change |
| `FileHeader.version_minor` | Backward-compatible extension |
| `FileHeader.generation` | Monotonically incremented on each `add`/`rm`/`repack` |
| `SectionDesc.version` | Per-section format version (e.g. shard v4 added checkpoints) |

---

## KMRX section

Stores L2-normalised k=4 canonical tetranucleotide frequency vectors (136 dimensions, not 256, because reverse-complement collapsing reduces unique k-mers). Layout:

```
KmrxHeader (32 B)
uint64_t genome_ids[n_genomes]       ← sorted ascending
float    profiles[n_genomes × 136]   ← parallel to genome_ids
```

Lookup is O(log n) binary search on the sorted `genome_ids` array. Profiles are stored uncompressed (float data compresses poorly).

---

## HNSW section

Embeds a serialised [hnswlib](https://github.com/nmslib/hnswlib) index blob with a label map that translates hnswlib internal labels back to `genome_id` values. Default build parameters: M=16, efConstruction=200.

```
HnswSectionHeader (64 B)
hnswlib serialised index blob
uint64_t label_map[n_elements]   ← label i → genome_id
```
