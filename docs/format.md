# Binary Format

A `.gpk` file is a seekable single-file container inspired by Parquet. The `TailLocator` fixed footer at EOF points to the current TOC, so the reader needs only two seeks to open any archive. Sections are appended across generations without rewriting existing data.

---

## File layout

```
┌──────────────────────────────────────────────────────────────────────┐
│  FileHeader          128 B                                           │
│  magic · version · uuid · created_ts · generation                   │
├──────────────────────────────────────────────────────────────────────┤
│  SHRD × N  (generation 1)    compressed genome shards               │
│  SHRD × N  (generation 2+)   appended; old shards untouched         │
│  ...                                                                 │
├──────────────────────────────────────────────────────────────────────┤
│  CATL  columnar genome metadata (SoA, sorted by oph_fingerprint)    │
│  GIDX  genome_id → (section_id, dir_index, catl_row)                │
│  ACCX  FNV-1a hash table: accession → genome_id                     │
│  CIDX  sorted (FNV-1a-64(contig), genome_id) array                  │
│  TAXN  FNV-1a hash table: accession → lineage string                │
│  TXDB  full taxonomy tree (taxid/parent/rank/name + acc→taxid)      │
│  KMRX  float[n × 136]  L2-normalised k=4 tetranucleotide profiles  │
│  HNSW  hnswlib serialised blob + label map                          │
│  TOMB  tombstone records for soft-deleted genomes                   │
├──────────────────────────────────────────────────────────────────────┤
│  TOC                         zstd-compressed SectionDesc[]          │
├──────────────────────────────────────────────────────────────────────┤
│  TailLocator          64 B   fixed footer at EOF, points to TOC     │
└──────────────────────────────────────────────────────────────────────┘
```

---

## Section types

| Magic | Name | Description |
|-------|------|-------------|
| `SHRD` | Shard | Compressed genome blobs and directory |
| `CATL` | Catalog | Columnar genome metadata (SoA, sorted by oph) |
| `GIDX` | Genome index | genome_id to (section_id, dir_index, catl_row) |
| `ACCX` | Accession index | FNV-1a hash table: accession string to genome_id |
| `CIDX` | Contig index | Sorted (FNV-1a-64(contig_acc), genome_id) array |
| `TAXN` | Taxonomy strings | FNV-1a hash table: accession to lineage string |
| `TXDB` | Taxonomy tree | Parsed taxid/parent/rank/name nodes + acc-to-taxid table |
| `KMRX` | K-mer profiles | float[n × 136] L2-normalised k=4 tetranucleotide frequencies |
| `HNSW` | HNSW index | hnswlib serialised blob + label map (genome_id per vector) |
| `TOMB` | Tombstone | Soft-deleted genome_id records |

---

## `FileHeader` — 128 bytes

| Offset | Size | Field | Description |
|--------|------|-------|-------------|
| 0 | 4 B | `magic` | `GPK\x01` |
| 4 | 2 B | `version_major` | Breaking format change |
| 6 | 2 B | `version_minor` | Backward-compatible extension |
| 8 | 16 B | `uuid` | Archive UUID (stable across generations) |
| 24 | 8 B | `created_ts` | Unix timestamp of initial build |
| 32 | 8 B | `generation` | Monotonically incremented on each `add`/`rm`/`repack` |
| 40 | 88 B | _reserved_ | Zero-padded |

---

## Shard section (`SHRD`)

```
┌──────────────────────────────────────────────────────────────────────┐
│  ShardHeader         128 B                                           │
│  magic · shard_id · n_genomes · codec · dict_size                   │
│  dir_offset · dict_offset · blob_area_offset · checkpoint_offset    │
├──────────────────────────────────────────────────────────────────────┤
│  GenomeDirEntry[n_genomes]   64 B each                               │
│  genome_id · oph_fingerprint · blob_offset · blob_len_cmp/raw       │
│  checkpoint_idx · n_checkpoints  —  sorted by oph_fingerprint       │
├──────────────────────────────────────────────────────────────────────┤
│  zstd dictionary   dict_size B   (optional)                         │
├──────────────────────────────────────────────────────────────────────┤
│  Blob area                                                           │
│  blob[0] · blob[1] · ...   (each independently decompressible)      │
├──────────────────────────────────────────────────────────────────────┤
│  CheckpointEntry[]   16 B each   (optional)                         │
│  symbol_offset · block_offset   (enables sub-genome slice)          │
└──────────────────────────────────────────────────────────────────────┘
```

Genomes are sorted by `oph_fingerprint` within each shard. Nearby OPH values indicate similar k-mer content, maximising zstd LDM reuse and shared dictionary effectiveness.

### `ShardHeader` — 128 bytes

| Offset | Size | Field | Description |
|--------|------|-------|-------------|
| 0 | 4 B | `magic` | `SHRD` |
| 4 | 2 B | `version` | |
| 6 | 2 B | `codec` | Compression codec (see table below) |
| 8 | 8 B | `shard_id` | Unique shard identifier |
| 16 | 4 B | `n_genomes` | Number of genomes in this shard |
| 20 | 4 B | `dict_size` | Dictionary size in bytes (0 if none) |
| 24 | 8 B | `genome_dir_offset` | Byte offset of `GenomeDirEntry[]` from section start |
| 32 | 8 B | `dict_offset` | Byte offset of zstd dictionary |
| 40 | 8 B | `blob_area_offset` | Byte offset of blob area |
| 48 | 8 B | `checkpoint_area_offset` | Byte offset of `CheckpointEntry[]` (0 if none) |
| 56 | 72 B | _reserved_ | Zero-padded |

### `GenomeDirEntry` — 64 bytes each

| Offset | Size | Field | Description |
|--------|------|-------|-------------|
| 0 | 8 B | `genome_id` | |
| 8 | 8 B | `oph_fingerprint` | Order-preserving hash; entries sorted by this value |
| 16 | 8 B | `blob_offset` | Byte offset of compressed blob from blob area start |
| 24 | 8 B | `blob_len_cmp` | Compressed size in bytes |
| 32 | 8 B | `blob_len_raw` | Uncompressed size in bytes |
| 40 | 4 B | `checkpoint_idx` | Index into `CheckpointEntry[]` for first checkpoint |
| 44 | 4 B | `n_checkpoints` | Number of checkpoints (0 if slice not needed) |
| 48 | 16 B | _reserved_ | |

### `CheckpointEntry` — 16 bytes each (optional)

| Offset | Size | Field | Description |
|--------|------|-------|-------------|
| 0 | 8 B | `symbol_offset` | Byte offset within the decompressed genome at this checkpoint |
| 8 | 8 B | `block_offset` | Byte offset of the corresponding zstd block within the blob |

### Codec values

| Value | Name | Description |
|-------|------|-------------|
| 0 | `PLAIN` | Each blob is independent zstd |
| 1 | `ZSTD_DICT` | Shared dictionary trained on first N genomes |
| 2 | `REF_DICT` | First genome used as reference content dictionary |
| 3 | `DELTA` | Non-reference blobs zstd-compressed with `refPrefix` from genome 0 |
| 4 | `MEM_DELTA` | Seed with k=31 k-mers; store MEM list + zstd verbatim residue |

---

## Catalog section (`CATL`)

```
┌──────────────────────────────────────────────────────────────────────┐
│  CatlHeader          32 B                                            │
│  magic · n_rows · n_groups · stats_offset · rows_offset             │
├──────────────────────────────────────────────────────────────────────┤
│  RowGroupStatsV2[n_groups]   72 B each                               │
│  min/max oph · min/max completeness · min/max genome_length          │
│  (enables predicate pushdown — skip entire groups without row scan)  │
├──────────────────────────────────────────────────────────────────────┤
│  GenomeMeta[n_rows]   72 B each                                      │
│  sorted by oph_fingerprint                                           │
└──────────────────────────────────────────────────────────────────────┘
```

Multiple CATL fragments (one per generation) are merged by `MergedCatalogReader` at read time; newer fragments take precedence on duplicate `genome_id`.

### `CatlHeader` — 32 bytes

| Offset | Size | Field | Description |
|--------|------|-------|-------------|
| 0 | 4 B | `magic` | `CATL` |
| 4 | 4 B | `n_rows` | Total number of `GenomeMeta` rows |
| 8 | 4 B | `n_groups` | Number of row groups |
| 12 | 4 B | _reserved_ | |
| 16 | 8 B | `stats_offset` | Byte offset of `RowGroupStatsV2[]` from section start |
| 24 | 8 B | `rows_offset` | Byte offset of `GenomeMeta[]` from section start |

### `GenomeMeta` — 72 bytes each

| Offset | Size | Field | Description |
|--------|------|-------|-------------|
| 0 | 8 B | `genome_id` | |
| 8 | 8 B | `oph_fingerprint` | |
| 16 | 4 B | `completeness` | CheckM completeness × 100 (fixed-point) |
| 20 | 4 B | `contamination` | CheckM contamination × 100 |
| 24 | 8 B | `genome_length` | Total assembly length in bp |
| 32 | 4 B | `n_contigs` | |
| 36 | 4 B | `shard_id` | Which shard holds this genome |
| 40 | 32 B | _reserved_ | |

---

## Genome index (`GIDX`)

A flat array of fixed-size records sorted by `genome_id`, enabling O(log n) binary search. Each record maps:

```
genome_id  ->  (section_id, dir_index, catl_row_index)
```

`section_id` is resolved via the TOC to find the shard's file offset. `dir_index` is the entry's position in `GenomeDirEntry[]`. `catl_row_index` is the row in the merged catalog.

---

## TOC and TailLocator

The TOC is a zstd-compressed array of `SectionDesc` records.

### `SectionDesc`

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

### Open sequence

1. `lseek(-64, SEEK_END)` — read the 64-byte `TailLocator`
2. `lseek(toc_offset, SEEK_SET)` — read and decompress the TOC
3. Parse `SectionDesc[]` — mmap the entire file

---

## KMRX section

Stores L2-normalised k=4 canonical tetranucleotide frequency vectors (136 dimensions; reverse-complement collapsing reduces unique k-mers from 256 to 136).

| Offset | Size | Field | Description |
|--------|------|-------|-------------|
| 0 | 32 B | `KmrxHeader` | magic, n_genomes, flags |
| 32 | n × 8 B | `genome_ids[n]` | Sorted ascending; binary search for O(log n) lookup |
| 32 + n×8 | n × 544 B | `profiles[n][136]` | Parallel to `genome_ids`; stored uncompressed |

---

## HNSW section

Embeds a serialised [hnswlib](https://github.com/nmslib/hnswlib) index blob. Default build parameters: M=16, efConstruction=200.

| Offset | Size | Field | Description |
|--------|------|-------|-------------|
| 0 | 64 B | `HnswSectionHeader` | magic, n_elements, M, efConstruction |
| 64 | variable | hnswlib blob | Serialised hnswlib index |
| 64 + blob | n × 8 B | `label_map[n]` | Translates hnswlib internal label i to `genome_id` |
