# Binary Format

A `.gpk` file is a **seekable single-file container** inspired by Parquet. Sections can be appended without rewriting existing data, and the TOC at the end is updated with each generation. A 64-byte `TailLocator` at EOF points to the current TOC, allowing the reader to open the archive with two seeks.

---

## File layout

<div class="fmt-box">
  <div class="fmt-row">
    <span class="fmt-name">FileHeader<span class="fmt-sub">magic · version · UUID · created_ts · generation</span></span>
    <span class="fmt-size">128 B</span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name">SHRD × N<span class="fmt-sub">compressed genome shards — each generation appended, old shards untouched</span></span>
    <span class="fmt-size">variable</span>
  </div>
  <div class="fmt-row divider">
    <span class="fmt-name">CATL<span class="fmt-sub">columnar genome metadata (SoA, sorted by oph_fingerprint)</span></span>
    <span class="fmt-desc"></span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name">GIDX<span class="fmt-sub">genome_id → (section_id, dir_index, catl_row)</span></span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name">ACCX<span class="fmt-sub">FNV-1a hash table: accession → genome_id</span></span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name">CIDX<span class="fmt-sub">sorted (FNV-1a-64(contig), genome_id) array</span></span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name">TAXN<span class="fmt-sub">FNV-1a hash table: accession → lineage string</span></span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name">TXDB<span class="fmt-sub">full taxonomy tree (taxid/parent/rank/name + acc→taxid)</span></span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name">KMRX<span class="fmt-sub">float[n × 136] L2-normalised k=4 tetranucleotide profiles</span></span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name">HNSW<span class="fmt-sub">hnswlib serialised blob + label map</span></span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name">TOMB<span class="fmt-sub">tombstone records for soft-deleted genomes</span></span>
  </div>
  <div class="fmt-row divider">
    <span class="fmt-name">TOC<span class="fmt-sub">zstd-compressed SectionDesc[] — one record per section</span></span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name">TailLocator<span class="fmt-sub">fixed footer at EOF, points to TOC offset</span></span>
    <span class="fmt-size">64 B</span>
  </div>
</div>

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

<div class="fmt-box">
  <div class="fmt-row">
    <span class="fmt-name">ShardHeader<span class="fmt-sub">magic · shard_id · n_genomes · codec · dict_size<br>dir_offset · dict_offset · blob_area_offset · checkpoint_offset</span></span>
    <span class="fmt-size">128 B</span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name">GenomeDirEntry[n_genomes]<span class="fmt-sub">genome_id · oph_fingerprint · blob_offset · blob_len_cmp · blob_len_raw<br>checkpoint_idx · n_checkpoints — sorted by oph_fingerprint</span></span>
    <span class="fmt-size">64 B × n</span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name" style="font-style:italic;font-weight:400">zstd dictionary (optional)<span class="fmt-sub">present only for ZSTD_DICT and REF_DICT codecs</span></span>
    <span class="fmt-size">dict_size B</span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name">Blob area<span class="fmt-sub">blob[0] · blob[1] · ... — each independently decompressible</span></span>
    <span class="fmt-size">variable</span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name" style="font-style:italic;font-weight:400">CheckpointEntry[] (optional)<span class="fmt-sub">symbol_offset · block_offset — enables sub-genome slice without full decompress</span></span>
    <span class="fmt-size">16 B × n</span>
  </div>
</div>

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

<div class="fmt-box">
  <div class="fmt-row">
    <span class="fmt-name">CatlHeader<span class="fmt-sub">magic · n_rows · n_groups · stats_offset · rows_offset</span></span>
    <span class="fmt-size">32 B</span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name">RowGroupStatsV2[n_groups]<span class="fmt-sub">min/max oph · min/max completeness · min/max genome_length<br>enables predicate pushdown — skip entire groups without row scan</span></span>
    <span class="fmt-size">72 B × n</span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name">GenomeMeta[n_rows]<span class="fmt-sub">sorted by oph_fingerprint</span></span>
    <span class="fmt-size">72 B × n</span>
  </div>
</div>

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

## Versioning

| Field | Description |
|-------|-------------|
| `FileHeader.version_major` | Breaking format change |
| `FileHeader.version_minor` | Backward-compatible extension |
| `FileHeader.generation` | Monotonically incremented on each `add`/`rm`/`repack` |
| `SectionDesc.version` | Per-section format version (e.g. shard v4 added checkpoints) |

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
