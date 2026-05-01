# Binary Format

A `.gpk` archive is a **directory of seekable section files** plus a `toc.bin` entry, inspired by Parquet. Sections can be appended without rewriting existing data, and the TOC is updated with each generation; a 64-byte `TailLocator` at the end of `toc.bin` points to the current TOC. Multipart sets (used for very large or distributed builds) are directories containing one or more `part_*.gpk` archives — readers like `ArchiveSetReader` open them transparently.

> **Single-file layout.** Earlier versions stored everything in one `.gpk` file with `FileHeader · sections · TOC · TailLocator`. The directory layout below is the current on-disk form; the section-level binary structures (CATL, GIDX, ACCX, …) are unchanged.

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
    <span class="fmt-name">SKCH<span class="fmt-sub">OPH (one-permutation hash) sketches — V4 seekable, dual-seed, multi-k</span></span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name">NTDB<span class="fmt-sub">embedded NCBI taxonomy (nodes.dmp + names.dmp) — optional</span></span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name" style="font-style:italic;font-weight:400">KMRX<span class="fmt-sub">optional — float[n × 136] L2-normalised k=4 tetranucleotide profiles (library-only)</span></span>
  </div>
  <div class="fmt-row">
    <span class="fmt-name" style="font-style:italic;font-weight:400">HNSW<span class="fmt-sub">optional — hnswlib serialised blob + label map (library-only, not built by default)</span></span>
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
| `SKCH` | OPH sketches | One-permutation hash signatures (V4 seekable; dual-seed; multi-k) |
| `NTDB` | NCBI taxdump | Embedded `nodes.dmp` + `names.dmp` (zstd) for offline taxid resolution |
| `KMRX` | K-mer profiles | *Library-only.* float[n × 136] L2-normalised k=4 tetranucleotide frequencies |
| `HNSW` | HNSW index | *Library-only, optional.* hnswlib serialised blob + label map |
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
  <div class="fmt-row parent">
    <div class="fmt-header">
      <span class="fmt-name">SHRD</span>
      <span class="fmt-size">variable</span>
    </div>
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
        <span class="fmt-name" style="font-style:italic;font-weight:400">zstd dictionary<span class="fmt-sub">optional — ZSTD_DICT and REF_DICT codecs only</span></span>
        <span class="fmt-size">dict_size B</span>
      </div>
      <div class="fmt-row">
        <span class="fmt-name">Blob area<span class="fmt-sub">blob[0] · blob[1] · ... — each independently decompressible</span></span>
        <span class="fmt-size">variable</span>
      </div>
      <div class="fmt-row">
        <span class="fmt-name" style="font-style:italic;font-weight:400">CheckpointEntry[]<span class="fmt-sub">optional — symbol_offset · block_offset · enables sub-genome slice</span></span>
        <span class="fmt-size">16 B × n</span>
      </div>
    </div>
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
  <div class="fmt-row parent">
    <div class="fmt-header">
      <span class="fmt-name">CATL</span>
      <span class="fmt-size">variable</span>
    </div>
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

> **Library-only / optional.** Not built by default. There is no `genopack` CLI to build or query HNSW; readers (e.g. `ArchiveReader::find_similar`) accelerate cosine-similarity search if a section is present and fall back to linear scan otherwise.

Embeds a serialised [hnswlib](https://github.com/nmslib/hnswlib) index blob. Default build parameters: M=16, efConstruction=200.

| Offset | Size | Field | Description |
|--------|------|-------|-------------|
| 0 | 64 B | `HnswSectionHeader` | magic, n_elements, M, efConstruction |
| 64 | variable | hnswlib blob | Serialised hnswlib index |
| 64 + blob | n × 8 B | `label_map[n]` | Translates hnswlib internal label `i` to `genome_id` |

---

## SKCH section (V4 seekable)

Stores one-permutation-hash (OPH) MinHash signatures used for fast Jaccard estimation. V4 is dual-seed (two parallel sketches per genome at independent seeds), supports multiple k-mer sizes in a single section, and is laid out as zstd-compressed seekable frames so a reader can stream a single frame without decompressing the whole section.

> Backward compatibility: V1/V2/V3 archives are **rejected** by current readers. Use `genopack reindex --skch --force` to rebuild SKCH on legacy archives.

### Header constants

| Constant | Value | Meaning |
|---|---|---|
| `SKCH_V4_MAGIC` | `'SKC4'` (`0x34434B53`) | Section magic |
| `SKCH_V4_FRAME_SIZE` | `16384` | Genomes per frame (last frame may be shorter) |
| `seed1` (default) | `42` | Primary OPH seed |
| `seed2` (default) | `43` | Dual seed (independent second sketch, **not** a densification seed) |

### `SkchSeekHdr` — 96 bytes

| Offset | Size | Field | Description |
|---|---|---|---|
| 0   | 4 B  | `magic`         | `SKCH_V4_MAGIC` |
| 4   | 4 B  | `n_frames`      | Number of zstd frames |
| 8   | 4 B  | `frame_size`    | Genomes per frame (16 384 by default) |
| 12  | 4 B  | `n_genomes`     | Total genomes covered |
| 16  | 4 B  | `sketch_size`   | OPH bins per signature |
| 20  | 4 B  | `n_kmer_sizes`  | Number of k values stored (≤ 8) |
| 24  | 32 B | `kmer_sizes[8]` | Sorted ascending; tail zero-padded |
| 56  | 4 B  | `syncmer_s`     | Open-syncmer prefilter `s` (0 = disabled) |
| 60  | 4 B  | `mask_words`    | `ceil(sketch_size / 64)` |
| 64  | 8 B  | `seed1`         | Primary seed |
| 72  | 8 B  | `seed2`         | Dual seed |
| 80  | 16 B | _reserved_      | Zero-padded |

### `SkchFrameDesc` — 16 bytes (× `n_frames`)

| Offset | Size | Field | Description |
|---|---|---|---|
| 0 | 8 B | `data_offset`     | Byte offset to compressed frame from section payload start |
| 8 | 4 B | `compressed_size` | On-disk size of the zstd-compressed frame |
| 12| 4 B | `n_genomes`       | Genomes covered by this frame (≤ `frame_size`) |

### Section payload layout

```
[SkchSeekHdr                                     (96 B)]
[SkchFrameDesc       × n_frames                  (16 B each)]
[uint64_t genome_ids[n_genomes]                          ] ← binary search w/o decompression
[uint64_t genome_lengths[n_genomes]                      ] ← parallel to genome_ids
[Frame 0 .. Frame n_frames-1   — independent zstd frames ]
```

Each decompressed frame is **planar by k-mer size**:

```
[uint32_t n_real_bins[n_kmer_sizes × frame_n]]
[uint16_t sigs1      [n_kmer_sizes × frame_n × sketch_size]]   seed = seed1
[uint16_t sigs2      [n_kmer_sizes × frame_n × sketch_size]]   seed = seed2
[uint64_t masks      [n_kmer_sizes × frame_n × mask_words]]
```

`genome_ids[]` is uncompressed, so a reader can binary-search a target genome, compute the frame index, then decompress that frame only.

---

## NTDB section

Optional. Embeds the raw NCBI taxdump (`nodes.dmp` + `names.dmp`) as two consecutive zstd-compressed blobs so an archive can resolve NCBI taxids fully offline (e.g. on a compute node without `/cvmfs` or network).

### `NtdbHeader` — 64 bytes

| Offset | Size | Field | Description |
|---|---|---|---|
| 0  | 4 B  | `magic`            | `SEC_NTDB` |
| 4  | 2 B  | `version`          | 1 |
| 6  | 2 B  | `flags`            | Reserved |
| 8  | 8 B  | `taxdump_date`     | `YYYYMMDD` (0 if unknown) |
| 16 | 8 B  | `nodes_raw_size`   | Uncompressed size of `nodes.dmp` |
| 24 | 8 B  | `nodes_zstd_size`  | Compressed size; blob immediately follows the header |
| 32 | 8 B  | `names_raw_size`   | Uncompressed size of `names.dmp` |
| 40 | 8 B  | `names_zstd_size`  | Compressed size; blob follows the nodes blob |
| 48 | 16 B | _reserved_         | Zero-padded |

Section layout: `[NtdbHeader (64 B)] [zstd(nodes.dmp)] [zstd(names.dmp)]`.

Built by `genopack coordinator --ntdb DIR` or via library `NtdbWriter`. Read with `NtdbReader::nodes_dmp()` / `names_dmp()` (decompresses on demand).

---

# `.gpd` — Geodesic Derep Archive Format v1

A `.gpd` file is the on-disk artefact produced by [geodesic](https://github.com/genomewalker/geodesic) when it dereplicates a genopack archive set. It is a **single file** (not a directory) and is consumed read-only by genopack via the `DerepView` API. The format is designed to (a) be cheap to mmap, (b) reuse genopack's accession universe verbatim, and (c) detect staleness against the source pack without re-running derep.

## File layout

```
[GpdFileHeader                64 B]
[Section blobs ...] (HDR · ASTR · ASOF · ARMP · RTBL · G2RM · EMBD · [CSTAT])
[TOC                                ]
[GpdTailLocator              16 B   ]
```

All multi-byte integers are little-endian; section payloads are 8-byte aligned (zero padding).

### `GpdFileHeader` — 64 bytes (offset 0)

| Offset | Size | Field | Description |
|---|---|---|---|
| 0  | 4 B  | `magic`        | `'GPDF'` (`0x46445047`) |
| 4  | 2 B  | `format_major` | 1 |
| 6  | 2 B  | `format_minor` | 0 |
| 8  | 8 B  | `toc_offset`   | Byte offset of TOC section |
| 16 | 8 B  | `toc_size`     | TOC byte size |
| 24 | 40 B | _reserved_     | Zero-padded |

### `GpdTailLocator` — 16 bytes (last bytes of file)

| Offset | Size | Field | Description |
|---|---|---|---|
| 0  | 8 B | `toc_offset` | Duplicate of header `toc_offset` (corruption check) |
| 8  | 4 B | `magic`      | `'GPDT'` (`0x54445047`) |
| 12 | 4 B | `crc32`      | crc32 of TOC bytes |

### Section magics

| Magic | Name | Description |
|---|---|---|
| `GHDR` | HDR  | Identity, params, source-part fingerprints (always uncompressed) |
| `GAST` | ASTR | Concatenated accession string pool, ASCIIbetically sorted |
| `GASO` | ASOF | `uint32_t offsets[n_genomes+1]` boundaries into ASTR |
| `GARM` | ARMP | Optional open-addressed accession→ordinal hash map |
| `GRTB` | RTBL | Rep table (one entry per representative) |
| `G2RM` | G2RM | `uint32_t rep_id[n_genomes]` (sentinels: `0xFFFFFFFE` unclustered, `0xFFFFFFFF` tombstoned) |
| `GEMB` | EMBD | Rep-only embedding matrix (default f16 × dim, typically 256) |
| `GCST` | CSTAT| Optional per-cluster QC statistics |
| `GTOC` | TOC  | Section descriptor table |

Sections may be zstd-compressed (`flags & 1`); the reader decompresses transparently and the TOC carries both compressed and uncompressed sizes.

### `GpdHeader` (HDR section)

```c
struct GpdHeader {
    uint32_t magic;            // 'GPDH'
    uint16_t format_major;     // 1
    uint16_t format_minor;     // 0
    uint64_t created_at_unix;
    uint8_t  run_id[16];       // UUID v4 — unique per derep run
    uint16_t n_parts;          // source pack parts at derep time
    uint16_t embedding_dim;    // typically 256
    uint8_t  embedding_dtype;  // 0=f32, 1=f16
    uint8_t  has_cstats;       // 1 if CSTAT present
    uint8_t  pad0[2];
    uint64_t n_genomes;        // total genomes covered
    uint64_t n_reps;
    uint64_t n_unclustered;
    // followed by:
    //   GpdSourcePart parts[n_parts];
    //   GpdDerepParams params;
};

struct GpdSourcePart {           // 48 bytes
    uint8_t  archive_uuid[16];   // from genopack header
    uint64_t generation;
    uint64_t n_genomes_total;
    uint64_t n_genomes_live;
    uint64_t accession_set_hash; // xxh3-64 of '\n'-joined sorted live accessions
};

struct GpdDerepParams {
    uint8_t  n_kmer_sizes;
    uint8_t  kmer_sizes[7];      // up to 7; tail zero-padded
    uint32_t sketch_size;
    uint64_t sig1_seed;
    uint64_t sig2_seed;
    float    jaccard_thresh;
    uint16_t geodesic_ver_len;
    uint8_t  pad1[2];
    char     geodesic_ver[];     // not null-terminated; padded to 8
};
```

### `accession_set_hash`

The operational identity of a source part: take all **live** (non-tombstoned) accessions, sort ASCIIbetically, join with `\n` (no trailing newline), hash with xxh3-64. Both the writer (geodesic, at derep time) and the reader (genopack `DerepView::check`) compute it the same way.

### Staleness levels

Returned by `DerepView::check(pack)`:

| Level | Trigger |
|---|---|
| `Valid`                      | All parts: same `accession_set_hash`, `archive_uuid`, `generation` |
| `LayoutChangedSameLiveSet`   | Same live accession set, but UUID/generation differs (e.g. repacked) |
| `StaleNewGenomes`            | Pack contains accessions absent from the .gpd |
| `StaleTombstones`            | .gpd contains accessions no longer live in the pack |
| `Mismatch`                   | Structural difference (missing part, etc.) |

### Versioning

`format_major` bumps for incompatible changes; `format_minor` bumps for additive ones (new optional sections). Readers must reject unknown `format_major` and accept higher `format_minor` (ignoring unknown sections).
