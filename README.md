# genopack

A high-performance single-file genome archive format for large-scale microbial genome collections. Stores millions of FASTA genomes with compressed shards, fast accession lookup, taxonomy, tetranucleotide profiles, and a full taxonomy tree — all in one seekable `.gpk` file.

## Features

- **Single-file archive** — one `.gpk` replaces millions of per-genome files on NFS
- **Fast random access** — O(1) genome lookup by accession (ACCX hash table) or genome ID (GIDX index)
- **zstd-compressed shards** — 3–4× compression, zero-copy mmap reads
- **Taxonomy storage** — full NCBI/GTDB taxonomy tree (TXDB) + per-accession lineage strings (TAXN)
- **k=4 tetranucleotide profiles** — 136-dim L2-normalised float vectors (KMRX) for fast cosine-similarity nearest-neighbour without decompressing FASTA
- **Distributed build** — split TSV across N nodes, build parts in parallel, merge with parallel pwrite
- **Append and tombstone** — add genomes or mark deleted without full rebuild

## Format overview

```
FileHeader (128 B)
SHRD × N          — zstd-compressed genome shards (sorted by OPH fingerprint)
CATL              — columnar genome metadata (genome_id, shard_id, GC%, length, OPH, quality)
GIDX              — genome_id → (section_id, dir_index, catl_row) for O(1) fetch
ACCX              — FNV-1a hash table: accession string → genome_id
TAXN              — FNV-1a hash table: accession string → full lineage string
TXDB              — full taxonomy tree (taxid/parent/rank/name nodes + acc→taxid table)
KMRX              — float[n_genomes × 136] L2-normalised k=4 tetranucleotide profiles
TOC               — section descriptor table
TailLocator (64B) — fixed-size footer at EOF, points to TOC offset (Parquet-style)
```

See the [wiki](../../wiki) for the full binary format specification.

## Build

**Dependencies:** `cmake ≥ 3.20`, `zstd`, `zlib`, `libdeflate` (optional, faster gzip), C++20 compiler

```bash
git clone https://github.com/genomewalker/genopack
cd genopack
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

> **Cluster note:** compile with `-mavx2 -mfma` (already the default). Do not use `-march=native` on a login node if compute nodes are a different microarchitecture.

## Quick start

### Build an archive

```bash
# Input TSV: accession  file_path  [completeness  contamination  taxonomy  ...]
genopack build -i genomes.tsv -o mydb.gpk -t 24 -z 3
```

### Archive statistics

```bash
genopack stat mydb.gpk
```

### Extract by accession

```bash
genopack extract mydb.gpk --accession GCA_000008085.1 -o out.fasta
genopack extract mydb.gpk --accessions-file accessions.txt -o out.fasta
```

### Find similar genomes (cosine similarity on KMRX profiles)

```bash
genopack similar mydb.gpk GCA_000008085.1 -k 20 --min-sim 0.95
```

### Taxonomy lookup

```bash
genopack taxonomy mydb.gpk --accession GCA_000008085.1
```

### Merge archives

```bash
# From file list (recommended for large merges)
ls parts/part_*.gpk > parts.txt
genopack merge -l parts.txt -o merged.gpk
```

### Append genomes

```bash
genopack add mydb.gpk -i new_genomes.tsv
```

### Remove genomes

```bash
genopack rm mydb.gpk GCA_000001405 GCA_000002655
```

## Distributed build

For collections too large for a single node, use the distributed build script:

```bash
.scripts/gpk-build-distributed.sh genomes.tsv output.gpk \
    -t 24 -z 3 \
    node1 node2 node3 node4
```

Each node builds its slice to local NVMe (`/scratch`), then rsyncs the part back. Merge uses parallel `pwrite` — one thread per part archive, reading shards sequentially for NFS readahead efficiency.

See [Distributed Build](../../wiki/Distributed-Build) in the wiki for details.

## Input TSV format

| Column | Required | Description |
|--------|----------|-------------|
| `accession` | ✓ | Unique genome identifier |
| `file_path` | ✓ | Path to `.fa`, `.fna`, `.fa.gz`, or `.fna.gz` |
| `completeness` | | CheckM completeness % (0–100) |
| `contamination` | | CheckM contamination % |
| `taxonomy` | | Full lineage string (`d__;p__;c__;o__;f__;g__;s__`) |
| any extras | | Stored verbatim in `meta.tsv` sidecar |

## CLI reference

| Command | Description |
|---------|-------------|
| `build` | Build new archive from TSV |
| `merge` | Merge multiple `.gpk` archives (parallel) |
| `stat` | Show archive statistics |
| `extract` | Extract genomes to FASTA |
| `add` | Append genomes to existing archive |
| `rm` | Tombstone (soft-delete) genomes |
| `dedup` | Remove sequence duplicates |
| `taxonomy` | Query taxonomy tree |
| `similar` | Find similar genomes by KMRX cosine similarity |
| `reindex` | Append missing GIDX section to archive |

Full option reference: [CLI Reference](../../wiki/CLI-Reference)

## Library usage

```cpp
#include <genopack/archive.hpp>

// Read
genopack::ArchiveReader reader;
reader.open("mydb.gpk");

auto genome = reader.fetch_by_accession("GCA_000008085.1");
if (genome) {
    // genome->fasta contains the decompressed FASTA string
}

auto taxon = reader.taxonomy_for_accession("GCA_000008085.1");

// Build
genopack::ArchiveBuilder builder("mydb.gpk");
builder.add_from_tsv("genomes.tsv");
builder.finalize();
```

## License

MIT
