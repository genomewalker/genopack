# genopack

A high-performance single-file genome archive format for large-scale microbial genome collections. Stores millions of FASTA genomes with compressed shards, fast accession lookup, taxonomy, tetranucleotide profiles, and a full taxonomy tree - all in one seekable `.gpk` file.

## Features

- **Single-file archive** - one `.gpk` replaces millions of per-genome files on NFS
- **Fast genome access** - O(1) accession and genome-id lookup via ACCX and GIDX, plus per-genome shard-local blob fetch
- **zstd-compressed shards** - dictionary, reference-delta, and MEM-delta modes with mmap-backed reads
- **Contig index (CIDX)** - FNV-1a sorted array mapping every contig accession → genome_id; batch lookup at ~150M queries/s
- **Taxonomy storage** - full NCBI/GTDB taxonomy tree (TXDB) + per-accession lineage strings (TAXN)
- **Taxonomy export** - NCBI taxdump (`names.dmp`/`nodes.dmp`/`acc2taxid.dmp`) and high-performance columnar binary (`acc2taxid.bin`/`taxnodes.bin`) plus TSV sidecars
- **k=4 tetranucleotide profiles** - 136-dim L2-normalised float vectors (KMRX) for fast cosine-similarity nearest-neighbour without decompressing FASTA
- **Taxonomy repack** - re-shard by genus/family for 10–13× faster per-taxon NFS access
- **Distributed build** - split TSV across N nodes, build parts in parallel, merge with parallel pwrite
- **Append and tombstone** - add genomes or mark deleted without full rebuild

## Format overview

```
FileHeader (128 B)
SHRD × N          - zstd-compressed genome shards (sorted by OPH fingerprint)
CATL              - columnar genome metadata (genome_id, shard_id, GC%, length, OPH, quality)
GIDX              - genome_id → (section_id, dir_index, catl_row) for O(1) fetch
ACCX              - FNV-1a hash table: accession string → genome_id
CIDX              - sorted (FNV-1a-64(contig_acc), genome_id) array for contig → genome lookup
TAXN              - FNV-1a hash table: accession string → full lineage string
TXDB              - full taxonomy tree (taxid/parent/rank/name nodes + acc→taxid table)
KMRX              - float[n_genomes × 136] L2-normalised k=4 tetranucleotide profiles
TOC               - section descriptor table
TailLocator (64B) - fixed-size footer at EOF, points to TOC offset (Parquet-style)
```

See the [documentation site](https://genomewalker.github.io/genopack/) for the full binary format specification, API reference, and CLI reference.

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

### Contig accession lookup

```bash
# Single lookup: contig accession → genome_id
genopack cidx mydb.gpk --accession NZ_JAVJIU010000001.1

# Batch lookup from file (one contig per line)
genopack cidx mydb.gpk --accessions-file contigs.txt --threads 8
```

### Taxonomy lookup

```bash
genopack taxonomy mydb.gpk --accession GCA_000008085.1
```

### Export taxonomy packages

```bash
# NCBI taxdump (names.dmp, nodes.dmp, acc2taxid.dmp - Kraken/Kaiju compatible)
genopack taxdump mydb.gpk -f taxdump -o ./taxdump/

# Columnar binary (acc2taxid.bin + taxnodes.bin) + TSV sidecars
genopack taxdump mydb.gpk -f columnar -o ./taxonomy/
```

The columnar binary format is designed for applications that need taxonomy without linking genopack:

| File | Size (GTDB r226 reps) | Description |
|------|----------------------|-------------|
| `acc2taxid.bin` | 2.2 MB | Sorted `(FNV-1a-64(acc), taxid)` pairs - O(log n) binary search |
| `taxnodes.bin` | 6.3 MB | Sorted node records + name pool - O(log n) taxid lookup |
| `acc2taxid.tsv` | 3.7 MB | `accession\ttaxid` - pandas/polars/R |
| `taxonomy.tsv` | 9.0 MB | `taxid\tparent_taxid\trank\tname\tis_synthetic` |

Binary file layout documented in [Taxonomy Export](../../wiki/Taxonomy-Export).

### Re-shard by taxonomy

```bash
# Re-shard an existing archive by genus (fast per-taxon NFS access)
genopack repack mydb.gpk mydb_taxon.gpk -t 24 -z 6
```

Reads only shard directory headers in a fast first pass (~minutes), then a single sequential decompression pass routes each genome to its taxonomy group. The result allows tools like geodesic to read only the shards for a target taxon instead of the entire archive.

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

Each node builds its slice to local NVMe (`/scratch`), then rsyncs the part back. Merge uses parallel `pwrite` - one thread per part archive, reading shards sequentially for NFS readahead efficiency.

See [Distributed Build](https://genomewalker.github.io/genopack/getting-started/#distributed-build) in the docs for details.

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
| `slice` | Extract a subsequence by accession and coordinates |
| `add` | Append genomes to existing archive |
| `rm` | Tombstone (soft-delete) genomes |
| `dedup` | Remove sequence duplicates |
| `taxonomy` | Query taxonomy tree |
| `taxdump` | Export taxonomy as NCBI taxdump or columnar binary |
| `similar` | Find similar genomes by KMRX cosine similarity |
| `repack` | Re-shard by taxonomy for fast per-taxon NFS access |
| `reindex` | Append missing GIDX/HNSW sections to archive |

Full option reference: [CLI Reference](https://genomewalker.github.io/genopack/cli/)

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

// Contig → genome lookup (single)
uint32_t genome_id = reader.find_contig_genome_id("NZ_JAVJIU010000001.1");

// Contig → genome lookup (batch, parallelised)
std::vector<std::string_view> contigs = { /* ... */ };
std::vector<uint32_t> genome_ids(contigs.size());
reader.batch_find_contig_genome_ids(contigs.data(), genome_ids.data(),
                                    contigs.size(), /*n_threads=*/8);

// Build
genopack::ArchiveBuilder builder("mydb.gpk");
builder.add_from_tsv("genomes.tsv");
builder.finalize();
```

## License

MIT
