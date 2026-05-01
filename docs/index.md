# genopack

**High-performance single-file genome archive for large-scale microbial collections.**

genopack stores millions of FASTA genomes with compressed shards, O(1) accession and genome-id lookup, full taxonomy, multi-k OPH sketches, and tetranucleotide profiles — all in one seekable `.gpk` archive (a directory containing `toc.bin` + section files) or a multipart set of `part_*.gpk` archives read transparently via `ArchiveSetReader`.

---

## Quick start

```bash
# Build an archive from a TSV (defaults: per-taxon shards, multi-stage codec, OPH sketches)
genopack build -i genomes.tsv -o mydb.gpk -t 24 -z 6

# Multi-k SKCH in one pass (e.g. for geodesic dual-k)
genopack build -i genomes.tsv -o mydb.gpk -t 24 -z 6 --sketch-kmers 16,21

# Extract by accession (works on a single archive or a multipart directory)
genopack extract mydb.gpk --accession GCA_000008085.1 -o genome.fasta

# Re-shard by taxonomy (10× faster per-taxon access on NFS)
genopack repack mydb.gpk mydb_taxon.gpk -t 24 -z 6

# Inspect SKCH layout and preload cost
genopack inspect mydb.gpk
```

```cpp
#include <genopack/archive.hpp>
#include <genopack/archive_set_reader.hpp>

genopack::ArchiveSetReader reader;
reader.open("mydb.gpk");        // archive directory; or directory of part_*.gpk

auto g = reader.fetch_by_accession("GCA_000008085.1");
if (g) std::cout << g->fasta;
```

---

## Features

| Feature | Description |
|---------|-------------|
| **Single archive** | One `.gpk` directory (`toc.bin` + section files) replaces millions of per-genome files on NFS |
| **Multipart sets** | Directories containing `part_*.gpk` are read transparently as one logical archive |
| **zstd shards** | Auto, shared-dict, reference-delta, MEM-delta, and 2-bit-pack codecs |
| **O(1) lookup** | FNV-1a hash table (ACCX) + direct index (GIDX) for genome fetch |
| **Contig index** | CIDX maps every contig accession → genome_id at ~150M queries/s |
| **Taxonomy** | Per-accession lineage strings (TAXN) + full tree (TXDB); NCBI taxdump + columnar binary export |
| **OPH sketches (SKCH v4)** | Dual-seed (42, 43) one-permutation-hash signatures over one or more k; seekable zstd frames (16 384 genomes/frame) |
| **k=4 profiles** | 136-dim L2-normalised tetranucleotide vectors (KMRX); cosine similarity via library API |
| **HNSW index** | Optional hnswlib serialised blob over KMRX (library only; not built by default) |
| **Taxonomy repack** | Re-shard by genus/family for fast per-taxon NFS access |
| **Distributed build** | Split TSV across N nodes; merge with parallel `pwrite` or coordinate via NFS manifest |
| **Append / tombstone** | Add genomes or soft-delete without full rebuild |
| **`.gpd` derep view** | Read derep state produced by [geodesic](https://github.com/genomewalker/geodesic) — `accession → rep_id`, `rep_id → embedding`, with staleness detection |

---

## Build

**Dependencies:** cmake ≥ 3.20, zstd, zlib, C++20 compiler. `libdeflate` optional (3× faster gzip input).

```bash
git clone https://github.com/genomewalker/genopack
cd genopack
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

!!! note "Cluster builds"
    The default flags include `-mavx2 -mfma`. Do **not** use `-march=native` on a login node if compute nodes have a different microarchitecture.

### Install Python docs dependencies

```bash
pip install mkdocs-material
mkdocs serve   # live preview at http://localhost:8000
mkdocs build   # static site → site/
```
