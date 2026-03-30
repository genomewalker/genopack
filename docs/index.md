# genopack

**High-performance single-file genome archive for large-scale microbial collections.**

genopack stores millions of FASTA genomes with compressed shards, O(1) accession and genome-id lookup, full taxonomy, tetranucleotide profiles, and an ANN index - all in one seekable `.gpk` file optimised for NFS and NVMe.

---

## Quick start

```bash
# Build an archive from a TSV
genopack build -i genomes.tsv -o mydb.gpk -t 24

# Extract by accession
genopack extract mydb.gpk --accession GCA_000008085.1 -o genome.fasta

# Find similar genomes by k=4 profile cosine similarity
genopack similar mydb.gpk GCA_000008085.1 -k 20 --min-sim 0.95

# Re-shard by taxonomy (10× faster per-taxon access on NFS)
genopack repack mydb.gpk mydb_taxon.gpk -t 24 -z 6
```

```cpp
#include <genopack/archive.hpp>

genopack::ArchiveReader reader;
reader.open("mydb.gpk");

// Fetch one genome
auto g = reader.fetch_by_accession("GCA_000008085.1");
if (g) std::cout << g->fasta;

// Stream all genomes shard-by-shard (peak memory = one shard)
reader.visit_shard_batches(accessions, [](auto& batch) {
    for (auto& [idx, genome] : batch) { /* ... */ }
});
```

---

## Features

| Feature | Description |
|---------|-------------|
| **Single file** | One `.gpk` replaces millions of per-genome files on NFS |
| **zstd shards** | Dictionary, reference-delta, and MEM-delta compression modes |
| **O(1) lookup** | FNV-1a hash table (ACCX) + direct index (GIDX) for genome fetch |
| **Contig index** | CIDX maps every contig accession → genome_id at ~150M queries/s |
| **Taxonomy** | Per-accession lineage strings (TAXN) + full tree (TXDB); NCBI taxdump export |
| **k=4 profiles** | 136-dim L2-normalised tetranucleotide vectors (KMRX) for ANN similarity |
| **HNSW index** | Approximate nearest-neighbour index on KMRX profiles (hnswlib, M=16) |
| **Taxonomy repack** | Re-shard by genus/family for fast per-taxon NFS access |
| **Distributed build** | Split TSV across N nodes, build parts, merge with parallel pwrite |
| **Append / tombstone** | Add genomes or soft-delete without full rebuild |

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
