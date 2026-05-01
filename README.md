# genopack

A high-performance single-file genome archive format for large-scale microbial genome collections. Stores millions of FASTA genomes with compressed shards, fast accession lookup, taxonomy, multi-k OPH sketches, and a full taxonomy tree — all in one seekable `.gpk` archive (a directory containing `toc.bin` plus section files) or a multipart set of `part_*.gpk` archives read transparently as one.

## Features

- **Single-file archive** — one `.gpk` directory replaces millions of per-genome files on NFS
- **Multipart sets** — directories containing `part_0.gpk`, `part_1.gpk`, … are read transparently via `ArchiveSetReader`
- **Fast genome access** — O(1) accession and genome-id lookup via ACCX and GIDX, plus per-genome shard-local blob fetch
- **zstd-compressed shards** — auto, shared-dict, reference-delta, MEM-delta, and 2-bit-pack codecs with mmap-backed reads
- **OPH sketches (SKCH v4)** — dual-seed (42, 43) one-permutation-hash signatures over one or more k-mer sizes, stored as seekable zstd frames (16 384 genomes per frame), so a `sketch_for_ids` query decompresses only the frames it touches
- **Contig index (CIDX)** — FNV-1a sorted array mapping every contig accession → genome_id; batch lookup at ~150M queries/s
- **Taxonomy storage** — full NCBI/GTDB taxonomy tree (TXDB) + per-accession lineage strings (TAXN)
- **Taxonomy export** — NCBI taxdump (`names.dmp`/`nodes.dmp`/`acc2taxid.dmp`) and high-performance columnar binary (`acc2taxid.bin`/`taxnodes.bin`) plus TSV sidecars
- **k=4 tetranucleotide profiles (KMRX)** — 136-dim L2-normalised float vectors for cosine similarity (library API; no CLI)
- **Taxonomy repack** — re-shard by genus/family for 10–13× faster per-taxon NFS access
- **Distributed build** — split TSV across N nodes, build parts in parallel, merge or coordinate via NFS manifest
- **Append and tombstone** — add genomes or mark deleted without full rebuild
- **`.gpd` derep archives** — read derep state produced by [geodesic](https://github.com/genomewalker/geodesic) via `DerepView`: O(1) `accession → rep_id`, O(1) `rep_id → embedding`, with staleness detection against the source pack

## Format overview

```
toc.bin             - sections + tail locator
SHRD × N            - zstd-compressed genome shards (sorted by OPH fingerprint)
CATL                - columnar genome metadata (genome_id, shard_id, GC%, length, OPH, quality)
GIDX                - genome_id → (section_id, dir_index, catl_row) for O(1) fetch
ACCX                - FNV-1a hash table: accession string → genome_id
CIDX                - sorted (FNV-1a-64(contig_acc), genome_id) array for contig → genome lookup
TAXN                - FNV-1a hash table: accession string → full lineage string
TXDB                - full taxonomy tree (taxid/parent/rank/name nodes + acc→taxid table)
SKCH × N            - OPH sketches: dual-seed sigs + occupancy masks, seekable zstd frames
KMRX (optional)     - float[n × 136] L2-normalised k=4 tetranucleotide profiles
HNSW (optional)     - hnswlib serialised blob for cosine ANN over KMRX (library only)
NTDB (optional)     - embedded NCBI nodes.dmp/names.dmp tree (set by `coordinator --ntdb`)
TOMB                - tombstones for soft-deleted genomes
TailLocator (64 B)  - fixed footer at end of toc.bin pointing to the TOC offset
```

See the [documentation site](https://genomewalker.github.io/genopack/) for the full binary format specification, API reference, and CLI reference. The `.gpd` derep archive format is documented in [`DEREP_FORMAT.md`](DEREP_FORMAT.md) and in the binary format docs.

## Build

**Dependencies** (resolved automatically via CMake FetchContent): CLI11, spdlog, BS::thread_pool, hnswlib, Catch2, rapidgzip, Eigen3, xxHash.

**System dependencies:** `cmake ≥ 3.20`, `zstd`, `zlib`, `libdeflate` (optional, faster gzip), C++20 compiler.

```bash
git clone https://github.com/genomewalker/genopack
cd genopack
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

> **Cluster note:** the default flags include `-mavx2 -mfma`. Do not use `-march=native` on a login node if compute nodes are a different microarchitecture.

## Quick start

### Build an archive

```bash
# Input TSV: accession  file_path  [completeness  contamination  taxonomy  ...]
genopack build -i genomes.tsv -o mydb.gpk -t 24 -z 6
```

The output `mydb.gpk` is a directory containing `toc.bin` and section files. Defaults: per-taxon shard grouping (`--taxon-group`, genus rank), kmer-NN sort within shards, OPH sketches (k=16, sketch size 10 000), CIDX contig index. Disable individual stages with `--no-taxon-group`, `--no-kmer-sort`, `--no-sketch`, `--no-cidx`. Use `--sketch-kmers 16,21,31` to build a multi-k SKCH section in one pass.

### Archive statistics

```bash
genopack stat mydb.gpk            # text
genopack stat mydb.gpk --json     # JSON
```

Works on a single archive directory or on a directory of `part_*.gpk` (multipart set).

### Inspect SKCH layout and preload cost

```bash
genopack inspect mydb.gpk           # per-archive sketch summary + estimated preload bytes
genopack inspect parts_dir/         # iterate all part_*.gpk under a directory
```

### Extract by accession

```bash
genopack extract mydb.gpk --accession GCA_000008085.1 -o out.fasta
genopack extract mydb.gpk --accessions-file accessions.txt -o out.fasta

# With quality filter
genopack extract mydb.gpk --min-completeness 95 --max-contamination 5 -o filtered.fasta

# One file per genome
genopack extract mydb.gpk --accessions-file list.txt --output-dir per_genome/
```

### Slice a sub-region

```bash
genopack slice mydb.gpk GCA_000008085.1 --start 100000 --length 5000 --fasta
```

Decompresses only the checkpoint blocks covering the requested range.

### Taxonomy

```bash
genopack taxonomy show mydb.gpk --accession GCA_000008085.1
genopack taxonomy show mydb.gpk --json

# Normalize a build TSV against a taxdump
genopack taxonomy normalize -i raw.tsv -o normalized.tsv --ncbi-taxdump taxdump/

# Partition a normalized TSV into N balanced parts for distributed build
genopack taxonomy partition -i normalized.tsv -n 8 -o parts/ -r f

# Assign stable taxids
genopack taxonomy assign-taxids -i normalized.tsv -o registry.tsv

# Diff against a new GTDB release
genopack taxonomy diff --current normalized.tsv --gtdb new_*.tsv -o diff/

# Apply a curated patch
genopack taxonomy patch --patch patch.tsv --archive mydb.gpk
```

### Export taxonomy packages

```bash
# NCBI taxdump (names.dmp, nodes.dmp, acc2taxid.dmp — Kraken/Kaiju compatible)
genopack taxdump mydb.gpk -f taxdump -o ./taxdump/

# Columnar binary (acc2taxid.bin + taxnodes.bin) + TSV sidecars
genopack taxdump mydb.gpk -f columnar -o ./taxonomy/
```

The columnar binary export is designed for applications that need taxonomy without linking genopack:

| File | Size (GTDB r226 reps) | Description |
|------|----------------------|-------------|
| `acc2taxid.bin` | 2.2 MB | Sorted `(FNV-1a-64(acc), taxid)` pairs — O(log n) binary search |
| `taxnodes.bin` | 6.3 MB | Sorted node records + name pool — O(log n) taxid lookup |
| `acc2taxid.tsv` | 3.7 MB | `accession\ttaxid` — pandas/polars/R |
| `taxonomy.tsv` | 9.0 MB | `taxid\tparent_taxid\trank\tname\tis_synthetic` |

### Re-shard by taxonomy

```bash
genopack repack mydb.gpk mydb_taxon.gpk -t 24 -z 6 -m 32
```

Reads only shard directory headers in a fast first pass (~minutes), then a single sequential decompression pass routes each genome to its taxonomy group. The result allows tools like `geodesic` to read only the shards for a target taxon instead of the entire archive.

### Merge archives

```bash
ls parts/part_*.gpk > parts.txt
genopack merge -l parts.txt -o merged.gpk
```

### Append, remove, and dedup

```bash
genopack add mydb.gpk -i new_genomes.tsv
genopack rm  mydb.gpk GCA_000001405 GCA_000002655
genopack dedup mydb.gpk --dry-run
```

### Add or rebuild indexes

```bash
# Add OPH sketches at multiple k (e.g. for geodesic dual-k)
genopack reindex mydb.gpk --skch --sketch-kmers 16,21 --skch-threads 16

# Build the full taxonomy tree (TXDB) from existing TAXN lineage strings
genopack reindex mydb.gpk --txdb

# Build CIDX from the original build TSV
genopack reindex mydb.gpk --cidx genomes.tsv --cidx-threads 16

# Force rebuild of any of the above
genopack reindex mydb.gpk --skch --force
```

## Distributed build

For collections too large for a single node, use the NFS manifest coordinator:

```bash
# 1. Start coordinator on any node with NFS access
genopack coordinator -o /nfs/output.gpk --nfs-dir /nfs/manifest/ --workers 4 \
    --ntdb /path/to/ncbi_taxdump/

# 2. On each worker (sbatch / parallel / pdsh)
genopack build -i part_N.tsv -o /scratch/part_N.gpk -t 24 -z 6 \
    --coordinator /nfs/manifest/:/nfs/output.gpk
```

Each worker builds its slice locally, then transfers sections into the shared output via `pwrite()` at coordinator-allocated offsets. Coordination is done through manifest files (`.pending` / `.alloc` / `.done`) on the shared filesystem; no TCP connections required.

`--ntdb` makes the coordinator embed the NCBI tree as an NTDB section in the final archive, so workers do not each need to parse `nodes.dmp`/`names.dmp`.

Alternatively, build parts independently and merge:

```bash
for i in 0 1 2 3; do
    genopack build -i part_$i.tsv -o parts/part_$i.gpk -t 24 -z 6
done
genopack merge -l <(ls parts/*.gpk) -o merged.gpk
```

For workflows that prefer reading parts directly without a final merge, every CLI command that takes an archive path also accepts a directory containing `part_*.gpk`.

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
| `build` | Build a new archive from TSV |
| `merge` | Merge multiple `.gpk` archives (parallel `pwrite`) |
| `stat` | Show archive statistics |
| `inspect` | Report SKCH layout and preload memory cost |
| `extract` | Extract genomes to FASTA (single file or per-accession files) |
| `slice` | Extract a subsequence by accession and coordinates |
| `add` | Append genomes to an existing archive |
| `rm` | Tombstone (soft-delete) genomes |
| `dedup` | Detect and tombstone duplicate genomes |
| `taxonomy show` | Look up taxonomy for an accession |
| `taxonomy normalize` | Normalize an input TSV against an NCBI taxdump |
| `taxonomy partition` | Split a TSV into N balanced parts for distributed build |
| `taxonomy assign-taxids` | Assign stable taxids to a normalized TSV |
| `taxonomy diff` | Diff against a new GTDB release |
| `taxonomy patch` | Apply a curated taxonomy patch to an archive |
| `taxdump` | Export NCBI taxdump or columnar binary taxonomy |
| `repack` | Re-shard by taxonomy for fast per-taxon NFS access |
| `reindex` | Append or rebuild GIDX / TXDB / CIDX / SKCH sections |
| `coordinator` | NFS manifest coordinator for distributed build |

Full option reference: [CLI Reference](https://genomewalker.github.io/genopack/cli/).

## Library usage

```cpp
#include <genopack/archive.hpp>
#include <genopack/archive_set_reader.hpp>

// Single archive or multipart set — same API
genopack::ArchiveSetReader reader;
reader.open("mydb.gpk");           // archive directory
// or
reader.open("path/to/parts_dir/"); // dir with part_*.gpk

auto genome = reader.fetch_by_accession("GCA_000008085.1");
if (genome) {
    // genome->fasta contains the decompressed FASTA string
}

auto taxon = reader.taxonomy_for_accession("GCA_000008085.1");

// Sub-region slice
auto region = reader.fetch_sequence_slice_by_accession(
    "GCA_000008085.1", /*start=*/100'000, /*length=*/5'000);
```

```cpp
// Single-archive reader exposes lower-level surfaces (KMRX, SKCH, contig index)
genopack::ArchiveReader reader;
reader.open("mydb.gpk");

// Contig → genome lookup (single)
uint32_t genome_id = reader.find_contig_genome_id("NZ_JAVJIU010000001.1");

// Contig → genome lookup (batch, parallelised)
std::vector<std::string_view> contigs = { /* ... */ };
std::vector<uint32_t> genome_ids(contigs.size());
reader.batch_find_contig_genome_ids(contigs.data(), genome_ids.data(),
                                    contigs.size(), /*n_threads=*/8);

// k=4 profile (136-dim L2-normalised); nullptr if KMRX absent
const float* p = reader.kmer_profile_by_accession("GCA_000008085.1");
```

```cpp
// Build
genopack::ArchiveBuilder builder("mydb.gpk");
builder.add_from_tsv("genomes.tsv");
builder.finalize();
```

```cpp
// Read derep state produced by geodesic
#include <genopack/derep_view.hpp>

auto derep = genopack::DerepView::open("run.gpd");

// Detect drift relative to the current pack
auto level = derep.check_against(reader);   // Valid / LayoutChangedSameLiveSet / Stale… / Mismatch

auto status = derep.rep_status("GCA_000008085.1");
if (status.kind == genopack::DerepView::RepStatus::Member) {
    auto rep = derep.rep(status.rep_id);
    const void* emb = derep.embedding(status.rep_id);   // dim × dtype bytes
}
```

## License

MIT
