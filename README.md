# genopack

A single-file genome archive format for large-scale microbial genome collections. Stores millions of FASTA genomes with compressed shards, fast accession lookup, taxonomy, multi-k OPH sketches, and a full taxonomy tree — all in one seekable `.gpk` archive (a directory containing `toc.bin` plus section files) or a multipart set of `part_*.gpk` archives read transparently as one.

## Features

- **Single-file archive** — one `.gpk` directory replaces millions of per-genome files on NFS
- **Multipart sets** — directories containing `part_0.gpk`, `part_1.gpk`, … are read transparently via `ArchiveSetReader`
- **Fast genome access** — O(1) accession and genome-id lookup via ACCX and GIDX, plus per-genome shard-local blob fetch
- **zstd-compressed shards** — auto, shared-dict, reference-delta, MEM-delta, and 2-bit-pack codecs with mmap-backed reads
- **OPH sketches (SKCH v4)** — dual-seed (42, 43) one-permutation-hash signatures over one or more k-mer sizes, stored as seekable zstd frames (16 384 genomes per frame), so a `sketch_for_ids` query decompresses only the frames it touches
- **Contig index (CIDX)** — FNV-1a sorted array mapping every contig accession → genome_id; batch lookup at ~150M queries/s
- **Taxonomy storage** — full NCBI/GTDB taxonomy tree (TXDB) + per-accession lineage strings (TAXN)
- **Taxonomy export** — NCBI taxdump (`names.dmp`/`nodes.dmp`/`acc2taxid.dmp`) and columnar binary (`acc2taxid.bin`/`taxnodes.bin`) plus TSV sidecars
- **k=4 tetranucleotide profiles (KMRX)** — 136-dim L2-normalised float vectors for cosine similarity (library API; no CLI)
- **Taxonomy repack** — re-shard by genus/family for 10–13× faster per-taxon NFS access
- **Distributed build** — split TSV across N nodes, build parts in parallel, merge or coordinate via NFS manifest
- **Append and tombstone** — add genomes or mark deleted without full rebuild
- **Quality scoring (`genopack check`)** — per-genome completeness and contamination signals (QUAL section): cluster-relative completeness, leakage, TNF excess, chromosome skew closure, contig-level Mahalanobis outlier (CCO), SPE, sibling outlier (family-vs-genus), marker-gene completeness/redundancy, FracMinHash minority fraction
- **Covariance sections (GCOV/FCOV)** — per-genus and per-family TNF covariance matrices built in one pass at build time or via `genopack gcov`; used by `check` for contamination detection
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
GSTX (optional)     - per-genus sketch stats: TNF centroid, p90 completeness, OPH consensus
GCOV (optional)     - per-genus biological covariance (Ledoit-Wolf TNF, eigenvectors, SPE thresholds)
FCOV (optional)     - per-family biological covariance (same layout as GCOV, keyed by family hash)
FMHR (optional)     - per-genus FracMinHash reference sketches (k=21, c=125)
QUAL (optional)     - per-genome quality records (56 B each): completeness, contamination, flags
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

The output `mydb.gpk` is a directory containing `toc.bin` and section files. Defaults: per-taxon shard grouping (`--taxon-group`, genus rank), kmer-NN sort within shards, OPH sketches (k=16, sketch size 10 000). The CIDX contig index is opt-in (`--cidx`). Disable default stages with `--no-taxon-group`, `--no-kmer-sort`, `--no-sketch`. Use `--sketch-kmers 16,21,31` to build a multi-k SKCH section in one pass.

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

The columnar binary export is for applications that need taxonomy without linking genopack:

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

Reads only shard directory headers in a fast first pass (~minutes), then a single sequential decompression pass routes each genome to its taxonomy group. Tools like `geodesic` then read only the shards for a target taxon, not the whole archive.

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

For large collections on a **shared NFS filesystem**, partition by genus, build parts in parallel, and query the directory directly — no merge or coordinator needed:

```bash
# 1. Partition by genus into N balanced parts (LPT bin-packing, genus stays whole)
genopack taxonomy partition -i genomes.tsv -n 6 -o parts/ -r g

# 2. Build each part on a separate node in parallel
genopack build -i parts/part_0.tsv -o parts/part_0.gpk -p 16 -t 8
genopack build -i parts/part_1.tsv -o parts/part_1.gpk -p 8  -t 8
# ...

# 3. All commands accept the parts directory directly
genopack check parts/ --genomes query.txt -o quality.tsv -t 16
genopack stat  parts/
```

Every CLI command that takes an archive path also accepts a directory containing `part_*.gpk`. Genus-partitioning keeps each genus whole so GSTX is built from complete genus membership per part.

For disconnected nodes with local scratch, use the NFS coordinator instead:

```bash
genopack coordinator -o /nfs/output.gpk --nfs-dir /nfs/manifest/ --workers 4
genopack build -i part_N.tsv -o /scratch/part_N.gpk -t 24 -z 6 \
    --coordinator /nfs/manifest/:/nfs/output.gpk
```

To produce a single merged archive:

```bash
genopack merge -l <(ls parts/*.gpk) -o merged.gpk
```

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
| `verify` | Verify XXH128 checksums for all sections |
| `check` | Compute per-genome quality scores; writes QUAL section and TSV |
| `gcov` | Build or rebuild GCOV + FCOV + FMHR covariance sections (one pass) |
| `calibrate` | Fit isotonic+OLS completeness calibration model vs CheckM2 |
| `markers build` | Build a `.mrk` single-copy marker panel from a GTDB-Tk database |
| `markers score` | Score a FASTA against a `.mrk` panel (completeness + redundancy) |
| `bench-grid` | Spike-fraction × taxonomic-distance contamination benchmark |

Full option reference: [CLI Reference](https://genomewalker.github.io/genopack/cli/).

## Quality scoring and contamination detection

`genopack check` computes per-genome quality and contamination signals and writes them
back to the archive as a QUAL section, and to a TSV for downstream use.

```bash
# Score all genomes in an archive
genopack check mydb.gpk -o quality.tsv -t 16

# With single-copy marker genes for completeness
genopack check mydb.gpk -o quality.tsv --markers markers_r232.mrk -t 16

# Force rescore even if QUAL section already exists
genopack check mydb.gpk -o quality.tsv --recompute -t 16
```

The GCOV and FCOV covariance sections must be present for contig-level outlier
scoring (CCO, SPE, sibling outlier). Build them with:

```bash
genopack gcov mydb.gpk -t 16
```

`gcov` runs a single shard-scan building GCOV (per-genus), FCOV (per-family),
and FMHR (per-genus FracMinHash references) simultaneously.

### Key quality metrics (TSV columns)

| Column | Description |
|--------|-------------|
| `quality_tier` | `HQ` / `MQ` / `LQ` — **completeness-only** (`comp_eff`); contamination is decoupled and never demotes to LQ. Only the calibrated duplication channel `D` caps HQ→MQ (`run_check.cpp:181-197`) |
| `completeness_effective` | Intrinsic completeness: `marker_completeness`, else `completeness_aamer_core`, else `completeness_post_decontam`. `completeness_cluster_relative` corroborates (geomean) only when intrinsic < 0.50 and sits > 0.30 above it |
| `completeness_cluster_relative` | Fraction of the genus pangenome covered (accessory breadth) — not intrinsic completeness |
| `contamination_leakage` | Minimizer mass leakage outside expected genus range |
| `contamination_tnf_excess` | TNF Mahalanobis distance from genus centroid |
| `contamination_contig_outlier` | Fraction bp where T² or SPE > 95th percentile (requires GCOV) |
| `contamination_spe` | SPE-based contig outlier fraction (requires GCOV) |
| `contamination_rho_outlier` | Genus-outlier AND family-inlier fraction (requires GCOV+FCOV) |
| `marker_completeness` | Single-copy marker gene completeness (requires `--markers`) |
| `marker_redundancy` | Single-copy marker gene redundancy (requires `--markers`) |

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

genopack::DerepView derep;
derep.open("run.gpd");

// Detect drift relative to the current pack
auto staleness = derep.check(reader);   // Valid / LayoutChangedSameLiveSet / Stale… / Mismatch

auto status = derep.status_for_accession("GCA_000008085.1");
if (status.kind == genopack::RepStatus::Kind::Member) {
    std::vector<float> emb(derep.stats().embedding_dim);
    derep.embedding_for_rep(status.rep_id, emb);
}
```

## License

MIT
