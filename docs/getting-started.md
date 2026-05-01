# Getting Started

## Input TSV format

`genopack build` and `genopack add` accept a tab-separated file:

| Column | Required | Description |
|--------|----------|-------------|
| `accession` | ✓ | Unique genome identifier (e.g. `GCA_000008085.1`) |
| `file_path` | ✓ | Path to `.fa`, `.fna`, `.fa.gz`, or `.fna.gz` |
| `completeness` | | CheckM completeness % (0–100) |
| `contamination` | | CheckM contamination % |
| `taxonomy` | | Full lineage string (`d__;p__;c__;o__;f__;g__;s__`) |
| extra columns | | Stored verbatim in `meta.tsv` sidecar |

Example:

```tsv
accession	file_path	completeness	contamination	taxonomy
GCA_000008085.1	/data/GCA_000008085.1.fna.gz	98.7	0.3	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;s__Salmonella enterica
GCA_000006945.2	/data/GCA_000006945.2.fna.gz	99.1	0.1	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
```

---

## Build an archive

```bash
genopack build -i genomes.tsv -o mydb.gpk -t 24 -z 6
```

The output `mydb.gpk` is a directory containing `toc.bin` and section files. Defaults: per-taxon shard grouping (genus rank), kmer-NN sort within each shard, OPH sketches (k=16, sketch size 10 000), CIDX contig index, auto codec.

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | required | Input TSV |
| `-o / --output` | required | Output archive directory (`.gpk`) |
| `-t / --threads` | 16 | I/O threads (decompression + compression) |
| `-z / --zstd-level` | 6 | zstd compression level (1–22) |
| `-p / --parallel` | 1 | Parallel build workers (auto-merge) |
| `--no-dict` | off | Disable shared dictionary training |
| `--ref-dict` | off | Use first genome in each shard as reference content dictionary |
| `--delta` | off | Compress non-reference blobs against first genome via zstd prefix |
| `--mem-delta` | off | k=31 k-mer seeded exact-match encoding for highly similar shard groups |
| `--2bit` | off | Pack nucleotides to 2 bits/base before zstd (~1.5–2× extra compression) |
| `--no-cidx` | off | Skip CIDX contig index (recommended for >1M genomes) |
| `--kmer-sort` / `--no-kmer-sort` | on | Sort genomes within each shard by kmer4 NN chain |
| `--taxon-group` / `--no-taxon-group` | on | Group genomes into per-taxon shards (requires taxonomy column) |
| `--taxon-rank` | `g` | Rank for grouping (`g` = genus, `f` = family) |
| `--sketch` / `--no-sketch` | on | Compute OPH sketches |
| `--sketch-kmer` | 16 | OPH sketch k-mer size |
| `--sketch-kmers` | unset | Comma list (e.g. `16,21,31`) → multi-k SKCH v2 in one pass |
| `--sketch-size` | 10000 | Number of OPH bins |
| `--sketch-syncmer` | 0 | Open syncmer prefilter `s` (0 disables) |
| `--coordinator` | unset | NFS manifest coordinator: `manifest_dir:/output.gpk` |
| `-v / --verbose` | off | Verbose progress |

---

## Archive statistics

```bash
genopack stat mydb.gpk             # text
genopack stat mydb.gpk --json      # JSON
```

`stat` works on a single archive directory or on a directory containing `part_*.gpk` (multipart set). Output: generation, shard count, genome count (total + live), total bp, compression ratio.

## Inspect SKCH layout

```bash
genopack inspect mydb.gpk
genopack inspect parts_dir/        # iterate all part_*.gpk
genopack inspect mydb.gpk --json
```

Reports per-archive sketch size, mask words, available k-mer sizes, bytes per genome, and total preload-cost estimate. Useful when planning memory budgets for downstream consumers (e.g. geodesic).

---

## Extract genomes

```bash
# Single accession
genopack extract mydb.gpk --accession GCA_000008085.1 -o genome.fasta

# From file (one accession per line)
genopack extract mydb.gpk --accessions-file list.txt -o genomes.fasta

# With quality filter
genopack extract mydb.gpk --min-completeness 95 --max-contamination 5 -o filtered.fasta
```

---

## Taxonomy

`genopack taxonomy` is a subcommand group:

```bash
# Look up one genome
genopack taxonomy show mydb.gpk --accession GCA_000008085.1
genopack taxonomy show mydb.gpk --json

# Normalize an input TSV against an NCBI taxdump
genopack taxonomy normalize -i raw.tsv -o normalized.tsv --ncbi-taxdump taxdump/

# Partition a normalized TSV into N balanced parts for distributed build
genopack taxonomy partition -i normalized.tsv -n 8 -o parts/ -r f

# Assign stable taxids to a normalized TSV
genopack taxonomy assign-taxids -i normalized.tsv -o registry.tsv

# Diff against a new GTDB release
genopack taxonomy diff --current normalized.tsv --gtdb new_*.tsv -o diff/

# Apply a curated taxonomy patch
genopack taxonomy patch --patch patch.tsv --archive mydb.gpk
```

Export the archive's taxonomy as NCBI taxdump or columnar binary:

```bash
genopack taxdump mydb.gpk -f taxdump -o ./taxdump/      # Kraken/Kaiju compatible
genopack taxdump mydb.gpk -f columnar -o ./taxonomy/    # fast offline lookup
```

The columnar binary export produces four files:

| File | Description |
|------|-------------|
| `acc2taxid.bin` | Sorted `(FNV-1a-64(acc), taxid)` pairs — O(log n) binary search |
| `taxnodes.bin` | Sorted node records + name pool — O(log n) taxid lookup |
| `acc2taxid.tsv` | `accession\ttaxid` — pandas/polars/R compatible |
| `taxonomy.tsv` | `taxid\tparent_taxid\trank\tname\tis_synthetic` |

---

## Similarity (library only)

There is no `genopack similar` CLI. KMRX k=4 cosine similarity and the optional HNSW ANN index are exposed through the C++ API:

```cpp
genopack::ArchiveReader reader; reader.open("mydb.gpk");
const float* p = reader.kmer_profile_by_accession("GCA_000008085.1");
auto hits = reader.find_similar_by_accession("GCA_000008085.1", 20);
```

`find_similar` uses the HNSW section if present; otherwise it falls back to a linear KMRX scan. HNSW is not built by default — use the C++ `ArchiveBuilder` API to opt in.

---

## Contig accession lookup (library only)

There is no `genopack cidx` CLI. CIDX is exposed via `ArchiveReader::find_contig_genome_id` / `batch_find_contig_genome_ids`:

```cpp
uint32_t gid = reader.find_contig_genome_id("NZ_JAVJIU010000001.1");
```

---

## Repack by taxonomy

After building a large archive with randomly distributed shards, re-shard by genus for fast per-taxon NFS access. Geodesic-style workflows that process one taxon at a time read only the relevant shards instead of the entire archive.

```bash
genopack repack mydb.gpk mydb_taxon.gpk -t 24 -z 6 -m 32
```

| Flag | Default | Description |
|------|---------|-------------|
| `-t / --threads` | 1 | OMP decompression threads |
| `-z / --zstd-level` | 6 | Output compression level |
| `-m / --max-memory` | 32 | Max buffered FASTA GB before eviction |
| `--taxonomy-rank` | `g` | `g` = genus, `f` = family |

The algorithm runs in three phases:

1. **Directory scan** - reads only shard headers + `GenomeDirEntry` arrays (~300 MB vs 3.1 TB for full blobs), builds the complete genome→taxonomy routing table.
2. **Sort** - sorts all genome records by `(taxonomy, oph_fingerprint)` in memory.
3. **Decompress + route + write** - single sequential pass; flushes the *largest* active shard writer when memory cap is hit (minimises shard fragmentation vs "flush all").

---

## Append genomes

```bash
genopack add mydb.gpk -i new_genomes.tsv
```

Appends a new shard generation. Existing shards are untouched; the catalog is updated with a new CATL fragment. Use `genopack repack` afterwards if taxonomy grouping is important.

---

## Remove genomes

```bash
genopack rm mydb.gpk GCA_000001405 GCA_000002655
```

Soft-deletes (tombstones) genomes. Physical space is not reclaimed; use `genopack repack` to compact.

---

## Distributed build

For collections too large for a single node:

## Distributed build

For collections too large for a single node, use the NFS manifest coordinator:

```bash
# Coordinator: allocates write offsets, assembles final TOC
genopack coordinator -o /nfs/output.gpk --nfs-dir /nfs/manifest/ --workers 4 \
    --ntdb /path/to/ncbi_taxdump/

# Workers: build to local scratch, transfer sections via NFS manifest
genopack build -i part_N.tsv -o /scratch/part_N.gpk -t 24 -z 6 \
    --coordinator /nfs/manifest/:/nfs/output.gpk
```

`--ntdb` makes the coordinator embed the NCBI tree as an NTDB section in the final archive.

Alternatively, build parts independently and merge:

```bash
for i in 0 1 2 3; do
    genopack build -i part_$i.tsv -o parts/part_$i.gpk -t 24 -z 6
done
genopack merge -l <(ls parts/*.gpk) -o merged.gpk
```

Or skip the merge and keep the parts as a multipart set — `ArchiveSetReader` (and every CLI command that takes an archive path) accepts a directory containing `part_*.gpk`.

---

## Deduplication

```bash
genopack dedup mydb.gpk              # tombstone duplicates in place
genopack dedup mydb.gpk --dry-run    # report duplicates without modifying the archive
```

Detects exact-sequence duplicates (same content under different accessions). The genome with the highest completeness wins; the rest are tombstoned in a new catalog fragment. Reclaim space afterwards with `genopack repack`.

---

## Reindex

Append or rebuild auxiliary sections on an existing archive:

```bash
# Add OPH sketches (single k or multi-k)
genopack reindex mydb.gpk --skch --sketch-kmer 16 --skch-threads 16
genopack reindex mydb.gpk --skch --sketch-kmers 16,21,31 --skch-threads 16

# Build TXDB from existing TAXN strings
genopack reindex mydb.gpk --txdb

# Build CIDX from the original build TSV
genopack reindex mydb.gpk --cidx genomes.tsv --cidx-threads 16

# Force rebuild
genopack reindex mydb.gpk --skch --force

# Skip GIDX (when only --skch is needed and GIDX is absent/unwanted)
genopack reindex mydb.gpk --skch --no-gidx
```
