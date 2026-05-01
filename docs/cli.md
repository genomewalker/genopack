# CLI Reference

## Synopsis

```
genopack <command> [options]
```

All commands accept `--help` / `-h` for option details. Wherever a flag below takes an *archive*, the path may be either a single archive directory (containing `toc.bin`) or a directory containing one or more `part_*.gpk` (multipart set).

---

## `build`

Build a new archive from a TSV.

```bash
genopack build -i genomes.tsv -o mydb.gpk [options]
```

The output `mydb.gpk` is a directory containing `toc.bin` plus section files. Defaults: per-taxon shard grouping (`--taxon-group`, genus rank), kmer-NN sort within shards, OPH sketches (`--sketch`, k=16, sketch size 10 000), CIDX contig index, auto codec.

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | required | Input TSV (`accession\tfile_path\t[completeness\tcontamination\ttaxonomy\t…]`) |
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
| `--sketch-kmers` | unset | Comma list (e.g. `16,21,31`) → multi-k SKCH v2 in a single pass |
| `--sketch-size` | 10000 | Number of OPH bins |
| `--sketch-syncmer` | 0 | Open syncmer prefilter `s` (0 disables) |
| `--coordinator` | unset | NFS manifest coordinator: `manifest_dir:/output.gpk` |
| `-v / --verbose` | off | Verbose progress |

---

## `merge`

Merge multiple `.gpk` archives into one. Uses parallel `pwrite` (one thread per part) for NFS efficiency.

```bash
genopack merge -l parts.txt -o merged.gpk
# or
genopack merge part1.gpk part2.gpk part3.gpk -o merged.gpk
```

| Flag | Default | Description |
|------|---------|-------------|
| `-l / --list` | | Text file with one `.gpk` path per line |
| `-o / --output` | required | Output path |
| `-t / --threads` | auto | Merge threads (one per input part) |

---

## `stat`

Print archive statistics.

```bash
genopack stat mydb.gpk [--json]
```

Output: generation, shard count, live/total genome count, total bp, compression ratio, per-section inventory. Accepts a single archive directory or a multipart set; multipart sets show an aggregated total plus per-part breakdown.

---

## `inspect`

Report SKCH layout and preload memory cost.

```bash
genopack inspect mydb.gpk [--json]
```

For each archive (single or each `part_*.gpk` in a multipart directory) prints: live genome count, `sketch_size` (bins), `mask_words` (`ceil(sketch_size/64)`), the list of `kmer_sizes` stored, bytes per sketch per k, bytes per genome, and total preload size. Use this to decide whether to mmap-preload sketches or stream them frame-by-frame on memory-tight nodes. `--json` emits machine-readable output.

---

## `extract`

Extract genomes as FASTA.

```bash
genopack extract mydb.gpk [filters] -o out.fasta
```

| Flag | Description |
|------|-------------|
| `--accession ACC` | Extract single genome |
| `--accessions-file FILE` | Extract list of accessions (one per line) |
| `--min-completeness FLOAT` | Completeness filter (0–100) |
| `--max-contamination FLOAT` | Contamination filter |
| `--min-length INT` | Minimum assembly length in bp |
| `--max-contigs INT` | Maximum contig count |
| `-o / --output` | Output FASTA (default: stdout) |

---

## `slice`

Extract a subsequence by accession and coordinates.

```bash
genopack slice mydb.gpk --accession GCA_000008085.1 --start 100000 --length 5000 -o region.fasta
```

Decompresses only the checkpoint blocks covering the requested region (sub-genome granularity).

---

## `add`

Append genomes to an existing archive (new shard generation).

```bash
genopack add mydb.gpk -i new_genomes.tsv [-t 16]
```

Existing shards are untouched. The catalog receives a new CATL fragment. Use `repack` afterwards if taxonomy grouping is required.

---

## `rm`

Soft-delete (tombstone) genomes.

```bash
genopack rm mydb.gpk GCA_000001405 GCA_000002655
```

Marks genomes as deleted in a new catalog fragment. Physical space is not reclaimed; use `repack` to compact.

---

## `dedup`

Tombstone duplicate genomes (same sequence, different accession) in place.

```bash
genopack dedup mydb.gpk [--dry-run]
```

Walks every shard, hashes each genome's canonical FASTA content, groups duplicates and tombstones all but one representative per group. Modifies the archive in place by appending a new CATL fragment with the tombstones; physical bytes are reclaimed only by `repack`. With `--dry-run`, prints the duplicate groups without writing.

---

## `repack`

Re-shard an archive by taxonomy for fast per-taxon NFS access.

```bash
genopack repack input.gpk output.gpk [options]
```

| Flag | Default | Description |
|------|---------|-------------|
| `-t / --threads` | 1 | OMP decompression threads |
| `-z / --zstd-level` | 6 | Output compression level |
| `-m / --max-memory` | 32 | Max buffered FASTA data in GB before eviction |
| `--taxonomy-rank` | `g` | `g` = genus, `f` = family |
| `-v / --verbose` | off | Log every source shard processed |

**When to use:** After building a large archive without `--taxonomy-group`, or after many `add` operations that scattered genomes across shards. A repacked archive allows geodesic-style tools to read only the ~1,900 shards belonging to Salmonella instead of all 24,000+ shards.

**Algorithm:** Three-phase two-pass design.

1. **Phase 1 (fast):** Reads only `GenomeDirEntry` arrays from each shard (~300 MB total for a 3.1 TB archive). Builds a full genome→taxonomy routing table in memory.
2. **Phase 2:** Sorts all genome records by `(taxonomy, oph_fingerprint)`.
3. **Phase 3:** Single sequential pass through the source archive; decompresses FASTAs with OMP parallelism; routes each genome to its pre-determined taxon writer; flushes the *largest* writer when memory cap is hit (minimises partial-shard fragmentation).

---

## `taxonomy`

Group of taxonomy utilities. Each operation is its own subcommand.

### `taxonomy show`

```bash
genopack taxonomy show mydb.gpk [--accession ACC] [--json]
```

Print the lineage for one accession (or every accession when `--accession` is omitted). With `--json`, emits one JSON object per line.

### `taxonomy normalize`

```bash
genopack taxonomy normalize -i raw.tsv -o normalized.tsv [--ncbi-taxdump DIR]
```

Take an `accession\ttaxonomy\tfile_path` TSV and produce a normalized 10-rank lineage TSV. With `--ncbi-taxdump`, fills missing ranks from NCBI `nodes.dmp` + `names.dmp`.

### `taxonomy partition`

```bash
genopack taxonomy partition -i normalized.tsv -n 8 -o parts/ [-r g]
```

Partition a normalized TSV into `N` balanced parts at a given rank (`g` = genus, `f` = family) for parallel/distributed builds. Writes `part_0.tsv` … `part_{N-1}.tsv` under the output directory.

### `taxonomy assign-taxids`

```bash
genopack taxonomy assign-taxids -i normalized.tsv -o registry.tsv [--acc-map acc2cid.tsv]
```

Assign canonical concept IDs to lineage paths and emit `canonical_path\tconcept_id` (sorted by path). Optionally writes a per-accession `accession\tconcept_id\ttaxonomy` map.

### `taxonomy diff`

```bash
genopack taxonomy diff --current current.tsv --gtdb bac120_taxonomy.tsv --gtdb ar53_taxonomy.tsv -o out/
```

Diff a current taxonomy TSV against a new GTDB release and write per-category TSVs (added, removed, reassigned) plus a `summary.txt` to `out/`. With `--write-unchanged`, also writes the (often huge) `unchanged.tsv`.

### `taxonomy patch`

Patch taxonomy assignments in place, either against a `.gpk` archive or a flat input TSV.

```bash
# Patch the archive directly
genopack taxonomy patch --archive mydb.gpk --patch reassignments.tsv

# Patch an input TSV before rebuilding
genopack taxonomy patch --tsv genomes.tsv --patch reassignments.tsv [--tsv-out patched.tsv]

# GTDB-Tk classify_summary input
genopack taxonomy patch --archive mydb.gpk --patch gtdbtk.summary.tsv --gtdbtk
```

Default patch format is `accession\tnew_taxonomy`. `--gtdbtk` accepts GTDB-Tk's `classify_summary` format directly. `--no-normalize` disables 7→10 rank normalization (default: on).

---

## `taxdump`

Export taxonomy in NCBI or columnar binary format.

```bash
genopack taxdump mydb.gpk -f taxdump -o ./taxdump/
genopack taxdump mydb.gpk -f columnar -o ./taxonomy/
```

| Format | Output files | Description |
|--------|-------------|-------------|
| `taxdump` | `names.dmp`, `nodes.dmp`, `acc2taxid.dmp` | NCBI taxdump - Kraken/Kaiju compatible |
| `columnar` | `acc2taxid.bin`, `taxnodes.bin`, `acc2taxid.tsv`, `taxonomy.tsv` | Fast offline lookup |

---

## Similarity search and contig lookup (library-only)

There is no `genopack similar` or `genopack cidx` CLI. KMRX/HNSW similarity search and CIDX contig→genome lookup are exposed through the C++ API only:

- `ArchiveReader::find_similar(...)` / `find_similar_by_accession(...)` — KMRX cosine similarity, HNSW-accelerated when an HNSW section is present, linear-scan fallback otherwise.
- `ArchiveReader::find_contig_genome_id(accession)` — CIDX binary search, ~150M queries/s/core.

See [API → Similarity & contig lookup](api.md#similarity--contig-lookup).

---

## `reindex`

Build or rebuild auxiliary index sections in place.

```bash
genopack reindex mydb.gpk [options]
```

| Flag | Default | Description |
|------|---------|-------------|
| `--force` | off | Rebuild indexes even if already present |
| `--no-gidx` | off | Skip GIDX (useful when only `--skch` is needed and GIDX is absent/unwanted) |
| `--txdb` | off | Build the taxonomy tree (TXDB) from existing TAXN lineage strings |
| `--cidx FILE` | unset | Build the contig accession index (CIDX) from a build TSV (`accession\ttaxonomy\tfile_path`) |
| `--cidx-threads` | 8 | Threads for parallel FASTA decompression while indexing contigs |
| `--skch` | off | Compute OPH sketches for genomes missing from existing SKCH sections |
| `--skch-threads` | 8 | Threads for parallel sketch computation |
| `--sketch-kmer` | inherit / 16 | OPH k-mer size for a single-k SKCH section |
| `--sketch-kmers` | unset | Comma list (e.g. `16,21,31`) → multi-k SKCH v2 in one pass |
| `--sketch-size` | inherit / 10000 | OPH sketch size |
| `--sketch-syncmer` | inherit / 0 | Open-syncmer prefilter `s` (0 disables) |

Typical uses: an archive built with `--no-cidx` later wants contig lookup (`--cidx genomes.tsv`); a TAXN-only archive needs the tree (`--txdb`); SKCH layout needs to be upgraded to V4 seekable (`--skch --force`); or a multi-k variant is needed (`--skch --sketch-kmers 16,21,31 --force`).

---

## `coordinator`

NFS-coordinated assembly mode for distributed builds. Workers run `genopack build` with `--coordinator <manifest_dir>:<output.gpk>`; the coordinator process waits for the expected number of worker manifests, then merges parts into a single archive.

```bash
genopack coordinator -o mydb.gpk --workers 64 --nfs-dir /shared/manifests/ \
    [--ntdb /path/to/ncbi/taxdump/]
```

| Flag | Default | Description |
|------|---------|-------------|
| `-o / --output` | required | Final merged archive path |
| `--workers` | required | Expected number of worker manifests |
| `--nfs-dir` | required | Shared directory where workers drop manifests |
| `--ntdb` | unset | NCBI `nodes.dmp` + `names.dmp` directory; embeds an NTDB section for offline taxid resolution |
