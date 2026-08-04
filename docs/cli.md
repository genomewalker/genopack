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

The output `mydb.gpk` is a directory containing `toc.bin` plus section files. Defaults: per-taxon shard grouping (`--taxon-group`, genus rank), kmer-NN sort within shards, OPH sketches (`--sketch`, k=16, sketch size 10 000), auto codec. The CIDX contig index is opt-in (`--cidx`).

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
| `--mem-delta` / `--no-mem-delta` | on | k=31 k-mer seeded exact-match encoding for highly similar shard groups |
| `--2bit` / `--no-2bit` | on | Pack nucleotides to 2 bits/base before zstd (~1.5–2× extra compression) |
| `--cidx` / `--no-cidx` | off (skip) | Build CIDX contig index (off by default; pass `--cidx` to build it) |
| `--kmer-sort` / `--no-kmer-sort` | on | Sort genomes within each shard by kmer4 NN chain |
| `--taxon-group` / `--no-taxon-group` | on | Group genomes into per-taxon shards (requires taxonomy column) |
| `--taxon-rank` | `g` | Rank for grouping (`g` = genus, `f` = family) |
| `--sketch` / `--no-sketch` | on | Compute OPH sketches |
| `--sketch-kmer` | 16 | OPH sketch k-mer size |
| `--sketch-kmers` | `16,21,31` | Comma list (e.g. `16,21,31`) → multi-k SKCH in a single pass |
| `--sketch-size` | 10000 | Number of OPH bins |
| `--sketch-syncmer` | `-1` (auto: `s=k/3`) | Open syncmer prefilter `s` (0 disables) |
| `--gstx` / `--no-gstx` | on | Build GSTX genus-stats index (`--no-gstx` to disable) |
| `--pcore` / `--no-pcore` | on | Build the dense PCORE per-genus small-contig reference (needed for small-contig contamination; dominant build cost and on-disk size). `--no-pcore` to disable |
| `--pcore-frac` | 100 | FMH subsampling factor `N` for PCORE aamer extraction (keep 1/N aamers). Lower = denser; 1 = keep all. Overridden by `GENOPACK_AAMER_FRAC` or a `--markers` panel |
| `--tier` | off | Emit a `.ptier` side-channel file for later `genopack tier merge` (records every `(aamer_hash, genus_hash)` pair for a global IDF tier table). Implies `--pcore` |
| `--micro-genus-threshold` | 0 (auto) | Min genomes for a genus to get its own shard + GSTX/GCOV/FMHR consensus model; smaller genera are bin-packed and unmodeled. Auto: 4 (≤50k genomes), 8 (≤500k), 16 (≤2M), else 32 |
| `--markers` | auto | Path to `markers.mrk` for build-time marker completeness scoring. Auto-resolved from `$GENOPACK_DATA/markers.mrk` or the install share dir |
| `--contam-panel` | auto | Path to contamination `.csp` panel for build-time intra-genome duplication scoring. Auto-resolved from `$GENOPACK_DATA/contamination.csp` or the install share dir; absent → axis skipped |
| `--from-gpk` | unset | Rebuild from an existing `.gpk`: stream decoded sequence from its shards instead of opening per-genome FASTA files. Mutually exclusive with `-i` |
| `--tmpdir` | `/scratch` or `/tmp` | Directory for PCORE/SKCH spill files; use a fast local NVMe path |
| `--thin` | off | Ingest-only preset: write sequences, sketches, TAXN/GIDX/CIDX only; skip compute sections (GSTX/GCOV/CORE/PCORE/GAMI). Do not combine with `--taxon-group` |
| `--coordinator` | unset | NFS manifest coordinator: `manifest_dir:/output.gpk` |
| `-v / --verbose` | off | Verbose progress |

---

## `enrich`

Add compute sections (GSTX/GCOV/CORE/PCORE/GAMI) to a thin archive produced by `build --thin`. Reads sequences locally from the thin `.gpk` (no NFS opens). Pass `--taxon-group` and `--markers` as for a full build.

```bash
genopack enrich -i thin.gpk -o enriched.gpk [options]
```

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | required | Thin source archive (`.gpk`) |
| `-o / --output` | required | Output enriched archive (`.gpk`) |
| `-t / --threads` | 16 | Worker threads |
| `--taxon-group` | off | Group genomes by genus for GSTX/GCOV computation (required for GCOV) |
| `--taxon-rank` | `g` | Taxonomy rank for grouping (`g` = genus, `f` = family) |
| `--markers` | unset | Marker panel (`.mrk`) for CORE/PCORE/GAMI sections |
| `--contam-panel` | unset | Contamination panel (`.csp`) for duplication scoring |
| `--sketch-kmers` | `16,21,31` | Comma-separated k-mer sizes for multi-k SKCH (default: 16,21,31) |
| `--sketch-size` | 10000 | OPH bins |
| `--sketch-syncmer` | auto | Syncmer `s` (default: auto) |
| `--pcore` | off | Also build PCORE/GAMI (use for thin archives built without `--markers`) |
| `--tmpdir` | `/scratch` or `/tmp` | Directory for PCORE/SKCH spill files |

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
| *(no thread option)* | — | No thread flag: merge copies each input part concurrently via `std::async` + `pwrite` (one task per part), so parallelism scales with the number of input archives |

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

## `neighbors`

A priori nearest-in-DB member per query accession from the derep `.gpd` cluster structure (cluster rep = nearest representative-DB member, no alignment). For ancient-metasim absent-taxon placement.

```bash
genopack neighbors derep.gpd --accessions-file queries.txt -o neighbors.tsv
```

| Flag | Default | Description |
|------|---------|-------------|
| `derep` | required | Derep `.gpd` file or directory of `.gpd` parts |
| `--accessions-file` | required | Query accessions, one per line |
| `-o / --output` | stdout | Output TSV |

---

## `extract`

Extract genomes as FASTA.

```bash
genopack extract mydb.gpk [filters] -o out.fasta
```

| Flag | Description |
|------|-------------|
| `--accession ACC` | Extract genome(s) by accession (repeatable: `--accession A --accession B`) |
| `--accessions-file FILE` | Extract list of accessions (one per line) |
| `--min-completeness FLOAT` | Completeness filter (0–100) |
| `--max-contamination FLOAT` | Contamination filter |
| `-o / --out` | Output FASTA (default: stdout) |

---

## `slice`

Extract a subsequence by accession and coordinates.

```bash
genopack slice mydb.gpk --accession GCA_000008085.1 --start 100000 --length 5000 --fasta
```

Decompresses only the checkpoint blocks covering the requested region (sub-genome granularity).

---

## `add`

Append genomes to an existing archive (new shard generation).

```bash
genopack add mydb.gpk -i new_genomes.tsv
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

## `subset`

Subset an archive to a set of accessions by **copying** the kept genomes' sections (sequences, sketches, QUAL, taxonomy, KMRX) directly — no source FASTA, no re-sketching, original QUAL preserved. For a multipart directory, subsets each `part_*.gpk`. Genome IDs are preserved (sparse); corpus-derived sections (GSTX/CIDX/HNSW) are dropped.

```bash
genopack subset input output --accessions-file accs.txt [options]
```

| Flag | Default | Description |
|------|---------|-------------|
| `input` | (required) | Source `.gpk` archive or multipart directory |
| `output` | (required) | Output `.gpk` (or directory, for a multipart source) |
| `--accessions-file` | (required) | File of accessions to keep: one per line, or a TSV with an `accession`/`name` column |
| `-t, --threads` | `1` | Threads for shard decompression |

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

## `check`

Compute per-genome quality scores (completeness, contamination) and write the QUAL section. See [Quality Scoring](quality.md) for the underlying formulas.

```bash
genopack check mydb.gpk -o quality.tsv [options]
```

| Flag | Default | Description |
|------|---------|-------------|
| `pack` | required | Path to `.gpk` archive or directory of parts |
| `-g / --genomes` | all in pack | Optional accession list (one per line) |
| `-o / --output` | | TSV output path for quality table |
| `-t / --threads` | 8 | Threads |
| `--min-genus-size` | 3 | Min genus members for saturated tier |
| `--leakage-threshold` | 0.05 | Containment leakage threshold |
| `--recompute` | off | Ignore existing QUAL section and force full rescan |
| `--scan-all` | off | Force every genome through FASTA-level analysis (intrinsic completeness for all, not just flagged); implied when `--markers` is given |
| `--markers` | unset | Path to markers `.mrk` DB; enables marker-based completeness/redundancy scoring |
| `--cross-genus-margin` | 0.0 | log-LR margin a foreign genus must beat the host by for `cross_genus` to flag a contig; 0 flags on any foreign LL > host LL |
| `--dup-restore` | unset | Re-inject the build-time `core_dup` axis from a quality TSV or `cladesplit score` TSV, when a prior run overwrote QUAL without it |

---

## `ingest`

Ingest external quality (CheckM2/anvi'o) into the archive as provenance-carrying XQAL columns.

```bash
genopack ingest mydb.gpk --checkm2 quality_report.tsv
```

| Flag | Default | Description |
|------|---------|-------------|
| `pack` | required | Path to `.gpk` archive or directory of parts |
| `--checkm2` | unset | CheckM2 `quality_report.tsv` (Name/Completeness/Contamination) |
| `--anvio` | unset | anvi'o completeness TSV (bin name/% completion/% redundancy) |

---

## `report`

Emit a unified per-genome quality table resolved through a named profile. Each axis is sourced from exactly one column (built-in rule or stored policy) and the report carries that column's provenance (tool/method) — fusion is explicit, never silent.

```bash
genopack report mydb.gpk -p best -o report.tsv
```

| Flag | Default | Description |
|------|---------|-------------|
| `pack` | required | Path to `.gpk` archive or directory of parts |
| `-p / --profile` | `best` | Profile name: built-in `{intrinsic,external,calibrated,best}` or a stored profile |
| `-o / --output` | stdout | TSV output path |
| `--list` | off | List built-in + stored profiles and the available provenance columns, then exit |

---

## `profile`

Manage named reporting/fusion profiles stored in the archive (SEC_PROF). Each operation is its own subcommand.

### `profile add`

```bash
genopack profile add mydb.gpk --name mypolicy -s completeness=checkm2:default -s contamination=intrinsic:default
```

Author a named profile pinning each axis to an exact column identity present in the archive.

| Flag | Default | Description |
|------|---------|-------------|
| `pack` | required | Path to `.gpk` archive or directory of parts |
| `--name` | required | Profile name (must not be a built-in name) |
| `--description` | unset | Human-facing description (cosmetic; excluded from policy hash) |
| `-s / --select` | required | Axis selector `axis=tool:method[@cal][/tool:method[@cal]]` (repeatable); the optional `/…` is a fallback |

### `profile list`

```bash
genopack profile list mydb.gpk
```

List stored profiles and available columns.

---

## `calibrate`

Calibrate intrinsic completeness (aamer genus-core) against ground truth and write a CQAL section of calibrated per-genome columns. Ground truth: `--checkm2` TSV, or the ingested XQAL if omitted. Also writes an isotonic+OLS JSON model and prints RMSE.

```bash
genopack calibrate mydb.gpk --checkm2 quality_report.tsv -o calibration.json
```

| Flag | Default | Description |
|------|---------|-------------|
| `archive` | required | Path to `.gpk` archive or directory of parts |
| `--checkm2` | unset | CheckM2 TSV ground truth (falls back to ingested XQAL CheckM2 columns) |
| `-o / --output` | `calibration.json` | Output JSON path for calibration model |
| `-t / --threads` | 8 | Threads |

---

## `qcontig`

Dump the per-contig quality overlay (SEC_QCONTIG): one row per (genome, contig) with offset/length/TNF/leakage, so you can see which contigs drive a genome's contamination.

```bash
genopack qcontig mydb.gpk --min-foreign 0.30 -o flagged_contigs.tsv
```

| Flag | Default | Description |
|------|---------|-------------|
| `pack` | required | Path to `.gpk` archive or directory of parts |
| `-a / --accession` | all | Restrict to one genome accession |
| `-o / --output` | stdout | TSV output path |
| `--min-outlier` | 0.0 | TNF channel (LONG contigs): emit contigs whose GCOV T² or SPE percentile ≥ this (e.g. 0.99) |
| `--min-foreign` | 0.0 | Protein-aamer channel (SMALL contigs): emit contigs whose `foreign_fraction` ≥ this (e.g. 0.30) |
| `--min-lr` | 0.0 | Protein-aamer channel: also require foreign log-LR ≥ this; pair with `--min-foreign` (e.g. 3.0) |
| `--min-leakage` | 0.0 | Also flag contigs with containment-leakage score ≥ this |

---

## `decontaminate`

Iteratively remove contaminated genomes (per-contig CCO above a threshold), rebuilding genus/family models from the survivors each round so the DB and its consensus stay clean.

```bash
genopack decontaminate mydb.gpk [--dry-run]
```

| Flag | Default | Description |
|------|---------|-------------|
| `archive` | required | Path to `.gpk` archive |
| `--max-fmh-z` | 5.0 | Remove genomes whose per-genus FMH-minority z (robust median+MAD delta from genus baseline) exceeds this — the genus-band signal |
| `--min-delta` | 0.02 | Also require `fmh_minority` to exceed the genus baseline by this absolute amount (guards tight-baseline genera from z blow-up) |
| `--max-cco` | 1.5 | Remove genomes whose CCO contamination % (per-contig T²∪SPE TNF-outlier bp vs the genus GCOV covariance) exceeds this — the distant-band signal (family/order/phylum); needs contigs ≥20kb |
| `--max-iters` | 5 | Max decontamination rounds |
| `-t / --threads` | 8 | Threads |
| `--dry-run` | off | Report what would be removed (with fmh_z/cco per axis) without modifying the archive |

---

## `gcov`

Build (or rebuild) the GCOV per-genus covariance section in an existing `.gpk` archive.

```bash
genopack gcov mydb.gpk
```

| Flag | Default | Description |
|------|---------|-------------|
| `archive` | required | Path to `.gpk` archive |
| `-t / --threads` | 8 | Parallel shard readers |

---

## `pcore`

Build SEC_PCORE: the unified dense per-genus aamer reference (every aamer + prevalence). Dense enough to detect SMALL foreign contigs (the conserved-only CORE is too sparse for 1-2kb).

```bash
genopack pcore mydb.gpk
```

| Flag | Default | Description |
|------|---------|-------------|
| `archive` | required | Path to `.gpk` archive or directory of parts |
| `-t / --threads` | 8 | Parallel shard readers |
| `--members` | all live | Reference accession list (one per line); only these genomes build the reference |

---

## `gami`

SEC_GAMI tools: precomputed global aamer multiplicity index.

### `gami build`

```bash
genopack gami build mydb.gpk
```

Build and append SEC_GAMI v2 to a `.gpk` archive. Decodes PCORE/CORE once, serialises the global multiplicity index (exact sorted pairs, zstd), eliminating the 10-30 min GMI rebuild at each `genopack check` run.

| Flag | Default | Description |
|------|---------|-------------|
| `archive` | required | Path to `.gpk` or directory of parts |
| `-t / --threads` | 8 | Parallel decode threads |

---

## `tier`

IDF tier table tools: merge per-part `.ptier` files into a global `.gtier`, then stamp tier bytes into PCORE sections producing PCORE v2.

### `tier merge`

```bash
genopack tier merge -i part1.ptier -i part2.ptier -o global.gtier
```

Merge `.ptier` side-channel files from multi-part builds into a global `.gtier` IDF table.

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | required | Input `.ptier` files (one per build part; repeatable) |
| `-o / --output` | required | Output global `.gtier` table |

### `tier stamp`

```bash
genopack tier stamp -i mydb.gpk --table global.gtier -o mydb.tier.gpk
```

Stamp per-aamer tier bytes from a global `.gtier` into a PCORE section, rewriting PCORE v1 → v2 in-place.

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | required | Input `.gpk` archive |
| `--table` | required | Global `.gtier` table from `genopack tier merge` |
| `-o / --output` | `<input>.tier.gpk` | Output `.gpk` |

---

## `cladesplit`

Protein-aamer clade-split contamination tier (GUNC-style, fast). Build a panel of genus-diagnostic aamers from clean genomes, then score genomes by lineage-split.

### `cladesplit build`

```bash
genopack cladesplit build -i genomes.tsv -o panel.csp
```

Build a `.csp` panel from clean reference genomes.

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | required | TSV (accession, file_path, taxonomy) |
| `-o / --output` | required | Output `.csp` panel |
| `-t / --threads` | all cores | Worker threads |
| `--frac-c` | 30 | FracMinHash density `1/c` on aamers |
| `--min-aa` | 8 | Min inter-stop AA segment length |
| `--primitive` | `aa` | Protein k-mer primitive: `aamer` \| `metamer` \| `strobemer` |

### `cladesplit score`

```bash
genopack cladesplit score --gpk mydb.gpk --panel panel.csp -o scores.tsv
```

Score genomes against a `.csp` panel (per-genome `minority_fraction`).

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | | TSV (accession, file_path) — or use `--gpk` |
| `--gpk` | | Score live genomes inside a `.gpk` archive directly (flag contamination in a genopack file) |
| `--panel` | required | Panel `.csp` |
| `-o / --output` | required | Output per-genome TSV (incl. `redundancy_fraction`) |
| `-t / --threads` | all cores | Worker threads |
| `--frac-c` | 30 | FracMinHash density `1/c` (must match build) |
| `--min-aa` | 8 | Min AA segment (must match build) |
| `--primitive` | `aa` | Protein k-mer primitive (must match build): `aamer` \| `metamer` \| `strobemer` |

### `cladesplit aamers`

```bash
genopack cladesplit aamers -i genomes.tsv -o dump.bin
```

Dump per-genome sorted-unique aamer sets (binary) for core/completeness R&D.

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | required | TSV (accession, file_path) |
| `-o / --output` | required | Output binary aamer dump |
| `-t / --threads` | all cores | Worker threads |
| `--frac-c` | 30 | FracMinHash density `1/c` |
| `--min-aa` | 8 | Min AA segment |
| `--primitive` | `aa` | Protein k-mer primitive: `aamer` \| `metamer` \| `strobemer` |

---

## `markers`

Build or manage the SCG marker aamer panel (`.mrk` sidecar). Each operation is its own subcommand.

### `markers build`

```bash
genopack markers build --gtdbtk-db /path/to/gtdbtk_db -o markers.mrk
```

Build `markers.mrk` from GTDB-Tk MSA files (bac120 + ar53).

| Flag | Default | Description |
|------|---------|-------------|
| `--gtdbtk-db` | required | GTDB-Tk reference database root |
| `-o / --output` | required | Output `.mrk` file |
| `--threshold` | 0.30 | Default presence threshold |
| `-t / --threads` | 1 | Thread count |
| `--scale` | 1 | FracMinHash scale factor: keep 1/N of hashes (ignored with `--dayhoff6`) |
| `--dayhoff6` | off | Build Dayhoff-6 k=12 syncmer profile pool (recommended; compact, robust to divergence) |
| `--ic-threshold` | 0.25 | Min per-column IC fraction to include k-mer positions (`--dayhoff6` only) |
| `--expected-min-frac` | 0.50 | A marker is counted as expected for a genus only if detectable (≥1 IC-passing syncmer) in at least this fraction of the genus's reference genomes; mirrors CheckM2's single-copy universality criterion |

### `markers remerge`

```bash
genopack markers remerge markers.mrk
```

Append the pre-merged pool to an existing `.mrk` panel (in place) so `MarkerReader` uses its zero-copy mmap fast-path instead of copying ~450MB per worker. No source genomes needed.

| Flag | Default | Description |
|------|---------|-------------|
| `panel` | required | Path to `.mrk` panel |

### `markers score`

```bash
genopack markers score --fasta genome.fna --markers markers.mrk --genus g__Escherichia
```

Score a FASTA file for SCG marker completeness.

| Flag | Default | Description |
|------|---------|-------------|
| `--fasta` | required | Input FASTA file |
| `--markers` | required | Markers `.mrk` file |
| `--genus` | required | Target genus key (e.g. `g__Escherichia`) |
| `--min-hits` | 1 | Min hits per marker to call present |

---

## `bench-grid`

Heterogeneous spike-fraction × ANI-distance benchmark from a manifest TSV.

```bash
genopack bench-grid mydb.gpk --manifest manifest.tsv -o results.tsv
```

| Flag | Default | Description |
|------|---------|-------------|
| `archive` | required | Archive path (`.gpk`) |
| `-o / --output` | required | TSV output path |
| `--manifest` | required | Manifest TSV (`host_genus, contam_genus, ani_label, intra_offset`) |
| `-t / --threads` | 4 | Thread count |
| `-r / --reps` | 5 | Replicates per cell |
| `--seed` | 42 | Random seed |
| `--completeness` | `1.0` | Comma-separated host completeness fractions 0.0-1.0 |

---

## `sim`

Generate synthetic fragmented/contaminated genomes for completeness/contamination benchmarking (CheckM2-compatible 20kb chunk fragmentation).

```bash
genopack sim --ref genome.fna --output-dir sim_out/
```

| Flag | Default | Description |
|------|---------|-------------|
| `--ref` | required | Reference genome FASTA (repeatable) |
| `--taxonomy` | unset | GTDB taxonomy; one per `--ref` in order |
| `--contam` | unset | Contamination source FASTA (repeatable; each adds a grid dimension) |
| `--contam-label` | unset | Label for each `--contam` source (e.g. genus,family,phylum); one per `--contam` |
| `--contam-taxonomy` | unset | GTDB taxonomy for each `--contam` source; one per `--contam` |
| `--contam-self` | off | Add each ref as its own contam source (tests `marker_redundancy`; label=`self`) |
| `--completeness` | `0.1,0.3,0.5,0.7,0.9,1.0` | Comma-separated completeness fractions 0.0-1.0 |
| `--contamination` | `0.0` | Comma-separated contamination fractions 0.0-0.5 |
| `--reps` | 3 | Replicates per combination |
| `--seed` | 42 | Base random seed |
| `--chunk-size` | 20000 | Fixed fragment size in bp (ignored if `--contig-n50` set) |
| `--min-chunk` | 1000 | Min chunk size to keep in bp |
| `--contig-n50` | 0 | N50 for lognormal contig length distribution (0 = fixed `--chunk-size`) |
| `--contig-sigma` | 1.2 | Lognormal sigma for contig lengths |
| `--output-dir` | required | Output directory for FASTA files |
| `--output-tsv` | `<output-dir>/sim_manifest.tsv` | Ground-truth TSV |
| `--manifest-tsv` | `<output-dir>/add_manifest.tsv` | `genopack add` manifest TSV |
| `-t / --threads` | 4 | Parallel workers |

---

## Similarity search and contig lookup (library-only)

There is no `genopack similar` or `genopack cidx` CLI. KMRX similarity search and CIDX contig→genome lookup are exposed through the C++ API only:

- `ArchiveReader::find_similar(...)` / `find_similar_by_accession(...)` — linear KMRX cosine similarity scan.
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
| `--sketch-kmers` | `16,21,31` | Comma list (e.g. `16,21,31`) → multi-k SKCH in one pass |
| `--sketch-size` | inherit / 10000 | OPH sketch size |
| `--sketch-syncmer` | `-1` (auto: `s=k/3`) | Open-syncmer prefilter `s` (0 disables) |

Typical uses: an archive built with `--no-cidx` later wants contig lookup (`--cidx genomes.tsv`); a TAXN-only archive needs the tree (`--txdb`); SKCH layout needs to be upgraded to V4 seekable (`--skch --force`); or a multi-k variant is needed (`--skch --sketch-kmers 16,21,31 --force`).

---

## `verify`

```bash
genopack verify mydb.gpk [--strict] [--coverage-only] [--verbose]
```

Runs two checks. First, **semantic coverage**: cross-checks per-genome/per-genus section row counts against the catalog from the TOC only (no data scan) — SKCH (per k, every k must cover the whole catalog), QUAL, GIDX, ACCX, CATL rows must equal the genome count; GSTX and other compute sections are reported as presence + coverage, not gated. Second, **byte-level checksums**: TailLocator, TocHeader, and per-section XXH128 content checksums. Exits 0 if all checks pass, non-zero on any failure.

| Flag | Default | Description |
|------|---------|-------------|
| `archive` | required | Archive directory to verify |
| `--strict` | off | Fail if any section lacks a content checksum |
| `--coverage-only` | off | Run only the section-coverage cross-check; skip byte-level checksums (seconds) |
| `--verbose` | off | Print OK lines in addition to failures |

---

## `coordinator`

NFS-coordinated assembly mode for distributed builds. Workers run `genopack build` with `--coordinator <manifest_dir>:<output.gpk>`; the coordinator process waits for the expected number of worker manifests, then merges parts into a single archive.

```bash
genopack coordinator -o mydb.gpk --workers 64 --nfs-dir /shared/manifests/
```

| Flag | Default | Description |
|------|---------|-------------|
| `-o / --output` | required | Final merged archive path |
| `--workers` | required | Expected number of worker manifests |
| `--nfs-dir` | required | Shared directory where workers drop manifests |
