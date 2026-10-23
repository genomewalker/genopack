# CLI Reference

> Status as of 2026-10-23: `dedup` describes the sequence-hash fix (2026-10-23); binaries built before it, including 2ac05b64, remove on fingerprint + length alone. Every subcommand and flag in `genopack --help` is documented here (checked by a per-command flag diff against the installed binary). The `--help` strings for `enrich --sketch-kmers`, `build --thin`, `repack` and `reindex --sketch-syncmer` were corrected to match the code on 2026-10-23; binaries built earlier, including 2ac05b64, still print the old text. Open `--help` bugs, where the code wins: the `enrich` description, `enrich --pcore` and `enrich --markers` name CORE/GAMI sections that enrich never writes.

## Synopsis

```
genopack <command> [options]
```

All commands accept `--help` / `-h` for option details. Wherever a flag below takes an *archive*, the path may be either a single `.gpk` file or a directory containing one or more `part_*.gpk` (multipart set).

---

## `build`

Build a new archive from a TSV.

```bash
genopack build -i genomes.tsv -o mydb.gpk [options]
```

The output `mydb.gpk` is a single file (`-o mydb` also writes `mydb.gpk`). Defaults: per-taxon shard grouping (`--taxon-group`, genus rank), kmer-NN sort within shards, OPH sketches (`--sketch`, k=16, sketch size 10 000), auto codec. The CIDX contig index is opt-in (`--cidx`).

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | required unless `--from-gpk` | Input TSV (`accession\tfile_path\t[completeness\tcontamination\ttaxonomy\t…]`); exactly one of `-i` and `--from-gpk` must be given |
| `-o / --output` | required | Output archive file (`.gpk`; appended if missing) |
| `-t / --threads` | 16 | I/O threads (decompression + compression) |
| `-z / --zstd-level` | 6 | zstd compression level (1-22) |
| `-p / --parallel` | 1 | Parallel build workers (auto-merge) |
| `--no-dict` | off | Disable shared dictionary training |
| `--ref-dict` | off | Use first genome in each shard as reference content dictionary |
| `--delta` | off | Compress non-reference blobs against first genome via zstd prefix |
| `--mem-delta` / `--no-mem-delta` | on | k=31 k-mer seeded exact-match encoding for highly similar shard groups |
| `--2bit` / `--no-2bit` | on | Pack nucleotides to 2 bits/base before zstd (~1.5-2× extra compression) |
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
| `--ncbi-taxdump` | unset | Path to an NCBI taxdump dir (`nodes.dmp`+`names.dmp`). Lineage taxa whose scientific name is unambiguous in NCBI keep their native NCBI taxid; taxa NCBI doesn't have get a synthetic concept_id. See [Unified Taxonomy](Taxonomy-Graft.md) |
| `--from-gpk` | unset | Rebuild from an existing `.gpk`: stream decoded sequence from its shards instead of opening per-genome FASTA files. Mutually exclusive with `-i` |
| `--tmpdir` | `/scratch` or `/tmp` | Directory for PCORE/SKCH spill files; use a fast local NVMe path |
| `--thin` | off | Ingest-only preset: skips GSTX and GCOV, which `enrich` builds later; PCORE and GAMI still follow `--pcore`. `--taxon-group` is on by default and `--thin` does not turn it off; pass `--no-taxon-group` so shards flush by size, not genus |
| `--coordinator` | unset | NFS manifest coordinator: `manifest_dir:/output.gpk` |
| `-v / --verbose` | off | Verbose progress |

---

## `enrich`

Rebuild a thin archive produced by `build --thin` into a new archive with GSTX and GCOV, reading sequences from the thin `.gpk` (no NFS opens). PCORE is added only with `--pcore`. Enrich never writes CORE or GAMI: the build path it uses returns before the GAMI step, and no build writes CORE. GSTX, and so GCOV, are accumulated only with taxon grouping and multi-k sketches (`src/builder.cpp:1350`), and both are off by default here, so pass `--taxon-group` and `--sketch-kmers 16,21,31`. The `--help` text still claims CORE/PCORE/GAMI.

```bash
genopack enrich -i thin.gpk -o enriched.gpk [options]
```

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | required | Thin source archive (`.gpk`) |
| `-o / --output` | required | Output enriched archive (`.gpk`) |
| `-t / --threads` | 16 | Worker threads |
| `--taxon-group` | off | Group genomes by genus; required for GSTX and therefore GCOV |
| `--taxon-rank` | `g` | Taxonomy rank for grouping (`g` = genus, `f` = family) |
| `--markers` | unset | Marker panel (`.mrk`) for CORE/PCORE/GAMI sections |
| `--contam-panel` | unset | Contamination panel (`.csp`) for duplication scoring |
| `--sketch-kmers` | unset (single-k, k=16) | Comma-separated k-mer sizes for multi-k SKCH; pass `16,21,31` explicitly for multi-k. The single-k size is fixed at 16; enrich exposes no option for it. |
| `--sketch-size` | 10000 | OPH bins |
| `--sketch-syncmer` | auto | Syncmer `s` (default: auto) |
| `--pcore` | off | Also build PCORE. GAMI is not built despite the `--help` text; run `genopack gami build` afterwards |
| `--tmpdir` | `/scratch` or `/tmp` | Directory for PCORE/SKCH spill files |

---

## `merge`

Merge multiple `.gpk` archives into one. Uses parallel `pwrite`, with one thread per part.

```bash
genopack merge -l parts.txt -o merged.gpk
# or
genopack merge part1.gpk part2.gpk part3.gpk -o merged.gpk
```

| Flag | Default | Description |
|------|---------|-------------|
| `-l / --list` | | Text file with one `.gpk` path per line |
| `-o / --output` | required | Output path |
| `--ncbi-taxdump` | unset | Reconcile the merged taxonomy against an NCBI taxdump dir: shared taxa keep native NCBI taxids |
| `--tmpdir` | `$GENOPACK_SKETCH_SPILL_DIR`, `$GENOPACK_SPILL_DIR`, `/scratch`, `/tmp` | Directory for SKCH spill files during merge |
| *(no thread option)* | - | No thread flag: merge copies each input part concurrently via `std::async` + `pwrite` (one task per part), so parallelism scales with the number of input archives |

---

## `stat`

Print archive statistics.

```bash
genopack stat mydb.gpk [--json]
```

Output: generation, shard count, live/total genome count, total bp, compression ratio, per-section inventory. Accepts a single `.gpk` file or a multipart set; multipart sets show an aggregated total plus per-part breakdown.

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
| `--min-completeness FLOAT` | Completeness filter (0-100) |
| `--max-contamination FLOAT` | Contamination filter |
| `-o / --out` | Output FASTA (default: stdout) |
| `--output-dir DIR` | Write one `{accession}.fa` per genome to this directory |
| `--concat-contigs` | One record per genome: contigs joined by N-runs, header = accession |
| `--concat-ns INT` | N-run length between contigs (default 300; must exceed the longest read) |
| `--nmask` | With `--concat-contigs` and `-o`: also write `<out>.nmask` (per-reference effective length + N-intervals, for kaiku) |
| `--manifest FILE` | Per-genome property TSV (bp, n_contigs, gc, completeness, contamination, oph_fingerprint, taxonomy) for the selected accessions, or the whole catalog when none are given. Metadata-only without `-o`; with `-o` it is written in the same pass as the FASTA |

---

## `slice`

Extract a subsequence by accession and coordinates.

```bash
genopack slice mydb.gpk GCA_000008085.1 --start 100000 --length 5000 --fasta
```

The accession is a positional argument. `--start` (0-based) and `--length` are required; `--fasta` adds a FASTA header.

For archives without 2-bit packing, only the checkpoint blocks covering the region are decompressed. Archives built with the default `build --2bit` packing decompress the whole genome and then slice (`src/shard.cpp`, 2-bit branch of `decompress_sequence_slice`).

---

## `add`

Append genomes to an existing archive (new shard generation).

```bash
genopack add mydb.gpk -i new_genomes.tsv
```

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | required | Input TSV of genomes to append |

Existing shards are untouched. CATL, ACCX, GIDX, TAXN, TXDB, CIDX and KMRX each receive a new fragment. Use `repack` afterwards if taxonomy grouping is required.

---

## `rm`

Soft-delete (tombstone) genomes.

```bash
genopack rm mydb.gpk GCA_000001405 GCA_000002655
```

Records the genome ids in a new tombstone (TOMB) section. Physical space is not reclaimed; use `repack` to compact.

---

## `dedup`

Tombstone genomes whose sequence equals another genome's, keeping one copy.

```bash
genopack dedup mydb.gpk --dry-run --groups groups.tsv   # review first
genopack dedup mydb.gpk --groups removed.tsv
```

| Flag | Default | Description |
|---|---|---|
| `--dry-run` | off | Write the groups TSV and change nothing |
| `--groups` | `dedup_groups.tsv` | Group-to-accession TSV, written on every run |

Genomes already tombstoned are skipped. Candidates are genomes with equal OPH fingerprint and genome length. Each candidate's FASTA is decoded and hashed (XXH3-128 over residues only: header lines and line breaks are ignored, contig boundaries are kept). Only candidates with equal hashes are duplicates. In each duplicate set the genome with the highest completeness is kept (ties: lowest genome_id) and the rest are tombstoned in a new TOMB section; physical bytes are reclaimed only by `repack`. Residue case is hashed as stored, so a soft-masked copy is not a duplicate.

The groups TSV has one row per candidate genome: `candidate_group`, `genome_id`, `accession`, `genome_length`, `oph_fingerprint`, `sequence_hash`, `completeness`, `action`. `action` is `keep`, `remove`, or `distinct` (shares fingerprint and length but no other candidate has its sequence).

---

## `subset`

Subset an archive to a set of accessions by **copying** the kept genomes' sections (sequences, sketches, QUAL, taxonomy, KMRX) directly. This requires no source FASTA or re-sketching and preserves the original QUAL. For a multipart directory, subsets each `part_*.gpk`. Genome IDs are preserved (sparse); corpus-derived sections (GSTX/CIDX/HNSW) are dropped.

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

**When to use:** After building a large archive with `build --no-taxon-group`, or after many `add` operations that scattered genomes across shards. A repacked archive allows geodesic-style tools to read only the ~1,900 shards belonging to Salmonella instead of all 24,000+ shards.

**Algorithm:** Three-phase two-pass design.

1. **Phase 1:** Reads only `GenomeDirEntry` arrays from each shard (~300 MB total for a 3.1 TB archive). Builds a full genome→taxonomy routing table in memory.
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

Take an `accession\ttaxonomy\tfile_path` TSV and produce a normalized 10-rank lineage TSV (`d,l,k,p,c,o,f,g,s,S`). GTDB lineages are normalized by rank propagation. With `--ncbi-taxdump`, eukaryote and virus lineages are resolved against NCBI `nodes.dmp` + `names.dmp`; the dump is downloaded if absent or older than 30 days.

Accessions reused by different source files are resolved by decompressed-sequence content hash: exact repeats are dropped and distinct genomes get a `~dupN` suffix.

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | required | Input TSV (`accession\ttaxonomy\tfile_path`) |
| `-o / --output` | required | Output normalized TSV |
| `--ncbi-taxdump` | unset | NCBI taxdump dir for eukaryote/virus resolution |
| `--dup-manifest` | unset | Write a TSV of resolved accession collisions (`accession, file_path, resolved_accession, action, content_hash`) |

### `taxonomy partition`

```bash
genopack taxonomy partition -i normalized.tsv -n 8 -o parts/ [-r g]
```

Partition a normalized TSV into `N` balanced parts at a given rank (`g` = genus, `f` = family) for parallel/distributed builds, using LPT bin-packing. Every genome of one taxon lands in the same part, sorted by taxonomy. Writes `part_0.tsv` … `part_{N-1}.tsv` under the output directory.

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | required | Normalized TSV (`accession\ttaxonomy\tfile_path`) |
| `-n / --parts` | required | Number of output parts |
| `-o / --output-dir` | required | Output directory |
| `-r / --rank` | `g` | Partition rank: `g` = genus, `f` = family |
| `--contiguous` | off | Fill parts in taxonomy order so sibling taxa land in the same part (lower read fan-out downstream) instead of size-only LPT. Use with `--cap` |
| `--cap` | `0` | Max genomes per taxon group; larger groups are split into contiguous chunks of at most `cap`. `0` = no cap |

### `taxonomy assign-taxids`

```bash
genopack taxonomy assign-taxids -i normalized.tsv -o registry.tsv [--acc-map acc2cid.tsv]
```

Assign canonical concept IDs (`FNV-64(canonical_path) | (1<<63)`, deterministic, no central counter) to lineage paths and emit `canonical_path\tconcept_id` (sorted by path). Optionally writes a per-accession `accession\tconcept_id\ttaxonomy` map.

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | required | Normalized TSV |
| `-o / --registry` | required | Output registry TSV |
| `--acc-map` | unset | Output per-accession TSV |

### `taxonomy diff`

```bash
genopack taxonomy diff --current current.tsv --gtdb bac120_taxonomy.tsv --gtdb ar53_taxonomy.tsv -o out/
```

Diff a current genome set against a new GTDB release. Each genome is classified as `unchanged`, `taxonomy_changed`, `needs_gtdbtk`, `dropped_from_gtdb` or `new_in_release`; one TSV per category plus `summary.txt` go to `-o / --output-dir`. `--current` accepts the original input TSV or the genopack meta TSV. `--gtdb` is repeatable (pass `bac120_taxonomy.tsv` and `ar53_taxonomy.tsv`, or a merged file). With `--write-unchanged`, also writes the (often huge) `unchanged.tsv`.

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

### `taxonomy unify`

Reconcile the taxonomies of several independent `.gpk` packs into one common taxdump without merging their sequence data. The command reads each pack's lineage strings (TAXN), builds one tree, and writes `nodes.dmp`/`names.dmp`/`acc2taxid.dmp` (or columnar binary + TSV) covering every accession. The result provides one taxid space for per-pack mapping references and downstream classification. See [Unified Taxonomy](Taxonomy-Graft.md).

```bash
genopack taxonomy unify prok.gpk organelle.gpk viral.gpk \
    --ncbi-taxdump ./taxdump -o ./common_taxonomy [-f taxdump|columnar]
```

| Flag | Default | Description |
|------|---------|-------------|
| `inputs` | required | Input `.gpk` packs (any number) |
| `-o / --output` | required | Output directory |
| `-f / --format` | `taxdump` | `taxdump` (NCBI `nodes`/`names`/`acc2taxid.dmp`) or `columnar` (binary + TSV) |
| `--ncbi-taxdump` | unset | NCBI taxdump dir: taxa with an unambiguous NCBI name keep their native taxid; GTDB/other-only taxa get a synthetic concept_id |
| `--preserve-taxids` | unset | TSV of `accession → NCBI taxid` (e.g. organelle host taxids). Those accessions keep the given taxid instead of resolving by name. Its NCBI lineage is imported down to a strain leaf. Stale ids remap via `merged.dmp` |
| `--ncbi-backbone` | off | Fold the entire NCBI tree in as backbone nodes, so NCBI taxa without a reference genome still exist in the scaffold |
| `--gtdb-replaces-prok` | off | With `--ncbi-backbone`: excise NCBI's `Bacteria`/`Archaea` descendants (anchors kept) so GTDB stays authoritative for prokaryotes |
| `--lineage-overrides` | unset | `accession\tlineage` TSV replacing those accessions' TAXN lineages at read time (GTDB 7-rank lineages normalized as in `taxonomy patch`). Duplicate or absent accessions are errors |
| `--stable-taxids-from` | unset | A prior unify's columnar `taxonomy.tsv`. Synthetic nodes whose full (rank, name) path exists there keep that taxid; new nodes are numbered above its maximum |

**Always pass `--stable-taxids-from` when reunifying a released taxonomy.** Without it, synthetic ids are numbered in concept-hash order, so any added, removed or renamed node shifts later ids and old taxids can name other taxa. Unify fails if two nodes would share one taxid. See [Unified Taxonomy](Taxonomy-Graft.md#reunifying-a-released-taxonomy).

Output is deterministic (identical inputs → identical output). The NCBI taxdump parse is cached to `<taxdump-dir>/.ncbi_taxdb.cache`.

---

## `taxdump`

Export taxonomy in NCBI or columnar binary format.

```bash
genopack taxdump mydb.gpk -f taxdump -o ./taxdump/
genopack taxdump mydb.gpk -f columnar -o ./taxonomy/
```

| Flag | Default | Description |
|------|---------|-------------|
| `-f / --format` | `columnar` | `taxdump` or `columnar` |
| `-o / --output` | required | Output directory |

| Format | Output files | Description |
|--------|-------------|-------------|
| `taxdump` | `names.dmp`, `nodes.dmp`, `acc2taxid.dmp`, empty `merged.dmp` and `delnodes.dmp` | NCBI taxdump - Kraken/Kaiju compatible |
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
| `--scan-all` | off | Force every genome through FASTA-level analysis (intrinsic completeness for all, rather than only flagged genomes); implied when `--markers` is given |
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

Emit a per-genome quality table resolved through a named profile. Each axis is sourced from one column (built-in rule or stored policy), and the report includes that column's provenance (tool/method). The selected fusion is recorded.

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

Iteratively remove contaminated genomes (per-contig CCO above a threshold), rebuilding genus/family models from the survivors after each round.

```bash
genopack decontaminate mydb.gpk [--dry-run]
```

| Flag | Default | Description |
|------|---------|-------------|
| `archive` | required | Path to `.gpk` archive |
| `--max-fmh-z` | 5.0 | Remove genomes whose per-genus FMH-minority z (median+MAD delta from genus baseline) exceeds this; the genus-band signal |
| `--min-delta` | 0.02 | Also require `fmh_minority` to exceed the genus baseline by this absolute amount (guards tight-baseline genera from z blow-up) |
| `--max-cco` | 1.5 | Remove genomes whose CCO contamination % (per-contig T²∪SPE TNF-outlier bp vs the genus GCOV covariance) exceeds this; the distant-band signal (family/order/phylum). Needs contigs ≥20kb |
| `--max-iters` | 5 | Max decontamination rounds |
| `-t / --threads` | 8 | Threads |
| `--dry-run` | off | Report what would be removed (with fmh_z/cco per axis) without modifying the archive |

---

## `gcov`

Build (or rebuild) the GCOV (per-genus), FCOV (per-family) and FMHR (per-genus FracMinHash) sections in an existing `.gpk` archive; all three are rewritten together.

```bash
genopack gcov mydb.gpk
```

| Flag | Default | Description |
|------|---------|-------------|
| `archive` | required | Path to `.gpk` archive |
| `-t / --threads` | 8 | Parallel shard readers |

---

## `pcore`

Build SEC_PCORE, the dense per-genus aamer reference containing every aamer and its prevalence. It detects foreign contigs of 1-2kb; the conserved-only CORE is too sparse at that length.

```bash
genopack pcore mydb.gpk
```

| Flag | Default | Description |
|------|---------|-------------|
| `archive` | required | Path to `.gpk` archive or directory of parts |
| `-t / --threads` | 8 | Parallel shard readers |
| `--members` | all live | Reference accession list (one per line); only these genomes build the reference |
| `--theta` | `0.90` | Conserved-core prevalence threshold; used only when the archive has no SEC_CORE |
| `--min-seg-aa` | `8` (= k) | Min inter-stop AA segment length; used only when no SEC_CORE |
| `--frac-max-hash` | `UINT64_MAX` (dense) | FracMinHash max hash; used only when no SEC_CORE |

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

Protein-aamer clade-split contamination tier (GUNC-style). Build a panel of genus-diagnostic aamers from reference genomes, then score genomes by lineage split.

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
| `-i / --input` | | TSV (accession, file_path), or use `--gpk` |
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
| `--dayhoff6` | off | Build Dayhoff-6 k=12 syncmer profile pool (compact and less sensitive to divergence) |
| `--ic-threshold` | 0.25 | Min per-column IC fraction to include k-mer positions (`--dayhoff6` only) |
| `--expected-min-frac` | 0.50 | A marker is counted as expected for a genus only if detectable (≥1 IC-passing syncmer) in at least this fraction of the genus's reference genomes; mirrors CheckM2's single-copy universality criterion |

### `markers remerge`

```bash
genopack markers remerge markers.mrk
```

Append the pre-merged pool to an existing `.mrk` panel (in place) so `MarkerReader` uses its zero-copy mmap path instead of copying ~450MB per worker. This does not require source genomes.

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
| `--min-hits` | 1 | Floor on hits per marker; the effective call threshold is `max(--min-hits, query_kmers / 20)` |

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

There is no `genopack similar` or `genopack cidx` CLI. KMRX profiles and CIDX contig→genome lookup are exposed through the C++ API only:

- `ArchiveReader::kmer_profile(...)` / `kmer_profile_by_accession(...)` plus `cosine_similarity(a, b)`: there is no built-in top-k search.
- `ArchiveReader::find_contig_genome_id(accession)`: CIDX binary search, ~150M queries/s/core.

See [API → Similarity & contig lookup](api.md#similarity-contig-lookup).

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
| `--no-gstx` | off | Skip the GSTX genus-stats index |
| `--no-qual` | off | Skip the QUAL per-genome quality build |
| `--gcov` | off | Rebuild GCOV/FCOV/FMHR (full shard rescan). When GSTX and GCOV both exist afterwards, reindex clears the stale-stats header flag that `merge` sets |
| `--txdb` | off | Build the taxonomy tree (TXDB) from existing TAXN lineage strings |
| `--cidx FILE` | unset | Build the contig accession index (CIDX) from a build TSV (`accession\ttaxonomy\tfile_path`) |
| `--cidx-threads` | 8 | Threads for parallel FASTA decompression while indexing contigs |
| `--skch` | off | Re-sketch **every** genome into one coalesced SKCH section; pre-existing SKCH sections are dropped |
| `--skch-threads` | 8 | Threads for parallel sketch computation |
| `--sketch-kmer` | inherit / 16 | OPH k-mer size for a single-k SKCH section |
| `--sketch-kmers` | unset (single-k) | Comma list (e.g. `16,21,31`) → multi-k SKCH v2; takes precedence over `--sketch-kmer` when it has more than one value |
| `--sketch-size` | inherit / 10000 | OPH sketch size |
| `--sketch-syncmer` | `-1` (auto) | Open-syncmer prefilter `s` (0 disables). Auto is `max(2, k/3)` with k from `--sketch-kmer` (or 16), or the smallest `--sketch-kmers` value; it is never inherited from existing SKCH sections (only k and size are). When k is inherited, auto `s` still uses `--sketch-kmer` or 16, not the inherited k |

Typical uses: an archive built with `build --no-cidx` later wants contig lookup (`--cidx genomes.tsv`); a TAXN-only archive needs the tree (`--txdb`); SKCH layout needs to be upgraded to V4 seekable (`--skch --force`); or a multi-k variant is needed (`--skch --sketch-kmers 16,21,31 --force`).

---

## `verify`

```bash
genopack verify mydb.gpk [--strict] [--coverage-only] [--verbose]
```

Runs two checks. **Semantic coverage** cross-checks per-genome and per-genus section row counts against the catalog from the TOC only, without a data scan. SKCH (per k, every k must cover the whole catalog), QUAL, GIDX, ACCX, and CATL rows must equal the genome count. GSTX and other compute sections are reported as presence + coverage, not gated. **Byte-level checksums** cover TailLocator, TocHeader, and per-section XXH128 content. The command exits 0 if all checks pass and non-zero on any failure.

| Flag | Default | Description |
|------|---------|-------------|
| `archive` | required | Archive directory to verify |
| `--strict` | off | Fail if any section lacks a content checksum |
| `--coverage-only` | off | Run only the section-coverage cross-check; skip byte-level checksums (seconds) |
| `--verbose` | off | Print OK lines in addition to failures |

---

## `coordinator`

NFS-coordinated assembly mode for distributed builds. Workers run `build --coordinator <manifest_dir>:<output.gpk>`; the coordinator process waits for the expected number of worker manifests, then merges parts into a single archive.

```bash
genopack coordinator -o mydb.gpk --workers 64 --nfs-dir /shared/manifests/
```

| Flag | Default | Description |
|------|---------|-------------|
| `-o / --output` | required | Final merged archive path |
| `--workers` | required | Expected number of worker manifests |
| `--nfs-dir` | required | Shared directory where workers drop manifests |
