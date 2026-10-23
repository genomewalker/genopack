# Changelog

> Status as of 2026-10-23: Unreleased covers main up to 43238fc. Entries for features merged after 1.1.0 that were missing (2026-07-19 to 2026-08-13) were added under "Commands and flags".

## Unreleased

### dedup

- `dedup` confirms fingerprint + length candidates by an XXH3-128 hash of the decoded residues before tombstoning, so genomes that only share a sketch minimum are no longer removed.
- Candidates are sorted by fingerprint then length; before, equal-fingerprint genomes of different lengths could separate true duplicates so they were never compared.
- Genomes already in a TOMB section are skipped; before, a second run could keep an already-tombstoned genome and tombstone the live copy.
- New `--groups` TSV (default `dedup_groups.tsv`) lists every candidate with its hash and action, also under `--dry-run`. Check: `tests/dedup_sequence_confirm.py`.

### Unified taxonomy: stable taxids

- **`taxonomy unify --stable-taxids-from <prior>/taxonomy.tsv`** (7887ad4): synthetic nodes whose full (rank, name) path exists in a prior unify's columnar `taxonomy.tsv` keep that taxid; new nodes are numbered above its maximum. Without it, export ids follow concept-hash order and a reunify after renaming 464,051 species reused 1,215,247 prior taxids for other taxa. Pass it on every reunify of a released taxonomy.
- **`taxonomy unify --lineage-overrides <tsv>`** (834573a): replaces listed accessions' pack TAXN lineages at read time. Duplicate or absent accessions are errors.
- **Unify fails on shared taxids** (365e976): export throws if two nodes would receive one stable or export taxid, instead of silently dropping the second.

### Quality

- **`quality_score` QCOL column** (f88d727, merged ad3be6f): the continuous completeness-only score is stored in QCOL (axis `support`, unit `Fraction01`, U16, 0 = not scored), readable via `ArchiveReader::qcol_reader()`.

### Commands and flags

- **`subset`** (490b34d): copy-only accession subset of an archive or multipart directory; no source FASTA or re-sketching.
- **`extract --concat-contigs`** (a9e61cd) and **`--nmask`** (d1e0bd4): one N-joined record per genome, with a kaiku `KNMASK01` side file.
- **`neighbors`** (eb4228c): a priori nearest-in-DB member per query from derep `.gpd` clusters.
- **`pcore` bootstrap** (7e6b310): archives without `SEC_CORE` build PCORE from `--theta` / `--min-seg-aa` / `--frac-max-hash` instead of being skipped.
- **`taxonomy partition --contiguous` / `--cap`** (3bd0840): taxonomy-contiguous part filling with a per-group size cap.

### Unified taxonomy

- **`taxonomy unify <packs…>`:** reconciles the taxonomies of several independent `.gpk` packs into one taxdump (`nodes`/`names`/`acc2taxid.dmp` or columnar binary + TSV) without merging sequence data. Per-part TAXN reads run in parallel; output is deterministic.
- **NCBI reconciliation** (`--ncbi-taxdump`, also on `build` and `merge`): taxa whose scientific name is unambiguous in the NCBI taxdump keep their native NCBI taxid (`bit63=0`); taxa NCBI does not name get a stable content-hash concept_id (`bit63=1`). The namespaces are disjoint, with no renumbering.
- **`--preserve-taxids`:** accessions with an authoritative NCBI taxid (e.g. organelle host taxids) keep it verbatim, importing the taxid's real NCBI lineage down to a strain leaf. Stale ids remap via `merged.dmp`. On the RefSeq organelle set, accessions matching their authoritative host taxid rose from 98.3% to 99.997%.
- **`--ncbi-backbone` / `--gtdb-replaces-prok`:** optionally folds the entire NCBI tree in as backbone nodes, excising NCBI's own `Bacteria`/`Archaea` descendants so GTDB stays authoritative for prokaryotes.
- **Performance:** caches the NCBI taxdump parse in `<dir>/.ncbi_taxdb.cache` (parallel dmp parse, staleness-checked). A 2.1M-accession, 1.6M-node unify runs in ~40 s with the cache.
- **Fixes:** accession lookup (`find_accession`) and `scan_taxonomy` skip tombstoned genomes, so in a multipart set `extract` reads the live copy in a later part instead of the deleted (empty) one, and taxonomy scans no longer list deleted accessions (`src/archive.cpp`).
- **Fixes:** `NcbiTaxdb` name parsing dropped every name because a `.dmp` record terminator was left on the last field. Taxdump export flattened all taxids to `1..N`; it now keeps native NCBI ids and renumbers only concept_ids above the maximum NCBI taxid.

## 1.1.0: 2026-07-18

### Quality tier: contamination decoupled

- **`quality_tier` is completeness-only.** The old NA-safe max over seven correlated contamination axes demoted 362,935 GTDB r232 genomes at 6.9% PPV against CheckM2; it is removed. LQ is now genuine incompleteness (`comp_eff < 0.50`) only. The single CheckM2-calibrated duplication channel `D` may cap HQ→MQ, and nothing forces LQ (`src/check/run_check.cpp:181-197`).
- **`quality_score` is completeness-only** (`compute_quality_score(comp_eff)`), replacing the multiplicative completeness×contamination score.
- GTDB r232 tiers: HQ 5,257,430 (55.16%), MQ 3,737,697 (39.22%), LQ 535,855 (5.62%); zero missing.

### Contamination: three reported channels

- `contamination_channels` reports D (duplication, the only CheckM2-calibrated axis), S (FracMinHash minority), and G (median of the correlated geometry axes). They are used for dereplication's D→S→G tiebreak, never as a discard gate.
- New TSV columns: `contam_D`, `contam_S`, `contam_G`, `contam_score` (noisy-OR, display only), `channels_fired`, `contamination_tnf_minor`.

### Completeness: marker scoring for all genomes

- **`--markers` scores marker completeness for every genome** (`score_all_completeness = scan_all || markers given`), not only pass-A-flagged ones. Marker completeness populated for ~99% of GTDB r232 (was ~34%); validated against CheckM2 at slope 0.81, pooled bias ~0.
- `completeness_effective` priority: `marker → aamer_core → post_decontam`, with `cluster_relative` admitted only as a soft corroborator.

### Pack write-safety

- `check` unions carried-forward QUAL records instead of overwriting them and refuses to drop build-time `core_dup_mass`. A subset rescore no longer deletes untouched genomes.
- New `--dup-restore` flag re-injects the build-time duplication axis from a quality or `cladesplit score` TSV.
- `--scan-all` forces every genome through pass-B (implied by `--markers`).

## 1.0.0: 2026-05-10

- **Directory archive layout**: `.gpk` is now a directory of section files plus `toc.bin`. Single-file `.gpk` archives are no longer produced; section binary structures (CATL, GIDX, ACCX, …) are unchanged. *(Corrected 2026-10-23: current builders write a single `.gpk` file and `ArchiveReader` maps one file; see [Format](format.md).)*
- **Multipart sets via `ArchiveSetReader`**: open a directory of `part_*.gpk` archives as one set. Per-accession routing across parts, aggregated stats, multipart-aware `extract` / `fetch_by_accession` / `batch_fetch` / `fetch_sequence_slice` / `taxonomy_for_accession`. `LocatedGenomeMeta` carries `(part_index, part_path)` so callers can disambiguate part-local `genome_id` values that may collide across parts.
- **SKCH V4 (seekable, dual-seed, multi-k)**: replaces V1/V2/V3 (now hard-rejected). 96-byte `SkchSeekHdr`, 16-byte `SkchFrameDesc`, 16 384 genomes per zstd-compressed frame, planar (n_real_bins / sigs1 / sigs2 / masks) with seeds default 42/43, up to 8 k-mer sizes (`--sketch-kmers 16,21,31`).
- **`.gpd` Geodesic Derep Archive Format v1** + `DerepView` reader: representatives, rep-only embedding matrix (default f16), `accession_set_hash` (xxh3-64 of sorted live accessions joined with `\n`), staleness levels (`Valid` / `LayoutChangedSameLiveSet` / `StaleNewGenomes` / `StaleTombstones` / `Mismatch`).
- **NTDB section**: optional embedded NCBI `nodes.dmp` + `names.dmp` (zstd) for offline taxid resolution; built by `coordinator --ntdb DIR`. *(Removed 2026-07-09 in be1b440; current archives have no NTDB section and `coordinator` has no `--ntdb`.)*
- **`inspect` subcommand**: reports SKCH layout (V3 vs V4 seekable, `n_kmers`, `sketch_size`, frame size, dual-seed presence) and preload memory cost. Accepts a single archive or a multipart directory.
- **`taxonomy` subcommand group**: replaces the single `genopack taxonomy` invocation with `show` / `normalize` / `partition` / `assign-taxids` / `diff` / `patch` (the last accepts GTDB-Tk classify_summary directly via `--gtdbtk`).
- **`coordinator --ntdb`**: NFS-coordinated distributed build can now embed an NTDB section in the merged archive. *(Removed in be1b440.)*
- **`reindex` extended**: `--force`, `--no-gidx`, `--txdb`, `--cidx FILE` (+ `--cidx-threads`), `--skch` with full `--sketch-kmer(s)` / `--sketch-size` / `--sketch-syncmer` (+ `--skch-threads`). Use `--skch --force` to upgrade legacy SKCH layouts to V4 in place.
- **`dedup` is now in-place**: tombstones duplicate genomes by content hash directly in the archive (`--dry-run` to preview); physical bytes reclaimed only by `repack`.
- **CLI surface trimmed**: `genopack similar` and `genopack cidx` removed. The functionality remains in the C++ API (`ArchiveReader::find_similar` / `find_similar_by_accession` / `find_contig_genome_id`); HNSW is no longer built by default and is library-only. *(Corrected 2026-10-23: `ArchiveReader` has no `find_similar`; only the HNSW section reader had one, and HNSW was removed in be1b440. Rank KMRX profiles with `cosine_similarity`.)*
- **Default build flags**: `--taxon-group` and `--taxon-rank` (genus) are on by default. `build` uses `taxon-rank`, while `repack` uses the separate `--taxonomy-rank` option.
- **`stat --json` and multipart support**: `stat` accepts either a single archive directory or a multipart set; `--json` emits machine-readable output.
- **`repack` command**: two-pass taxonomy-grouped reshard for fast per-taxon NFS access
  - Phase 1 reads only `GenomeDirEntry` arrays (~300 MB) instead of full blobs; builds complete routing table in minutes
  - Phase 2 sorts all genome records by `(taxonomy, oph_fingerprint)` in memory
  - Phase 3 single-pass decompress + route with OMP parallelism; cap eviction of the largest writer only minimises shard fragmentation
  - `--taxonomy-rank g|f` flag (genus default, family fallback)
  - `-m / --max-memory` flag to control eviction threshold
- OpenMP parallel decompression in repack Phase 3
- **Fix: `MmapFileReader::advise()` page alignment:** `madvise()` requires a page-aligned address. The raw `data_ + offset` pointer was not aligned, so every `MADV_DONTNEED` call on SKCH sections returned `EINVAL`. The address is now rounded down to the nearest 4 KiB page boundary, with the length extended accordingly. On large runs (e.g. GTDB r232, 5.2 M genomes), this caused sketch mmap pages to accumulate across waves and increased peak RSS by ~200-300 GiB.

## 0.1.0

- Initial release: single-file `.gpk` format
- `build`, `merge`, `stat`, `extract`, `slice`, `add`, `rm`, `taxonomy`, `taxdump`, `similar`, `cidx`, `reindex`, `dedup` commands
- SHRD / CATL / GIDX / ACCX / CIDX / TAXN / TXDB / KMRX / HNSW / TOMB sections
- Columnar catalog with row-group predicate pushdown
- HNSW approximate nearest-neighbour index on KMRX profiles
- NCBI taxdump export (`names.dmp` / `nodes.dmp` / `acc2taxid.dmp`)
- Columnar binary taxonomy export (`acc2taxid.bin` / `taxnodes.bin`)
- Distributed build scripts
- `ScanEngine` with I/O / worker thread separation (NFS and NVMe modes)
- MEM-delta codec (k=31 seed + zstd verbatim residue)
