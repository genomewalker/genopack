# Changelog

## 1.1.0 — 2026-07-18

### Quality tier — contamination decoupled

- **`quality_tier` is completeness-only.** The old NA-safe max over seven correlated contamination axes demoted 362,935 GTDB r232 genomes at 6.9% PPV against CheckM2; it is removed. LQ is now genuine incompleteness (`comp_eff < 0.50`) only. The single CheckM2-calibrated duplication channel `D` may cap HQ→MQ, and nothing forces LQ (`src/check/run_check.cpp:181-197`).
- **`quality_score` is completeness-only** (`compute_quality_score(comp_eff)`), replacing the multiplicative completeness×contamination score.
- GTDB r232 tiers: HQ 5,257,430 (55.16%), MQ 3,737,697 (39.22%), LQ 535,855 (5.62%); zero missing.

### Contamination — three reported channels

- `contamination_channels` reports D (duplication, the only CheckM2-calibrated axis), S (FracMinHash minority), G (median of the correlated geometry axes) — reported for dereplication's D→S→G tiebreak, never a discard gate.
- New TSV columns: `contam_D`, `contam_S`, `contam_G`, `contam_score` (noisy-OR, display only), `channels_fired`, `contamination_tnf_minor`.

### Completeness — marker scoring for all genomes

- **`--markers` scores marker completeness for every genome** (`score_all_completeness = scan_all || markers given`), not only pass-A-flagged ones. Marker completeness populated for ~99% of GTDB r232 (was ~34%); validated against CheckM2 at slope 0.81, pooled bias ~0.
- `completeness_effective` priority: `marker → aamer_core → post_decontam`, with `cluster_relative` admitted only as a soft corroborator.

### Pack write-safety

- `check` unions carried-forward QUAL records instead of overwriting them, and refuses to drop build-time `core_dup_mass` — a subset rescore can no longer silently delete untouched genomes.
- New `--dup-restore` flag re-injects the build-time duplication axis from a quality or `cladesplit score` TSV.
- `--scan-all` forces every genome through pass-B (implied by `--markers`).

## 1.0.0 — 2026-05-10

- **Directory archive layout**: `.gpk` is now a directory of section files plus `toc.bin`. Single-file `.gpk` archives are no longer produced; section binary structures (CATL, GIDX, ACCX, …) are unchanged.
- **Multipart sets via `ArchiveSetReader`**: open a directory of `part_*.gpk` archives transparently. Per-accession routing across parts, aggregated stats, multipart-aware `extract` / `fetch_by_accession` / `batch_fetch` / `fetch_sequence_slice` / `taxonomy_for_accession`. `LocatedGenomeMeta` carries `(part_index, part_path)` so callers can disambiguate part-local `genome_id` values that may collide across parts.
- **SKCH V4 (seekable, dual-seed, multi-k)**: replaces V1/V2/V3 (now hard-rejected). 96-byte `SkchSeekHdr`, 16-byte `SkchFrameDesc`, 16 384 genomes per zstd-compressed frame, planar (n_real_bins / sigs1 / sigs2 / masks) with seeds default 42/43, up to 8 k-mer sizes (`--sketch-kmers 16,21,31`).
- **`.gpd` Geodesic Derep Archive Format v1** + `DerepView` reader: representatives, rep-only embedding matrix (default f16), `accession_set_hash` (xxh3-64 of sorted live accessions joined with `\n`), staleness levels (`Valid` / `LayoutChangedSameLiveSet` / `StaleNewGenomes` / `StaleTombstones` / `Mismatch`).
- **NTDB section**: optional embedded NCBI `nodes.dmp` + `names.dmp` (zstd) for offline taxid resolution; built by `coordinator --ntdb DIR`.
- **`inspect` subcommand**: reports SKCH layout (V3 vs V4 seekable, `n_kmers`, `sketch_size`, frame size, dual-seed presence) and preload memory cost. Accepts a single archive or a multipart directory.
- **`taxonomy` subcommand group**: replaces the single `genopack taxonomy` invocation with `show` / `normalize` / `partition` / `assign-taxids` / `diff` / `patch` (the last accepts GTDB-Tk classify_summary directly via `--gtdbtk`).
- **`coordinator --ntdb`**: NFS-coordinated distributed build can now embed an NTDB section in the merged archive.
- **`reindex` extended**: `--force`, `--no-gidx`, `--txdb`, `--cidx FILE` (+ `--cidx-threads`), `--skch` with full `--sketch-kmer(s)` / `--sketch-size` / `--sketch-syncmer` (+ `--skch-threads`). Use `--skch --force` to upgrade legacy SKCH layouts to V4 in place.
- **`dedup` is now in-place**: tombstones duplicate genomes by content hash directly in the archive (`--dry-run` to preview); physical bytes reclaimed only by `repack`.
- **CLI surface trimmed**: `genopack similar` and `genopack cidx` removed. The functionality remains in the C++ API (`ArchiveReader::find_similar` / `find_similar_by_accession` / `find_contig_genome_id`); HNSW is no longer built by default and is library-only.
- **Default build flags**: `--taxon-group` and `--taxon-rank` (genus) are on by default. Note: `taxon-rank` for `build` vs `--taxonomy-rank` for `repack` (intentional, separate options).
- **`stat --json` and multipart support**: `stat` accepts either a single archive directory or a multipart set; `--json` emits machine-readable output.
- **`repack` command**: two-pass taxonomy-grouped reshard for fast per-taxon NFS access
  - Phase 1 reads only `GenomeDirEntry` arrays (~300 MB) instead of full blobs; builds complete routing table in minutes
  - Phase 2 sorts all genome records by `(taxonomy, oph_fingerprint)` in memory
  - Phase 3 single-pass decompress + route with OMP parallelism; smart cap eviction (largest writer only) minimises shard fragmentation
  - `--taxonomy-rank g|f` flag (genus default, family fallback)
  - `-m / --max-memory` flag to control eviction threshold
- OpenMP parallel decompression in repack Phase 3
- **Fix: `MmapFileReader::advise()` page-alignment** — `madvise()` requires a page-aligned address; the raw `data_ + offset` pointer was not aligned, so every `MADV_DONTNEED` call on SKCH sections silently returned `EINVAL`. Fixed by rounding the address down to the nearest 4 KiB page boundary and extending the length accordingly. On large runs (e.g. GTDB r232, 5.2 M genomes) this caused sketch mmap pages to accumulate across waves, inflating peak RSS by ~200–300 GiB.

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
