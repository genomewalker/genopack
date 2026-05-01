# Changelog

## Unreleased

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
