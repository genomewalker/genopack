# Changelog

## Unreleased

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
