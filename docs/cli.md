# CLI Reference

## Synopsis

```
genopack <command> [options]
```

All commands accept `--help` / `-h` for option details.

---

## `build`

Build a new archive from a TSV.

```bash
genopack build -i genomes.tsv -o mydb.gpk [options]
```

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | required | Input TSV (`accession  file_path  [completeness  contamination  taxonomy ...]`) |
| `-o / --output` | required | Output `.gpk` path |
| `-t / --threads` | 1 | Parallel FASTA decompression threads |
| `-z / --zstd-level` | 6 | zstd compression level (1–22) |
| `--no-hnsw` | off | Skip HNSW index (recommended for > 1M genomes) |
| `--no-cidx` | off | Skip CIDX contig index |
| `--taxonomy-group` | off | Bucket genomes by genus before shard formation |
| `--taxonomy-rank` | `g` | `g` = genus, `f` = family (requires `--taxonomy-group`) |

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
genopack stat mydb.gpk
```

Output: generation, shard count, live/total genome count, total bp, compression ratio, per-section inventory.

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

Remove near-identical sequences.

```bash
genopack dedup mydb.gpk -o mydb_dedup.gpk [--min-ani 0.999]
```

Uses KMRX cosine similarity as pre-filter, then exact comparison. Keeps the genome with the highest completeness.

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

Query taxonomy for one genome.

```bash
genopack taxonomy mydb.gpk --accession GCA_000008085.1
```

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

## `similar`

Find similar genomes by KMRX cosine similarity.

```bash
genopack similar mydb.gpk GCA_000008085.1 -k 20 --min-sim 0.90
```

Uses the HNSW approximate nearest-neighbour index (if present) or falls back to linear scan.

| Flag | Default | Description |
|------|---------|-------------|
| `-k / --top-k` | 10 | Number of results |
| `--min-sim` | 0.0 | Minimum cosine similarity threshold |

---

## `cidx`

Look up genome_id from contig accession via the CIDX index.

```bash
# Single
genopack cidx mydb.gpk --accession NZ_JAVJIU010000001.1

# Batch
genopack cidx mydb.gpk --accessions-file contigs.txt --threads 8 -o results.tsv
```

Throughput: ~150M queries/s on a single core (binary search on sorted FNV-1a-64 hash array).

---

## `reindex`

Append a missing GIDX or HNSW section to an existing archive.

```bash
genopack reindex mydb.gpk [--hnsw] [--gidx]
```

Useful when an archive was built with `--no-hnsw` and the index is needed later.
