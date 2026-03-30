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

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | required | Input TSV |
| `-o / --output` | required | Output `.gpk` path |
| `-t / --threads` | 1 | Parallel decompression threads |
| `-z / --zstd-level` | 6 | zstd compression level (1–22) |
| `--no-hnsw` | off | Skip HNSW index (recommended for > 1M genomes) |
| `--no-cidx` | off | Skip contig index (saves memory for large archives) |
| `--taxonomy-group` | off | Group genomes by genus before shard formation |

---

## Archive statistics

```bash
genopack stat mydb.gpk
```

Output includes: generation, shard count, genome count, total bp, compression ratio, section inventory.

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

```bash
# Lookup one genome
genopack taxonomy mydb.gpk --accession GCA_000008085.1

# Export as NCBI taxdump (Kraken/Kaiju compatible)
genopack taxdump mydb.gpk -f taxdump -o ./taxdump/

# Export columnar binary (fast offline lookup)
genopack taxdump mydb.gpk -f columnar -o ./taxonomy/
```

The columnar binary export produces four files:

| File | Description |
|------|-------------|
| `acc2taxid.bin` | Sorted `(FNV-1a-64(acc), taxid)` pairs - O(log n) binary search |
| `taxnodes.bin` | Sorted node records + name pool - O(log n) taxid lookup |
| `acc2taxid.tsv` | `accession\ttaxid` - pandas/polars/R compatible |
| `taxonomy.tsv` | `taxid\tparent_taxid\trank\tname\tis_synthetic` |

---

## Find similar genomes

```bash
genopack similar mydb.gpk GCA_000008085.1 -k 20 --min-sim 0.90
```

Uses HNSW approximate nearest-neighbour on KMRX k=4 tetranucleotide profiles. Results are sorted by cosine similarity descending.

---

## Contig accession lookup

```bash
# Single contig → genome_id
genopack cidx mydb.gpk --accession NZ_JAVJIU010000001.1

# Batch (one contig per line)
genopack cidx mydb.gpk --accessions-file contigs.txt --threads 8
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

```bash
.scripts/gpk-build-distributed.sh genomes.tsv output.gpk \
    -t 24 -z 3 \
    node1 node2 node3 node4
```

Each node builds its slice to local NVMe (`/scratch`), rsyncs the part back, then a parallel merge (`pwrite` one thread per part) assembles the final archive. NFS readahead efficiency is maintained by reading each part archive sequentially.

---

## Deduplication

```bash
genopack dedup mydb.gpk -o mydb_dedup.gpk --min-ani 0.999
```

Removes near-identical sequences using KMRX cosine similarity as a pre-filter, then exact sequence comparison. Keeps the genome with the highest completeness.
