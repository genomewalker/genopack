# Contig Index (CIDX)

The CIDX section maps contig accession strings (e.g. `NZ_JAVJIU010000001.1`) to the `genome_id` of the containing assembly. This enables efficient lookup from any reference sequence seen in metagenomic alignments to the genome record in the archive.

## Section layout

```
CidxHeader (32 bytes)
CidxEntry[] sorted ascending by acc_hash (16 bytes each)
```

```c
struct CidxHeader {
    uint32_t magic;        // SEC_CIDX = 0x58444943 ("CIDX")
    uint32_t version;      // 1
    uint64_t n_entries;
    uint64_t batch_id;     // monotonically increasing; used for multi-section merge order
    uint8_t  reserved[8];
};

struct CidxEntry {
    uint64_t acc_hash;     // FNV-1a-64(contig_accession_string)
    uint32_t genome_id;
    uint32_t _pad;
};
```

**Hash function:** FNV-1a-64 — same as ACCX and the columnar taxonomy binary, so lookup code is shared across all index types.

**Multiple sections:** incremental `add` operations append a new CIDX section per batch. `MergedCidxReader` searches newest→oldest and returns the first match. `genopack merge` consolidates all CIDX sections into one.

## Build integration

CIDX is built automatically during `genopack build`. For each genome, the builder parses every `>header` line in the FASTA and records the first whitespace-delimited token as the contig accession:

```
>NZ_JAVJIU010000001.1 Christiangramia sp. SM2212 ...
 ^^^^^^^^^^^^^^^^^^^^
 contig accession stored in CIDX
```

## Library API

```cpp
#include <genopack/archive.hpp>

genopack::ArchiveReader ar;
ar.open("mydb.gpk");

// Single lookup — O(log n) binary search across all CIDX sections
uint32_t gid = ar.find_contig_genome_id("NZ_JAVJIU010000001.1");
if (gid != UINT32_MAX) {
    std::string acc = ar.accession_for_genome_id(gid);
    auto tax = ar.taxonomy_for_accession(acc);
}

// Batch lookup — parallel hash + sort + merge-join
std::vector<std::string_view> contigs = { "NZ_ABC.1", "NZ_DEF.1", /* ... */ };
std::vector<uint32_t> genome_ids(contigs.size());
ar.batch_find_contig_genome_ids(contigs.data(), genome_ids.data(),
                                contigs.size(), /*n_threads=*/8);
// genome_ids[i] == UINT32_MAX if not found
```

## Batch lookup algorithm

`batch_find_contig_genome_ids` runs in three phases:

1. **Hash** — compute FNV-1a-64 for all N accessions in parallel across T threads — O(N/T)
2. **Sort** — sort queries by hash (pdqsort) — O(N log N)
3. **Merge-join** — for each CIDX section (newest→oldest), merge-join sorted queries against the sorted section entries — O(N + M) per section

Within the merge-join, two optimisations avoid scanning the full entry array:

- **Range clip:** binary-search the entry array to `[ei_start, ei_end)` covering only the hash range of the current query batch — free for clustered queries (e.g. one organism)
- **Hardware prefetch:** `__builtin_prefetch` 32 entries (8 cache lines) ahead on both the query and entry arrays

After each section pass, only unresolved queries enter the next pass.

## Performance (GTDB r226 reps, 143k genomes, 23.7M contig entries)

| Batch size | Threads | Per-query |
|-----------|---------|-----------|
| 1,000 | 1 | 68 ns |
| 10,000 | 1 | 6.8 ns |
| 10,000 | 4 | 6.9 ns |
| single lookup | — | ~250 µs (cold mmap) |

Thread scaling is flat at this archive size because the bottleneck is the merge-join sequential pass, not hashing. Parallelism helps at N > 100k with multi-section archives.
