# `.gpd` — Geodesic Derep Archive Format v1

A `.gpd` file is the canonical artifact of a geodesic dereplication run.
It is a self-contained binary archive that answers, for every genome in
the input pack: "is this a representative? if not, who is? what is the
rep's embedding?" — and detects when the source pack has drifted.

## Goals

1. Single source of truth for derep state — no TSVs, no `meta.tsv`
   `representative` column.
2. Versioned over a frozen source-pack generation set; staleness is
   detectable on every open.
3. Survives genopack repack (accession-keyed; physical locators are
   advisory and self-heal).
4. Append-aware — genomes added to the source after derep return
   `UnknownSinceGeneration(N)`.
5. Self-contained — readable without the source pack for archival
   transport.
6. Fast — O(1) accession→rep_status, O(1) rep_id→embedding.

## File layout

A `.gpd` is a **single file**, not a directory:

```
[FileHeader  64B]
[Section blobs  ...] (HDR, ASTR, ASOF, ARMP, RTBL, G2RM, EMBD, [CSTAT])
[TOC               ]
[TailLocator 16B   ]
```

All multi-byte integers are little-endian. All offsets are byte offsets
from the start of the file. All sections are 8-byte aligned; padding is
zero.

## File header (offset 0, fixed 64 bytes)

```c
struct GpdFileHeader {
    uint32_t magic;          // 'GPDF' = 0x46445047  ('F','D','P','G')
    uint16_t format_major;   // 1
    uint16_t format_minor;   // 0
    uint64_t toc_offset;     // bytes to start of TOC section
    uint64_t toc_size;       // bytes
    uint64_t reserved[6];    // zero
};
```

## TailLocator (last 16 bytes of file)

```c
struct GpdTailLocator {
    uint64_t toc_offset;     // duplicate of header.toc_offset (corruption check)
    uint32_t magic;          // 'GPDT' = 0x54445047
    uint32_t crc32;          // crc32 of TOC bytes
};
```

## Section identifiers

```c
constexpr uint32_t GPD_SEC_HDR  = 0x52444847u; // 'GHDR'
constexpr uint32_t GPD_SEC_AST  = 0x54534147u; // 'GAST'  accession string pool
constexpr uint32_t GPD_SEC_ASO  = 0x4F534147u; // 'GASO'  accession offsets
constexpr uint32_t GPD_SEC_ARM  = 0x4D524147u; // 'GARM'  accession hash map (optional)
constexpr uint32_t GPD_SEC_RTB  = 0x42545247u; // 'GRTB'  rep table
constexpr uint32_t GPD_SEC_G2R  = 0x4D523247u; // 'G2RM'  genome→rep map
constexpr uint32_t GPD_SEC_EMB  = 0x424D4547u; // 'GEMB'  embeddings
constexpr uint32_t GPD_SEC_CST  = 0x54534347u; // 'GCST'  cluster stats (optional)
constexpr uint32_t GPD_SEC_TOC  = 0x434F5447u; // 'GTOC'
```

## TOC section

```c
struct GpdSectionDesc {
    uint32_t type;             // GPD_SEC_*
    uint32_t flags;            // bit0 = zstd-compressed payload
    uint64_t file_offset;      // absolute byte offset to section payload
    uint64_t compressed_size;  // bytes on disk (== uncompressed_size if not compressed)
    uint64_t uncompressed_size;
    uint64_t section_id;       // unique within file (monotonic, starts at 1)
    uint64_t reserved[2];
};

// TOC payload layout:
//   uint32_t magic = 'GPDT'
//   uint32_t n_sections
//   uint32_t crc32_of_following_descs
//   uint32_t pad
//   GpdSectionDesc descs[n_sections]
```

Sections may be compressed (`flags & 1`). The reader transparently
decompresses on first access; the file header above gives uncompressed
sizes so the reader can pre-size buffers.

## HDR section (always present, always uncompressed)

The fingerprint that establishes identity and detects staleness.

```c
struct GpdHeader {
    uint32_t magic;            // 'GPDH' = 0x48445047
    uint16_t format_major;     // 1
    uint16_t format_minor;     // 0
    uint64_t created_at_unix;
    uint8_t  run_id[16];       // UUID v4 — unique per derep run
    uint16_t n_parts;          // source pack parts at derep time
    uint16_t embedding_dim;    // typically 256
    uint8_t  embedding_dtype;  // 0=float32, 1=float16
    uint8_t  has_cstats;       // 1 if CSTAT section is present
    uint8_t  pad0[2];
    uint64_t n_genomes;        // total genomes covered (= sum of part live_counts)
    uint64_t n_reps;
    uint64_t n_unclustered;
    // Followed by:
    //   GpdSourcePart parts[n_parts];
    //   GpdDerepParams params;
};

struct GpdSourcePart {           // 48 bytes
    uint8_t  archive_uuid[16];   // from genopack header
    uint64_t generation;
    uint64_t n_genomes_total;
    uint64_t n_genomes_live;
    uint64_t accession_set_hash; // xxhash3-64 of '\n'-joined sorted live accessions
};

struct GpdDerepParams {
    uint8_t  n_kmer_sizes;
    uint8_t  kmer_sizes[7];      // up to 7; tail zero-padded
    uint32_t sketch_size;
    uint64_t sig1_seed;
    uint64_t sig2_seed;
    float    jaccard_thresh;
    uint16_t geodesic_ver_len;
    uint8_t  pad1[2];
    char     geodesic_ver[];     // not null-terminated; padded to 8
};
```

`accession_set_hash` is the operational identity. Computation: take all
**live** accessions of the part (excluding tombstoned), UTF-8, sort
ASCIIbetically, join with `\n` (no trailing), hash with xxh3-64. The
writer must compute this from the source pack at derep time; the reader
recomputes it from the current pack to validate.

## ASTR section — accession string pool

Concatenated accession strings, sorted ASCIIbetically, no separators
(offsets give boundaries). Includes **all** genomes from the source set
(reps + members + unclustered + tombstoned-at-derep-time? **no, only
live at derep time**).

Compressed with zstd if `flags & 1`. Payload after decompression is the
raw concatenated bytes.

## ASOF section — accession offsets

```c
uint32_t magic = 'GASO';
uint32_t n_genomes;
uint64_t pad;
uint32_t offsets[n_genomes + 1];   // offsets[i+1] - offsets[i] = length of accession i
```

For 5M genomes × 4 bytes = 20 MB. Compressed.

`acc_string(ord)` = `string_view(astr + offsets[ord], offsets[ord+1] - offsets[ord])`.

`acc_to_ord(s)`: binary search over `[0, n_genomes)` comparing against
`acc_string(mid)`. ~22 comparisons for 5M, ~300 ns NFS-cold.

## ARMP section — accession hash map (optional but recommended)

Open-addressed hash table for O(1) `acc_to_ord`.

```c
uint32_t magic = 'GARM';
uint32_t n_buckets;       // power of two, load factor target 0.7
uint32_t hash_seed;       // xxh3-64 seed
uint32_t pad;
struct GpdArmpEntry {
    uint64_t hash;        // upper 64 bits of xxh3-128 of accession
    uint32_t ordinal;     // index into ASOF; 0xFFFFFFFF = empty bucket
    uint32_t pad;
} entries[n_buckets];
```

Lookup: bucket = `hash & (n_buckets - 1)`. Linear probe on collision.
Verify by comparing accession string at the candidate ordinal.

For 5M genomes × 7M buckets × 16 bytes = ~110 MB. Compressed to ~40 MB
with zstd.

If absent (older writers, or `--no-armp` flag), reader falls back to
ASOF binary search.

## RTBL section — rep table

```c
uint32_t magic = 'GRTB';
uint32_t n_reps;
uint64_t pad;
struct GpdRepEntry {           // 24 bytes
    uint32_t rep_acc_ord;      // index into ASOF (the rep is itself a genome)
    uint32_t cluster_size;     // including the rep itself; ≥ 1
    uint64_t source_locator;   // (part_idx<<48)|local_genome_id at derep time, advisory
    uint16_t sketch_kmer;      // which k-mer size produced the winning sketch
    uint8_t  flags;            // bit0=has_embedding(must be 1 in v1)
    uint8_t  pad;
    uint32_t cstat_offset;     // index into CSTAT, or 0xFFFFFFFF if absent
} entries[n_reps];
```

Rep ordering convention: sorted by `rep_acc_ord` ascending. This makes
`rep_id ↔ rep_acc_ord` strictly monotonic and makes embedding-row order
match rep-id order for cache locality in cosine-similarity sweeps.

## G2RM section — genome→rep map

```c
uint32_t magic = 'G2RM';
uint32_t n_genomes;
uint64_t pad;
uint32_t rep_id[n_genomes];    // sentinels: 0xFFFFFFFE = unclustered
                               //            0xFFFFFFFF = tombstoned-at-derep-time (rare)
```

Indexed by genome ordinal (i.e., index into ASOF). For genome i,
`rep_id[i]` is the index into RTBL; for a rep, `rep_id[rep.rep_acc_ord] == self_rep_id`.

5M × 4 bytes = 20 MB raw, ~5 MB zstd.

## EMBD section — embeddings (rep-only)

```c
uint32_t magic = 'GEMB';
uint16_t dim;                  // typically 256
uint8_t  dtype;                // 0=f32, 1=f16
uint8_t  pad0;
uint32_t n_reps;
uint32_t pad1;
// followed by: dtype_bytes(dtype) * dim * n_reps  raw matrix
//   row i corresponds to rep_id i (RTBL ordering)
```

Default dtype is **f16** (270 MB for 546k×256). Loss is acceptable for
cosine search; mapping pipelines that need f32 can request via
`derep --embedding-dtype f32` at write time.

## CSTAT section — cluster stats (optional)

```c
uint32_t magic = 'GCST';
uint32_t n_entries;            // ≤ n_reps; sparse if cstat_offset != UINT32_MAX in RTBL
uint64_t pad;
struct GpdClusterStat {        // 12 bytes
    float    mean_jaccard;     // mean pairwise jaccard within cluster
    float    max_distance;     // max member-to-rep jaccard distance
    uint16_t n_members_used;   // sample size for the stats above
    uint8_t  flags;            // bit0=outlier
    uint8_t  pad;
} stats[n_entries];
```

Optional. Used for QC reports; mapping pipelines ignore.

## Public API (genopack)

```cpp
namespace genopack {

class DerepView {
public:
    static DerepView open(const std::filesystem::path& gpd_file);
    void close();
    bool is_open() const;

    enum class StalenessLevel {
        Valid,
        LayoutChangedSameLiveSet,
        StaleNewGenomes,
        StaleTombstones,
        Mismatch
    };
    StalenessLevel check_against(const ArchiveSetReader& pack) const;

    struct RepStatus {
        enum Kind {
            Representative, Member,
            Unclustered, UnknownSinceGeneration, Absent, Tombstoned
        };
        Kind     kind;
        uint32_t rep_id;            // valid for Rep/Member
        uint64_t since_generation;  // valid for UnknownSinceGeneration
    };
    RepStatus rep_status(std::string_view accession) const;

    struct RepInfo {
        std::string_view accession;
        uint32_t         cluster_size;
        uint64_t         source_locator;
        uint16_t         sketch_kmer;
    };
    std::optional<RepInfo> rep(uint32_t rep_id) const;

    const void* embedding(uint32_t rep_id) const;
    uint16_t    embedding_dim() const;
    uint8_t     embedding_dtype() const;

    std::vector<std::string_view> members_of(uint32_t rep_id) const;

    uint32_t n_reps() const;
    uint64_t n_genomes() const;
    uint64_t n_unclustered() const;
    void scan_reps(const std::function<void(uint32_t rep_id, RepInfo)>& cb) const;

    // Iterate all (accession, rep_status) pairs in genome-ordinal order
    void scan_genomes(const std::function<void(std::string_view acc, RepStatus s)>& cb) const;

    // Header introspection
    std::array<uint8_t, 16> run_id() const;
    uint64_t created_at_unix() const;
    std::string_view geodesic_version() const;
    uint16_t source_n_parts() const;
};

DerepView open_derep(const std::filesystem::path& gpd_file);

} // namespace genopack
```

## Public API (geodesic — writer)

```cpp
namespace geodesic {

struct DerepArchiveBuilderConfig {
    std::filesystem::path output_path;     // ".../my_run.gpd"
    uint16_t              embedding_dim;
    uint8_t               embedding_dtype; // 0=f32, 1=f16
    bool                  emit_armp = true;
    bool                  emit_cstat = false;
    int                   zstd_level = 19;
    std::string           geodesic_version;
};

class DerepArchiveBuilder {
public:
    explicit DerepArchiveBuilder(DerepArchiveBuilderConfig cfg);
    ~DerepArchiveBuilder();

    // Set source pack fingerprint info BEFORE adding genomes.
    void set_source_pack(const genopack::ArchiveSetReader& pack);

    // Set derep params for HDR.
    void set_params(const std::vector<uint8_t>& kmer_sizes,
                    uint32_t sketch_size,
                    uint64_t sig1_seed, uint64_t sig2_seed,
                    float jaccard_thresh);

    // Add one genome with its assigned status.
    // For Representative: rep_accession must equal the genome's own accession;
    //                     embedding ptr is required.
    // For Member: rep_accession is the rep this genome maps to;
    //             embedding ptr is ignored.
    // For Unclustered: rep_accession empty; embedding ignored.
    enum class Kind { Representative, Member, Unclustered };
    void add(std::string_view accession,
             Kind kind,
             std::string_view rep_accession,        // self if Representative
             uint64_t source_locator,               // (part_idx<<48)|local_id
             uint16_t sketch_kmer,                  // 0 if Unclustered
             uint32_t cluster_size,                 // ≥1 for Rep
             const void* embedding_or_null);        // dim×dtype bytes if Rep

    // Finalize: writes all sections and TOC. After this the archive is ready.
    void finalize();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace geodesic
```

## Computing `accession_set_hash`

```cpp
// Given an ArchiveSetReader (open) and a part_index:
uint64_t compute_accession_set_hash(const genopack::ArchiveSetReader& pack,
                                    size_t part_idx) {
    std::vector<std::string> accs;
    accs.reserve(part->n_genomes_live());
    part->scan_genome_accessions([&](std::string_view a, GenomeId id) {
        // Skip tombstoned: caller must check live set
        accs.emplace_back(a);
    });
    std::sort(accs.begin(), accs.end());
    // Stream into xxh3-64
    XXH3_state_t* st = XXH3_createState();
    XXH3_64bits_reset(st);
    for (size_t i = 0; i < accs.size(); ++i) {
        if (i > 0) {
            char nl = '\n';
            XXH3_64bits_update(st, &nl, 1);
        }
        XXH3_64bits_update(st, accs[i].data(), accs[i].size());
    }
    uint64_t h = XXH3_64bits_digest(st);
    XXH3_freeState(st);
    return h;
}
```

## Staleness detection

```
For each part_idx in 0..n_parts:
    If pack does not have part_idx → Mismatch (archive structurally different)
    p = pack.part(part_idx)
    h = compute_accession_set_hash(pack, part_idx)
    if h == HDR.parts[i].accession_set_hash:
        if p.archive_uuid == HDR.parts[i].archive_uuid
           and p.generation == HDR.parts[i].generation:
            continue (Valid for this part)
        else:
            level = max(level, LayoutChangedSameLiveSet)
    else:
        n_added = (count of accs in pack not in .gpd)
        n_dropped = (count of accs in .gpd not in pack)
        if n_added > 0: level = max(level, StaleNewGenomes)
        elif n_dropped > 0: level = max(level, StaleTombstones)
        else: level = max(level, Mismatch)
```

For efficiency, the simple version recomputes the hash. A future
optimization can compare incremental Bloom filters or a Merkle tree.

## Versioning

`format_major` bumps for incompatible changes. `format_minor` bumps for
additive changes (e.g., new optional sections). Readers must reject
unknown `format_major`; readers must accept higher `format_minor` and
ignore unknown sections.

## Open questions for v1

- Should embeddings be stored normalized (unit L2) or raw? **Raw.**
  Mapping pipelines often want raw and normalize on demand.
- Should we store per-genome sketch hashes for verification? **No.**
  That's what genopack SKCH is for; .gpd is derep state, not sketch state.
- Should we support derep-of-derep (hierarchical)? **No in v1.** A v2
  could add a parent-gpd reference.
