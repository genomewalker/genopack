#pragma once
#include <cstddef>
#include <cstdint>

namespace genopack {

static constexpr uint32_t GPK_MAGIC = 0x334B5047u; // "GPK3" — container magic (GPK3 format break)
static constexpr uint32_t TOCB_MAGIC = 0x42434F54u; // "TOCB"
static constexpr uint32_t GPKT_MAGIC = 0x544B5047u; // "GPKT"
static constexpr uint16_t FORMAT_MAJOR = 3;
static constexpr uint16_t FORMAT_MINOR = 0;

// Byte-order / ABI canary. A correctly-byte-ordered reader reads this back
// verbatim; a mismatch means the file was written on an incompatible ABI (P26).
static constexpr uint64_t ENDIAN_ABI_TAG = 0x0102040810204080ull;

static constexpr uint32_t SEC_SHRD = 0x44524853u; // "SHRD"
static constexpr uint32_t SEC_CATL = 0x4C544143u; // "CATL"
static constexpr uint32_t SEC_ACCX = 0x58434341u; // "ACCX"
static constexpr uint32_t SEC_TAXN = 0x4E584154u; // "TAXN"
static constexpr uint32_t SEC_TXDB = 0x42445854u; // "TXDB"
static constexpr uint32_t SEC_TOMB = 0x424D4F54u; // "TOMB"
static constexpr uint32_t SEC_TOCB = 0x42434F54u; // "TOCB"
static constexpr uint32_t SEC_KMRX = 0x58524D4Bu; // "KMRX"
static constexpr uint32_t SEC_GIDX = 0x58444947u; // "GIDX" - genome_id index
static constexpr uint32_t SEC_CIDX = 0x58444943u; // "CIDX" - contig accession index
static constexpr uint32_t SEC_SKCH    = 0x48434B53u; // "SKCH" - OPH sketch store (format detected by payload magic)
static constexpr uint32_t SEC_GSTX = 0x58545347u; // "GSTX" - genus sketch stats (consensus + p90 + TNF centroid)
static constexpr uint32_t SEC_QUAL = 0x4C415551u; // "QUAL" - per-genome quality scores
static constexpr uint32_t SEC_GCOV = 0x564F4347u; // "GCOV" - per-genus biological covariance (whitening basis + calibrated quantiles)
static constexpr uint32_t SEC_FCOV = 0x564F4346u; // "FCOV" - per-family biological covariance (same layout as GCOV, keyed by family hash)
static constexpr uint32_t SEC_FMHR = 0x52484D46u; // "FMHR" - per-genus FracMinHash reference sketches (k=21,c=125)
static constexpr uint32_t SEC_BPRM = 0x4D525042u; // "BPRM" - self-describing build parameters (one per archive)
// ── Provenance-first columnar quality store (replaces flat SEC_QUAL) ──
// Cell = (genome, axis, source/tool, method). Intrinsic and external measures
// live in separate sections and are never silently fused; CORE holds the
// per-genus prevalence models the intrinsic columns reference.
static constexpr uint32_t SEC_CORE = 0x45524F43u; // "CORE" - per-genus prevalence cores (content-hashed models)
static constexpr uint32_t SEC_PCORE = 0x524F4350u;// "PCOR" - unified per-genus aamer reference: ALL aamers + u8 prevalence (dense; supersedes CORE)
static constexpr uint32_t SEC_GAMI  = 0x494D4147u;// "GAMI" - global aamer multiplicity index: BlockedBloom of rare aamers (genus-count ≤ K), precomputed at build time
static constexpr uint32_t SEC_QCOL = 0x4C4F4351u; // "QCOL" - intrinsic quality, columnar (genopack-computed)
static constexpr uint32_t SEC_XQAL = 0x4C415158u; // "XQAL" - external quality, columnar (CheckM2/anvi'o/...)
static constexpr uint32_t SEC_CQAL = 0x4C415143u; // "CQAL" - calibrated quality, columnar (excluded from derivation)
static constexpr uint32_t SEC_PROF = 0x464F5250u; // "PROF" - reporting profiles (named, content-hashed fusion policy)
static constexpr uint32_t SEC_QCONTIG = 0x47544351u; // "QCTG" - per-contig quality overlay (which contigs drive a score)

// Model/pipeline versions folded into BPRM.params_hash (and GPK3 derivation_hash).
// Bump when a section's computation semantics change in a byte-incompatible way.
static constexpr uint32_t GSTX_MODEL_VERSION = 1;
static constexpr uint32_t GCOV_MODEL_VERSION = 1;
static constexpr uint32_t TNF_ORDER          = 4;

// ── FileHeader: 256 bytes (GPK3) ──────────────────────────────────────────────
// Offsets 0..39 are byte-identical to GPK2 (magic/version/uuid/created_at/flags);
// the GPK2 reserved[88] region is now named GPK3 fields + a smaller reserved tail.
// All writers/readers use sizeof(FileHeader), so the growth is layout-transparent.
struct FileHeader {
    uint32_t magic;             //   0  GPK_MAGIC (GPK3)
    uint16_t version_major;     //   4
    uint16_t version_minor;     //   6
    uint64_t file_uuid_lo;      //   8  real random UUID (P23)
    uint64_t file_uuid_hi;      //  16
    uint64_t created_at_unix;   //  24
    uint64_t flags;             //  32
    // ── GPK3 fields (were reserved[88]) ──
    uint64_t endian_abi_tag;    //  40  == ENDIAN_ABI_TAG (P26)
    uint64_t build_params_hash; //  48  BPRM params_hash (cross-section param link)
    uint64_t generation;        //  56  durable, bumped on each seal
    uint64_t shard_set_uuid_lo; //  64  shard-set identity
    uint64_t shard_set_uuid_hi; //  72
    uint32_t shard_id;          //  80  durable, assigned at first seal
    uint32_t _pad0;             //  84
    uint64_t dir_back_offset;   //  88  header-side anchor → directory (T2 torn-tail recovery)
    uint64_t dir_back_size;     //  96
    uint8_t  dir_back_xxh128[16];// 104
    uint64_t dir_back_generation;//120  front/back agreement check
    uint8_t  reserved[128];     // 128  pad to 256
};
static_assert(sizeof(FileHeader) == 256);
static_assert(offsetof(FileHeader, endian_abi_tag) == 40);
static_assert(offsetof(FileHeader, dir_back_offset) == 88);

// FileHeader.flags bits.
// Set by merge when it drops the derived genus/quality sections
// (GSTX/GCOV/FCOV/FMHR/QUAL): those are corpus-wide aggregates / per-genome
// scores that cannot be concatenated across inputs (a genus can span inputs, so
// each input's stats are partial). A merged archive carries this marker until
// the sections are re-derived (genopack reindex --gcov). Readers may warn before
// using genus stats while it is set.
static constexpr uint64_t FH_FLAG_STATS_STALE = 0x1;

// ── SectionDesc: 112 bytes (GPK3) ─────────────────────────────────────────────
// Offsets 0..79 are byte-identical to GPK2; the three GPK3 fields are appended.
// `checksum` is the per-section content_xxh128 (P1). derivation_hash /
// semantic_schema_hash live in the directory so --from-gpk can plan the whole
// reuse set from one directory read without touching any section body (§5).
struct SectionDesc {
    uint32_t type;                  //   0
    uint16_t version;               //   4  per-section schema version (enforced)
    uint16_t flags;                 //   6
    uint64_t section_id;            //   8
    uint64_t file_offset;           //  16  bounds-checked at directory load (P5)
    uint64_t compressed_size;       //  24
    uint64_t uncompressed_size;     //  32
    uint64_t item_count;            //  40
    uint64_t aux0;                  //  48
    uint64_t aux1;                  //  56
    uint8_t  checksum[16];          //  64  content_xxh128 (P1)
    // ── GPK3 fields ──
    uint64_t derivation_hash;       //  80  content-addressed reuse key (§5); 0 = SOURCE/unset
    uint64_t semantic_schema_hash;  //  88  ordering/hash/collation identity (§4); 0 = unset
    uint16_t header_size;           //  96  section preamble size (tolerate growth); 0 = legacy
    uint8_t  reserved2[14];         //  98  pad to 112
};
static_assert(sizeof(SectionDesc) == 112);
static_assert(offsetof(SectionDesc, checksum) == 64);
static_assert(offsetof(SectionDesc, derivation_hash) == 80);

// ── TocHeader: 128 bytes ──────────────────────────────────────────────────────
// 4+2+2+8+8+8+8+8+8+8+8+16+40 = 128
struct TocHeader {
    uint32_t magic;                      // TOCB_MAGIC
    uint16_t version;
    uint16_t flags;
    uint64_t generation;
    uint64_t prev_toc_offset;
    uint64_t section_count;
    uint64_t live_genome_count;
    uint64_t total_genome_count;
    uint64_t catalog_root_section_id;
    uint64_t accession_root_section_id;
    uint64_t tombstone_root_section_id;
    uint8_t  checksum[16];
    uint8_t  reserved[40];
};
static_assert(sizeof(TocHeader) == 128);

// ── TailLocator: 128 bytes (GPK3) ─────────────────────────────────────────────
// Offsets 0..63 are byte-identical to GPK2; GPK3 appends file_uuid + directory hash.
// Discovery is a fixed read of the final sizeof(TailLocator) bytes at EOF.
// reserved[0..7]  = meta_bundle_offset (u64, 0 = not present)
// reserved[8..15] = meta_bundle_size   (u64, 0 = not present)
// When both are non-zero, a single pread of [meta_bundle_offset, +meta_bundle_size)
// covers all hot sections (CATL, ACCX, QUAL, GSTX, GIDX, CIDX, TAXN) contiguously.
struct TailLocator {
    uint32_t magic;             //   0  GPKT_MAGIC
    uint16_t version;           //   4
    uint16_t flags;             //   6
    uint64_t toc_offset;        //   8
    uint64_t toc_size;          //  16
    uint64_t generation;        //  24
    uint8_t  toc_checksum[16];  //  32
    uint8_t  reserved[16];      //  48  meta_bundle_offset/size (via accessors below)
    // ── GPK3 fields ──
    uint64_t file_uuid_lo;      //  64  cross-check vs FileHeader.file_uuid (P16/P23)
    uint64_t file_uuid_hi;      //  72
    uint8_t  directory_xxh128[16];// 80  aggregate hash of the directory
    uint8_t  reserved2[32];     //  96  pad to 128

    uint64_t meta_bundle_offset() const {
        uint64_t v; __builtin_memcpy(&v, reserved,     8); return v;
    }
    uint64_t meta_bundle_size() const {
        uint64_t v; __builtin_memcpy(&v, reserved + 8, 8); return v;
    }
    void set_meta_bundle(uint64_t offset, uint64_t size) {
        __builtin_memcpy(reserved,     &offset, 8);
        __builtin_memcpy(reserved + 8, &size,   8);
    }
};
static_assert(sizeof(TailLocator) == 128);
static_assert(offsetof(TailLocator, file_uuid_lo) == 64);

// ── MetaBundle: embedded section directory ────────────────────────────────────
// Written inside the .gpk just before the TOC; pointed to by TailLocator.reserved.
// One pread of this region returns the complete section directory so the reader
// never needs to page-fault near EOF of a multi-TB mmap just to find sections.

static constexpr uint32_t GMET_MAGIC       = 0x54454D47u; // "GMET"
static constexpr uint32_t GMET_TRAIL_MAGIC = 0x6C617274u; // "tral"

struct MetaHeader {                              // 128 bytes
    uint32_t magic;                              //   4  GMET_MAGIC
    uint16_t version;                            //   2  1
    uint16_t flags;                              //   2
    uint32_t n_sections;                         //   4
    uint32_t _pad;                               //   4
    uint64_t generation;                         //   8
    uint64_t gpk_uuid_lo;                        //   8
    uint64_t gpk_uuid_hi;                        //   8
    uint64_t live_genome_count;                  //   8
    uint64_t total_genome_count;                 //   8
    uint64_t catalog_root_section_id;            //   8
    uint64_t accession_root_section_id;          //   8
    uint64_t tombstone_root_section_id;          //   8
    uint8_t  checksum[16];                       //  16  XXH128 of bundle, this field zeroed
    uint8_t  _reserved[32];                      //  32
};
static_assert(sizeof(MetaHeader) == 128);

struct MetaTrailer {                             // 32 bytes — mirror anchor
    uint32_t magic;                              //   4  GMET_TRAIL_MAGIC
    uint32_t n_sections;                         //   4
    uint64_t generation;                         //   8
    uint8_t  _reserved[16];                      //  16
};
static_assert(sizeof(MetaTrailer) == 32);

} // namespace genopack
