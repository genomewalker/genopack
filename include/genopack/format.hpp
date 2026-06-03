#pragma once
#include <cstdint>

namespace genopack {

static constexpr uint32_t GPK2_MAGIC = 0x324B5047u; // "GPK2"
static constexpr uint32_t TOCB_MAGIC = 0x42434F54u; // "TOCB"
static constexpr uint32_t GPKT_MAGIC = 0x544B5047u; // "GPKT"
static constexpr uint16_t FORMAT_MAJOR = 2;
static constexpr uint16_t FORMAT_MINOR = 0;

static constexpr uint32_t SEC_SHRD = 0x44524853u; // "SHRD"
static constexpr uint32_t SEC_CATL = 0x4C544143u; // "CATL"
static constexpr uint32_t SEC_ACCX = 0x58434341u; // "ACCX"
static constexpr uint32_t SEC_TAXN = 0x4E584154u; // "TAXN"
static constexpr uint32_t SEC_TXDB = 0x42445854u; // "TXDB"
static constexpr uint32_t SEC_TOMB = 0x424D4F54u; // "TOMB"
static constexpr uint32_t SEC_TOCB = 0x42434F54u; // "TOCB"
static constexpr uint32_t SEC_KMRX = 0x58524D4Bu; // "KMRX"
static constexpr uint32_t SEC_GIDX = 0x58444947u; // "GIDX" - genome_id index
static constexpr uint32_t SEC_TIDX = 0x58444954u; // "TIDX" - taxonomy chunk bitmaps (future)
static constexpr uint32_t SEC_HNSW = 0x574E5348u; // "HNSW" - ANN index (future)
static constexpr uint32_t SEC_CIDX = 0x58444943u; // "CIDX" - contig accession index
static constexpr uint32_t SEC_SKCH    = 0x48434B53u; // "SKCH" - OPH sketch store (format detected by payload magic)
static constexpr uint32_t SEC_NTDB = 0x4244544Eu; // "NTDB" - NCBI taxonomy database (complete tree)
static constexpr uint32_t SEC_GTAX = 0x58415447u; // "GTAX" - taxonomy alias/redirect table (merges, renames, splits)
static constexpr uint32_t SEC_GSTX = 0x58545347u; // "GSTX" - genus sketch stats (consensus + p90 + TNF centroid)
static constexpr uint32_t SEC_QUAL = 0x4C415551u; // "QUAL" - per-genome quality scores
static constexpr uint32_t SEC_GCOV = 0x564F4347u; // "GCOV" - per-genus biological covariance (whitening basis + calibrated quantiles)
static constexpr uint32_t SEC_FCOV = 0x564F4346u; // "FCOV" - per-family biological covariance (same layout as GCOV, keyed by family hash)
static constexpr uint32_t SEC_FMHR = 0x52484D46u; // "FMHR" - per-genus FracMinHash reference sketches (k=21,c=125)
static constexpr uint32_t SEC_BPRM = 0x4D525042u; // "BPRM" - self-describing build parameters (one per archive)

// Model/pipeline versions folded into BPRM.params_hash (and GPK3 derivation_hash).
// Bump when a section's computation semantics change in a byte-incompatible way.
static constexpr uint32_t GSTX_MODEL_VERSION = 1;
static constexpr uint32_t GCOV_MODEL_VERSION = 1;
static constexpr uint32_t TNF_ORDER          = 4;

// ── FileHeader: 128 bytes ─────────────────────────────────────────────────────
// 4+2+2+8+8+8+8+88 = 128
struct FileHeader {
    uint32_t magic;             // GPK2_MAGIC
    uint16_t version_major;
    uint16_t version_minor;
    uint64_t file_uuid_lo;
    uint64_t file_uuid_hi;
    uint64_t created_at_unix;
    uint64_t flags;
    uint8_t  reserved[88];
};
static_assert(sizeof(FileHeader) == 128);

// ── SectionDesc: 80 bytes ─────────────────────────────────────────────────────
// 4+2+2+8+8+8+8+8+8+8+16 = 80
struct SectionDesc {
    uint32_t type;
    uint16_t version;
    uint16_t flags;
    uint64_t section_id;
    uint64_t file_offset;
    uint64_t compressed_size;
    uint64_t uncompressed_size;
    uint64_t item_count;
    uint64_t aux0;
    uint64_t aux1;
    uint8_t  checksum[16];
};
static_assert(sizeof(SectionDesc) == 80);

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

// ── TailLocator: 64 bytes ─────────────────────────────────────────────────────
// 4+2+2+8+8+8+16+16 = 64
// reserved[0..7]  = meta_bundle_offset (u64, 0 = not present)
// reserved[8..15] = meta_bundle_size   (u64, 0 = not present)
// When both are non-zero, a single pread of [meta_bundle_offset, +meta_bundle_size)
// covers all hot sections (CATL, ACCX, QUAL, GSTX, GIDX, CIDX, TAXN) contiguously.
struct TailLocator {
    uint32_t magic;             // GPKT_MAGIC
    uint16_t version;
    uint16_t flags;
    uint64_t toc_offset;
    uint64_t toc_size;
    uint64_t generation;
    uint8_t  toc_checksum[16];
    uint8_t  reserved[16];

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
static_assert(sizeof(TailLocator) == 64);

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
