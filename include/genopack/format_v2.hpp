#pragma once
#include <cstdint>

namespace genopack {

static constexpr uint32_t GPK2_MAGIC = 0x324B5047u; // "GPK2"
static constexpr uint32_t TOCB_MAGIC = 0x42434F54u; // "TOCB"
static constexpr uint32_t GPKT_MAGIC = 0x544B5047u; // "GPKT"
static constexpr uint16_t FORMAT_V2_MAJOR = 2;
static constexpr uint16_t FORMAT_V2_MINOR = 0;

static constexpr uint32_t SEC_SHRD = 0x44524853u; // "SHRD"
static constexpr uint32_t SEC_CATL = 0x4C544143u; // "CATL"
static constexpr uint32_t SEC_ACCX = 0x58434341u; // "ACCX"
static constexpr uint32_t SEC_TOMB = 0x424D4F54u; // "TOMB"
static constexpr uint32_t SEC_TOCB = 0x42434F54u; // "TOCB"

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
struct TailLocator {
    uint32_t magic;             // GPKT_MAGIC
    uint16_t version;
    uint16_t flags;
    uint64_t toc_offset;
    uint64_t toc_size;
    uint64_t generation;
    uint8_t  toc_checksum[16];
    uint8_t  reserved[16];
};
static_assert(sizeof(TailLocator) == 64);

} // namespace genopack
