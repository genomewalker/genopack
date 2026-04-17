#pragma once
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>

namespace genopack {

// ── NtdbHeader: 64 bytes ──────────────────────────────────────────────────────
// The NTDB section embeds the raw NCBI taxdump (nodes.dmp + names.dmp) as two
// consecutive zstd-compressed blobs. This preserves the full NCBI tree for
// provenance and offline NCBI taxid resolution without any lossy transformation.
struct NtdbHeader {
    uint32_t magic;             // SEC_NTDB
    uint16_t version;           // 1
    uint16_t flags;
    uint64_t taxdump_date;      // YYYYMMDD as uint64 (0 if unknown)
    uint64_t nodes_raw_size;    // uncompressed size of nodes.dmp
    uint64_t nodes_zstd_size;   // compressed size (immediately follows header)
    uint64_t names_raw_size;    // uncompressed size of names.dmp
    uint64_t names_zstd_size;   // compressed size (follows nodes blob)
    uint8_t  _pad[16];
};
static_assert(sizeof(NtdbHeader) == 64);

// Layout: [NtdbHeader(64)] [zstd(nodes.dmp)] [zstd(names.dmp)]

// ── NtdbWriter ────────────────────────────────────────────────────────────────
// Reads nodes.dmp and names.dmp from a taxdump directory, zstd-compresses them,
// and writes a NTDB section.
class NtdbWriter {
public:
    // taxdump_dir must contain nodes.dmp and names.dmp.
    // taxdump_date: YYYYMMDD (e.g. 20240401); 0 = unknown.
    void load(const std::filesystem::path& taxdump_dir, uint64_t taxdump_date = 0);

    SectionDesc finalize(AppendWriter& writer, uint64_t section_id);

private:
    std::string nodes_raw_;
    std::string names_raw_;
    uint64_t    taxdump_date_ = 0;
};

// ── NtdbReader ────────────────────────────────────────────────────────────────
// Provides access to the embedded NCBI taxdump bytes (decompressed on demand).
class NtdbReader {
public:
    void open(const uint8_t* base, uint64_t offset, uint64_t size);

    // Decompress and return nodes.dmp content
    std::string nodes_dmp() const;

    // Decompress and return names.dmp content
    std::string names_dmp() const;

    uint64_t taxdump_date() const { return hdr_ ? hdr_->taxdump_date : 0; }
    bool     valid()        const { return hdr_ != nullptr; }

private:
    const NtdbHeader* hdr_         = nullptr;
    const uint8_t*    nodes_blob_  = nullptr;  // points into mmap
    const uint8_t*    names_blob_  = nullptr;
};

} // namespace genopack
