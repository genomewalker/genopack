#pragma once
#include "format.hpp"
#include "mmap_file.hpp"
#include "types.hpp"
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace genopack {

// SEC_FMHR — per-genus FracMinHash reference section.
// Built by `genopack gcov`; read at check time for FMH minority scoring.
// Uses k=21, c=125 (hardcoded at build time; recorded in header for verification).

// ── On-disk layout ────────────────────────────────────────────────────────────

struct FmhrHeader {           // 64 bytes
    uint32_t magic;           // SEC_FMHR
    uint32_t n_genera;
    uint32_t n_buckets;       // open-addressing table, power of 2
    uint32_t k;               // k-mer size used at build
    uint32_t c;               // density = 1/c used at build
    uint32_t _pad[3];
    uint64_t entries_offset;  // bytes from section start → FmhrEntry[n_genera]
    uint64_t buckets_offset;  // bytes from section start → uint32_t[n_buckets]
    uint64_t hashes_offset;   // bytes from section start → concatenated uint64_t hash arrays
    uint64_t _pad2;
};
static_assert(sizeof(FmhrHeader) == 64);

struct FmhrEntry {            // 24 bytes
    uint64_t genus_hash;
    uint64_t hashes_offset;   // bytes from section start to this genus's hash array
    uint32_t n_hashes;
    uint32_t _pad;
};
static_assert(sizeof(FmhrEntry) == 24);

static constexpr uint32_t FMHR_EMPTY = UINT32_MAX;

// Zero-copy view into the mmap'd FMHR section for one genus.
struct FmhrView {
    uint64_t        genus_hash = 0;
    const uint64_t* hashes     = nullptr;
    uint32_t        n_hashes   = 0;
    bool valid() const noexcept { return hashes != nullptr && n_hashes > 0; }
    std::span<const uint64_t> span() const noexcept { return {hashes, n_hashes}; }
};

// ── Writer ────────────────────────────────────────────────────────────────────

class FmhrWriter {
public:
    // hashes must already be sorted and deduplicated.
    void add(uint64_t genus_hash, std::vector<uint64_t> hashes) {
        entries_.push_back({genus_hash, std::move(hashes)});
    }
    SectionDesc finalize(AppendWriter& w, uint64_t section_id);
    size_t n_genera() const { return entries_.size(); }

private:
    struct Entry { uint64_t genus_hash; std::vector<uint64_t> hashes; };
    std::vector<Entry> entries_;

    static uint32_t next_pow2(uint32_t v) noexcept {
        if (v == 0) return 1;
        --v;
        v |= v>>1; v |= v>>2; v |= v>>4; v |= v>>8; v |= v>>16;
        return v + 1;
    }
};

// ── Reader ────────────────────────────────────────────────────────────────────

class FmhrReader {
public:
    // data points to the start of the mmap'd file; offset is section start within it.
    void open(const uint8_t* data, uint64_t offset, uint64_t size) {
        if (size < sizeof(FmhrHeader))
            throw std::runtime_error("FMHR section too small");
        data_    = data + offset;
        header_  = reinterpret_cast<const FmhrHeader*>(data_);
        if (header_->magic != SEC_FMHR)
            throw std::runtime_error("FMHR: bad magic");
        entries_ = reinterpret_cast<const FmhrEntry*>(data_ + header_->entries_offset);
        buckets_ = reinterpret_cast<const uint32_t*>(data_ + header_->buckets_offset);
    }

    FmhrView lookup(uint64_t genus_hash) const noexcept {
        if (!header_ || header_->n_buckets == 0) return {};
        const uint32_t mask = header_->n_buckets - 1;
        uint32_t slot = static_cast<uint32_t>(genus_hash) & mask;
        for (uint32_t probe = 0; probe < header_->n_buckets; ++probe) {
            const uint32_t idx = buckets_[slot];
            if (idx == FMHR_EMPTY) return {};
            if (entries_[idx].genus_hash == genus_hash) {
                const auto& e = entries_[idx];
                return {genus_hash,
                        reinterpret_cast<const uint64_t*>(data_ + e.hashes_offset),
                        e.n_hashes};
            }
            slot = (slot + 1) & mask;
        }
        return {};
    }

    int k() const noexcept { return header_ ? static_cast<int>(header_->k) : 0; }
    int c() const noexcept { return header_ ? static_cast<int>(header_->c) : 0; }
    uint32_t n_genera() const noexcept { return header_ ? header_->n_genera : 0; }

private:
    const uint8_t*   data_    = nullptr;
    const FmhrHeader* header_  = nullptr;
    const FmhrEntry*  entries_ = nullptr;
    const uint32_t*   buckets_ = nullptr;
};

} // namespace genopack
