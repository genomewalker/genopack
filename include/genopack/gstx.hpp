#pragma once
#include "format.hpp"
#include "mmap_file.hpp"
#include "types.hpp"
#include <cstdint>
#include <functional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace genopack {

// ── Constants ─────────────────────────────────────────────────────────────────

static constexpr uint32_t GSTX_MAX_K    = 3;    // max k-values stored per genus
static constexpr uint32_t GSTX_BINS     = 10000; // OPH bins per sketch

// ── On-disk layout ────────────────────────────────────────────────────────────

struct GstxHeader {         // 72 bytes, at section start
    uint32_t magic;         // SEC_GSTX
    uint32_t n_genera;
    uint32_t n_buckets;     // open-addressing table size, power of 2
    uint32_t n_k;           // k-values stored (≤ GSTX_MAX_K)
    uint32_t sketch_size;   // bins per sketch (must equal GSTX_BINS)
    uint32_t kmer_sizes[4]; // k values e.g. [16,21,31,0]
    uint32_t _pad0;         // pad to 8-byte align entries_offset
    uint64_t entries_offset;// bytes from section start → GstxEntry array
    uint64_t buckets_offset;// bytes from section start → uint32_t bucket array
    uint64_t entry_stride;  // sizeof(GstxEntry) for forward-compat check
    uint8_t  reserved[4];
    uint32_t _pad1;
};
static_assert(sizeof(GstxHeader) == 72);

// Fixed-stride entry; the open-addressing table stores uint32_t indices into this array.
// genus_hash == 0 means empty slot (genus "" has hash remapped to 1 to avoid clash).
struct GstxEntry {          // 60576 bytes
    uint64_t genus_hash;    // FNV-1a of genus string; 0 = empty slot
    uint32_t n_members;
    uint16_t flags;         // bit 0: has_tnf
    uint8_t  n_k_stored;    // actual k-values in consensus[] (≤ GSTX_MAX_K)
    uint8_t  _pad0;
    float    p90_containment[GSTX_MAX_K]; // [k0,k1,k2]
    float    nrb_p90;                     // p90 of n_real_bins at k=0 across genus members
    float    tnf_mu[136];                 // L2-normalised TNF centroid
    uint16_t consensus[GSTX_MAX_K][GSTX_BINS]; // majority-vote sig per bin
};
static_assert(sizeof(GstxEntry) ==
    8 + 4 + 2 + 1 + 1 +
    GSTX_MAX_K * 4 + 4 +
    136 * 4 +
    GSTX_MAX_K * GSTX_BINS * 2);

// ── Writer ────────────────────────────────────────────────────────────────────

class GstxWriter {
public:
    // Call once per genus after flush_staging_buf finalises its consensus/stats.
    // consensus[ki][bin] = majority uint16 value for bin at k-index ki.
    // p90[ki] = p90 containment fraction at k-index ki.
    // tnf_mu[136] = per-genus TNF centroid (may be nullptr if n_members < 2).
    void add_genus(std::string_view genus,
                   uint32_t         n_members,
                   uint32_t         n_k,
                   const std::vector<std::vector<uint16_t>>& consensus, // [n_k][GSTX_BINS]
                   const float*     p90,          // float[n_k]
                   const float*     tnf_mu,       // float[136] or nullptr
                   const uint32_t*  kmer_sizes,   // uint32_t[n_k]
                   float            nrb_p90 = 0.0f);

    // Write GSTX section to w; returns the SectionDesc to add to the TOC.
    SectionDesc finalize(AppendWriter& w, uint64_t section_id);

    size_t n_genera() const { return entries_.size(); }

private:
    static uint64_t hash_genus(std::string_view s) noexcept;
    static uint32_t next_pow2(uint32_t v) noexcept;

    std::vector<GstxEntry>   entries_;
    std::vector<uint32_t>    kmer_sizes_hdr_{0, 0, 0, 0};
    uint32_t                 n_k_hdr_ = 0;
};

// ── Reader ────────────────────────────────────────────────────────────────────

class GstxReader {
public:
    void open(const uint8_t* data, uint64_t offset, uint64_t size);

    bool is_open()   const { return data_ != nullptr; }
    uint32_t n_genera()   const { return header_ ? header_->n_genera  : 0; }
    uint32_t n_k()        const { return header_ ? header_->n_k       : 0; }
    uint32_t sketch_size()const { return header_ ? header_->sketch_size: 0; }
    uint32_t kmer_size(int ki) const {
        return (header_ && ki < 4) ? header_->kmer_sizes[ki] : 0;
    }

    // Returns nullptr if genus not found or reader not open.
    const GstxEntry* lookup(std::string_view genus) const;

    // Iterate all entries (valid and empty slots filtered out).
    void scan(const std::function<void(const GstxEntry&)>& cb) const {
        if (!data_) return;
        for (uint32_t i = 0; i < header_->n_genera; ++i)
            if (entries_[i].genus_hash != 0) cb(entries_[i]);
    }

private:
    static uint64_t hash_genus(std::string_view s) noexcept;

    const uint8_t*    data_      = nullptr;
    const GstxHeader* header_    = nullptr;
    const GstxEntry*  entries_   = nullptr;
    const uint32_t*   buckets_   = nullptr;
    uint32_t          n_buckets_ = 0;
};

} // namespace genopack
