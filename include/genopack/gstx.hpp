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
// Members with k=0 containment < GSTX_P90_COMPLETE_FRAC * cluster_max are excluded
// from p90 computation — this prevents incomplete genomes from deflating the p90
// reference and causing systematic CCR underestimates at high completeness.
static constexpr float    GSTX_P90_COMPLETE_FRAC = 0.5f;

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
//
// FORMAT NOTE (dispersion bump): median_containment/mad_containment are APPENDED at the
// end of the entry. Every field before them keeps its byte offset, so an old-layout entry
// (written before the bump) reinterpreted through this struct is valid for all prefix
// fields; only median/mad are absent. The GstxHeader::entry_stride discriminates the two
// layouts (GSTX_STRIDE_LEGACY vs sizeof(GstxEntry)); the reader indexes entries by the
// on-disk stride and exposes has_dispersion() so callers never read median/mad from an
// old pack. This mirrors the .csp v2→v3 dual-magic reader.
struct GstxEntry {
    uint64_t genus_hash;    // FNV-1a of genus string; 0 = empty slot
    uint32_t n_members;
    uint16_t flags;         // bit 0: has_tnf
    uint8_t  n_k_stored;    // actual k-values in consensus[] (≤ GSTX_MAX_K)
    uint8_t  _pad0;
    float    p90_containment[GSTX_MAX_K]; // [k0,k1,k2]
    float    nrb_p90;                     // p90 of n_real_bins at k=0 across genus members
    float    tnf_mu[136];                 // L2-normalised TNF centroid
    uint16_t consensus[GSTX_MAX_K][GSTX_BINS]; // majority-vote sig per bin
    // ── dispersion bump (appended; absent in GSTX_STRIDE_LEGACY packs) ──
    float    median_containment[GSTX_MAX_K]; // median member OPH-containment per k
    float    mad_containment[GSTX_MAX_K];    // 1.4826·MAD of member OPH-containment per k
};

// Byte size of the entry as written before the dispersion bump (median/mad absent).
static constexpr uint64_t GSTX_STRIDE_LEGACY =
    8 + 4 + 2 + 1 + 1 +
    GSTX_MAX_K * 4 + 4 +
    136 * 4 +
    static_cast<uint64_t>(GSTX_MAX_K) * GSTX_BINS * 2;
static_assert(GSTX_STRIDE_LEGACY == 60576);
static_assert(sizeof(GstxEntry) ==
    GSTX_STRIDE_LEGACY + GSTX_MAX_K * 4 + GSTX_MAX_K * 4);

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
                   float            nrb_p90 = 0.0f,
                   const float*     median = nullptr, // float[n_k] or nullptr → NaN
                   const float*     mad    = nullptr);// float[n_k] or nullptr → NaN

    // Write GSTX section to w; returns the SectionDesc to add to the TOC.
    SectionDesc finalize(AppendWriter& w, uint64_t section_id);

    size_t n_genera() const { return entries_.size(); }

private:
    static uint64_t hash_genus(std::string_view s) noexcept;
    static uint32_t next_pow2(uint32_t v) noexcept;

    std::vector<GstxEntry>   entries_;
    std::vector<uint32_t>    entry_order_; // arrival indices — sort by these in finalize() for reproducibility
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

    // True when the pack was written with the dispersion bump (median/mad present).
    // Old packs (GSTX_STRIDE_LEGACY) return false → median/mad are absent; callers must
    // not read GstxEntry::median_containment/mad_containment from them.
    bool has_dispersion() const { return has_dispersion_; }

    // Returns nullptr if genus not found or reader not open.
    const GstxEntry* lookup(std::string_view genus) const;

    // Iterate all entries (valid and empty slots filtered out).
    void scan(const std::function<void(const GstxEntry&)>& cb) const {
        if (!data_) return;
        for (uint32_t i = 0; i < header_->n_genera; ++i) {
            const GstxEntry& e = entry_at(i);
            if (e.genus_hash != 0) cb(e);
        }
    }

private:
    static uint64_t hash_genus(std::string_view s) noexcept;

    // Stride-aware entry access: on old packs entry_stride_ < sizeof(GstxEntry), so we
    // index by the on-disk stride. All fields up to consensus[] share offsets across
    // layouts; median/mad are only valid when has_dispersion_ is true.
    const GstxEntry& entry_at(uint32_t i) const {
        return *reinterpret_cast<const GstxEntry*>(
            entries_base_ + static_cast<uint64_t>(i) * entry_stride_);
    }

    const uint8_t*    data_          = nullptr;
    const GstxHeader* header_        = nullptr;
    const uint8_t*    entries_base_  = nullptr;
    uint64_t          entry_stride_  = 0;
    bool              has_dispersion_ = false;
    const uint32_t*   buckets_       = nullptr;
    uint32_t          n_buckets_     = 0;
};

} // namespace genopack
