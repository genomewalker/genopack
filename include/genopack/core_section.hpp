#pragma once
// ── SEC_CORE — per-genus prevalence cores ─────────────────────────────────────
// A genus prevalence core is the set of aamers present in >= ceil(theta * N) of a
// genus's N reference members. Genus-core coverage = |genome ∩ core| / |core| is
// the intrinsic completeness signal (validated: RMSE ~5.1% vs known level, ties
// CheckM2). This section makes the per-genus cores a FIRST-CLASS, content-hashed
// model so the intrinsic completeness columns are a pure function of build inputs
// (they reference core_model_hash, not the raw sibling set). Structurally this is
// FMHR (per-genus sorted u64 sets + open-addressing bucket table + concatenated
// pool) plus prevalence params, contributing member ids, and the model hash.
#include "format.hpp"
#include "mmap_file.hpp"
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace genopack {

struct CoreHeader {              // 64 bytes
    uint32_t magic;              //  0  SEC_CORE
    uint32_t n_genera;           //  4
    uint32_t n_buckets;          //  8  open-addressing table, power of 2
    uint32_t k;                  // 12  aamer k (AAMER_K = 8)
    uint32_t min_seg_aa;         // 16  min inter-stop AA segment at extraction
    float    theta;              // 20  prevalence threshold (e.g. 0.90)
    uint64_t frac_max_hash;      // 24  FracMinHash max hash at extraction (UINT64_MAX = none)
    uint64_t core_model_hash;    // 32  content hash of (params + per-genus cores); referenced by QCOL
    uint64_t entries_offset;     // 40  -> CoreEntry[n_genera]
    uint64_t buckets_offset;     // 48  -> uint32_t[n_buckets]
    uint64_t _reserved;          // 56
};
static_assert(sizeof(CoreHeader) == 64);

struct CoreEntry {               // 32 bytes
    uint64_t genus_hash;         //  0
    uint64_t aamers_offset;      //  8  bytes from section start -> sorted u64 core aamers
    uint64_t ids_offset;         // 16  bytes from section start -> sorted u64 contributing genome ids
    uint32_t n_aamers;           // 24
    uint32_t n_members;          // 28  count of contributing genome ids
};
static_assert(sizeof(CoreEntry) == 32);

static constexpr uint32_t CORE_EMPTY = UINT32_MAX;

// Zero-copy view of one genus's core.
struct CoreView {
    uint64_t        genus_hash = 0;
    const uint64_t* aamers     = nullptr;
    uint32_t        n_aamers   = 0;
    const uint64_t* members    = nullptr;
    uint32_t        n_members  = 0;
    bool valid() const noexcept { return aamers != nullptr && n_aamers > 0; }
    std::span<const uint64_t> aamer_span() const noexcept { return {aamers, n_aamers}; }
};

// ── Writer ────────────────────────────────────────────────────────────────────
class CoreWriter {
public:
    // section_type is SEC_CORE (per-genus, default) or SEC_FCORE (per-family). It
    // only changes the on-disk magic + TOC type; the layout is identical, so one
    // CoreReader reads either (mirrors GcovWriter emitting GCOV and FCOV).
    CoreWriter(uint32_t k, uint32_t min_seg_aa, float theta, uint32_t section_type = SEC_CORE)
        : k_(k), min_seg_aa_(min_seg_aa), theta_(theta), section_type_(section_type) {}

    // Build a prevalence core from per-member sorted-unique aamer sets and add it.
    // Keeps aamers present in >= ceil(theta * n_members) members. No-op if the
    // resulting core is empty. member_ids is stored (sorted) for provenance.
    void add_from_members(uint64_t genus_hash,
                          const std::vector<std::vector<uint64_t>>& member_qmers,
                          std::vector<uint64_t> member_ids);

    // Direct add of an already-computed sorted-unique core.
    void add(uint64_t genus_hash, std::vector<uint64_t> core_aamers,
             std::vector<uint64_t> member_ids);

    // Deterministic content hash of (k, min_seg_aa, theta, per-genus cores sorted
    // by genus_hash). Independent of insertion order and of member_ids. Valid any
    // time after the adds; QCOL columns reference this value.
    uint64_t core_model_hash() const;

    SectionDesc finalize(AppendWriter& w, uint64_t section_id, uint64_t frac_max_hash);

    size_t n_genera() const { return entries_.size(); }

private:
    struct Entry {
        uint64_t genus_hash;
        std::vector<uint64_t> aamers;   // sorted-unique core
        std::vector<uint64_t> members;  // sorted contributing genome ids
    };
    uint32_t k_, min_seg_aa_;
    float    theta_;
    uint32_t section_type_ = SEC_CORE;
    std::vector<Entry> entries_;

    static uint32_t next_pow2(uint32_t v) noexcept {
        if (v == 0) return 1;
        --v; v |= v>>1; v |= v>>2; v |= v>>4; v |= v>>8; v |= v>>16;
        return v + 1;
    }
};

// ── Reader ────────────────────────────────────────────────────────────────────
class CoreReader {
public:
    void open(const uint8_t* data, uint64_t offset, uint64_t size) {
        if (size < sizeof(CoreHeader))
            throw std::runtime_error("CORE section too small");
        data_   = data + offset;
        header_ = reinterpret_cast<const CoreHeader*>(data_);
        if (header_->magic != SEC_CORE && header_->magic != SEC_FCORE)
            throw std::runtime_error("CORE: bad magic");
        const uint64_t ent_end = header_->entries_offset
            + static_cast<uint64_t>(header_->n_genera) * sizeof(CoreEntry);
        const uint64_t buk_end = header_->buckets_offset
            + static_cast<uint64_t>(header_->n_buckets) * sizeof(uint32_t);
        if (ent_end > size || buk_end > size)
            throw std::runtime_error("CORE section truncated");
        entries_ = reinterpret_cast<const CoreEntry*>(data_ + header_->entries_offset);
        buckets_ = reinterpret_cast<const uint32_t*>(data_ + header_->buckets_offset);
        size_    = size;
    }

    CoreView lookup(uint64_t genus_hash) const noexcept {
        if (!header_ || header_->n_buckets == 0) return {};
        const uint32_t mask = header_->n_buckets - 1;
        uint32_t slot = static_cast<uint32_t>(genus_hash) & mask;
        for (uint32_t probe = 0; probe < header_->n_buckets; ++probe) {
            const uint32_t idx = buckets_[slot];
            if (idx == CORE_EMPTY) return {};
            if (entries_[idx].genus_hash == genus_hash) {
                const auto& e = entries_[idx];
                return {genus_hash,
                        reinterpret_cast<const uint64_t*>(data_ + e.aamers_offset), e.n_aamers,
                        reinterpret_cast<const uint64_t*>(data_ + e.ids_offset),    e.n_members};
            }
            slot = (slot + 1) & mask;
        }
        return {};
    }

    bool     is_open() const noexcept { return data_ != nullptr; }
    uint32_t n_genera() const noexcept { return header_ ? header_->n_genera : 0; }
    // Enumerate the stored genus hashes (for deterministic background-panel sampling).
    uint64_t genus_hash_at(uint32_t i) const noexcept {
        return (header_ && i < header_->n_genera) ? entries_[i].genus_hash : 0;
    }
    int      k() const noexcept { return header_ ? static_cast<int>(header_->k) : 0; }
    int      min_seg_aa() const noexcept { return header_ ? static_cast<int>(header_->min_seg_aa) : 0; }
    float    theta() const noexcept { return header_ ? header_->theta : 0.0f; }
    uint64_t frac_max_hash() const noexcept { return header_ ? header_->frac_max_hash : 0; }
    uint64_t core_model_hash() const noexcept { return header_ ? header_->core_model_hash : 0; }

private:
    const uint8_t*    data_    = nullptr;
    const CoreHeader* header_  = nullptr;
    const CoreEntry*  entries_ = nullptr;
    const uint32_t*   buckets_ = nullptr;
    uint64_t          size_    = 0;
};

} // namespace genopack
