#pragma once
// ── SEC_PCORE — unified prevalence-annotated per-genus aamer reference ─────────
// One per-genus reference that supersedes CORE (prevalence core) and FCORE: it
// stores EVERY aamer seen in ≥1 member (the dense union) together with a u8
// prevalence (percent of members carrying it). From this single section a consumer
// derives, at query time:
//   prevalence ≥ θ%  → the conserved completeness core  (== old CORE; a filter)
//   prevalence ≥ 1   → the DENSE foreign reference       (small-contig detection)
//   prevalence       → the per-aamer weight for a real log-likelihood ratio
// Layout mirrors CORE (open-addressing hash over key_hash → sorted u64 aamers),
// plus a parallel u8 prevalence array per entry. Keyed by genus hash (a family
// tier may be added with the same key space).
#include "format.hpp"
#include "mmap_file.hpp"
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace genopack {

struct PcoreHeader {             // 64 bytes
    uint32_t magic;              //  0  SEC_PCORE
    uint32_t n_entries;          //  4
    uint32_t n_buckets;          //  8  power of 2
    uint32_t k;                  // 12  aamer k
    uint32_t min_seg_aa;         // 16
    float    theta;              // 20  completeness-core prevalence threshold (== old CORE theta)
    uint64_t frac_max_hash;      // 24
    uint64_t model_hash;         // 32
    uint64_t entries_offset;     // 40  -> PcoreEntry[n_entries]
    uint64_t buckets_offset;     // 48  -> uint32_t[n_buckets]
    uint64_t _reserved;          // 56
};
static_assert(sizeof(PcoreHeader) == 64);

struct PcoreEntry {              // 32 bytes
    uint64_t key_hash;           //  0  genus (or family) hash
    uint64_t aamers_offset;      //  8  -> sorted u64 aamers (8-aligned)
    uint64_t prev_offset;        // 16  -> u8 prevalence percent[n_aamers]
    uint32_t n_aamers;           // 24
    uint32_t n_members;          // 28
};
static_assert(sizeof(PcoreEntry) == 32);

static constexpr uint32_t PCORE_EMPTY = UINT32_MAX;

struct PcoreView {
    uint64_t        key_hash  = 0;
    const uint64_t* aamers    = nullptr;   // sorted
    const uint8_t*  prev      = nullptr;   // member COUNT carrying each aamer (cap 255), parallel
    uint32_t        n_aamers  = 0;
    uint32_t        n_members = 0;
    bool valid() const noexcept { return aamers != nullptr && n_aamers > 0; }
};

// ⌈theta · n_members⌉ — the member-count threshold for the conserved completeness
// core (identical rule to CORE), so the prev>=this slice reproduces CORE bit-for-bit.
inline uint32_t pcore_core_threshold(uint32_t n_members, float theta) {
    return static_cast<uint32_t>(std::ceil(static_cast<double>(theta) * n_members));
}

// Completeness coverage = |query ∩ {aamers with count >= ⌈theta·n⌉}| / |{...}|.
// query_sorted must be sorted-unique. Returns NaN if the core slice is empty.
inline float pcore_core_coverage(const PcoreView& v, const std::vector<uint64_t>& query_sorted,
                                 float theta) {
    if (!v.valid() || v.n_members == 0) return std::numeric_limits<float>::quiet_NaN();
    const uint32_t thr = pcore_core_threshold(v.n_members, theta);
    uint32_t core_size = 0, inter = 0;
    size_t qi = 0;
    for (uint32_t i = 0; i < v.n_aamers; ++i) {
        if (v.prev[i] < thr) continue;          // not in the conserved core
        ++core_size;
        while (qi < query_sorted.size() && query_sorted[qi] < v.aamers[i]) ++qi;
        if (qi < query_sorted.size() && query_sorted[qi] == v.aamers[i]) ++inter;
    }
    if (core_size == 0) return std::numeric_limits<float>::quiet_NaN();
    return static_cast<float>(inter) / static_cast<float>(core_size);
}

// ── Writer ────────────────────────────────────────────────────────────────────
class PcoreWriter {
public:
    PcoreWriter(uint32_t k, uint32_t min_seg_aa, float theta = 0.90f)
        : k_(k), min_seg_aa_(min_seg_aa), theta_(theta) {}

    // member_qmers = per-member sorted-unique aamer sets. Stores every aamer present
    // in ≥1 member with prev[i] = member count carrying it (capped at 255). The
    // prev>=⌈theta·n⌉ slice is the conserved completeness core (== old CORE).
    void add_from_members(uint64_t key_hash,
                          const std::vector<std::vector<uint64_t>>& member_qmers);

    uint64_t model_hash() const;
    SectionDesc finalize(AppendWriter& w, uint64_t section_id, uint64_t frac_max_hash);
    size_t n_entries() const { return entries_.size(); }

private:
    struct Entry {
        uint64_t key_hash;
        std::vector<uint64_t> aamers;   // sorted
        std::vector<uint8_t>  prev;     // prevalence percent, parallel
        uint32_t n_members;
    };
    static uint32_t next_pow2(uint32_t v) noexcept {
        if (v == 0) return 1;
        --v; v |= v>>1; v |= v>>2; v |= v>>4; v |= v>>8; v |= v>>16; return v + 1;
    }
    uint32_t k_, min_seg_aa_;
    float    theta_;
    std::vector<Entry> entries_;
};

// ── Reader ────────────────────────────────────────────────────────────────────
class PcoreReader {
public:
    void open(const uint8_t* data, uint64_t offset, uint64_t size) {
        if (size < sizeof(PcoreHeader)) throw std::runtime_error("PCORE section too small");
        data_   = data + offset;
        header_ = reinterpret_cast<const PcoreHeader*>(data_);
        if (header_->magic != SEC_PCORE) throw std::runtime_error("PCORE: bad magic");
        const uint64_t ent_end = header_->entries_offset
            + static_cast<uint64_t>(header_->n_entries) * sizeof(PcoreEntry);
        const uint64_t buk_end = header_->buckets_offset
            + static_cast<uint64_t>(header_->n_buckets) * sizeof(uint32_t);
        if (ent_end > size || buk_end > size) throw std::runtime_error("PCORE section truncated");
        entries_ = reinterpret_cast<const PcoreEntry*>(data_ + header_->entries_offset);
        buckets_ = reinterpret_cast<const uint32_t*>(data_ + header_->buckets_offset);
    }

    PcoreView lookup(uint64_t key_hash) const noexcept {
        if (!header_ || header_->n_buckets == 0) return {};
        const uint32_t mask = header_->n_buckets - 1;
        uint32_t slot = static_cast<uint32_t>(key_hash) & mask;
        for (uint32_t probe = 0; probe < header_->n_buckets; ++probe) {
            const uint32_t idx = buckets_[slot];
            if (idx == PCORE_EMPTY) return {};
            if (entries_[idx].key_hash == key_hash) {
                const auto& e = entries_[idx];
                return {key_hash,
                        reinterpret_cast<const uint64_t*>(data_ + e.aamers_offset),
                        reinterpret_cast<const uint8_t*>(data_ + e.prev_offset),
                        e.n_aamers, e.n_members};
            }
            slot = (slot + 1) & mask;
        }
        return {};
    }

    bool     is_open() const noexcept { return data_ != nullptr; }
    uint32_t n_entries() const noexcept { return header_ ? header_->n_entries : 0; }
    uint64_t key_hash_at(uint32_t i) const noexcept {
        return (header_ && i < header_->n_entries) ? entries_[i].key_hash : 0;
    }
    int      k() const noexcept { return header_ ? static_cast<int>(header_->k) : 0; }
    int      min_seg_aa() const noexcept { return header_ ? static_cast<int>(header_->min_seg_aa) : 0; }
    float    theta() const noexcept { return header_ ? header_->theta : 0.90f; }
    uint64_t frac_max_hash() const noexcept { return header_ ? header_->frac_max_hash : 0; }
    uint64_t model_hash() const noexcept { return header_ ? header_->model_hash : 0; }

private:
    const uint8_t*     data_    = nullptr;
    const PcoreHeader* header_  = nullptr;
    const PcoreEntry*  entries_ = nullptr;
    const uint32_t*    buckets_ = nullptr;
};

} // namespace genopack
