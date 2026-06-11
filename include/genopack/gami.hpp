#pragma once
// ── GAMI section: build-time per-aamer genus-multiplicity index ───────────────
// Stores a BlockedBloom filter over "rare" aamers (appearing in ≥1 but ≤K
// distinct genera across the whole archive). Written once by the builder after
// all genera are processed; loaded at check time via mmap — zero decode cost.
//
// Why Bloom, not the full GMI flat hash table:
//   The full GMI for a dense-PCORE archive can exceed 90 GB (5B unique aamers).
//   The Bloom filter over rare aamers is ~10× smaller: for K=3 and ~3B rare
//   aamers at 10 bits/element ≈ 3.75 GB — loadable as a single mmap at check time.
//   False positives (non-rare aamers classified as rare): conservative, not harmful.
#include <genopack/markers.hpp>   // BlockedBloom
#include <genopack/format.hpp>    // SEC_GAMI
#include <cstdint>
#include <vector>

namespace genopack {

// SEC_GAMI is defined in format.hpp
static constexpr uint64_t GAMI_MAGIC = 0x494D414700000001ULL;
static constexpr uint32_t GAMI_SIZE_CAP_GB = 20; // skip writing if Bloom > this

// ── On-disk section header (32 bytes) ────────────────────────────────────────
struct GamiHeader {
    uint64_t magic;      // GAMI_MAGIC
    uint32_t n_blocks;   // BlockedBloom n_blocks
    uint32_t K;          // rare threshold: count ∈ [1, K] → foreign-specific
    uint64_t n_rare;     // number of rare aamers inserted
    uint64_t frac_max;   // PCORE frac_max_hash the archive was built with
};
static_assert(sizeof(GamiHeader) == 32);

// ── Check-time view (mmap-backed, zero-copy) ──────────────────────────────────
struct GamiView {
    const uint64_t* bloom_data = nullptr;
    uint32_t        n_blocks   = 0;
    uint32_t        K          = 3;
    uint64_t        frac_max   = 0;

    bool valid() const noexcept { return bloom_data != nullptr && n_blocks > 0; }

    // Returns true if h is possibly a rare aamer (genus-count ∈ [1,K]).
    // False positives: aamers with count=0 may return true (novel organisms,
    // small fraction in practice). False negatives: never (Bloom guarantee).
    bool is_rare(uint64_t h) const noexcept {
        if (!bloom_data) return false;
        const uint32_t b = static_cast<uint32_t>(h % n_blocks) * 8;
        uint64_t mix = h;
        for (int i = 0; i < 4; ++i) {
            mix = mix * 0x9e3779b97f4a7c15ULL ^ (mix >> 32);
            if (!((bloom_data[b + ((mix >> 9) & 7)] >> (mix & 63)) & 1)) return false;
        }
        return true;
    }
};

// ── Build-time flat GMI (global multiplicity index) ───────────────────────────
// Open-addressing hash table: aamer_hash → genus_count (saturated at 255).
// Built incrementally during genus finalization; call maybe_resize() after
// each genus to stay below 50% load factor.
struct GlobalMultiplicityIndex {
    static constexpr uint64_t kEmpty = ~uint64_t(0);

    std::vector<uint64_t> keys_;
    std::vector<uint8_t>  vals_;
    uint64_t              mask_  = 0;
    size_t                count_ = 0;

    void reserve(size_t expected) {
        size_t cap = 1;
        while (cap < expected * 2) cap <<= 1;
        keys_.assign(cap, kEmpty);
        vals_.assign(cap, 0);
        mask_  = cap - 1;
        count_ = 0;
    }

    void increment(uint64_t h) noexcept {
        uint64_t i = h & mask_;
        while (keys_[i] != kEmpty && keys_[i] != h) i = (i + 1) & mask_;
        if (keys_[i] == kEmpty) { keys_[i] = h; vals_[i] = 1; ++count_; }
        else if (vals_[i] < 255) ++vals_[i];
    }

    void maybe_resize() {
        if (count_ * 2 <= keys_.size()) return;
        const size_t new_cap = keys_.size() * 4;
        std::vector<uint64_t> nk(new_cap, kEmpty);
        std::vector<uint8_t>  nv(new_cap, 0);
        const uint64_t new_mask = static_cast<uint64_t>(new_cap - 1);
        for (size_t i = 0; i < keys_.size(); ++i) {
            if (keys_[i] == kEmpty) continue;
            uint64_t j = keys_[i] & new_mask;
            while (nk[j] != kEmpty) j = (j + 1) & new_mask;
            nk[j] = keys_[i]; nv[j] = vals_[i];
        }
        keys_ = std::move(nk); vals_ = std::move(nv); mask_ = new_mask;
    }

    uint8_t lookup(uint64_t h) const noexcept {
        if (keys_.empty()) return 0;
        uint64_t i = h & mask_;
        while (keys_[i] != kEmpty && keys_[i] != h) i = (i + 1) & mask_;
        return keys_[i] == h ? vals_[i] : 0;
    }

    bool   empty() const noexcept { return count_ == 0; }
    size_t count() const noexcept { return count_; }
    size_t bytes() const noexcept { return keys_.size() * 9; }
};

} // namespace genopack
