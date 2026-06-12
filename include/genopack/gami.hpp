#pragma once
// ── GAMI section: precomputed per-aamer genus-multiplicity index ──────────────
// V1 (legacy, removed): BlockedBloom filter over rare aamers.
// V2 (current): exact sorted (hash, count8) pairs, zstd-compressed.
//   Eliminates the 10–30 min GMI rebuild cost at each `genopack check` run.
//   Built offline via `genopack gami build`; loaded in seconds at check start.
#include <genopack/format.hpp>    // SEC_GAMI, SectionDesc
#include <genopack/mmap_file.hpp> // AppendWriter
#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace genopack {

// SEC_GAMI is defined in format.hpp
static constexpr uint64_t GAMI_MAGIC_V1 = 0x494D414700000001ULL;  // legacy Bloom (unused)
static constexpr uint64_t GAMI_MAGIC_V2 = 0x494D414700000002ULL;  // exact sorted pairs

// ── V2 on-disk section header (64 bytes) ─────────────────────────────────────
// Immediately followed by a zstd-compressed payload of sorted
// (uint64_t hash, uint8_t count) pairs for all aamers with count ≥ 1.
struct GamiHeaderV2 {
    uint64_t magic          = GAMI_MAGIC_V2;
    uint64_t n_entries      = 0;     // total (hash, count) pairs stored
    uint64_t frac_max       = 0;     // PCORE frac_max_hash for consistency check
    uint64_t payload_bytes  = 0;     // compressed payload size in bytes
    uint64_t _pad[4]        = {};
};
static_assert(sizeof(GamiHeaderV2) == 64);

// ── Legacy V1 placeholder (kept for call-site compat in foreign_contam.hpp) ──
struct GamiView {
    bool valid() const noexcept { return false; }
    bool is_rare(uint64_t) const noexcept { return false; }
};

// ── Elias-Fano sorted integer set used by SEC_GAMI v2 load path ──────────────
// Stores sorted uint64_t hashes without changing the on-disk GAMI payload.
struct EliasFano {
    size_t n = 0;

    void build(const uint64_t* vals, size_t count, uint64_t universe) {
        n = count;
        universe_ = universe;
        l_ = 0;
        max_high_ = 0;
        high_bit_count_ = 0;
        low_.clear();
        high_.clear();
        high_rank_.clear();
        if (n == 0) return;

        const uint64_t ratio = universe_ / static_cast<uint64_t>(n);
        l_ = static_cast<uint8_t>(ratio == 0 ? 0 : floor_log2(ratio));
        if (l_ > 63) l_ = 63;

        const uint64_t mask = low_mask();
        const uint64_t low_bit_count = static_cast<uint64_t>(n) * l_;
        low_.assign(static_cast<size_t>((low_bit_count + 63) / 64), 0);

        max_high_ = vals[n - 1] >> l_;
        high_bit_count_ = max_high_ + static_cast<uint64_t>(n) + 1;
        high_.assign(static_cast<size_t>((high_bit_count_ + 63) / 64), 0);

        for (size_t i = 0; i < n; ++i) {
            write_low(i, vals[i] & mask);
            set_high_bit((vals[i] >> l_) + static_cast<uint64_t>(i));
        }

        high_rank_.assign(high_.size() + 1, 0);
        for (size_t i = 0; i < high_.size(); ++i)
            high_rank_[i + 1] = high_rank_[i] + std::popcount(high_[i]);
    }

    size_t find(uint64_t v) const noexcept {
        if (n == 0 || v >= universe_) return n;
        const uint64_t high = v >> l_;
        if (high > max_high_) return n;

        const uint64_t begin = high == 0 ? 0 : rank1_before(select0(high - 1));
        const uint64_t end = rank1_before(select0(high));
        const uint64_t low = v & low_mask();

        size_t lo = static_cast<size_t>(begin);
        size_t hi = static_cast<size_t>(end);
        while (lo < hi) {
            const size_t mid = lo + (hi - lo) / 2;
            if (read_low(mid) < low) lo = mid + 1;
            else hi = mid;
        }
        return lo < static_cast<size_t>(end) && read_low(lo) == low ? lo : n;
    }

    uint64_t access(size_t i) const noexcept {
        if (i >= n) return 0;
        const uint64_t pos = select1(static_cast<uint64_t>(i));
        const uint64_t high = pos - static_cast<uint64_t>(i);
        return (high << l_) | read_low(i);
    }

    size_t bytes() const noexcept {
        return low_.size() * sizeof(uint64_t)
             + high_.size() * sizeof(uint64_t)
             + high_rank_.size() * sizeof(uint64_t);
    }

private:
    uint64_t universe_ = 0;
    uint64_t max_high_ = 0;
    uint64_t high_bit_count_ = 0;
    uint8_t  l_ = 0;
    std::vector<uint64_t> low_;
    std::vector<uint64_t> high_;
    std::vector<uint64_t> high_rank_;

    static uint8_t floor_log2(uint64_t v) noexcept {
        uint8_t r = 0;
        while (v >>= 1) ++r;
        return r;
    }

    uint64_t low_mask() const noexcept {
        return l_ == 0 ? 0 : ((uint64_t{1} << l_) - 1);
    }

    void write_low(size_t i, uint64_t low) noexcept {
        if (l_ == 0) return;
        const uint64_t bit = static_cast<uint64_t>(i) * l_;
        const size_t word = static_cast<size_t>(bit / 64);
        const uint32_t shift = static_cast<uint32_t>(bit % 64);
        low_[word] |= low << shift;
        if (shift + l_ > 64)
            low_[word + 1] |= low >> (64 - shift);
    }

    uint64_t read_low(size_t i) const noexcept {
        if (l_ == 0) return 0;
        const uint64_t bit = static_cast<uint64_t>(i) * l_;
        const size_t word = static_cast<size_t>(bit / 64);
        const uint32_t shift = static_cast<uint32_t>(bit % 64);
        uint64_t low = low_[word] >> shift;
        if (shift + l_ > 64)
            low |= low_[word + 1] << (64 - shift);
        return low & low_mask();
    }

    void set_high_bit(uint64_t pos) noexcept {
        high_[static_cast<size_t>(pos / 64)] |= uint64_t{1} << (pos % 64);
    }

    uint64_t rank1_before(uint64_t pos) const noexcept {
        if (pos >= high_bit_count_) return static_cast<uint64_t>(n);
        const size_t word = static_cast<size_t>(pos / 64);
        const uint32_t bit = static_cast<uint32_t>(pos % 64);
        const uint64_t mask = bit == 0 ? 0 : ((uint64_t{1} << bit) - 1);
        return high_rank_[word] + std::popcount(high_[word] & mask);
    }

    uint64_t select1(uint64_t k) const noexcept {
        size_t lo = 0, hi = high_.size();
        while (lo < hi) {
            const size_t mid = lo + (hi - lo) / 2;
            if (high_rank_[mid + 1] > k) hi = mid;
            else lo = mid + 1;
        }
        uint64_t word = high_[lo];
        uint64_t in_word = k - high_rank_[lo];
        while (in_word-- > 0) word &= word - 1;
        return static_cast<uint64_t>(lo) * 64 + std::countr_zero(word);
    }

    uint64_t select0(uint64_t k) const noexcept {
        size_t lo = 0, hi = high_.size();
        while (lo < hi) {
            const size_t mid = lo + (hi - lo) / 2;
            const uint64_t valid_bits_end =
                std::min<uint64_t>((static_cast<uint64_t>(mid) + 1) * 64, high_bit_count_);
            const uint64_t zeros_end = valid_bits_end - high_rank_[mid + 1];
            if (zeros_end > k) hi = mid;
            else lo = mid + 1;
        }

        const uint64_t valid_bits_before =
            std::min<uint64_t>(static_cast<uint64_t>(lo) * 64, high_bit_count_);
        const uint64_t zeros_before = valid_bits_before - high_rank_[lo];
        uint64_t word = high_[lo];
        if (lo + 1 == high_.size() && (high_bit_count_ % 64) != 0) {
            const uint32_t valid = static_cast<uint32_t>(high_bit_count_ % 64);
            word |= ~((uint64_t{1} << valid) - 1);
        }
        word = ~word;
        uint64_t in_word = k - zeros_before;
        while (in_word-- > 0) word &= word - 1;
        return static_cast<uint64_t>(lo) * 64 + std::countr_zero(word);
    }
};

// ── Build-time flat GMI (global multiplicity index) ───────────────────────────
// Open-addressing hash table: aamer_hash → genus_count (saturated at 255).
// Built incrementally during genus finalization; call maybe_resize() after
// each genus to stay below 50% load factor.
struct GlobalMultiplicityIndex {
    static constexpr uint64_t kEmpty = ~uint64_t(0);

    // Hash-table storage: used by the runtime GMI build path (increment/maybe_resize).
    std::vector<uint64_t> keys_;
    std::vector<uint8_t>  vals_;
    uint64_t              mask_  = 0;
    size_t                count_ = 0;

    // Elias-Fano key storage plus flat counts: used by the SEC_GAMI load path.
    std::vector<uint8_t>  sorted_vals_;
    EliasFano             ef_;

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
        if (count_ * 10 <= keys_.size() * 7) return;  // trigger at 70% load
        const size_t new_cap = keys_.size() * 2;       // grow 2× (was 4× at 50%)
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
        if (ef_.n != 0) {
            const size_t i = ef_.find(h);
            return i < ef_.n ? sorted_vals_[i] : 0;
        }
        if (keys_.empty()) return 0;
        uint64_t i = h & mask_;
        while (keys_[i] != kEmpty && keys_[i] != h) i = (i + 1) & mask_;
        return keys_[i] == h ? vals_[i] : 0;
    }

    bool   empty() const noexcept { return count_ == 0 && ef_.n == 0; }
    size_t count() const noexcept { return count_ + ef_.n; }
    size_t bytes() const noexcept {
        return keys_.size() * 9
             + ef_.bytes() + sorted_vals_.size();
    }
};

// ── V2 writer: serializes a GlobalMultiplicityIndex → SEC_GAMI v2 ────────────
// Call finalize() once; the writer sorts all (hash, count) pairs, zstd-compresses
// them, and appends a SEC_GAMI section to the archive.
struct GamiV2Writer {
    SectionDesc finalize(AppendWriter& w, uint64_t section_id,
                         const GlobalMultiplicityIndex& gmi, uint64_t frac_max);
};

// ── V2 loader: decompresses SEC_GAMI v2 payload → GlobalMultiplicityIndex ────
// section_data must point to the start of the SEC_GAMI section (GamiHeaderV2).
// section_size is the number of bytes in the section.
void gami_v2_load(const uint8_t* section_data, uint64_t section_size,
                  GlobalMultiplicityIndex& out);

} // namespace genopack
