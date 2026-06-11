#pragma once
// ── GAMI section: precomputed per-aamer genus-multiplicity index ──────────────
// V1 (legacy, removed): BlockedBloom filter over rare aamers.
// V2 (current): exact sorted (hash, count8) pairs, zstd-compressed.
//   Eliminates the 10–30 min GMI rebuild cost at each `genopack check` run.
//   Built offline via `genopack gami build`; loaded in seconds at check start.
#include <genopack/format.hpp>    // SEC_GAMI, SectionDesc
#include <genopack/mmap_file.hpp> // AppendWriter
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
