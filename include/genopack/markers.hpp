#pragma once
#include "aamer.hpp"
#include "mmap_file.hpp"
#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <numeric>
#include <queue>
#include <span>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace genopack {

// SEC_MRKR sidecar (.mrk file): global deduplicated marker aamer pool +
// per-lineage calibration. Built once per GTDB release by `genopack markers build`.
// Read at check time via mmap; the pool (sorted uint64 arrays) supports O(log n)
// membership test; a per-genome in-RAM hash map is built for pass-B scoring.

static constexpr uint32_t MRKR_MAGIC       = 0x504B524Du; // "MRKP"
static constexpr uint8_t  MRKR_VERSION     = 2;
static constexpr uint8_t  MRKR_VERSION_MIN = 1; // lowest readable version
static constexpr uint8_t  MRKR_DOMAIN_BAC = 0;
static constexpr uint8_t  MRKR_DOMAIN_ARC = 1;
static constexpr uint8_t  MRKR_ALPHABET_FULL_AA   = 0;
static constexpr uint8_t  MRKR_ALPHABET_MURPHY10  = 1; // Murphy 10-class reduced alphabet
static constexpr uint8_t  MRKR_ALPHABET_DAYHOFF6  = 2; // Dayhoff-6 groups, k=12 syncmers

// ── On-disk structs (little-endian POD) ──────────────────────────────────────

struct MarkerHeader {         // 40 bytes (v2+; v1 was 32 bytes)
    uint32_t magic;           // MRKR_MAGIC
    uint16_t version;         // MRKR_VERSION
    uint8_t  k;               // aamer k
    uint8_t  alphabet;        // MRKR_ALPHABET_*
    uint32_t n_lineages;      // entries in lookup table
    uint8_t  n_bac_markers;   // 120 for GTDB bac120
    uint8_t  n_arc_markers;   // 53 for GTDB ar53
    uint16_t frac_scale;      // FracMinHash divisor (0 or 1 = keep all hashes)
    uint32_t lookup_off;      // byte offset → MarkerLookupEntry[n_lineages]
    uint32_t pool_idx_off;    // byte offset → MarkerPoolEntry[n_bac+n_arc]
    uint32_t calib_off;       // byte offset → packed CalibEntry data
    uint32_t pool_off;        // byte offset → concatenated sorted uint64 hash arrays
    // Merged-pool fast-path: pre-sorted (hash,id) pairs written at build time.
    // If non-zero, build_merged_pool() uses mmap pointers instead of k-way merge.
    // Layout: N_bac uint64 hashes | N_bac uint8 ids | (arc equivalent at merged_arc_off).
    // N_bac/N_arc are the respective sums of per-marker pool sizes.
    uint32_t merged_bac_off;  // byte offset → pre-merged bac pool (0 if absent)
    uint32_t merged_arc_off;  // byte offset → pre-merged arc pool (0 if absent)
};
static_assert(sizeof(MarkerHeader) == 40);

// Sorted by lineage_hash for O(log n) binary search over ~20k lineages.
struct MarkerLookupEntry {    // 16 bytes
    uint64_t lineage_hash;
    uint32_t calib_data_off;  // byte offset from calib_off to CalibEntry for this lineage
    uint8_t  domain;          // MRKR_DOMAIN_BAC or MRKR_DOMAIN_ARC
    uint8_t  _pad[3];
};
static_assert(sizeof(MarkerLookupEntry) == 16);

// One per marker family (bac markers 0..n_bac-1, then arc markers 0..n_arc-1).
struct MarkerPoolEntry {      // 16 bytes
    uint64_t hashes_off;      // byte offset from pool_off to first uint64 hash
    uint32_t n_hashes;
    uint32_t _pad;
};
static_assert(sizeof(MarkerPoolEntry) == 16);

// Variable-size on-disk per-lineage calibration block, located at calib_off + calib_data_off.
// Layout: CalibEntryHeader (8 bytes)
//       + CalibSlot[n_markers] (4 bytes each)
//       + expected_bitset (ceil(n_markers/8) bytes, rounded to 8)
struct CalibEntryHeader {     // 8 bytes
    uint16_t n_markers;       // 120 (bac) or 53 (arc)
    uint16_t ref_count;       // number of reference genomes used
    uint8_t  domain;
    uint8_t  flags;
    uint16_t _pad;
};
static_assert(sizeof(CalibEntryHeader) == 8);

struct CalibSlot {            // 4 bytes per marker
    uint16_t null_floor_u16;  // background containment fraction × 65535
    uint16_t threshold_u16;   // presence call threshold × 65535
};
static_assert(sizeof(CalibSlot) == 4);

// Per-genus redundancy calibration entry (v2+).
// Sorted by genus_hash for O(log n) lookup.
// When n has REDUN_BORROWED_FLAG set, median/mad are from the family-level distribution.
static constexpr uint16_t REDUN_BORROWED_FLAG = 0x8000u;
struct RedunCalibEntry {      // 24 bytes
    uint64_t genus_hash;      //  8  FNV-1a of "g__Name" (same key as gcov/markers lookup)
    float    median;          //  4  median marker_redundancy of reference genomes in this genus
    float    mad;             //  4  median absolute deviation (floored at max(0.01*median, 1e-4))
    uint16_t n;               //  2  reference count; high bit = REDUN_BORROWED_FLAG (family fallback)
    uint8_t  _pad[6];         //  6
};
static_assert(sizeof(RedunCalibEntry) == 24);

// ── In-RAM acceleration structures ────────────────────────────────────────────

// Flat open-addressing hash map: uint64 hash → uint8 marker_id.
// ~220k entries per genus, L2/L3-resident. Much faster than the 12M merged pool.
struct FlatHitMap {
    static constexpr uint64_t EMPTY = UINT64_MAX;
    std::vector<uint64_t> keys;
    std::vector<uint8_t>  vals;
    uint32_t              mask = 0;
    uint32_t              count = 0;

    void build(const std::vector<std::pair<uint64_t,uint8_t>>& entries) {
        uint32_t cap = 1;
        while (cap < entries.size() * 2) cap <<= 1;
        keys.assign(cap, EMPTY);
        vals.resize(cap, 0);
        mask = cap - 1;
        count = 0;
        for (auto& [k, v] : entries) {
            uint32_t i = (uint32_t)(k * 0x9e3779b97f4a7c15ULL >> 32) & mask;
            while (keys[i] != EMPTY) i = (i + 1) & mask;
            keys[i] = k; vals[i] = v; ++count;
        }
    }

    // Returns UINT8_MAX if not found.
    uint8_t lookup(uint64_t k) const noexcept {
        if (keys.empty()) return UINT8_MAX;
        uint32_t i = (uint32_t)(k * 0x9e3779b97f4a7c15ULL >> 32) & mask;
        while (keys[i] != EMPTY) {
            if (keys[i] == k) return vals[i];
            i = (i + 1) & mask;
        }
        return UINT8_MAX;
    }

    bool empty() const { return count == 0; }
    size_t memory_bytes() const { return keys.size() * 9; }
};

// Blocked Bloom filter: one 64-byte cache-line block per bucket.
// k=4 independent bit probes per element. ~1% FP at 10 bits/element.
struct BlockedBloom {
    static constexpr int BITS_PER_BLOCK = 512; // 64 bytes
    static constexpr int K = 4;

    std::vector<uint64_t> data; // n_blocks * 8 uint64s
    uint32_t n_blocks = 0;

    void build(const FlatHitMap& map, float bits_per_elem = 10.0f) {
        if (map.empty()) return;
        n_blocks = std::max(1u, (uint32_t)(map.count * bits_per_elem / BITS_PER_BLOCK));
        data.assign((size_t)n_blocks * 8, 0);
        for (size_t i = 0; i < map.keys.size(); ++i)
            if (map.keys[i] != FlatHitMap::EMPTY) insert(map.keys[i]);
    }

    // Build directly from a sorted/unsorted hash array — no FlatHitMap needed.
    void build_from_hashes(const uint64_t* hashes, size_t n, float bits_per_elem = 10.0f) {
        if (n == 0) return;
        n_blocks = std::max(1u, (uint32_t)(n * bits_per_elem / BITS_PER_BLOCK));
        data.assign((size_t)n_blocks * 8, 0);
        for (size_t i = 0; i < n; ++i) insert(hashes[i]);
    }

    void insert(uint64_t h) noexcept {
        const uint32_t b = (uint32_t)(h % n_blocks) * 8;
        uint64_t mix = h;
        for (int i = 0; i < K; ++i) {
            mix = mix * 0x9e3779b97f4a7c15ULL ^ (mix >> 32);
            data[b + ((mix >> 9) & 7)] |= (uint64_t)1 << (mix & 63);
        }
    }

    bool might_contain(uint64_t h) const noexcept {
        if (data.empty()) return true;
        const uint32_t b = (uint32_t)(h % n_blocks) * 8;
        uint64_t mix = h;
        for (int i = 0; i < K; ++i) {
            mix = mix * 0x9e3779b97f4a7c15ULL ^ (mix >> 32);
            if (!((data[b + ((mix >> 9) & 7)] >> (mix & 63)) & 1)) return false;
        }
        return true;
    }
};

// ── Writer ────────────────────────────────────────────────────────────────────

class MarkerWriter {
public:
    // Register a marker family's global hash set (sorted, deduplicated).
    // bac markers: call with marker_id 0..n_bac-1.
    // arc markers: call with marker_id n_bac..n_bac+n_arc-1.
    void set_pool_entry(uint8_t marker_id, std::vector<uint64_t> sorted_hashes) {
        if (marker_id >= pool_.size()) pool_.resize(marker_id + 1);
        pool_[marker_id] = std::move(sorted_hashes);
    }

    // Add a per-lineage calibration entry.
    void add_lineage(uint64_t lineage_hash, uint8_t domain, uint16_t ref_count,
                     const std::vector<CalibSlot>& slots,
                     const std::vector<bool>& expected_markers) {
        LineageEntry e;
        e.hash       = lineage_hash;
        e.domain     = domain;
        e.ref_count  = ref_count;
        e.slots      = slots;
        e.expected   = expected_markers;
        lineages_.push_back(std::move(e));
    }

    // Write the complete .mrk file. n_bac/n_arc are the schema sizes (e.g., 120, 53).
    void finalize(const std::filesystem::path& path,
                  uint8_t n_bac, uint8_t n_arc, uint8_t k = AAMER_K,
                  uint16_t frac_scale = 1,
                  uint8_t alphabet = MRKR_ALPHABET_FULL_AA) const;

    // Set the per-genus redundancy calibration table (sorted by genus_hash).
    void set_redun_calib(std::vector<RedunCalibEntry> entries) {
        redun_calib_ = std::move(entries);
    }

    size_t n_lineages()    const { return lineages_.size(); }
    size_t n_pool_entries() const { return pool_.size(); }

private:
    struct LineageEntry {
        uint64_t             hash;
        uint8_t              domain;
        uint16_t             ref_count;
        std::vector<CalibSlot> slots;
        std::vector<bool>    expected;
    };
    std::vector<LineageEntry>            lineages_;
    std::vector<std::vector<uint64_t>>   pool_;        // indexed by marker_id
    std::vector<RedunCalibEntry>         redun_calib_; // sorted by genus_hash
};

// ── Reader ────────────────────────────────────────────────────────────────────

class MarkerReader {
public:
    void open(const std::filesystem::path& path) {
        mmap_.open(path);
        init(mmap_.data(), mmap_.size());
    }

    // Open from a raw buffer (e.g., for testing).
    void open(const uint8_t* data, uint64_t size) { init(data, size); }

    bool is_open()       const { return hdr_ != nullptr; }
    uint32_t n_lineages() const { return hdr_ ? hdr_->n_lineages : 0; }
    uint8_t  n_bac()      const { return hdr_ ? hdr_->n_bac_markers : 0; }
    uint8_t  n_arc()      const { return hdr_ ? hdr_->n_arc_markers : 0; }
    uint8_t  k()          const { return hdr_ ? hdr_->k : 0; }
    uint8_t  alphabet()   const { return hdr_ ? hdr_->alphabet : MRKR_ALPHABET_FULL_AA; }
    bool     is_murphy10()  const { return alphabet() == MRKR_ALPHABET_MURPHY10; }
    bool     is_dayhoff6()  const { return alphabet() == MRKR_ALPHABET_DAYHOFF6; }
    uint16_t frac_scale() const { return hdr_ ? (hdr_->frac_scale ? hdr_->frac_scale : 1) : 1; }
    uint64_t frac_max_hash() const {
        const uint32_t s = frac_scale();
        return (s <= 1) ? UINT64_MAX : UINT64_MAX / s;
    }

    // Find the calibration entry for a lineage (O(log n) binary search).
    // Returns {calib_header_ptr, calib_slots_ptr, expected_bitset_ptr} or {nullptr,...} if absent.
    struct CalibView {
        const CalibEntryHeader* header    = nullptr;
        const CalibSlot*        slots     = nullptr;
        const uint8_t*          expected  = nullptr;  // ceil(n_markers/8) bytes
        bool valid() const { return header != nullptr; }
        bool marker_expected(int idx) const {
            return expected && (expected[idx / 8] >> (idx % 8)) & 1;
        }
    };

    CalibView lookup_lineage(uint64_t lineage_hash) const {
        if (!hdr_) return {};
        const auto* lo = lookup_;
        const auto* hi = lo + hdr_->n_lineages;
        auto it = std::lower_bound(lo, hi, lineage_hash,
            [](const MarkerLookupEntry& e, uint64_t h) { return e.lineage_hash < h; });
        if (it == hi || it->lineage_hash != lineage_hash) return {};
        return calib_view_at(it->calib_data_off);
    }

    // Sorted hash array for one marker family (direct mmap view, no copy).
    std::span<const uint64_t> pool_hashes(uint8_t marker_id) const {
        if (!hdr_ || marker_id >= hdr_->n_bac_markers + hdr_->n_arc_markers) return {};
        const auto& pe = pool_idx_[marker_id];
        const auto* base = reinterpret_cast<const uint64_t*>(data_ + hdr_->pool_off + pe.hashes_off);
        return {base, pe.n_hashes};
    }

    // O(log n) membership test: is `hash` in marker_id's pool?
    bool marker_contains(uint8_t marker_id, uint64_t hash) const {
        auto sp = pool_hashes(marker_id);
        return std::binary_search(sp.begin(), sp.end(), hash);
    }

    // Build an in-RAM aamer→marker_idx map for the expected markers of one lineage.
    // marker_idx is local (0..n_markers-1 within the domain schema).
    // For bac domain: pool_id = marker_idx; for arc domain: pool_id = n_bac + marker_idx.
    // Returns empty map if calib is invalid.
    std::unordered_map<uint64_t, uint8_t>
    build_hit_map(const CalibView& calib) const {
        std::unordered_map<uint64_t, uint8_t> map;
        if (!calib.valid() || !hdr_) return map;
        const uint8_t n  = calib.header->n_markers;
        const uint8_t base_id = (calib.header->domain == MRKR_DOMAIN_ARC) ? hdr_->n_bac_markers : 0;
        map.reserve(static_cast<size_t>(n) * 600); // ~500 hashes/marker average
        for (uint8_t i = 0; i < n; ++i) {
            if (!calib.marker_expected(i)) continue;
            auto hashes = pool_hashes(base_id + i);
            for (uint64_t h : hashes)
                map.emplace(h, i);  // first writer wins on collision across markers
        }
        return map;
    }

    // Build a per-genus FlatHitMap from the expected markers of one lineage.
    // Caller should cache this and reuse across all genomes of the same genus.
    FlatHitMap build_flat_hit_map(const CalibView& calib) const {
        FlatHitMap fm;
        if (!calib.valid() || !hdr_) return fm;
        const uint8_t n  = calib.header->n_markers;
        const uint8_t base_id = (calib.header->domain == MRKR_DOMAIN_ARC) ? hdr_->n_bac_markers : 0;
        std::vector<std::pair<uint64_t,uint8_t>> entries;
        entries.reserve(static_cast<size_t>(n) * 600); // ~500 hashes/marker average
        for (uint8_t i = 0; i < n; ++i) {
            if (!calib.marker_expected(i)) continue;
            auto hashes = pool_hashes(base_id + i);
            for (uint64_t h : hashes) entries.emplace_back(h, i);
        }
        fm.build(entries);
        return fm;
    }

    // Build merged sorted pools (bac and arc separately) for fast single-pass scoring.
    // k-way merge of all per-marker pools — call once after open().
    // After this, merged_hashes_bac/arc() return contiguous sorted arrays where
    // merged_ids_bac/arc()[i] gives the local marker index for merged_hashes[i].
    // Cross-family filtering guarantees each hash appears in at most one marker pool.
    void build_merged_pool() {
        if (!hdr_) return;
        const uint8_t n_bac = hdr_->n_bac_markers;
        const uint8_t n_arc = hdr_->n_arc_markers;

        if (hdr_->merged_bac_off != 0) {
            // Fast path: pre-merged pool stored in file — just mmap pointers, no copy/merge.
            const auto load_pre = [&](uint32_t sec_off, uint8_t base, uint8_t count,
                                      std::vector<uint64_t>& out_h,
                                      std::vector<uint8_t>&  out_id) {
                size_t total = 0;
                for (uint8_t i = 0; i < count; ++i) total += pool_idx_[base + i].n_hashes;
                if (total == 0) return;
                const auto* h_ptr  = reinterpret_cast<const uint64_t*>(data_ + sec_off);
                const auto* id_ptr = reinterpret_cast<const uint8_t*>(h_ptr + total);
                out_h.assign(h_ptr,  h_ptr  + total);
                out_id.assign(id_ptr, id_ptr + total);
            };
            load_pre(hdr_->merged_bac_off, 0,     n_bac, merged_hashes_bac_, merged_ids_bac_);
            if (hdr_->merged_arc_off != 0)
                load_pre(hdr_->merged_arc_off, n_bac, n_arc, merged_hashes_arc_, merged_ids_arc_);
        } else {
            // Fallback: k-way merge for old pools without pre-merged section.
            build_merged(0,     n_bac, merged_hashes_bac_, merged_ids_bac_);
            build_merged(n_bac, n_arc, merged_hashes_arc_, merged_ids_arc_);
        }

        merged_bloom_bac_.build_from_hashes(merged_hashes_bac_.data(),
                                            merged_hashes_bac_.size());
        merged_bloom_arc_.build_from_hashes(merged_hashes_arc_.data(),
                                            merged_hashes_arc_.size());
    }

    bool has_merged_pool() const { return !merged_hashes_bac_.empty(); }
    uint32_t pool_n_hashes(uint8_t mi) const {
        return pool_idx_ ? pool_idx_[mi].n_hashes : 0;
    }


    bool has_redun_calib() const { return false; }

    // Redundancy calibration removed — always returns nullptr.
    const RedunCalibEntry* lookup_redun_calib(uint64_t) const { return nullptr; }

    // Compute z-score for observed redundancy given a calibration entry.
    // Returns NaN if calibration is absent or invalid.
    static float redun_zscore(float observed, const RedunCalibEntry* ce) {
        if (!ce) return NAN;
        const float mad_eff = std::max({ce->mad, 0.01f * ce->median, 1e-4f});
        return (observed - ce->median) / (1.4826f * mad_eff);
    }

    std::span<const uint64_t>  merged_hashes_bac() const { return merged_hashes_bac_; }
    std::span<const uint8_t>   merged_ids_bac()    const { return merged_ids_bac_; }
    std::span<const uint64_t>  merged_hashes_arc() const { return merged_hashes_arc_; }
    std::span<const uint8_t>   merged_ids_arc()    const { return merged_ids_arc_; }
    const BlockedBloom& merged_bloom_bac()          const { return merged_bloom_bac_; }
    const BlockedBloom& merged_bloom_arc()          const { return merged_bloom_arc_; }

private:
    MmapFileReader           mmap_;
    const uint8_t*           data_        = nullptr;
    const MarkerHeader*      hdr_         = nullptr;
    const MarkerLookupEntry* lookup_      = nullptr;
    const MarkerPoolEntry*   pool_idx_    = nullptr;
    const RedunCalibEntry*   redun_calib_ = nullptr;

    std::vector<uint64_t> merged_hashes_bac_;
    std::vector<uint8_t>  merged_ids_bac_;
    std::vector<uint64_t> merged_hashes_arc_;
    std::vector<uint8_t>  merged_ids_arc_;
    BlockedBloom          merged_bloom_bac_;
    BlockedBloom          merged_bloom_arc_;

    void build_merged(uint8_t base, uint8_t count,
                      std::vector<uint64_t>& out_hashes,
                      std::vector<uint8_t>&  out_ids) {
        // Preload all pools sequentially into RAM to avoid mmap page-fault storms
        // during the priority-queue merge (which has semi-random access patterns).
        std::vector<std::vector<uint64_t>> ram(count);
        size_t total = 0;
        for (uint8_t i = 0; i < count; ++i) {
            auto sp = pool_hashes(base + i);
            ram[i].assign(sp.begin(), sp.end());
            total += ram[i].size();
        }
        out_hashes.resize(total);
        out_ids.resize(total);

        // k-way merge: priority queue of (hash, local_marker_id, position_in_pool)
        using E = std::tuple<uint64_t, uint8_t, uint32_t>;
        std::priority_queue<E, std::vector<E>, std::greater<E>> pq;
        for (uint8_t i = 0; i < count; ++i)
            if (!ram[i].empty()) pq.push({ram[i][0], i, 0u});

        size_t out = 0;
        while (!pq.empty()) {
            auto [h, mi, pos] = pq.top(); pq.pop();
            out_hashes[out] = h;
            out_ids[out]    = mi;
            ++out;
            if (static_cast<size_t>(pos) + 1 < ram[mi].size())
                pq.push({ram[mi][pos + 1], mi, pos + 1});
        }
    }

    void init(const uint8_t* data, uint64_t size) {
        if (size < 32)  // v1 header minimum
            throw std::runtime_error("markers.mrk: file too small");
        data_ = data;
        hdr_  = reinterpret_cast<const MarkerHeader*>(data);
        if (hdr_->magic != MRKR_MAGIC)
            throw std::runtime_error("markers.mrk: bad magic");
        if (hdr_->version < MRKR_VERSION_MIN || hdr_->version > MRKR_VERSION)
            throw std::runtime_error("markers.mrk: unsupported version — rebuild required");
        const bool k_ok = (hdr_->k == AAMER_K) ||
                          (hdr_->k == AAMER_K_D6 &&
                           hdr_->alphabet == MRKR_ALPHABET_DAYHOFF6);
        if (!k_ok)
            throw std::runtime_error("markers.mrk: k mismatch — rebuild required");
        lookup_   = reinterpret_cast<const MarkerLookupEntry*>(data + hdr_->lookup_off);
        pool_idx_ = reinterpret_cast<const MarkerPoolEntry*>(data + hdr_->pool_idx_off);
        // redun_calib_ always null — field repurposed for merged_bac/arc_off in v2+
    }

    CalibView calib_view_at(uint32_t off) const {
        const uint8_t* p = data_ + hdr_->calib_off + off;
        CalibView v;
        v.header   = reinterpret_cast<const CalibEntryHeader*>(p);
        v.slots    = reinterpret_cast<const CalibSlot*>(p + sizeof(CalibEntryHeader));
        const uint8_t n = v.header->n_markers;
        v.expected = p + sizeof(CalibEntryHeader) + n * sizeof(CalibSlot);
        return v;
    }
};

// Combined per-genus marker acceleration index: flat hit map + bloom prefilter.
// Built once per genus and cached (see GenusIndexCache); reused across all genomes.
struct GenusMarkerIndex {
    FlatHitMap    map;
    BlockedBloom  bloom;

    void build(const MarkerReader& mrk_rd, const MarkerReader::CalibView& calib) {
        map = mrk_rd.build_flat_hit_map(calib);
        bloom.build(map);
    }
    size_t memory_bytes() const { return map.memory_bytes() + bloom.data.size() * 8; }
};

} // namespace genopack
