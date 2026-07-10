#pragma once
// ── SEC_PCORE — unified prevalence-annotated per-genus aamer reference ─────────
// One per-genus reference that supersedes CORE (prevalence core): it
// stores EVERY aamer seen in ≥1 member (the dense union) together with a per-aamer
// member prevalence. From this single section a consumer derives, at query time:
//   prevalence ≥ θ   → the conserved completeness core   (== old CORE; exact)
//   prevalence ≥ 1   → the DENSE foreign reference        (small-contig detection)
//   prevalence       → the per-aamer weight for a real log-likelihood ratio
//
// CODEC 0 (legacy): open-addressing hash key_hash → sorted u64 aamers + u8 *member
//   count* (capped 255 — BUG on genera with >255 members; read-only back-compat).
// CODEC 1 (v1): three contiguous PFOR-delta runs per genus —
//   [singleton: count==1, prev implicit] [multi: 1<count<C] [core: count≥C=⌈θ·n⌉]
//   with quantized (log-normalized) prevalence for multi+core; the core run is
//   decided from TRUE u32 counts at build, so the completeness slice is exact and
//   independent of prevalence quantization. Consumed by sequential scan only.
#include "format.hpp"
#include "mmap_file.hpp"
#include "pfor_codec.hpp"
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace genopack {

static constexpr uint8_t PCORE_CODEC_LEGACY = 0;   // dense u64 + u8 member count
static constexpr uint8_t PCORE_CODEC_V1     = 1;   // 3-run stratified PFOR + quant prev
static constexpr uint8_t PCORE_CODEC_V2     = 2;   // v1 + per-aamer IDF tier bytes
static constexpr uint8_t PCORE_PREVQ_NONE   = 0;   // raw count (legacy)
static constexpr uint8_t PCORE_PREVQ_LOGN   = 1;   // q = round(255·log(c)/log(n))

struct PcoreHeader {             // 64 bytes — legacy-compatible prefix
    uint32_t magic;              //  0  SEC_PCORE
    uint32_t n_entries;          //  4
    uint32_t n_buckets;          //  8  power of 2
    uint32_t k;                  // 12  aamer k
    uint32_t min_seg_aa;         // 16
    float    theta;              // 20  completeness-core prevalence threshold
    uint64_t frac_max_hash;      // 24  FMH threshold (pinned)
    uint64_t model_hash;         // 32
    uint64_t entries_offset;     // 40  -> PcoreEntry[]  (legacy) or PcoreEntryV1[]
    uint64_t buckets_offset;     // 48  -> uint32_t[n_buckets]
    uint8_t  codec;              // 56  PCORE_CODEC_*
    uint8_t  prev_quant;         // 57  PCORE_PREVQ_*
    uint8_t  frame_table_id;     // 58  translation-table id (pinned)
    uint8_t  _flags;             // 59
    uint32_t header_bytes;       // 60  64 (legacy) or 128 (v1, ext follows)
};
static_assert(sizeof(PcoreHeader) == 64);

// v1 extension, written at offset 64 (when codec==1). Pins the join invariants.
struct PcoreHeaderExtV1 {        // 64 bytes
    uint64_t fmh_seed;           //  0
    uint64_t taxonomy_hash;      //  8
    uint64_t k_frame_pin;        // 16  (k<<8)|frame_table_id — redundant guard
    uint64_t total_genera_G = 0; // 24  total genus count across the corpus (for IDF weighting)
    uint64_t tier_table_hash = 0;// 32  XXH3 of the global tier table (zero if not stamped)
    uint64_t _reserved[3];       // 40..64
};
static_assert(sizeof(PcoreHeaderExtV1) == 64);

struct PcoreEntry {              // 32 bytes — legacy
    uint64_t key_hash;
    uint64_t aamers_offset;
    uint64_t prev_offset;
    uint32_t n_aamers;
    uint32_t n_members;
};
static_assert(sizeof(PcoreEntry) == 32);

struct PcoreEntryV1 {            // 48 bytes — v1 per-genus run directory
    uint64_t key_hash;           //  0
    uint64_t run_offset;         //  8  section-relative start of this genus's runs (8-aligned)
    uint32_t n_singleton;        // 16
    uint32_t n_multi;            // 20
    uint32_t n_core;             // 24
    uint32_t n_members;          // 28
    uint32_t enc_singleton;      // 32  byte length of singleton PFOR stream
    uint32_t enc_multi;          // 36
    uint32_t enc_core;           // 40
    uint32_t _pad;               // 44
};
static_assert(sizeof(PcoreEntryV1) == 48);

static constexpr uint32_t PCORE_EMPTY = UINT32_MAX;

// ── prevalence quantization (codec v1) ────────────────────────────────────────
// q = round(255·log(count)/log(n)); reserve nothing — core membership is a run, not
// a quant value. dequant prevalence p = n^(q/255 − 1). Exact at count==n (q=255→p=1).
inline uint8_t pcore_quant_prev(uint32_t count, uint32_t n_members) noexcept {
    if (n_members <= 1 || count >= n_members) return 255;
    if (count <= 1) return 0;
    const double q = 255.0 * std::log(static_cast<double>(count)) / std::log(static_cast<double>(n_members));
    int qi = static_cast<int>(q + 0.5);
    if (qi < 1) qi = 1;
    if (qi > 255) qi = 255;
    return static_cast<uint8_t>(qi);
}
inline float pcore_dequant_prev(uint8_t q, uint32_t n_members) noexcept {
    if (n_members <= 1) return 1.0f;
    return static_cast<float>(std::pow(static_cast<double>(n_members), q / 255.0 - 1.0));
}

// IDF tier: tier4 = floor(log2(genus_count)), capped at 15.
// genus_count == 0 or 1 → tier 0 (ubiquitous within one genus → not informative inter-genus).
inline uint8_t pcore_tier_from_count(uint32_t genus_count) noexcept {
    if (genus_count <= 1) return 0;
    uint32_t t = 0, c = genus_count;
    while (c > 1 && t < 15) { c >>= 1; ++t; }
    return static_cast<uint8_t>(t);
}

// IDF weight from a 4-bit tier and corpus genus count G.
// w = 1 - tier4/log2(G), clamped to [0,1]. Returns 1 if G<=1 (single-genus corpus).
inline float pcore_idf_weight(uint8_t tier4, uint32_t G) noexcept {
    if (G <= 1) return 1.0f;
    const float denom = std::log2(static_cast<float>(G));
    const float w = 1.0f - static_cast<float>(tier4 & 0x0Fu) / denom;
    return w > 0.0f ? w : 0.0f;
}

struct PcoreView {
    uint64_t key_hash  = 0;
    uint32_t n_members = 0;
    float    theta     = 0.90f;
    bool     is_v1     = false;
    // legacy (codec 0): raw mmap pointers
    const uint64_t* aamers   = nullptr;   // sorted
    const uint8_t*  prev     = nullptr;   // member count (cap 255)
    uint32_t        n_aamers = 0;
    // v1 (codec 1): encoded runs (into mmap) + quant prev bytes
    const uint8_t* enc_s = nullptr; uint32_t n_singleton = 0; uint32_t enc_s_bytes = 0;
    const uint8_t* enc_m = nullptr; uint32_t n_multi     = 0; uint32_t enc_m_bytes = 0; const uint8_t* prevq_m = nullptr;
    const uint8_t* enc_c = nullptr; uint32_t n_core      = 0; uint32_t enc_c_bytes = 0; const uint8_t* prevq_c = nullptr;
    // v2 (codec 2): per-aamer IDF tier bytes (parallel to run order: singleton|multi|core)
    const uint8_t* tier_s = nullptr;   // tier byte per singleton aamer
    const uint8_t* tier_m = nullptr;   // tier byte per multi aamer
    const uint8_t* tier_c = nullptr;   // tier byte per core aamer
    bool     has_tier           = false;
    uint32_t total_genera_G     = 0;    // corpus-wide genus count (from ext header)

    uint32_t total() const noexcept { return is_v1 ? (n_singleton + n_multi + n_core) : n_aamers; }
    bool valid() const noexcept { return total() > 0; }

    // Decode the full union into aamers_out with dequantized prevalence prevf_out.
    void materialize(std::vector<uint64_t>& aamers_out, std::vector<float>& prevf_out) const {
        aamers_out.clear(); prevf_out.clear();
        if (!is_v1) {
            aamers_out.assign(aamers, aamers + n_aamers);
            prevf_out.resize(n_aamers);
            for (uint32_t i = 0; i < n_aamers; ++i)
                prevf_out[i] = (prev && n_members) ? static_cast<float>(prev[i]) / static_cast<float>(n_members) : 1.0f;
            return;
        }
        aamers_out.reserve(total()); prevf_out.reserve(total());
        std::vector<uint64_t> run;
        const float p_single = n_members ? 1.0f / static_cast<float>(n_members) : 1.0f;
        if (n_singleton) { pfor::decode_into(enc_s, enc_s_bytes, n_singleton, run);
            for (uint64_t h : run) { aamers_out.push_back(h); prevf_out.push_back(p_single); } }
        if (n_multi) { pfor::decode_into(enc_m, enc_m_bytes, n_multi, run);
            for (uint32_t i = 0; i < n_multi; ++i) { aamers_out.push_back(run[i]); prevf_out.push_back(pcore_dequant_prev(prevq_m[i], n_members)); } }
        if (n_core) { pfor::decode_into(enc_c, enc_c_bytes, n_core, run);
            for (uint32_t i = 0; i < n_core; ++i) { aamers_out.push_back(run[i]); prevf_out.push_back(pcore_dequant_prev(prevq_c[i], n_members)); } }
    }

    // Decode just the conserved core run (sorted) — the exact completeness slice.
    void core_into(std::vector<uint64_t>& out) const {
        if (is_v1) { pfor::decode_into(enc_c, enc_c_bytes, n_core, out); return; }
        // legacy: filter by count ≥ ⌈θ·n⌉ (capped-count BUG inherited from old archives)
        out.clear();
        const uint32_t thr = static_cast<uint32_t>(std::ceil(static_cast<double>(theta) * n_members));
        for (uint32_t i = 0; i < n_aamers; ++i) if (prev && prev[i] >= thr) out.push_back(aamers[i]);
    }
};

// ⌈theta · n_members⌉ — member-count threshold for the conserved completeness core.
inline uint32_t pcore_core_threshold(uint32_t n_members, float theta) {
    return static_cast<uint32_t>(std::ceil(static_cast<double>(theta) * n_members));
}

// Completeness coverage = |query ∩ core| / |core|. query_sorted must be sorted-unique.
inline float pcore_core_coverage(const PcoreView& v, const std::vector<uint64_t>& query_sorted, float theta) {
    (void)theta;                                  // core slice is intrinsic to v (v.theta for legacy)
    if (!v.valid() || v.n_members == 0) return std::numeric_limits<float>::quiet_NaN();
    std::vector<uint64_t> core;
    v.core_into(core);                            // exact for v1; legacy filters by count
    if (core.empty()) return std::numeric_limits<float>::quiet_NaN();
    uint32_t inter = 0; size_t qi = 0;
    for (uint64_t a : core) {                     // core is sorted ascending
        while (qi < query_sorted.size() && query_sorted[qi] < a) ++qi;
        if (qi < query_sorted.size() && query_sorted[qi] == a) ++inter;
    }
    return static_cast<float>(inter) / static_cast<float>(core.size());
}

// ── Writer (streaming) ────────────────────────────────────────────────────────
// Bounded memory: each genus's encoded runs are spilled to a temp file as added
// (only ~64 bytes/genus metadata stays in RAM). $GENOPACK_SPILL_DIR or system temp.
// Model hash folded incrementally and order-independently (XOR of per-entry hashes).
class PcoreWriter {
public:
    PcoreWriter(uint32_t k, uint32_t min_seg_aa, float theta = 0.90f,
                std::string spill_dir = {});
    ~PcoreWriter();
    PcoreWriter(const PcoreWriter&) = delete;
    PcoreWriter& operator=(const PcoreWriter&) = delete;

    // Pin the join invariants into the header (call before finalize; defaults if unset).
    void set_pins(uint64_t fmh_seed, uint64_t taxonomy_hash, uint8_t frame_table_id) {
        fmh_seed_ = fmh_seed; taxonomy_hash_ = taxonomy_hash; frame_table_id_ = frame_table_id;
    }

    // Open a side-channel .ptier file that records (aamer_hash, genus_hash) pairs
    // for every aamer in every genus. Used post-build by `genopack tier merge`.
    // Call before add_sorted/add_from_members/add_from_counts; no-op if path is empty.
    void set_tier_sidechannel(const std::string& path);

    // member_qmers = per-member sorted-unique aamer sets. Builds the dense union with
    // true u32 member counts, stratifies into singleton/multi/core runs, PFOR-encodes.
    void add_from_members(uint64_t key_hash, const std::vector<std::vector<uint64_t>>& member_qmers);

    // Streaming hash-aggregation entry point: caller maintains one per-genus
    // aamer→member-count map (one member's aamers resident at a time), so peak
    // memory is a single genus's union — bounded even for 600k-member genera.
    void add_from_counts(uint64_t key_hash, const std::unordered_map<uint64_t, uint32_t>& cnt,
                         uint32_t n_members);
    // Fast path: caller supplies the genus union as globally-sorted (aamer, count)
    // arrays (from GenusAamerCounter::finalize_sorted) — no internal re-sort.
    void add_sorted(uint64_t key_hash, const std::vector<uint64_t>& aamers,
                    const std::vector<uint32_t>& counts, uint32_t n_members);

    uint64_t model_hash() const;
    SectionDesc finalize(AppendWriter& w, uint64_t section_id, uint64_t frac_max_hash);
    size_t n_entries() const { return meta_.size(); }

private:
    struct EntryMeta {
        uint64_t key_hash;
        uint64_t run_off;            // offset within the spilled pool
        uint32_t n_singleton, n_multi, n_core, n_members;
        uint32_t enc_singleton, enc_multi, enc_core;
    };
    static uint32_t next_pow2(uint32_t v) noexcept {
        if (v == 0) return 1;
        --v; v |= v>>1; v |= v>>2; v |= v>>4; v |= v>>8; v |= v>>16; return v + 1;
    }
    void open_spill_();

    uint32_t k_, min_seg_aa_;
    float    theta_;
    uint64_t fmh_seed_ = 0, taxonomy_hash_ = 0;
    uint8_t  frame_table_id_ = 0;
    std::vector<EntryMeta> meta_;
    std::FILE*  spill_       = nullptr;
    std::string spill_dir_;   // explicit override; empty = use env/auto-detect
    std::string spill_path_;
    uint64_t    pool_cursor_ = 0;
    uint64_t    hash_fold_   = 0;
    // Side-channel .ptier file state
    std::FILE*  tier_sc_      = nullptr;
    std::string tier_sc_path_;
    uint64_t    tier_sc_n_pairs_ = 0;  // running count of emitted (aamer, genus) pairs
};

// ── Reader ────────────────────────────────────────────────────────────────────
class PcoreReader {
public:
    void open(const uint8_t* data, uint64_t offset, uint64_t size) {
        if (size < sizeof(PcoreHeader)) throw std::runtime_error("PCORE section too small");
        data_   = data + offset;
        header_ = reinterpret_cast<const PcoreHeader*>(data_);
        if (header_->magic != SEC_PCORE) throw std::runtime_error("PCORE: bad magic");
        codec_ = header_->codec;
        buckets_ = reinterpret_cast<const uint32_t*>(data_ + header_->buckets_offset);
        if (codec_ == PCORE_CODEC_V1) {
            if (size < 128) throw std::runtime_error("PCORE v1 header truncated");
            ext_      = reinterpret_cast<const PcoreHeaderExtV1*>(data_ + sizeof(PcoreHeader));
            entries1_ = reinterpret_cast<const PcoreEntryV1*>(data_ + header_->entries_offset);
        } else {
            entries0_ = reinterpret_cast<const PcoreEntry*>(data_ + header_->entries_offset);
        }
    }

    PcoreView lookup(uint64_t key_hash) const noexcept {
        if (!header_ || header_->n_buckets == 0) return {};
        const uint32_t mask = header_->n_buckets - 1;
        uint32_t slot = static_cast<uint32_t>(key_hash) & mask;
        for (uint32_t probe = 0; probe < header_->n_buckets; ++probe) {
            const uint32_t idx = buckets_[slot];
            if (idx == PCORE_EMPTY) return {};
            const uint64_t kh = (codec_ == PCORE_CODEC_V1) ? entries1_[idx].key_hash : entries0_[idx].key_hash;
            if (kh == key_hash) return view_at(idx);
            slot = (slot + 1) & mask;
        }
        return {};
    }

    bool     is_open() const noexcept { return data_ != nullptr; }
    uint32_t n_entries() const noexcept { return header_ ? header_->n_entries : 0; }
    uint64_t key_hash_at(uint32_t i) const noexcept {
        if (!header_ || i >= header_->n_entries) return 0;
        return (codec_ == PCORE_CODEC_V1) ? entries1_[i].key_hash : entries0_[i].key_hash;
    }
    uint64_t total_union_aamers() const noexcept {
        if (!header_) return 0;
        uint64_t total = 0;
        if (codec_ == PCORE_CODEC_V1 || codec_ == PCORE_CODEC_V2) {
            for (uint32_t i = 0; i < header_->n_entries; ++i)
                total += uint64_t(entries1_[i].n_singleton) + entries1_[i].n_multi + entries1_[i].n_core;
        } else {
            for (uint32_t i = 0; i < header_->n_entries; ++i)
                total += entries0_[i].n_aamers;
        }
        return total;
    }
    int      k() const noexcept { return header_ ? static_cast<int>(header_->k) : 0; }
    int      min_seg_aa() const noexcept { return header_ ? static_cast<int>(header_->min_seg_aa) : 0; }
    float    theta() const noexcept { return header_ ? header_->theta : 0.90f; }
    uint64_t frac_max_hash() const noexcept { return header_ ? header_->frac_max_hash : 0; }
    uint64_t fmh_seed() const noexcept { return ext_ ? ext_->fmh_seed : 0; }
    uint64_t taxonomy_hash() const noexcept { return ext_ ? ext_->taxonomy_hash : 0; }
    uint8_t  codec() const noexcept { return codec_; }
    uint64_t model_hash() const noexcept { return header_ ? header_->model_hash : 0; }

private:
    PcoreView view_at(uint32_t idx) const noexcept {
        PcoreView v;
        v.theta = header_->theta;
        if (codec_ == PCORE_CODEC_V1 || codec_ == PCORE_CODEC_V2) {
            const auto& e = entries1_[idx];
            v.is_v1 = true; v.key_hash = e.key_hash; v.n_members = e.n_members;
            v.n_singleton = e.n_singleton; v.n_multi = e.n_multi; v.n_core = e.n_core;
            v.enc_s_bytes = e.enc_singleton; v.enc_m_bytes = e.enc_multi; v.enc_c_bytes = e.enc_core;
            const uint8_t* p = data_ + e.run_offset;
            v.enc_s = p;                       p += e.enc_singleton;
            v.enc_m = p;                       p += e.enc_multi;
            v.enc_c = p;                       p += e.enc_core;
            v.prevq_m = p;                     p += e.n_multi;
            v.prevq_c = p;                     p += e.n_core;
            if (codec_ == PCORE_CODEC_V2) {
                // tier bytes appended after prevq_c in run order: singleton|multi|core
                v.tier_s = p;                  p += e.n_singleton;
                v.tier_m = p;                  p += e.n_multi;
                v.tier_c = p;
                v.has_tier = true;
                v.total_genera_G = ext_ ? static_cast<uint32_t>(ext_->total_genera_G) : 0;
            }
        } else {
            const auto& e = entries0_[idx];
            v.is_v1 = false; v.key_hash = e.key_hash; v.n_members = e.n_members; v.n_aamers = e.n_aamers;
            v.aamers = reinterpret_cast<const uint64_t*>(data_ + e.aamers_offset);
            v.prev   = reinterpret_cast<const uint8_t*>(data_ + e.prev_offset);
        }
        return v;
    }

    const uint8_t*          data_     = nullptr;
    const PcoreHeader*      header_   = nullptr;
    const PcoreHeaderExtV1* ext_      = nullptr;
    const PcoreEntry*       entries0_ = nullptr;
    const PcoreEntryV1*     entries1_ = nullptr;
    const uint32_t*         buckets_  = nullptr;
    uint8_t                 codec_    = PCORE_CODEC_LEGACY;
};

} // namespace genopack
