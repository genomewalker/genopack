#pragma once
// ── Per-contig foreign-aamer containment (small-contig contamination) ─────────
#include <genopack/gami.hpp>   // GlobalMultiplicityIndex, GamiView
// The TNF/GCOV per-contig outlier signal needs ≥5–20kb contigs; below that
// tetranucleotide composition is too noisy. But most real contamination lives in
// SMALL contigs. This module adds a composition-INDEPENDENT per-contig signal that
// works down to ~1–2kb: protein-aamer (k=8, 6-frame) CONTAINMENT against the dense
// PCORE reference (every per-genus aamer + member count).
//
// For each contig, unique 6-frame aamers are classified against a precomputed
// reference set { host genus ∪ host family ∪ same-family siblings ∪ background }:
//   host-specific    : in the host genus reference
//   foreign-specific : in some foreign reference, NOT host, NOT family-conserved,
//                      present in ≤3 foreign references (taxon-specific, not universal)
//   family/conserved : in host family or >3 foreign references → excluded
//   foreign_fraction = foreign_specific / (host_specific + foreign_specific)   (DETECTION)
//   prot_loglr       = Σ_specific log((p_foreign+ε)/(p_host+ε))                (prevalence-weighted)
// where p_* are the per-aamer member-prevalences from PCORE (1.0 for the legacy
// presence-only CORE fallback). Detection stays presence-based; the log-LR is the
// richer statistic the design room wanted for the intra-family case.
#include <genopack/core_section.hpp>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace genopack::check {

// A reference aamer span. `prev`/`n_members` carry PCORE member counts (per-aamer
// prevalence); for the sparse CORE fallback they are null → prevalence defaults 1.0.
struct AamerSpan {
    const uint64_t* aamers = nullptr;
    uint32_t        n      = 0;
    const float*    prevf  = nullptr;   // dequantized prevalence in (0,1], parallel; null => 1.0
    bool  valid() const noexcept { return aamers != nullptr && n > 0; }
    float prevalence(uint32_t i) const noexcept { return prevf ? prevf[i] : 1.0f; }
};

// Owns a decoded per-genus union (PCORE v1 runs decode here; legacy/CORE copy here)
// so an AamerSpan view stays valid for the lifetime of the refset build.
struct OwnedSpan {
    std::vector<uint64_t> aamers;
    std::vector<float>    prevf;
    bool valid() const noexcept { return !aamers.empty(); }
    AamerSpan view() const noexcept {
        return { aamers.data(), static_cast<uint32_t>(aamers.size()),
                 prevf.empty() ? nullptr : prevf.data() };
    }
};

// 16-bit prefix index over a sorted OwnedSpan.
// offsets[b] = first entry index where h>>48 == b; offsets[65536] = n (sentinel).
// Avg bucket width ≈ n/65536. For a 10M-entry host union: ~153 entries/bucket.
// contains() costs ~7 comparisons in a contiguous 1.2 KB span vs a full 80 MB scan.
struct IndexedSpan {
    const uint64_t* data  = nullptr;
    const float*    prevf = nullptr;  // parallel prevalence; null → 1.0
    uint32_t        n     = 0;
    std::vector<uint32_t> offsets;   // 65537 entries, 256 KB

    bool valid() const noexcept { return data != nullptr && n > 0; }
    float prevalence(uint32_t i) const noexcept { return prevf ? prevf[i] : 1.0f; }

    // Returns index into data if found, n if absent.
    uint32_t find(uint64_t h) const noexcept {
        if (offsets.empty()) return n;
        const uint32_t b     = static_cast<uint32_t>(h >> 48);
        const uint64_t* lo   = data + offsets[b];
        const uint64_t* hi   = data + offsets[b + 1];
        const uint64_t* it   = std::lower_bound(lo, hi, h);
        return (it != hi && *it == h) ? static_cast<uint32_t>(it - data) : n;
    }
    bool contains(uint64_t h) const noexcept { return find(h) < n; }
};

inline IndexedSpan build_index(const OwnedSpan& os) {
    IndexedSpan ix;
    ix.data  = os.aamers.data();
    ix.prevf = os.prevf.empty() ? nullptr : os.prevf.data();
    ix.n     = static_cast<uint32_t>(os.aamers.size());
    ix.offsets.assign(65537, 0);
    for (uint64_t h : os.aamers) ++ix.offsets[static_cast<uint32_t>(h >> 48)];
    uint32_t acc = 0;
    for (int b = 0; b <= 65535; ++b) { uint32_t c = ix.offsets[b]; ix.offsets[b] = acc; acc += c; }
    ix.offsets[65536] = ix.n;
    return ix;
}

// Packed per-aamer reference info.
struct ForeignAamerInfo {
    float    host_prev    = 0.0f;   // prevalence in host genus reference (0 = absent)
    float    foreign_prev = 0.0f;   // best foreign reference prevalence (0 = absent)
    uint8_t  fcount       = 0;      // # foreign references containing it (cap 255)
    uint8_t  in_family    = 0;      // present in host family reference
    uint32_t genus_tag    = 0;      // dominant foreign genus tag (22-bit)
};

struct ForeignRefSet {
    std::unordered_map<uint64_t, ForeignAamerInfo> m;
    bool     valid          = false;
    uint32_t host_core_size = 0;
};

// Build the reference map once per host genus (single-threaded; read-only in the
// parallel contig loop). foreign = (genus_hash, span) for siblings + background.
inline ForeignRefSet build_foreign_refset(const AamerSpan& host, const AamerSpan& family,
                                          const std::vector<std::pair<uint64_t, AamerSpan>>& foreign) {
    ForeignRefSet rs;
    rs.valid = host.valid();
    rs.host_core_size = host.n;
    if (!rs.valid) return rs;

    rs.m.reserve(static_cast<size_t>(host.n) * 2 + 64);
    for (uint32_t i = 0; i < host.n; ++i) rs.m[host.aamers[i]].host_prev = host.prevalence(i);
    if (family.valid())
        for (uint32_t i = 0; i < family.n; ++i) rs.m[family.aamers[i]].in_family = 1;

    for (const auto& [ghash, fc] : foreign) {
        if (!fc.valid()) continue;
        const uint32_t tag = static_cast<uint32_t>(ghash & 0x3FFFFFu);
        for (uint32_t i = 0; i < fc.n; ++i) {
            ForeignAamerInfo& v = rs.m[fc.aamers[i]];
            if (v.fcount < 0xFFu) ++v.fcount;
            const float p = fc.prevalence(i);
            if (p > v.foreign_prev) { v.foreign_prev = p; v.genus_tag = tag; }
        }
    }
    return rs;
}

struct ForeignScore {
    uint32_t classifiable      = 0;   // host_specific + foreign_specific
    uint32_t host_specific     = 0;
    uint32_t foreign_specific  = 0;
    uint32_t family_hits       = 0;
    uint32_t best_genus_tag    = 0;
    float    host_prev_sum     = 0.0f; // sum of host aamer prevalences (weighted denominator)
    float    foreign_fraction  = NAN;
    float    score             = NAN; // prevalence-weighted log-LR
};

// Score one contig's sorted-unique aamer set against a reference map.
inline ForeignScore score_contig_foreign(const ForeignRefSet& rs,
                                         const std::vector<uint64_t>& contig_aamers) {
    ForeignScore s;
    if (!rs.valid) return s;
    constexpr float eps = 0.01f;
    double loglr = 0.0;
    std::unordered_map<uint32_t, uint32_t> genus_tally;
    for (uint64_t h : contig_aamers) {
        auto it = rs.m.find(h);
        if (it == rs.m.end()) continue;             // not in any reference → not classifiable
        const ForeignAamerInfo& v = it->second;
        const bool in_host = v.host_prev > 0.0f;
        if (v.in_family) ++s.family_hits;
        if (in_host) {
            ++s.host_specific;
            loglr += std::log((v.foreign_prev + eps) / (v.host_prev + eps));
        } else if (v.fcount >= 1 && !v.in_family && v.fcount <= 3) {   // taxon-specific foreign
            ++s.foreign_specific;
            ++genus_tally[v.genus_tag];
            loglr += std::log((v.foreign_prev + eps) / (0.0f + eps));
        }
    }
    s.classifiable = s.host_specific + s.foreign_specific;
    if (s.classifiable > 0) {
        s.foreign_fraction = static_cast<float>(s.foreign_specific) / static_cast<float>(s.classifiable);
        s.score = static_cast<float>(loglr / s.classifiable);   // mean per-aamer log-LR (length-robust)
    }
    uint32_t best = 0, best_n = 0;
    for (const auto& [g, n] : genus_tally) if (n > best_n) { best_n = n; best = g; }
    s.best_genus_tag = best;
    return s;
}

// ── TierView — IDF tier lookup from a sorted (aamer_hash, tier) array ─────────
// Loaded from the per-pack V2 PCORE tier_s/m/c bytes or from a global .gtier
// side-table. When valid(), replaces the GAMI/GMI rare-aamer classification with
// a soft IDF weight: w = pcore_idf_weight(tier4, G) ∈ (0,1].
struct TierView {
    // Sorted (aamer_hash, tier) pairs (8+1 bytes each, packed).
    struct Entry { uint64_t h; uint8_t t; } __attribute__((packed));
    const Entry* data = nullptr;
    uint64_t     n    = 0;
    uint32_t     G    = 0;   // total corpus genera (for IDF weight denominator)

    bool valid() const noexcept { return data != nullptr && n > 0 && G > 1; }

    // Returns IDF weight in (0,1], or 1.0 if aamer not found (treat as rare).
    float weight(uint64_t aamer_hash) const noexcept {
        if (!valid()) return 1.0f;
        uint64_t lo = 0, hi = n;
        while (lo < hi) {
            const uint64_t mid = lo + (hi - lo) / 2;
            const uint64_t mh  = data[mid].h;
            if (mh == aamer_hash) return genopack::pcore_idf_weight(data[mid].t, G);
            if (mh  < aamer_hash) lo = mid + 1;
            else                  hi = mid;
        }
        return 1.0f;  // absent → treat as fully specific
    }
};

// GlobalMultiplicityIndex is defined in <genopack/gami.hpp> (namespace genopack).
using genopack::GlobalMultiplicityIndex;
using genopack::GamiView;

// Score one contig's sorted-unique aamer set using the GMI.
// host/family are sorted AamerSpans from the host's PCORE union; merge-scan
// avoids per-aamer binary search. K = max genus-count for an aamer to be
// considered taxon-specific foreign (equivalent to the old fcount ≤ 3).
inline ForeignScore score_contig_foreign_v2(const GlobalMultiplicityIndex& gmi,
                                            const AamerSpan& host,
                                            const AamerSpan& family,
                                            const std::vector<uint64_t>& contig_aamers,
                                            uint8_t K = 3) {
    ForeignScore s;
    if (!host.valid() || gmi.empty()) return s;
    constexpr float eps = 0.01f;
    double loglr = 0.0;
    uint32_t hi = 0, fi = 0;
    for (uint64_t q : contig_aamers) {
        while (hi < host.n && host.aamers[hi] < q) ++hi;
        const bool in_host = (hi < host.n && host.aamers[hi] == q);
        bool in_fam = false;
        if (family.valid()) {
            while (fi < family.n && family.aamers[fi] < q) ++fi;
            in_fam = (fi < family.n && family.aamers[fi] == q);
        }
        if (in_fam) ++s.family_hits;
        if (in_host) {
            ++s.host_specific;
            const float hp = host.prevalence(hi);
            loglr += std::log(eps / (hp + eps));
        } else if (!in_fam) {
            const uint8_t gc = gmi.lookup(q);
            if (gc >= 1 && gc <= K) {
                ++s.foreign_specific;
                loglr += std::log(1.0f / eps);  // simplified: no per-aamer foreign_prev in GMI
            }
        }
    }
    s.classifiable = s.host_specific + s.foreign_specific;
    if (s.classifiable > 0) {
        s.foreign_fraction = static_cast<float>(s.foreign_specific) /
                             static_cast<float>(s.classifiable);
        s.score = static_cast<float>(loglr / static_cast<double>(s.classifiable));
    }
    return s;
}

// Unified scorer using IndexedSpan for O(log 153) host/family lookup instead of
// O(10M) merge-scan. Handles TierView (IDF weights), GAMI (Bloom), and GMI
// (hash table) confirmation — in preference order.
inline ForeignScore score_contig_foreign_indexed(
        const IndexedSpan& host,
        const IndexedSpan& family,
        const GlobalMultiplicityIndex& gmi,
        const GamiView& gv,
        uint8_t K,
        const std::vector<uint64_t>& contig_aamers,
        const TierView& tier_view = TierView{}) {
    ForeignScore s;
    if (!host.valid()) return s;
    constexpr float eps = 0.01f;
    double loglr = 0.0;
    for (uint64_t q : contig_aamers) {
        const uint32_t hpos    = host.find(q);
        const bool     in_host = hpos < host.n;
        const bool     in_fam  = family.valid() && family.contains(q);
        if (in_fam) ++s.family_hits;
        if (in_host) {
            ++s.host_specific;
            const float hp = host.prevalence(hpos);
            s.host_prev_sum += hp;
            loglr += std::log(eps / (hp + eps));
        } else if (!in_fam) {
            if (tier_view.valid()) {
                // IDF path: soft weight replaces binary is_rare.
                const float w = tier_view.weight(q);
                if (w > 0.0f) {
                    ++s.foreign_specific;
                    loglr += w * std::log(1.0f / eps);
                }
            } else if (gv.valid()) {
                // Legacy GAMI path (prebuilt Bloom).
                if (gv.is_rare(q)) { ++s.foreign_specific; loglr += std::log(1.0f / eps); }
            } else if (!gmi.empty()) {
                // Legacy runtime-GMI path.
                const uint8_t gc = gmi.lookup(q);
                if (gc >= 1 && gc <= K) { ++s.foreign_specific; loglr += std::log(1.0f / eps); }
            } else {
                // No tier / GAMI / GMI: unweighted baseline — every foreign aamer counts.
                ++s.foreign_specific;
                loglr += std::log(1.0f / eps);
            }
        }
    }
    s.classifiable = s.host_specific + s.foreign_specific;
    if (s.classifiable >= 16) {
        // Prevalence-weighted fraction: host side weighted by within-genus prevalence,
        // foreign side uniform (unknown foreign prevalence). Amplifies signal for small
        // or noisy genera where sampled host aamers have low prevalence.
        const float denom = s.host_prev_sum + static_cast<float>(s.foreign_specific);
        s.foreign_fraction = denom > 0.0f
            ? static_cast<float>(s.foreign_specific) / denom
            : static_cast<float>(s.foreign_specific) / static_cast<float>(s.classifiable);
        s.score = static_cast<float>(loglr / static_cast<double>(s.classifiable));
    } else if (family.valid() && !contig_aamers.empty()) {
        // Sub-floor (classifiable < 16): too few host/foreign aamers for the standard scorer.
        // Use family-hit density as contamination proxy: a family member should hit roughly
        // expected = |contig| * |family| / hash_space aamers; far fewer hits → foreign.
        const float expected = static_cast<float>(contig_aamers.size()) *
                               static_cast<float>(family.n) / 2e6f;
        if (expected >= 1.0f) {
            s.foreign_fraction = 1.0f - std::min(1.0f,
                static_cast<float>(s.family_hits) / expected);
        }
    }
    return s;
}

} // namespace genopack::check
