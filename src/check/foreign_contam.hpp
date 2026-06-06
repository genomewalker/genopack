#pragma once
// ── Per-contig foreign-aamer containment (small-contig contamination) ─────────
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

} // namespace genopack::check
