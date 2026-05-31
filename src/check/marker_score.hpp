#pragma once
#include <genopack/markers.hpp>
#include <genopack/metamer.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

namespace genopack::check {

struct MarkerScoreResult {
    float completeness = NAN;
    float redundancy   = NAN;
    int   n_present    = -1;
    int   n_expected   = -1;
};

// Per-contig SCG marker scoring.
// Each contig votes once per marker (argmax of best ORF score in that contig).
// contig_votes[mi] = number of contigs where marker mi was the best argmax match.
// redundancy = markers with >=2 contig votes (two separate genomic copies).
// On a clean single-contig genome: max 1 vote per marker -> 0% redundancy.
//
// Lookup path selection:
//   - genus_idx != nullptr: bloom prefilter + per-genus FlatHitMap (Phase 1/2, fast).
//   - genus_idx == nullptr: binary search over the merged sorted pool (legacy fallback).
//
// Templated on the contig container so it works with the caller's local contig type
// (must expose a `.seq` member convertible to std::string_view).
template <typename Contigs>
MarkerScoreResult score_markers_for_genome(
    const Contigs&                        contigs,
    const MarkerReader&                   mrk_rd,
    const MarkerReader::CalibView&        calib,
    uint32_t                              min_hits,
    uint64_t                              mrk_frac_max,
    const GenusMarkerIndex*               genus_idx)
{
    MarkerScoreResult res;
    if (!calib.valid()) return res;

    const bool is_arc = (calib.header->domain == MRKR_DOMAIN_ARC);
    // Only fetch merged pool spans if falling back to binary search (genus_idx == nullptr).
    auto mh  = genus_idx ? std::span<const uint64_t>{} : (is_arc ? mrk_rd.merged_hashes_arc() : mrk_rd.merged_hashes_bac());
    auto mid = genus_idx ? std::span<const uint8_t>{}  : (is_arc ? mrk_rd.merged_ids_arc()    : mrk_rd.merged_ids_bac());
    const uint8_t n_markers = calib.header->n_markers;

    constexpr int MRK_MIN_ORF_AA  = 50;
    constexpr int MRK_MIN_FRAC_PC = 3;

    thread_local std::vector<uint64_t> orf_mers;
    uint32_t contig_votes[173] = {};

    for (const auto& contig : contigs) {
        // Per-contig best score per marker (across all ORFs in this contig).
        // Each marker votes independently: one vote per contig per marker that
        // clears threshold. Avoids argmax suppression where a contig containing
        // RpoB + RecA would vote for only one of them.
        uint32_t contig_best[173] = {};

        translate_6frame(contig.seq, MRK_MIN_ORF_AA,
            [&](int /*frame*/, const uint8_t* seg, int len, int, int) {
                orf_mers.clear();
                for (int i = 0; i + METAMER_K <= len; ++i) {
                    if (metamer_is_low_complexity(seg + i, METAMER_K)) continue;
                    const uint64_t h = mrk_rd.is_murphy10()
                        ? metamer_hash_murphy10(seg + i)
                        : metamer_hash(seg + i);
                    if (h <= mrk_frac_max) orf_mers.push_back(h);
                }
                if (orf_mers.empty()) return;
                std::sort(orf_mers.begin(), orf_mers.end());
                orf_mers.erase(std::unique(orf_mers.begin(), orf_mers.end()),
                               orf_mers.end());

                uint32_t local[173] = {};
                const uint64_t* mhp = mh.data();
                const uint8_t*  mip = mid.data();
                const uint64_t* mhe = mhp + mh.size();
                for (uint64_t h : orf_mers) {
                    if (genus_idx) {
                        if (!genus_idx->bloom.might_contain(h)) continue;
                        const uint8_t mi = genus_idx->map.lookup(h);
                        if (mi != UINT8_MAX) local[mi]++;
                    } else {
                        auto it = std::lower_bound(mhp, mhe, h);
                        if (it != mhe && *it == h) local[mip[it - mhp]]++;
                    }
                }
                // Update per-marker best across ORFs (no argmax winner-take-all).
                const uint32_t thr = std::max(min_hits,
                    (uint32_t)(orf_mers.size() * MRK_MIN_FRAC_PC / 100));
                for (uint8_t mi = 0; mi < n_markers; ++mi)
                    if (local[mi] >= thr && local[mi] > contig_best[mi])
                        contig_best[mi] = local[mi];
            });

        // Vote for every marker that this contig clears threshold for.
        for (uint8_t mi = 0; mi < n_markers; ++mi)
            if (contig_best[mi] > 0) contig_votes[mi]++;
    }

    res.n_present  = 0;
    res.n_expected = 0;
    int redundant_n = 0;
    for (uint8_t mi = 0; mi < n_markers; ++mi) {
        if (!calib.marker_expected(mi)) continue;
        ++res.n_expected;
        if (contig_votes[mi] >= 1) ++res.n_present;
        if (contig_votes[mi] >= 2) ++redundant_n;
    }
    if (res.n_expected > 0) {
        res.completeness = static_cast<float>(res.n_present)
                         / static_cast<float>(res.n_expected);
        res.redundancy   = static_cast<float>(redundant_n)
                         / static_cast<float>(res.n_expected);
    }
    return res;
}

} // namespace genopack::check
