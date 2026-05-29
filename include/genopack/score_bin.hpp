#pragma once
#include "fmhr.hpp"
#include "gcov.hpp"
#include "gstx.hpp"
#include "skani.hpp"
#include <cmath>
#include <cstdint>
#include <string_view>
#include <vector>

namespace genopack {

struct BinScoreConfig {
    float    outlier_percentile = 0.99f;
    uint32_t min_long_bp        = 20000;
};

struct BinScores {
    enum class Status : uint8_t { OK, NoGcovEntry, TooFewContigs };

    Status   status            = Status::NoGcovEntry;
    float    cco_fraction      = NAN;  // (T²|SPE outlier bp) / scored_bp
    float    sibling_fraction  = NAN;  // genus-outlier ∩ family-inlier bp / scored_bp
    float    rho_fraction      = NAN;  // ρ* Mahalanobis outlier bp / scored_bp
    uint32_t scored_bp         = 0;
};

// Per-contig containment split: for each contig build an OPH sketch and find its
// best-matching genus by comparing against GSTX consensus sketches.
// Returns minority_fraction = bp(best_match != primary_genus) / scored_bp.
// NaN if candidate_genera has < 2 entries or no contig meets min_contig_bp.
// Thread-safe: no shared mutable state.
struct ContainmentSplitResult {
    float minority_fraction     = NAN;  // bp(best-match genus != primary) / scored_bp
    float self_outlier_fraction = NAN;  // bp of self-containment outliers (<mean-2σ) / scored_bp
    float fiedler_oph           = NAN;  // sketch Fiedler: 1−λ₂/2 of OPH-Jaccard normalized Laplacian
    float fiedler_tnf_bimod     = NAN;  // TNF-kernel: normalized max-gap bimodality of Fiedler v₂
    float fiedler_tnf_gap       = NAN;  // TNF-kernel: λ₃−λ₂ eigengap of normalized Laplacian
};

ContainmentSplitResult score_bin_containment(
    const std::vector<std::string_view>& seqs,
    uint64_t primary_genus_hash,
    const std::vector<const GstxEntry*>& candidate_genera,
    const GstxReader& gstx_rd,
    uint32_t min_contig_bp = 10000);

// ── FracMinHash containment sweep ─────────────────────────────────────────────

// Merged FMH reference sketch for one genus (sorted union of all member hashes).
struct FmhGenusRef {
    uint64_t              genus_hash;
    std::vector<uint64_t> hashes;  // sorted, deduplicated
};

// Build a merged FMH genus reference from raw genome sequences (not FASTA — plain seq).
// Call once per (genus, k, c) combination before scoring bins.
FmhGenusRef build_fmh_genus_ref(uint64_t genus_hash,
                                  const std::vector<std::string_view>& genome_seqs,
                                  int k, int c);

// For each contig >= min_contig_bp, build an FMH sketch and compute containment vs
// each genus ref. Returns fraction of bin bp where best-match genus != primary_genus_hash.
// Returns NaN if no contig is long enough or refs has < 2 entries.
float score_bin_fmh_containment(
    const std::vector<std::string_view>& bin_seqs,
    uint64_t primary_genus_hash,
    const std::vector<FmhGenusRef>& refs,
    int k, int c,
    uint32_t min_contig_bp = 10000);

// ── Multi-(k,c) sweep — 2 sequence scans instead of 6 ────────────────────────
// Fixed sweep: k∈{21,31} × c∈{30,125,500}, index = ki*3+ci.

// Build 6 FMH genus refs in 2 sequence scans (one per k).
std::array<FmhGenusRef, 6> build_fmh_genus_refs_sweep(
    uint64_t genus_hash,
    const std::vector<std::string_view>& genome_seqs);

// Score bin against host_refs vs contam_refs across all 6 combos in 2 contig scans.
// Returns minority fraction [0,1] per combo (NaN if no scoreable contigs).
std::array<float, 6> score_bin_fmh_sweep(
    const std::vector<std::string_view>& bin_seqs,
    const std::array<FmhGenusRef, 6>& host_refs,
    const std::array<FmhGenusRef, 6>& contam_refs,
    uint32_t min_contig_bp = 10000);

// Same as score_bin_fmh_containment but uses zero-copy FmhrView refs (mmap'd FMHR section).
// primary_genus_hash identifies which ref is "host"; remaining refs are candidates.
float score_bin_fmh_containment(
    const std::vector<std::string_view>& bin_seqs,
    uint64_t primary_genus_hash,
    const std::vector<FmhrView>& refs,
    int k, int c,
    uint32_t min_contig_bp = 10000);

// Score a bin (set of contig sequences) against GCOV genus + FCOV family models.
// gcov_entry/fcov_entry are const ptrs into memory-mapped archive data — no I/O.
// fcov_entry/fcov_rd may be nullptr; sibling_fraction will be NaN in that case.
// Thread-safe: pure computation, no mutable state.
BinScores score_bin(
    const std::vector<std::string_view>& seqs,
    const GcovEntry*      gcov_entry,
    const GcovReader*     gcov_rd,
    const GcovEntry*      fcov_entry,
    const GcovReader*     fcov_rd,
    const BinScoreConfig& cfg = {});

} // namespace genopack
