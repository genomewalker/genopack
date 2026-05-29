#pragma once
#include "format.hpp"
#include "mmap_file.hpp"
#include <cstdint>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace genopack {

// ── Constants ─────────────────────────────────────────────────────────────────

static constexpr uint32_t GCOV_N_EIGVECS    = 15;   // top eigenvectors of Σ_bio stored
static constexpr uint32_t GCOV_RHO_DIM     = 16;   // dinucleotide pairs (A/C/G/T × A/C/G/T)
static constexpr uint32_t GCOV_RHO_PREC_N  = GCOV_RHO_DIM*(GCOV_RHO_DIM+1)/2; // 136 lower-triangle entries
static constexpr uint32_t GCOV_N_QUANTILES  = 128;  // percentile points of Mahalanobis dist
static constexpr uint32_t GCOV_MIN_LONG_BP  = 20000; // min contig length for Σ_bio estimation
static constexpr uint32_t GCOV_MIN_LONG_N   = 20;   // min long contigs per genus

// ── On-disk layout ────────────────────────────────────────────────────────────

struct GcovHeader {             // 64 bytes
    uint32_t magic;             // SEC_GCOV
    uint32_t n_genera;
    uint32_t n_buckets;         // open-addressing table, power of 2
    uint32_t n_eigvecs;         // K (≤ GCOV_N_EIGVECS)
    uint32_t n_quantiles;       // Q (≤ GCOV_N_QUANTILES)
    uint32_t min_long_bp;       // contig length threshold used at build time
    uint32_t min_long_n;        // min long contigs required per genus
    uint32_t _pad;
    uint64_t entries_offset;    // bytes from section start → GcovEntry array
    uint64_t buckets_offset;    // bytes from section start → uint32_t bucket array
    uint64_t entry_stride;      // sizeof(GcovEntry)
    uint64_t _pad2;
};
static_assert(sizeof(GcovHeader) == 64);

// Fixed-stride per-genus entry.
// Σ_bio estimated from long contigs (≥ min_long_bp) across all clean genus members.
// Stored as top-K eigenvectors (rows of eigvecs[][136]) + eigenvalues.
// Mahalanobis distances of those long contigs are stored as Q sorted quantile points.
// flags bit 0: has valid covariance (n_long_contigs >= GCOV_MIN_LONG_N)
struct GcovEntry {
    uint64_t genus_hash;                              //  8
    uint32_t n_members;                               //  4
    uint32_t n_long_contigs;                          //  4  long contigs used for Σ_bio
    uint32_t flags;                                   //  4  bit0 = valid
    uint32_t _pad;                                    //  4
    float    eigenvalues[GCOV_N_EIGVECS];             //   60  top-K eigenvalues of LW-shrunk Σ_bio (descending)
    float    eigvecs[GCOV_N_EIGVECS][136];            // 8160  top-K eigenvectors row-major (descending)
    float    quantiles[GCOV_N_QUANTILES];             //  512  sorted Mahalanobis distances at Q equally-spaced percentile points
    float    p_bar[136];                              //  544  raw L1-normalised genus mean frequency
    float    sigma2_resid;                            //    4  PPCA isotropic floor = mean(λ_{k>K})
    float    n_eff_alpha;                             //    4  N_eff = n_eff_alpha * n_kmers (default 0.25)
    float    mu[136];                                 //  544  L2-normalised genus centroid (mean of long-contig profiles)
    uint32_t _pad2;                                   //    4
    float    spe_quantiles[GCOV_N_QUANTILES];         //  512  sorted SPE (‖residual‖²) values at Q percentile points
    float    rho_mean[GCOV_RHO_DIM];                  //   64  genus mean dinucleotide ρ* vector
    float    rho_prec_lower[GCOV_RHO_PREC_N];         //  544  precision matrix Σ_ρ⁻¹ lower triangle
    float    rho_quantiles[GCOV_N_QUANTILES];          //  512  sorted ρ* Mahalanobis distances for percentile lookup
};
static_assert(sizeof(GcovEntry) ==
    8 + 4 + 4 + 4 + 4 +
    GCOV_N_EIGVECS * 4 +
    GCOV_N_EIGVECS * 136 * 4 +
    GCOV_N_QUANTILES * 4 +
    136 * 4 + 4 + 4 + 136 * 4 + 4 +
    GCOV_N_QUANTILES * 4 +
    GCOV_RHO_DIM * 4 + GCOV_RHO_PREC_N * 4 + GCOV_N_QUANTILES * 4);  // rho fields

static constexpr uint32_t GCOV_FLAG_VALID = 0x1u;

// ── Writer ────────────────────────────────────────────────────────────────────

class GcovWriter {
public:
    // Add one genus. eigvecs[K][136] row-major; eigenvalues[K]; quantiles[Q] sorted ascending.
    // Pass flags=0 (invalid) for genera below the long-contig threshold.
    void add(uint64_t    genus_hash,
             uint32_t    n_members,
             uint32_t    n_long_contigs,
             uint32_t    flags,
             const float eigenvalues[GCOV_N_EIGVECS],
             const float eigvecs[GCOV_N_EIGVECS][136],
             const float quantiles[GCOV_N_QUANTILES],
             const float p_bar[136],
             float       sigma2_resid,
             float       n_eff_alpha,
             const float mu[136],
             const float spe_quantiles[GCOV_N_QUANTILES],
             const float rho_mean[GCOV_RHO_DIM],
             const float rho_prec_lower[GCOV_RHO_PREC_N],
             const float rho_quantiles[GCOV_N_QUANTILES])
    {
        GcovEntry e{};
        e.genus_hash      = genus_hash;
        e.n_members       = n_members;
        e.n_long_contigs  = n_long_contigs;
        e.flags           = flags;
        for (uint32_t k = 0; k < GCOV_N_EIGVECS; ++k) {
            e.eigenvalues[k] = eigenvalues[k];
            for (int d = 0; d < 136; ++d) e.eigvecs[k][d] = eigvecs[k][d];
        }
        for (uint32_t q = 0; q < GCOV_N_QUANTILES; ++q)
            e.quantiles[q] = quantiles[q];
        for (int d = 0; d < 136; ++d) e.p_bar[d] = p_bar[d];
        e.sigma2_resid = sigma2_resid;
        e.n_eff_alpha  = n_eff_alpha;
        for (int d = 0; d < 136; ++d) e.mu[d] = mu[d];
        for (uint32_t q = 0; q < GCOV_N_QUANTILES; ++q)
            e.spe_quantiles[q] = spe_quantiles[q];
        for (uint32_t i = 0; i < GCOV_RHO_DIM; ++i) e.rho_mean[i] = rho_mean[i];
        for (uint32_t i = 0; i < GCOV_RHO_PREC_N; ++i) e.rho_prec_lower[i] = rho_prec_lower[i];
        for (uint32_t q = 0; q < GCOV_N_QUANTILES; ++q) e.rho_quantiles[q] = rho_quantiles[q];
        entries_.push_back(e);
    }

    SectionDesc finalize(AppendWriter& w, uint64_t section_id, uint32_t section_type = SEC_GCOV);
    size_t n_genera() const { return entries_.size(); }
    const GcovEntry& last_entry() const { return entries_.back(); }

    // FNV-1a 64-bit hash of genus name; 0 remapped to 1 (0 = empty sentinel).
    static uint64_t hash_genus(std::string_view s) noexcept {
        uint64_t h = 0xcbf29ce484222325ULL;
        for (unsigned char c : s) { h ^= c; h *= 0x100000001b3ULL; }
        return h == 0 ? 1 : h;
    }

private:
    static uint32_t next_pow2(uint32_t v) noexcept {
        if (v == 0) return 1;
        --v;
        v |= v>>1; v |= v>>2; v |= v>>4; v |= v>>8; v |= v>>16;
        return v + 1;
    }
    std::vector<GcovEntry> entries_;
};

// ── Reader ────────────────────────────────────────────────────────────────────

class GcovReader {
public:
    void open(const uint8_t* data, uint64_t offset, uint64_t size) {
        if (size < sizeof(GcovHeader))
            throw std::runtime_error("GCOV section too small");
        data_   = data + offset;
        header_ = reinterpret_cast<const GcovHeader*>(data_);
        if (header_->magic != SEC_GCOV && header_->magic != SEC_FCOV)
            throw std::runtime_error("GCOV/FCOV: bad magic");
        if (header_->entry_stride != sizeof(GcovEntry))
            throw std::runtime_error("GCOV: unknown entry_stride — rebuild required");
        entries_ = reinterpret_cast<const GcovEntry*>(data_ + header_->entries_offset);
        buckets_ = reinterpret_cast<const uint32_t*>(data_ + header_->buckets_offset);
        n_buckets_ = header_->n_buckets;
    }

    bool     is_open()    const { return data_ != nullptr; }
    uint32_t n_genera()   const { return header_ ? header_->n_genera    : 0; }
    uint32_t n_eigvecs()  const { return header_ ? header_->n_eigvecs   : GCOV_N_EIGVECS; }
    uint32_t n_quantiles()const { return header_ ? header_->n_quantiles : GCOV_N_QUANTILES; }
    uint32_t min_long_bp()const { return header_ ? header_->min_long_bp : GCOV_MIN_LONG_BP; }

    // Returns nullptr if genus not found or section not open.
    const GcovEntry* lookup(uint64_t genus_hash) const {
        if (!data_ || n_buckets_ == 0) return nullptr;
        const uint32_t mask = n_buckets_ - 1;
        uint32_t slot = static_cast<uint32_t>(genus_hash) & mask;
        static constexpr uint32_t EMPTY = UINT32_MAX;
        for (uint32_t tries = 0; tries < n_buckets_; ++tries) {
            uint32_t idx = buckets_[slot];
            if (idx == EMPTY) return nullptr;
            if (entries_[idx].genus_hash == genus_hash) return &entries_[idx];
            slot = (slot + 1) & mask;
        }
        return nullptr;
    }

    // Iterate all valid entries.
    void scan(const std::function<void(const GcovEntry&)>& cb) const {
        if (!data_) return;
        for (uint32_t i = 0; i < header_->n_genera; ++i)
            cb(entries_[i]);
    }

    // Fraction of calibration contigs with Mahalanobis distance <= d.
    float percentile(const GcovEntry& e, float d) const {
        const uint32_t Q = header_->n_quantiles;
        if (Q == 0) return NAN;
        uint32_t lo = 0, hi = Q;
        while (lo < hi) {
            uint32_t mid = (lo + hi) / 2;
            if (e.quantiles[mid] <= d) lo = mid + 1; else hi = mid;
        }
        return static_cast<float>(lo) / static_cast<float>(Q);
    }

    // Fraction of calibration contigs with ρ* Mahalanobis distance <= d.
    float rho_percentile(const GcovEntry& e, float d) const {
        const uint32_t Q = header_->n_quantiles;
        if (Q == 0) return NAN;
        uint32_t lo = 0, hi = Q;
        while (lo < hi) {
            uint32_t mid = (lo + hi) / 2;
            if (e.rho_quantiles[mid] <= d) lo = mid + 1; else hi = mid;
        }
        return static_cast<float>(lo) / static_cast<float>(Q);
    }

    // Fraction of calibration contigs with SPE (‖residual‖²) <= spe.
    float spe_percentile(const GcovEntry& e, float spe) const {
        const uint32_t Q = header_->n_quantiles;
        if (Q == 0) return NAN;
        uint32_t lo = 0, hi = Q;
        while (lo < hi) {
            uint32_t mid = (lo + hi) / 2;
            if (e.spe_quantiles[mid] <= spe) lo = mid + 1; else hi = mid;
        }
        return static_cast<float>(lo) / static_cast<float>(Q);
    }

private:
    const uint8_t*    data_      = nullptr;
    const GcovHeader* header_    = nullptr;
    const GcovEntry*  entries_   = nullptr;
    const uint32_t*   buckets_   = nullptr;
    uint32_t          n_buckets_ = 0;
};

// PPCA+Δ-method distance (length-adjusted; use for statistical p-value estimation).
// x_minus_mu = contig_tnf[136] - e.mu[136]; n_kmers = contig_bp - 3.
// Standalone percentile lookups — no GcovReader required; use compile-time GCOV_N_QUANTILES.
inline float gcov_percentile(const GcovEntry& e, float d) noexcept {
    constexpr uint32_t Q = GCOV_N_QUANTILES;
    uint32_t lo = 0, hi = Q;
    while (lo < hi) { uint32_t mid=(lo+hi)/2; if (e.quantiles[mid]<=d) lo=mid+1; else hi=mid; }
    return static_cast<float>(lo) / static_cast<float>(Q);
}
inline float gcov_spe_percentile(const GcovEntry& e, float spe) noexcept {
    constexpr uint32_t Q = GCOV_N_QUANTILES;
    uint32_t lo = 0, hi = Q;
    while (lo < hi) { uint32_t mid=(lo+hi)/2; if (e.spe_quantiles[mid]<=spe) lo=mid+1; else hi=mid; }
    return static_cast<float>(lo) / static_cast<float>(Q);
}
inline float gcov_rho_percentile(const GcovEntry& e, float d) noexcept {
    constexpr uint32_t Q = GCOV_N_QUANTILES;
    uint32_t lo = 0, hi = Q;
    while (lo < hi) { uint32_t mid=(lo+hi)/2; if (e.rho_quantiles[mid]<=d) lo=mid+1; else hi=mid; }
    return static_cast<float>(lo) / static_cast<float>(Q);
}

float gcov_ppca_distance(const GcovEntry& e, const float* x_minus_mu, uint32_t n_kmers) noexcept;

// Pure Mahalanobis distance in top-K eigenbasis (no length correction).
float gcov_mahalanobis(const GcovEntry& e, const float* x_minus_mu) noexcept;

// SPE (squared prediction error): ‖x_minus_mu − Σ_k proj_k v_k‖²
// Catches contamination whose direction lies outside the top-K biological subspace.
float gcov_spe(const GcovEntry& e, const float* x_minus_mu) noexcept;

// Mahalanobis distance in 16-dim ρ* space using stored precision matrix.
// rho_minus_mean = contig_rho[16] - e.rho_mean[16].
float gcov_rho_distance(const GcovEntry& e, const float* rho_minus_mean) noexcept;

} // namespace genopack
