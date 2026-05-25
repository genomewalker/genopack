#include <genopack/build_genus_stats.hpp>
#include <genopack/archive.hpp>
#include <genopack/util.hpp>
#include <Eigen/Eigenvalues>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace genopack {

static constexpr int   D           = 136;
static constexpr int   K           = static_cast<int>(GCOV_N_EIGVECS);
static constexpr float N_EFF_ALPHA = 0.25f;

using Mat136 = Eigen::Matrix<float, D, D>;
using Vec136 = Eigen::Matrix<float, D, 1>;

// ── Genus accumulator ─────────────────────────────────────────────────────────

struct ProfileEntry {
    std::array<float, D> p;
    float freq_norm;
    uint32_t bp;
};

struct GenusAccum {
    std::vector<std::vector<ProfileEntry>> genomes;  // one inner vec per genome

    void add_genome(const std::vector<ContigProfile>& cps) {
        if (cps.empty()) return;
        auto& g = genomes.emplace_back();
        g.reserve(cps.size());
        for (const auto& cp : cps) {
            std::array<float, D> a;
            std::memcpy(a.data(), cp.p, D * sizeof(float));
            g.push_back({a, cp.freq_norm, cp.bp});
        }
    }

    uint32_t n_genomes_with_long() const { return static_cast<uint32_t>(genomes.size()); }

    uint32_t n_long_contigs() const {
        uint32_t n = 0;
        for (const auto& g : genomes) n += static_cast<uint32_t>(g.size());
        return n;
    }
};

// ── Full eigendecomposition with Ledoit-Wolf shrinkage ────────────────────────
// Builds scatter S from centred profiles, shrinks with LW oracle approximation,
// decomposes via SelfAdjointEigenSolver (exact, sub-ms for 136×136).
// Outputs top-K eigenvalues/vectors descending; sigma2_resid = mean(λ_{0..D-K-1}).

static void eigen_lw(const std::vector<std::array<float, D>>& X_c,
                     const std::vector<float>& weights,
                     int n_eff,
                     float eigenvalues_out[K],
                     float eigvecs_out[K][D],
                     float& sigma2_resid_out)
{
    Mat136 S = Mat136::Zero();
    for (size_t i = 0; i < X_c.size(); ++i) {
        Vec136 xi = Eigen::Map<const Vec136>(X_c[i].data());
        S += weights[i] * xi * xi.transpose();
    }

    // LW oracle using n_eff = n_genomes (independent units, not total contigs)
    const float alpha = static_cast<float>(D + 2) /
                        static_cast<float>((n_eff - 1) * (D + 2) + 2);
    const float tr    = S.trace();
    S = (1.0f - alpha) * S + (alpha * tr / static_cast<float>(D)) * Mat136::Identity();

    Eigen::SelfAdjointEigenSolver<Mat136> eig(S, Eigen::ComputeEigenvectors);

    float resid_sum = 0.0f;
    for (int k = 0; k < D - K; ++k)
        resid_sum += eig.eigenvalues()(k);
    sigma2_resid_out = (D - K > 0) ? resid_sum / static_cast<float>(D - K) : 0.0f;

    for (int k = 0; k < K; ++k) {
        eigenvalues_out[k] = eig.eigenvalues()(D - 1 - k);
        for (int d = 0; d < D; ++d)
            eigvecs_out[k][d] = eig.eigenvectors()(d, D - 1 - k);
    }
}

// ── Per-genus finalization ────────────────────────────────────────────────────

static void finalize_genus(const std::string& genus_key,
                            uint32_t           n_members,
                            GenusAccum&        acc,
                            GcovWriter&        writer)
{
    const uint32_t n_genomes = acc.n_genomes_with_long();
    const uint32_t n_long    = acc.n_long_contigs();
    const uint64_t h         = GcovWriter::hash_genus(genus_key);

    float eigenvalues[GCOV_N_EIGVECS] = {};
    float eigvecs[GCOV_N_EIGVECS][D]  = {};
    float quantiles[GCOV_N_QUANTILES] = {};
    float spe_quantiles[GCOV_N_QUANTILES] = {};
    float p_bar[D]                    = {};
    float sigma2_resid                = 0.0f;

    if (n_genomes < GCOV_MIN_LONG_N) {
        float zero_mu[D] = {};
        writer.add(h, n_members, n_long, 0u,
                   eigenvalues, eigvecs, quantiles, p_bar, sigma2_resid, N_EFF_ALPHA, zero_mu,
                   spe_quantiles);
        return;
    }

    // p_bar: contig-weighted mean raw frequency (for delta-method sampling variance)
    float total_fn = 0.0f;
    for (const auto& gp : acc.genomes)
        for (const auto& pe : gp) {
            for (int d = 0; d < D; ++d) p_bar[d] += pe.p[d] * pe.freq_norm;
            total_fn += pe.freq_norm;
        }
    if (total_fn > 0.0f)
        for (int d = 0; d < D; ++d) p_bar[d] /= total_fn;

    // Genome-weighted centroid: μ = (1/G) × Σ_g mean_g
    float mean[D] = {};
    const float inv_G = 1.0f / static_cast<float>(n_genomes);
    for (const auto& gp : acc.genomes) {
        const float inv_ng = 1.0f / static_cast<float>(gp.size());
        for (const auto& pe : gp)
            for (int d = 0; d < D; ++d) mean[d] += pe.p[d] * inv_ng * inv_G;
    }

    // Genome-balanced contig deviations: w_ij = 1/(n_g × G)
    // Each genome contributes total weight 1/G regardless of fragmentation.
    std::vector<std::array<float, D>> X_c;
    std::vector<float> weights;
    X_c.reserve(n_long);
    weights.reserve(n_long);
    for (const auto& gp : acc.genomes) {
        const float w = inv_G / static_cast<float>(gp.size());
        for (const auto& pe : gp) {
            std::array<float, D> xc{};
            for (int d = 0; d < D; ++d) xc[d] = pe.p[d] - mean[d];
            X_c.push_back(xc);
            weights.push_back(w);
        }
    }

    float evecs_local[K][D];
    float evals_local[K];
    eigen_lw(X_c, weights, static_cast<int>(n_genomes), evals_local, evecs_local, sigma2_resid);

    for (int k = 0; k < K; ++k) {
        eigenvalues[k] = evals_local[k];
        for (int d = 0; d < D; ++d) eigvecs[k][d] = evecs_local[k][d];
    }

    // Build tmp_e for calibration distance functions
    GcovEntry tmp_e{};
    for (int k = 0; k < K; ++k) {
        tmp_e.eigenvalues[k] = evals_local[k];
        for (int d = 0; d < D; ++d) tmp_e.eigvecs[k][d] = evecs_local[k][d];
    }
    for (int d = 0; d < D; ++d) { tmp_e.p_bar[d] = p_bar[d]; tmp_e.mu[d] = mean[d]; }
    tmp_e.sigma2_resid = sigma2_resid;

    // T² and SPE calibration — both on individual contigs (matches query time)
    std::vector<float> dists(n_long);
    std::vector<float> spe_dists(n_long);
    uint32_t ci = 0;
    for (const auto& gp : acc.genomes) {
        for (const auto& pe : gp) {
            float xmu[D];
            for (int d = 0; d < D; ++d) xmu[d] = pe.p[d] - mean[d];
            dists[ci]     = gcov_mahalanobis(tmp_e, xmu);
            spe_dists[ci] = gcov_spe(tmp_e, xmu);
            ++ci;
        }
    }
    std::sort(dists.begin(), dists.end());
    std::sort(spe_dists.begin(), spe_dists.end());

    for (uint32_t q = 0; q < GCOV_N_QUANTILES; ++q) {
        const size_t idx = std::min(
            static_cast<size_t>(
                static_cast<float>(q) / static_cast<float>(GCOV_N_QUANTILES - 1)
                * static_cast<float>(n_long - 1) + 0.5f),
            static_cast<size_t>(n_long - 1));
        quantiles[q]     = dists[idx];
        spe_quantiles[q] = spe_dists[idx];
    }

    writer.add(h, n_members, n_long, GCOV_FLAG_VALID,
               eigenvalues, eigvecs, quantiles, p_bar, sigma2_resid, N_EFF_ALPHA, mean,
               spe_quantiles);
}

// ── Public API ────────────────────────────────────────────────────────────────

GcovWriter build_genus_covariance(const ArchiveReader& ar, int threads) {
    std::unordered_map<std::string, std::vector<std::string>> genus_to_accs;
    std::unordered_map<std::string, uint32_t> genus_member_count;

    ar.scan_taxonomy([&](std::string_view acc, std::string_view tax) {
        auto pos = tax.find("g__");
        if (pos == std::string_view::npos) return;
        auto end = tax.find(';', pos + 3);
        std::string genus(tax.substr(pos, end == std::string_view::npos
                                        ? tax.size() - pos
                                        : end - pos));
        if (genus == "g__") return;
        genus_to_accs[genus].emplace_back(acc);
    });

    if (genus_to_accs.empty()) {
        spdlog::warn("build_genus_covariance: no taxonomy data found in archive");
        return {};
    }

    spdlog::info("build_genus_covariance: {} genera to process", genus_to_accs.size());

    for (const auto& [g, accs] : genus_to_accs)
        genus_member_count[g] = static_cast<uint32_t>(accs.size());

    GcovWriter writer;

    std::vector<std::pair<std::string, std::vector<std::string>>> genera(
        genus_to_accs.begin(), genus_to_accs.end());
    std::sort(genera.begin(), genera.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    uint32_t n_valid = 0, n_invalid = 0;
    const uint32_t min_long_bp = GCOV_MIN_LONG_BP;

    for (auto& [genus, accs] : genera) {
        GenusAccum acc;

        if (threads > 1) {
            std::mutex mu;
            ar.visit_shard_batches_parallel(accs, threads,
                [&](ArchiveReader::ShardBatch& batch) {
                    std::vector<std::vector<ContigProfile>> local;
                    local.reserve(batch.size());
                    for (auto& [idx, eg] : batch)
                        local.push_back(compute_long_contig_profiles(eg.fasta, min_long_bp));
                    std::lock_guard<std::mutex> lk(mu);
                    for (auto& cps : local) acc.add_genome(cps);
                });
        } else {
            ar.visit_by_shard(accs, [&](size_t, ExtractedGenome eg) {
                auto cps = compute_long_contig_profiles(eg.fasta, min_long_bp);
                acc.add_genome(cps);
            });
        }

        const uint32_t n_members = genus_member_count[genus];
        if (acc.n_long_contigs() >= GCOV_MIN_LONG_N) ++n_valid;
        else ++n_invalid;

        finalize_genus(genus, n_members, acc, writer);
    }

    spdlog::info("build_genus_covariance: {} valid, {} invalid (below min_long_n={}) genera",
                 n_valid, n_invalid, GCOV_MIN_LONG_N);
    return writer;
}

} // namespace genopack
