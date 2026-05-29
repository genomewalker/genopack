#include "gmm.hpp"
#include <algorithm>
#include <cmath>

namespace genopack {

static constexpr int GK = static_cast<int>(genopack::GCOV_N_EIGVECS);

// Log-likelihood of one point under a diagonal Gaussian (omits the constant -GK/2 * log(2π)
// which cancels in BIC delta comparisons between k=1 and k=2).
static float log_gauss_diag(const float* x, const float* mu, const float* var) noexcept {
    float ll = 0.f;
    for (int k = 0; k < GK; ++k) {
        float d = x[k] - mu[k];
        ll -= 0.5f * (std::log(var[k]) + d * d / var[k]);
    }
    return ll;
}

GmmResult fit_gmm2(
    const std::vector<std::array<float, genopack::GCOV_N_EIGVECS>>& pts,
    const std::vector<float>& weights,
    int   max_iter,
    float tol)
{
    const int N = static_cast<int>(pts.size());
    if (N < 10) return {};

    float W = 0.f;
    for (float w : weights) W += w;
    if (W < 1e-12f) return {};

    // ── k=1 fit ───────────────────────────────────────────────────────────────
    float mu1[GK] = {}, var1[GK] = {};
    for (int i = 0; i < N; ++i) {
        float w = weights[i] / W;
        for (int k = 0; k < GK; ++k) mu1[k] += w * pts[i][k];
    }
    for (int i = 0; i < N; ++i) {
        float w = weights[i] / W;
        for (int k = 0; k < GK; ++k) {
            float d = pts[i][k] - mu1[k];
            var1[k] += w * d * d;
        }
    }
    for (int k = 0; k < GK; ++k) var1[k] = std::max(var1[k], 1e-8f);

    float ll1 = 0.f;
    for (int i = 0; i < N; ++i)
        ll1 += weights[i] * log_gauss_diag(pts[i].data(), mu1, var1);

    // ── k=2 initialisation: split along highest-variance dimension ─────────────
    int sdim = 0;
    for (int k = 1; k < GK; ++k)
        if (var1[k] > var1[sdim]) sdim = k;
    float halfstep = 0.5f * std::sqrt(var1[sdim]);

    float muA[GK], muB[GK], varA[GK], varB[GK];
    for (int k = 0; k < GK; ++k) {
        muA[k]  = mu1[k] + (k == sdim ?  halfstep : 0.f);
        muB[k]  = mu1[k] + (k == sdim ? -halfstep : 0.f);
        varA[k] = varB[k] = var1[k];
    }
    float piA = 0.5f, piB = 0.5f;

    std::vector<float> rA(N);
    float prev_ll2 = -1e30f;

    // ── EM iterations ──────────────────────────────────────────────────────────
    for (int iter = 0; iter < max_iter; ++iter) {
        float ll2 = 0.f;
        for (int i = 0; i < N; ++i) {
            float la = std::log(piA) + log_gauss_diag(pts[i].data(), muA, varA);
            float lb = std::log(piB) + log_gauss_diag(pts[i].data(), muB, varB);
            float lmax  = std::max(la, lb);
            float denom = lmax + std::log(std::exp(la - lmax) + std::exp(lb - lmax));
            rA[i] = std::exp(la - denom);
            ll2  += weights[i] * denom;
        }

        float WA = 0.f, WB = 0.f;
        for (int i = 0; i < N; ++i) {
            WA += weights[i] * rA[i];
            WB += weights[i] * (1.f - rA[i]);
        }
        // Prevent degenerate collapse (< 5% in either component)
        if (WA / W < 0.05f || WB / W < 0.05f) { prev_ll2 = ll2; break; }

        for (int k = 0; k < GK; ++k) muA[k] = muB[k] = 0.f;
        for (int i = 0; i < N; ++i) {
            float wa = weights[i] * rA[i]          / WA;
            float wb = weights[i] * (1.f - rA[i])  / WB;
            for (int k = 0; k < GK; ++k) {
                muA[k] += wa * pts[i][k];
                muB[k] += wb * pts[i][k];
            }
        }
        for (int k = 0; k < GK; ++k) varA[k] = varB[k] = 0.f;
        for (int i = 0; i < N; ++i) {
            float wa = weights[i] * rA[i]          / WA;
            float wb = weights[i] * (1.f - rA[i])  / WB;
            for (int k = 0; k < GK; ++k) {
                float da = pts[i][k] - muA[k];
                float db = pts[i][k] - muB[k];
                varA[k] += wa * da * da;
                varB[k] += wb * db * db;
            }
        }
        for (int k = 0; k < GK; ++k) {
            varA[k] = std::max(varA[k], 1e-8f);
            varB[k] = std::max(varB[k], 1e-8f);
        }
        piA = WA / W;
        piB = WB / W;

        if (ll2 - prev_ll2 < tol) { prev_ll2 = ll2; break; }
        prev_ll2 = ll2;
    }

    // ── BIC comparison ─────────────────────────────────────────────────────────
    // n_params: k=1 → 2*GK (means + variances); k=2 → 4*GK + 1 (+ mixing proportion)
    float logN      = std::log(static_cast<float>(N));
    float bic1      = -2.f * ll1      + 2.f * GK           * logN;
    float bic2      = -2.f * prev_ll2 + (4.f * GK + 1.f)   * logN;
    float bic_delta = bic1 - bic2;

    bool  bimodal  = bic_delta > 0.f;
    float minority = bimodal ? std::min(piA, piB) : 0.f;
    return {minority, bic_delta, bimodal};
}

} // namespace genopack
