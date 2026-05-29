#pragma once
#include <array>
#include <vector>
#include <genopack/gcov.hpp>

namespace genopack {

struct GmmResult {
    float minority_fraction = 0.0f;  // min(piA, piB) when bimodal; 0 otherwise
    float bic_delta         = 0.0f;  // BIC(k=1) - BIC(k=2); positive → bimodal favored
    bool  bimodal           = false;
};

// 2-component diagonal-covariance GMM on whitened GCOV PCA projections.
// pts: per-contig whitened projections (xmu @ eigvecs / sqrt(lambda)).
// weights: contig lengths (bp).
// Returns minority fraction and BIC delta; bimodal = true when k=2 beats k=1 by BIC.
GmmResult fit_gmm2(
    const std::vector<std::array<float, genopack::GCOV_N_EIGVECS>>& pts,
    const std::vector<float>& weights,
    int   max_iter = 50,
    float tol      = 1e-4f);

} // namespace genopack
