#include <genopack/gcov.hpp>
#include <cmath>
#include <cstring>
#include <stdexcept>

namespace genopack {

SectionDesc GcovWriter::finalize(AppendWriter& w, uint64_t section_id, uint32_t section_type) {
    const uint32_t n = static_cast<uint32_t>(entries_.size());
    const uint32_t min_buckets = (n == 0) ? 1
        : static_cast<uint32_t>(static_cast<double>(n) / 0.7 + 1.0);
    const uint32_t n_buckets = next_pow2(min_buckets);
    const uint32_t mask      = n_buckets - 1;

    static constexpr uint32_t EMPTY = UINT32_MAX;
    std::vector<uint32_t> buckets(n_buckets, EMPTY);
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t slot = static_cast<uint32_t>(entries_[i].genus_hash) & mask;
        while (buckets[slot] != EMPTY) slot = (slot + 1) & mask;
        buckets[slot] = i;
    }

    w.align(8);
    const uint64_t section_start = w.current_offset();

    const uint64_t entries_off = sizeof(GcovHeader);
    const uint64_t buckets_off = entries_off + static_cast<uint64_t>(n) * sizeof(GcovEntry);

    GcovHeader hdr{};
    hdr.magic          = section_type;
    hdr.n_genera       = n;
    hdr.n_buckets      = n_buckets;
    hdr.n_eigvecs      = GCOV_N_EIGVECS;
    hdr.n_quantiles    = GCOV_N_QUANTILES;
    hdr.min_long_bp    = GCOV_MIN_LONG_BP;
    hdr.min_long_n     = GCOV_MIN_LONG_N;
    hdr._pad           = 0;
    hdr.entries_offset = entries_off;
    hdr.buckets_offset = buckets_off;
    hdr.entry_stride   = sizeof(GcovEntry);

    w.append(&hdr, sizeof(hdr));
    if (n > 0)
        w.append(entries_.data(), static_cast<uint64_t>(n) * sizeof(GcovEntry));
    w.append(buckets.data(), static_cast<uint64_t>(n_buckets) * sizeof(uint32_t));

    const uint64_t section_end = w.current_offset();

    SectionDesc sd{};
    sd.type              = section_type;
    sd.version           = 1;
    sd.flags             = 0;
    sd.section_id        = section_id;
    sd.file_offset       = section_start;
    sd.compressed_size   = section_end - section_start;
    sd.uncompressed_size = section_end - section_start;
    sd.item_count        = n;
    return sd;
}

float gcov_ppca_distance(const GcovEntry& e, const float* x_minus_mu, uint32_t n_kmers) noexcept {
    static constexpr int D = 136;
    static constexpr int K = static_cast<int>(GCOV_N_EIGVECS);

    if (n_kmers == 0) return 0.0f;

    // û = unit direction of genus centroid
    float mu_norm2 = 0.0f;
    for (int d = 0; d < D; ++d) mu_norm2 += e.mu[d] * e.mu[d];
    float u_hat[D] = {};
    if (mu_norm2 > 1e-12f) {
        float inv = 1.0f / std::sqrt(mu_norm2);
        for (int d = 0; d < D; ++d) u_hat[d] = e.mu[d] * inv;
    }

    float pbar_norm2 = 0.0f;
    for (int d = 0; d < D; ++d) pbar_norm2 += e.p_bar[d] * e.p_bar[d];
    if (pbar_norm2 < 1e-12f) return 0.0f;

    const float inv_neff = 1.0f / (pbar_norm2 * e.n_eff_alpha * static_cast<float>(n_kmers));

    float resid[D];
    std::memcpy(resid, x_minus_mu, D * sizeof(float));

    float d2_top = 0.0f;
    for (int k = 0; k < K; ++k) {
        if (e.eigenvalues[k] < 1e-9f) continue;

        float vk_dot_u = 0.0f;
        for (int d = 0; d < D; ++d) vk_dot_u += e.eigvecs[k][d] * u_hat[d];

        float sigma_samp_k = 0.0f;
        for (int d = 0; d < D; ++d) {
            float vk_perp_d = e.eigvecs[k][d] - vk_dot_u * u_hat[d];
            sigma_samp_k += e.p_bar[d] * vk_perp_d * vk_perp_d;
        }
        sigma_samp_k *= inv_neff;

        float proj = 0.0f;
        for (int d = 0; d < D; ++d) proj += x_minus_mu[d] * e.eigvecs[k][d];
        d2_top += (proj * proj) / (e.eigenvalues[k] + sigma_samp_k);
        for (int d = 0; d < D; ++d) resid[d] -= proj * e.eigvecs[k][d];
    }

    float tr_proj = 0.0f;
    for (int d = 0; d < D; ++d) tr_proj += e.p_bar[d] * (1.0f - u_hat[d] * u_hat[d]);
    const float sigma_samp_resid = (D - K > 0)
        ? tr_proj * inv_neff / static_cast<float>(D - K)
        : 0.0f;

    float resid_norm2 = 0.0f;
    for (int d = 0; d < D; ++d) resid_norm2 += resid[d] * resid[d];
    const float d2_resid = resid_norm2 / (e.sigma2_resid + sigma_samp_resid + 1e-12f);

    return std::sqrt(d2_top + d2_resid);
}

float gcov_mahalanobis(const GcovEntry& e, const float* x_minus_mu) noexcept {
    static constexpr int K = static_cast<int>(GCOV_N_EIGVECS);

    // Top-K only: project onto biological eigenvectors, divide by eigenvalue.
    // The residual (D-K) directions have near-zero biological variance but large
    // noise variance, so including them masks the biological signal.
    float d2 = 0.0f;
    for (int k = 0; k < K; ++k) {
        if (e.eigenvalues[k] < 1e-9f) continue;
        float proj = 0.0f;
        for (int d = 0; d < 136; ++d) proj += x_minus_mu[d] * e.eigvecs[k][d];
        d2 += (proj * proj) / e.eigenvalues[k];
    }
    return std::sqrt(d2);
}

float gcov_spe(const GcovEntry& e, const float* x_minus_mu) noexcept {
    static constexpr int K = static_cast<int>(GCOV_N_EIGVECS);
    float resid[136];
    std::memcpy(resid, x_minus_mu, 136 * sizeof(float));
    for (int k = 0; k < K; ++k) {
        float proj = 0.0f;
        for (int d = 0; d < 136; ++d) proj += x_minus_mu[d] * e.eigvecs[k][d];
        for (int d = 0; d < 136; ++d) resid[d] -= proj * e.eigvecs[k][d];
    }
    float spe = 0.0f;
    for (int d = 0; d < 136; ++d) spe += resid[d] * resid[d];
    return spe;
}

float gcov_rho_distance(const GcovEntry& e, const float* rho_minus_mean) noexcept {
    static constexpr int R = static_cast<int>(GCOV_RHO_DIM);
    // d² = z^T P z, P stored as lower triangle: P[i][j] = rho_prec_lower[i*(i+1)/2+j] for i>=j
    float d2 = 0.0f;
    for (int i = 0; i < R; ++i) {
        float row = 0.0f;
        for (int j = 0; j <= i; ++j)
            row += e.rho_prec_lower[i*(i+1)/2 + j] * rho_minus_mean[j];
        for (int j = i+1; j < R; ++j)
            row += e.rho_prec_lower[j*(j+1)/2 + i] * rho_minus_mean[j];
        d2 += rho_minus_mean[i] * row;
    }
    return std::sqrt(d2 > 0.f ? d2 : 0.f);
}

} // namespace genopack
