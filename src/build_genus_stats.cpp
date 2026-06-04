#include <genopack/build_genus_stats.hpp>
#include <genopack/archive.hpp>
#include <genopack/skani.hpp>
#include <genopack/util.hpp>
#include <Eigen/Eigenvalues>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <mutex>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>

namespace genopack {

static constexpr int   D           = 136;
static constexpr int   K           = static_cast<int>(GCOV_N_EIGVECS);
static constexpr float N_EFF_ALPHA = 0.25f;

using Mat136 = Eigen::Matrix<float, D, D>;
using Vec136 = Eigen::Matrix<float, D, 1>;

// GenusAccum / ProfileEntry defined in build_genus_stats.hpp (exposed for builder one-pass use).

// ── Full eigendecomposition with Ledoit-Wolf shrinkage ────────────────────────
// Builds scatter S from centred profiles, shrinks with LW oracle approximation,
// decomposes via SelfAdjointEigenSolver (exact, sub-ms for 136×136).
// Outputs top-K eigenvalues/vectors descending; sigma2_resid = mean(λ_{0..D-K-1}).

// Xw: [N×D] row-major, row i = sqrt(w_i) * (x_i - mu). S = Xw^T Xw via SYRK (BLAS-3).
static void eigen_lw(const Eigen::MatrixXf& Xw,
                     int n_eff,
                     float eigenvalues_out[K],
                     float eigvecs_out[K][D],
                     float& sigma2_resid_out)
{
    Mat136 S = Mat136::Zero();
    S.selfadjointView<Eigen::Lower>().rankUpdate(Xw.transpose());
    S = S.selfadjointView<Eigen::Lower>();  // mirror lower → full (needed for LW scalar add)

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
        float zero_rho[GCOV_RHO_DIM] = {}, zero_prec[GCOV_RHO_PREC_N] = {}, zero_rhoq[GCOV_N_QUANTILES] = {};
        writer.add(h, n_members, n_long, 0u,
                   eigenvalues, eigvecs, quantiles, p_bar, sigma2_resid, N_EFF_ALPHA, zero_mu,
                   spe_quantiles, zero_rho, zero_prec, zero_rhoq);
        return;
    }

    // p_bar: contig-weighted mean raw frequency (for delta-method sampling variance)
    float total_fn = 0.0f;
    for (const auto& pe : acc.flat) {
        for (int d = 0; d < D; ++d) p_bar[d] += pe.p[d] * pe.freq_norm;
        total_fn += pe.freq_norm;
    }
    if (total_fn > 0.0f)
        for (int d = 0; d < D; ++d) p_bar[d] /= total_fn;

    // Genome-weighted centroid: μ = (1/G) × Σ_g mean_g  (CSR layout)
    float mean[D] = {};
    const float inv_G = 1.0f / static_cast<float>(n_genomes);
    for (uint32_t g = 0; g < n_genomes; ++g) {
        const uint32_t g0 = acc.genome_offsets[g], g1 = acc.genome_offsets[g + 1];
        const float inv_ng = 1.0f / static_cast<float>(g1 - g0);
        for (uint32_t ci = g0; ci < g1; ++ci)
            for (int d = 0; d < D; ++d) mean[d] += acc.flat[ci].p[d] * inv_ng * inv_G;
    }

    // Genome-balanced SYRK: build Xw [n_long × D] with row = sqrt(w_ij) * (x_ij - mu).
    // S = Xw^T Xw  via selfadjointView::rankUpdate — BLAS-3 vs 90k rank-1 updates.
    Eigen::MatrixXf Xw(n_long, D);
    {
        uint32_t row = 0;
        for (uint32_t g = 0; g < n_genomes; ++g) {
            const uint32_t g0 = acc.genome_offsets[g], g1 = acc.genome_offsets[g + 1];
            const float sw = std::sqrt(inv_G / static_cast<float>(g1 - g0));
            for (uint32_t ci = g0; ci < g1; ++ci, ++row)
                for (int d = 0; d < D; ++d)
                    Xw(row, d) = sw * (acc.flat[ci].p[d] - mean[d]);
        }
    }

    float evecs_local[K][D];
    float evals_local[K];
    eigen_lw(Xw, static_cast<int>(n_genomes), evals_local, evecs_local, sigma2_resid);

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
    for (uint32_t ci = 0; ci < n_long; ++ci) {
        float xmu[D];
        for (int d = 0; d < D; ++d) xmu[d] = acc.flat[ci].p[d] - mean[d];
        dists[ci]     = gcov_mahalanobis(tmp_e, xmu);
        spe_dists[ci] = gcov_spe(tmp_e, xmu);
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

    // ── ρ* (Karlin genomic signature) ─────────────────────────────────────────
    // Genome-balanced mean ρ*: per-genome bp-weighted mean, then average over genomes.
    using Vec16 = Eigen::Matrix<float, 16, 1>;
    using Mat16 = Eigen::Matrix<float, 16, 16>;

    std::vector<Vec16> genome_rho(n_genomes);
    for (uint32_t g = 0; g < n_genomes; ++g) {
        genome_rho[g] = Vec16::Zero();
        uint64_t total_bp = 0;
        for (uint32_t ci = acc.genome_offsets[g]; ci < acc.genome_offsets[g + 1]; ++ci) {
            for (int i = 0; i < 16; ++i) genome_rho[g](i) += acc.flat[ci].rho[i] * acc.flat[ci].bp;
            total_bp += acc.flat[ci].bp;
        }
        if (total_bp > 0) genome_rho[g] /= static_cast<float>(total_bp);
    }

    Vec16 rho_mu = Vec16::Zero();
    for (const auto& gr : genome_rho) rho_mu += gr;
    rho_mu /= static_cast<float>(n_genomes);

    Mat16 Sigma_rho = Mat16::Zero();
    for (const auto& gr : genome_rho) {
        Vec16 z = gr - rho_mu;
        Sigma_rho += z * z.transpose();
    }
    Sigma_rho /= static_cast<float>(n_genomes);

    // LW shrinkage for 16×16
    const float alpha_rho = static_cast<float>(18) /
                            static_cast<float>((n_genomes - 1) * 18 + 2);
    Sigma_rho = (1.f - alpha_rho) * Sigma_rho
              + (alpha_rho * Sigma_rho.trace() / 16.f) * Mat16::Identity();

    Mat16 Prec_rho = Sigma_rho.inverse();

    float rho_mean_out[GCOV_RHO_DIM];
    float rho_prec_lower_out[GCOV_RHO_PREC_N];
    for (int i = 0; i < 16; ++i) rho_mean_out[i] = rho_mu(i);
    for (int i = 0; i < 16; ++i)
        for (int j = 0; j <= i; ++j)
            rho_prec_lower_out[i*(i+1)/2+j] = Prec_rho(i,j);

    // Calibration: per-contig ρ* Mahalanobis distances
    std::vector<float> rho_dists;
    rho_dists.reserve(n_long);
    for (const auto& pe : acc.flat) {
        Vec16 z;
        for (int i = 0; i < 16; ++i) z(i) = pe.rho[i] - rho_mu(i);
        float d2 = z.dot(Prec_rho * z);
        rho_dists.push_back(std::sqrt(d2 > 0.f ? d2 : 0.f));
    }
    std::sort(rho_dists.begin(), rho_dists.end());

    float rho_quantiles_out[GCOV_N_QUANTILES] = {};
    for (uint32_t q = 0; q < GCOV_N_QUANTILES; ++q) {
        const size_t idx = std::min(
            static_cast<size_t>(
                static_cast<float>(q) / static_cast<float>(GCOV_N_QUANTILES - 1)
                * static_cast<float>(rho_dists.size() - 1) + 0.5f),
            rho_dists.size() - 1);
        rho_quantiles_out[q] = rho_dists[idx];
    }

    writer.add(h, n_members, n_long, GCOV_FLAG_VALID,
               eigenvalues, eigvecs, quantiles, p_bar, sigma2_resid, N_EFF_ALPHA, mean,
               spe_quantiles, rho_mean_out, rho_prec_lower_out, rho_quantiles_out);
}

// ── Public API ────────────────────────────────────────────────────────────────

GcovEntry finalize_and_add_genus(const std::string& key, uint32_t n_members,
                                  GenusAccum& acc, GcovWriter& writer) {
    finalize_genus(key, n_members, acc, writer);
    return writer.last_entry();
}

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

GcovWriter build_family_covariance(const ArchiveReader& ar, int threads) {
    std::unordered_map<std::string, std::vector<std::string>> family_to_accs;
    std::unordered_map<std::string, uint32_t> family_member_count;

    ar.scan_taxonomy([&](std::string_view acc, std::string_view tax) {
        auto pos = tax.find("f__");
        if (pos == std::string_view::npos) return;
        auto end = tax.find(';', pos + 3);
        std::string family(tax.substr(pos, end == std::string_view::npos
                                          ? tax.size() - pos
                                          : end - pos));
        if (family == "f__") return;
        family_to_accs[family].emplace_back(acc);
    });

    if (family_to_accs.empty()) {
        spdlog::warn("build_family_covariance: no taxonomy data found in archive");
        return {};
    }

    spdlog::info("build_family_covariance: {} families to process", family_to_accs.size());

    for (const auto& [f, accs] : family_to_accs)
        family_member_count[f] = static_cast<uint32_t>(accs.size());

    GcovWriter writer;

    std::vector<std::pair<std::string, std::vector<std::string>>> families(
        family_to_accs.begin(), family_to_accs.end());
    std::sort(families.begin(), families.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    uint32_t n_valid = 0, n_invalid = 0;
    const uint32_t min_long_bp = GCOV_MIN_LONG_BP;

    for (auto& [family, accs] : families) {
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

        const uint32_t n_members = family_member_count[family];
        if (acc.n_long_contigs() >= GCOV_MIN_LONG_N) ++n_valid;
        else ++n_invalid;

        finalize_genus(family, n_members, acc, writer);
    }

    spdlog::info("build_family_covariance: {} valid, {} invalid (below min_long_n={}) families",
                 n_valid, n_invalid, GCOV_MIN_LONG_N);
    return writer;
}

GcovFcovFmhrResult build_gcov_fcov_fmhr(const ArchiveReader& ar, int threads) {
    static constexpr int K_FMH = 21, C_FMH = 125;
    static constexpr uint32_t FMH_MIN_BP  = 1000;

    // ── Step 1: scan taxonomy → per-genome (genus, family) ────────────────────
    struct Taxon { std::string genus; std::string family; };
    std::unordered_map<std::string, Taxon> acc_taxon;
    std::unordered_map<std::string, uint32_t> genus_count, family_count;

    // Exclude tombstoned genomes so the consensus is built from live genomes only.
    std::unordered_set<std::string> deleted_accs;
    ar.scan_genome_accessions([&](std::string_view acc, GenomeId gid) {
        if (ar.is_deleted(gid)) deleted_accs.insert(std::string(acc));
    });

    ar.scan_taxonomy([&](std::string_view acc, std::string_view tax) {
        if (deleted_accs.count(std::string(acc))) return;   // skip tombstoned
        Taxon t;
        auto gp = tax.find("g__");
        if (gp != std::string_view::npos) {
            auto ge = tax.find(';', gp + 3);
            std::string g(tax.substr(gp, ge == std::string_view::npos ? tax.size()-gp : ge-gp));
            if (g != "g__") { t.genus = std::move(g); ++genus_count[t.genus]; }
        }
        auto fp = tax.find("f__");
        if (fp != std::string_view::npos) {
            auto fe = tax.find(';', fp + 3);
            std::string f(tax.substr(fp, fe == std::string_view::npos ? tax.size()-fp : fe-fp));
            if (f != "f__") { t.family = std::move(f); ++family_count[t.family]; }
        }
        if (!t.genus.empty())
            acc_taxon.emplace(std::string(acc), std::move(t));
    });

    if (acc_taxon.empty()) {
        spdlog::warn("build_gcov_fcov_fmhr: no taxonomy data found");
        return {};
    }
    spdlog::info("build_gcov_fcov_fmhr: {} genera, {} families, {} genomes",
                 genus_count.size(), family_count.size(), acc_taxon.size());

    // ── Step 2: assign genus/family indices for lock-indexed vectors ───────────
    std::vector<std::string> genus_names, family_names;
    std::unordered_map<std::string, size_t> genus_idx, family_idx;
    for (const auto& [g, _] : genus_count) {
        genus_idx[g] = genus_names.size();
        genus_names.push_back(g);
    }
    for (const auto& [f, _] : family_count) {
        family_idx[f] = family_names.size();
        family_names.push_back(f);
    }
    const size_t ng = genus_names.size(), nf = family_names.size();

    std::vector<GenusAccum>             genus_accum(ng), family_accum(nf);
    std::vector<std::mutex>             genus_mu(ng),    family_mu(nf);
    std::vector<std::vector<uint64_t>>  fmhr_hashes(ng);
    std::vector<std::mutex>             fmhr_mu(ng);

    // ── Step 3: collect all accessions for a single shard scan ────────────────
    std::vector<std::string> all_accs;
    all_accs.reserve(acc_taxon.size());
    for (const auto& [acc, _] : acc_taxon) all_accs.push_back(acc);

    // ── Step 4: single shard-major pass ───────────────────────────────────────
    const int c_vals[1] = {C_FMH};
    ar.visit_shard_batches_parallel(all_accs, threads,
        [&](ArchiveReader::ShardBatch& batch) {
            for (const auto& [idx, eg] : batch) {
                if (eg.fasta.empty()) continue;
                auto tit = acc_taxon.find(all_accs[idx]);
                if (tit == acc_taxon.end()) continue;
                const auto& [genus, family] = tit->second;
                const size_t gi = genus_idx.at(genus);

                // TNF profiles for long contigs (≥ GCOV_MIN_LONG_BP)
                auto cps = compute_long_contig_profiles(eg.fasta, GCOV_MIN_LONG_BP);
                if (!cps.empty()) {
                    { std::lock_guard<std::mutex> lk(genus_mu[gi]);
                      genus_accum[gi].add_genome(cps); }
                    if (!family.empty()) {
                        const size_t fi = family_idx.at(family);
                        std::lock_guard<std::mutex> lk(family_mu[fi]);
                        family_accum[fi].add_genome(cps);
                    }
                }

                // FMH hashes for all contigs ≥ FMH_MIN_BP
                // Parse FASTA in-place, skip the multi-alloc in build_fmh_genus_section
                {
                    const char* p = eg.fasta.data(), *end = p + eg.fasta.size();
                    std::vector<uint64_t> local_hashes;
                    while (p < end) {
                        while (p < end && *p != '>') ++p;
                        if (p >= end) break;
                        while (p < end && *p != '\n') ++p;
                        if (p < end) ++p;
                        const char* seq_start = p;
                        size_t seq_len = 0;
                        while (p < end && *p != '>') {
                            while (p < end && *p != '\n' && *p != '\r') { ++seq_len; ++p; }
                            while (p < end && (*p == '\n' || *p == '\r')) ++p;
                        }
                        if (seq_len < FMH_MIN_BP) continue;
                        std::string seq;
                        seq.reserve(seq_len);
                        for (const char* sp = seq_start; sp < p; ++sp)
                            if (*sp != '\n' && *sp != '\r') seq.push_back(*sp);
                        auto vecs = fmh_multi_c(seq, K_FMH, c_vals, 1);
                        if (!vecs.empty())
                            local_hashes.insert(local_hashes.end(),
                                                vecs[0].begin(), vecs[0].end());
                    }
                    if (!local_hashes.empty()) {
                        std::lock_guard<std::mutex> lk(fmhr_mu[gi]);
                        fmhr_hashes[gi].insert(fmhr_hashes[gi].end(),
                                               local_hashes.begin(), local_hashes.end());
                    }
                }
            }
        });

    // ── Step 5: finalize all three writers ────────────────────────────────────
    GcovFcovFmhrResult result;
    uint32_t gcov_valid = 0, gcov_invalid = 0;
    for (size_t i = 0; i < ng; ++i) {
        if (genus_accum[i].n_long_contigs() >= GCOV_MIN_LONG_N) ++gcov_valid; else ++gcov_invalid;
        finalize_genus(genus_names[i], genus_count.at(genus_names[i]), genus_accum[i], result.gcov);
        auto& h = fmhr_hashes[i];
        if (!h.empty()) {
            std::sort(h.begin(), h.end());
            h.erase(std::unique(h.begin(), h.end()), h.end());
            result.fmhr.add(GcovWriter::hash_genus(genus_names[i]), std::move(h));
        }
    }
    uint32_t fcov_valid = 0, fcov_invalid = 0;
    for (size_t i = 0; i < nf; ++i) {
        if (family_accum[i].n_long_contigs() >= GCOV_MIN_LONG_N) ++fcov_valid; else ++fcov_invalid;
        finalize_genus(family_names[i], family_count.at(family_names[i]), family_accum[i], result.fcov);
    }
    spdlog::info("build_gcov_fcov_fmhr: gcov {}/{} valid genera, fcov {}/{} valid families, fmhr {} genera",
                 gcov_valid, ng, fcov_valid, nf, result.fmhr.n_genera());
    return result;
}

FmhrWriter build_fmh_genus_section(const ArchiveReader& ar, int threads) {
    static constexpr int K_FMH = 21;
    static constexpr int C_FMH = 125;
    static constexpr uint32_t MIN_CONTIG_BP = 1000;

    std::unordered_map<std::string, std::vector<std::string>> genus_to_accs;
    ar.scan_taxonomy([&](std::string_view acc, std::string_view tax) {
        auto pos = tax.find("g__");
        if (pos == std::string_view::npos) return;
        auto end = tax.find(';', pos + 3);
        std::string genus(tax.substr(pos, end == std::string_view::npos
                                        ? tax.size() - pos : end - pos));
        if (genus == "g__") return;
        genus_to_accs[genus].emplace_back(acc);
    });

    if (genus_to_accs.empty()) {
        spdlog::warn("build_fmh_genus_section: no taxonomy data found");
        return {};
    }
    spdlog::info("build_fmh_genus_section: {} genera to process", genus_to_accs.size());

    std::vector<std::pair<std::string, std::vector<std::string>>> genera(
        genus_to_accs.begin(), genus_to_accs.end());
    std::sort(genera.begin(), genera.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    FmhrWriter writer;
    const int c_vals[1] = {C_FMH};

    for (auto& [genus, accs] : genera) {
        std::vector<uint64_t> all_hashes;
        std::mutex mu;

        auto collect = [&](std::string_view fasta) {
            const char* p = fasta.data(), *end = p + fasta.size();
            // Iterate contigs in-place without storing strings; each contig is parsed
            // sequentially and hashed immediately.
            while (p < end) {
                while (p < end && *p != '>') ++p;
                if (p >= end) break;
                while (p < end && *p != '\n') ++p;
                if (p < end) ++p;
                // collect seq
                const char* seq_start = p;
                size_t seq_len = 0;
                while (p < end && *p != '>') {
                    while (p < end && *p != '\n' && *p != '\r') { ++seq_len; ++p; }
                    while (p < end && (*p == '\n' || *p == '\r')) ++p;
                }
                if (seq_len < MIN_CONTIG_BP) continue;
                // Build a contiguous seq string (needed for ntHash)
                std::string seq;
                seq.reserve(seq_len);
                const char* sp = seq_start;
                while (sp < p) {
                    if (*sp != '\n' && *sp != '\r') seq.push_back(*sp);
                    ++sp;
                }
                auto result = fmh_multi_c(seq, K_FMH, c_vals, 1);
                if (!result.empty() && !result[0].empty()) {
                    std::lock_guard<std::mutex> lk(mu);
                    all_hashes.insert(all_hashes.end(), result[0].begin(), result[0].end());
                }
            }
        };

        if (threads > 1) {
            ar.visit_shard_batches_parallel(accs, threads,
                [&](ArchiveReader::ShardBatch& batch) {
                    for (auto& [idx, eg] : batch)
                        if (!eg.fasta.empty()) collect(eg.fasta);
                });
        } else {
            ar.visit_by_shard(accs, [&](size_t, ExtractedGenome eg) {
                if (!eg.fasta.empty()) collect(eg.fasta);
            });
        }

        if (all_hashes.empty()) continue;
        std::sort(all_hashes.begin(), all_hashes.end());
        all_hashes.erase(std::unique(all_hashes.begin(), all_hashes.end()), all_hashes.end());
        writer.add(GcovWriter::hash_genus(genus), std::move(all_hashes));
    }

    spdlog::info("build_fmh_genus_section: {} genera written", writer.n_genera());
    return writer;
}

} // namespace genopack
