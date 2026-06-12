#include <genopack/score_bin.hpp>
#include <genopack/gcov.hpp>
#include <genopack/oph_sketch.hpp>
#include <genopack/skani.hpp>
#include "tnf.hpp"
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <numeric>
#include <vector>

namespace genopack {

// OPH Jaccard between two contig sketches, using union of real bins as denominator.
static float oph_jaccard_contigs(const OPHSketchResult& a, const OPHSketchResult& b,
                                  uint32_t sketch_sz) noexcept {
    uint32_t match = 0, union_cnt = 0;
    for (uint32_t b_ = 0; b_ < sketch_sz; ++b_) {
        const bool ra = !!(a.real_bins_bitmask[b_/64] & (1ULL << (b_%64)));
        const bool rb = !!(b.real_bins_bitmask[b_/64] & (1ULL << (b_%64)));
        if (!ra && !rb) continue;
        ++union_cnt;
        if (ra && rb &&
            a.signature[b_] == b.signature[b_]) ++match;
    }
    return union_cnt ? static_cast<float>(match) / static_cast<float>(union_cnt) : 0.f;
}

ContainmentSplitResult score_bin_containment(
    const std::vector<std::string_view>& seqs,
    uint64_t primary_genus_hash,
    const std::vector<const GstxEntry*>& candidate_genera,
    const GstxReader& gstx_rd,
    uint32_t min_contig_bp)
{
    if (!gstx_rd.is_open()) return {};

    // Use highest stored k for best species discrimination.
    const uint32_t n_k      = gstx_rd.n_k();
    const uint32_t kmer_size = gstx_rd.kmer_size(n_k > 0 ? n_k - 1 : 0);
    // GSTX consensus[ki] index for the chosen k (last slot).
    const uint32_t ki_ref    = n_k > 0 ? n_k - 1 : 0;
    const uint32_t sketch_sz = gstx_rd.sketch_size();
    if (kmer_size == 0 || sketch_sz == 0) return {};

    OPHSketchConfig cfg;
    cfg.kmer_size   = static_cast<int>(kmer_size);
    cfg.sketch_size = static_cast<int>(sketch_sz);
    cfg.syncmer_s   = 0;
    cfg.seed        = 42;

    // Build per-contig sketches + TNF vectors for all long-enough sequences.
    struct ContigSketch { OPHSketchResult sk; Eigen::VectorXf tnf; uint32_t len; bool has_tnf = false; };
    std::vector<ContigSketch> csk;
    csk.reserve(seqs.size());
    for (const auto& seq : seqs) {
        if (seq.size() < min_contig_bp) continue;
        auto sk = sketch_oph_from_buffer(seq.data(), seq.size(), cfg);
        if (sk.n_real_bins == 0) continue;
        Eigen::VectorXf tnf;
        bool ht = compute_tnf(seq, tnf);
        csk.push_back({std::move(sk), std::move(tnf), static_cast<uint32_t>(seq.size()), ht});
    }
    if (csk.empty()) return {};

    ContainmentSplitResult r;
    const uint32_t N = static_cast<uint32_t>(csk.size());

    // ── 1. Reference minority split ───────────────────────────────────────────
    if (candidate_genera.size() >= 2) {
        uint64_t minority_bp = 0, scored_bp = 0;
        for (const auto& cs : csk) {
            scored_bp += cs.len;
            uint64_t best_hash = 0;
            float    best_cont = -1.f;
            for (const auto* e : candidate_genera) {
                if (!e || e->genus_hash == 0 || e->n_k_stored == 0) continue;
                const uint32_t ki = std::min(ki_ref, static_cast<uint32_t>(e->n_k_stored - 1));
                uint32_t matches = 0;
                for (uint32_t b = 0; b < sketch_sz; ++b) {
                    if (!(cs.sk.real_bins_bitmask[b/64] & (1ULL << (b%64)))) continue;
                    if (cs.sk.signature[b] == e->consensus[ki][b]) ++matches;
                }
                float cont = static_cast<float>(matches) / static_cast<float>(cs.sk.n_real_bins);
                if (cont > best_cont) { best_cont = cont; best_hash = e->genus_hash; }
            }
            if (best_hash != 0 && best_hash != primary_genus_hash)
                minority_bp += cs.len;
        }
        if (scored_bp > 0)
            r.minority_fraction = static_cast<float>(minority_bp) / static_cast<float>(scored_bp);
    }

    // ── 2. Self-containment null (reference-free) ────────────────────────────
    // Bin consensus via Boyer-Moore majority vote on real bins.
    std::vector<uint16_t> consensus(sketch_sz, 0);
    std::vector<int32_t>  votes(sketch_sz, 0);
    for (const auto& cs : csk) {
        for (uint32_t b = 0; b < sketch_sz; ++b) {
            if (!(cs.sk.real_bins_bitmask[b/64] & (1ULL << (b%64)))) continue;
            const uint16_t v = cs.sk.signature[b];
            if (votes[b] == 0) { consensus[b] = v; votes[b] = 1; }
            else if (v == consensus[b]) ++votes[b];
            else                        --votes[b];
        }
    }
    // Per-contig containment vs bin consensus.
    std::vector<float> self_cont(N);
    for (uint32_t i = 0; i < N; ++i) {
        uint32_t match = 0;
        for (uint32_t b = 0; b < sketch_sz; ++b) {
            if (!(csk[i].sk.real_bins_bitmask[b/64] & (1ULL << (b%64)))) continue;
            if (csk[i].sk.signature[b] == consensus[b]) ++match;
        }
        self_cont[i] = static_cast<float>(match) / static_cast<float>(csk[i].sk.n_real_bins);
    }
    if (N >= 3) {
        float mu = 0.f;
        for (float v : self_cont) mu += v;
        mu /= static_cast<float>(N);
        float var = 0.f;
        for (float v : self_cont) var += (v - mu) * (v - mu);
        const float sd = std::sqrt(var / static_cast<float>(N));
        const float thresh = mu - 2.f * sd;
        uint64_t out_bp = 0, tot_bp = 0;
        for (uint32_t i = 0; i < N; ++i) {
            tot_bp += csk[i].len;
            if (self_cont[i] < thresh) out_bp += csk[i].len;
        }
        if (tot_bp > 0)
            r.self_outlier_fraction = static_cast<float>(out_bp) / static_cast<float>(tot_bp);
    }

    // ── 3. Sketch Fiedler (spectral gap of pairwise Jaccard graph) ───────────
    static constexpr uint32_t k_max_fiedler = 256;
    if (N >= 3) {
        // Cap to longest k_max_fiedler contigs to bound O(N³) eigensolver cost.
        std::vector<uint32_t> idx(N);
        std::iota(idx.begin(), idx.end(), 0u);
        const uint32_t Nf = std::min(N, k_max_fiedler);
        if (N > Nf) {
            std::partial_sort(idx.begin(), idx.begin() + Nf, idx.end(),
                              [&](uint32_t a, uint32_t b){ return csk[a].len > csk[b].len; });
            idx.resize(Nf);
        }
        Eigen::MatrixXf W(Nf, Nf);
        for (uint32_t i = 0; i < Nf; ++i) {
            W(i, i) = 0.f;
            for (uint32_t j = i + 1; j < Nf; ++j) {
                float J = oph_jaccard_contigs(csk[idx[i]].sk, csk[idx[j]].sk, sketch_sz);
                W(i, j) = W(j, i) = J;
            }
        }
        Eigen::VectorXf deg = W.rowwise().sum();
        Eigen::MatrixXf L   = Eigen::MatrixXf::Identity(Nf, Nf);
        for (uint32_t i = 0; i < Nf; ++i) {
            if (deg(i) < 1e-9f) continue;
            for (uint32_t j = 0; j < Nf; ++j) {
                if (i == j) continue;
                const float dij = deg(i) * deg(j);
                if (dij > 0.f) L(i, j) = -W(i, j) / std::sqrt(dij);
            }
        }
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(L, Eigen::EigenvaluesOnly);
        if (es.info() == Eigen::Success) {
            const float lam2 = es.eigenvalues()(1);
            r.fiedler_oph = 1.f - std::min(1.f, std::max(0.f, lam2 / 2.f));
        }
    }

    // ── 4. TNF-kernel Fiedler: v₂ bimodality + λ₃−λ₂ eigengap ───────────────
    // Self-tuning affinity on 136-dim TNF vectors (Zelnik-Manor σᵢ = dist to 7-NN).
    // v₂ bimodality = normalized max-gap in sorted v₂ entries (high → two clusters).
    if (N >= 4) {
        std::vector<uint32_t> valid;
        valid.reserve(N);
        for (uint32_t i = 0; i < N; ++i)
            if (csk[i].has_tnf) valid.push_back(i);
        if (valid.size() > k_max_fiedler) {
            std::partial_sort(valid.begin(), valid.begin() + k_max_fiedler, valid.end(),
                              [&](uint32_t a, uint32_t b){ return csk[a].len > csk[b].len; });
            valid.resize(k_max_fiedler);
        }
        const uint32_t M = static_cast<uint32_t>(valid.size());
        if (M >= 4) {
            Eigen::MatrixXf D2(M, M);
            for (uint32_t a = 0; a < M; ++a) {
                D2(a, a) = 0.f;
                for (uint32_t b = a + 1; b < M; ++b) {
                    float d = (csk[valid[a]].tnf - csk[valid[b]].tnf).squaredNorm();
                    D2(a, b) = D2(b, a) = d;
                }
            }
            const uint32_t knn = std::min(7u, M - 1);
            Eigen::VectorXf sigma(M);
            for (uint32_t i = 0; i < M; ++i) {
                std::vector<float> row(D2.row(i).data(), D2.row(i).data() + M);
                std::nth_element(row.begin(), row.begin() + knn, row.end());
                sigma(i) = std::sqrt(std::max(0.f, row[knn]));
                if (sigma(i) < 1e-9f) sigma(i) = 1e-9f;
            }
            Eigen::MatrixXf Wt(M, M);
            for (uint32_t a = 0; a < M; ++a) {
                Wt(a, a) = 0.f;
                for (uint32_t b = a + 1; b < M; ++b) {
                    float w = std::exp(-D2(a, b) / (sigma(a) * sigma(b)));
                    Wt(a, b) = Wt(b, a) = w;
                }
            }
            Eigen::VectorXf degt = Wt.rowwise().sum();
            Eigen::MatrixXf Lt  = Eigen::MatrixXf::Identity(M, M);
            for (uint32_t a = 0; a < M; ++a) {
                if (degt(a) < 1e-9f) continue;
                for (uint32_t b = 0; b < M; ++b) {
                    if (a == b) continue;
                    const float dab = degt(a) * degt(b);
                    if (dab > 0.f) Lt(a, b) = -Wt(a, b) / std::sqrt(dab);
                }
            }
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> est(Lt);
            if (est.info() == Eigen::Success) {
                r.fiedler_tnf_gap = est.eigenvalues()(2) - est.eigenvalues()(1);
                Eigen::VectorXf v2 = est.eigenvectors().col(1);
                std::vector<float> v2s(v2.data(), v2.data() + M);
                std::sort(v2s.begin(), v2s.end());
                const float rng = v2s.back() - v2s.front();
                float max_gap = 0.f;
                for (uint32_t i = 0; i + 1 < M; ++i)
                    max_gap = std::max(max_gap, v2s[i + 1] - v2s[i]);
                r.fiedler_tnf_bimod = (rng > 1e-9f) ? max_gap / rng : 0.f;
            }
        }
    }

    return r;
}

FmhGenusRef build_fmh_genus_ref(uint64_t genus_hash,
                                  const std::vector<std::string_view>& genome_seqs,
                                  int k, int c)
{
    FmhGenusRef ref;
    ref.genus_hash = genus_hash;
    for (const auto& seq : genome_seqs) {
        auto sk = build_sketch("", seq, k, c);
        ref.hashes.insert(ref.hashes.end(), sk.hashes.begin(), sk.hashes.end());
    }
    std::sort(ref.hashes.begin(), ref.hashes.end());
    ref.hashes.erase(std::unique(ref.hashes.begin(), ref.hashes.end()), ref.hashes.end());
    return ref;
}

float score_bin_fmh_containment(
    const std::vector<std::string_view>& bin_seqs,
    uint64_t primary_genus_hash,
    const std::vector<FmhGenusRef>& refs,
    int k, int c,
    uint32_t min_contig_bp)
{
    if (refs.size() < 2) return NAN;
    uint64_t minority_bp = 0, scored_bp = 0;
    for (const auto& seq : bin_seqs) {
        if (seq.size() < min_contig_bp) continue;
        const uint32_t clen = static_cast<uint32_t>(seq.size());
        scored_bp += clen;
        auto sk = build_sketch("", seq, k, c);
        if (sk.hashes.empty()) continue;
        uint64_t best_hash = 0;
        double   best_cont = -1.0;
        for (const auto& ref : refs) {
            double cont = fmh_containment(sk.hashes, ref.hashes);
            if (cont > best_cont) { best_cont = cont; best_hash = ref.genus_hash; }
        }
        if (best_hash != 0 && best_hash != primary_genus_hash)
            minority_bp += clen;
    }
    if (scored_bp == 0) return NAN;
    return static_cast<float>(minority_bp) / static_cast<float>(scored_bp);
}

float score_bin_fmh_containment(
    const std::vector<std::string_view>& bin_seqs,
    uint64_t primary_genus_hash,
    const std::vector<FmhrView>& refs,
    int k, int c,
    uint32_t min_contig_bp)
{
    if (refs.size() < 2) return NAN;
    uint64_t minority_bp = 0, scored_bp = 0;
    for (const auto& seq : bin_seqs) {
        if (seq.size() < min_contig_bp) continue;
        const uint32_t clen = static_cast<uint32_t>(seq.size());
        scored_bp += clen;
        auto sk = build_sketch("", seq, k, c);
        if (sk.hashes.empty()) continue;
        uint64_t best_hash = 0;
        double   best_cont = -1.0;
        for (const auto& ref : refs) {
            if (!ref.valid()) continue;
            double cont = fmh_containment(sk.hashes, ref.hashes, ref.n_hashes);
            if (cont > best_cont) { best_cont = cont; best_hash = ref.genus_hash; }
        }
        if (best_hash != 0 && best_hash != primary_genus_hash)
            minority_bp += clen;
    }
    if (scored_bp == 0) return NAN;
    return static_cast<float>(minority_bp) / static_cast<float>(scored_bp);
}

// ── Multi-(k,c) sweep implementations ────────────────────────────────────────

static constexpr int SWEEP_K[2] = {21, 31};
static constexpr int SWEEP_C[3] = {30, 125, 500};

std::array<FmhGenusRef, 6> build_fmh_genus_refs_sweep(
    uint64_t genus_hash,
    const std::vector<std::string_view>& genome_seqs)
{
    std::array<FmhGenusRef, 6> refs;
    for (auto& r : refs) r.genus_hash = genus_hash;

    for (int ki = 0; ki < 2; ++ki) {
        std::array<std::vector<uint64_t>, 3> all;
        for (const auto& seq : genome_seqs) {
            auto vecs = fmh_multi_c(seq, SWEEP_K[ki], SWEEP_C, 3);
            for (int ci = 0; ci < 3; ++ci)
                all[ci].insert(all[ci].end(), vecs[ci].begin(), vecs[ci].end());
        }
        for (int ci = 0; ci < 3; ++ci) {
            auto& h = all[ci];
            std::sort(h.begin(), h.end());
            h.erase(std::unique(h.begin(), h.end()), h.end());
            refs[ki * 3 + ci].hashes = std::move(h);
        }
    }
    return refs;
}

std::array<float, 6> score_bin_fmh_sweep(
    const std::vector<std::string_view>& bin_seqs,
    const std::array<FmhGenusRef, 6>& host_refs,
    const std::array<FmhGenusRef, 6>& contam_refs,
    uint32_t min_contig_bp)
{
    uint64_t minority_bp[6] = {}, scored_bp[6] = {};

    for (const auto& seq : bin_seqs) {
        if (seq.size() < min_contig_bp) continue;
        const uint64_t clen = seq.size();

        for (int ki = 0; ki < 2; ++ki) {
            auto vecs = fmh_multi_c(seq, SWEEP_K[ki], SWEEP_C, 3);
            for (int ci = 0; ci < 3; ++ci) {
                const int idx = ki * 3 + ci;
                const auto& hv   = host_refs[idx].hashes;
                const auto& cv   = contam_refs[idx].hashes;
                if (hv.empty() || cv.empty()) continue;
                scored_bp[idx] += clen;
                const auto& qv = vecs[ci];
                if (qv.empty()) continue;
                if (fmh_containment(qv, cv) > fmh_containment(qv, hv))
                    minority_bp[idx] += clen;
            }
        }
    }

    std::array<float, 6> result;
    for (int i = 0; i < 6; ++i)
        result[i] = scored_bp[i] > 0
                    ? static_cast<float>(minority_bp[i]) / static_cast<float>(scored_bp[i])
                    : NAN;
    return result;
}

BinScores score_bin(
    const std::vector<std::string_view>& seqs,
    const GcovEntry*      gcov_entry,
    const GcovReader*     gcov_rd,
    const GcovEntry*      fcov_entry,
    const GcovReader*     fcov_rd,
    const BinScoreConfig& cfg)
{
    if (!gcov_entry || !gcov_rd)
        return BinScores{BinScores::Status::NoGcovEntry};

    uint32_t gcov_out_bp = 0, spe_out_bp = 0, sibling_out_bp = 0, rho_out_bp = 0, scored_bp = 0;

    for (const auto& seq : seqs) {
        if (seq.size() < cfg.min_long_bp) continue;

        Eigen::VectorXf tnf;
        if (!compute_tnf(seq, tnf)) continue;

        float xmu[136];
        for (int d = 0; d < 136; ++d) xmu[d] = tnf[d] - gcov_entry->mu[d];

        float dist, spe;
        gcov_mahalanobis_spe(*gcov_entry, xmu, &dist, &spe);
        const float pct     = gcov_rd->percentile(*gcov_entry, dist);
        const float spe_pct = gcov_rd->spe_percentile(*gcov_entry, spe);
        const bool  t2_out  = (pct     >= cfg.outlier_percentile);
        const bool  spe_flg = (spe_pct >= cfg.outlier_percentile);

        const uint32_t clen = static_cast<uint32_t>(seq.size());
        if (t2_out || spe_flg) gcov_out_bp += clen;
        if (spe_flg)           spe_out_bp  += clen;
        scored_bp += clen;

        // ρ* Mahalanobis outlier (independent of TNF covariance)
        float rho[16];
        if (compute_rho(seq, rho)) {
            float rho_z[16];
            for (int d = 0; d < 16; ++d) rho_z[d] = rho[d] - gcov_entry->rho_mean[d];
            const float rho_dist = gcov_rho_distance(*gcov_entry, rho_z);
            const float rho_pct  = gcov_rd->rho_percentile(*gcov_entry, rho_dist);
            if (rho_pct >= cfg.outlier_percentile) rho_out_bp += clen;
        }

        // Sibling: genus-outlier but family-inlier
        if ((t2_out || spe_flg) && fcov_entry && fcov_rd) {
            float xmu_fam[136];
            for (int d = 0; d < 136; ++d) xmu_fam[d] = tnf[d] - fcov_entry->mu[d];
            float f_dist, f_spe;
            gcov_mahalanobis_spe(*fcov_entry, xmu_fam, &f_dist, &f_spe);
            const float f_pct  = fcov_rd->percentile(*fcov_entry, f_dist);
            const float f_spct = fcov_rd->spe_percentile(*fcov_entry, f_spe);
            const bool  fam_out = (f_pct  >= cfg.outlier_percentile) ||
                                  (f_spct >= cfg.outlier_percentile);
            if (!fam_out) sibling_out_bp += clen;
        }
    }

    if (scored_bp == 0)
        return BinScores{BinScores::Status::TooFewContigs};

    BinScores r;
    r.status   = BinScores::Status::OK;
    r.scored_bp = scored_bp;
    r.cco_fraction     = static_cast<float>(gcov_out_bp) / static_cast<float>(scored_bp);
    r.sibling_fraction = (fcov_entry && fcov_rd)
        ? static_cast<float>(sibling_out_bp) / static_cast<float>(scored_bp)
        : NAN;
    r.rho_fraction     = static_cast<float>(rho_out_bp) / static_cast<float>(scored_bp);
    return r;
}

} // namespace genopack
