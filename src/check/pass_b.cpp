#include "pass_b.hpp"
#include "tnf.hpp"
#include <genopack/markers.hpp>
#include <genopack/metamer.hpp>
#include <genopack/score_bin.hpp>
#include <genopack/gcov.hpp>
#include <genopack/qual.hpp>
#include <genopack/util.hpp>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <mutex>
#include <omp.h>
#include <spdlog/spdlog.h>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace genopack::check {

// Defined here; declared in pass_b.hpp for use by pass_marker.cpp.
std::vector<ContigRecord> parse_fasta(std::string_view fasta) {
    std::vector<ContigRecord> result;
    const char* p   = fasta.data();
    const char* end = fasta.data() + fasta.size();

    while (p < end) {
        while (p < end && *p != '>') ++p;
        if (p >= end) break;

        const char* hstart = p + 1;
        while (p < end && *p != '\n') ++p;
        std::string header(hstart, p - hstart);
        if (!header.empty() && header.back() == '\r') header.pop_back();
        if (p < end) ++p;

        std::string seq;
        while (p < end && *p != '>') {
            const char* lstart = p;
            while (p < end && *p != '\n' && *p != '\r') ++p;
            seq.append(lstart, p - lstart);
            while (p < end && (*p == '\n' || *p == '\r')) ++p;
        }
        result.push_back({std::move(header), std::move(seq)});
    }
    return result;
}

// Canonical index and 2-bit encoding tables — built once at startup.
struct TnfTables {
    std::array<uint8_t, 256> base2bit;  // ASCII → 2-bit (0-3), 0xFF for non-ACGT
    std::array<uint8_t, 256> canon;     // 8-bit 4-mer index → canonical index (0-135)

    TnfTables() {
        base2bit.fill(0xFF);
        for (char c : {'A','a'}) base2bit[static_cast<uint8_t>(c)] = 0;
        for (char c : {'C','c'}) base2bit[static_cast<uint8_t>(c)] = 1;
        for (char c : {'G','g'}) base2bit[static_cast<uint8_t>(c)] = 2;
        for (char c : {'T','t'}) base2bit[static_cast<uint8_t>(c)] = 3;

        std::array<int, 256> canonical{};
        canonical.fill(-1);
        int next = 0;
        for (int i = 0; i < 256; ++i) {
            if (canonical[i] >= 0) continue;
            int b0 = (i >> 6) & 3, b1 = (i >> 4) & 3,
                b2 = (i >> 2) & 3, b3 = i & 3;
            int rc = (3-b3)*64 + (3-b2)*16 + (3-b1)*4 + (3-b0);
            canonical[i] = canonical[rc] = next++;
            canon[i] = canon[rc] = static_cast<uint8_t>(canonical[i]);
        }
    }
};
static const TnfTables kTnfTables;

// Stack-array overload: writes 136 normalized floats into caller-provided storage.
bool compute_tnf(std::string_view seq, float out[136]) {
    if (seq.size() < 5000) return false;
    std::fill(out, out + 136, 0.0f);

    const uint8_t* b2b = kTnfTables.base2bit.data();
    uint32_t counts[256] = {};

    uint32_t idx = 0, valid = 0;
    for (size_t i = 0; i < seq.size(); ++i) {
        uint8_t bits = b2b[static_cast<uint8_t>(seq[i])];
        if (bits == 0xFF) { valid = 0; idx = 0; continue; }
        idx = ((idx << 2) | bits) & 0xFF;
        if (++valid >= 4)
            ++counts[idx];
    }

    const uint8_t* canon = kTnfTables.canon.data();
    for (int i = 0; i < 256; ++i)
        out[canon[i]] += static_cast<float>(counts[i]);

    float norm_sq = 0.0f;
    for (int i = 0; i < 136; ++i) norm_sq += out[i] * out[i];
    const float norm = std::sqrt(norm_sq);
    if (norm < 1e-8f) return false;
    for (int i = 0; i < 136; ++i) out[i] /= norm;
    return true;
}

void run_pass_b(ICheckReader& pack,
                const PassAResult& pass_a,
                std::unordered_map<std::string, GenomeQuality>& quality,
                const PassBConfig& cfg,
                int threads,
                const std::unordered_map<uint64_t, QualRecord>* qual_cache,
                const std::unordered_map<uint64_t, QualRecord>* baseline_cache)
{
    // needs_full_analysis: completeness or leakage triggered — run mixture+Fiedler window model.
    // tnf-only flagged genomes skip the expensive PCA/GMM pass (skip_mixture=true).
    std::vector<std::string> flagged;
    std::unordered_set<std::string> needs_full_set;
    flagged.reserve(pass_a.accessions.size() / 10);
    for (const auto& acc : pass_a.accessions) {
        auto it = quality.find(acc);
        if (it == quality.end()) continue;
        const auto& q = it->second;
        const bool by_leakage     = q.contamination_leakage > cfg.contamination_flag_threshold;
        const bool by_tnf         = q.contamination_tnf_excess > cfg.tnf_flag_threshold;
        const bool by_completeness = !std::isnan(q.completeness_cluster_relative) &&
                                     q.completeness_cluster_relative < cfg.completeness_flag_threshold;
        if (by_leakage || by_tnf || by_completeness) {
            flagged.push_back(acc);
            if (by_leakage || by_completeness)
                needs_full_set.insert(acc);
        }
    }

    if (flagged.empty()) return;
    spdlog::info("check pass-B: {} genomes flagged for FASTA-level analysis", flagged.size());

    const FmhrReader* fmhr_rd = pack.fmhr_reader();

    // Serve from QUAL cache where both skew and post-decontam scores are precomputed.
    std::vector<std::string> to_scan;
    to_scan.reserve(flagged.size());
    if (qual_cache) {
        for (const auto& acc : flagged) {
            auto meta = pack.genome_meta_by_accession(acc);
            if (!meta) { to_scan.push_back(acc); continue; }
            auto it = qual_cache->find(meta->genome_id);
            if (it == qual_cache->end()) { to_scan.push_back(acc); continue; }
            const QualRecord& r = it->second;
            const bool has_build_signals = (r.qual_flags & QualRecord::QUAL_FLAG_BUILD_SIGNALS) != 0;
            const bool has_gcov_scored   = (r.qual_flags & QualRecord::QUAL_FLAG_GCOV_SCORED)   != 0;
            if (std::isnan(r.chromosome_skew_closure) || std::isnan(r.completeness_post_decontam) ||
                (!has_build_signals && r.self_coherence == 0.0f) ||
                (fmhr_rd && r.fmh_minority_u8 == 0 && !std::isnan(r.chromosome_skew_closure)) ||
                (!cfg.markers_path.empty() && r.marker_completeness_u8 == 0) ||
                !has_gcov_scored) {
                to_scan.push_back(acc); continue;
            }
            auto& q = quality.at(acc);
            q.chromosome_skew_closure    = r.chromosome_skew_closure;
            q.completeness_post_decontam = r.completeness_post_decontam;
            q.self_coherence             = r.self_coherence;
            q.chargaff_parity            = r.chargaff_parity;
            q.spectral_gap               = r.spectral_gap;
            q.scale_kink                 = r.scale_kink;
            // No FMHR section -> the FMH axis is unavailable; degrade to NaN so
            // cache-served genomes match the FASTA-scan path (which leaves it NaN)
            // and the axis-fallback flag fires consistently.
            q.fmh_minority_fraction      = fmhr_rd ? (r.fmh_minority_u8 / 255.0f) : NAN;
            if (r.marker_completeness_u8 > 0)
                q.marker_completeness = (r.marker_completeness_u8 - 1) / 254.0f;
            q.contamination_contig_outlier  = r.contig_outlier_u8  / 255.0f;
            q.contamination_spe             = r.spe_outlier_u8     / 255.0f;
            q.contamination_sibling_outlier = r.sibling_outlier_u8 / 255.0f;
            q.contamination_rho_outlier     = r.rho_outlier_u8     / 255.0f;
            q.contamination_cross_genus     = r.cross_genus_u8     / 255.0f;
            if (r.sketch_fill_u8 > 0)
                q.completeness_sketch_fill = r.sketch_fill_u8 / 200.0f;
        }
        if (to_scan.empty()) {
            spdlog::info("check pass-B: complete — all scores from QUAL cache");
            return;
        }
        spdlog::info("check pass-B: {}/{} genomes need FASTA scan ({} from cache)",
                     to_scan.size(), flagged.size(), flagged.size() - to_scan.size());
    } else {
        to_scan = std::move(flagged);
    }

    // Build acc→genus map once; reused for centroid accumulation and per-genome lookup.
    std::unordered_map<std::string, std::string> flagged_genus;
    for (const auto& [genus, members] : pass_a.genus_members)
        for (const auto& acc : members)
            flagged_genus[acc] = genus;

    std::unordered_map<std::string, Eigen::VectorXf> genus_centroid;
    {
        std::unordered_map<std::string, std::pair<Eigen::VectorXf, int>> sums;
        for (const auto& acc : pass_a.accessions) {
            const float* p = pack.kmer_profile_by_accession(acc);
            if (!p) continue;
            auto git = flagged_genus.find(acc);
            if (git == flagged_genus.end()) continue;
            const std::string& genus = git->second;
            auto& [sum, cnt] = sums[genus];
            if (cnt == 0) sum = Eigen::VectorXf::Zero(136);
            sum += Eigen::Map<const Eigen::VectorXf>(p, 136);
            ++cnt;
        }
        for (auto& [g, sc] : sums) {
            if (sc.second > 0) {
                sc.first /= static_cast<float>(sc.second);
                genus_centroid[g] = sc.first;
            }
        }
    }

    std::unordered_map<std::string, std::string> flagged_family;
    for (const auto& [family, members] : pass_a.family_members)
        for (const auto& acc : members)
            flagged_family[acc] = family;

    // Pre-build family → GSTX candidate entries for containment split.
    const GstxReader* gstx_rd = pack.gstx_reader();
    std::unordered_map<std::string, std::vector<const GstxEntry*>> family_gstx_candidates;
    if (gstx_rd) {
        for (const auto& [family, genera] : pass_a.family_to_genera) {
            auto& candidates = family_gstx_candidates[family];
            for (const auto& genus_name : genera) {
                const GstxEntry* ge = pack.gstx_for_genus(genus_name);
                if (ge && ge->genus_hash != 0 && ge->n_k_stored > 0)
                    candidates.push_back(ge);
            }
        }
    }

    // Pre-cache FmhrView for each unique genus (zero-copy pointers into mmap'd FMHR section).
    // Also cache candidate refs per family = all genera in that family with valid FMHR entries.
    std::unordered_map<std::string, FmhrView> fmhr_host_cache;
    std::unordered_map<std::string, std::vector<FmhrView>> fmhr_family_candidates;
    if (fmhr_rd) {
        for (const auto& [acc, genus] : flagged_genus) {
            if (fmhr_host_cache.count(genus)) continue;
            auto v = pack.fmhr_for_genus(genus);
            if (v.valid()) fmhr_host_cache[genus] = v;
        }
        for (const auto& [family, genera] : pass_a.family_to_genera) {
            auto& cands = fmhr_family_candidates[family];
            for (const auto& g : genera) {
                auto v = pack.fmhr_for_genus(g);
                if (v.valid()) cands.push_back(v);
            }
        }
        spdlog::info("check pass-B: FMHR cache: {} host genera, {} families with candidates",
                     fmhr_host_cache.size(), fmhr_family_candidates.size());
    }

    // Open marker panel (read-only, mmap — safe for concurrent reads across all threads).
    std::unique_ptr<MarkerReader> mrk_rd;
    if (!cfg.markers_path.empty()) {
        mrk_rd = std::make_unique<MarkerReader>();
        try {
            mrk_rd->open(cfg.markers_path);
            spdlog::info("check pass-B: marker panel: {} bac + {} arc, {} lineages",
                         mrk_rd->n_bac(), mrk_rd->n_arc(), mrk_rd->n_lineages());
        } catch (const std::exception& ex) {
            spdlog::warn("check pass-B: cannot open markers ({}); skipping marker scoring", ex.what());
            mrk_rd.reset();
        }
    }
    if (mrk_rd && mrk_rd->is_open()) {
        mrk_rd->build_merged_pool();
        spdlog::info("check pass-B: merged pool: {} M bac + {} M arc hashes",
                     mrk_rd->merged_hashes_bac().size() / 1'000'000,
                     mrk_rd->merged_hashes_arc().size() / 1'000'000);
    }

    const uint64_t mrk_frac_max = mrk_rd ? mrk_rd->frac_max_hash() : UINT64_MAX;

    // N parallel reader threads each own an independent fd and sequential shard band.
    // Callback is invoked concurrently; quality_mtx guards shared state.
    constexpr int k_shard_readers = 8;
    const int inner_threads = std::max(1, threads / k_shard_readers);
    std::mutex quality_mtx;
    size_t n_done = 0;
    const size_t n_flagged = to_scan.size();

    pack.visit_shard_batches_parallel(to_scan, k_shard_readers,
        [&](ArchiveReader::ShardBatch& batch) {
            const int n = static_cast<int>(batch.size());
            #pragma omp parallel for schedule(dynamic, 1) num_threads(inner_threads)
            for (int j = 0; j < n; ++j) {
                const auto& [idx, eg] = batch[static_cast<size_t>(j)];
                if (eg.fasta.empty()) continue;
                const std::string& acc = to_scan[idx];

                auto contigs = parse_fasta(eg.fasta);
                if (contigs.empty()) continue;

                // GC skew inline — no full concatenation needed.
                float skew_score = NAN;
                {
                    size_t total_len = 0;
                    for (const auto& c : contigs) total_len += c.seq.size();
                    if (total_len >= 10000) {
                        int64_t cum = 0, max_abs = 0;
                        for (const auto& c : contigs) {
                            for (char ch : c.seq) {
                                switch (ch | 0x20) {
                                case 'g': ++cum; break;
                                case 'c': --cum; break;
                                default:  break;
                                }
                                if (std::abs(cum) > max_abs) max_abs = std::abs(cum);
                            }
                        }
                        if (max_abs > 0)
                            skew_score = 1.0f - static_cast<float>(std::abs(cum)) /
                                                static_cast<float>(max_abs);
                    }
                }

                auto git = flagged_genus.find(acc);
                const Eigen::VectorXf* centroid = nullptr;
                const GcovEntry* gcov_entry = nullptr;
                const GcovReader* gcov_rd   = pack.gcov_reader();
                if (git != flagged_genus.end()) {
                    auto cit = genus_centroid.find(git->second);
                    if (cit != genus_centroid.end()) centroid = &cit->second;
                    const GcovEntry* ge = pack.gcov_for_genus(git->second);
                    if (ge && (ge->flags & GCOV_FLAG_VALID)) gcov_entry = ge;
                }

                const GcovReader* fcov_rd    = pack.fcov_reader();
                const GcovEntry*  fcov_entry = nullptr;
                auto fit = flagged_family.find(acc);
                if (fit != flagged_family.end() && fcov_rd) {
                    const GcovEntry* fe = pack.fcov_for_family(fit->second);
                    if (fe && (fe->flags & GCOV_FLAG_VALID)) fcov_entry = fe;
                }

                std::vector<ContigFlag> flags;
                // ctg_cross_genus[i] = true if contigs[i] was flagged by cross_genus log-LR test.
                // Used below for the joint marker+GCOV contamination signal.
                std::vector<bool> ctg_cross_genus(contigs.size(), false);
                uint32_t byte_offset    = 0;
                uint32_t clean_len      = 0;
                uint32_t total_len      = 0;
                uint32_t gcov_out_bp    = 0;
                uint32_t gcov_scored_bp = 0;
                uint32_t spe_out_bp     = 0;
                uint32_t sibling_out_bp = 0;
                uint32_t rho_out_bp     = 0;
                uint32_t cross_genus_bp = 0;

                // Per-genome TNF cache: keeps computed TNF alive for compute_all_signals overload.
                // Thread-local to avoid per-genome allocation.
                thread_local std::vector<std::array<float,136>> tnf_cache;
                tnf_cache.clear();
                tnf_cache.reserve(contigs.size());

                // Qualifying contigs for the progressive-bound cross-genus scan.
                // τ = host_ll + margin is the threshold a foreign genus must beat.
                // Clean contigs (host fits well) prune ~99% of genera at the L2/λ_max step.
                // Contaminated contigs early-exit the scan on the first winner.
                struct QualContig {
                    const float* tnf;
                    float tau;          // host_ll + cross_genus_lr_margin
                    bool  family_inlier;
                    uint32_t length;
                    size_t   ci;
                    bool     found;     // true once flagged — skip in subsequent genera
                };
                thread_local std::vector<QualContig> qual_contigs;
                qual_contigs.clear();
                // Build ContigAccum for fused compute_all_signals (eliminates raw FASTA re-scan).
                thread_local std::vector<ContigAccum> accum_contigs;
                accum_contigs.clear();
                accum_contigs.reserve(contigs.size());

                // When no genus centroid is available, all tnf_scores are NaN →
                // every contig counts as clean → completeness_post_decontam = 1.0.
                // Skip the per-contig TNF computation entirely in that case.
                if (!centroid) {
                    for (const auto& contig : contigs) {
                        uint32_t clen = static_cast<uint32_t>(contig.seq.size());
                        flags.push_back({byte_offset, clen, NAN, NAN});
                        byte_offset += clen;
                        clean_len   += clen;
                        total_len   += clen;
                        accum_contigs.push_back({contig.seq, nullptr, clen});
                    }
                } else {
                    // Hoist centroid norm — constant across all contigs for this genome.
                    const float norm_c = centroid->norm();
                    flags.reserve(contigs.size());
                    size_t contig_ci = 0;
                    for (const auto& contig : contigs) {
                        const size_t ci = contig_ci++;
                        ContigFlag cf;
                        cf.contig_offset = byte_offset;
                        cf.contig_length = static_cast<uint32_t>(contig.seq.size());
                        byte_offset += cf.contig_length;
                        total_len   += cf.contig_length;

                        if (contig.seq.size() >= 5000) {
                            // Stack-allocated TNF stored in per-genome cache for compute_all_signals.
                            tnf_cache.push_back({});
                            float* tnf = tnf_cache.back().data();
                            if (compute_tnf(contig.seq, tnf)) {
                                accum_contigs.push_back({contig.seq, tnf, cf.contig_length});
                                // Inline d.norm() — avoids Eigen heap vector for d = tnf - centroid.
                                float d_norm_sq = 0.0f;
                                for (int di = 0; di < 136; ++di) {
                                    float diff = tnf[di] - (*centroid)[di];
                                    d_norm_sq += diff * diff;
                                }
                                cf.tnf_score = (norm_c > 1e-8f)
                                    ? std::sqrt(d_norm_sq) / norm_c : NAN;

                                if (gcov_entry && gcov_rd && cf.contig_length >= gcov_rd->min_long_bp()) {
                                    float xmu[136];
                                    for (int di = 0; di < 136; ++di) xmu[di] = tnf[di] - gcov_entry->mu[di];
                                    const float dist    = gcov_mahalanobis(*gcov_entry, xmu);
                                    const float pct     = gcov_rd->percentile(*gcov_entry, dist);
                                    const float spe     = gcov_spe(*gcov_entry, xmu);
                                    const float spe_pct = gcov_rd->spe_percentile(*gcov_entry, spe);
                                    const bool t2_out   = (pct     >= cfg.gcov_outlier_percentile);
                                    const bool spe_flag = (spe_pct >= cfg.gcov_outlier_percentile);
                                    if (t2_out || spe_flag) gcov_out_bp += cf.contig_length;
                                    if (spe_flag)           spe_out_bp  += cf.contig_length;
                                    gcov_scored_bp += cf.contig_length;

                                    // ρ* outlier (independent of TNF covariance model)
                                    float rho_q[16];
                                    if (compute_rho(contig.seq, rho_q)) {
                                        float rho_z[16];
                                        for (int di2=0;di2<16;++di2) rho_z[di2]=rho_q[di2]-gcov_entry->rho_mean[di2];
                                        const float rd  = gcov_rho_distance(*gcov_entry, rho_z);
                                        const float rpc = gcov_rd->rho_percentile(*gcov_entry, rd);
                                        if (rpc >= cfg.gcov_outlier_percentile) rho_out_bp += cf.contig_length;
                                    }

                                    // Collect for progressive-bound cross-genus scan.
                                    {
                                        const float host_ll = gcov_log_likelihood(
                                            *gcov_entry, dist * dist, spe);
                                        const float tau = host_ll + cfg.cross_genus_lr_margin;
                                        bool family_inlier = false;
                                        if (fcov_entry && fcov_rd) {
                                            float xmu_fam[136];
                                            for (int di = 0; di < 136; ++di)
                                                xmu_fam[di] = tnf[di] - fcov_entry->mu[di];
                                            const float f_dist = gcov_mahalanobis(*fcov_entry, xmu_fam);
                                            const float f_pct  = fcov_rd->percentile(*fcov_entry, f_dist);
                                            const float f_spe  = gcov_spe(*fcov_entry, xmu_fam);
                                            const float f_spct = fcov_rd->spe_percentile(*fcov_entry, f_spe);
                                            family_inlier = (f_pct  < cfg.gcov_outlier_percentile) &&
                                                            (f_spct < cfg.gcov_outlier_percentile);
                                        }
                                        qual_contigs.push_back({tnf, tau,
                                                                family_inlier, cf.contig_length, ci, false});
                                    }
                                }
                            } else {
                                // TNF compute failed (degenerate contig) — still include in
                                // accum_contigs so transitions/chargaff walk sees the sequence.
                                tnf_cache.pop_back();
                                accum_contigs.push_back({contig.seq, nullptr, cf.contig_length});
                            }
                        }

                        bool contig_clean = !std::isnan(cf.tnf_score) &&
                                            cf.tnf_score < cfg.contig_tnf_threshold;
                        if (contig_clean || std::isnan(cf.tnf_score))
                            clean_len += cf.contig_length;

                        flags.push_back(cf);
                    }
                }

                // Progressive-bound cross-genus scan.
                //
                // For each (genus, contig) pair we maintain a running admissible lower
                // bound on (d² + spe/σ²) after j eigenvector projections:
                //   B_j = Σ_{k≤j} proj_k²/λ_k  +  (‖xf‖² − Σ_{k≤j} proj_k²) / max(λ_{j+1}, σ²)
                //
                // B_j is non-decreasing in j; at j=−1 it collapses to ‖xf‖²/λ_max (cheap L2 guard).
                // If −0.5·(B_j + log_det) ≤ τ at any point, the genus can't beat the threshold
                // and we abandon. At j=K−1 the bound equals the exact LL.
                //
                // τ = host_ll + margin is pre-computed per contig. Clean contigs (host fits well)
                // prune ~99% of genera at j=−1. Contaminated contigs exit the scan on first hit.
                if (!qual_contigs.empty() && gcov_rd) {
                    uint32_t n_unfound = static_cast<uint32_t>(qual_contigs.size());
                    constexpr int K = static_cast<int>(GCOV_N_EIGVECS);

                    gcov_rd->scan_early([&](const GcovEntry& fe) -> bool {
                        if (&fe == gcov_entry || !(fe.flags & GCOV_FLAG_VALID)) return true;

                        const float ld  = gcov_log_det(fe);
                        const float vmax = std::max(fe.eigenvalues[0], fe.sigma2_resid);

                        for (auto& qc : qual_contigs) {
                            if (qc.found) continue;

                            // j=−1: L2/λ_max admissible upper bound (136 multiply-adds)
                            float xf[136], l2 = 0.f;
                            for (int d = 0; d < 136; ++d) {
                                const float v = qc.tnf[d] - fe.mu[d];
                                xf[d] = v; l2 += v * v;
                            }
                            if (-0.5f * (l2 / vmax + ld) <= qc.tau) continue;

                            // j=0..K−1: tighten bound one eigenvector at a time
                            float sum_p2 = 0.f, d2 = 0.f;
                            bool pruned = false;
                            for (int k = 0; k < K; ++k) {
                                float dot = 0.f;
                                for (int d = 0; d < 136; ++d) dot += xf[d] * fe.eigvecs[k][d];
                                sum_p2 += dot * dot;
                                d2     += dot * dot / fe.eigenvalues[k];
                                const float vsub = (k + 1 < K)
                                    ? std::max(fe.eigenvalues[k + 1], fe.sigma2_resid)
                                    : std::max(fe.sigma2_resid, 1e-12f);
                                const float spe_lb = std::max(0.f, l2 - sum_p2);
                                if (-0.5f * (d2 + spe_lb / vsub + ld) <= qc.tau) {
                                    pruned = true; break;
                                }
                            }
                            if (pruned) continue;

                            // Survived all K projections → LL > τ (bound is exact at j=K−1)
                            qc.found = true;
                            cross_genus_bp      += qc.length;
                            ctg_cross_genus[qc.ci] = true;
                            if (qc.family_inlier) sibling_out_bp += qc.length;
                            if (--n_unfound == 0) return false; // stop scan
                        }
                        return true; // continue scan
                    });
                }

                float completeness_post = total_len > 0
                    ? static_cast<float>(clean_len) / static_cast<float>(total_len)
                    : NAN;

                // Per-contig outlier detection.
                // Primary: GCOV archive-calibrated percentile (no self-reference).
                // Fallback: within-genome median/MAD z-score (when no GCOV entry).
                float contamination_contig_outlier  = NAN;
                float contamination_spe             = NAN;
                float contamination_sibling_outlier = NAN;
                float contamination_rho_outlier     = NAN;
                float contamination_cross_genus     = NAN;
                if (gcov_entry && gcov_rd && gcov_scored_bp > 0) {
                    contamination_contig_outlier =
                        static_cast<float>(gcov_out_bp) / static_cast<float>(gcov_scored_bp);
                    contamination_spe =
                        static_cast<float>(spe_out_bp) / static_cast<float>(gcov_scored_bp);
                    contamination_rho_outlier =
                        static_cast<float>(rho_out_bp) / static_cast<float>(gcov_scored_bp);
                    contamination_cross_genus =
                        static_cast<float>(cross_genus_bp) / static_cast<float>(gcov_scored_bp);
                    if (fcov_entry)
                        contamination_sibling_outlier =
                            static_cast<float>(sibling_out_bp) / static_cast<float>(gcov_scored_bp);
                } else if (centroid) {
                    std::vector<std::pair<float,uint32_t>> scored;
                    for (const auto& cf : flags) {
                        if (!std::isnan(cf.tnf_score) && cf.contig_length >= 5000)
                            scored.emplace_back(cf.tnf_score, cf.contig_length);
                    }
                    if (scored.size() >= 3) {
                        std::vector<float> sv;
                        sv.reserve(scored.size());
                        for (const auto& [s, l] : scored) sv.push_back(s);
                        std::sort(sv.begin(), sv.end());
                        const float null_med = sv[sv.size() / 2];
                        std::vector<float> dv;
                        dv.reserve(sv.size());
                        for (float s : sv) dv.push_back(std::abs(s - null_med));
                        std::sort(dv.begin(), dv.end());
                        float null_mad = 1.4826f * dv[dv.size() / 2];
                        if (null_mad < 1e-6f) null_mad = 1e-6f;
                        uint32_t outlier_bp = 0, scored_bp = 0;
                        for (const auto& [score, clen] : scored) {
                            const float n_win = std::max(1.0f, static_cast<float>(clen) / 16384.0f);
                            const float z = (score - null_med) / (null_mad / std::sqrt(n_win));
                            if (z > 3.0f) outlier_bp += clen;
                            scored_bp += clen;
                        }
                        if (scored_bp > 0)
                            contamination_contig_outlier =
                                static_cast<float>(outlier_bp) / static_cast<float>(scored_bp);
                    }
                }

                // ── Marker-panel completeness + redundancy ────────────────────────────────
                // Per-contig voting: each contig independently scores all 6-frame ORFs
                // against the marker pool. A contig "votes" for marker M if any ORF
                // clears the containment threshold. Redundancy = markers with ≥2 votes
                // (two separate genomic loci → contamination signal).
                // Works for both full-AA k=8 (legacy) and Dayhoff-6 k=12 syncmer pools.
                float marker_completeness        = NAN;
                float marker_redundancy          = NAN;
                float marker_joint_contamination = NAN;
                int   marker_n_present = -1, marker_n_expected = -1;
                if (mrk_rd && mrk_rd->has_merged_pool() && git != flagged_genus.end()) {
                    const uint64_t genus_hash = GcovWriter::hash_genus(git->second);
                    auto calib = mrk_rd->lookup_lineage(genus_hash);
                    if (calib.valid()) {
                        const bool is_arc      = (calib.header->domain == MRKR_DOMAIN_ARC);
                        const bool is_d6       = mrk_rd->is_dayhoff6();
                        auto mh   = is_arc ? mrk_rd->merged_hashes_arc() : mrk_rd->merged_hashes_bac();
                        auto mid  = is_arc ? mrk_rd->merged_ids_arc()    : mrk_rd->merged_ids_bac();
                        const BlockedBloom& bloom =
                            is_arc ? mrk_rd->merged_bloom_arc() : mrk_rd->merged_bloom_bac();
                        const uint8_t  n_markers = calib.header->n_markers;
                        // min_seg=50 for Dayhoff-6: ORFs <50 AA can never accumulate
                        // min_hits=5 syncmer hits (max ~5 syncmers at 11% rate → needs 100% hit
                        // rate at any divergence — impossible). Skip them with zero sensitivity loss.
                        const int min_seg = is_d6 ? 50 : METAMER_K;
                        const uint32_t min_hits  = static_cast<uint32_t>(cfg.marker_min_hits);

                        // contig_votes_native[mi]  = host-origin contigs voting for marker mi.
                        // contig_votes_foreign[mi] = cross_genus-flagged contigs voting for marker mi.
                        // Joint signal: both vote → marker present in two different organisms.
                        uint32_t contig_votes_native[173]  = {};
                        uint32_t contig_votes_foreign[173] = {};

                        thread_local std::vector<uint64_t> orf_mers;
                        for (size_t ci2 = 0; ci2 < contigs.size(); ++ci2) {
                            const auto& contig = contigs[ci2];
                            const bool is_foreign = (ci2 < ctg_cross_genus.size()) && ctg_cross_genus[ci2];
                            uint32_t best[173] = {};

                            if (is_d6) {
                                // Dayhoff-6: per-ORF scoring — threshold normalised to each
                                // individual ORF's syncmer count, not the whole contig.
                                // A 300AA ORF has ~32 syncmers; threshold = max(min_hits, 32/20)=2.
                                translate_6frame(contig.seq, min_seg,
                                    [&](int, const uint8_t* seg, int len, int, int) {
                                        orf_mers.clear();
                                        extract_d6_orf_syncmers(seg, len, orf_mers);
                                        if (orf_mers.empty()) return;
                                        std::sort(orf_mers.begin(), orf_mers.end());
                                        orf_mers.erase(std::unique(orf_mers.begin(), orf_mers.end()),
                                                       orf_mers.end());

                                        uint32_t local[173] = {};
                                        const uint64_t* mhp = mh.data();
                                        const uint64_t* mhe = mhp + mh.size();
                                        const uint8_t*  mip = mid.data();
                                        for (uint64_t h : orf_mers) {
                                            if (!bloom.might_contain(h)) continue;
                                            auto it = std::lower_bound(mhp, mhe, h);
                                            if (it != mhe && *it == h)
                                                local[mip[it - mhp]]++;
                                        }
                                        const uint32_t q_sz = static_cast<uint32_t>(orf_mers.size());
                                        const uint32_t thr  = std::max(min_hits,
                                                                        std::max(1u, q_sz / 20u));
                                        for (uint8_t mi = 0; mi < n_markers; ++mi)
                                            if (local[mi] >= thr && local[mi] > best[mi])
                                                best[mi] = local[mi];
                                    });
                            } else {
                                // Full-AA k=8: whole-contig accumulation (legacy).
                                orf_mers.clear();
                                extract_metamers_dna_into(contig.seq, METAMER_K, min_seg,
                                                         mrk_frac_max, orf_mers);
                                if (!orf_mers.empty()) {
                                    std::sort(orf_mers.begin(), orf_mers.end());
                                    orf_mers.erase(std::unique(orf_mers.begin(), orf_mers.end()),
                                                   orf_mers.end());
                                    uint32_t hits[173] = {};
                                    const uint64_t* mhp = mh.data();
                                    const uint64_t* mhe = mhp + mh.size();
                                    const uint8_t*  mip = mid.data();
                                    for (uint64_t h : orf_mers) {
                                        if (!bloom.might_contain(h)) continue;
                                        auto it = std::lower_bound(mhp, mhe, h);
                                        if (it != mhe && *it == h) hits[mip[it - mhp]]++;
                                    }
                                    const uint32_t q_sz = static_cast<uint32_t>(orf_mers.size());
                                    const uint32_t thr  = std::max(min_hits,
                                                                    std::max(1u, q_sz / 20u));
                                    for (uint8_t mi = 0; mi < n_markers; ++mi)
                                        if (hits[mi] >= thr && hits[mi] > best[mi])
                                            best[mi] = hits[mi];
                                }
                            }

                            for (uint8_t mi = 0; mi < n_markers; ++mi) {
                                if (best[mi] > 0) {
                                    if (is_foreign) contig_votes_foreign[mi]++;
                                    else            contig_votes_native[mi]++;
                                }
                            }
                        }

                        marker_n_present  = 0;
                        marker_n_expected = 0;
                        int redundant_n = 0, joint_n = 0;
                        for (uint8_t mi = 0; mi < n_markers; ++mi) {
                            if (!calib.marker_expected(mi)) continue;
                            ++marker_n_expected;
                            const uint32_t tot = contig_votes_native[mi] + contig_votes_foreign[mi];
                            if (tot >= 1) ++marker_n_present;
                            if (tot >= 2) ++redundant_n;
                            // Joint: confirmed duplication across organisms (native + foreign vote)
                            if (contig_votes_native[mi] >= 1 && contig_votes_foreign[mi] >= 1)
                                ++joint_n;
                        }
                        if (marker_n_expected > 0) {
                            const float ne = static_cast<float>(marker_n_expected);
                            marker_completeness        = static_cast<float>(marker_n_present) / ne;
                            marker_redundancy          = static_cast<float>(redundant_n) / ne;
                            marker_joint_contamination = static_cast<float>(joint_n) / ne;
                        }
                    }
                }

                // Joint: marker_joint_contamination already assigned inside scoring block above.

                // Z-score marker_redundancy against per-genus calibration.
                float marker_redundancy_z = NAN;
                if (!std::isnan(marker_redundancy) && mrk_rd && git != flagged_genus.end()) {
                    const uint64_t genus_hash = GcovWriter::hash_genus(git->second);
                    const auto* ce = mrk_rd->lookup_redun_calib(genus_hash);
                    marker_redundancy_z = MarkerReader::redun_zscore(marker_redundancy, ce);
                }

                // Per-contig containment split against GSTX family-member consensus sketches.
                float contamination_contig_split = NAN;
                float contamination_self_outlier = NAN;
                float fiedler_oph_val            = NAN;
                float fiedler_tnf_bimod_val      = NAN;
                float fiedler_tnf_gap_val        = NAN;
                if (gstx_rd) {
                    // Gather family candidates (may be empty; minority_fraction needs ≥2).
                    const std::vector<const GstxEntry*>* candidates_ptr = nullptr;
                    auto ffit2 = flagged_family.find(acc);
                    if (ffit2 != flagged_family.end()) {
                        auto cit = family_gstx_candidates.find(ffit2->second);
                        if (cit != family_gstx_candidates.end())
                            candidates_ptr = &cit->second;
                    }
                    // Always run for reference-free metrics (self_outlier, Fiedler).
                    // minority_fraction / contig_split are set only when ≥2 candidates
                    // (enforced inside score_bin_containment section 1).
                    {
                        static const std::vector<const GstxEntry*> empty_cands;
                        const auto& cands = candidates_ptr ? *candidates_ptr : empty_cands;
                        std::vector<std::string_view> seqs;
                        seqs.reserve(contigs.size());
                        for (const auto& c : contigs) seqs.push_back(c.seq);
                        uint64_t pgh = 0;
                        if (git != flagged_genus.end())
                            pgh = GcovWriter::hash_genus(git->second);
                        auto csr = score_bin_containment(seqs, pgh, cands, *gstx_rd);
                        contamination_contig_split = csr.minority_fraction;
                        contamination_self_outlier = csr.self_outlier_fraction;
                        fiedler_oph_val            = csr.fiedler_oph;
                        fiedler_tnf_bimod_val      = csr.fiedler_tnf_bimod;
                        fiedler_tnf_gap_val        = csr.fiedler_tnf_gap;
                    }
                }

                // FMH minority fraction: zero-copy mmap'd FMHR refs.
                float fmh_minority = NAN;
                if (fmhr_rd && git != flagged_genus.end()) {
                    auto hit = fmhr_host_cache.find(git->second);
                    if (hit != fmhr_host_cache.end()) {
                        std::vector<FmhrView> refs;
                        refs.push_back(hit->second);
                        auto ffit3 = flagged_family.find(acc);
                        if (ffit3 != flagged_family.end()) {
                            auto cit = fmhr_family_candidates.find(ffit3->second);
                            if (cit != fmhr_family_candidates.end()) {
                                for (const auto& v : cit->second)
                                    if (v.genus_hash != hit->second.genus_hash)
                                        refs.push_back(v);
                            }
                        }
                        if (refs.size() >= 2) {
                            std::vector<std::string_view> seqs;
                            seqs.reserve(contigs.size());
                            for (const auto& c : contigs) seqs.push_back(c.seq);
                            fmh_minority = score_bin_fmh_containment(
                                seqs, hit->second.genus_hash, refs, 21, 125);
                        }
                    }
                }

                {
                    const bool skip_mix = (needs_full_set.find(acc) == needs_full_set.end());
                    AllSignals sig = compute_all_signals(
                        std::span<const ContigAccum>(accum_contigs), 5000, 16384, 5, skip_mix);

                    std::lock_guard<std::mutex> lk(quality_mtx);
                    auto& q = quality.at(acc);
                    q.chromosome_skew_closure    = skew_score;
                    q.completeness_post_decontam = completeness_post;
                    q.self_coherence             = sig.self_coherence;
                    q.chargaff_parity            = sig.chargaff_parity;
                    q.spectral_gap               = sig.spectral_gap;
                    q.scale_kink                 = sig.scale_kink;
                    q.contamination_mixture      = sig.contamination_mixture;
                    q.mixture_sources            = sig.mixture_sources;
                    q.n_mix_windows              = sig.n_mix_windows;
                    q.fiedler_value                   = sig.fiedler_value;
                    q.contamination_contig_outlier    = contamination_contig_outlier;
                    q.contamination_spe               = contamination_spe;
                    q.contamination_sibling_outlier   = contamination_sibling_outlier;
                    q.contamination_rho_outlier       = contamination_rho_outlier;
                    q.contamination_cross_genus       = contamination_cross_genus;
                    q.contamination_contig_split      = contamination_contig_split;
                    q.contamination_self_outlier      = contamination_self_outlier;
                    q.fiedler_oph_split               = fiedler_oph_val;
                    q.fiedler_tnf_bimod               = fiedler_tnf_bimod_val;
                    q.fiedler_tnf_gap                 = fiedler_tnf_gap_val;
                    q.fmh_minority_fraction           = fmh_minority;
                    q.marker_completeness             = marker_completeness;
                    q.marker_redundancy               = marker_redundancy;
                    q.marker_redundancy_z             = marker_redundancy_z;
                    q.marker_joint_contamination      = marker_joint_contamination;
                    q.marker_n_present                = marker_n_present;
                    q.marker_n_expected               = marker_n_expected;

                    if (sig.mix_no_data)
                        q.qual_flags |= QualRecord::QUAL_FLAG_MIX_NO_DATA;
                    q.contig_flags               = std::move(flags);
                    ++n_done;
                    if (n_done % 50000 == 0)
                        spdlog::info("check pass-B: {}/{} genomes done", n_done, n_flagged);
                }
            }
        });

    spdlog::info("check pass-B: complete ({} genomes)", n_flagged);

    // ── CCO baseline subtraction ─────────────────────────────────────────────
    // Per-genus 20th-percentile CCO of reference genomes (from qual_cache) gives
    // the natural genomic-heterogeneity floor. Subtract it from target CCO and
    // clamp to [0,1] to produce contamination_contig_outlier_adj.
    // ── CCO self-calibrating baseline ───────────────────────────────────────
    // Use the 20th percentile of ALL non-NaN CCO values across current targets as
    // the natural-heterogeneity floor. Under the assumption that most genomes in the
    // query set are clean, the low percentile captures the background noise level.
    // This avoids dependency on cached reference QUAL records (which are overwritten
    // by each check run) while still removing the systematic 6% FP offset.
    {
        std::unordered_map<std::string, std::vector<float>> genus_cco;
        for (const auto& [acc, q] : quality) {
            if (std::isnan(q.contamination_contig_outlier)) continue;
            for (const auto& [genus, members] : pass_a.genus_members) {
                bool found = false;
                for (const auto& m : members) { if (m == acc) { found = true; break; } }
                if (!found) continue;
                genus_cco[genus].push_back(q.contamination_contig_outlier);
                break;
            }
        }

        std::unordered_map<std::string, float> genus_baseline;
        for (auto& [genus, vec] : genus_cco) {
            if (vec.size() < 5) continue;
            size_t p20_idx = static_cast<size_t>(0.2f * static_cast<float>(vec.size()));
            if (p20_idx >= vec.size()) p20_idx = 0;
            std::nth_element(vec.begin(), vec.begin() + static_cast<ptrdiff_t>(p20_idx), vec.end());
            genus_baseline[genus] = vec[p20_idx];
        }

        if (!genus_baseline.empty()) {
            spdlog::info("CCO baseline: {} genera with self-calibrated floor", genus_baseline.size());
            for (auto& [acc, q] : quality) {
                if (std::isnan(q.contamination_contig_outlier)) continue;
                for (const auto& [genus, members] : pass_a.genus_members) {
                    auto bit = genus_baseline.find(genus);
                    if (bit == genus_baseline.end()) continue;
                    for (const auto& m : members) {
                        if (m == acc) {
                            q.contamination_contig_outlier_adj = std::max(
                                0.0f, q.contamination_contig_outlier - bit->second);
                            goto next_acc;
                        }
                    }
                }
                next_acc:;
            }
        }
    }
}

} // namespace genopack::check
