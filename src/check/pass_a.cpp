#include "pass_a.hpp"
#include <genopack/gstx.hpp>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <mutex>
#include <omp.h>
#include <spdlog/spdlog.h>
#include <unordered_set>

namespace genopack::check {

namespace {

float percentile_90(std::vector<uint64_t> lengths) {
    if (lengths.empty()) return 0.0f;
    std::sort(lengths.begin(), lengths.end());
    size_t idx = static_cast<size_t>(std::ceil(0.9f * lengths.size()));
    if (idx > 0) --idx;
    return static_cast<float>(lengths[std::min(idx, lengths.size() - 1)]);
}

float fragmentation_score(uint32_t n) {
    if (n <= 1) return 1.0f;
    return 1.0f / (1.0f + 0.333f * std::log2f(static_cast<float>(n)));
}

float oas_alpha(float trS, float trS2, int n, int p) {
    float num = (1.0f - 2.0f / p) * trS2 + trS * trS;
    float den = (static_cast<float>(n) + 1.0f - 2.0f / p) * (trS2 - trS * trS / p);
    if (std::abs(den) < 1e-12f) return 1.0f;
    return std::clamp(num / den, 0.0f, 1.0f);
}

struct GenusSignatureVote {
    std::vector<uint16_t> candidate;
    std::vector<int32_t>  count;

    explicit GenusSignatureVote(uint32_t sketch_sz)
        : candidate(sketch_sz, 0), count(sketch_sz, 0) {}

    GenusSignatureVote(const uint16_t* prebuilt, uint32_t sketch_sz)
        : candidate(prebuilt, prebuilt + sketch_sz), count(sketch_sz, 1) {}

    void update(const uint16_t* sig, uint32_t sketch_sz) {
        for (uint32_t b = 0; b < sketch_sz; ++b) {
            if (count[b] == 0)             { candidate[b] = sig[b]; count[b] = 1; }
            else if (sig[b] == candidate[b]) { ++count[b]; }
            else                             { --count[b]; }
        }
    }

    float leakage(const uint16_t* sig, uint32_t sketch_sz) const {
        uint32_t mismatch = 0;
        for (uint32_t b = 0; b < sketch_sz; ++b)
            if (sig[b] != candidate[b]) ++mismatch;
        return static_cast<float>(mismatch) / static_cast<float>(sketch_sz);
    }
};

// Score contamination/consistency from per-k leakage values.
void apply_leakage_scores(GenomeQuality& q,
                          const float* lk, int n_k,
                          const uint32_t* avail_k,
                          float p90_c0)
{
    if (p90_c0 > 0.01f && !std::isnan(lk[0]))
        q.completeness_cluster_relative = std::clamp(
            (1.0f - lk[0]) / p90_c0, 0.0f, 1.5f);

    // Slope-only OLS: log(C_k) = k * β (α=1 forced), β = Σ(k_i*log(c_i))/Σ(k_i²)
    if (n_k >= 3 && !std::isnan(lk[0]) && !std::isnan(lk[1]) && !std::isnan(lk[2])) {
        const float c0 = 1.0f - lk[0], c1 = 1.0f - lk[1], c2 = 1.0f - lk[2];
        if (c0 >= 0.01f && c1 >= 0.01f && c2 >= 0.01f) {
            const float k0 = static_cast<float>(avail_k[0]);
            const float k1 = static_cast<float>(avail_k[1]);
            const float k2 = static_cast<float>(avail_k[2]);
            const float lc0 = std::log(c0), lc1 = std::log(c1), lc2 = std::log(c2);
            const float beta    = (k0*lc0 + k1*lc1 + k2*lc2) / (k0*k0 + k1*k1 + k2*k2);
            const float c2_pred = std::exp(beta * k2);
            q.contamination_leakage = std::max(0.0f, c2_pred - c2) / std::max(c2_pred, 0.01f);
            const float r0=lc0-beta*k0, r1=lc1-beta*k1, r2=lc2-beta*k2;
            q.leakage_residual = std::sqrt((r0*r0 + r1*r1 + r2*r2) / 2.0f);
        }
    } else if (n_k == 2 && !std::isnan(lk[0]) && !std::isnan(lk[1])) {
        const float c0 = 1.0f - lk[0], c1 = 1.0f - lk[1];
        if (c0 >= 0.01f && c1 >= 0.01f) {
            const float k0 = static_cast<float>(avail_k[0]);
            const float k1 = static_cast<float>(avail_k[1]);
            const float beta    = (k0*std::log(c0) + k1*std::log(c1)) / (k0*k0 + k1*k1);
            const float c1_pred = std::exp(beta * k1);
            q.contamination_leakage = std::max(0.0f, c1_pred - c1) / std::max(c1_pred, 0.01f);
        }
    } else if (n_k == 1 && !std::isnan(lk[0])) {
        q.contamination_leakage = lk[0];
    }
}

} // namespace

PassAResult run_pass_a(ICheckReader& pack,
                       const std::vector<std::string>& accessions,
                       int threads,
                       int min_genus_size,
                       float /*leakage_threshold*/,
                       const std::unordered_map<uint64_t, QualRecord>* qual_cache)
{
    struct MemberMeta {
        std::string acc;
        GenomeId    gid   = 0;
        uint64_t    len   = 0;
        uint32_t    nctg  = 0;
    };
    std::unordered_map<std::string, std::vector<MemberMeta>> genus_all;
    genus_all.reserve(8192);
    std::unordered_map<std::string, std::string> acc_to_family;

    struct SvHash {
        using is_transparent = void;
        size_t operator()(std::string_view s) const noexcept
            { return std::hash<std::string_view>{}(s); }
        size_t operator()(const std::string& s) const noexcept
            { return std::hash<std::string_view>{}(s); }
    };
    const std::unordered_set<std::string, SvHash, std::equal_to<>> target_set(
        accessions.begin(), accessions.end());

    // Optional genus pre-filter for small target sets
    std::unordered_set<std::string, SvHash, std::equal_to<>> target_genera_filter;
    if (accessions.size() < 100000) {
        pack.scan_taxonomy_with_id([&](std::string_view acc, std::string_view tax, GenomeId) {
            if (!target_set.count(acc)) return;
            constexpr std::string_view needle = ";g__";
            auto pos = tax.find(needle);
            std::string_view genus_sv;
            if (pos != std::string_view::npos) {
                auto start = pos + 1;
                genus_sv   = tax.substr(start, tax.find(';', start) - start);
            } else if (tax.starts_with("g__")) {
                auto end = tax.find(';', 3);
                genus_sv  = tax.substr(0, end);
            }
            if (!genus_sv.empty()) target_genera_filter.emplace(genus_sv);
        });
    }

    pack.scan_taxonomy_with_id([&](std::string_view acc, std::string_view tax, GenomeId gid) {
        constexpr std::string_view needle = ";g__";
        auto pos = tax.find(needle);
        std::string_view genus_sv;
        if (pos != std::string_view::npos) {
            auto start = pos + 1;
            auto end   = tax.find(';', start);
            genus_sv   = tax.substr(start, end == std::string_view::npos ? end : end - start);
        } else if (tax.starts_with("g__")) {
            auto end = tax.find(';', 3);
            genus_sv  = tax.substr(0, end);
        }
        if (!target_genera_filter.empty() && !target_genera_filter.count(genus_sv)) return;
        MemberMeta m{std::string(acc), gid, 0, 0};
        if (gid != 0) {
            if (auto mo = pack.genome_meta_by_accession(acc)) {
                m.len  = mo->genome_length;
                m.nctg = mo->n_contigs;
            }
        }
        genus_all[std::string(genus_sv)].push_back(std::move(m));

        // Extract family for sibling contamination scoring in pass-B
        constexpr std::string_view fneedle = ";f__";
        auto fpos = tax.find(fneedle);
        std::string_view family_sv;
        if (fpos != std::string_view::npos) {
            auto fstart = fpos + 1;
            auto fend   = tax.find(';', fstart);
            family_sv   = tax.substr(fstart, fend == std::string_view::npos ? fend : fend - fstart);
        } else if (tax.starts_with("f__")) {
            auto fend = tax.find(';', 3);
            family_sv = tax.substr(0, fend);
        }
        if (!family_sv.empty() && family_sv != "f__")
            acc_to_family.emplace(std::string(acc), std::string(family_sv));
    });

    std::unordered_map<std::string, std::vector<const MemberMeta*>> genus_targets;
    for (const auto& [g, members] : genus_all)
        for (const auto& m : members)
            if (target_set.count(m.acc))
                genus_targets[g].push_back(&m);

    std::vector<std::string> active_genera;
    active_genera.reserve(genus_targets.size());
    for (const auto& [g, _] : genus_targets) active_genera.push_back(g);
    const int n_active = static_cast<int>(active_genera.size());

    const bool have_tnf    = !accessions.empty() &&
                             pack.kmer_profile_by_accession(accessions[0]) != nullptr;
    const bool have_sketch = pack.has_sketches();
    const std::vector<uint32_t> avail_k =
        have_sketch ? pack.available_kmer_sizes() : std::vector<uint32_t>{};
    const uint32_t sketch_sz = have_sketch ? pack.sketch_sketch_size() : 0u;
    const int n_k = static_cast<int>(avail_k.size());

    {
        std::string ks;
        for (size_t i = 0; i < avail_k.size(); ++i) {
            if (i) ks += ',';
            ks += std::to_string(avail_k[i]);
        }
        spdlog::info("check pass-A: sketch_size={} k=[{}] n_k={}", sketch_sz, ks, n_k);
    }
    spdlog::info("check pass-A: {} target genomes, {} relevant genera (of {} total)",
                 accessions.size(), active_genera.size(), genus_all.size());

    constexpr float kExpMahal = 11.63f;

    // ── Phase 1: TNF + genus metadata (CPU-bound, parallel, no SKCH I/O) ──────
    struct GenusSlot {
        SupportTier      tier;
        float            interval_width;
        float            S_c;
        const GstxEntry* gstx_e = nullptr;
        bool             has_tnf = false;
        Eigen::VectorXf  tnf_mu;
        Eigen::MatrixXf  tnf_L;
    };
    std::vector<GenusSlot> slots(n_active);

    #pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (int gi = 0; gi < n_active; ++gi) {
        const std::string& genus = active_genera[gi];
        const auto& all_members  = genus_all.at(genus);
        const int n              = static_cast<int>(all_members.size());
        auto& sl = slots[gi];

        sl.tier = (genus.empty() || genus == "g__") ? SupportTier::Singleton
                : (n < min_genus_size)              ? SupportTier::GenusSparse
                                                    : SupportTier::GenusSaturated;
        sl.interval_width = (sl.tier == SupportTier::GenusSaturated) ? 0.05f
                          : (sl.tier == SupportTier::GenusSparse)    ? 0.20f
                                                                      : 0.50f;

        std::vector<uint64_t> lengths;
        lengths.reserve(n);
        for (const auto& m : all_members) if (m.len > 0) lengths.push_back(m.len);
        sl.S_c = percentile_90(lengths);

        if (have_sketch) {
            sl.gstx_e = pack.gstx_for_genus(genus);
            if (sl.gstx_e && sl.gstx_e->n_k_stored < static_cast<uint8_t>(n_k))
                sl.gstx_e = nullptr;
        }

        static constexpr int kMaxTnfMembers = 5000;
        if (have_tnf && n >= 2) {
            std::vector<const float*> ptrs;
            ptrs.reserve(std::min(n, kMaxTnfMembers));
            int step = std::max(1, n / kMaxTnfMembers);
            for (int mi = 0; mi < n; mi += step) {
                const float* p = pack.kmer_profile_by_accession(all_members[mi].acc);
                if (p) ptrs.push_back(p);
            }
            const int nv = static_cast<int>(ptrs.size());
            if (nv >= 2) {
                Eigen::MatrixXf X(nv, 136);
                for (int r = 0; r < nv; ++r)
                    X.row(r) = Eigen::Map<const Eigen::VectorXf>(ptrs[r], 136).transpose();
                sl.tnf_mu = X.colwise().mean();
                if (nv >= 3) {
                    Eigen::MatrixXf Xc = X.rowwise() - sl.tnf_mu.transpose();
                    Eigen::MatrixXf S  = (Xc.transpose() * Xc) / static_cast<float>(nv - 1);
                    float trS  = S.trace();
                    float trS2 = (S * S).trace();
                    float alpha  = oas_alpha(trS, trS2, nv, 136);
                    float sigma2 = trS / 136.0f;
                    Eigen::MatrixXf Sig = (1.0f - alpha) * S
                                       + alpha * sigma2 * Eigen::MatrixXf::Identity(136, 136);
                    Eigen::LLT<Eigen::MatrixXf> llt(Sig);
                    if (llt.info() == Eigen::Success) {
                        sl.tnf_L   = llt.matrixL();
                        sl.has_tnf = true;
                    }
                }
            }
        }
    }

    // ── Phase 2: Build per-genus TargetQual lists (TNF + fragmentation scores) ──
    struct TargetQual {
        const MemberMeta* m;
        GenomeQuality     q;
    };
    std::vector<std::vector<TargetQual>> genus_tq(n_active);

    for (int gi = 0; gi < n_active; ++gi) {
        const std::string& genus = active_genera[gi];
        const auto& targets      = genus_targets.at(genus);
        const auto& sl           = slots[gi];
        auto& tq                 = genus_tq[gi];
        tq.reserve(targets.size());
        for (const MemberMeta* m : targets) {
            GenomeQuality q;
            q.support_tier   = sl.tier;
            q.interval_width = sl.interval_width;

            // When a QUAL cache is provided, load sketch-derived scores from it.
            // TNF and fragmentation are always recomputed (cheap, no I/O).
            if (qual_cache && m->gid > 0) {
                auto it = qual_cache->find(m->gid);
                if (it != qual_cache->end()) {
                    const QualRecord& r = it->second;
                    q.completeness_cluster_relative = r.completeness_cluster_relative;
                    q.contamination_leakage         = r.contamination_leakage;
                    q.leakage_residual              = r.leakage_residual;
                }
            } else if (sl.S_c > 0.0f) {
                q.completeness_cluster_relative =
                    std::clamp(static_cast<float>(m->len) / sl.S_c, 0.0f, 1.5f);
            }

            q.completeness_fragmentation = fragmentation_score(m->nctg);
            if (sl.has_tnf) {
                const float* p = pack.kmer_profile_by_accession(m->acc);
                if (p) {
                    Eigen::VectorXf x = Eigen::Map<const Eigen::VectorXf>(p, 136);
                    float dist = sl.tnf_L.triangularView<Eigen::Lower>()
                                         .solve(x - sl.tnf_mu).norm();
                    q.contamination_tnf_excess =
                        std::max(0.0f, (dist - kExpMahal) / kExpMahal);
                }
            }
            tq.push_back({m, std::move(q)});
        }
    }

    // ── Phase 3: Single sequential SKCH pass (GSTX genera) ─────────────────────
    // Skipped when a qual_cache is provided — sketch-derived scores already loaded.
    if (qual_cache) {
        spdlog::info("check pass-A: using cached QUAL scores, SKCH scan skipped");
    } else
    // Genera with GSTX: consensus + p90 already known. One global sorted-gid scan
    // to compute leakage for each target, then derive all scores in one post-pass.
    // Genera without GSTX: legacy per-genus path (serial, expected rare/empty).

    if (have_sketch && n_k > 0) {
        // Build global sorted target gid list, split by GSTX availability
        struct GlobalTarget {
            GenomeId gid;
            int      gi;      // index into active_genera / slots / genus_tq
            int      tq_idx;  // index within genus_tq[gi]
        };
        std::vector<GlobalTarget> gstx_targets, nogstx_targets;

        for (int gi = 0; gi < n_active; ++gi) {
            const auto& tq = genus_tq[gi];
            auto& dest = slots[gi].gstx_e ? gstx_targets : nogstx_targets;
            for (int ti = 0; ti < static_cast<int>(tq.size()); ++ti)
                if (tq[ti].m->gid > 0)
                    dest.push_back({tq[ti].m->gid, gi, ti});
        }

        auto sort_by_gid = [](std::vector<GlobalTarget>& v) {
            std::sort(v.begin(), v.end(),
                      [](const GlobalTarget& a, const GlobalTarget& b) { return a.gid < b.gid; });
        };

        // ── 3a: GSTX single-pass ────────────────────────────────────────────────
        if (!gstx_targets.empty()) {
            sort_by_gid(gstx_targets);

            // Per-target leakage storage: lk[target_idx * GSTX_MAX_K + ki]
            // visit_sketch_gids delivers ki as index into available_kmer_sizes() (0,1,2…),
            // NOT as the raw k-value — use directly as array index.
            const int n_tgt = static_cast<int>(gstx_targets.size());
            std::vector<float> lk(static_cast<size_t>(n_tgt) * GSTX_MAX_K, NAN);

            std::vector<GenomeId> sorted_gids;
            sorted_gids.reserve(gstx_targets.size());
            for (const auto& gt : gstx_targets) sorted_gids.push_back(gt.gid);

            pack.visit_sketch_gids(sorted_gids, sketch_sz,
                [&](size_t bidx, uint32_t ki, const SketchResult& sk) {
                    if (ki >= static_cast<uint32_t>(n_k)) return;
                    const GlobalTarget& gt  = gstx_targets[bidx];
                    const GstxEntry*    gst = slots[gt.gi].gstx_e;
                    const uint16_t* cons    = gst->consensus[ki];
                    uint32_t mismatch = 0;
                    for (uint32_t b = 0; b < sk.sketch_size; ++b)
                        if (sk.sig[b] != cons[b]) ++mismatch;
                    lk[static_cast<size_t>(bidx) * GSTX_MAX_K + ki] =
                        static_cast<float>(mismatch) / static_cast<float>(sk.sketch_size);
                });

            for (int tidx = 0; tidx < n_tgt; ++tidx) {
                const GlobalTarget& gt = gstx_targets[tidx];
                const GenusSlot&    sl = slots[gt.gi];
                const float* tlk       = &lk[static_cast<size_t>(tidx) * GSTX_MAX_K];
                apply_leakage_scores(genus_tq[gt.gi][gt.tq_idx].q,
                                     tlk, n_k, avail_k.data(),
                                     sl.gstx_e->p90_containment[0]);
            }
            spdlog::info("check pass-A: {} targets scored via single GSTX pass",
                         gstx_targets.size());
        }

        // ── 3b: Legacy per-genus path for no-GSTX genera (serial) ───────────────
        if (!nogstx_targets.empty()) {
            spdlog::info("check pass-A: {} targets in no-GSTX genera (legacy path)",
                         nogstx_targets.size());

            std::unordered_map<uint32_t, int> k_to_ki;
            for (int ki = 0; ki < n_k; ++ki) k_to_ki[avail_k[ki]] = ki;

            // Collect unique genus indices with no GSTX
            std::vector<int> nogstx_gi;
            {
                std::unordered_set<int> seen;
                for (const auto& gt : nogstx_targets)
                    if (seen.insert(gt.gi).second) nogstx_gi.push_back(gt.gi);
            }

            for (int gi : nogstx_gi) {
                const std::string& genus  = active_genera[gi];
                const auto& all_members   = genus_all.at(genus);
                auto& tq                  = genus_tq[gi];
                const int n               = static_cast<int>(all_members.size());

                // Sort all member gids for sequential frame access
                std::vector<std::pair<GenomeId, std::string>> id_acc;
                id_acc.reserve(n);
                for (const auto& m : all_members)
                    if (m.gid > 0) id_acc.emplace_back(m.gid, m.acc);
                std::sort(id_acc.begin(), id_acc.end(),
                          [](const auto& a, const auto& b) { return a.first < b.first; });

                std::vector<GenomeId> all_gids;
                all_gids.reserve(id_acc.size());
                for (const auto& [gid, _] : id_acc) all_gids.push_back(gid);

                // BM consensus pass
                std::vector<GenusSignatureVote> votes;
                votes.reserve(n_k);
                for (int ki = 0; ki < n_k; ++ki) votes.emplace_back(sketch_sz);

                pack.visit_sketch_gids(all_gids, sketch_sz,
                    [&](size_t, uint32_t ki_raw, const SketchResult& sk) {
                        auto it = k_to_ki.find(ki_raw);
                        if (it == k_to_ki.end()) return;
                        votes[it->second].update(sk.sig, sk.sketch_size);
                    });

                // Sample for scoring (targets + uniform sample up to 2000)
                static constexpr int kMaxSample = 2000;
                std::unordered_map<GenomeId, int> tq_by_gid;
                for (int ti = 0; ti < static_cast<int>(tq.size()); ++ti)
                    if (tq[ti].m->gid > 0) tq_by_gid[tq[ti].m->gid] = ti;

                std::vector<GenomeId> sample_gids;
                for (const auto& [gid, _] : id_acc)
                    if (tq_by_gid.count(gid)) sample_gids.push_back(gid);
                const int budget = kMaxSample - static_cast<int>(sample_gids.size());
                if (budget > 0) {
                    int step = std::max(1, static_cast<int>(id_acc.size()) / budget);
                    for (int i = 0; i < static_cast<int>(id_acc.size()) &&
                             static_cast<int>(sample_gids.size()) < kMaxSample; i += step)
                        if (!tq_by_gid.count(id_acc[i].first))
                            sample_gids.push_back(id_acc[i].first);
                    std::sort(sample_gids.begin(), sample_gids.end());
                }
                const int n_sample = static_cast<int>(sample_gids.size());

                // Leakage pass over sample
                std::vector<float> sample_lk(
                    static_cast<size_t>(n_sample) * GSTX_MAX_K, NAN);

                pack.visit_sketch_gids(sample_gids, sketch_sz,
                    [&](size_t sidx, uint32_t ki_raw, const SketchResult& sk) {
                        auto it = k_to_ki.find(ki_raw);
                        if (it == k_to_ki.end()) return;
                        const int ki = it->second;
                        sample_lk[sidx * GSTX_MAX_K + ki] =
                            votes[ki].leakage(sk.sig, sk.sketch_size);
                    });

                // p90 from k=0 containment across full sample
                std::vector<float> c0_all;
                c0_all.reserve(n_sample);
                for (int si = 0; si < n_sample; ++si) {
                    float v = sample_lk[static_cast<size_t>(si) * GSTX_MAX_K + 0];
                    if (!std::isnan(v)) c0_all.push_back(1.0f - v);
                }
                float p90_c0 = 0.0f;
                if (!c0_all.empty()) {
                    size_t p90_idx = static_cast<size_t>(
                        std::ceil(0.9f * static_cast<float>(c0_all.size())) - 1);
                    std::nth_element(c0_all.begin(),
                                     c0_all.begin() + static_cast<ptrdiff_t>(p90_idx),
                                     c0_all.end());
                    p90_c0 = c0_all[p90_idx];
                }

                // Apply scores to target genomes
                std::unordered_map<GenomeId, int> sample_pos;
                for (int si = 0; si < n_sample; ++si)
                    sample_pos[sample_gids[si]] = si;

                for (auto& tqr : tq) {
                    if (!tqr.m->gid) continue;
                    auto it = sample_pos.find(tqr.m->gid);
                    if (it == sample_pos.end()) continue;
                    const float* tlk = &sample_lk[static_cast<size_t>(it->second) * GSTX_MAX_K];
                    apply_leakage_scores(tqr.q, tlk, n_k, avail_k.data(), p90_c0);
                }
            }
        }
    }

    // ── Phase 4: Assemble quality map ───────────────────────────────────────────
    std::unordered_map<std::string, GenomeQuality> quality;
    quality.reserve(accessions.size());
    for (int gi = 0; gi < n_active; ++gi)
        for (auto& tqr : genus_tq[gi])
            quality[tqr.m->acc] = std::move(tqr.q);

    PassAResult result;
    result.accessions = accessions;
    result.quality    = std::move(quality);
    for (auto& [g, members] : genus_all) {
        auto& vec = result.genus_members[g];
        vec.reserve(members.size());
        for (const auto& m : members) vec.push_back(m.acc);
    }
    for (auto& [acc, family] : acc_to_family)
        result.family_members[family].push_back(acc);

    // Build family → genera mapping for containment split in pass-B.
    {
        std::unordered_map<std::string, std::string> acc_to_genus;
        for (const auto& [g, members] : genus_all)
            for (const auto& m : members)
                acc_to_genus[m.acc] = g;
        std::unordered_map<std::string, std::unordered_set<std::string>> fam_genera;
        for (const auto& [acc, fam] : acc_to_family) {
            auto git = acc_to_genus.find(acc);
            if (git != acc_to_genus.end())
                fam_genera[fam].insert(git->second);
        }
        for (auto& [fam, genera_set] : fam_genera)
            result.family_to_genera[fam] = {genera_set.begin(), genera_set.end()};
    }

    spdlog::info("check pass-A: complete — {} genera processed", active_genera.size());
    return result;
}

} // namespace genopack::check
