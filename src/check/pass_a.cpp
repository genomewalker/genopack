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

    if (n_k == 1) {
        if (!std::isnan(lk[0])) q.contamination_leakage = lk[0];
        return;
    }

    // Excess containment decay at the largest k, relative to a slope-only OLS
    // (log C_k = k*beta, alpha=1 forced) fit on the SMALLER k only. The held-out k must not
    // enter its own fit: at k=[16,21,31] the top k carries 961/1658 = 58% of the weight, so a
    // fit including it predicts itself, the residual collapses, and max(0,.) clamps the rest to
    // zero. That is what this was doing -- exactly 0.0 on 99.97% of 9.53 M genomes, global max
    // 0.005, four times under its own 0.02 firing threshold. The axis could not fire, ever.
    const int last = n_k - 1;
    for (int i = 0; i <= last; ++i)
        if (std::isnan(lk[i])) return;

    const float c_last = 1.0f - lk[last];
    if (c_last < 0.01f) return;

    float num = 0.0f, den = 0.0f;
    for (int i = 0; i < last; ++i) {
        const float c = 1.0f - lk[i];
        if (c < 0.01f) return;
        const float k = static_cast<float>(avail_k[i]);
        num += k * std::log(c);
        den += k * k;
    }
    if (den <= 0.0f) return;

    const float c_pred = std::exp((num / den) * static_cast<float>(avail_k[last]));
    q.contamination_leakage = std::max(0.0f, c_pred - c_last) / std::max(c_pred, 0.01f);
}

} // namespace

PassAResult run_pass_a(ICheckReader& pack,
                       const std::vector<std::string>& accessions,
                       int threads,
                       int min_genus_size,
                       float /*leakage_threshold*/,
                       int min_gstx_members,
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
                    // S is symmetric, so trace(S*S) == ||S||_F^2 — avoid the 136x136 matmul.
                    float trS2 = S.squaredNorm();
                    float alpha  = oas_alpha(trS, trS2, nv, 136);
                    float sigma2 = trS / 136.0f;
                    Eigen::MatrixXf Sig = (1.0f - alpha) * S;
                    Sig.diagonal().array() += alpha * sigma2;  // shrink toward sigma2*I
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

            // When a QUAL cache is provided, load sketch-derived scores from it.
            // TNF and fragmentation are always recomputed (cheap, no I/O).
            bool qual_cache_hit = false;
            if (qual_cache && m->gid > 0) {
                auto it = qual_cache->find(m->gid);
                if (it != qual_cache->end()) {
                    const QualRecord& r = it->second;
                    q.completeness_cluster_relative = r.completeness_cluster_relative;
                    q.contamination_leakage         = r.contamination_leakage;
                    q.contamination_duplication     = QualRecord::decode_dup(r.contamination_duplication_u16);
                    // Phase-2 estimators (raw float; NAN for old-stride packs via reader upcast).
                    q.contamination_core_dup_mass   = r.contamination_core_dup_mass;
                    q.accessory_ratio               = r.accessory_ratio;
                    qual_cache_hit = true;
                }
            }
            // Cache miss (gid absent from QUAL section, e.g. GenusSparse genera not scored at
            // build time) falls through to the length-ratio proxy so completeness_cluster_relative
            // is non-NaN. Without this, pass_b's by_completeness gate never fires and these
            // genomes skip FASTA scoring entirely, leaving all quality fields NaN → stuck LQ.
            if (!qual_cache_hit && sl.S_c > 0.0f) {
                q.completeness_cluster_relative =
                    std::clamp(static_cast<float>(m->len) / sl.S_c, 0.0f, 1.5f);
            }

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
            GenomeId         gid;
            int              gi;       // index into active_genera / slots / genus_tq
            int              tq_idx;   // index within genus_tq[gi]
            const GstxEntry* gstx_e;  // species-level GSTX (or genus fallback)
        };
        std::vector<GlobalTarget> gstx_targets, nogstx_targets;

        // Build acc→species map for target genomes (for species-level GSTX lookup).
        std::unordered_map<std::string, std::string> target_species;
        if (have_sketch) {
            pack.scan_taxonomy_with_id([&](std::string_view acc, std::string_view tax, GenomeId) {
                if (!target_set.count(acc)) return;
                constexpr std::string_view needle = ";s__";
                std::string_view sv;
                auto pos = tax.find(needle);
                if (pos != std::string_view::npos) {
                    auto s = pos + 1, e = tax.find(';', s);
                    sv = tax.substr(s, e == std::string_view::npos ? e : e - s);
                } else if (tax.starts_with("s__")) {
                    sv = tax.substr(0, tax.find(';', 3));
                }
                if (!sv.empty() && sv != "s__")
                    target_species.emplace(std::string(acc), std::string(sv));
            });
        }

        for (int gi = 0; gi < n_active; ++gi) {
            const auto& tq = genus_tq[gi];
            for (int ti = 0; ti < static_cast<int>(tq.size()); ++ti) {
                if (tq[ti].m->gid == 0) continue;
                const GstxEntry* tgt_gstx = nullptr;
                if (have_sketch) {
                    // Prefer species-level GSTX with sufficient members; fall back to genus.
                    auto sit = target_species.find(tq[ti].m->acc);
                    if (sit != target_species.end()) {
                        const GstxEntry* sp = pack.gstx_for_genus(sit->second);
                        if (sp && sp->n_members >= static_cast<uint32_t>(min_gstx_members))
                            tgt_gstx = sp;
                    }
                    if (!tgt_gstx) {
                        const GstxEntry* ge = slots[gi].gstx_e;
                        if (ge && ge->n_members >= static_cast<uint32_t>(min_gstx_members))
                            tgt_gstx = ge;
                    }
                }
                auto& dest = tgt_gstx ? gstx_targets : nogstx_targets;
                dest.push_back({tq[ti].m->gid, gi, ti, tgt_gstx});
            }
        }

        auto sort_by_gid = [](std::vector<GlobalTarget>& v) {
            std::sort(v.begin(), v.end(),
                      [](const GlobalTarget& a, const GlobalTarget& b) { return a.gid < b.gid; });
        };

        // ── 3a: GSTX single-pass ────────────────────────────────────────────────
        if (!gstx_targets.empty()) {
            sort_by_gid(gstx_targets);

            const int n_tgt = static_cast<int>(gstx_targets.size());
            std::vector<float>    lk(static_cast<size_t>(n_tgt) * GSTX_MAX_K, NAN);
            std::vector<uint32_t> nrb(static_cast<size_t>(n_tgt), 0); // n_real_bins at k=0

            std::vector<GenomeId> sorted_gids;
            sorted_gids.reserve(gstx_targets.size());
            for (const auto& gt : gstx_targets) sorted_gids.push_back(gt.gid);

            pack.visit_sketch_gids(sorted_gids, sketch_sz,
                [&](size_t bidx, uint32_t ki, const SketchResult& sk) {
                    if (ki >= static_cast<uint32_t>(n_k)) return;
                    const GlobalTarget& gt  = gstx_targets[bidx];
                    const GstxEntry*    gst = gt.gstx_e;
                    const uint16_t* cons    = gst->consensus[ki];
                    uint32_t mismatch = 0;
                    for (uint32_t b = 0; b < sk.sketch_size; ++b)
                        if (sk.sig[b] != cons[b]) ++mismatch;
                    lk[static_cast<size_t>(bidx) * GSTX_MAX_K + ki] =
                        static_cast<float>(mismatch) / static_cast<float>(sk.sketch_size);
                    if (ki == 0) nrb[bidx] = sk.n_real_bins;
                });

            // Dispersion bump: median/mad are only present in new-format GSTX packs.
            const GstxReader* gstx_rd   = pack.gstx_reader();
            const bool gstx_dispersion  = gstx_rd && gstx_rd->has_dispersion();

            for (int tidx = 0; tidx < n_tgt; ++tidx) {
                const GlobalTarget& gt  = gstx_targets[tidx];
                const GstxEntry*    gst = gt.gstx_e;
                const float* tlk        = &lk[static_cast<size_t>(tidx) * GSTX_MAX_K];
                auto& q = genus_tq[gt.gi][gt.tq_idx].q;
                apply_leakage_scores(q, tlk, n_k, avail_k.data(), gst->p90_containment[0]);

                // Phase-1 relative-conspecific-containment via pre-baked GSTX dispersion.
                // Mirrors the no-GSTX path (accessory_ratio/accessory_z from c0 median+MAD),
                // but reads the median/mad baked into the entry at build time. Gated on
                // GenusSaturated, present dispersion, mad>0, ≥10 members, and a real leakage.
                if (gstx_dispersion &&
                    q.support_tier == SupportTier::GenusSaturated &&
                    gst->n_members >= 10 && !std::isnan(tlk[0])) {
                    const float med = gst->median_containment[0];
                    const float mad = gst->mad_containment[0];
                    if (!std::isnan(med) && !std::isnan(mad) && mad > 0.0f) {
                        const float c0_query = 1.0f - tlk[0];
                        if (med > 0.0f)
                            q.accessory_ratio = c0_query / med;
                        q.accessory_z = (c0_query - med) / mad;
                    }
                }
            }
            spdlog::info("check pass-A: {} targets scored via single GSTX pass",
                         gstx_targets.size());
        }

        // ── 3b: No-GSTX genera — two global sorted passes then parallel post-process
        // Old: 2 × N_genera sequential visit_sketch_gids calls (N≈500 → ~1000 NFS scans)
        // New: 2 global sorted scans (one consensus, one leakage) + parallel per-genus scoring
        if (!nogstx_targets.empty()) {
            spdlog::info("check pass-A: {} targets in no-GSTX genera (legacy path)",
                         nogstx_targets.size());

            static constexpr int kMaxSample = 2000;

            std::unordered_map<uint32_t, int> k_to_ki;
            for (int ki = 0; ki < n_k; ++ki) k_to_ki[avail_k[ki]] = ki;

            std::vector<int> nogstx_gi;
            {
                std::unordered_set<int> seen;
                for (const auto& gt : nogstx_targets)
                    if (seen.insert(gt.gi).second) nogstx_gi.push_back(gt.gi);
            }
            const int n_nogstx = static_cast<int>(nogstx_gi.size());

            struct GenusState {
                std::vector<GenusSignatureVote>   votes;
                std::vector<GenomeId>             sample_gids;
                std::unordered_map<GenomeId, int> sample_pos;
                std::vector<float>                sample_lk;   // [n_sample × GSTX_MAX_K]
                std::vector<uint32_t>             sample_nrb;
                std::unordered_set<GenomeId>      tq_gids;
                std::unordered_map<GenomeId, int> tq_by_gid;
            };
            std::vector<GenusState> gstates(n_nogstx);

            // gid → gstates index for the two global passes
            std::unordered_map<GenomeId, int> gid_to_gsidx;
            std::unordered_map<GenomeId, int> sgid_to_gsidx;

            for (int gii = 0; gii < n_nogstx; ++gii) {
                const int gi             = nogstx_gi[gii];
                const auto& all_members  = genus_all.at(active_genera[gi]);
                const auto& tq           = genus_tq[gi];
                auto& gs                 = gstates[gii];
                const int n              = static_cast<int>(all_members.size());

                gs.votes.reserve(n_k);
                for (int ki = 0; ki < n_k; ++ki) gs.votes.emplace_back(sketch_sz);

                for (int ti = 0; ti < static_cast<int>(tq.size()); ++ti) {
                    if (!tq[ti].m->gid) continue;
                    gs.tq_by_gid[tq[ti].m->gid] = ti;
                    gs.tq_gids.insert(tq[ti].m->gid);
                }

                std::vector<std::pair<GenomeId, std::string>> id_acc;
                id_acc.reserve(n);
                for (const auto& m : all_members)
                    if (m.gid > 0) id_acc.emplace_back(m.gid, m.acc);
                std::sort(id_acc.begin(), id_acc.end(),
                          [](const auto& a, const auto& b) { return a.first < b.first; });

                for (const auto& [gid, _] : id_acc) gid_to_gsidx[gid] = gii;

                for (const auto& [gid, _] : id_acc)
                    if (gs.tq_by_gid.count(gid)) gs.sample_gids.push_back(gid);
                const int budget = kMaxSample - static_cast<int>(gs.sample_gids.size());
                if (budget > 0) {
                    int step = std::max(1, n / budget);
                    for (int i = 0; i < n &&
                             static_cast<int>(gs.sample_gids.size()) < kMaxSample; i += step)
                        if (!gs.tq_by_gid.count(id_acc[i].first))
                            gs.sample_gids.push_back(id_acc[i].first);
                    std::sort(gs.sample_gids.begin(), gs.sample_gids.end());
                }
                const int ns = static_cast<int>(gs.sample_gids.size());
                for (int si = 0; si < ns; ++si) {
                    gs.sample_pos[gs.sample_gids[si]] = si;
                    sgid_to_gsidx[gs.sample_gids[si]] = gii;
                }
                gs.sample_lk.assign(static_cast<size_t>(ns) * GSTX_MAX_K, NAN);
                gs.sample_nrb.assign(static_cast<size_t>(ns), 0);
            }

            // Pass 1: one global sorted scan — update BM consensus for every genus
            {
                std::vector<GenomeId> all_gids;
                all_gids.reserve(gid_to_gsidx.size());
                for (const auto& [gid, _] : gid_to_gsidx) all_gids.push_back(gid);
                std::sort(all_gids.begin(), all_gids.end());

                pack.visit_sketch_gids(all_gids, sketch_sz,
                    [&](size_t bidx, uint32_t ki_raw, const SketchResult& sk) {
                        auto kit = k_to_ki.find(ki_raw);
                        if (kit == k_to_ki.end()) return;
                        auto git = gid_to_gsidx.find(all_gids[bidx]);
                        if (git == gid_to_gsidx.end()) return;
                        gstates[git->second].votes[kit->second].update(sk.sig, sk.sketch_size);
                    });
            }

            // Pass 2: one global sorted scan — compute leakage against per-genus consensus
            {
                std::vector<GenomeId> all_sample;
                all_sample.reserve(sgid_to_gsidx.size());
                for (const auto& [gid, _] : sgid_to_gsidx) all_sample.push_back(gid);
                std::sort(all_sample.begin(), all_sample.end());

                pack.visit_sketch_gids(all_sample, sketch_sz,
                    [&](size_t bidx, uint32_t ki_raw, const SketchResult& sk) {
                        auto kit = k_to_ki.find(ki_raw);
                        if (kit == k_to_ki.end()) return;
                        const int ki   = kit->second;
                        const GenomeId gid = all_sample[bidx];
                        auto git = sgid_to_gsidx.find(gid);
                        if (git == sgid_to_gsidx.end()) return;
                        auto& gs   = gstates[git->second];
                        auto sit   = gs.sample_pos.find(gid);
                        if (sit == gs.sample_pos.end()) return;
                        const int si = sit->second;
                        gs.sample_lk[static_cast<size_t>(si) * GSTX_MAX_K + ki] =
                            gs.votes[ki].leakage(sk.sig, sk.sketch_size);
                        if (ki == 0) gs.sample_nrb[si] = sk.n_real_bins;
                    });
            }

            // Post-process: p90 + score application — no I/O, parallel across genera
            #pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
            for (int gii = 0; gii < n_nogstx; ++gii) {
                auto& gs     = gstates[gii];
                auto& tq     = genus_tq[nogstx_gi[gii]];
                const int ns = static_cast<int>(gs.sample_gids.size());

                std::vector<float> c0_all, nrb_ref;
                c0_all.reserve(ns);
                nrb_ref.reserve(ns);
                for (int si = 0; si < ns; ++si) {
                    const float v = gs.sample_lk[static_cast<size_t>(si) * GSTX_MAX_K + 0];
                    if (!std::isnan(v)) c0_all.push_back(1.0f - v);
                    if (!gs.tq_gids.count(gs.sample_gids[si]) && gs.sample_nrb[si] > 0)
                        nrb_ref.push_back(static_cast<float>(gs.sample_nrb[si]));
                }
                float p90_c0 = 0.0f;
                if (!c0_all.empty()) {
                    size_t idx = static_cast<size_t>(
                        std::ceil(0.9f * static_cast<float>(c0_all.size())) - 1);
                    std::nth_element(c0_all.begin(),
                                     c0_all.begin() + static_cast<ptrdiff_t>(idx),
                                     c0_all.end());
                    p90_c0 = c0_all[idx];
                }
                float nrb_p90 = 0.0f;
                if (!nrb_ref.empty()) {
                    size_t idx = static_cast<size_t>(
                        0.9f * static_cast<float>(nrb_ref.size()));
                    if (idx >= nrb_ref.size()) idx = nrb_ref.size() - 1;
                    std::nth_element(nrb_ref.begin(),
                                     nrb_ref.begin() + static_cast<ptrdiff_t>(idx),
                                     nrb_ref.end());
                    nrb_p90 = nrb_ref[idx];
                }

                // Phase-1 relative-conspecific-containment: median + MAD of the genus
                // OPH-containment distribution (c0_all). Only when the genus is saturated
                // with ≥10 members and MAD>0; else abstain (NaN). c0_all was mutated by the
                // p90 nth_element above, so recompute the order statistics from scratch here.
                float c0_median = NAN, c0_mad = 0.0f;
                const bool acc_ok =
                    !tq.empty() && tq.front().q.support_tier == SupportTier::GenusSaturated &&
                    c0_all.size() >= 10;
                if (acc_ok) {
                    std::vector<float> tmp = c0_all;
                    const size_t mid = tmp.size() / 2;
                    std::nth_element(tmp.begin(), tmp.begin() + mid, tmp.end());
                    c0_median = tmp[mid];
                    for (float& v : tmp) v = std::fabs(v - c0_median);
                    std::nth_element(tmp.begin(), tmp.begin() + mid, tmp.end());
                    c0_mad = 1.4826f * tmp[mid];
                }

                for (auto& tqr : tq) {
                    if (!tqr.m->gid) continue;
                    auto it = gs.sample_pos.find(tqr.m->gid);
                    if (it == gs.sample_pos.end()) continue;
                    const float* tlk = &gs.sample_lk[static_cast<size_t>(it->second) * GSTX_MAX_K];
                    apply_leakage_scores(tqr.q, tlk, n_k, avail_k.data(), p90_c0);
                    if (acc_ok && c0_mad > 0.0f && !std::isnan(tlk[0])) {
                        const float c0_query = 1.0f - tlk[0];
                        if (c0_median > 0.0f)
                            tqr.q.accessory_ratio = c0_query / c0_median;
                        tqr.q.accessory_z = (c0_query - c0_median) / c0_mad;
                    }
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
        for (auto& [fam, genera_set] : fam_genera) {
            auto& v = result.family_to_genera[fam];
            v.assign(genera_set.begin(), genera_set.end());
            std::sort(v.begin(), v.end());
        }
    }

    spdlog::info("check pass-A: complete — {} genera processed", active_genera.size());
    return result;
}

} // namespace genopack::check
