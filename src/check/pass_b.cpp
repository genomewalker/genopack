#include "pass_b.hpp"
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
#include <vector>

namespace genopack::check {

namespace {


struct ContigRecord {
    std::string header;
    std::string seq;
};

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

bool compute_tnf(std::string_view seq, Eigen::VectorXf& out) {
    if (seq.size() < 5000) return false;
    out = Eigen::VectorXf::Zero(136);

    const uint8_t* b2b = kTnfTables.base2bit.data();
    uint32_t counts[256] = {};

    // Rolling 2-bit index: one lookup + shift/or/mask per base instead of 4 lookups.
    // Reset window on any N/ambiguous base.
    uint32_t idx   = 0;
    uint32_t valid = 0; // consecutive valid bases in current window
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

    float norm = out.norm();
    if (norm < 1e-8f) return false;
    out /= norm;
    return true;
}

} // namespace

void run_pass_b(ICheckReader& pack,
                const PassAResult& pass_a,
                std::unordered_map<std::string, GenomeQuality>& quality,
                const PassBConfig& cfg,
                int threads,
                const std::unordered_map<uint64_t, QualRecord>* qual_cache)
{
    std::vector<std::string> flagged;
    flagged.reserve(pass_a.accessions.size() / 10);
    for (const auto& acc : pass_a.accessions) {
        auto it = quality.find(acc);
        if (it == quality.end()) continue;
        const auto& q = it->second;
        bool needs_b = (q.contamination_leakage > cfg.contamination_flag_threshold) ||
                       (q.contamination_tnf_excess  > cfg.tnf_flag_threshold)        ||
                       (!std::isnan(q.completeness_cluster_relative) &&
                        q.completeness_cluster_relative < cfg.completeness_flag_threshold);
        if (needs_b) flagged.push_back(acc);
    }

    if (flagged.empty()) return;
    spdlog::info("check pass-B: {} genomes flagged for FASTA-level analysis", flagged.size());

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
            if (std::isnan(r.chromosome_skew_closure) || std::isnan(r.completeness_post_decontam) ||
                r.self_coherence == 0.0f) {
                to_scan.push_back(acc); continue;
            }
            auto& q = quality.at(acc);
            q.chromosome_skew_closure    = r.chromosome_skew_closure;
            q.completeness_post_decontam = r.completeness_post_decontam;
            q.self_coherence             = r.self_coherence;
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

    std::unordered_map<std::string, Eigen::VectorXf> genus_centroid;
    {
        std::unordered_map<std::string, std::string> acc_to_genus;
        for (const auto& [genus, members] : pass_a.genus_members)
            for (const auto& acc : members)
                acc_to_genus[acc] = genus;

        std::unordered_map<std::string, std::pair<Eigen::VectorXf, int>> sums;
        for (const auto& acc : pass_a.accessions) {
            const float* p = pack.kmer_profile_by_accession(acc);
            if (!p) continue;
            auto git = acc_to_genus.find(acc);
            if (git == acc_to_genus.end()) continue;
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

    std::unordered_map<std::string, std::string> flagged_genus;
    for (const auto& [genus, members] : pass_a.genus_members)
        for (const auto& acc : members)
            flagged_genus[acc] = genus;

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

                std::vector<ContigFlag> flags;
                uint32_t byte_offset    = 0;
                uint32_t clean_len      = 0;
                uint32_t total_len      = 0;
                uint32_t gcov_out_bp    = 0;
                uint32_t gcov_scored_bp = 0;
                uint32_t spe_out_bp     = 0;

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
                    }
                } else {
                    flags.reserve(contigs.size());
                    for (const auto& contig : contigs) {
                        ContigFlag cf;
                        cf.contig_offset = byte_offset;
                        cf.contig_length = static_cast<uint32_t>(contig.seq.size());
                        byte_offset += cf.contig_length;
                        total_len   += cf.contig_length;

                        if (contig.seq.size() >= 5000) {
                            Eigen::VectorXf tnf;
                            if (compute_tnf(contig.seq, tnf)) {
                                Eigen::VectorXf d = tnf - *centroid;
                                float norm_c = centroid->norm();
                                cf.tnf_score = (norm_c > 1e-8f) ? d.norm() / norm_c : NAN;

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
                                }
                            }
                        }

                        bool contig_clean = !std::isnan(cf.tnf_score) &&
                                            cf.tnf_score < cfg.contig_tnf_threshold;
                        if (contig_clean || std::isnan(cf.tnf_score))
                            clean_len += cf.contig_length;

                        flags.push_back(cf);
                    }
                }

                float completeness_post = total_len > 0
                    ? static_cast<float>(clean_len) / static_cast<float>(total_len)
                    : NAN;

                // Per-contig outlier detection.
                // Primary: GCOV archive-calibrated percentile (no self-reference).
                // Fallback: within-genome median/MAD z-score (when no GCOV entry).
                float contamination_contig_outlier = NAN;
                float contamination_spe            = NAN;
                if (gcov_entry && gcov_rd && gcov_scored_bp > 0) {
                    contamination_contig_outlier =
                        static_cast<float>(gcov_out_bp) / static_cast<float>(gcov_scored_bp);
                    contamination_spe =
                        static_cast<float>(spe_out_bp) / static_cast<float>(gcov_scored_bp);
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

                {
                    AllSignals sig = compute_all_signals(eg.fasta, 5000, 16384, 5);

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
}

} // namespace genopack::check
