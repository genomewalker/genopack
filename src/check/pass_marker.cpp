#include "pass_marker.hpp"
#include "pass_b.hpp"
#include <genopack/archive.hpp>
#include <genopack/markers.hpp>
#include <genopack/metamer.hpp>
#include <genopack/gcov.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cmath>
#include <mutex>

namespace genopack::check {

void run_pass_marker(ICheckReader&                                   pack,
                     const PassAResult&                              pass_a,
                     std::unordered_map<std::string, GenomeQuality>& quality,
                     const PassMarkerConfig&                         cfg,
                     int                                             threads)
{
    if (cfg.markers_path.empty()) return;

    // acc → genus (same build as pass_b's flagged_genus)
    std::unordered_map<std::string, std::string> acc_genus;
    for (const auto& [genus, members] : pass_a.genus_members)
        for (const auto& acc : members)
            acc_genus[acc] = genus;

    // Only score genomes still missing marker_completeness after pass_b
    std::vector<std::string> to_score;
    to_score.reserve(pass_a.accessions.size());
    for (const auto& acc : pass_a.accessions) {
        if (!acc_genus.count(acc)) continue;
        auto it = quality.find(acc);
        if (it == quality.end() || std::isnan(it->second.marker_completeness))
            to_score.push_back(acc);
    }
    if (to_score.empty()) return;
    spdlog::info("check pass-Marker: {} genomes to score", to_score.size());

    MarkerReader mrk_rd;
    try {
        mrk_rd.open(cfg.markers_path);
    } catch (const std::exception& ex) {
        spdlog::warn("check pass-Marker: cannot open markers ({}); skipping", ex.what());
        return;
    }
    mrk_rd.build_merged_pool();
    if (!mrk_rd.is_open() || !mrk_rd.has_merged_pool()) return;

    spdlog::info("check pass-Marker: merged pool {} M bac hashes, D6={}",
                 mrk_rd.merged_hashes_bac().size() / 1'000'000,
                 mrk_rd.is_dayhoff6() ? "yes" : "no");

    const bool     is_d6   = mrk_rd.is_dayhoff6();
    const uint64_t frmax   = mrk_rd.frac_max_hash();
    const uint32_t minhit  = static_cast<uint32_t>(cfg.marker_min_hits);
    const int      min_seg = is_d6 ? 50 : METAMER_K;

    std::mutex quality_mtx;
    constexpr int k_readers    = 8;
    const int     inner_threads = std::max(1, threads / k_readers);

    pack.visit_shard_batches_parallel(to_score, k_readers,
        [&](ArchiveReader::ShardBatch& batch) {
            const int n = static_cast<int>(batch.size());
            #pragma omp parallel for schedule(dynamic,1) num_threads(inner_threads)
            for (int j = 0; j < n; ++j) {
                const auto& [idx, eg] = batch[static_cast<size_t>(j)];
                if (eg.fasta.empty()) continue;
                const std::string& acc = to_score[idx];

                auto git = acc_genus.find(acc);
                if (git == acc_genus.end()) continue;

                const uint64_t gh = GcovWriter::hash_genus(git->second);
                auto calib = mrk_rd.lookup_lineage(gh);
                if (!calib.valid()) continue;

                auto contigs = parse_fasta(eg.fasta);
                if (contigs.empty()) continue;

                const bool     is_arc = (calib.header->domain == MRKR_DOMAIN_ARC);
                auto mh  = is_arc ? mrk_rd.merged_hashes_arc() : mrk_rd.merged_hashes_bac();
                auto mid = is_arc ? mrk_rd.merged_ids_arc()    : mrk_rd.merged_ids_bac();
                const BlockedBloom& bloom =
                    is_arc ? mrk_rd.merged_bloom_arc() : mrk_rd.merged_bloom_bac();
                const uint8_t n_markers = calib.header->n_markers;
                const uint8_t base_id   = is_arc ? mrk_rd.n_bac() : 0;

                std::vector<uint32_t> pool_sizes(n_markers);
                for (uint8_t mi = 0; mi < n_markers; ++mi)
                    pool_sizes[mi] = static_cast<uint32_t>(
                        mrk_rd.pool_hashes(static_cast<uint8_t>(base_id + mi)).size());

                // null_floor_u16 / threshold_u16 are calibrated for complete-genome presence
                // calls (threshold ~47% median containment). For completeness *estimation* on
                // fragments we need high sensitivity, so we keep the original raw-count gate
                // (>= minhit) and only use pool_sizes for future fractional scoring.
                // TODO: apply per-ORF null_floor correction as max(0, hits - nf * q_sz) >= minhit.

                std::vector<uint32_t> votes(n_markers, 0);
                thread_local std::vector<uint64_t> orf_mers;

                for (const auto& contig : contigs) {
                    std::vector<uint32_t> best(n_markers, 0);

                    if (is_d6) {
                        translate_6frame(contig.seq, min_seg,
                            [&](int, const uint8_t* seg, int len, int, int) {
                                orf_mers.clear();
                                extract_d6_orf_syncmers(seg, len, orf_mers);
                                if (orf_mers.empty()) return;
                                std::sort(orf_mers.begin(), orf_mers.end());
                                orf_mers.erase(std::unique(orf_mers.begin(), orf_mers.end()),
                                               orf_mers.end());
                                std::vector<uint32_t> local(n_markers, 0);
                                const uint64_t* mhp = mh.data(), *mhe = mhp + mh.size();
                                const uint8_t*  mip = mid.data();
                                for (uint64_t h : orf_mers) {
                                    if (!bloom.might_contain(h)) continue;
                                    auto it2 = std::lower_bound(mhp, mhe, h);
                                    if (it2 != mhe && *it2 == h) {
                                        const uint8_t id = mip[it2 - mhp];
                                        if (id < n_markers) local[id]++;
                                    }
                                }
                                for (uint8_t mi = 0; mi < n_markers; ++mi)
                                    if (local[mi] > best[mi]) best[mi] = local[mi];
                            });
                    } else {
                        orf_mers.clear();
                        extract_metamers_dna_into(contig.seq, METAMER_K, min_seg, frmax, orf_mers);
                        if (!orf_mers.empty()) {
                            std::sort(orf_mers.begin(), orf_mers.end());
                            orf_mers.erase(std::unique(orf_mers.begin(), orf_mers.end()),
                                           orf_mers.end());
                            std::vector<uint32_t> hits(n_markers, 0);
                            const uint64_t* mhp = mh.data(), *mhe = mhp + mh.size();
                            const uint8_t*  mip = mid.data();
                            for (uint64_t h : orf_mers) {
                                if (!bloom.might_contain(h)) continue;
                                auto it2 = std::lower_bound(mhp, mhe, h);
                                if (it2 != mhe && *it2 == h) {
                                    const uint8_t id = mip[it2 - mhp];
                                    if (id < n_markers) hits[id]++;
                                }
                            }
                            for (uint8_t mi = 0; mi < n_markers; ++mi)
                                if (hits[mi] > best[mi]) best[mi] = hits[mi];
                        }
                    }

                    for (uint8_t mi = 0; mi < n_markers; ++mi)
                        if (best[mi] >= minhit) votes[mi]++;
                }

                // Aggregate + per-SCG bitmask
                int n_present = 0, n_expected = 0, n_redundant = 0;
                const size_t nbytes = (static_cast<size_t>(n_markers) + 7) / 8;
                std::vector<uint8_t> bits(nbytes, 0);
                for (uint8_t mi = 0; mi < n_markers; ++mi) {
                    if (!calib.marker_expected(mi)) continue;
                    ++n_expected;
                    if (votes[mi] >= 1) {
                        ++n_present;
                        bits[mi / 8] |= uint8_t(1u << (mi % 8));
                    }
                    if (votes[mi] >= 2) ++n_redundant;
                }
                if (n_expected == 0) continue;

                const float ne   = static_cast<float>(n_expected);
                const float comp = static_cast<float>(n_present)  / ne;
                const float redun = static_cast<float>(n_redundant) / ne;

                std::lock_guard<std::mutex> lk(quality_mtx);
                auto& q = quality[acc];
                q.marker_completeness  = comp;
                q.marker_redundancy    = redun;
                q.marker_n_present     = n_present;
                q.marker_n_expected    = n_expected;
                q.marker_present_bits  = std::move(bits);
            }
        });

    spdlog::info("check pass-Marker: done");
}

} // namespace genopack::check
