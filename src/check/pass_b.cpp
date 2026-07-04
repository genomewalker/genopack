#include "pass_b.hpp"
#include "tnf.hpp"
#include "foreign_contam.hpp"
#include <genopack/markers.hpp>
#include <genopack/aamer.hpp>
#include <genopack/score_bin.hpp>
#include <genopack/gcov.hpp>
#include <genopack/qual.hpp>
#include <genopack/util.hpp>
#include <Eigen/Dense>
#include <algorithm>
#include <atomic>
#include <chrono>
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
                const std::unordered_map<uint64_t, QualRecord>* baseline_cache,
                genopack::GlobalMultiplicityIndex* preloaded_gmi,
                std::future<void>* gmi_future)
{
    std::atomic<int64_t> s_extract_ns{0}, s_score_ns{0}, s_genome_ns{0}, s_mutex_ns{0};
    std::atomic<int64_t> s_parse_ns{0}, s_completeness_ns{0}, s_contig_ns{0},
                         s_gcov_ns{0}, s_marker_ns{0}, s_mixture_ns{0};
    std::atomic<int64_t> s_gstx_ns{0}, s_fmh_ns{0}, s_lock_ns{0};
    using Clock = std::chrono::steady_clock;
    auto rss_kb = []() -> long {
        long v = 0;
        if (FILE* f = fopen("/proc/self/status", "r")) {
            char line[128];
            while (fgets(line, sizeof(line), f))
                if (sscanf(line, "VmRSS: %ld", &v) == 1) break;
            fclose(f);
        }
        return v;
    };
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
        if (cfg.scan_all || by_leakage || by_tnf || by_completeness) {
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
                (fmhr_rd && r.fmh_minority_u8 == 0 && r.fmh_minority_u16 == 0 && !std::isnan(r.chromosome_skew_closure)) ||
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
            if (fmhr_rd)
                q.fmh_minority_fraction = r.fmh_minority_u16 > 0
                    ? r.fmh_minority_u16 / 65535.0f
                    : r.fmh_minority_u8  / 255.0f;
            if (r.marker_completeness_u8 > 0)
                q.marker_completeness = (r.marker_completeness_u8 - 1) / 254.0f;
            q.contamination_contig_outlier  = r.contig_outlier_u16 > 0
                ? r.contig_outlier_u16 / 65535.0f
                : r.contig_outlier_u8  / 255.0f;
            q.contamination_spe             = r.spe_outlier_u8     / 255.0f;
            q.contamination_sibling_outlier = r.sibling_outlier_u8 / 255.0f;
            q.contamination_rho_outlier     = r.rho_outlier_u8     / 255.0f;
            q.contamination_cross_genus     = r.cross_genus_u8     / 255.0f;
            if (r.sketch_fill_u8 > 0)
                q.completeness_sketch_fill = r.sketch_fill_u8 / 200.0f;
            q.contamination_mixture = r.contamination_mixture;
            q.mixture_sources       = r.mixture_sources;
            q.n_mix_windows         = r.n_mix_windows;
            if (!std::isnan(r.completeness_aamer_core))
                q.completeness_aamer_core        = r.completeness_aamer_core;
            if (!std::isnan(r.completeness_aamer_family_core))
                q.completeness_aamer_family_core = r.completeness_aamer_family_core;
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

    // ── Foreign-aamer containment via Global Multiplicity Index ─────────────
    // Build a global aamer→genus_count map once, then score with only the host
    // genus union per genome. Avoids the 523GB per-genus refset blowup: peak
    // memory is GMI (~5–15 GB) + host/family unions for flagged genera (~20 GB).
    int      foreign_k = 0, foreign_min_seg = 0;
    uint64_t foreign_frac_max = UINT64_MAX;
    const bool use_pcore = pack.has_pcore();
    if (use_pcore) {
        const PcoreReader* pr = pack.pcore_reader();
        foreign_k = pr->k(); foreign_min_seg = pr->min_seg_aa(); foreign_frac_max = pr->frac_max_hash();
    } else if (pack.has_core()) {
        const CoreReader* cr = pack.core_reader();
        foreign_k = cr->k(); foreign_min_seg = cr->min_seg_aa(); foreign_frac_max = cr->frac_max_hash();
    }

    // Materialize helpers — defined outside the CORE/PCORE guard so both the
    // per-genus refset build loop and the scan callback can capture them via [&].
    // Return an empty OwnedSpan (valid()==false) when no section is present.
    auto materialize_genus = [&](const std::string& g) -> OwnedSpan {
        OwnedSpan o;
        if (use_pcore) { PcoreView v = pack.pcore_for_genus(g); if (v.valid()) v.materialize(o.aamers, o.prevf); }
        else if (pack.has_core()) { CoreView v = pack.core_for_genus(g); if (v.valid()) o.aamers.assign(v.aamers, v.aamers + v.n_aamers); }
        return o;
    };
    auto materialize_family = [&](const std::string& f) -> OwnedSpan {
        OwnedSpan o;
        if (use_pcore) { PcoreView v = pack.pcore_for_family(f); if (v.valid()) v.materialize(o.aamers, o.prevf); }
        else if (pack.has_core()) { CoreView v = pack.core_for_family(f); if (v.valid()) o.aamers.assign(v.aamers, v.aamers + v.n_aamers); }
        return o;
    };
    const uint32_t n_entries = use_pcore ? pack.pcore_reader()->n_entries()
                             : (pack.has_core() ? pack.core_reader()->n_genera() : 0u);
    auto entry_hash_at = [&](uint32_t i) -> uint64_t {
        if (use_pcore) return pack.pcore_reader()->key_hash_at(i);
        if (pack.has_core()) return pack.core_reader()->genus_hash_at(i);
        return 0u;
    };

    // GMI: the multiplicity index that drives foreign-aamer scoring.
    // Runtime GMI build (~10-30 min) from PCORE/CORE union when present.
    const genopack::GamiView gami_view{};   // legacy Bloom path (unused)
    GlobalMultiplicityIndex local_gmi;
    GlobalMultiplicityIndex& gmi = preloaded_gmi ? *preloaded_gmi : local_gmi;
    // Preloaded GMI (started concurrently with pass-A) is joined after host unions below.
    // Fallback: load synchronously here only when no preloaded GMI is provided.
    if (!preloaded_gmi && pack.has_gami_v2()) {
        pack.load_gami_into(gmi);
        spdlog::info("check pass-B: GMI loaded from SEC_GAMI v2 ({} M entries, {:.1f} GB)",
                     gmi.count() / 1'000'000, gmi.bytes() / 1e9);
        spdlog::info("pass-B RSS after GMI load: {:.1f} GB", rss_kb() / 1e6);
    }
    std::unordered_map<std::string, OwnedSpan>           genus_host_union;
    std::unordered_map<std::string, OwnedSpan>           family_union_map;
    std::unordered_map<std::string, IndexedSpan>         genus_host_index;
    std::unordered_map<std::string, IndexedSpan>         family_index_map;
    std::unordered_map<std::string, std::vector<uint64_t>> genus_core_cache;
    std::unordered_map<std::string, std::vector<uint64_t>> family_core_cache;
    std::mutex eviction_mtx;
    std::unordered_map<std::string, std::atomic<int>> genus_refcount, family_refcount;
    const uint8_t foreign_K = [&]() -> uint8_t {
        if (const char* e = std::getenv("GENOPACK_FOREIGN_K")) {
            int v = std::atoi(e); if (v > 0 && v < 256) return static_cast<uint8_t>(v);
        }
        const uint32_t k = std::max(2u, n_entries / 500u);
        return static_cast<uint8_t>(k < 255u ? k : 255u);
    }();
    if (!preloaded_gmi && gmi.empty() && (use_pcore || pack.has_core()) && !flagged_genus.empty()) {
        const int gmi_par = std::min(threads, 8);
        // Exact prescan eliminates all rehash transients.
        gmi.reserve(use_pcore ? pack.pcore_reader()->total_union_aamers() : 64'000'000ULL);
        // Batched parallel decode → sequential insert avoids atomic/mutex on GMI.
        std::vector<OwnedSpan> buf(gmi_par);
        for (uint32_t bs = 0; bs < n_entries; bs += static_cast<uint32_t>(gmi_par)) {
            const uint32_t be = std::min(bs + static_cast<uint32_t>(gmi_par), n_entries);
            const int bsz = static_cast<int>(be - bs);
            #pragma omp parallel for schedule(static) num_threads(bsz)
            for (int bi = 0; bi < bsz; ++bi) {
                buf[bi] = OwnedSpan{};
                const uint64_t gh = entry_hash_at(bs + static_cast<uint32_t>(bi));
                if (use_pcore) { PcoreView v = pack.pcore_reader()->lookup(gh); if (v.valid()) v.materialize(buf[bi].aamers, buf[bi].prevf); }
                else if (pack.has_core()) { CoreView v = pack.core_reader()->lookup(gh); if (v.valid()) buf[bi].aamers.assign(v.aamers, v.aamers + v.n_aamers); }
            }
            for (int bi = 0; bi < bsz; ++bi) {
                for (uint64_t h : buf[bi].aamers) gmi.increment(h);
                buf[bi].aamers.clear(); buf[bi].aamers.shrink_to_fit();
                gmi.maybe_resize();
            }
        }
        spdlog::info("check pass-B: GMI built: {} M unique aamers, {} entries, {:.1f} GB",
                     gmi.count() / 1'000'000, n_entries, gmi.bytes() / 1e9);
    }
    // Host and family unions: built only when PCORE/CORE present and genera are flagged.
    // Required for host-specific aamer classification in score_contig_foreign_indexed.
    if ((use_pcore || pack.has_core()) && !flagged_genus.empty()) {
        const int gmi_par = std::min(threads, 8);
        std::unordered_set<std::string> unique_genera, unique_families;
        for (const auto& acc : to_scan) {
            auto git = flagged_genus.find(acc);  if (git != flagged_genus.end()) unique_genera.insert(git->second);
            auto fit = flagged_family.find(acc);  if (fit != flagged_family.end()) unique_families.insert(fit->second);
        }
        std::vector<std::string> gvec(unique_genera.begin(), unique_genera.end());
        std::vector<std::string> fvec(unique_families.begin(), unique_families.end());
        genus_refcount.reserve(unique_genera.size());
        family_refcount.reserve(unique_families.size());
        for (const auto& acc : to_scan) {
            if (auto it = flagged_genus.find(acc); it != flagged_genus.end()) {
                auto [rit, inserted] = genus_refcount.try_emplace(it->second, 0);
                (void)inserted;
                rit->second.fetch_add(1, std::memory_order_relaxed);
            }
            if (auto it = flagged_family.find(acc); it != flagged_family.end()) {
                auto [rit, inserted] = family_refcount.try_emplace(it->second, 0);
                (void)inserted;
                rit->second.fetch_add(1, std::memory_order_relaxed);
            }
        }
        for (const auto& g : gvec) genus_host_union[g];
        for (const auto& f : fvec) family_union_map[f];
        #pragma omp parallel for schedule(dynamic, 1) num_threads(gmi_par)
        for (size_t i = 0; i < gvec.size(); ++i) genus_host_union[gvec[i]] = materialize_genus(gvec[i]);
        #pragma omp parallel for schedule(dynamic, 1) num_threads(gmi_par)
        for (size_t i = 0; i < fvec.size(); ++i) family_union_map[fvec[i]] = materialize_family(fvec[i]);
        spdlog::info("check pass-B: host/family unions: {} genera + {} families, K={}",
                     genus_host_union.size(), family_union_map.size(), foreign_K);

        // Pre-decode PCORE core slices per unique genus/family for completeness scoring.
        // Avoids O(G_flagged × T_per_genus) PFor re-decodes in the scan loop; reduces to O(G_flagged).
        if (use_pcore) {
            std::vector<std::vector<uint64_t>> gcore_buf(gvec.size());
            std::vector<std::vector<uint64_t>> fcore_buf(fvec.size());
            #pragma omp parallel for schedule(dynamic,1) num_threads(gmi_par)
            for (size_t i = 0; i < gvec.size(); ++i) {
                PcoreView v = pack.pcore_for_genus(gvec[i]);
                if (v.valid()) v.core_into(gcore_buf[i]);
            }
            #pragma omp parallel for schedule(dynamic,1) num_threads(gmi_par)
            for (size_t i = 0; i < fvec.size(); ++i) {
                PcoreView v = pack.pcore_for_family(fvec[i]);
                if (v.valid()) v.core_into(fcore_buf[i]);
            }
            genus_core_cache.reserve(gvec.size());
            family_core_cache.reserve(fvec.size());
            for (size_t i = 0; i < gvec.size(); ++i) genus_core_cache[gvec[i]] = std::move(gcore_buf[i]);
            for (size_t i = 0; i < fvec.size(); ++i) family_core_cache[fvec[i]] = std::move(fcore_buf[i]);
        }
        // Build 16-bit prefix indexes — replaces O(10M) merge-scan with O(~153) binary search.
        genus_host_index.reserve(genus_host_union.size());
        for (const auto& [g, os] : genus_host_union)
            if (os.valid()) genus_host_index.emplace(g, build_index(os));
        family_index_map.reserve(family_union_map.size());
        for (const auto& [f, os] : family_union_map)
            if (os.valid()) family_index_map.emplace(f, build_index(os));
        spdlog::info("pass-B RSS after host unions: {:.1f} GB", rss_kb() / 1e6);
        spdlog::info("check pass-B: host index: {} genera, {} families built",
                     genus_host_index.size(), family_index_map.size());
    }
    // Join preloaded GMI (started concurrently with pass-A, ~37s).
    // Host unions ran concurrently; by now GMI is almost certainly done.
    if (gmi_future && gmi_future->valid()) {
        gmi_future->get();
        spdlog::info("check pass-B: GMI loaded from SEC_GAMI v2 ({} M entries, {:.1f} GB)",
                     gmi.count() / 1'000'000, gmi.bytes() / 1e9);
        spdlog::info("pass-B RSS after GMI load: {:.1f} GB", rss_kb() / 1e6);
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
    // TODO(D3): load a mobile-element marker index when GENOPACK_MOBILE_INDEX is set.
    // const char* mob_idx_path = std::getenv("GENOPACK_MOBILE_INDEX");
    // if (mob_idx_path) { ... open MobileMarkerIndex ... }
    if (mrk_rd && mrk_rd->is_open()) {
        mrk_rd->build_merged_pool();
        spdlog::info("check pass-B: merged pool: {} M bac + {} M arc hashes",
                     mrk_rd->merged_hashes_bac().size() / 1'000'000,
                     mrk_rd->merged_hashes_arc().size() / 1'000'000);
    }

    const uint64_t mrk_frac_max = mrk_rd ? mrk_rd->frac_max_hash() : UINT64_MAX;

    // Pre-build per-genus GCOV candidate lists (family siblings + ≤50 background).
    // Replaces the O(n_all_genera) scan_early per contig with O(50) direct iteration —
    // critical for 9.2M+ archives where n_genera can be in the thousands.
    std::unordered_map<std::string, std::vector<const GcovEntry*>> genus_gcov_cands;
    {
        const GcovReader* gr = pack.gcov_reader();
        if (gr && gr->is_open() && gr->n_genera() > 0) {
            // Collect all valid entries once for background sampling
            std::vector<const GcovEntry*> all_valid;
            all_valid.reserve(gr->n_genera());
            gr->scan([&](const GcovEntry& fe) { if (fe.flags & GCOV_FLAG_VALID) all_valid.push_back(&fe); });

            std::unordered_map<std::string, std::string> genus_to_family;
            for (const auto& acc : to_scan) {
                auto git = flagged_genus.find(acc); if (git == flagged_genus.end()) continue;
                auto fit = flagged_family.find(acc);
                genus_to_family.emplace(git->second, fit != flagged_family.end() ? fit->second : std::string{});
            }

            for (const auto& [genus, family] : genus_to_family) {
                auto& cands = genus_gcov_cands[genus];
                std::unordered_set<uint64_t> used;
                const uint64_t self_h = GcovWriter::hash_genus(genus);
                used.insert(self_h);
                // Family siblings
                auto fmit = pass_a.family_to_genera.find(family);
                if (fmit != pass_a.family_to_genera.end())
                    for (const auto& sib : fmit->second) {
                        const uint64_t sh = GcovWriter::hash_genus(sib);
                        if (!used.insert(sh).second) continue;
                        const GcovEntry* ge = gr->lookup(sh);
                        if (ge && (ge->flags & GCOV_FLAG_VALID)) cands.push_back(ge);
                    }
                // Background: deterministic sample from all_valid
                if (!all_valid.empty()) {
                    const size_t n = all_valid.size();
                    const size_t start = self_h % n;
                    uint32_t added = 0;
                    for (size_t s = 0; s < n && added < 50; ++s) {
                        const GcovEntry* fe = all_valid[(start + s) % n];
                        if (used.insert(fe->genus_hash).second) { cands.push_back(fe); ++added; }
                    }
                }
            }
            size_t total_cands = 0;
            for (auto& kv : genus_gcov_cands) total_cands += kv.second.size();
            if (!genus_gcov_cands.empty())
                spdlog::info("check pass-B: bounded GCOV scan: {} genera, avg {:.0f} candidates (vs {} total)",
                             genus_gcov_cands.size(),
                             static_cast<double>(total_cands) / genus_gcov_cands.size(),
                             gr->n_genera());
        }
    }

    // N parallel reader threads each own an independent fd and sequential shard band.
    // Callback is invoked concurrently; quality_mtx guards shared state.
    const int k_shard_readers = std::min(16, std::max(4, threads / 2));
    const int inner_threads = std::max(1, threads / k_shard_readers);
    std::mutex quality_mtx;
    size_t n_done = 0;
    const size_t n_flagged = to_scan.size();
    auto release_host_refs = [&](const std::string& acc) {
        if (auto git = flagged_genus.find(acc); git != flagged_genus.end()) {
            const std::string& g = git->second;
            if (auto rit = genus_refcount.find(g); rit != genus_refcount.end()) {
                if (rit->second.fetch_sub(1, std::memory_order_acq_rel) == 1) {
                    std::lock_guard<std::mutex> lk(eviction_mtx);
                    genus_host_index.erase(g);
                    genus_host_union.erase(g);
                }
            }
        }
        if (auto fit = flagged_family.find(acc); fit != flagged_family.end()) {
            const std::string& f = fit->second;
            if (auto rit = family_refcount.find(f); rit != family_refcount.end()) {
                if (rit->second.fetch_sub(1, std::memory_order_acq_rel) == 1) {
                    std::lock_guard<std::mutex> lk(eviction_mtx);
                    family_index_map.erase(f);
                    family_union_map.erase(f);
                }
            }
        }
    };

    pack.visit_shard_batches_parallel(to_scan, k_shard_readers,
        [&](ArchiveReader::ShardBatch& batch) {
            const int n = static_cast<int>(batch.size());
            #pragma omp parallel for schedule(dynamic, 1) num_threads(inner_threads)
            for (int j = 0; j < n; ++j) {
                const auto& [idx, eg] = batch[static_cast<size_t>(j)];
                const std::string& acc = to_scan[idx];
                if (eg.fasta.empty()) {
                    release_host_refs(acc);
                    continue;
                }

                auto contigs = parse_fasta(eg.fasta);
                if (contigs.empty()) {
                    release_host_refs(acc);
                    continue;
                }

                const auto t_genome_start = Clock::now();
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

                // Indexed host/family spans for foreign-aamer scoring (read-only, zero-copy).
                // Use pointers into the map entries — valid for the genome's lifetime since
                // release_host_refs() only evicts AFTER processing completes (line ~1178).
                const IndexedSpan* host_ix_ptr = nullptr;
                const IndexedSpan* fam_ix_ptr  = nullptr;
                if (git != flagged_genus.end()) {
                    auto t_m = Clock::now();
                    {
                    std::lock_guard<std::mutex> lk(eviction_mtx);
                    if (!genus_host_index.empty()) {
                        auto hit = genus_host_index.find(git->second);
                        if (hit != genus_host_index.end()) host_ix_ptr = &hit->second;
                        auto famit = flagged_family.find(acc);
                        if (famit != flagged_family.end()) {
                            auto fit = family_index_map.find(famit->second);
                            if (fit != family_index_map.end()) fam_ix_ptr = &fit->second;
                        }
                    }
                    }
                    s_mutex_ns.fetch_add((Clock::now()-t_m).count(), std::memory_order_relaxed);
                }

                const auto t_c1 = Clock::now(); s_parse_ns.fetch_add((t_c1-t_genome_start).count(), std::memory_order_relaxed);
                // ── Intrinsic completeness: genus prevalence-core coverage ──────
                // coverage = |genome aamers ∩ genus core| / |core|, extracting the
                // query genome's aamers with the exact params the CORE section was
                // built with (k / min_seg_aa / FracMinHash max). Validated to tie
                // CheckM2 (RMSE ~5.1%) and ~280x faster than cluster-relative.
                // Genus-core coverage is the primary signal; family-core coverage is
                // the fallback for novel/sparse genera (SEC_FCORE). Both are recorded
                // as DISTINCT columns — the reporting profile prefers genus, falls back
                // to family, never silently substitutes. The query aamer set is the
                // same for both (FCORE built with CORE's k/min_seg_aa/frac_max), so it
                // is extracted once and intersected against each core.
                float completeness_aamer_core        = NAN;  // genus
                float completeness_aamer_family_core = NAN;  // family fallback
                {
                    auto famit = flagged_family.find(acc);
                    const std::string* genus_s  = (git   != flagged_genus.end())  ? &git->second   : nullptr;
                    const std::string* family_s = (famit != flagged_family.end()) ? &famit->second : nullptr;
                    if (pack.has_pcore()) {
                        auto gcit = genus_s  ? genus_core_cache.find(*genus_s)   : genus_core_cache.end();
                        auto fcit = family_s ? family_core_cache.find(*family_s) : family_core_cache.end();
                        const std::vector<uint64_t>* gcore = (gcit != genus_core_cache.end()  && !gcit->second.empty()) ? &gcit->second  : nullptr;
                        const std::vector<uint64_t>* fcore = (fcit != family_core_cache.end() && !fcit->second.empty()) ? &fcit->second : nullptr;
                        if (gcore || fcore) {
                            const PcoreReader* pr = pack.pcore_reader();
                            thread_local std::vector<uint64_t> core_qmers;
                            core_qmers.clear();
                            for (const auto& c : contigs)
                                extract_aamers_dna_into(c.seq, pr->k(), pr->min_seg_aa(),
                                                        pr->frac_max_hash(), core_qmers);
                            std::sort(core_qmers.begin(), core_qmers.end());
                            core_qmers.erase(std::unique(core_qmers.begin(), core_qmers.end()),
                                             core_qmers.end());
                            auto cached_coverage = [](const std::vector<uint64_t>* core,
                                                      const std::vector<uint64_t>& qm) -> float {
                                if (!core || core->empty()) return NAN;
                                uint32_t inter = 0; size_t qi = 0;
                                for (uint64_t a : *core) {
                                    while (qi < qm.size() && qm[qi] < a) ++qi;
                                    if (qi < qm.size() && qm[qi] == a) ++inter;
                                }
                                return static_cast<float>(inter) / static_cast<float>(core->size());
                            };
                            completeness_aamer_core        = cached_coverage(gcore, core_qmers);
                            completeness_aamer_family_core = cached_coverage(fcore, core_qmers);
                        }
                    } else {
                        CoreView gcore = (pack.has_core() && genus_s)  ? pack.core_for_genus(*genus_s)   : CoreView{};
                        CoreView fcore = (pack.has_fcore() && family_s) ? pack.core_for_family(*family_s) : CoreView{};
                        if (gcore.valid() || fcore.valid()) {
                            const CoreReader* cr = pack.has_core() ? pack.core_reader() : pack.fcore_reader();
                            thread_local std::vector<uint64_t> core_qmers;
                            core_qmers.clear();
                            for (const auto& c : contigs)
                                extract_aamers_dna_into(c.seq, cr->k(), cr->min_seg_aa(),
                                                        cr->frac_max_hash(), core_qmers);
                            std::sort(core_qmers.begin(), core_qmers.end());
                            core_qmers.erase(std::unique(core_qmers.begin(), core_qmers.end()),
                                             core_qmers.end());
                            auto coverage = [&](const CoreView& cv) -> float {
                                if (!cv.valid()) return NAN;
                                uint32_t inter = 0;
                                size_t qi = 0, kj = 0;
                                while (qi < core_qmers.size() && kj < cv.n_aamers) {
                                    if      (core_qmers[qi] < cv.aamers[kj]) ++qi;
                                    else if (core_qmers[qi] > cv.aamers[kj]) ++kj;
                                    else { ++inter; ++qi; ++kj; }
                                }
                                return static_cast<float>(inter) / static_cast<float>(cv.n_aamers);
                            };
                            completeness_aamer_core        = coverage(gcore);
                            completeness_aamer_family_core = coverage(fcore);
                        }
                    }
                }

                const GcovReader* fcov_rd    = pack.fcov_reader();
                const GcovEntry*  fcov_entry = nullptr;
                auto fit = flagged_family.find(acc);
                if (fit != flagged_family.end() && fcov_rd) {
                    const GcovEntry* fe = pack.fcov_for_family(fit->second);
                    if (fe && (fe->flags & GCOV_FLAG_VALID)) fcov_entry = fe;
                }

                const auto t_c2 = Clock::now(); s_completeness_ns.fetch_add((t_c2-t_c1).count(), std::memory_order_relaxed);
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
                                    float dist, spe;
                                    gcov_mahalanobis_spe(*gcov_entry, xmu, &dist, &spe);
                                    const float pct     = gcov_rd->percentile(*gcov_entry, dist);
                                    const float spe_pct = gcov_rd->spe_percentile(*gcov_entry, spe);
                                    cf.gcov_t2_pct  = pct;       // persisted into QCONTIG
                                    cf.gcov_spe_pct = spe_pct;   // the cross-genus / foreign-contig signal
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
                                            float f_dist, f_spe;
                                            gcov_mahalanobis_spe(*fcov_entry, xmu_fam, &f_dist, &f_spe);
                                            const float f_pct  = fcov_rd->percentile(*fcov_entry, f_dist);
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

                        // ── Foreign-aamer containment: the SMALL-contig channel ──
                        // Composition-independent; runs for every contig >=1kb (where TNF
                        // abstains). Classifies the contig's 6-frame aamers as host- vs
                        // foreign-specific against the precomputed reference set.
                        // Skip extraction for large (>=5kb) contigs that are TNF inliers:
                        // the composition signal already confirms host origin.
                        if (host_ix_ptr && host_ix_ptr->valid() && contig.seq.size() >= 1000) {
                            const bool skip_aamer = contig.seq.size() >= 5000
                                && !std::isnan(cf.tnf_score)
                                && cf.tnf_score < cfg.contig_tnf_threshold;
                            thread_local std::vector<uint64_t> prot_qm;
                            if (!skip_aamer) {
                            prot_qm.clear();
                            { auto t0 = Clock::now();
                            extract_aamers_dna_into(contig.seq, foreign_k, foreign_min_seg,
                                                    foreign_frac_max, prot_qm);
                            std::sort(prot_qm.begin(), prot_qm.end());
                            prot_qm.erase(std::unique(prot_qm.begin(), prot_qm.end()), prot_qm.end());
                            s_extract_ns.fetch_add((Clock::now()-t0).count(), std::memory_order_relaxed); }
                            static const IndexedSpan empty_ix;
                            const IndexedSpan& fam_ref = fam_ix_ptr ? *fam_ix_ptr : empty_ix;
                            { auto t0 = Clock::now();
                            ForeignScore fs = score_contig_foreign_indexed(
                                *host_ix_ptr, fam_ref, gmi, gami_view, foreign_K, prot_qm);
                            cf.prot_classifiable     = static_cast<uint16_t>(std::min<uint32_t>(fs.classifiable, 65535));
                            cf.prot_host_specific    = static_cast<uint16_t>(std::min<uint32_t>(fs.host_specific, 65535));
                            cf.prot_foreign_specific = static_cast<uint16_t>(std::min<uint32_t>(fs.foreign_specific, 65535));
                            cf.prot_family_hits      = static_cast<uint16_t>(std::min<uint32_t>(fs.family_hits, 65535));
                            cf.prot_best_genus       = fs.best_genus_tag;
                            cf.prot_foreign_frac     = fs.foreign_fraction;
                            cf.prot_loglr            = fs.score;
                            uint8_t pf = 0;
                            // Abstain floor = 16 classifiable aamers (validated: 5kb contigs
                            // carry ~25 against the prevalence core; below 16 the fraction is
                            // too noisy). 1-2kb contigs fall below this against CORE — they need
                            // a denser per-genus proteome reference (a v2 reference section).
                            if (fs.classifiable < 16) pf |= PROT_ABSTAIN_LOW_N;
                            if (fs.foreign_specific >= 12 && fs.family_hits >= fs.foreign_specific)
                                pf |= PROT_MOBILE_NATIVE;
                            // TODO(D3): score against mobile-element index when loaded;
                            // set PROT_MOBILE_ELEMENT flag for contigs matching a mobile MGE.
                            cf.prot_flags = pf;
                            s_score_ns.fetch_add((Clock::now()-t0).count(), std::memory_order_relaxed); }
                            } // if (!skip_aamer)
                        }

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

                    // score_entry: evaluate one foreign genus GCOV entry against all pending
                    // qual_contigs using the progressive admissible-bound early-exit.
                    // Returns false when all contigs are accounted for (stop signal).
                    auto score_entry = [&](const GcovEntry& fe) -> bool {
                        if (&fe == gcov_entry || !(fe.flags & GCOV_FLAG_VALID)) return true;
                        const float ld   = gcov_log_det(fe);
                        const float vmax = std::max(fe.eigenvalues[0], fe.sigma2_resid);
                        for (auto& qc : qual_contigs) {
                            if (qc.found) continue;
                            float xf[136], l2 = 0.f;
                            for (int d = 0; d < 136; ++d) {
                                const float v = qc.tnf[d] - fe.mu[d];
                                xf[d] = v; l2 += v * v;
                            }
                            if (-0.5f * (l2 / vmax + ld) <= qc.tau) continue;
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
                                if (-0.5f * (d2 + std::max(0.f, l2 - sum_p2) / vsub + ld) <= qc.tau) {
                                    pruned = true; break;
                                }
                            }
                            if (pruned) continue;
                            qc.found = true;
                            cross_genus_bp         += qc.length;
                            ctg_cross_genus[qc.ci]  = true;
                            if (qc.family_inlier) sibling_out_bp += qc.length;
                            if (--n_unfound == 0) return false;
                        }
                        return true;
                    };

                    // Use pre-built candidate list (O(50)) if present; fall back to
                    // scan_early over all genera for archives without candidate data.
                    const std::vector<const GcovEntry*>* gcov_cands = nullptr;
                    if (git != flagged_genus.end()) {
                        auto cit = genus_gcov_cands.find(git->second);
                        if (cit != genus_gcov_cands.end()) gcov_cands = &cit->second;
                    }
                    if (gcov_cands) {
                        for (const GcovEntry* fe : *gcov_cands)
                            if (!score_entry(*fe)) break;
                    } else {
                        gcov_rd->scan_early(score_entry);
                    }
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
                    static thread_local std::vector<std::pair<float,uint32_t>> scored;
                    scored.clear();
                    for (const auto& cf : flags) {
                        if (!std::isnan(cf.tnf_score) && cf.contig_length >= 5000)
                            scored.emplace_back(cf.tnf_score, cf.contig_length);
                    }
                    if (scored.size() >= 3) {
                        static thread_local std::vector<float> sv;
                        sv.clear();
                        sv.reserve(scored.size());
                        for (const auto& [s, l] : scored) sv.push_back(s);
                        std::sort(sv.begin(), sv.end());
                        const float null_med = sv[sv.size() / 2];
                        static thread_local std::vector<float> dv;
                        dv.clear();
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

                const auto t_c3 = Clock::now(); s_contig_ns.fetch_add((t_c3-t_c2).count(), std::memory_order_relaxed);
                const auto t_c4 = Clock::now(); s_gcov_ns.fetch_add((t_c4-t_c3).count(), std::memory_order_relaxed);
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
                        const int min_seg = is_d6 ? 50 : AAMER_K;
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
                                extract_aamers_dna_into(contig.seq, AAMER_K, min_seg,
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
                            // TODO(D4): adjust denominator when mobile-element marker sets are
                            // loaded — subtract mobile-element markers from ne before dividing
                            // so completeness is not penalised for absent MGE-associated markers.
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

                const auto t_c5 = Clock::now(); s_marker_ns.fetch_add((t_c5-t_c4).count(), std::memory_order_relaxed);
                // Per-contig containment split against GSTX family-member consensus sketches.
                float contamination_contig_split = NAN;
                float contamination_self_outlier = NAN;
                float fiedler_oph_val            = NAN;
                float fiedler_tnf_bimod_val      = NAN;
                float fiedler_tnf_gap_val        = NAN;
                {
                    const auto t_gstx0 = Clock::now();
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
                    s_gstx_ns.fetch_add((Clock::now()-t_gstx0).count(), std::memory_order_relaxed);
                }

                // FMH minority fraction: zero-copy mmap'd FMHR refs.
                float fmh_minority = NAN;
                {
                    const auto t_fmh0 = Clock::now();
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
                                // Build seq list sorted descending by length; cap at 30.
                                // Contamination signal concentrates in long contigs; FMH
                                // is bp-weighted so this is unbiased for long-contig contamination.
                                constexpr size_t kMaxFmhContigs = 30;
                                std::vector<std::string_view> seqs;
                                seqs.reserve(std::min(contigs.size(), kMaxFmhContigs));
                                for (const auto& c : contigs) seqs.push_back(c.seq);
                                if (seqs.size() > kMaxFmhContigs) {
                                    std::partial_sort(seqs.begin(), seqs.begin() + kMaxFmhContigs,
                                                      seqs.end(),
                                                      [](auto& a, auto& b){ return a.size() > b.size(); });
                                    seqs.resize(kMaxFmhContigs);
                                }
                                fmh_minority = score_bin_fmh_containment(
                                    seqs, hit->second.genus_hash, refs, 21, 125);
                            }
                        }
                    }
                    s_fmh_ns.fetch_add((Clock::now()-t_fmh0).count(), std::memory_order_relaxed);
                }

                {
                    const bool skip_mix = (needs_full_set.find(acc) == needs_full_set.end());
                    AllSignals sig = compute_all_signals(
                        std::span<const ContigAccum>(accum_contigs), 5000, 16384, 5, skip_mix);
                    release_host_refs(acc);

                    const auto t_lock0 = Clock::now();
                    std::lock_guard<std::mutex> lk(quality_mtx);
                    s_lock_ns.fetch_add((Clock::now()-t_lock0).count(), std::memory_order_relaxed);
                    auto& q = quality.at(acc);
                    q.chromosome_skew_closure    = skew_score;
                    q.completeness_post_decontam = completeness_post;
                    q.completeness_aamer_core    = completeness_aamer_core;
                    q.completeness_aamer_family_core = completeness_aamer_family_core;
                    q.self_coherence             = sig.self_coherence;
                    q.chargaff_parity            = sig.chargaff_parity;
                    q.spectral_gap               = sig.spectral_gap;
                    q.scale_kink                 = sig.scale_kink;
                    q.contamination_mixture      = sig.contamination_mixture;
                    q.contamination_tnf_minor    = sig.contamination_tnf_minor;
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
                    s_mixture_ns.fetch_add((Clock::now()-t_c5).count(), std::memory_order_relaxed);
                    s_genome_ns.fetch_add((Clock::now()-t_genome_start).count(), std::memory_order_relaxed);
                    ++n_done;
                    if (n_done % 50000 == 0)
                        spdlog::info("check pass-B: {}/{} genomes done", n_done, n_flagged);
                }
            }
        });

    spdlog::info("check pass-B: complete ({} genomes)", n_flagged);
    spdlog::info("pass-B RSS after scan: {:.1f} GB", rss_kb() / 1e6);
    {
        double gn  = s_genome_ns.load()       / 1e9;
        double ps  = s_parse_ns.load()        / 1e9;
        double cm  = s_completeness_ns.load() / 1e9;
        double ct  = s_contig_ns.load()       / 1e9;
        double ex  = s_extract_ns.load()      / 1e9;
        double sc  = s_score_ns.load()        / 1e9;
        double gc  = s_gcov_ns.load()         / 1e9;
        double mk  = s_marker_ns.load()       / 1e9;
        double mx  = s_mixture_ns.load()      / 1e9;
        double mu  = s_mutex_ns.load()        / 1e9;
        auto pct = [&](double v){ return gn>0 ? 100*v/gn : 0.0; };
        spdlog::info("pass-B timing (CPU-s, {} genomes): genome={:.1f}s", n_flagged, gn);
        spdlog::info("  parse+skew={:.1f}s ({:.0f}%)  completeness={:.1f}s ({:.0f}%)  contig-loop={:.1f}s ({:.0f}%)",
                     ps, pct(ps), cm, pct(cm), ct, pct(ct));
        spdlog::info("    └─ extract={:.1f}s ({:.0f}%)  score={:.1f}s ({:.0f}%)", ex, pct(ex), sc, pct(sc));
        spdlog::info("  gcov={:.1f}s ({:.0f}%)  markers={:.1f}s ({:.0f}%)  mixture={:.1f}s ({:.0f}%)  mutex={:.1f}s ({:.0f}%)",
                     gc, pct(gc), mk, pct(mk), mx, pct(mx), mu, pct(mu));
        auto cas = genopack::get_cas_timings();
        double cas_tot = static_cast<double>(cas.scan_ns + cas.coh_ns + cas.chargaff_ns + cas.spectral_ns + cas.mix_ns);
        if (cas_tot > 0) {
            auto cpct = [&](int64_t ns){ return ns * 100.0 / cas_tot; };
            auto csec = [](int64_t ns){ return ns * 1e-9; };
            spdlog::info("  compute_all_signals: scan={:.1f}s ({:.0f}%)  coh={:.1f}s ({:.0f}%)  chargaff={:.1f}s ({:.0f}%)  spectral={:.1f}s ({:.0f}%)  mix={:.1f}s ({:.0f}%)",
                         csec(cas.scan_ns),    cpct(cas.scan_ns),
                         csec(cas.coh_ns),     cpct(cas.coh_ns),
                         csec(cas.chargaff_ns),cpct(cas.chargaff_ns),
                         csec(cas.spectral_ns),cpct(cas.spectral_ns),
                         csec(cas.mix_ns),     cpct(cas.mix_ns));
        }
        spdlog::info("  gstx={:.1f}s ({:.0f}%)  fmh={:.1f}s ({:.0f}%)  lock={:.1f}s ({:.0f}%)",
                     s_gstx_ns.load()/1e9, pct(s_gstx_ns.load()/1e9),
                     s_fmh_ns.load()/1e9,  pct(s_fmh_ns.load()/1e9),
                     s_lock_ns.load()/1e9, pct(s_lock_ns.load()/1e9));
    }

    // ── CCO self-calibrating baseline ───────────────────────────────────────
    // Use the 20th percentile of CCO values per genus as the natural-heterogeneity
    // floor; subtract from target CCO and clamp to [0,1].
    // Build acc→genus once (O(n)) to avoid the prior O(n²) nested-loop.
    {
        std::unordered_map<std::string, const std::string*> acc_to_genus;
        acc_to_genus.reserve(quality.size() * 2);
        for (const auto& [genus, members] : pass_a.genus_members)
            for (const auto& m : members)
                acc_to_genus.emplace(m, &genus);

        std::unordered_map<std::string, std::vector<float>> genus_cco;
        for (const auto& [acc, q] : quality) {
            if (std::isnan(q.contamination_contig_outlier)) continue;
            auto it = acc_to_genus.find(acc);
            if (it == acc_to_genus.end()) continue;
            genus_cco[*it->second].push_back(q.contamination_contig_outlier);
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
                auto it = acc_to_genus.find(acc);
                if (it == acc_to_genus.end()) continue;
                auto bit = genus_baseline.find(*it->second);
                if (bit == genus_baseline.end()) continue;
                q.contamination_contig_outlier_adj = std::max(
                    0.0f, q.contamination_contig_outlier - bit->second);
            }
        }
    }
}

} // namespace genopack::check
