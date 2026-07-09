#include "pass_aamer.hpp"
#include "pass_b.hpp"
#include <genopack/aamer.hpp>
#include <algorithm>
#include <cmath>
#include <omp.h>
#include <spdlog/spdlog.h>
#include <unordered_set>
#include <vector>

namespace genopack::check {

void run_pass_aamer(ICheckReader& pack,
                    const PassAResult& pass_a,
                    std::unordered_map<std::string, GenomeQuality>& quality,
                    int threads) {
    if (!pack.has_pcore()) return;

    std::vector<std::string> candidates;
    for (const auto& acc : pass_a.accessions) {
        auto it = quality.find(acc);
        if (it != quality.end() && std::isnan(it->second.completeness_cluster_relative))
            candidates.push_back(acc);
    }
    if (candidates.empty()) return;

    std::unordered_map<std::string, std::string> acc_to_genus, acc_to_family;
    acc_to_genus.reserve(candidates.size());
    acc_to_family.reserve(candidates.size());
    for (const auto& [genus, members] : pass_a.genus_members)
        for (const auto& acc : members)
            acc_to_genus[acc] = genus;
    for (const auto& [family, members] : pass_a.family_members)
        for (const auto& acc : members)
            acc_to_family[acc] = family;

    std::unordered_set<std::string> unique_genera, unique_families;
    for (const auto& acc : candidates) {
        if (auto it = acc_to_genus.find(acc);  it != acc_to_genus.end())  unique_genera.insert(it->second);
        if (auto it = acc_to_family.find(acc); it != acc_to_family.end()) unique_families.insert(it->second);
    }

    std::vector<std::string> gvec(unique_genera.begin(), unique_genera.end());
    std::vector<std::string> fvec(unique_families.begin(), unique_families.end());
    const int load_par = std::min(threads, 8);
    std::vector<std::vector<uint64_t>> gcore_buf(gvec.size()), fcore_buf(fvec.size());
    #pragma omp parallel for schedule(dynamic, 1) num_threads(load_par)
    for (size_t i = 0; i < gvec.size(); ++i) {
        PcoreView v = pack.pcore_for_genus(gvec[i]);
        if (v.valid()) v.core_into(gcore_buf[i]);
    }
    // Family cores: prefer the inline SEC_FCORE section (a CoreReader whose entries are
    // already the prevalence core — count >= ceil(theta*n)). Fall back to the PCORE
    // family tier (added post-hoc by `genopack pcore`) for packs built without inline
    // FCORE. Both yield a sorted core-aamer vector for coverage().
    const bool use_fcore_section = pack.has_fcore();
    #pragma omp parallel for schedule(dynamic, 1) num_threads(load_par)
    for (size_t i = 0; i < fvec.size(); ++i) {
        if (use_fcore_section) {
            CoreView v = pack.core_for_family(fvec[i]);
            if (v.valid()) fcore_buf[i].assign(v.aamers, v.aamers + v.n_aamers);
        } else {
            PcoreView v = pack.pcore_for_family(fvec[i]);
            if (v.valid()) v.core_into(fcore_buf[i]);
        }
    }
    std::unordered_map<std::string, const std::vector<uint64_t>*> genus_core_ptr, family_core_ptr;
    for (size_t i = 0; i < gvec.size(); ++i)
        if (!gcore_buf[i].empty()) genus_core_ptr[gvec[i]] = &gcore_buf[i];
    for (size_t i = 0; i < fvec.size(); ++i)
        if (!fcore_buf[i].empty()) family_core_ptr[fvec[i]] = &fcore_buf[i];

    spdlog::info("check pass-aamer: {} candidates, {} genera, {} families",
                 candidates.size(), genus_core_ptr.size(), family_core_ptr.size());

    const PcoreReader* pr = pack.pcore_reader();
    const int k_shard_readers = std::min(16, std::max(4, threads / 2));
    const int inner_threads   = std::max(1, threads / k_shard_readers);

    struct AamerResult { float genus_cov = NAN; float family_cov = NAN; };
    std::vector<AamerResult> results(candidates.size());

    auto coverage = [](const std::vector<uint64_t>* core, const std::vector<uint64_t>& qm) -> float {
        if (!core || core->empty()) return NAN;
        uint32_t inter = 0; size_t qi = 0;
        for (uint64_t a : *core) {
            while (qi < qm.size() && qm[qi] < a) ++qi;
            if (qi < qm.size() && qm[qi] == a) ++inter;
        }
        return static_cast<float>(inter) / static_cast<float>(core->size());
    };

    pack.visit_shard_batches_parallel(candidates, k_shard_readers,
        [&](ArchiveReader::ShardBatch& batch) {
            const int n = static_cast<int>(batch.size());
            #pragma omp parallel for schedule(dynamic, 1) num_threads(inner_threads)
            for (int j = 0; j < n; ++j) {
                const auto& [idx, eg] = batch[static_cast<size_t>(j)];
                if (eg.fasta.empty()) continue;
                auto contigs = parse_fasta(eg.fasta);
                if (contigs.empty()) continue;
                const std::string& acc = candidates[idx];
                const std::vector<uint64_t>* gcore = nullptr;
                const std::vector<uint64_t>* fcore = nullptr;
                if (auto it = acc_to_genus.find(acc); it != acc_to_genus.end())
                    if (auto pit = genus_core_ptr.find(it->second); pit != genus_core_ptr.end())
                        gcore = pit->second;
                if (auto it = acc_to_family.find(acc); it != acc_to_family.end())
                    if (auto pit = family_core_ptr.find(it->second); pit != family_core_ptr.end())
                        fcore = pit->second;
                if (!gcore && !fcore) continue;
                thread_local std::vector<uint64_t> core_qmers;
                core_qmers.clear();
                for (const auto& c : contigs)
                    extract_aamers_dna_into(c.seq, pr->k(), pr->min_seg_aa(), pr->frac_max_hash(), core_qmers);
                std::sort(core_qmers.begin(), core_qmers.end());
                core_qmers.erase(std::unique(core_qmers.begin(), core_qmers.end()), core_qmers.end());
                results[idx] = {coverage(gcore, core_qmers), coverage(fcore, core_qmers)};
            }
        });

    size_t scored = 0;
    for (size_t i = 0; i < candidates.size(); ++i) {
        auto it = quality.find(candidates[i]);
        if (it == quality.end()) continue;
        it->second.completeness_aamer_core        = results[i].genus_cov;
        it->second.completeness_aamer_family_core = results[i].family_cov;
        if (!std::isnan(results[i].genus_cov) || !std::isnan(results[i].family_cov)) ++scored;
    }
    spdlog::info("check pass-aamer: complete — {}/{} genomes scored", scored, candidates.size());
}

} // namespace genopack::check
