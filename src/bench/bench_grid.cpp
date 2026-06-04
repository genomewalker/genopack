#include <genopack/archive.hpp>
#include <genopack/score_bin.hpp>
#include <genopack/skani.hpp>
#include <check/pack_iface.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <sstream>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace genopack::bench {

namespace {



static constexpr std::array<float, 6> SPIKE_FRACS  = {0.0f, 0.01f, 0.02f, 0.05f, 0.10f, 0.20f};
static constexpr uint32_t N_HOST       = 8;
static constexpr uint32_t N_SPIKE_POOL = 8;
static constexpr uint32_t MIN_LONG_BP  = 20000;

static constexpr int N_FMH = 6; // k∈{21,31} × c∈{30,125,500}

// ── Manifest ──────────────────────────────────────────────────────────────────

struct ManifestRow {
    std::string host_genus;
    std::string contam_genus;
    std::string ani_label;
    int         intra_offset = 0;  // used only when host_genus == contam_genus
};

std::vector<ManifestRow> load_manifest(const std::filesystem::path& p) {
    std::ifstream f(p);
    if (!f) throw std::runtime_error("bench-grid: cannot open manifest " + p.string());
    std::string line;
    std::getline(f, line); // skip header
    std::vector<ManifestRow> rows;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        ManifestRow r;
        std::string off;
        std::getline(ss, r.host_genus,   '\t');
        std::getline(ss, r.contam_genus, '\t');
        std::getline(ss, r.ani_label,    '\t');
        if (std::getline(ss, off, '\t') && !off.empty())
            r.intra_offset = std::stoi(off);
        if (!r.host_genus.empty() && !r.contam_genus.empty() && !r.ani_label.empty())
            rows.push_back(std::move(r));
    }
    return rows;
}

// ── Per-host context ──────────────────────────────────────────────────────────

struct HostCtx {
    const GcovEntry*              gcov_e    = nullptr;
    const GcovReader*             gcov_rd   = nullptr;
    const GcovEntry*              fcov_e    = nullptr;
    const GcovReader*             fcov_rd   = nullptr;
    const GstxReader*             gstx_rd   = nullptr;
    std::vector<const GstxEntry*> gstx_cands;
    uint64_t                      genus_hash = 0;
    std::vector<std::string>      host_accs;  // first N_HOST
};

// ── FASTA helpers ─────────────────────────────────────────────────────────────

std::vector<std::string> parse_fasta_seqs(const std::string& fasta) {
    std::vector<std::string> seqs;
    std::string cur;
    for (size_t i = 0; i < fasta.size(); ) {
        if (fasta[i] == '>') {
            if (!cur.empty()) { seqs.push_back(std::move(cur)); cur.clear(); }
            while (i < fasta.size() && fasta[i] != '\n') ++i;
            if (i < fasta.size()) ++i;
        } else {
            size_t end = fasta.find('\n', i);
            if (end == std::string::npos) end = fasta.size();
            cur.append(fasta, i, end - i);
            i = end + 1;
        }
    }
    if (!cur.empty()) seqs.push_back(std::move(cur));
    return seqs;
}

std::vector<std::string> fragment_genome(
    const std::vector<std::string>& seqs,
    std::mt19937& rng,
    uint32_t min_bp = MIN_LONG_BP)
{
    std::lognormal_distribution<double> dist(10.5, 1.0);
    std::vector<std::string> frags;
    for (const auto& s : seqs) {
        size_t pos = 0;
        while (pos < s.size()) {
            auto raw  = static_cast<size_t>(dist(rng));
            size_t fl = std::max<size_t>(min_bp, std::min(raw, s.size() - pos));
            if (fl >= min_bp) frags.push_back(s.substr(pos, fl));
            pos += fl;
        }
    }
    return frags;
}

std::unordered_map<std::string, std::vector<std::string>>
extract_genomes(check::ICheckReader& pack,
                const std::vector<std::string>& accs,
                int threads)
{
    std::unordered_map<std::string, std::vector<std::string>> result;
    result.reserve(accs.size());
    pack.visit_shard_batches_parallel(accs, std::max(1, threads / 2),
        [&](ArchiveReader::ShardBatch& batch) {
            for (const auto& [idx, eg] : batch) {
                if (eg.fasta.empty()) continue;
                result[accs[idx]] = parse_fasta_seqs(eg.fasta);
            }
        });
    return result;
}

} // namespace

// ── Public entry point ────────────────────────────────────────────────────────

int cmd_bench_grid(
    const std::filesystem::path& archive_path,
    const std::filesystem::path& manifest_path,
    const std::filesystem::path& output,
    int      threads,
    int      reps,
    uint32_t seed,
    const std::vector<float>&    completeness_fracs)
{
    auto manifest = load_manifest(manifest_path);
    if (manifest.empty()) throw std::runtime_error("bench-grid: empty manifest");
    spdlog::info("bench-grid: {} manifest rows", manifest.size());

    ArchiveReader ar;
    ar.open(archive_path);
    check::SingleArchiveCheckReader pack(ar);


    // ── Taxonomy scan ──────────────────────────────────────────────────────────
    std::unordered_map<std::string, std::vector<std::string>> genus_accs;
    std::unordered_map<std::string, std::string>              acc_family;

    pack.scan_taxonomy_with_id([&](std::string_view acc, std::string_view tax, GenomeId) {
        auto extract_tag = [&](const char* prefix) -> std::string {
            auto pos = tax.find(prefix);
            if (pos == std::string_view::npos) return {};
            pos += std::strlen(prefix);
            auto end = tax.find(';', pos);
            auto tag = std::string(tax.substr(pos, end == std::string_view::npos
                                                   ? tax.size() - pos : end - pos));
            return tag.empty() ? std::string{} : std::string(prefix) + tag;
        };
        std::string g = extract_tag("g__");
        std::string f = extract_tag("f__");
        if (!g.empty()) {
            genus_accs[g].emplace_back(acc);
            if (!f.empty()) acc_family[std::string(acc)] = f;
        }
    });

    // ── Build per-host context cache ───────────────────────────────────────────
    std::unordered_map<std::string, HostCtx> host_ctxs;
    for (const auto& row : manifest) {
        if (host_ctxs.count(row.host_genus)) continue;
        HostCtx ctx;
        ctx.gcov_rd = pack.gcov_reader();
        ctx.gcov_e  = pack.gcov_for_genus(row.host_genus);
        if (!ctx.gcov_e || !(ctx.gcov_e->flags & GCOV_FLAG_VALID)) {
            spdlog::warn("bench-grid: no valid GCOV for {} — skipping", row.host_genus);
            host_ctxs[row.host_genus] = ctx;
            continue;
        }
        ctx.fcov_rd = pack.fcov_reader();
        auto& all_accs = genus_accs[row.host_genus];
        if (!all_accs.empty()) {
            auto fit = acc_family.find(all_accs[0]);
            if (fit != acc_family.end() && ctx.fcov_rd) {
                ctx.fcov_e = pack.fcov_for_family(fit->second);
                if (ctx.fcov_e && !(ctx.fcov_e->flags & GCOV_FLAG_VALID)) ctx.fcov_e = nullptr;
            }
        }
        ctx.gstx_rd = pack.gstx_reader();
        if (ctx.gstx_rd && !all_accs.empty()) {
            ctx.genus_hash = GcovWriter::hash_genus(row.host_genus);
            std::string host_fam;
            auto ffit = acc_family.find(all_accs[0]);
            if (ffit != acc_family.end()) host_fam = ffit->second;
            if (!host_fam.empty()) {
                std::unordered_set<std::string> fam_genera;
                for (const auto& [acc, fam] : acc_family)
                    if (fam == host_fam) {
                        for (const auto& [g, gaccs] : genus_accs)
                            for (const auto& a : gaccs)
                                if (a == acc) { fam_genera.insert(g); break; }
                    }
                for (const auto& g : fam_genera) {
                    const GstxEntry* ge = pack.gstx_for_genus(g);
                    if (ge && ge->genus_hash != 0 && ge->n_k_stored > 0)
                        ctx.gstx_cands.push_back(ge);
                }
            }
        }
        if (all_accs.size() < N_HOST)
            spdlog::warn("bench-grid: only {} accs for host {} (need {})",
                         all_accs.size(), row.host_genus, N_HOST);
        for (uint32_t i = 0; i < std::min<uint32_t>(N_HOST, all_accs.size()); ++i)
            ctx.host_accs.push_back(all_accs[i]);
        spdlog::info("bench-grid: host={} fcov={} gstx_cands={}",
                     row.host_genus, ctx.fcov_e ? "yes" : "no", ctx.gstx_cands.size());
        host_ctxs[row.host_genus] = std::move(ctx);
    }

    // ── Collect all accessions to extract ─────────────────────────────────────
    std::unordered_set<std::string> to_extract_set;
    for (const auto& [g, ctx] : host_ctxs)
        for (const auto& a : ctx.host_accs) to_extract_set.insert(a);

    for (const auto& row : manifest) {
        auto& all = genus_accs[row.contam_genus];
        if (row.host_genus == row.contam_genus) {
            int off = row.intra_offset;
            for (int i = 0; i < static_cast<int>(N_SPIKE_POOL) &&
                            off + i < static_cast<int>(all.size()); ++i)
                to_extract_set.insert(all[off + i]);
        } else {
            for (size_t i = 0; i < std::min<size_t>(N_SPIKE_POOL, all.size()); ++i)
                to_extract_set.insert(all[i]);
        }
    }
    std::vector<std::string> to_extract(to_extract_set.begin(), to_extract_set.end());
    std::sort(to_extract.begin(), to_extract.end());

    spdlog::info("bench-grid: extracting {} genomes…", to_extract.size());
    auto genome_seqs = extract_genomes(pack, to_extract, threads);
    spdlog::info("bench-grid: extraction done");

    // ── Grid run ───────────────────────────────────────────────────────────────
    struct Row {
        std::string host_genus, contam_genus, ani_label;
        int         rep;
        float       spike_frac_pct, actual_frac_pct;
        float       cco, sibling, rho, contig_split, self_outlier, fiedler_oph;
        float       fiedler_tnf_bimod, fiedler_tnf_gap;
        float       fmh_minority;
        std::array<float, N_FMH> fmh;  // [k21c30, k21c125, k21c500, k31c30, k31c125, k31c500]
        float       completeness_target;
        uint32_t    scored_bp;
    };
    std::vector<Row> rows;

    const BinScoreConfig score_cfg;
    const int n_fracs   = static_cast<int>(SPIKE_FRACS.size());
    const int n_total   = static_cast<int>(manifest.size()) * reps * n_fracs
                          * static_cast<int>(completeness_fracs.size());
    int done = 0;

    for (const auto& mrow : manifest) {
        auto& ctx = host_ctxs[mrow.host_genus];
        if (!ctx.gcov_e) continue;

        auto& contam_all = genus_accs[mrow.contam_genus];
        std::vector<std::string> spike_accs;
        if (mrow.host_genus == mrow.contam_genus) {
            int off = mrow.intra_offset;
            for (int i = 0; i < static_cast<int>(N_SPIKE_POOL) &&
                            off + i < static_cast<int>(contam_all.size()); ++i)
                spike_accs.push_back(contam_all[off + i]);
        } else {
            for (size_t i = 0; i < std::min<size_t>(N_SPIKE_POOL, contam_all.size()); ++i)
                spike_accs.push_back(contam_all[i]);
        }
        if (spike_accs.empty()) {
            spdlog::warn("bench-grid: no spike accs for {} ({})", mrow.contam_genus, mrow.ani_label);
            continue;
        }

        // ── Pre-build FMH genus refs (2 sequence scans via sweep) ─────────────
        uint64_t contam_genus_hash = GcovWriter::hash_genus(mrow.contam_genus);
        bool fmh_ok = (mrow.host_genus != mrow.contam_genus);
        std::array<FmhGenusRef, 6> host_fmh_refs, contam_fmh_refs;
        if (fmh_ok) {
            std::vector<std::string_view> host_sv, spike_sv;
            for (const auto& acc : ctx.host_accs) {
                auto it = genome_seqs.find(acc);
                if (it != genome_seqs.end())
                    for (const auto& s : it->second) host_sv.push_back(s);
            }
            for (const auto& acc : spike_accs) {
                auto it = genome_seqs.find(acc);
                if (it != genome_seqs.end())
                    for (const auto& s : it->second) spike_sv.push_back(s);
            }
            fmh_ok = !host_sv.empty() && !spike_sv.empty();
            if (fmh_ok) {
                host_fmh_refs   = build_fmh_genus_refs_sweep(ctx.genus_hash,   host_sv);
                contam_fmh_refs = build_fmh_genus_refs_sweep(contam_genus_hash, spike_sv);
            }
        }

        for (float comp : completeness_fracs) {
        for (int rep = 0; rep < reps; ++rep) {
            std::mt19937 host_rng(seed + static_cast<uint32_t>(rep) * 1000);
            std::vector<std::string> host_frags_all;
            for (const auto& acc : ctx.host_accs) {
                auto it = genome_seqs.find(acc);
                if (it == genome_seqs.end()) continue;
                auto f = fragment_genome(it->second, host_rng);
                host_frags_all.insert(host_frags_all.end(), f.begin(), f.end());
            }
            // Subsample host fragments to target completeness
            std::vector<std::string> host_frags;
            if (comp >= 1.0f) {
                host_frags = host_frags_all;
            } else {
                uint64_t total_bp = 0;
                for (const auto& f : host_frags_all) total_bp += f.size();
                uint64_t target_bp = static_cast<uint64_t>(static_cast<double>(total_bp) * comp);
                std::mt19937 comp_rng(seed + static_cast<uint32_t>(rep) * 777 +
                                      static_cast<uint32_t>(comp * 1000));
                std::shuffle(host_frags_all.begin(), host_frags_all.end(), comp_rng);
                uint64_t acc = 0;
                for (auto& f : host_frags_all) {
                    if (acc >= target_bp) break;
                    acc += f.size();
                    host_frags.push_back(std::move(f));
                }
            }
            uint64_t host_bp = 0;
            for (const auto& f : host_frags) host_bp += f.size();

            std::mt19937 spike_rng(seed + static_cast<uint32_t>(rep) * 9999 +
                                   static_cast<uint32_t>(std::hash<std::string>{}(mrow.ani_label)));
            std::vector<std::string> spike_frags;
            for (const auto& acc : spike_accs) {
                auto it = genome_seqs.find(acc);
                if (it == genome_seqs.end()) continue;
                auto f = fragment_genome(it->second, spike_rng);
                spike_frags.insert(spike_frags.end(), f.begin(), f.end());
            }
            std::shuffle(spike_frags.begin(), spike_frags.end(), spike_rng);

            for (float frac : SPIKE_FRACS) {
                std::vector<std::string_view> bin_seqs;
                bin_seqs.reserve(host_frags.size() + spike_frags.size());
                for (const auto& f : host_frags) bin_seqs.push_back(f);

                float actual_frac_pct = 0.f;
                if (frac > 0.f && !spike_frags.empty()) {
                    uint64_t target_bp = static_cast<uint64_t>(
                        static_cast<double>(host_bp) * frac / (1.0 - frac));
                    uint64_t added = 0;
                    for (const auto& f : spike_frags) {
                        if (added >= target_bp) break;
                        bin_seqs.push_back(f);
                        added += f.size();
                    }
                    actual_frac_pct = 100.f * static_cast<float>(added) /
                                      static_cast<float>(host_bp + added);
                }

                auto sc = score_bin(bin_seqs, ctx.gcov_e, ctx.gcov_rd,
                                    ctx.fcov_e, ctx.fcov_rd, score_cfg);
                // Fiedler signals are dead (confirmed) — skip score_bin_containment
                float cs_pct = NAN, cs_self = NAN, cs_fiedler = NAN;
                float cs_tnf_bimod = NAN, cs_tnf_gap = NAN;

                // ── FMH sweep (2 contig scans for all 6 combos) ───────────────
                std::array<float, N_FMH> fmh_scores;
                fmh_scores.fill(NAN);
                if (fmh_ok) {
                    auto raw = score_bin_fmh_sweep(bin_seqs,
                                                   host_fmh_refs, contam_fmh_refs);
                    for (int i = 0; i < N_FMH; ++i)
                        fmh_scores[i] = std::isnan(raw[i]) ? NAN : raw[i] * 100.f;
                }

                rows.push_back({mrow.host_genus, mrow.contam_genus, mrow.ani_label,
                                rep, frac * 100.f, actual_frac_pct,
                                std::isnan(sc.cco_fraction)      ? NAN : sc.cco_fraction      * 100.f,
                                std::isnan(sc.sibling_fraction)  ? NAN : sc.sibling_fraction  * 100.f,
                                std::isnan(sc.rho_fraction)      ? NAN : sc.rho_fraction      * 100.f,
                                cs_pct, cs_self, cs_fiedler, cs_tnf_bimod, cs_tnf_gap,
                                fmh_ok ? fmh_scores[1] : NAN,
                                fmh_scores,
                                comp * 100.f,
                                sc.scored_bp});

                if (++done % 20 == 0)
                    spdlog::info("bench-grid: {}/{}", done, n_total);
            }
        }
        } // comp loop
    }

    // ── TSV output ─────────────────────────────────────────────────────────────
    std::ofstream tsv(output);
    if (!tsv) throw std::runtime_error("bench-grid: cannot open " + output.string());

    tsv << "host_genus\tcontam_genus\tani_label\trep\tspike_frac\tactual_frac"
           "\tcco\tsibling_outlier\trho_outlier\tcontig_split\tself_outlier"
           "\tfiedler_oph\tfiedler_tnf_bimod\tfiedler_tnf_gap"
           "\tfmh_minority"
           "\tfmh_k21_c30\tfmh_k21_c125\tfmh_k21_c500"
           "\tfmh_k31_c30\tfmh_k31_c125\tfmh_k31_c500"
           "\tcompleteness_target"
           "\tscored_bp\n";
    auto fv = [](float v) { return std::isnan(v) ? std::string("NA") : std::to_string(v); };
    for (const auto& r : rows) {
        tsv << r.host_genus   << '\t' << r.contam_genus  << '\t' << r.ani_label << '\t'
            << r.rep          << '\t'
            << fv(r.spike_frac_pct)  << '\t' << fv(r.actual_frac_pct) << '\t'
            << fv(r.cco)             << '\t' << fv(r.sibling)         << '\t'
            << fv(r.rho)             << '\t' << fv(r.contig_split)    << '\t'
            << fv(r.self_outlier)    << '\t' << fv(r.fiedler_oph)     << '\t'
            << fv(r.fiedler_tnf_bimod) << '\t' << fv(r.fiedler_tnf_gap) << '\t'
            << fv(r.fmh_minority) << '\t'
            << fv(r.fmh[0]) << '\t' << fv(r.fmh[1]) << '\t' << fv(r.fmh[2]) << '\t'
            << fv(r.fmh[3]) << '\t' << fv(r.fmh[4]) << '\t' << fv(r.fmh[5]) << '\t'
            << fv(r.completeness_target) << '\t'
            << r.scored_bp    << '\n';
    }
    spdlog::info("bench-grid: results → {}", output.string());
    return 0;
}

} // namespace genopack::bench
