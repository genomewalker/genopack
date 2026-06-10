#include <genopack/markers_build.hpp>
#include <genopack/gcov.hpp>     // GcovWriter::hash_genus (FNV-1a)
#include <genopack/aamer.hpp>
#include <genopack/markers.hpp>  // RedunCalibEntry, REDUN_BORROWED_FLAG

#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <spdlog/spdlog.h>

namespace genopack {

// ── Marker column ranges (GTDB r232, PFAM then TIGRFAM order) ────────────────

struct MarkerRange { const char* name; int col_start; int col_end; };

static constexpr MarkerRange k_bac120_ranges[] = {
    {"PF00380.20", 0, 121},
    {"PF00410.20", 121, 247},
    {"PF00466.21", 247, 347},
    {"PF01025.20", 347, 513},
    {"PF02576.18", 513, 586},
    {"PF03726.15", 586, 669},
    {"TIGR00006", 669, 979},
    {"TIGR00019", 979, 1340},
    {"TIGR00020", 1340, 1705},
    {"TIGR00029", 1705, 1792},
    {"TIGR00043", 1792, 1903},
    {"TIGR00054", 1903, 2324},
    {"TIGR00059", 2324, 2436},
    {"TIGR00061", 2436, 2537},
    {"TIGR00064", 2537, 2816},
    {"TIGR00065", 2816, 3169},
    {"TIGR00082", 3169, 3284},
    {"TIGR00083", 3284, 3574},
    {"TIGR00084", 3574, 3766},
    {"TIGR00086", 3766, 3910},
    {"TIGR00088", 3910, 4143},
    {"TIGR00090", 4143, 4242},
    {"TIGR00092", 4242, 4610},
    {"TIGR00095", 4610, 4804},
    {"TIGR00115", 4804, 5214},
    {"TIGR00116", 5214, 5507},
    {"TIGR00138", 5507, 5690},
    {"TIGR00158", 5690, 5838},
    {"TIGR00166", 5838, 5933},
    {"TIGR00168", 5933, 6098},
    {"TIGR00186", 6098, 6338},
    {"TIGR00194", 6338, 6912},
    {"TIGR00250", 6912, 7042},
    {"TIGR00337", 7042, 7568},
    {"TIGR00344", 7568, 8415},
    {"TIGR00362", 8415, 8852},
    {"TIGR00382", 8852, 9266},
    {"TIGR00392", 9266, 10127},
    {"TIGR00396", 10127, 10970},
    {"TIGR00398", 10970, 11500},
    {"TIGR00414", 11500, 11918},
    {"TIGR00416", 11918, 12372},
    {"TIGR00420", 12372, 12723},
    {"TIGR00431", 12723, 12933},
    {"TIGR00435", 12933, 13399},
    {"TIGR00436", 13399, 13669},
    {"TIGR00442", 13669, 14075},
    {"TIGR00445", 14075, 14396},
    {"TIGR00456", 14396, 14965},
    {"TIGR00459", 14965, 15551},
    {"TIGR00460", 15551, 15866},
    {"TIGR00468", 15866, 16190},
    {"TIGR00472", 16190, 16988},
    {"TIGR00487", 16988, 17575},
    {"TIGR00496", 17575, 17751},
    {"TIGR00539", 17751, 18112},
    {"TIGR00580", 18112, 19035},
    {"TIGR00593", 19035, 19925},
    {"TIGR00615", 19925, 20121},
    {"TIGR00631", 20121, 20779},
    {"TIGR00634", 20779, 21342},
    {"TIGR00635", 21342, 21647},
    {"TIGR00643", 21647, 22276},
    {"TIGR00663", 22276, 22643},
    {"TIGR00717", 22643, 23159},
    {"TIGR00755", 23159, 23415},
    {"TIGR00810", 23415, 23488},
    {"TIGR00922", 23488, 23660},
    {"TIGR00928", 23660, 24096},
    {"TIGR00959", 24096, 24524},
    {"TIGR00963", 24524, 25311},
    {"TIGR00964", 25311, 25368},
    {"TIGR00967", 25368, 25782},
    {"TIGR01009", 25782, 25994},
    {"TIGR01011", 25994, 26219},
    {"TIGR01017", 26219, 26419},
    {"TIGR01021", 26419, 26575},
    {"TIGR01029", 26575, 26729},
    {"TIGR01032", 26729, 26843},
    {"TIGR01039", 26843, 27305},
    {"TIGR01044", 27305, 27408},
    {"TIGR01059", 27408, 28047},
    {"TIGR01063", 28047, 28847},
    {"TIGR01066", 28847, 28988},
    {"TIGR01071", 28988, 29132},
    {"TIGR01079", 29132, 29236},
    {"TIGR01082", 29236, 29685},
    {"TIGR01087", 29685, 30126},
    {"TIGR01128", 30126, 30440},
    {"TIGR01146", 30440, 30726},
    {"TIGR01164", 30726, 30852},
    {"TIGR01169", 30852, 31079},
    {"TIGR01171", 31079, 31354},
    {"TIGR01302", 31354, 31804},
    {"TIGR01391", 31804, 32218},
    {"TIGR01393", 32218, 32813},
    {"TIGR01394", 32813, 33407},
    {"TIGR01510", 33407, 33562},
    {"TIGR01632", 33562, 33702},
    {"TIGR01951", 33702, 33833},
    {"TIGR01953", 33833, 34173},
    {"TIGR02012", 34173, 34494},
    {"TIGR02013", 34494, 35732},
    {"TIGR02027", 35732, 36030},
    {"TIGR02075", 36030, 36263},
    {"TIGR02191", 36263, 36482},
    {"TIGR02273", 36482, 36648},
    {"TIGR02350", 36648, 37244},
    {"TIGR02386", 37244, 38391},
    {"TIGR02397", 38391, 38746},
    {"TIGR02432", 38746, 38935},
    {"TIGR02729", 38935, 39264},
    {"TIGR03263", 39264, 39444},
    {"TIGR03594", 39444, 39876},
    {"TIGR03625", 39876, 40078},
    {"TIGR03632", 40078, 40195},
    {"TIGR03654", 40195, 40370},
    {"TIGR03723", 40370, 40684},
    {"TIGR03725", 40684, 40896},
    {"TIGR03953", 40896, 41084},
};
static_assert(std::size(k_bac120_ranges) == 120);

static constexpr MarkerRange k_ar53_ranges[] = {
    {"PF04919.13", 0, 181},
    {"PF07541.13", 181, 295},
    {"PF01000.27", 295, 407},
    {"PF00687.22", 407, 611},
    {"PF00466.21", 611, 711},
    {"PF00827.18", 711, 900},
    {"PF01280.21", 900, 1044},
    {"PF01090.20", 1044, 1181},
    {"PF01200.19", 1181, 1245},
    {"PF01015.19", 1245, 1440},
    {"PF00900.21", 1440, 1515},
    {"PF00410.20", 1515, 1641},
    {"TIGR00037", 1641, 1771},
    {"TIGR00064", 1771, 2050},
    {"TIGR00111", 2050, 2401},
    {"TIGR00134", 2401, 3023},
    {"TIGR00279", 3023, 3195},
    {"TIGR00291", 3195, 3426},
    {"TIGR00323", 3426, 3641},
    {"TIGR00335", 3641, 3965},
    {"TIGR00373", 3965, 4127},
    {"TIGR00405", 4127, 4272},
    {"TIGR00448", 4272, 4451},
    {"TIGR00483", 4451, 4878},
    {"TIGR00491", 4878, 5472},
    {"TIGR00522", 5472, 5730},
    {"TIGR00967", 5730, 6144},
    {"TIGR00982", 6144, 6283},
    {"TIGR01008", 6283, 6478},
    {"TIGR01012", 6478, 6674},
    {"TIGR01018", 6674, 6836},
    {"TIGR01020", 6836, 7048},
    {"TIGR01028", 7048, 7234},
    {"TIGR01046", 7234, 7333},
    {"TIGR01052", 7333, 7821},
    {"TIGR01171", 7821, 8096},
    {"TIGR01213", 8096, 8483},
    {"TIGR01952", 8483, 8624},
    {"TIGR02236", 8624, 8935},
    {"TIGR02338", 8935, 9045},
    {"TIGR02389", 9045, 9412},
    {"TIGR02390", 9412, 10279},
    {"TIGR03626", 10279, 10610},
    {"TIGR03627", 10610, 10740},
    {"TIGR03628", 10740, 10857},
    {"TIGR03629", 10857, 11001},
    {"TIGR03670", 11001, 11600},
    {"TIGR03671", 11600, 12010},
    {"TIGR03672", 12010, 12261},
    {"TIGR03673", 12261, 12392},
    {"TIGR03674", 12392, 12730},
    {"TIGR03676", 12730, 13133},
    {"TIGR03680", 13133, 13540},
};
static_assert(std::size(k_ar53_ranges) == 53);

// ── Helpers ───────────────────────────────────────────────────────────────────

// Extract genus name from GTDB taxonomy string (e.g. "d__Bac;...;g__Escherichia;s__...").
static std::string_view extract_genus(std::string_view tax) {
    auto pos = tax.find("g__");
    if (pos == std::string_view::npos) return {};
    auto end = tax.find(';', pos + 3);
    return tax.substr(pos, end == std::string_view::npos ? tax.size() - pos : end - pos);
}

static std::string_view extract_family(std::string_view tax) {
    auto pos = tax.find("f__");
    if (pos == std::string_view::npos) return {};
    auto end = tax.find(';', pos + 3);
    return tax.substr(pos, end == std::string_view::npos ? tax.size() - pos : end - pos);
}

// Load acc→taxonomy map from GTDB-Tk taxonomy TSV (tab-separated: acc \t taxonomy).
static std::unordered_map<std::string, std::string>
load_taxonomy(const std::filesystem::path& tsv) {
    std::unordered_map<std::string, std::string> map;
    std::ifstream f(tsv);
    if (!f) throw std::runtime_error("markers: cannot open taxonomy: " + tsv.string());
    std::string line;
    while (std::getline(f, line)) {
        auto tab = line.find('\t');
        if (tab == std::string::npos) continue;
        map.emplace(line.substr(0, tab), line.substr(tab + 1));
    }
    return map;
}

// Lineage accumulator: tracks per-genus ref count and per-marker detection rate.
struct LineageAcc {
    uint64_t genus_hash;
    uint8_t  domain;
    uint32_t ref_count = 0;
    std::string genus_name;
    std::vector<uint32_t> marker_detected; // [mi] = #ref genomes with ≥1 IC-passing syncmer
};

// ── Per-panel scan ────────────────────────────────────────────────────────────
//
// Strategy: load all genome sequences into a single contiguous buffer in one
// I/O pass, then do 120 cheap in-RAM per-marker loops using sort+unique.
// Peak memory = MSA size (~7.8 GB bac120) + one marker's raw hashes (~500 MB).
// This avoids holding 120 large hash sets simultaneously (~50 GB with
// std::unordered_set).

struct GenomeInfo {
    uint64_t genus_hash;
    uint32_t lineage_idx;
    uint32_t seq_offset; // byte offset into seqbuf (UINT32_MAX = invalid)
};

struct PanelResult {
    std::vector<std::vector<uint64_t>> pool_hashes; // [marker_id] sorted unique
    std::vector<LineageAcc> lineages;
    std::vector<GenomeInfo> ginfo;   // per-genome genus/offset (including invalid)
    std::vector<char>       seqbuf;  // all valid genome sequences flat
    uint8_t n_markers;
    uint8_t domain;
};

static PanelResult scan_panel(
    const std::filesystem::path& msa_path,
    const std::unordered_map<std::string, std::string>& taxonomy,
    const MarkerRange* ranges,
    int n_markers,
    uint8_t domain,
    int threads = 1)
{
    PanelResult res;
    res.n_markers = static_cast<uint8_t>(n_markers);
    res.domain    = domain;
    res.pool_hashes.resize(n_markers);

    const int expected_cols = ranges[n_markers - 1].col_end;

    // ── Pass 1: stream MSA, validate genomes, store sequences flat ───────────
    // seqbuf holds all valid genome sequences back-to-back (each expected_cols bytes).
    // seq_index[i] = byte offset into seqbuf for valid genome i.

    std::vector<GenomeInfo>& ginfo  = res.ginfo;
    std::vector<char>&       seqbuf = res.seqbuf;
    ginfo.reserve(200000);
    seqbuf.reserve((size_t)200000 * expected_cols);

    std::unordered_map<uint64_t, size_t> genus_idx_map;
    size_t n_valid = 0, n_skipped = 0;

    {
        std::ifstream f(msa_path);
        if (!f) throw std::runtime_error("markers: cannot open MSA: " + msa_path.string());

        std::string acc, seq, line;

        auto commit = [&]() {
            if (acc.empty()) return;
            if ((int)seq.size() != expected_cols) { ++n_skipped; ginfo.push_back({0, 0, UINT32_MAX}); return; }
            auto tax_it = taxonomy.find(acc);
            if (tax_it == taxonomy.end())           { ++n_skipped; ginfo.push_back({0, 0, UINT32_MAX}); return; }
            auto genus = extract_genus(tax_it->second);
            if (genus.empty() || genus.size() <= 3)  { ++n_skipped; ginfo.push_back({0, 0, UINT32_MAX}); return; }

            const uint64_t gh = GcovWriter::hash_genus(genus);
            auto [it, inserted] = genus_idx_map.emplace(gh, res.lineages.size());
            if (inserted) {
                res.lineages.push_back({gh, domain, 0, std::string(genus),
                                        std::vector<uint32_t>(n_markers, 0)});
            }
            res.lineages[it->second].ref_count++;

            const uint32_t off = static_cast<uint32_t>(seqbuf.size());
            seqbuf.insert(seqbuf.end(), seq.begin(), seq.end());
            ginfo.push_back({gh, static_cast<uint32_t>(it->second), off});
            ++n_valid;
        };

        while (std::getline(f, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                commit();
                acc = line.substr(1);
                if (auto sp = acc.find(' '); sp != std::string::npos) acc.resize(sp);
                seq.clear();
                seq.reserve(expected_cols + 4);
            } else {
                seq += line;
            }
        }
        commit();
    }

    spdlog::info("markers: {} panel: {} genomes loaded ({} skipped), {} genera, seqbuf={:.1f} GB",
                 domain == MRKR_DOMAIN_BAC ? "bac120" : "ar53",
                 n_valid, n_skipped, res.lineages.size(),
                 static_cast<double>(seqbuf.size()) / 1e9);

    // ── Pass 2+: per-marker in-RAM extraction ────────────────────────────────
    // Each marker is independent — parallelise over markers with OpenMP.
    // Each thread has its own `raw` buffer (declared inside the loop body).

    std::atomic<int> done_count{0};
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(threads)
#endif
    for (int mi = 0; mi < n_markers; ++mi) {
        std::vector<uint64_t> raw;
        raw.reserve(1 << 22); // 4M initial

        for (const auto& g : ginfo) {
            if (g.seq_offset == UINT32_MAX) continue;
            std::string_view cols(seqbuf.data() + g.seq_offset + ranges[mi].col_start,
                                   ranges[mi].col_end - ranges[mi].col_start);
            auto hashes = extract_aamers_aa(cols);
            raw.insert(raw.end(), hashes.begin(), hashes.end());
        }

        std::sort(raw.begin(), raw.end());
        raw.erase(std::unique(raw.begin(), raw.end()), raw.end());
        res.pool_hashes[mi] = std::move(raw);

        int d = ++done_count;
        if (d % 20 == 0 || d == n_markers)
            spdlog::info("markers: {}/{} markers processed, pool[{}]={} hashes",
                         d, n_markers, mi, res.pool_hashes[mi].size());
    }

    return res;
}

// ── Dayhoff-6 profile pool builder ───────────────────────────────────────────
// Builds a compact discriminative pool using Dayhoff-6 k=12 syncmers filtered
// by per-column information content (IC). Only k-mers whose weakest MSA column
// has IC ≥ ic_threshold are retained, giving a sparse, high-quality pool.
// Syncmers (closed, t=0) replace FracMinHash: consistent ~11% density without
// the Poisson empty-window variance that destabilizes per-ORF containment scores.

static PanelResult scan_panel_d6(
    const std::filesystem::path& msa_path,
    const std::unordered_map<std::string, std::string>& taxonomy,
    const MarkerRange* ranges,
    int n_markers,
    uint8_t domain,
    float ic_threshold,
    int threads = 1)
{
    // ── Pass 1: load MSA (identical to scan_panel) ────────────────────────────
    PanelResult res;
    res.n_markers = static_cast<uint8_t>(n_markers);
    res.domain    = domain;
    res.pool_hashes.resize(n_markers);

    const int expected_cols = ranges[n_markers - 1].col_end;
    std::vector<GenomeInfo>& ginfo  = res.ginfo;
    std::vector<char>&       seqbuf = res.seqbuf;
    ginfo.reserve(200000);
    seqbuf.reserve((size_t)200000 * expected_cols);

    std::unordered_map<uint64_t, size_t> genus_idx_map;
    size_t n_valid = 0, n_skipped = 0;
    {
        std::ifstream f(msa_path);
        if (!f) throw std::runtime_error("markers: cannot open MSA: " + msa_path.string());
        std::string acc, seq, line;
        auto commit = [&]() {
            if (acc.empty()) return;
            if ((int)seq.size() != expected_cols) { ++n_skipped; ginfo.push_back({0, 0, UINT32_MAX}); return; }
            auto tax_it = taxonomy.find(acc);
            if (tax_it == taxonomy.end())           { ++n_skipped; ginfo.push_back({0, 0, UINT32_MAX}); return; }
            auto genus = extract_genus(tax_it->second);
            if (genus.empty() || genus.size() <= 3)  { ++n_skipped; ginfo.push_back({0, 0, UINT32_MAX}); return; }
            const uint64_t gh = GcovWriter::hash_genus(genus);
            auto [it, inserted] = genus_idx_map.emplace(gh, res.lineages.size());
            if (inserted) {
                res.lineages.push_back({gh, domain, 0, std::string(genus),
                                        std::vector<uint32_t>(n_markers, 0)});
            }
            res.lineages[it->second].ref_count++;
            const uint32_t off = static_cast<uint32_t>(seqbuf.size());
            seqbuf.insert(seqbuf.end(), seq.begin(), seq.end());
            ginfo.push_back({gh, static_cast<uint32_t>(it->second), off});
            ++n_valid;
        };
        while (std::getline(f, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                commit();
                acc = line.substr(1);
                if (auto sp = acc.find(' '); sp != std::string::npos) acc.resize(sp);
                seq.clear(); seq.reserve(expected_cols + 4);
            } else { seq += line; }
        }
        commit();
    }
    spdlog::info("markers: d6 {}: {} genomes ({} skipped), {} genera",
                 domain == MRKR_DOMAIN_BAC ? "bac120" : "ar53",
                 n_valid, n_skipped, res.lineages.size());

    // ── Pass 2: per-marker Dayhoff-6 profile extraction with IC filter ────────
    // One representative genome per genus — identical for every marker, so build
    // it once here instead of rebuilding it inside each marker's loop body.
    std::unordered_map<uint64_t, uint32_t> genus_rep; // genus_hash → ginfo index
    genus_rep.reserve(res.lineages.size());
    for (uint32_t gi = 0; gi < static_cast<uint32_t>(ginfo.size()); ++gi)
        if (ginfo[gi].seq_offset != UINT32_MAX)
            genus_rep.emplace(ginfo[gi].genus_hash, gi); // first-seen wins

    std::atomic<int> done_count{0};
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(threads)
#endif
    for (int mi = 0; mi < n_markers; ++mi) {
        const int col_start = ranges[mi].col_start;
        const int col_len   = ranges[mi].col_end - col_start;

        // Per-column IC using one representative per genus (phylogenetically de-biased).
        // Raw frequency counting over 189K genomes overweights dense clades (e.g. many
        // Escherichia strains dominate → "conserved" means conserved in Proteobacteria).
        // Using one genome per genus gives ~35K near-independent samples and makes IC
        // reflect cross-clade conservation rather than within-clade redundancy.
        std::vector<float> col_ic(col_len, 0.0f);
        {
            // genus_rep built once above the parallel region; read-only here.
            for (int j = 0; j < col_len; ++j) {
                int cnt[6] = {};
                int n_obs  = 0;
                for (const auto& [gh, gi] : genus_rep) {
                    const uint8_t aa = AA_ENC[static_cast<unsigned char>(
                        seqbuf[ginfo[gi].seq_offset + col_start + j])];
                    if (aa >= 20) continue;
                    cnt[AA_DAYHOFF6[aa]]++;
                    ++n_obs;
                }
                if (n_obs == 0) continue;
                float H = 0.0f;
                for (int g = 0; g < 6; ++g) {
                    if (!cnt[g]) continue;
                    const float p = static_cast<float>(cnt[g]) / n_obs;
                    H -= p * std::log(p);
                }
                col_ic[j] = std::max(0.0f, (std::log(6.0f) - H) / std::log(6.0f));
            }
        }

        // Extract k=12 syncmer k-mers: recode → IC filter → dedup.
        std::vector<uint64_t> raw;
        raw.reserve(1 << 20);

        for (const auto& g : ginfo) {
            if (g.seq_offset == UINT32_MAX) continue;

            // Build gap-skipped Dayhoff-6 segment + MSA column index mapping.
            static thread_local std::vector<uint8_t> d6;
            static thread_local std::vector<int>     cpos;
            d6.clear();
            cpos.clear();
            d6.reserve(col_len);
            cpos.reserve(col_len);
            for (int j = 0; j < col_len; ++j) {
                const uint8_t aa = AA_ENC[static_cast<unsigned char>(
                    seqbuf[g.seq_offset + col_start + j])];
                if (aa >= 20) continue;
                d6.push_back(AA_DAYHOFF6[aa]);
                cpos.push_back(j);
            }

            bool any_hit = false;
            for (int i = 0; i + AAMER_K_D6 <= (int)d6.size(); ++i) {
                if (!aamer_is_syncmer_d6(d6.data() + i)) continue;
                if (aamer_is_low_complexity(d6.data() + i, AAMER_K_D6)) continue;
                // Weakest-link IC over the k-mer's MSA columns.
                float min_ic = 1.0f;
                for (int j = 0; j < AAMER_K_D6; ++j)
                    min_ic = std::min(min_ic, col_ic[cpos[i + j]]);
                if (min_ic < ic_threshold) continue;
                raw.push_back(aamer_hash_d6(d6.data() + i));
                any_hit = true;
            }
            // Track how many reference genomes in this lineage detected marker mi.
            // Safe: each OpenMP thread owns a unique mi, so no data race on [mi].
            if (any_hit) res.lineages[g.lineage_idx].marker_detected[mi]++;
        }

        std::sort(raw.begin(), raw.end());
        raw.erase(std::unique(raw.begin(), raw.end()), raw.end());
        res.pool_hashes[mi] = std::move(raw);

        int d = ++done_count;
        if (d % 20 == 0 || d == n_markers)
            spdlog::info("markers: d6 {}/{} done, pool[{}]={} hashes",
                         d, n_markers, mi, res.pool_hashes[mi].size());
    }
    return res;
}

// ── Cross-family filter ───────────────────────────────────────────────────────
// Remove any hash appearing in ≥2 distinct marker families within the same panel.
//
// Pools are sorted uint64 vectors. We partition the uint64 space into N_PARTS
// slices (by high bits) and process one slice at a time so peak memory is
// O(total_hashes / N_PARTS) instead of O(total_hashes).

static void cross_family_filter(std::vector<std::vector<uint64_t>>& pool_hashes,
                                int threads = 1) {
    constexpr int N_PARTS = 256;

    // Phase 1 (parallel over partitions, read-only on pool_hashes):
    // Each partition scans all pools for the [lo,hi) hash-value slice and
    // collects hashes that appear in ≥2 distinct marker families.
    // Partitions write to their own slot — no contention.
    std::vector<std::vector<uint64_t>> shared_per_part(N_PARTS);

    // Hash top-byte selects the partition, so each partition sees ~1/N_PARTS of
    // the total hashes — reserve up front to avoid repeated rehashing.
    size_t total_hashes = 0;
    for (const auto& ph : pool_hashes) total_hashes += ph.size();
    const size_t est_per_part = total_hashes / N_PARTS + 1;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(threads)
#endif
    for (int part = 0; part < N_PARTS; ++part) {
        const uint64_t lo = (uint64_t)part << 56;
        const uint64_t hi = (part < N_PARTS - 1) ? (((uint64_t)(part + 1) << 56) - 1)
                                                  : UINT64_MAX;

        std::unordered_map<uint64_t, int> seen;
        seen.reserve(est_per_part);
        for (int mi = 0; mi < (int)pool_hashes.size(); ++mi) {
            auto beg = std::lower_bound(pool_hashes[mi].begin(), pool_hashes[mi].end(), lo);
            auto end = std::upper_bound(pool_hashes[mi].begin(), pool_hashes[mi].end(), hi);
            for (auto it = beg; it != end; ++it) {
                auto [sit, ins] = seen.emplace(*it, mi);
                if (!ins && sit->second != mi)
                    shared_per_part[part].push_back(*it);
            }
        }
        auto& sp = shared_per_part[part];
        std::sort(sp.begin(), sp.end());
        sp.erase(std::unique(sp.begin(), sp.end()), sp.end());
    }

    // Merge: partitions cover non-overlapping hash ranges → just concatenate.
    std::vector<uint64_t> shared;
    for (auto& sp : shared_per_part)
        shared.insert(shared.end(), sp.begin(), sp.end());
    shared_per_part.clear();
    shared_per_part.shrink_to_fit();

    if (shared.empty()) return;

    // Phase 2 (parallel over markers, each pool independent):
    // Remove shared hashes using linear merge (both vectors sorted).
    std::atomic<size_t> removed{0};

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(threads)
#endif
    for (int mi = 0; mi < (int)pool_hashes.size(); ++mi) {
        auto& hashes = pool_hashes[mi];
        const size_t old_sz = hashes.size();
        auto out = hashes.begin();
        auto si  = shared.cbegin();
        for (uint64_t h : hashes) {
            while (si != shared.cend() && *si < h) ++si;
            if (si == shared.cend() || *si != h) *out++ = h;
        }
        removed += old_sz - static_cast<size_t>(std::distance(hashes.begin(), out));
        hashes.erase(out, hashes.end());
    }

    spdlog::info("markers: cross-family filter removed {} hashes shared across ≥2 markers",
                 removed.load());
}

// ── Redundancy calibration helpers ────────────────────────────────────────────

static float vec_median(std::vector<float>& v) {
    if (v.empty()) return 0.0f;
    const size_t mid = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + static_cast<ptrdiff_t>(mid), v.end());
    return v[mid];
}

static float vec_mad(std::vector<float>& v, float median) {
    std::vector<float> devs(v.size());
    for (size_t i = 0; i < v.size(); ++i)
        devs[i] = std::abs(v[i] - median);
    const size_t mid = devs.size() / 2;
    std::nth_element(devs.begin(), devs.begin() + static_cast<ptrdiff_t>(mid), devs.end());
    return devs[mid];
}

// Score all reference genomes against the finalized pool and return sorted
// RedunCalibEntry array. Uses FracMinHash-filtered pool + per-genome seqbuf.
static std::vector<RedunCalibEntry>
compute_redun_calib(const PanelResult& panel,
                    const MarkerRange* ranges,
                    int n_markers,
                    uint64_t max_hash,
                    const std::unordered_map<std::string, std::string>& taxonomy)
{
    const int n_m = n_markers;

    // Pool sizes per marker.
    std::vector<uint32_t> pool_sz(n_m, 0);
    for (int mi = 0; mi < n_m; ++mi)
        pool_sz[mi] = static_cast<uint32_t>(panel.pool_hashes[mi].size());

    // Flat merged pool: (hash, marker_id) sorted by hash.
    struct ME { uint64_t hash; uint8_t mi; };
    std::vector<ME> merged;
    {
        size_t total = 0;
        for (int mi = 0; mi < n_m; ++mi) total += panel.pool_hashes[mi].size();
        merged.reserve(total);
        for (int mi = 0; mi < n_m; ++mi)
            for (uint64_t h : panel.pool_hashes[mi])
                merged.push_back({h, static_cast<uint8_t>(mi)});
        std::sort(merged.begin(), merged.end(),
                  [](const ME& a, const ME& b) { return a.hash < b.hash; });
    }

    // Build genus→family hash from taxonomy.
    std::unordered_map<uint64_t, uint64_t> genus_to_family;
    for (const auto& [acc_str, tax] : taxonomy) {
        auto genus  = extract_genus(tax);
        auto family = extract_family(tax);
        if (genus.empty() || genus.size() <= 3) continue;
        const uint64_t gh = GcovWriter::hash_genus(genus);
        if (!family.empty() && family.size() > 3 && !genus_to_family.count(gh))
            genus_to_family[gh] = GcovWriter::hash_genus(family);
    }

    // Accumulate redundancy per genus.
    std::unordered_map<uint64_t, std::vector<float>> genus_vals;
    genus_vals.reserve(panel.lineages.size());

    for (const auto& g : panel.ginfo) {
        if (g.seq_offset == UINT32_MAX) continue;
        const char* row = panel.seqbuf.data() + g.seq_offset;

        std::vector<uint64_t> qmers;
        qmers.reserve(512);
        for (int mi = 0; mi < n_m; ++mi) {
            std::string_view cols(row + ranges[mi].col_start,
                                  static_cast<size_t>(ranges[mi].col_end - ranges[mi].col_start));
            auto hv = extract_aamers_aa(cols);
            for (uint64_t h : hv)
                if (h <= max_hash) qmers.push_back(h);
        }
        std::sort(qmers.begin(), qmers.end());
        qmers.erase(std::unique(qmers.begin(), qmers.end()), qmers.end());
        if (qmers.empty()) continue;

        uint32_t hits[173] = {};
        size_t qi = 0, pi = 0;
        while (qi < qmers.size() && pi < merged.size()) {
            if      (qmers[qi] < merged[pi].hash) ++qi;
            else if (merged[pi].hash < qmers[qi]) ++pi;
            else { hits[merged[pi].mi]++; ++qi; ++pi; }
        }

        float sum = 0.0f; int cnt = 0;
        for (int mi = 0; mi < n_m; ++mi)
            if (pool_sz[mi] > 0 && hits[mi] > 0) {
                sum += static_cast<float>(hits[mi]) / static_cast<float>(pool_sz[mi]);
                ++cnt;
            }
        if (cnt > 0)
            genus_vals[g.genus_hash].push_back(sum / static_cast<float>(cnt));
    }

    // Build family-level distributions for fallback.
    std::unordered_map<uint64_t, std::vector<float>> family_vals;
    for (auto& [gh, vals] : genus_vals) {
        auto fit = genus_to_family.find(gh);
        if (fit != genus_to_family.end())
            for (float v : vals)
                family_vals[fit->second].push_back(v);
    }
    // Compute family median+MAD.
    std::unordered_map<uint64_t, std::pair<float,float>> family_stats;
    for (auto& [fh, vals] : family_vals) {
        float med = vec_median(vals);
        float mad = vec_mad(vals, med);
        family_stats[fh] = {med, mad};
    }

    // Build output entries.
    std::vector<RedunCalibEntry> out;
    out.reserve(genus_vals.size());
    for (auto& [gh, vals] : genus_vals) {
        RedunCalibEntry e{};
        e.genus_hash = gh;
        const uint16_t n_raw = static_cast<uint16_t>(std::min<size_t>(vals.size(), 0x7FFFu));
        bool borrowed = false;

        if (vals.size() < 10) {
            // Fall back to family stats.
            auto fit = genus_to_family.find(gh);
            if (fit != genus_to_family.end()) {
                auto sit = family_stats.find(fit->second);
                if (sit != family_stats.end()) {
                    e.median = sit->second.first;
                    e.mad    = sit->second.second;
                    borrowed = true;
                }
            }
        }
        if (!borrowed) {
            e.median = vec_median(vals);
            e.mad    = vec_mad(vals, e.median);
        }
        // MAD floor.
        e.mad    = std::max({e.mad, 0.01f * e.median, 1e-4f});
        e.n      = n_raw | (borrowed ? REDUN_BORROWED_FLAG : 0u);
        out.push_back(e);
    }

    std::sort(out.begin(), out.end(),
              [](const RedunCalibEntry& a, const RedunCalibEntry& b) {
                  return a.genus_hash < b.genus_hash;
              });
    return out;
}

// ── build_markers_panel ───────────────────────────────────────────────────────

void build_markers_panel(const MarkersBuildConfig& cfg) {
    const auto db = cfg.gtdbtk_db;

    // Load taxonomy.
    spdlog::info("markers: loading taxonomy");
    auto bac_tax = load_taxonomy(db / "taxonomy/bac120_taxonomy_r232_reps.tsv");
    auto arc_tax = load_taxonomy(db / "taxonomy/ar53_taxonomy_r232_reps.tsv");
    spdlog::info("markers: bac120 taxonomy {} entries, ar53 {} entries",
                 bac_tax.size(), arc_tax.size());

    // Scan bac120 MSA.
    spdlog::info("markers: scanning bac120 MSA ({}) [{}]",
                 (db / "msa/gtdb_r232_bac120.faa").string(),
                 cfg.dayhoff6 ? "dayhoff6+syncmer+IC" : "full-AA");
    auto bac = cfg.dayhoff6
        ? scan_panel_d6(db / "msa/gtdb_r232_bac120.faa",
                        bac_tax, k_bac120_ranges, 120, MRKR_DOMAIN_BAC,
                        cfg.ic_threshold, cfg.threads)
        : scan_panel(db / "msa/gtdb_r232_bac120.faa",
                     bac_tax, k_bac120_ranges, 120, MRKR_DOMAIN_BAC, cfg.threads);

    // Scan ar53 MSA.
    spdlog::info("markers: scanning ar53 MSA");
    auto arc = cfg.dayhoff6
        ? scan_panel_d6(db / "msa/gtdb_r232_ar53.faa",
                        arc_tax, k_ar53_ranges, 53, MRKR_DOMAIN_ARC,
                        cfg.ic_threshold, cfg.threads)
        : scan_panel(db / "msa/gtdb_r232_ar53.faa",
                     arc_tax, k_ar53_ranges, 53, MRKR_DOMAIN_ARC, cfg.threads);

    // Cross-family filter: removes k-mers shared across ≥2 marker families.
    // Essential for Murphy10 — without it, first-writer collision in FlatHitMap
    // makes most markers undetectable. Murphy10 (100M space) survives better than Dayhoff6.
    spdlog::info("markers: cross-family filter (bac120)");
    cross_family_filter(bac.pool_hashes, cfg.threads);
    spdlog::info("markers: cross-family filter (ar53)");
    cross_family_filter(arc.pool_hashes, cfg.threads);

    // FracMinHash: discard hashes above threshold (full-AA mode only;
    // dayhoff6 pools are already sparse via IC + syncmer filtering).
    if (!cfg.dayhoff6 && cfg.frac_scale > 1) {
        const uint64_t max_hash = UINT64_MAX / static_cast<uint64_t>(cfg.frac_scale);
        for (auto& pool : bac.pool_hashes)
            pool.erase(std::upper_bound(pool.begin(), pool.end(), max_hash), pool.end());
        for (auto& pool : arc.pool_hashes)
            pool.erase(std::upper_bound(pool.begin(), pool.end(), max_hash), pool.end());
    }

    // Log pool sizes.
    size_t bac_total = 0, arc_total = 0;
    for (auto& v : bac.pool_hashes) bac_total += v.size();
    for (auto& v : arc.pool_hashes) arc_total += v.size();
    spdlog::info("markers: bac120 pool: {} total hashes across {} markers",
                 bac_total, bac.pool_hashes.size());
    spdlog::info("markers: ar53 pool: {} total hashes across {} markers",
                 arc_total, arc.pool_hashes.size());

    // Redundancy calibration skipped: per-ORF argmax scoring gives ~0 baseline
    // redundancy on clean genomes by construction, so z-score normalization adds nothing.
    std::vector<RedunCalibEntry> all_calib;

    bac.seqbuf.clear(); bac.seqbuf.shrink_to_fit();
    arc.seqbuf.clear(); arc.seqbuf.shrink_to_fit();
    bac.ginfo.clear();  bac.ginfo.shrink_to_fit();
    arc.ginfo.clear();  arc.ginfo.shrink_to_fit();

    // Build MarkerWriter.
    MarkerWriter writer;

    // Register pool entries (bac markers first, then arc).
    for (int mi = 0; mi < 120; ++mi)
        writer.set_pool_entry(static_cast<uint8_t>(mi), std::move(bac.pool_hashes[mi]));
    for (int mi = 0; mi < 53; ++mi)
        writer.set_pool_entry(static_cast<uint8_t>(120 + mi), std::move(arc.pool_hashes[mi]));

    // Add bac lineages.
    const uint16_t thresh_u16 = static_cast<uint16_t>(cfg.default_threshold * 65535.0f);
    size_t bac_expected_set = 0, bac_expected_total = 0;
    for (const auto& lin : bac.lineages) {
        std::vector<CalibSlot> slots(120, CalibSlot{0, thresh_u16});
        std::vector<bool> expected(120);
        for (int mi = 0; mi < 120; ++mi) {
            const bool det = lin.ref_count > 0 &&
                static_cast<float>(lin.marker_detected[mi]) / static_cast<float>(lin.ref_count)
                    >= cfg.expected_min_frac;
            expected[mi] = det;
            bac_expected_set    += det ? 1 : 0;
            bac_expected_total  += 1;
        }
        writer.add_lineage(lin.genus_hash, MRKR_DOMAIN_BAC,
                           static_cast<uint16_t>(std::min<uint32_t>(lin.ref_count, 65535)),
                           slots, expected);
    }
    spdlog::info("markers: bac expected bits set: {}/{} ({:.1f}% of 120 per genus on avg)",
                 bac_expected_set, bac_expected_total,
                 bac.lineages.empty() ? 0.0
                     : 100.0 * static_cast<double>(bac_expected_set)
                             / static_cast<double>(bac.lineages.size() * 120));

    // Add arc lineages.
    for (const auto& lin : arc.lineages) {
        std::vector<CalibSlot> slots(53, CalibSlot{0, thresh_u16});
        std::vector<bool> expected(53);
        for (int mi = 0; mi < 53; ++mi)
            expected[mi] = lin.ref_count > 0 &&
                static_cast<float>(lin.marker_detected[mi]) / static_cast<float>(lin.ref_count)
                    >= cfg.expected_min_frac;
        writer.add_lineage(lin.genus_hash, MRKR_DOMAIN_ARC,
                           static_cast<uint16_t>(std::min<uint32_t>(lin.ref_count, 65535)),
                           slots, expected);
    }

    writer.set_redun_calib(std::move(all_calib));

    spdlog::info("markers: {} bac + {} arc lineages, writing {}",
                 bac.lineages.size(), arc.lineages.size(), cfg.output.string());

    const uint8_t  pool_k      = cfg.dayhoff6 ? AAMER_K_D6 : AAMER_K;
    const uint16_t pool_scale  = cfg.dayhoff6 ? 1            : cfg.frac_scale;
    const uint8_t  pool_alpha  = cfg.dayhoff6 ? MRKR_ALPHABET_DAYHOFF6 : MRKR_ALPHABET_FULL_AA;
    writer.finalize(cfg.output, 120, 53, pool_k, pool_scale, pool_alpha);
    spdlog::info("markers: done");
}

} // namespace genopack
