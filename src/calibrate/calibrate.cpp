#include "calibrate.hpp"
#include <genopack/archive.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace genopack::calibrate {

namespace {

// ── CheckM2 record ────────────────────────────────────────────────────────────

struct CheckM2Record {
    float completeness   = NAN;
    float contamination  = NAN;
};

// Parse CheckM2 TSV. Accepts any column order; header detection for
// "accession"/"name", "completeness", "contamination" (case-insensitive).
std::unordered_map<std::string, CheckM2Record>
load_checkm2(const std::string& path)
{
    std::ifstream f(path);
    if (!f) throw std::runtime_error("calibrate: cannot open CheckM2 TSV: " + path);

    auto to_lower = [](std::string s) {
        for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return s;
    };

    // Parse header
    std::string header_line;
    if (!std::getline(f, header_line))
        throw std::runtime_error("calibrate: empty CheckM2 TSV: " + path);
    if (!header_line.empty() && header_line.back() == '\r') header_line.pop_back();

    std::vector<std::string> col_names;
    {
        std::istringstream ss(header_line);
        std::string tok;
        while (std::getline(ss, tok, '\t'))
            col_names.push_back(to_lower(tok));
    }

    int acc_col  = -1;
    int comp_col = -1;
    int cont_col = -1;
    for (int i = 0; i < static_cast<int>(col_names.size()); ++i) {
        const auto& c = col_names[i];
        if (acc_col  < 0 && (c == "accession" || c == "name"))    acc_col  = i;
        if (comp_col < 0 && c == "completeness")                   comp_col = i;
        if (cont_col < 0 && c == "contamination")                  cont_col = i;
    }
    if (acc_col  < 0) throw std::runtime_error("calibrate: no accession/name column in CheckM2 TSV");
    if (comp_col < 0) throw std::runtime_error("calibrate: no completeness column in CheckM2 TSV");
    if (cont_col < 0) throw std::runtime_error("calibrate: no contamination column in CheckM2 TSV");

    std::unordered_map<std::string, CheckM2Record> result;
    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;

        std::vector<std::string> fields;
        {
            std::istringstream ss(line);
            std::string tok;
            while (std::getline(ss, tok, '\t'))
                fields.push_back(tok);
        }

        int n = static_cast<int>(fields.size());
        if (acc_col  >= n || comp_col >= n || cont_col >= n) continue;

        CheckM2Record r;
        try { r.completeness  = std::stof(fields[comp_col]); } catch (...) { continue; }
        try { r.contamination = std::stof(fields[cont_col]); } catch (...) { r.contamination = NAN; }
        result.emplace(fields[acc_col], r);
    }
    spdlog::info("calibrate: loaded {} records from CheckM2 TSV", result.size());
    return result;
}

// ── Isotonic regression (PAV algorithm) ──────────────────────────────────────

struct IsotonicModel {
    std::vector<float> x;   // breakpoints (sorted)
    std::vector<float> y;   // fitted values at breakpoints
};

// Fit isotonic regression: sort by x, pool adjacent violators.
// xy must be non-empty.
IsotonicModel fit_isotonic(std::vector<std::pair<float, float>> xy)
{
    // Filter NaNs
    xy.erase(std::remove_if(xy.begin(), xy.end(),
                             [](const auto& p) {
                                 return std::isnan(p.first) || std::isnan(p.second);
                             }),
             xy.end());
    if (xy.empty()) throw std::runtime_error("calibrate: no valid points for isotonic fit");

    std::sort(xy.begin(), xy.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

    // PAV: maintain stack of pools {sum_x, sum_y, count}
    struct Pool { double sx = 0, sy = 0; int n = 0; };
    std::vector<Pool> pools;
    pools.reserve(xy.size());

    for (const auto& [xi, yi] : xy) {
        Pool p;
        p.sx = xi; p.sy = yi; p.n = 1;
        while (pools.size() >= 2) {
            auto& prev = pools[pools.size() - 2];
            auto& cur  = pools[pools.size() - 1];
            if (cur.sy / cur.n < prev.sy / prev.n) {
                prev.sx += cur.sx; prev.sy += cur.sy; prev.n += cur.n;
                pools.pop_back();
            } else break;
        }
        // Check last pool vs new point
        if (!pools.empty()) {
            auto& last = pools.back();
            if (p.sy / p.n < last.sy / last.n) {
                last.sx += p.sx; last.sy += p.sy; last.n += p.n;
                continue;
            }
        }
        pools.push_back(p);
    }

    IsotonicModel m;
    m.x.reserve(pools.size());
    m.y.reserve(pools.size());
    for (const auto& pool : pools) {
        m.x.push_back(static_cast<float>(pool.sx / pool.n));
        m.y.push_back(static_cast<float>(pool.sy / pool.n));
    }
    return m;
}

// Linear interpolation using isotonic model breakpoints.
float apply_isotonic(const IsotonicModel& m, float x)
{
    if (m.x.empty()) return NAN;
    if (x <= m.x.front()) return m.y.front();
    if (x >= m.x.back())  return m.y.back();

    // Binary search for interval
    auto it = std::lower_bound(m.x.begin(), m.x.end(), x);
    size_t hi = static_cast<size_t>(it - m.x.begin());
    size_t lo = hi - 1;

    float dx = m.x[hi] - m.x[lo];
    if (dx < 1e-12f) return m.y[lo];
    float t = (x - m.x[lo]) / dx;
    return m.y[lo] + t * (m.y[hi] - m.y[lo]);
}

// ── OLS with ridge regularisation (Gaussian elimination) ─────────────────────
// X: n×d matrix (row-major), y: n-vector. Returns d-vector of weights.
std::vector<float> fit_ols(const std::vector<std::vector<float>>& X,
                            const std::vector<float>& y,
                            float lambda = 1e-4f)
{
    if (X.empty()) return {};
    const int n = static_cast<int>(X.size());
    const int d = static_cast<int>(X[0].size());

    // Build XTX (d×d) and XTy (d×1)
    std::vector<double> XTX(d * d, 0.0);
    std::vector<double> XTy(d, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < d; ++j) {
            XTy[j] += static_cast<double>(X[i][j]) * static_cast<double>(y[i]);
            for (int k = 0; k < d; ++k)
                XTX[j * d + k] += static_cast<double>(X[i][j]) * static_cast<double>(X[i][k]);
        }
    }

    // Ridge: add lambda to diagonal
    for (int j = 0; j < d; ++j)
        XTX[j * d + j] += static_cast<double>(lambda);

    // Gaussian elimination with partial pivoting
    // Augmented matrix [XTX | XTy]
    std::vector<double> A(d * (d + 1));
    for (int r = 0; r < d; ++r) {
        for (int c = 0; c < d; ++c)
            A[r * (d + 1) + c] = XTX[r * d + c];
        A[r * (d + 1) + d] = XTy[r];
    }

    for (int col = 0; col < d; ++col) {
        // Find pivot
        int pivot = col;
        double best = std::abs(A[col * (d + 1) + col]);
        for (int row = col + 1; row < d; ++row) {
            double v = std::abs(A[row * (d + 1) + col]);
            if (v > best) { best = v; pivot = row; }
        }
        if (pivot != col) {
            for (int c = 0; c <= d; ++c)
                std::swap(A[col * (d + 1) + c], A[pivot * (d + 1) + c]);
        }
        double diag = A[col * (d + 1) + col];
        if (std::abs(diag) < 1e-15) continue;
        for (int row = 0; row < d; ++row) {
            if (row == col) continue;
            double factor = A[row * (d + 1) + col] / diag;
            for (int c = col; c <= d; ++c)
                A[row * (d + 1) + c] -= factor * A[col * (d + 1) + c];
        }
    }

    std::vector<float> w(d);
    for (int j = 0; j < d; ++j) {
        double diag = A[j * (d + 1) + j];
        w[j] = (std::abs(diag) < 1e-15) ? 0.0f
                                         : static_cast<float>(A[j * (d + 1) + d] / diag);
    }
    return w;
}

// ── Feature names (must match order used in run_calibrate) ───────────────────
// Isotonic axis is features[0] — completeness_cluster_relative (always stored
// in QualRecord). self_coherence added as OLS feature when available.
static const std::vector<std::string> k_completeness_features = {
    "completeness_cluster_relative",   // [0] isotonic axis
    "self_coherence",
    "contamination_leakage",
    "gc_pct",
    "log1p_genome_length",
    "log1p_n_contigs",
    "bias"
};

static const std::vector<std::string> k_contamination_features = {
    "contamination_leakage",
    "contamination_tnf_excess",
    "gc_pct",
    "log1p_genome_length",
    "log1p_n_contigs",
    "bias"
};

} // namespace

// ── Main calibration entry point ──────────────────────────────────────────────

int run_calibrate(const std::string& archive_path,
                  const std::string& checkm2_path,
                  const std::string& output_path,
                  int /*threads*/)
{
    // 1. Load CheckM2 ground truth
    auto checkm2 = load_checkm2(checkm2_path);
    if (checkm2.empty())
        throw std::runtime_error("calibrate: CheckM2 TSV is empty");

    // 2. Open archive — reads only index, no FASTA decompression
    ArchiveReader ar;
    ar.open(archive_path);

    if (!ar.has_qual())
        throw std::runtime_error("calibrate: archive has no QUAL section — run 'check' first");

    // 3. Build genome_id → accession map (sequential index scan, O(n))
    std::unordered_map<uint64_t, std::string> id_to_acc;
    id_to_acc.reserve(ar.archive_stats().n_genomes_total);
    ar.scan_genome_accessions([&](std::string_view acc, GenomeId gid) {
        id_to_acc.emplace(static_cast<uint64_t>(gid), acc);
    });
    spdlog::info("calibrate: {} genomes in archive", id_to_acc.size());

    // 4. Scan QUAL section — flat mmapped array, instant
    const int dc = static_cast<int>(k_completeness_features.size());
    const int dp = static_cast<int>(k_contamination_features.size());

    struct Sample {
        float completeness   = NAN;
        float contamination  = NAN;
        std::vector<float> comp_feat;   // dc features
        std::vector<float> cont_feat;   // dp features
    };
    std::vector<Sample> samples;
    samples.reserve(checkm2.size());

    ar.scan_qual([&](const QualRecord& r) {
        auto ait = id_to_acc.find(r.genome_id);
        if (ait == id_to_acc.end()) return;
        auto cit = checkm2.find(ait->second);
        if (cit == checkm2.end()) return;

        auto meta = ar.genome_meta(r.genome_id);
        if (!meta) return;

        if (std::isnan(r.completeness_cluster_relative)) return; // need isotonic axis

        const float ccr = r.completeness_cluster_relative / 100.0f;
        const float sc  = (std::isnan(r.self_coherence) || r.self_coherence == 0.0f)
                          ? 0.0f : r.self_coherence;
        const float cl  = std::isnan(r.contamination_leakage) ? 0.0f
                          : r.contamination_leakage;
        const float cte = std::isnan(r.contamination_tnf_excess) ? 0.0f
                          : r.contamination_tnf_excess;
        const float gc  = static_cast<float>(meta->gc_pct_x100) / 10000.0f;
        const float lgl = std::log1p(static_cast<float>(meta->genome_length));
        const float lnc = std::log1p(static_cast<float>(meta->n_contigs));

        Sample s;
        s.completeness  = cit->second.completeness  / 100.0f;
        s.contamination = cit->second.contamination / 100.0f;

        // completeness features — axis is completeness_cluster_relative [0]
        s.comp_feat.resize(dc);
        s.comp_feat[0] = ccr;   // isotonic axis
        s.comp_feat[1] = sc;
        s.comp_feat[2] = cl;
        s.comp_feat[3] = gc;
        s.comp_feat[4] = lgl;
        s.comp_feat[5] = lnc;
        s.comp_feat[6] = 1.0f;

        // contamination features — require contamination_leakage
        if (!std::isnan(r.contamination_leakage)) {
            s.cont_feat.resize(dp);
            s.cont_feat[0] = cl;
            s.cont_feat[1] = cte;
            s.cont_feat[2] = gc;
            s.cont_feat[3] = lgl;
            s.cont_feat[4] = lnc;
            s.cont_feat[5] = 1.0f;
        }

        if (!s.comp_feat.empty() || !s.cont_feat.empty())
            samples.push_back(std::move(s));
    });

    ar.close();
    spdlog::info("calibrate: {} samples matched (QUAL + CheckM2)", samples.size());

    // Helper: fit + evaluate one isotonic+OLS model, log results, return {iso, ols_w, rmse}
    auto fit_model = [&](const char* label,
                         std::vector<std::pair<float,float>> iso_xy,
                         std::vector<std::vector<float>> X,
                         std::vector<float> y,
                         int feat_dim)
        -> std::tuple<IsotonicModel, std::vector<float>, float>
    {
        spdlog::info("calibrate: {} — {} samples", label, iso_xy.size());
        if (iso_xy.empty())
            throw std::runtime_error(std::string("calibrate: no usable samples for ") + label);

        IsotonicModel iso = fit_isotonic(iso_xy);

        std::vector<float> ols_y_resid = y;
        for (size_t i = 0; i < ols_y_resid.size(); ++i)
            ols_y_resid[i] -= apply_isotonic(iso, X[i][0]);
        std::vector<float> ols_w = fit_ols(X, ols_y_resid, 1e-4f);

        // RMSE + per-decile
        double sse = 0.0; int n = 0;
        struct Bin { int n = 0; double sse = 0; };
        std::array<Bin, 10> bins{};
        for (size_t i = 0; i < X.size(); ++i) {
            float iso_p = apply_isotonic(iso, X[i][0]);
            float ols_p = 0.0f;
            for (int j = 0; j < feat_dim; ++j) ols_p += ols_w[j] * X[i][j];
            double e = static_cast<double>((iso_p + ols_p) - y[i]);
            sse += e * e; ++n;
            int b = std::clamp(static_cast<int>(y[i] * 10.0f), 0, 9);
            bins[b].n++; bins[b].sse += e * e;
        }
        float rmse = n > 0 ? static_cast<float>(std::sqrt(sse / n)) : NAN;
        spdlog::info("calibrate: {} RMSE = {:.4f} (n={})", label, rmse, n);
        spdlog::info("calibrate: {} per-decile RMSE:", label);
        for (int b = 0; b < 10; ++b) {
            float rb = bins[b].n > 0 ? static_cast<float>(std::sqrt(bins[b].sse / bins[b].n)) : NAN;
            spdlog::info("  [{:.0f}%-{:.0f}%]: n={} rmse={:.4f}",
                         b * 10.0f, (b+1) * 10.0f, bins[b].n, rb);
        }
        return {std::move(iso), std::move(ols_w), rmse};
    };

    // ── Completeness model (isotonic axis: self_coherence) ────────────────────
    std::vector<std::pair<float,float>> comp_iso_xy;
    std::vector<std::vector<float>>     comp_X;
    std::vector<float>                  comp_y;
    for (const auto& s : samples) {
        if (s.comp_feat.empty() || std::isnan(s.completeness)) continue;
        comp_iso_xy.emplace_back(s.comp_feat[0], s.completeness);
        comp_X.push_back(s.comp_feat);
        comp_y.push_back(s.completeness);
    }
    auto [comp_iso, comp_w, comp_rmse] =
        fit_model("completeness", comp_iso_xy, comp_X, comp_y, dc);

    // ── Contamination model (isotonic axis: contamination_leakage) ────────────
    std::vector<std::pair<float,float>> cont_iso_xy;
    std::vector<std::vector<float>>     cont_X;
    std::vector<float>                  cont_y;
    for (const auto& s : samples) {
        if (s.cont_feat.empty() || std::isnan(s.contamination)) continue;
        cont_iso_xy.emplace_back(s.cont_feat[0], s.contamination);
        cont_X.push_back(s.cont_feat);
        cont_y.push_back(s.contamination);
    }
    auto [cont_iso, cont_w, cont_rmse] =
        fit_model("contamination", cont_iso_xy, cont_X, cont_y, dp);

    // ── Write JSON ────────────────────────────────────────────────────────────
    auto write_vec = [](std::ofstream& f, const auto& v) {
        f << "[";
        for (size_t i = 0; i < v.size(); ++i) { if (i) f << ", "; f << v[i]; }
        f << "]";
    };
    auto write_str_vec = [](std::ofstream& f, const std::vector<std::string>& v) {
        f << "[";
        for (size_t i = 0; i < v.size(); ++i) { if (i) f << ", "; f << "\"" << v[i] << "\""; }
        f << "]";
    };

    std::ofstream out(output_path);
    if (!out) throw std::runtime_error("calibrate: cannot write output: " + output_path);
    out << "{\n";
    out << "  \"completeness\": {\n";
    out << "    \"isotonic_x\": "; write_vec(out, comp_iso.x); out << ",\n";
    out << "    \"isotonic_y\": "; write_vec(out, comp_iso.y); out << ",\n";
    out << "    \"ols_weights\": "; write_vec(out, comp_w); out << ",\n";
    out << "    \"feature_names\": "; write_str_vec(out, k_completeness_features); out << ",\n";
    out << "    \"training_rmse\": " << comp_rmse << ",\n";
    out << "    \"n_samples\": " << comp_X.size() << "\n";
    out << "  },\n";
    out << "  \"contamination\": {\n";
    out << "    \"isotonic_x\": "; write_vec(out, cont_iso.x); out << ",\n";
    out << "    \"isotonic_y\": "; write_vec(out, cont_iso.y); out << ",\n";
    out << "    \"ols_weights\": "; write_vec(out, cont_w); out << ",\n";
    out << "    \"feature_names\": "; write_str_vec(out, k_contamination_features); out << ",\n";
    out << "    \"training_rmse\": " << cont_rmse << ",\n";
    out << "    \"n_samples\": " << cont_X.size() << "\n";
    out << "  }\n";
    out << "}\n";

    spdlog::info("calibrate: model written to {}", output_path);
    return 0;
}

} // namespace genopack::calibrate
