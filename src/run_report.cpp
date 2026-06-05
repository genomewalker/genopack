#include "run_report.hpp"
#include <genopack/archive.hpp>
#include <genopack/colstore.hpp>
#include <genopack/prof.hpp>
#include <genopack/qcontig.hpp>
#include <genopack/quality_schema.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <climits>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace genopack {

namespace {

// ── archive enumeration ───────────────────────────────────────────────────────
std::vector<std::filesystem::path> collect_gpk_paths(const std::filesystem::path& pack_path) {
    std::vector<std::filesystem::path> paths;
    if (pack_path.extension() == ".gpk") { paths.push_back(pack_path); return paths; }
    for (const auto& e : std::filesystem::directory_iterator(pack_path))
        if (e.path().extension() == ".gpk") paths.push_back(e.path());
    std::sort(paths.begin(), paths.end());
    return paths;
}

// ── available columns (with provenance + a cached zero-copy view) ─────────────
struct AvailColumn {
    ColumnKey                    key;
    uint64_t                     identity = 0;
    int                          store    = 0;       // 0=qcol, 1=xqual, 2=cqal
    const ColStoreReader*        reader   = nullptr;
    ColStoreReader::ColumnView   view{};
};

const char* store_name(int s) { return s == 0 ? "qcol" : s == 1 ? "xqual" : "cqal"; }

void gather_columns(const ArchiveReader& ar, std::vector<AvailColumn>& out) {
    auto add = [&](const ColStoreReader* r, int store) {
        if (!r) return;
        for (uint32_t c = 0; c < r->n_columns(); ++c) {
            ColStoreReader::ColumnView v = r->column(c);
            AvailColumn a;
            a.key.axis             = std::string(v.axis);
            a.key.tool             = std::string(v.tool);
            a.key.method           = std::string(v.method);
            a.key.unit             = v.unit;
            a.key.version          = v.version;
            a.key.ref_db_hash      = v.ref_db_hash;
            a.key.core_model_hash  = v.core_model_hash;
            a.key.calibration_hash = v.calibration_hash;
            a.key.label            = std::string(v.label);
            a.identity = v.identity_hash;
            a.store    = store;
            a.reader   = r;
            a.view     = v;
            out.push_back(std::move(a));
        }
    };
    add(ar.qcol_reader(), 0);
    add(ar.xqual_reader(), 1);
    add(ar.cqal_reader(), 2);
}

// ── built-in resolution policy ────────────────────────────────────────────────
enum class Policy { Intrinsic, External, Calibrated, Best };

bool parse_policy(const std::string& name, Policy& out) {
    if (name == "intrinsic")  { out = Policy::Intrinsic;  return true; }
    if (name == "external")   { out = Policy::External;   return true; }
    if (name == "calibrated") { out = Policy::Calibrated; return true; }
    if (name == "best")       { out = Policy::Best;       return true; }
    return false;
}

// Preferred method order within an axis (lower index = preferred). Unlisted
// methods rank after all listed ones; ties break toward the lower store index.
int method_rank(const std::string& axis, const std::string& method) {
    static const std::vector<std::string> compl_order = {
        qual_method::AAMER_GENUS_CORE, qual_method::AAMER_FAMILY_CORE,
        qual_method::CLUSTER_RELATIVE, qual_method::SKETCH_FILL, qual_method::MARKER_SCG,
        qual_method::ML, qual_method::FRAGMENTATION, qual_method::POST_DECONTAM};
    static const std::vector<std::string> contam_order = {
        qual_method::LEAKAGE, qual_method::TNF_EXCESS, qual_method::DUPLICATION,
        qual_method::MARKER_REDUNDANCY, qual_method::MIXTURE};
    const std::vector<std::string>* order = nullptr;
    if (axis == qual_axis::COMPLETENESS)       order = &compl_order;
    else if (axis == qual_axis::CONTAMINATION) order = &contam_order;
    if (order)
        for (size_t i = 0; i < order->size(); ++i)
            if (method == (*order)[i]) return static_cast<int>(i);
    return 100;
}

// Tool preference for the external policy (lower = preferred).
int tool_rank(const std::string& tool) {
    if (tool == qual_tool::CHECKM2) return 0;
    if (tool == qual_tool::ANVIO)   return 1;
    return 50;
}

const AvailColumn* pick_intrinsic(const std::vector<AvailColumn>& av, const std::string& axis) {
    const AvailColumn* best = nullptr; int br = INT_MAX, bs = INT_MAX;
    for (const auto& a : av) {
        if (a.key.axis != axis || a.key.tool != qual_tool::GENOPACK || a.key.calibration_hash != 0) continue;
        int r = method_rank(axis, a.key.method);
        if (r < br || (r == br && a.store < bs)) { best = &a; br = r; bs = a.store; }
    }
    return best;
}

const AvailColumn* pick_external(const std::vector<AvailColumn>& av, const std::string& axis) {
    const AvailColumn* best = nullptr; int btr = INT_MAX, bmr = INT_MAX;
    for (const auto& a : av) {
        if (a.key.axis != axis) continue;
        int tr = tool_rank(a.key.tool);
        if (tr >= 50) continue;  // not an external tool
        int mr = method_rank(axis, a.key.method);
        if (tr < btr || (tr == btr && mr < bmr)) { best = &a; btr = tr; bmr = mr; }
    }
    return best;
}

const AvailColumn* pick_calibrated(const std::vector<AvailColumn>& av, const std::string& axis) {
    const AvailColumn* best = nullptr; int br = INT_MAX, bs = INT_MAX;
    for (const auto& a : av) {
        if (a.key.axis != axis || a.key.calibration_hash == 0) continue;
        int r = method_rank(axis, a.key.method);
        if (r < br || (r == br && a.store < bs)) { best = &a; br = r; bs = a.store; }
    }
    return best;
}

const AvailColumn* resolve_policy(const std::vector<AvailColumn>& av, const std::string& axis, Policy p) {
    switch (p) {
        case Policy::Intrinsic:  return pick_intrinsic(av, axis);
        case Policy::External:   return pick_external(av, axis);
        case Policy::Calibrated: return pick_calibrated(av, axis);
        case Policy::Best: {
            if (const AvailColumn* c = pick_calibrated(av, axis)) return c;
            if (const AvailColumn* e = pick_external(av, axis))   return e;
            return pick_intrinsic(av, axis);
        }
    }
    return nullptr;
}

const AvailColumn* by_identity(const std::vector<AvailColumn>& av, uint64_t id) {
    if (id == 0) return nullptr;
    for (const auto& a : av) if (a.identity == id) return &a;
    return nullptr;
}

// ── presentation ──────────────────────────────────────────────────────────────
const char* AXIS_ORDER[] = {
    qual_axis::COMPLETENESS, qual_axis::CONTAMINATION, qual_axis::COHERENCE,
    qual_axis::CHROM_SKEW,   qual_axis::SUPPORT,       qual_axis::INTERVAL};

std::string provenance(const AvailColumn& a) {
    std::string s = a.key.tool + "/" + a.key.method;
    if (a.key.calibration_hash != 0) s += "+cal";
    return s;
}

std::string fmt_value(const AvailColumn* a, double raw) {
    if (!a || std::isnan(raw)) return "NA";
    double v = raw;
    if (a->key.unit == Unit::Fraction01) v = raw * 100.0;  // report completeness/contam as percent
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.4f", v);
    return buf;
}

char hex_lower[] = "0123456789abcdef";
std::string hex64(uint64_t v) {
    std::string s(18, '0'); s[0] = '0'; s[1] = 'x';
    for (int i = 0; i < 16; ++i) s[2 + i] = hex_lower[(v >> ((15 - i) * 4)) & 0xF];
    return s;
}

// One resolved axis: which column (if any) backs it under the active profile.
struct ResolvedAxis {
    std::string         axis;
    const AvailColumn*  col = nullptr;
};

// Emit per-genome rows for one open archive given the resolved axis columns.
void emit_rows(std::ostream& os, const ArchiveReader& ar, const std::vector<ResolvedAxis>& axes) {
    ar.scan_genome_accessions([&](std::string_view acc, GenomeId gid) {
        if (ar.is_deleted(gid)) return;
        std::string tax;
        if (auto t = ar.taxonomy_for_accession(acc)) tax = *t;
        os << acc << '\t' << tax;
        for (const auto& ra : axes) {
            double raw = std::nan("");
            if (ra.col) {
                int64_t row = ra.col->reader->row_index(static_cast<uint64_t>(gid));
                if (row >= 0) raw = ra.col->view.get_f64(static_cast<uint32_t>(row));
            }
            os << '\t' << fmt_value(ra.col, raw)
               << '\t' << (ra.col && !std::isnan(raw) ? provenance(*ra.col) : "-");
        }
        os << '\n';
    });
}

// ── --list ────────────────────────────────────────────────────────────────────
int do_list(const std::vector<std::filesystem::path>& gpks) {
    ArchiveReader ar;
    ar.open(gpks.front());
    std::vector<AvailColumn> av;
    gather_columns(ar, av);

    std::cout << "Built-in profiles (resolved by rule at report time):\n"
                 "  intrinsic   per axis, the genopack-computed (uncalibrated) column\n"
                 "  external    per axis, an external tool (CheckM2 preferred, then anvi'o)\n"
                 "  calibrated  per axis, a column carrying a calibration model\n"
                 "  best        calibrated > external > intrinsic (most trustworthy available)\n\n";

    if (ar.has_prof()) {
        const ProfReader* pr = ar.prof_reader();
        std::cout << "Stored profiles in " << gpks.front().filename().string() << ":\n";
        for (uint32_t i = 0; i < pr->n_profiles(); ++i) {
            Profile p = pr->profile(i);
            std::cout << "  " << p.name << "  (" << p.description << ")  policy_hash="
                      << hex64(pr->policy_hash(i)) << '\n';
            for (const auto& s : p.selections) {
                const AvailColumn* c = by_identity(av, s.column_identity);
                std::cout << "      " << s.axis << " -> "
                          << (c ? provenance(*c) : "MISSING(" + hex64(s.column_identity) + ")");
                if (s.fallback_identity) {
                    const AvailColumn* f = by_identity(av, s.fallback_identity);
                    std::cout << "  | fallback " << (f ? provenance(*f) : "MISSING");
                }
                std::cout << '\n';
            }
        }
        std::cout << '\n';
    } else {
        std::cout << "Stored profiles: none\n\n";
    }

    std::cout << "Available columns (axis  tool/method  store  cal  identity):\n";
    std::sort(av.begin(), av.end(), [](const AvailColumn& a, const AvailColumn& b) {
        if (a.key.axis != b.key.axis) return a.key.axis < b.key.axis;
        return a.store < b.store;
    });
    for (const auto& a : av)
        std::cout << "  " << a.key.axis << '\t' << a.key.tool << '/' << a.key.method
                  << '\t' << store_name(a.store)
                  << '\t' << (a.key.calibration_hash ? "cal" : "raw")
                  << '\t' << hex64(a.identity) << '\n';
    ar.close();
    return 0;
}

} // namespace

// ── report ────────────────────────────────────────────────────────────────────
int cmd_report(const std::filesystem::path& pack_path,
               const std::string& profile_name,
               const std::filesystem::path& output,
               bool list) {
    auto gpks = collect_gpk_paths(pack_path);
    if (gpks.empty()) throw std::runtime_error("report: no .gpk found at " + pack_path.string());
    if (list) return do_list(gpks);

    Policy policy{};
    const bool builtin = parse_policy(profile_name, policy);

    // Determine the ordered axis list + (for stored profiles) the pinned selections.
    std::vector<std::string> axis_order;          // header axes, in order
    std::vector<ProfileSelection> stored_sel;     // only for stored profiles
    if (builtin) {
        // Union of axes that resolve under the policy across all parts (canonical order).
        std::vector<bool> seen(sizeof(AXIS_ORDER) / sizeof(AXIS_ORDER[0]), false);
        for (const auto& gp : gpks) {
            ArchiveReader ar; ar.open(gp);
            std::vector<AvailColumn> av; gather_columns(ar, av);
            for (size_t i = 0; i < seen.size(); ++i)
                if (!seen[i] && resolve_policy(av, AXIS_ORDER[i], policy)) seen[i] = true;
            ar.close();
        }
        for (size_t i = 0; i < seen.size(); ++i)
            if (seen[i]) axis_order.emplace_back(AXIS_ORDER[i]);
        if (axis_order.empty())
            throw std::runtime_error("report: profile '" + profile_name +
                                     "' resolves no axis (no matching quality columns)");
    } else {
        // Stored profile: read it from the first part that carries it.
        bool found = false;
        for (const auto& gp : gpks) {
            ArchiveReader ar; ar.open(gp);
            if (ar.has_prof()) {
                int pi = ar.prof_reader()->find(profile_name);
                if (pi >= 0) {
                    Profile p = ar.prof_reader()->profile(static_cast<uint32_t>(pi));
                    for (const auto& s : p.selections) { axis_order.push_back(s.axis); stored_sel.push_back(s); }
                    found = true;
                }
            }
            ar.close();
            if (found) break;
        }
        if (!found)
            throw std::runtime_error("report: unknown profile '" + profile_name +
                                     "' (not a built-in {intrinsic,external,calibrated,best} "
                                     "and no stored profile by that name; see `report --list`)");
    }

    std::ofstream fout;
    std::ostream* os = &std::cout;
    if (!output.empty()) {
        fout.open(output);
        if (!fout) throw std::runtime_error("report: cannot write " + output.string());
        os = &fout;
    }

    *os << "accession\ttaxonomy";
    for (const auto& a : axis_order) *os << '\t' << a << '\t' << a << "_source";
    *os << '\n';

    size_t n_parts = 0;
    for (const auto& gp : gpks) {
        ArchiveReader ar; ar.open(gp);
        std::vector<AvailColumn> av; gather_columns(ar, av);

        std::vector<ResolvedAxis> axes;
        axes.reserve(axis_order.size());
        for (size_t i = 0; i < axis_order.size(); ++i) {
            ResolvedAxis ra; ra.axis = axis_order[i];
            if (builtin) {
                ra.col = resolve_policy(av, axis_order[i], policy);
            } else {
                ra.col = by_identity(av, stored_sel[i].column_identity);
                if (!ra.col) ra.col = by_identity(av, stored_sel[i].fallback_identity);
            }
            axes.push_back(ra);
        }
        emit_rows(*os, ar, axes);
        ar.close();
        ++n_parts;
    }

    spdlog::info("report: profile '{}' over {} part(s), {} axis column(s): {}",
                 profile_name, n_parts, axis_order.size(), [&] {
                     std::string s; for (const auto& a : axis_order) { if (!s.empty()) s += ","; s += a; } return s;
                 }());
    return 0;
}

// ── profile authoring ─────────────────────────────────────────────────────────
namespace {

struct Selector {
    std::string axis;
    std::string tool, method;   bool cal = false;
    std::string f_tool, f_method; bool f_cal = false; bool has_fallback = false;
};

// Parse "tool:method[@cal]" → (tool, method, cal).
bool parse_tm(const std::string& s, std::string& tool, std::string& method, bool& cal) {
    auto colon = s.find(':');
    if (colon == std::string::npos) return false;
    tool = s.substr(0, colon);
    std::string m = s.substr(colon + 1);
    cal = false;
    const std::string at = "@cal";
    if (m.size() > at.size() && m.compare(m.size() - at.size(), at.size(), at) == 0) {
        cal = true; m.erase(m.size() - at.size());
    }
    method = m;
    return !tool.empty() && !method.empty();
}

// Parse "axis=tool:method[@cal][/tool2:method2[@cal]]".
bool parse_selector(const std::string& spec, Selector& out) {
    auto eq = spec.find('=');
    if (eq == std::string::npos) return false;
    out.axis = spec.substr(0, eq);
    std::string rhs = spec.substr(eq + 1);
    auto slash = rhs.find('/');
    std::string primary = (slash == std::string::npos) ? rhs : rhs.substr(0, slash);
    if (!parse_tm(primary, out.tool, out.method, out.cal)) return false;
    if (slash != std::string::npos) {
        out.has_fallback = parse_tm(rhs.substr(slash + 1), out.f_tool, out.f_method, out.f_cal);
        if (!out.has_fallback) return false;
    }
    return !out.axis.empty();
}

// Resolve (axis,tool,method,cal) to an on-disk column identity; 0 if not present.
uint64_t resolve_identity(const std::vector<AvailColumn>& av, const std::string& axis,
                          const std::string& tool, const std::string& method, bool cal) {
    const AvailColumn* hit = nullptr; int n = 0;
    for (const auto& a : av) {
        if (a.key.axis == axis && a.key.tool == tool && a.key.method == method &&
            ((a.key.calibration_hash != 0) == cal)) { hit = &a; ++n; }
    }
    if (n > 1)
        spdlog::warn("profile: {}={}:{}{} matches {} columns; pinning the first",
                     axis, tool, method, cal ? "@cal" : "", n);
    return hit ? hit->identity : 0;
}

} // namespace

int cmd_profile_add(const std::filesystem::path& pack_path,
                    const std::string& name,
                    const std::string& description,
                    const std::vector<std::string>& selectors) {
    if (name.empty()) throw std::runtime_error("profile add: --name is required");
    {
        Policy tmp;
        if (parse_policy(name, tmp))
            throw std::runtime_error("profile add: '" + name +
                                     "' is a reserved built-in profile name");
    }
    if (selectors.empty()) throw std::runtime_error("profile add: at least one --select is required");

    std::vector<Selector> sels;
    for (const auto& s : selectors) {
        Selector sel;
        if (!parse_selector(s, sel))
            throw std::runtime_error("profile add: bad selector '" + s +
                                     "' (want axis=tool:method[@cal][/tool:method[@cal]])");
        sels.push_back(std::move(sel));
    }

    auto gpks = collect_gpk_paths(pack_path);
    if (gpks.empty()) throw std::runtime_error("profile add: no .gpk found at " + pack_path.string());

    for (const auto& gp : gpks) {
        std::vector<Profile> profiles;
        Profile newp;
        newp.name = name; newp.description = description;
        size_t resolved = 0;
        {
            ArchiveReader ar; ar.open(gp);
            std::vector<AvailColumn> av; gather_columns(ar, av);

            for (const auto& sel : sels) {
                ProfileSelection ps;
                ps.axis            = sel.axis;
                ps.column_identity = resolve_identity(av, sel.axis, sel.tool, sel.method, sel.cal);
                if (ps.column_identity == 0)
                    spdlog::warn("profile add: {}: {}={}:{}{} not present in {} (pinned as 0)",
                                 name, sel.axis, sel.tool, sel.method, sel.cal ? "@cal" : "",
                                 gp.filename().string());
                else ++resolved;
                if (sel.has_fallback)
                    ps.fallback_identity = resolve_identity(av, sel.axis, sel.f_tool, sel.f_method, sel.f_cal);
                newp.selections.push_back(std::move(ps));
            }
            // Preserve existing stored profiles, replacing any with the same name.
            if (ar.has_prof()) {
                const ProfReader* pr = ar.prof_reader();
                for (uint32_t i = 0; i < pr->n_profiles(); ++i) {
                    Profile ex = pr->profile(i);
                    if (ex.name != name) profiles.push_back(std::move(ex));
                }
            }
            ar.close();
        }
        profiles.push_back(std::move(newp));
        write_prof_to_archive(gp, profiles);
        spdlog::info("profile add: '{}' ({} selections, {} resolved) -> {}",
                     name, sels.size(), resolved, gp.filename().string());
    }
    return 0;
}

int cmd_profile_list(const std::filesystem::path& pack_path) {
    auto gpks = collect_gpk_paths(pack_path);
    if (gpks.empty()) throw std::runtime_error("profile list: no .gpk found at " + pack_path.string());
    return do_list(gpks);
}

// ── qcontig viewer ────────────────────────────────────────────────────────────
int cmd_qcontig(const std::filesystem::path& pack_path,
                const std::string& accession_filter,
                const std::filesystem::path& output,
                float min_outlier, float min_leakage,
                float min_foreign, float min_lr) {
    auto gpks = collect_gpk_paths(pack_path);
    if (gpks.empty()) throw std::runtime_error("qcontig: no .gpk found at " + pack_path.string());

    std::ofstream fout;
    std::ostream* os = &std::cout;
    if (!output.empty()) {
        fout.open(output);
        if (!fout) throw std::runtime_error("qcontig: cannot write " + output.string());
        os = &fout;
    }
    *os << "accession\tcontig_index\tcontig_offset\tcontig_length"
           "\tgcov_t2_pct\tgcov_spe_pct\tprot_foreign_frac\tprot_loglr"
           "\tprot_foreign_specific\tprot_classifiable\tprot_flags\tchannel\tsuspicious\n";

    auto fmt = [](float v) { if (std::isnan(v)) return std::string("NA");
                             char b[32]; std::snprintf(b, sizeof(b), "%.6g", v); return std::string(b); };
    // TWO complementary channels. TNF/SPE (composition) catches LONG foreign contigs;
    // protein foreign-containment catches SMALL ones. A contig is suspicious if either
    // fires (or leakage). NaN/abstain never fires.
    auto tnf_sus = [&](const QcontigRecord& r) {
        return min_outlier > 0.0f &&
               ((!std::isnan(r.gcov_t2_pct)  && r.gcov_t2_pct  >= min_outlier) ||
                (!std::isnan(r.gcov_spe_pct) && r.gcov_spe_pct >= min_outlier));
    };
    auto foreign_sus = [&](const QcontigRecord& r) {
        return min_foreign > 0.0f && !(r.prot_flags & 0x1 /*ABSTAIN_LOW_N*/) &&
               !(r.prot_flags & 0x2 /*MOBILE_NATIVE*/) &&
               !std::isnan(r.prot_foreign_frac) && r.prot_foreign_frac >= min_foreign &&
               !std::isnan(r.prot_loglr) && r.prot_loglr >= min_lr;
    };
    auto leak_sus = [&](const QcontigRecord& r) {
        return min_leakage > 0.0f && !std::isnan(r.leakage_score) && r.leakage_score >= min_leakage;
    };
    const bool filtering = (min_outlier > 0.0f || min_leakage > 0.0f || min_foreign > 0.0f);

    size_t total = 0, flagged = 0;
    for (const auto& gp : gpks) {
        ArchiveReader ar; ar.open(gp);
        if (!ar.has_qcontig()) { ar.close(); continue; }
        const QcontigReader* qr = ar.qcontig_reader();
        ar.scan_genome_accessions([&](std::string_view acc, GenomeId gid) {
            if (ar.is_deleted(gid)) return;
            if (!accession_filter.empty() && acc != accession_filter) return;
            for (const auto& r : qr->contigs_for(static_cast<uint64_t>(gid))) {
                const bool t = tnf_sus(r), f = foreign_sus(r), l = leak_sus(r);
                const bool sus = t || f || l;
                if (sus) ++flagged;
                if (filtering && !sus) continue;   // restrict output to suspicious contigs
                std::string ch = t && f ? "both" : f ? "foreign" : t ? "tnf" : l ? "leakage" : "-";
                *os << acc << '\t' << r.contig_index << '\t' << r.contig_offset << '\t'
                    << r.contig_length << '\t' << fmt(r.gcov_t2_pct) << '\t' << fmt(r.gcov_spe_pct)
                    << '\t' << fmt(r.prot_foreign_frac) << '\t' << fmt(r.prot_loglr)
                    << '\t' << r.prot_foreign_specific << '\t' << r.prot_classifiable
                    << '\t' << static_cast<int>(r.prot_flags) << '\t' << ch << '\t' << (sus ? 1 : 0) << '\n';
                ++total;
            }
        });
        ar.close();
    }
    spdlog::info("qcontig: {} rows{}{}", total,
                 accession_filter.empty() ? "" : (" for " + accession_filter),
                 filtering ? (" (" + std::to_string(flagged) + " suspicious)") : "");
    return 0;
}

} // namespace genopack
