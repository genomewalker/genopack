#include "genopack/ncbi_taxdb.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_set>
#include <thread>
#include <cstdint>

#include <spdlog/spdlog.h>

namespace genopack {

namespace {

static constexpr std::array<std::string_view, 10> kRanks = {
    "d__", "l__", "k__", "p__", "c__", "o__", "f__", "g__", "s__", "S__"
};

static const std::unordered_map<std::string, std::string> kNcbiToPrefix = {
    {"superkingdom", "d__"},
    {"kingdom",      "k__"},
    {"phylum",       "p__"},
    {"class",        "c__"},
    {"order",        "o__"},
    {"family",       "f__"},
    {"genus",        "g__"},
    {"species",      "s__"},
    {"subspecies",   "S__"},
};

// Clade names that serve as the l__ (lineage) slot for eukaryotes.
static const std::array<std::string_view, 8> kLineageClades = {
    "Opisthokonta",   // fungi, invertebrates, vertebrates
    "Viridiplantae",  // plants
    "SAR",            // stramenopiles, alveolates, rhizaria
    "Discoba",        // excavates
    "Amoebozoa",
    "Haptista",
    "Cryptista",
    "Rhodophyta",     // red algae (no kingdom rank in NCBI)
};

// Generic tokens that should not be used for taxid lookup.
static const std::unordered_set<std::string> kSkipTokens = {
    "unclassified", "environmental samples", "other sequences",
    "artificial sequences", "cellular organisms", "root",
    "metagenomes", "unidentified",
};

bool is_lineage_clade(const std::string& name) {
    for (auto c : kLineageClades)
        if (name == c) return true;
    return false;
}

std::string trim(const std::string& s) {
    auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return {};
    auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

std::vector<std::string> split_dmp(const std::string& line_in) {
    // NCBI .dmp records are "\t|\t"-separated and terminate with a "\t|" (or "|")
    // sentinel. Drop that terminator so the final field isn't left as "value\t|".
    std::string_view lv(line_in);
    if (lv.size() >= 2 && lv.substr(lv.size() - 2) == "\t|") lv.remove_suffix(2);
    else if (!lv.empty() && lv.back() == '|')                lv.remove_suffix(1);
    std::string line(lv);
    std::vector<std::string> fields;
    std::size_t pos = 0;
    while (true) {
        auto sep = line.find("\t|\t", pos);
        if (sep == std::string::npos) {
            fields.push_back(trim(line.substr(pos)));
            break;
        }
        fields.push_back(trim(line.substr(pos, sep - pos)));
        pos = sep + 3;
    }
    return fields;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Static utilities
// ---------------------------------------------------------------------------

std::optional<std::chrono::system_clock::time_point>
NcbiTaxdb::dump_timestamp(const fs::path& dir) {
    auto ts_file = dir / ".timestamp";
    std::ifstream in(ts_file);
    if (!in) return std::nullopt;
    std::time_t t = 0;
    in >> t;
    if (!in) return std::nullopt;
    return std::chrono::system_clock::from_time_t(t);
}

void NcbiTaxdb::ensure_fresh(const fs::path& dir, int max_age_days) {
    auto nodes = dir / "nodes.dmp";
    auto names = dir / "names.dmp";

    bool needs_download = !fs::exists(nodes) || !fs::exists(names);
    if (!needs_download) {
        auto opt_ts = dump_timestamp(dir);
        if (!opt_ts) {
            needs_download = true;
        } else {
            auto age = std::chrono::duration_cast<std::chrono::hours>(
                std::chrono::system_clock::now() - *opt_ts).count();
            needs_download = (age > max_age_days * 24);
        }
    }

    if (!needs_download) {
        spdlog::info("NCBI taxdump is up to date in {}", dir.string());
        return;
    }

    spdlog::info("Downloading NCBI taxdump → {}", dir.string());
    fs::create_directories(dir);

    auto archive = dir / "new_taxdump.tar.gz";
    std::string cmd = "curl -fsSL --retry 5 --retry-wait 10 -o " +
                      archive.string() + " " + kDownloadUrl;
    if (std::system(cmd.c_str()) != 0)
        throw std::runtime_error("Failed to download NCBI taxdump from " +
                                  std::string(kDownloadUrl));

    cmd = "tar -xzf " + archive.string() + " -C " + dir.string();
    if (std::system(cmd.c_str()) != 0)
        throw std::runtime_error("Failed to extract " + archive.string());

    fs::remove(archive);

    std::ofstream out(dir / ".timestamp");
    out << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    spdlog::info("NCBI taxdump downloaded and extracted ({})",
                 fs::exists(nodes) ? "ok" : "MISSING nodes.dmp");
}

// ---------------------------------------------------------------------------
// Load
// ---------------------------------------------------------------------------

namespace {
constexpr uint32_t kCacheMagic = 0x4E434258u;  // "NCBX"
constexpr uint32_t kCacheVer   = 1u;

template<class T> void wr(std::ostream& o, const T& v) { o.write(reinterpret_cast<const char*>(&v), sizeof(T)); }
template<class T> bool rd(std::istream& i, T& v) { return static_cast<bool>(i.read(reinterpret_cast<char*>(&v), sizeof(T))); }
void wrs(std::ostream& o, const std::string& s) { uint16_t n = (uint16_t)std::min<size_t>(s.size(), 65535); wr(o, n); o.write(s.data(), n); }
bool rds(std::istream& i, std::string& s) { uint16_t n; if (!rd(i, n)) return false; s.resize(n); return n == 0 || static_cast<bool>(i.read(&s[0], n)); }

// Cache is valid only if it exists and is newer than both source dumps.
bool cache_fresh(const fs::path& cache, const fs::path& dir) {
    std::error_code ec;
    if (!fs::exists(cache, ec)) return false;
    auto ct = fs::last_write_time(cache, ec); if (ec) return false;
    for (const char* f : {"nodes.dmp", "names.dmp"}) {
        auto st = fs::last_write_time(dir / f, ec); if (ec) return false;
        if (st > ct) return false;
    }
    return true;
}
}  // namespace

bool NcbiTaxdb::load_cache(const fs::path& cache) {
    std::ifstream in(cache, std::ios::binary);
    if (!in) return false;
    uint32_t magic, ver; if (!rd(in, magic) || !rd(in, ver) || magic != kCacheMagic || ver != kCacheVer) return false;
    uint64_t n;
    if (!rd(in, n)) return false;
    nodes_.reserve(n);
    for (uint64_t k = 0; k < n; ++k) { int32_t t, p; std::string r; if (!rd(in, t) || !rd(in, p) || !rds(in, r)) return false; nodes_[t] = {p, std::move(r)}; }
    if (!rd(in, n)) return false;
    name_.reserve(n);
    for (uint64_t k = 0; k < n; ++k) { int32_t t; std::string s; if (!rd(in, t) || !rds(in, s)) return false; name_[t] = std::move(s); }
    if (!rd(in, n)) return false;
    name2tid_.reserve(n);
    for (uint64_t k = 0; k < n; ++k) { std::string s; int32_t t; if (!rds(in, s) || !rd(in, t)) return false; name2tid_[std::move(s)] = t; }
    if (!rd(in, n)) return false;
    merged_.reserve(n);
    for (uint64_t k = 0; k < n; ++k) { int32_t a, b; if (!rd(in, a) || !rd(in, b)) return false; merged_[a] = b; }
    return static_cast<bool>(in);
}

void NcbiTaxdb::save_cache(const fs::path& cache) const {
    fs::path tmp = cache; tmp += ".tmp";
    { std::ofstream o(tmp, std::ios::binary); if (!o) return;
      wr(o, kCacheMagic); wr(o, kCacheVer);
      wr(o, (uint64_t)nodes_.size());     for (auto& [t, nd] : nodes_)   { wr(o, (int32_t)t); wr(o, (int32_t)nd.parent_taxid); wrs(o, nd.rank); }
      wr(o, (uint64_t)name_.size());      for (auto& [t, s]  : name_)    { wr(o, (int32_t)t); wrs(o, s); }
      wr(o, (uint64_t)name2tid_.size());  for (auto& [s, t]  : name2tid_){ wrs(o, s); wr(o, (int32_t)t); }
      wr(o, (uint64_t)merged_.size());    for (auto& [a, b]  : merged_)  { wr(o, (int32_t)a); wr(o, (int32_t)b); }
      if (!o) { std::error_code ec; fs::remove(tmp, ec); return; } }
    std::error_code ec; fs::rename(tmp, cache, ec);  // atomic publish; best-effort
}

NcbiTaxdb NcbiTaxdb::load(const fs::path& dir) {
    NcbiTaxdb db;

    // Fast path: a fresh binary cache skips the ~15-20s dmp parse. Any read error
    // falls through to a full parse, so a stale/corrupt cache never breaks a run.
    fs::path cache = dir / ".ncbi_taxdb.cache";
    if (cache_fresh(cache, dir) && db.load_cache(cache)) {
        spdlog::info("NcbiTaxdb loaded (cache): {} nodes, {} unique names, {} merged",
                     db.nodes_.size(), db.name2tid_.size(), db.merged_.size());
        return db;
    }
    db = NcbiTaxdb{};  // discard any partial cache read

    // Parse nodes.dmp and names.dmp concurrently (disjoint maps) + merged.dmp.
    auto parse_nodes = [&] {
        std::ifstream in(dir / "nodes.dmp");
        if (!in) throw std::runtime_error("Cannot open nodes.dmp in " + dir.string());
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            auto f = split_dmp(line);
            if (f.size() < 3) continue;
            db.nodes_[std::stoi(f[0])] = {std::stoi(f[1]), f[2]};
        }
    };
    auto parse_names = [&] {
        std::ifstream in(dir / "names.dmp");
        if (!in) throw std::runtime_error("Cannot open names.dmp in " + dir.string());
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            auto f = split_dmp(line);
            if (f.size() < 4 || f[3] != "scientific name") continue;
            int taxid = std::stoi(f[0]);
            const std::string& sci = f[1];
            db.name_[taxid] = sci;
            auto [it, inserted] = db.name2tid_.emplace(sci, taxid);
            if (!inserted) it->second = -1;  // ambiguous — mark unusable
        }
        std::erase_if(db.name2tid_, [](const auto& kv){ return kv.second == -1; });
    };
    std::exception_ptr err;
    std::thread tn([&]{ try { parse_names(); } catch (...) { err = std::current_exception(); } });
    parse_nodes();
    tn.join();
    if (err) std::rethrow_exception(err);

    {
        std::ifstream in(dir / "merged.dmp");
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            auto f = split_dmp(line);
            if (f.size() < 2) continue;
            try { db.merged_[std::stoi(f[0])] = std::stoi(f[1]); } catch (...) {}
        }
    }

    db.save_cache(cache);  // best-effort; next load is fast
    spdlog::info("NcbiTaxdb loaded: {} nodes, {} unique names, {} merged", db.nodes_.size(), db.name2tid_.size(), db.merged_.size());
    return db;
}

int NcbiTaxdb::resolve_merged(int taxid) const {
    // follow merged chain (bounded) to the current taxid
    for (int i = 0; i < 16; ++i) {
        if (nodes_.count(taxid)) return taxid;
        auto it = merged_.find(taxid);
        if (it == merged_.end()) return taxid;
        taxid = it->second;
    }
    return taxid;
}

std::vector<std::pair<std::string, uint64_t>>
NcbiTaxdb::canonical_path_taxids(int taxid) const {
    static const std::unordered_map<std::string, std::string> kPfx = {
        {"domain","d__"},{"superkingdom","d__"},{"realm","d__"},{"kingdom","k__"},
        {"phylum","p__"},{"class","c__"},{"order","o__"},{"family","f__"},
        {"genus","g__"},{"species","s__"}};
    taxid = resolve_merged(taxid);
    // collect ancestors root->taxon with (rank, name, taxid)
    std::vector<std::tuple<std::string,std::string,int>> chain;
    std::unordered_set<int> seen;
    int cur = taxid;
    while (cur != 1 && cur != 0 && !seen.count(cur)) {
        seen.insert(cur);
        auto nit = nodes_.find(cur);
        if (nit == nodes_.end()) break;
        auto nme = name_.find(cur);
        chain.emplace_back(nit->second.rank, nme != name_.end() ? nme->second : "", cur);
        cur = nit->second.parent_taxid;
    }
    std::reverse(chain.begin(), chain.end());
    std::vector<std::pair<std::string, uint64_t>> out;
    std::string path;
    int deepest = 0;
    for (auto& [rank, name, tid] : chain) {
        auto it = kPfx.find(rank);
        if (it == kPfx.end() || name.empty()) continue;
        std::string tok = it->second + name;
        path = path.empty() ? tok : path + ";" + tok;
        out.emplace_back(path, static_cast<uint64_t>(tid));
        deepest = tid;
    }
    // taxid below the deepest canonical rank (strain/subspecies): keep it as a leaf
    // under the species so the accession preserves its authoritative id.
    if (!out.empty() && deepest != taxid) {
        auto nme = name_.find(taxid);
        std::string sname = nme != name_.end() ? nme->second : ("taxon " + std::to_string(taxid));
        path += ";s__" + sname;
        out.emplace_back(path, static_cast<uint64_t>(taxid));
    }
    return out;
}

// ---------------------------------------------------------------------------
// Name lookup
// ---------------------------------------------------------------------------

int NcbiTaxdb::taxid_for_name(const std::string& name) const {
    auto it = name2tid_.find(name);
    return (it != name2tid_.end()) ? it->second : -1;
}

// ---------------------------------------------------------------------------
// Taxonomy string resolution
// ---------------------------------------------------------------------------

std::string NcbiTaxdb::taxonomy_for_string(const std::string& tax_str,
                                            const std::string& accession) const {
    // Split on ';', try tokens from most-specific (rightmost) to least-specific.
    // Skip known generic tokens. Return the first successful taxid resolution.
    std::vector<std::string> tokens;
    std::string_view sv(tax_str);
    while (!sv.empty()) {
        auto sep = sv.find(';');
        auto tok = trim(std::string(sv.substr(0, sep)));
        if (!tok.empty()) tokens.push_back(std::move(tok));
        sv = (sep == std::string_view::npos) ? "" : sv.substr(sep + 1);
    }

    for (auto it = tokens.rbegin(); it != tokens.rend(); ++it) {
        const std::string& tok = *it;
        // Skip prefixed GTDB-style tokens (caller should use GTDB path for d__)
        if (tok.size() >= 3 && tok[1] == '_' && tok[2] == '_') continue;
        // Skip generic noise tokens
        std::string lower = tok;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        if (kSkipTokens.count(lower)) continue;

        int taxid = taxid_for_name(tok);
        if (taxid > 0) {
            std::string result = taxonomy_for_taxid(taxid, accession);
            if (!result.empty()) return result;
        }
    }
    return {};
}

// ---------------------------------------------------------------------------
// Lineage builder
// ---------------------------------------------------------------------------

std::vector<std::pair<std::string, std::string>>
NcbiTaxdb::lineage(int taxid) const {
    std::vector<std::pair<std::string, std::string>> result;
    static std::unordered_set<int> visited;
    visited.clear();
    int cur = taxid;
    while (cur != 1 && cur != 0 && !visited.count(cur)) {
        visited.insert(cur);
        auto nit = nodes_.find(cur);
        if (nit == nodes_.end()) break;
        auto nme = name_.find(cur);
        result.push_back({nit->second.rank,
                          nme != name_.end() ? nme->second : ""});
        cur = nit->second.parent_taxid;
    }
    std::reverse(result.begin(), result.end());
    return result;
}

std::string NcbiTaxdb::lineage_to_10rank(
    const std::vector<std::pair<std::string, std::string>>& lin,
    const std::string& accession) const {

    std::unordered_map<std::string, std::string> rank_map;

    for (const auto& [ncbi_rank, name] : lin) {
        auto it = kNcbiToPrefix.find(ncbi_rank);
        if (it != kNcbiToPrefix.end()) {
            rank_map[it->second] = it->second + name;
        }
        if (!rank_map.count("l__") && is_lineage_clade(name)) {
            rank_map["l__"] = "l__" + name;
        }
    }

    auto propagate = [&](const std::string& child, const std::string& parent) {
        if (!rank_map.count(child) || rank_map[child] == child) {
            if (rank_map.count(parent))
                rank_map[child] = child + rank_map[parent].substr(3);
        }
    };
    propagate("l__", "d__");
    propagate("k__", "l__");
    propagate("p__", "k__");
    propagate("c__", "p__");
    propagate("o__", "c__");
    propagate("f__", "o__");
    propagate("g__", "f__");

    if (!rank_map.count("s__") || rank_map["s__"] == "s__") {
        if (rank_map.count("g__") && rank_map["g__"] != "g__")
            rank_map["s__"] = "s__" + rank_map["g__"].substr(3) + " unclassified";
        else {
            // No genus resolved — stamp accession at every rank below d__ so
            // each genome is its own taxon with no shared nodes
            for (const char* r : {"l__","k__","p__","c__","o__","f__","g__","s__"})
                rank_map[r] = std::string(r) + accession;
        }
    }

    rank_map["S__"] = "S__" + rank_map["s__"].substr(3);

    std::string result;
    result.reserve(200);
    for (auto pfx : kRanks) {
        if (!result.empty()) result += ';';
        auto it = rank_map.find(std::string(pfx));
        result += (it != rank_map.end()) ? it->second : std::string(pfx);
    }
    return result;
}

std::string NcbiTaxdb::taxonomy_for_taxid(int taxid,
                                           const std::string& accession) const {
    if (taxid <= 0) return {};
    auto lin = lineage(taxid);
    if (lin.empty()) return {};
    return lineage_to_10rank(lin, accession);
}

void NcbiTaxdb::for_each_node(
    const std::function<void(int, int, const std::string&, const std::string&)>& cb) const {
    static const std::string kEmpty;
    for (const auto& [taxid, nd] : nodes_) {
        auto nm = name_.find(taxid);
        cb(taxid, nd.parent_taxid, nd.rank, nm != name_.end() ? nm->second : kEmpty);
    }
}

std::unordered_set<int> NcbiTaxdb::descendants_including(const std::vector<int>& roots) const {
    std::unordered_map<int, std::vector<int>> children;
    children.reserve(nodes_.size());
    for (const auto& [taxid, nd] : nodes_)
        if (nd.parent_taxid != taxid) children[nd.parent_taxid].push_back(taxid);
    std::unordered_set<int> out;
    std::vector<int> stack;
    for (int r : roots) { if (out.insert(r).second) stack.push_back(r); }
    while (!stack.empty()) {
        int t = stack.back(); stack.pop_back();
        auto it = children.find(t);
        if (it == children.end()) continue;
        for (int c : it->second) if (out.insert(c).second) stack.push_back(c);
    }
    return out;
}

} // namespace genopack
