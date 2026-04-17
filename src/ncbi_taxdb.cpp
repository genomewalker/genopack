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

std::vector<std::string> split_dmp(const std::string& line) {
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

NcbiTaxdb NcbiTaxdb::load(const fs::path& dir) {
    NcbiTaxdb db;

    {
        std::ifstream in(dir / "nodes.dmp");
        if (!in) throw std::runtime_error("Cannot open nodes.dmp in " + dir.string());
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            auto f = split_dmp(line);
            if (f.size() < 3) continue;
            int taxid  = std::stoi(f[0]);
            int parent = std::stoi(f[1]);
            db.nodes_[taxid] = {parent, f[2]};
        }
    }

    // Build both taxid→name and name→taxid (mark ambiguous names with -1).
    {
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
            if (!inserted) it->second = -1; // ambiguous — mark unusable
        }
    }
    // Remove ambiguous entries so taxid_for_name() never returns a wrong hit.
    std::erase_if(db.name2tid_, [](const auto& kv){ return kv.second == -1; });

    spdlog::info("NcbiTaxdb loaded: {} nodes, {} unique names", db.nodes_.size(), db.name2tid_.size());
    return db;
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
    std::unordered_set<int> visited;
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

} // namespace genopack
