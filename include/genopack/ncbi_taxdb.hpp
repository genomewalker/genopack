#pragma once

#include <chrono>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace genopack {

namespace fs = std::filesystem;

// NCBI taxonomy database loader and 10-rank normaliser.
//
// Downloads new_taxdump.tar.gz from NCBI FTP on first use (or when stale).
// Parses nodes.dmp + names.dmp into a fast in-memory lookup.
// Converts NCBI taxonomy strings to the canonical 10-rank format:
//   d__ (superkingdom) | l__ (lineage/clade) | k__ (kingdom) |
//   p__ (phylum) | c__ (class) | o__ (order) | f__ (family) |
//   g__ (genus) | s__ (species) | S__ (subspecies)
//
// For GTDB prokaryotes (d__-prefixed), use the GTDB rank-propagation path.
// This class handles Eukaryota, Viruses, and other NCBI-only domains by
// matching names from the input taxonomy string against the NCBI tree.

class NcbiTaxdb {
public:
    static constexpr int kMaxAgeDays = 30;
    static constexpr const char* kDownloadUrl =
        "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz";

    // Load from an existing taxdump directory (must contain nodes.dmp + names.dmp).
    static NcbiTaxdb load(const fs::path& dir);

    // Ensure taxdump is present and not older than max_age_days.
    // Downloads and extracts new_taxdump.tar.gz if needed.
    static void ensure_fresh(const fs::path& dir, int max_age_days = kMaxAgeDays);

    // Return the timestamp of the current dump, or nullopt if not present.
    static std::optional<std::chrono::system_clock::time_point>
        dump_timestamp(const fs::path& dir);

    // Normalise a taxonomy string (e.g. "Eukaryota;Metazoa;Chordata;...") to
    // canonical 10-rank format. Tokens are matched from most-specific to
    // least-specific against the NCBI scientific-name index. Returns empty
    // string if no token could be resolved.
    // accession is used as the species leaf when s__ cannot be determined.
    std::string taxonomy_for_string(const std::string& tax_str,
                                    const std::string& accession) const;

    // Build a canonical 10-rank taxonomy string for a given NCBI taxid.
    // Returns empty string if taxid is not found.
    std::string taxonomy_for_taxid(int taxid, const std::string& accession) const;

    // Look up a taxid by scientific name. Returns -1 if not found or ambiguous.
    int taxid_for_name(const std::string& name) const;

    bool empty() const { return name_.empty(); }
    std::size_t size() const { return name_.size(); }

private:
    struct Node {
        int parent_taxid = 0;
        std::string rank;
    };

    std::unordered_map<int, Node>        nodes_;    // taxid → {parent, rank}
    std::unordered_map<int, std::string> name_;     // taxid → scientific name
    std::unordered_map<std::string, int> name2tid_; // scientific name → taxid (unique only)

    std::vector<std::pair<std::string, std::string>>
        lineage(int taxid) const;

    std::string lineage_to_10rank(
        const std::vector<std::pair<std::string, std::string>>& lineage,
        const std::string& accession) const;
};

} // namespace genopack
