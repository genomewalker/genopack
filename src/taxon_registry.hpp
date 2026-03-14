#pragma once
#include <genopack/types.hpp>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace genopack {

// Builds the TaxonId registry from GTDB lineage strings.
// Thread-safe for concurrent lookups after build.
class TaxonRegistry {
public:
    TaxonRegistry();

    // Parse GTDB string "d__Bacteria;…;s__Escherichia coli" and register all
    // nodes. Returns the species TaxonId (leaf node).
    TaxonId insert_lineage(std::string_view lineage);

    // Lookup by name at any rank ("Escherichia coli", "g__Salmonella", etc.)
    TaxonId find_by_name(std::string_view name) const;

    // Lookup by full lineage string
    TaxonId find_by_lineage(std::string_view lineage) const;

    // Subtree range [lo, hi) for a parent taxon_id (genus → all species)
    std::pair<TaxonId, TaxonId> subtree_range(TaxonId parent) const;

    // Info access
    const TaxonInfo& info(TaxonId id) const;
    size_t size() const { return nodes_.size(); }

    // Export sorted node list (for writing taxon.idx)
    const std::vector<TaxonInfo>& sorted_nodes() const;

    // Assign stable DFS order for subtree range queries
    void finalize();

private:
    // Tries: "d__Bacteria" → id, "g__Salmonella" → id, "s__..." → id
    std::unordered_map<std::string, TaxonId> by_qualified_name_;
    std::vector<TaxonInfo>                   nodes_;
    std::vector<TaxonId>                     children_[256];  // parent → children
    bool                                     finalized_ = false;

    TaxonId next_id_ = 1;  // 0 = INVALID

    TaxonId get_or_create(std::string_view qualified_name,
                          TaxonId parent_id, TaxonRank rank);
    void dfs_assign(TaxonId id, TaxonId& counter,
                    std::vector<TaxonInfo>& ordered);
};

// Parse GTDB lineage into rank-qualified tokens:
// "d__Bacteria;l__Bacteria;…;s__Escherichia coli"
// → [{"d__Bacteria", Domain}, {"l__Bacteria", Life}, …]
struct RankedToken { std::string token; TaxonRank rank; };
std::vector<RankedToken> parse_gtdb_lineage(std::string_view lineage);

} // namespace genopack
