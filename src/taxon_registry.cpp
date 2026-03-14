#include "taxon_registry.hpp"
#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace genopack {

static TaxonRank rank_from_prefix(char p) {
    switch (p) {
        case 'd': return TaxonRank::Domain;
        case 'l': return TaxonRank::Life;
        case 'k': return TaxonRank::Kingdom;
        case 'p': return TaxonRank::Phylum;
        case 'c': return TaxonRank::Class;
        case 'o': return TaxonRank::Order;
        case 'f': return TaxonRank::Family;
        case 'g': return TaxonRank::Genus;
        case 's': return TaxonRank::Species;
        default:  return TaxonRank::Domain;
    }
}

std::vector<RankedToken> parse_gtdb_lineage(std::string_view lineage) {
    std::vector<RankedToken> tokens;
    size_t start = 0;
    while (start < lineage.size()) {
        size_t end = lineage.find(';', start);
        if (end == std::string_view::npos) end = lineage.size();
        std::string_view tok = lineage.substr(start, end - start);
        if (tok.size() >= 3 && tok[1] == '_' && tok[2] == '_')
            tokens.push_back({std::string(tok), rank_from_prefix(tok[0])});
        start = end + 1;
    }
    return tokens;
}

TaxonRegistry::TaxonRegistry() {
    // Reserve ID 0 = INVALID, start at 1
    nodes_.push_back({INVALID_TAXON_ID, INVALID_TAXON_ID, 0, "", ""});
}

TaxonId TaxonRegistry::get_or_create(std::string_view qualified_name,
                                      TaxonId parent_id, TaxonRank rank) {
    std::string key(qualified_name);
    auto it = by_qualified_name_.find(key);
    if (it != by_qualified_name_.end()) return it->second;

    TaxonId id = next_id_++;
    by_qualified_name_[key] = id;

    TaxonInfo info;
    info.taxon_id  = id;
    info.parent_id = parent_id;
    info.rank      = static_cast<uint8_t>(rank);
    info.name      = std::string(qualified_name);
    nodes_.push_back(info);
    return id;
}

TaxonId TaxonRegistry::insert_lineage(std::string_view lineage) {
    if (finalized_)
        throw std::logic_error("TaxonRegistry: cannot insert after finalize()");

    auto tokens = parse_gtdb_lineage(lineage);
    if (tokens.empty()) return INVALID_TAXON_ID;

    TaxonId parent = INVALID_TAXON_ID;
    std::string cumulative;
    TaxonId leaf = INVALID_TAXON_ID;

    for (const auto& [tok, rank] : tokens) {
        if (!cumulative.empty()) cumulative += ';';
        cumulative += tok;
        TaxonId id = get_or_create(tok, parent, rank);
        nodes_[id].full_lineage = cumulative;
        parent = id;
        leaf   = id;
    }
    return leaf;
}

TaxonId TaxonRegistry::find_by_name(std::string_view name) const {
    auto it = by_qualified_name_.find(std::string(name));
    return (it != by_qualified_name_.end()) ? it->second : INVALID_TAXON_ID;
}

TaxonId TaxonRegistry::find_by_lineage(std::string_view lineage) const {
    // Try full lineage (last token is the key)
    auto tokens = parse_gtdb_lineage(lineage);
    if (tokens.empty()) return INVALID_TAXON_ID;
    return find_by_name(tokens.back().token);
}

const TaxonInfo& TaxonRegistry::info(TaxonId id) const {
    if (id >= nodes_.size())
        throw std::out_of_range("TaxonRegistry::info: invalid id");
    return nodes_[id];
}

const std::vector<TaxonInfo>& TaxonRegistry::sorted_nodes() const {
    return nodes_;
}

void TaxonRegistry::finalize() {
    // Sort nodes by full_lineage for stable DFS/taxonomy order.
    // After sort, node positions change — rebuild by_qualified_name_ accordingly.
    // We assign contiguous DFS-ordered IDs so that subtree_range is a simple [lo, hi).

    // Collect (lineage, old_id) pairs
    std::vector<std::pair<std::string, TaxonId>> sorted;
    sorted.reserve(nodes_.size() - 1);
    for (TaxonId id = 1; id < static_cast<TaxonId>(nodes_.size()); ++id)
        sorted.emplace_back(nodes_[id].full_lineage, id);

    std::sort(sorted.begin(), sorted.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    // Remap old IDs → new IDs (sorted order)
    std::vector<TaxonId> old_to_new(nodes_.size(), INVALID_TAXON_ID);
    std::vector<TaxonInfo> new_nodes;
    new_nodes.push_back(nodes_[0]);  // preserve INVALID at index 0

    for (TaxonId new_id = 1; const auto& [lineage, old_id] : sorted) {
        old_to_new[old_id] = new_id;
        TaxonInfo n = nodes_[old_id];
        n.taxon_id = new_id;
        new_nodes.push_back(std::move(n));
        ++new_id;
    }

    // Fix parent_ids
    for (auto& n : new_nodes) {
        if (n.parent_id != INVALID_TAXON_ID)
            n.parent_id = old_to_new[n.parent_id];
    }

    // Rebuild index
    by_qualified_name_.clear();
    for (const auto& n : new_nodes) {
        if (n.taxon_id != INVALID_TAXON_ID)
            by_qualified_name_[n.name] = n.taxon_id;
    }

    nodes_ = std::move(new_nodes);
    next_id_ = static_cast<TaxonId>(nodes_.size());
    finalized_ = true;
}

std::pair<TaxonId, TaxonId> TaxonRegistry::subtree_range(TaxonId parent) const {
    if (!finalized_)
        throw std::logic_error("TaxonRegistry: call finalize() first");
    if (parent == INVALID_TAXON_ID || parent >= static_cast<TaxonId>(nodes_.size()))
        return {1, static_cast<TaxonId>(nodes_.size())};

    // Since nodes are sorted by full_lineage, all children have lineages that
    // extend the parent's lineage (prefix match). Find the range.
    const std::string& prefix = nodes_[parent].full_lineage;
    TaxonId lo = parent;
    TaxonId hi = parent + 1;
    for (TaxonId id = parent + 1; id < static_cast<TaxonId>(nodes_.size()); ++id) {
        const std::string& lin = nodes_[id].full_lineage;
        if (lin.size() > prefix.size() &&
            lin.substr(0, prefix.size()) == prefix &&
            lin[prefix.size()] == ';') {
            hi = id + 1;
        } else {
            break;
        }
    }
    return {lo, hi};
}

} // namespace genopack
