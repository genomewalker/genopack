#include <genopack/txdb.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace genopack {

// ── rank helpers ─────────────────────────────────────────────────────────────

static TaxRank prefix_to_rank(std::string_view prefix) {
    if (prefix == "d") return TaxRank::DOMAIN;
    if (prefix == "l") return TaxRank::LIFE;
    if (prefix == "k") return TaxRank::KINGDOM;
    if (prefix == "p") return TaxRank::PHYLUM;
    if (prefix == "c") return TaxRank::CLASS;
    if (prefix == "o") return TaxRank::ORDER;
    if (prefix == "f") return TaxRank::FAMILY;
    if (prefix == "g") return TaxRank::GENUS;
    if (prefix == "s") return TaxRank::SPECIES;
    return TaxRank::NO_RANK;
}

static std::string_view rank_prefix(TaxRank r) {
    switch (r) {
        case TaxRank::DOMAIN:  return "d";
        case TaxRank::LIFE:    return "l";
        case TaxRank::KINGDOM: return "k";
        case TaxRank::PHYLUM:  return "p";
        case TaxRank::CLASS:   return "c";
        case TaxRank::ORDER:   return "o";
        case TaxRank::FAMILY:  return "f";
        case TaxRank::GENUS:   return "g";
        case TaxRank::SPECIES: return "s";
        default:               return "";
    }
}

static int rank_depth(TaxRank r) {
    return static_cast<int>(static_cast<uint8_t>(r));
}

static const TaxRank STD_FILL_ORDER[] = {
    TaxRank::DOMAIN, TaxRank::PHYLUM, TaxRank::CLASS, TaxRank::ORDER,
    TaxRank::FAMILY, TaxRank::GENUS,  TaxRank::SPECIES
};
static const TaxRank FULL_FILL_ORDER[] = {
    TaxRank::DOMAIN, TaxRank::LIFE, TaxRank::KINGDOM, TaxRank::PHYLUM,
    TaxRank::CLASS,  TaxRank::ORDER, TaxRank::FAMILY,  TaxRank::GENUS, TaxRank::SPECIES
};

// ── lineage parsing ──────────────────────────────────────────────────────────

struct ParsedRank {
    TaxRank     rank;
    std::string name;
};

static std::vector<ParsedRank> parse_lineage(std::string_view lineage) {
    std::vector<ParsedRank> result;
    size_t pos = 0;
    while (pos < lineage.size()) {
        size_t sep = lineage.find(';', pos);
        std::string_view token = (sep == std::string_view::npos)
            ? lineage.substr(pos)
            : lineage.substr(pos, sep - pos);

        size_t dunder = token.find("__");
        if (dunder != std::string_view::npos) {
            std::string_view prefix = token.substr(0, dunder);
            std::string_view name   = token.substr(dunder + 2);
            TaxRank r = prefix_to_rank(prefix);
            if (r != TaxRank::NO_RANK)
                result.push_back({r, std::string(name)});
        }

        pos = (sep == std::string_view::npos) ? lineage.size() : sep + 1;
    }
    return result;
}

// ── hash helpers ──────────────────────────────────────────────────────────────

uint32_t TaxonomyTree::fnv1a(std::string_view s) {
    uint32_t h = 2166136261u;
    for (unsigned char c : s) h = (h ^ c) * 16777619u;
    return h;
}

// FNV-1a 64-bit hash of canonical path, with GTDB namespace bit set.
// Returns a uint64_t with bit 63 = 1 (GTDB synthetic concept_id).
static uint64_t path_to_concept_id(std::string_view full_path) {
    uint64_t h = 14695981039346656037ULL;
    for (unsigned char c : full_path) h = (h ^ c) * 1099511628211ULL;
    h |= TAXID_GTDB_BIT;
    // Avoid the degenerate value where only the namespace bit is set
    if ((h & ~TAXID_GTDB_BIT) == 0) h |= 1;
    return h;
}

uint64_t TaxonomyTree::concept_id_for_path(std::string_view path) {
    return path_to_concept_id(path);
}

// ── build_internal ─────────────────────────────────────────────────────────--

struct NodeInfo {
    uint64_t    taxid;
    uint64_t    parent_taxid;
    TaxRank     rank;
    uint8_t     flags;      // bit 0: synthetic, bit 1: per_accession
    std::string name;
    std::string full_path;
};

struct BuiltTree {
    std::vector<NodeInfo>                         nodes_ordered;
    std::vector<std::pair<std::string, uint64_t>> acc_taxids;
};

static float dot136(const float* a, const float* b) {
    float s = 0.0f;
    for (int i = 0; i < 136; ++i) s += a[i] * b[i];
    return s;
}

static std::string cluster_group_name(const std::string& anchor,
                                       const std::string& rep_acc) {
    uint32_t h = TaxonomyTree::fnv1a(rep_acc);
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%08x", h);
    return "s__[" + anchor + " sp. group" + buf + "]";
}

static BuiltTree build_internal(
    const std::vector<std::pair<std::string, std::string>>& entries,
    const std::unordered_map<std::string, std::array<float, 136>>* kmer_profiles = nullptr,
    float similarity_threshold = 0.97f)
{
    std::unordered_map<std::string, std::string> acc_to_cluster_rep;

    if (kmer_profiles && !entries.empty()) {
        std::unordered_map<std::string, std::string> acc_to_anchor;
        for (const auto& [acc, lineage] : entries) {
            auto parsed = parse_lineage(lineage);
            std::string cur_path;
            std::string anchor_path;
            for (auto& pr : parsed) {
                if (pr.name.empty()) continue;
                std::string token = std::string(rank_prefix(pr.rank)) + "__" + pr.name;
                cur_path = cur_path.empty() ? token : cur_path + ";" + token;
                anchor_path = cur_path;
            }
            acc_to_anchor[acc] = anchor_path;
        }

        std::unordered_map<std::string, std::vector<std::string>> anchor_to_accs;
        for (const auto& [acc, anchor] : acc_to_anchor)
            anchor_to_accs[anchor].push_back(acc);

        for (auto& [anchor, accs] : anchor_to_accs) {
            std::vector<int> rep_idx(accs.size(), -1);
            std::vector<size_t> reps;

            for (size_t i = 0; i < accs.size(); ++i) {
                auto it_i = kmer_profiles->find(accs[i]);
                if (it_i == kmer_profiles->end()) {
                    rep_idx[i] = static_cast<int>(i);
                    reps.push_back(i);
                    continue;
                }
                const float* pi = it_i->second.data();
                bool found = false;
                for (size_t r : reps) {
                    auto it_r = kmer_profiles->find(accs[r]);
                    if (it_r == kmer_profiles->end()) continue;
                    if (dot136(pi, it_r->second.data()) >= similarity_threshold) {
                        rep_idx[i] = static_cast<int>(r);
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    rep_idx[i] = static_cast<int>(i);
                    reps.push_back(i);
                }
            }

            for (size_t i = 0; i < accs.size(); ++i) {
                int r = rep_idx[i];
                acc_to_cluster_rep[accs[i]] = accs[static_cast<size_t>(r)];
            }
        }
    }

    std::unordered_map<std::string, uint64_t> path_to_tid;
    std::unordered_map<uint64_t, std::string> tid_to_path;
    std::unordered_map<std::string, NodeInfo> path_to_node;
    std::vector<std::string> insertion_order;

    bool has_life = false, has_kingdom = false;
    for (const auto& [acc, lin] : entries) {
        auto parsed = parse_lineage(lin);
        for (auto& p : parsed) {
            if (p.rank == TaxRank::LIFE)    has_life    = true;
            if (p.rank == TaxRank::KINGDOM) has_kingdom = true;
        }
    }

    std::vector<TaxRank> fill_order;
    if (has_life || has_kingdom) {
        for (auto r : FULL_FILL_ORDER) fill_order.push_back(r);
    } else {
        for (auto r : STD_FILL_ORDER)  fill_order.push_back(r);
    }

    // Assign collision-resolved concept_id for a canonical path.
    // Root (path="root") gets taxid=1 without the GTDB bit.
    // All other nodes get FNV-64 with GTDB bit set.
    auto assign_taxid = [&](const std::string& path) -> uint64_t {
        auto it = path_to_tid.find(path);
        if (it != path_to_tid.end()) return it->second;

        uint64_t tid = path_to_concept_id(path);
        // Resolve collisions: increment within GTDB namespace
        while (true) {
            auto existing = tid_to_path.find(tid);
            if (existing == tid_to_path.end()) break;
            if (existing->second == path) break;
            ++tid;
            // Keep bit 63 set and avoid wrapping to 0 in low bits
            tid |= TAXID_GTDB_BIT;
            if ((tid & ~TAXID_GTDB_BIT) == 0) tid |= 1;
        }
        path_to_tid[path] = tid;
        tid_to_path[tid]  = path;
        return tid;
    };

    // Root node: taxid=1 (NCBI-style common root, no GTDB bit)
    {
        const std::string root_path = "root";
        uint64_t root_tid = 1;
        path_to_tid[root_path] = root_tid;
        tid_to_path[root_tid]  = root_path;
        if (path_to_node.find(root_path) == path_to_node.end()) {
            path_to_node[root_path] = {root_tid, root_tid, TaxRank::NO_RANK, 0, "root", root_path};
            insertion_order.push_back(root_path);
        }
    }

    std::vector<std::pair<std::string, uint64_t>> acc_taxids;
    acc_taxids.reserve(entries.size());

    for (const auto& [acc, lineage] : entries) {
        auto parsed = parse_lineage(lineage);

        int deepest_idx = -1;
        for (int i = static_cast<int>(parsed.size()) - 1; i >= 0; --i) {
            if (!parsed[i].name.empty()) { deepest_idx = i; break; }
        }

        std::string cur_path;
        uint64_t    prev_tid     = 1;
        std::string anchor_name;
        std::string anchor_path;
        TaxRank     deepest_rank = TaxRank::DOMAIN;

        for (auto& pr : parsed) {
            if (pr.name.empty()) continue;

            std::string token = std::string(rank_prefix(pr.rank)) + "__" + pr.name;
            if (cur_path.empty()) cur_path = token;
            else                  cur_path += ";" + token;

            uint64_t tid = assign_taxid(cur_path);

            if (path_to_node.find(cur_path) == path_to_node.end()) {
                path_to_node[cur_path] = {tid, prev_tid, pr.rank, 0, pr.name, cur_path};
                insertion_order.push_back(cur_path);
            }

            prev_tid     = tid;
            anchor_name  = pr.name;
            anchor_path  = cur_path;
            deepest_rank = pr.rank;
        }

        (void)deepest_idx;

        for (TaxRank fill : fill_order) {
            if (rank_depth(fill) <= rank_depth(deepest_rank)) continue;

            bool is_species = (fill == TaxRank::SPECIES);
            std::string syn_name, syn_path;

            if (is_species) {
                auto cluster_it = acc_to_cluster_rep.find(acc);
                if (cluster_it != acc_to_cluster_rep.end()) {
                    syn_name = cluster_group_name(anchor_name, cluster_it->second);
                } else {
                    syn_name = "s__[" + anchor_name + " sp. " + acc + "]";
                }
                syn_path = anchor_path + ";" + syn_name;
            } else {
                std::string pfx = std::string(rank_prefix(fill));
                syn_name = pfx + "__[unclassified " + anchor_name + "]";
                syn_path = anchor_path + ";" + syn_name;
            }

            uint64_t syn_tid = assign_taxid(syn_path);

            if (path_to_node.find(syn_path) == path_to_node.end()) {
                uint8_t syn_flags = is_species ? 0x02u : 0x01u;
                path_to_node[syn_path] = {syn_tid, prev_tid, fill, syn_flags, syn_name, syn_path};
                insertion_order.push_back(syn_path);
            }

            prev_tid = syn_tid;
        }

        acc_taxids.emplace_back(acc, prev_tid);
    }

    BuiltTree bt;
    bt.nodes_ordered.reserve(insertion_order.size());
    for (const auto& path : insertion_order)
        bt.nodes_ordered.push_back(std::move(path_to_node.at(path)));
    bt.acc_taxids = std::move(acc_taxids);
    return bt;
}

// ── TxdbWriter ────────────────────────────────────────────────────────────────

void TxdbWriter::add(std::string_view accession, std::string_view lineage_string) {
    entries_.emplace_back(std::string(accession), std::string(lineage_string));
}

void TxdbWriter::set_kmer_profiles(
    const std::unordered_map<std::string, std::array<float, 136>>& profiles,
    float similarity_threshold)
{
    kmer_profiles_        = &profiles;
    similarity_threshold_ = similarity_threshold;
}

TaxonomyTree TxdbWriter::build_tree() const {
    BuiltTree bt = build_internal(entries_, kmer_profiles_, similarity_threshold_);

    std::sort(bt.nodes_ordered.begin(), bt.nodes_ordered.end(),
              [](const NodeInfo& a, const NodeInfo& b) { return a.taxid < b.taxid; });

    std::string name_pool;
    std::vector<uint32_t> name_offsets;
    name_offsets.reserve(bt.nodes_ordered.size());
    for (const auto& n : bt.nodes_ordered) {
        name_offsets.push_back(static_cast<uint32_t>(name_pool.size()));
        name_pool.append(n.name);
        name_pool.push_back('\0');
    }

    std::vector<TaxNodeRecord> records;
    records.reserve(bt.nodes_ordered.size());
    for (size_t i = 0; i < bt.nodes_ordered.size(); ++i) {
        const auto& ni = bt.nodes_ordered[i];
        TaxNodeRecord rec{};
        rec.taxid        = ni.taxid;
        rec.parent_taxid = ni.parent_taxid;
        rec.name_offset  = name_offsets[i];
        rec.name_len     = static_cast<uint16_t>(ni.name.size());
        rec.rank         = static_cast<uint8_t>(ni.rank);
        rec.taxnamespace = static_cast<uint8_t>(
            taxid_is_gtdb(ni.taxid) ? TaxNamespace::GTDB : TaxNamespace::NCBI);
        rec.status       = static_cast<uint8_t>(TaxNodeStatus::ACTIVE);
        rec.flags        = ni.flags;
        records.push_back(rec);
    }

    auto& at = bt.acc_taxids;
    std::sort(at.begin(), at.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    uint32_t n_acc = static_cast<uint32_t>(at.size());

    std::vector<std::string> acc_sorted;
    std::vector<uint64_t>    acc_taxids_vec;
    acc_sorted.reserve(n_acc);
    acc_taxids_vec.reserve(n_acc);
    for (const auto& [acc, tid] : at) {
        acc_sorted.push_back(acc);
        acc_taxids_vec.push_back(tid);
    }

    uint32_t n_buckets = 1;
    while (n_buckets < 2 * n_acc) n_buckets <<= 1;
    std::vector<uint32_t> buckets(n_buckets, UINT32_MAX);
    for (uint32_t i = 0; i < n_acc; ++i) {
        uint32_t h = TaxonomyTree::fnv1a(acc_sorted[i]) & (n_buckets - 1);
        while (buckets[h] != UINT32_MAX) h = (h + 1) & (n_buckets - 1);
        buckets[h] = i;
    }

    TaxonomyTree::OwnedData od;
    od.nodes      = std::move(records);
    od.name_pool  = std::move(name_pool);
    od.acc_sorted   = std::move(acc_sorted);
    od.acc_taxids   = std::move(acc_taxids_vec);
    od.acc_buckets  = std::move(buckets);
    od.n_buckets    = n_buckets;

    return TaxonomyTree(std::move(od));
}

SectionDesc TxdbWriter::finalize(AppendWriter& writer, uint64_t section_id) {
    TaxonomyTree tree = build_tree();
    const TaxonomyTree::OwnedData& od = *tree.owned_data();

    uint32_t n_nodes = static_cast<uint32_t>(od.nodes.size());
    uint32_t n_acc   = static_cast<uint32_t>(od.acc_sorted.size());
    uint32_t n_buck  = od.n_buckets;

    std::string acc_str_pool;
    std::vector<uint32_t> acc_offsets;
    acc_offsets.reserve(n_acc);
    for (const auto& s : od.acc_sorted) {
        acc_offsets.push_back(static_cast<uint32_t>(acc_str_pool.size()));
        acc_str_pool.append(s);
        acc_str_pool.push_back('\0');
    }

    // Layout (all offsets from section start):
    //   [TxdbHeader(80)]
    //   [nodes: 32*n_nodes]    ← sizeof(TaxNodeRecord)=32 in v2
    //   [name_pool]
    //   [acc_offsets: 4*n_acc]
    //   [acc_taxids:  8*n_acc] ← uint64_t in v2
    //   [acc_buckets: 4*n_buckets]
    //   [acc_strings pool]

    uint64_t nodes_off     = sizeof(TxdbHeader);
    uint64_t name_pool_off = nodes_off     + sizeof(TaxNodeRecord) * n_nodes;
    uint64_t acc_off_off   = name_pool_off + od.name_pool.size();
    uint64_t acc_tid_off   = acc_off_off   + 4ULL * n_acc;
    uint64_t acc_buck_off  = acc_tid_off   + 8ULL * n_acc;  // uint64_t per entry
    uint64_t acc_str_off   = acc_buck_off  + 4ULL * n_buck;

    TxdbHeader hdr{};
    hdr.magic              = SEC_TXDB;
    hdr.version            = 2;
    hdr.flags              = 1; // DERIVED_FROM_LINEAGE
    hdr.n_nodes            = n_nodes;
    hdr.n_accessions       = n_acc;
    hdr.n_buckets          = n_buck;
    hdr.nodes_offset       = nodes_off;
    hdr.name_pool_offset   = name_pool_off;
    hdr.acc_offsets_offset = acc_off_off;
    hdr.acc_taxids_offset  = acc_tid_off;
    hdr.acc_buckets_offset = acc_buck_off;
    hdr.acc_strings_offset = acc_str_off;

    uint64_t section_start = writer.current_offset();

    writer.append(&hdr,                sizeof(hdr));
    writer.append(od.nodes.data(),     sizeof(TaxNodeRecord) * n_nodes);
    writer.append(od.name_pool.data(), od.name_pool.size());
    writer.append(acc_offsets.data(),  4ULL * n_acc);
    writer.append(od.acc_taxids.data(),8ULL * n_acc);
    writer.append(od.acc_buckets.data(),4ULL * n_buck);
    writer.append(acc_str_pool.data(), acc_str_pool.size());

    uint64_t section_end  = writer.current_offset();
    uint64_t payload_size = section_end - section_start;

    SectionDesc desc{};
    desc.type              = SEC_TXDB;
    desc.version           = 2;
    desc.flags             = 0;
    desc.section_id        = section_id;
    desc.file_offset       = section_start;
    desc.compressed_size   = payload_size;
    desc.uncompressed_size = payload_size;
    desc.item_count        = n_nodes;
    desc.aux0              = n_acc;
    desc.aux1              = 0;
    return desc;
}

// ── TxdbReader ────────────────────────────────────────────────────────────────

void TxdbReader::open(const uint8_t* base, uint64_t offset, uint64_t size) {
    if (size < sizeof(TxdbHeader))
        throw std::runtime_error("TxdbReader: section too small");

    const uint8_t* payload = base + offset;
    const auto* hdr = reinterpret_cast<const TxdbHeader*>(payload);

    if (hdr->magic != SEC_TXDB)
        throw std::runtime_error("TxdbReader: bad magic");
    if (hdr->version != 2)
        throw std::runtime_error("TxdbReader: unsupported version " +
                                 std::to_string(hdr->version) + " (expected 2)");

    TaxonomyTree::MmapView v;
    v.base          = base;
    v.size          = size;
    v.hdr           = hdr;
    v.n_nodes       = hdr->n_nodes;
    v.n_accessions  = hdr->n_accessions;
    v.n_buckets     = hdr->n_buckets;
    v.nodes         = reinterpret_cast<const TaxNodeRecord*>(payload + hdr->nodes_offset);
    v.name_pool     = reinterpret_cast<const char*>(payload + hdr->name_pool_offset);
    v.acc_offsets   = reinterpret_cast<const uint32_t*>(payload + hdr->acc_offsets_offset);
    v.acc_taxids    = reinterpret_cast<const uint64_t*>(payload + hdr->acc_taxids_offset);
    v.acc_buckets   = reinterpret_cast<const uint32_t*>(payload + hdr->acc_buckets_offset);
    v.acc_strings   = reinterpret_cast<const char*>(payload + hdr->acc_strings_offset);

    tree_ = TaxonomyTree(std::move(v));
}

// ── TaxonomyTree — unified accessors ─────────────────────────────────────────

const TaxNodeRecord* TaxonomyTree::nodes_ptr() const {
    if (owned_) return owned_->nodes.data();
    if (mmap_)  return mmap_->nodes;
    return nullptr;
}

uint32_t TaxonomyTree::nodes_count() const {
    if (owned_) return static_cast<uint32_t>(owned_->nodes.size());
    if (mmap_)  return mmap_->n_nodes;
    return 0;
}

const char* TaxonomyTree::name_pool_ptr() const {
    if (owned_) return owned_->name_pool.data();
    if (mmap_)  return mmap_->name_pool;
    return nullptr;
}

const uint32_t* TaxonomyTree::acc_offsets_ptr() const {
    if (mmap_) return mmap_->acc_offsets;
    return nullptr;
}

const uint64_t* TaxonomyTree::acc_taxids_ptr() const {
    if (owned_) return owned_->acc_taxids.data();
    if (mmap_)  return mmap_->acc_taxids;
    return nullptr;
}

const uint32_t* TaxonomyTree::acc_buckets_ptr() const {
    if (owned_) return owned_->acc_buckets.data();
    if (mmap_)  return mmap_->acc_buckets;
    return nullptr;
}

const char* TaxonomyTree::acc_strings_ptr() const {
    if (mmap_) return mmap_->acc_strings;
    return nullptr;
}

uint32_t TaxonomyTree::acc_buckets_count() const {
    if (owned_) return owned_->n_buckets;
    if (mmap_)  return mmap_->n_buckets;
    return 0;
}

uint32_t TaxonomyTree::acc_count() const {
    if (owned_) return static_cast<uint32_t>(owned_->acc_sorted.size());
    if (mmap_)  return mmap_->n_accessions;
    return 0;
}

size_t TaxonomyTree::n_nodes() const { return nodes_count(); }
size_t TaxonomyTree::n_accessions() const { return acc_count(); }

// Binary search by taxid (nodes sorted by taxid)
const TaxNodeRecord* TaxonomyTree::find_node(uint64_t taxid) const {
    const TaxNodeRecord* arr = nodes_ptr();
    uint32_t n = nodes_count();
    if (!arr || n == 0) return nullptr;

    uint32_t lo = 0, hi = n;
    while (lo < hi) {
        uint32_t mid = lo + (hi - lo) / 2;
        if (arr[mid].taxid == taxid) return &arr[mid];
        if (arr[mid].taxid < taxid)  lo = mid + 1;
        else                          hi = mid;
    }
    return nullptr;
}

// ── TaxonomyTree — query API ──────────────────────────────────────────────────

uint64_t TaxonomyTree::taxid_for_accession(std::string_view acc) const {
    uint32_t n_buck = acc_buckets_count();
    if (n_buck == 0) return 0;

    const uint32_t* bkts = acc_buckets_ptr();
    if (!bkts) return 0;

    uint32_t h = fnv1a(acc) & (n_buck - 1);

    if (owned_) {
        const auto& sorted = owned_->acc_sorted;
        const auto& taxids = owned_->acc_taxids;
        for (uint32_t probe = 0; probe < n_buck; ++probe) {
            uint32_t idx = bkts[(h + probe) & (n_buck - 1)];
            if (idx == UINT32_MAX) return 0;
            if (sorted[idx] == acc) return taxids[idx];
        }
    } else if (mmap_) {
        const uint32_t* offsets = mmap_->acc_offsets;
        const uint64_t* tid_arr = mmap_->acc_taxids;
        const char*     strs    = mmap_->acc_strings;
        uint32_t n_acc = mmap_->n_accessions;
        for (uint32_t probe = 0; probe < n_buck; ++probe) {
            uint32_t idx = bkts[(h + probe) & (n_buck - 1)];
            if (idx == UINT32_MAX) return 0;
            if (idx >= n_acc) return 0;
            if (acc == std::string_view(strs + offsets[idx])) return tid_arr[idx];
        }
    }
    return 0;
}

uint64_t TaxonomyTree::parent(uint64_t taxid) const {
    const TaxNodeRecord* n = find_node(taxid);
    if (!n) return 0;
    return n->parent_taxid;
}

TaxRank TaxonomyTree::rank(uint64_t taxid) const {
    const TaxNodeRecord* n = find_node(taxid);
    if (!n) return TaxRank::NO_RANK;
    return static_cast<TaxRank>(n->rank);
}

std::string_view TaxonomyTree::name(uint64_t taxid) const {
    const TaxNodeRecord* n = find_node(taxid);
    if (!n) return {};
    const char* pool = name_pool_ptr();
    if (!pool) return {};
    return std::string_view(pool + n->name_offset, n->name_len);
}

TaxNamespace TaxonomyTree::taxnamespace(uint64_t taxid) const {
    const TaxNodeRecord* n = find_node(taxid);
    if (!n) return TaxNamespace::NCBI;
    return static_cast<TaxNamespace>(n->taxnamespace);
}

TaxNodeStatus TaxonomyTree::status(uint64_t taxid) const {
    const TaxNodeRecord* n = find_node(taxid);
    if (!n) return TaxNodeStatus::TOMBSTONED;
    return static_cast<TaxNodeStatus>(n->status);
}

bool TaxonomyTree::is_synthetic(uint64_t taxid) const {
    const TaxNodeRecord* n = find_node(taxid);
    if (!n) return false;
    return (n->flags & 0x01u) != 0;
}

uint64_t TaxonomyTree::ancestor_at_rank(uint64_t taxid, TaxRank target) const {
    uint64_t cur = taxid;
    for (int i = 0; i < 40; ++i) {
        if (cur == 0 || cur == 1) return 0;
        const TaxNodeRecord* n = find_node(cur);
        if (!n) return 0;
        if (static_cast<TaxRank>(n->rank) == target) return cur;
        if (n->parent_taxid == cur) return 0;
        cur = n->parent_taxid;
    }
    return 0;
}

uint64_t TaxonomyTree::lca(uint64_t a, uint64_t b) const {
    std::vector<uint64_t> chain_a;
    uint64_t cur = a;
    for (int i = 0; i < 40 && cur > 1; ++i) {
        chain_a.push_back(cur);
        const TaxNodeRecord* n = find_node(cur);
        if (!n || n->parent_taxid == cur) break;
        cur = n->parent_taxid;
    }
    chain_a.push_back(1);

    std::unordered_set<uint64_t> anc_a(chain_a.begin(), chain_a.end());

    cur = b;
    for (int i = 0; i < 40; ++i) {
        if (anc_a.count(cur)) return cur;
        const TaxNodeRecord* n = find_node(cur);
        if (!n || n->parent_taxid == cur) return 1;
        cur = n->parent_taxid;
    }
    return 1;
}

TaxRank TaxonomyTree::lca_rank(uint64_t a, uint64_t b) const {
    return rank(lca(a, b));
}

void TaxonomyTree::scan_nodes(
    const std::function<void(uint64_t, uint64_t, TaxRank, std::string_view, bool)>& cb) const
{
    if (owned_) {
        for (const auto& rec : owned_->nodes) {
            std::string_view nm(owned_->name_pool.data() + rec.name_offset, rec.name_len);
            cb(rec.taxid, rec.parent_taxid,
               static_cast<TaxRank>(rec.rank),
               nm,
               (rec.flags & 0x01) != 0);
        }
    } else if (mmap_) {
        const auto& v = *mmap_;
        for (uint32_t i = 0; i < v.n_nodes; i++) {
            const auto& rec = v.nodes[i];
            std::string_view nm(v.name_pool + rec.name_offset, rec.name_len);
            cb(rec.taxid, rec.parent_taxid,
               static_cast<TaxRank>(rec.rank),
               nm,
               (rec.flags & 0x01) != 0);
        }
    }
}

void TaxonomyTree::scan_accessions(
    const std::function<void(std::string_view, uint64_t)>& cb) const
{
    if (owned_) {
        for (size_t i = 0; i < owned_->acc_sorted.size(); i++)
            cb(owned_->acc_sorted[i], owned_->acc_taxids[i]);
    } else if (mmap_) {
        const auto& v = *mmap_;
        for (uint32_t i = 0; i < v.n_accessions; i++) {
            const char* s = v.acc_strings + v.acc_offsets[i];
            cb(s, v.acc_taxids[i]);
        }
    }
}

} // namespace genopack
