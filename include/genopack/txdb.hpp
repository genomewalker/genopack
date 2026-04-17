#pragma once
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace genopack {

// ── Taxid namespace convention ────────────────────────────────────────────────
// Stored as uint64_t throughout.
// bit 63 = 0 : NCBI taxid (native 32-bit value; max ~3.5M, fits easily)
// bit 63 = 1 : GTDB synthetic concept_id (FNV-64 of canonical path + namespace bit)
// Special: taxid=1 is the common root (no namespace bit).
static constexpr uint64_t TAXID_GTDB_BIT   = (1ULL << 63);
static constexpr uint64_t TAXID_NCBI_MASK  = ~TAXID_GTDB_BIT;

inline bool taxid_is_gtdb(uint64_t t) { return (t & TAXID_GTDB_BIT) != 0; }
inline bool taxid_is_ncbi(uint64_t t) { return (t & TAXID_GTDB_BIT) == 0; }

enum class TaxRank : uint8_t {
    NO_RANK = 0,
    DOMAIN  = 1,  // d__
    LIFE    = 2,  // l__ (mOTUs extended format)
    KINGDOM = 3,  // k__ (mOTUs extended format)
    PHYLUM  = 4,  // p__
    CLASS   = 5,  // c__
    ORDER   = 6,  // o__
    FAMILY  = 7,  // f__
    GENUS   = 8,  // g__
    SPECIES = 9   // s__
};

enum class TaxNamespace : uint8_t {
    NCBI = 0,
    GTDB = 1,
};

enum class TaxNodeStatus : uint8_t {
    ACTIVE         = 0,
    MERGED         = 1,  // auto-follow via GTAX
    RENAMED        = 2,  // auto-follow, old name remains queryable
    RETIRED_SPLIT  = 3,  // AMBIGUOUS — never auto-follow; see GTAX for successors
    TOMBSTONED     = 4,  // removed, no successor
    PROVISIONAL    = 5,  // part-local only; must not appear in finalized archive
};

// ── TxdbHeader: 80 bytes (unchanged size, version bumped to 2) ────────────────
// Version 2 changes vs version 1:
//   - TaxNodeRecord is 32 bytes (was 16)
//   - acc_taxids_offset points to uint64_t[n_accessions] (was uint32_t[])
struct TxdbHeader {
    uint32_t magic;              // SEC_TXDB
    uint16_t version;            // 2
    uint16_t flags;              // bit 0: DERIVED_FROM_LINEAGE
    uint32_t n_nodes;
    uint32_t n_accessions;
    uint32_t n_buckets;          // power-of-2 hash table size
    uint32_t _pad;
    uint64_t _pad2;
    uint64_t nodes_offset;       // from section start → TaxNodeRecord array
    uint64_t name_pool_offset;   // from section start → null-terminated name strings
    uint64_t acc_offsets_offset; // from section start → uint32_t[n_accessions] into acc_strings
    uint64_t acc_taxids_offset;  // from section start → uint64_t[n_accessions] deepest taxid
    uint64_t acc_buckets_offset; // from section start → uint32_t[n_buckets] FNV hash table
    uint64_t acc_strings_offset; // from section start → null-terminated accession strings
};
static_assert(sizeof(TxdbHeader) == 80);

// ── TaxNodeRecord: 32 bytes (was 16 in version 1), sorted by taxid ───────────
struct TaxNodeRecord {
    uint64_t taxid;         // 8
    uint64_t parent_taxid;  // 8
    uint32_t name_offset;   // 4 — into name_pool
    uint16_t name_len;      // 2
    uint8_t  rank;          // 1 — TaxRank as uint8_t
    uint8_t  taxnamespace;  // 1 — TaxNamespace as uint8_t
    uint8_t  status;        // 1 — TaxNodeStatus as uint8_t
    uint8_t  flags;         // 1 — bit 0: synthetic, bit 1: per_accession
    uint8_t  _pad[6];       // 6
};
static_assert(sizeof(TaxNodeRecord) == 32);

// ── TaxonomyTree ──────────────────────────────────────────────────────────────
class TaxonomyTree {
public:
    TaxonomyTree() = default;

    // Query interface (all taxids are uint64_t)
    uint64_t taxid_for_accession(std::string_view acc) const;
    uint64_t parent(uint64_t taxid) const;
    TaxRank  rank(uint64_t taxid) const;
    std::string_view name(uint64_t taxid) const;
    TaxNamespace  taxnamespace(uint64_t taxid) const;
    TaxNodeStatus status(uint64_t taxid) const;
    bool     is_synthetic(uint64_t taxid) const;

    // Walk up tree to find ancestor at exactly this rank; returns 0 if not found
    uint64_t ancestor_at_rank(uint64_t taxid, TaxRank rank) const;

    // LCA of two taxids
    uint64_t lca(uint64_t a, uint64_t b) const;
    TaxRank  lca_rank(uint64_t a, uint64_t b) const;

    size_t n_nodes() const;
    size_t n_accessions() const;
    bool empty() const { return n_nodes() == 0; }

    static uint32_t fnv1a(std::string_view s);

    // Compute the stable GTDB concept_id for a canonical path string.
    // Returns FNV-64(path) | TAXID_GTDB_BIT. Collision-free up to ~2 billion unique paths.
    // This is the same function used internally by TxdbWriter::build_tree().
    static uint64_t concept_id_for_path(std::string_view canonical_path);

    // Iterate all nodes: cb(taxid, parent_taxid, rank, name, is_synthetic)
    void scan_nodes(const std::function<void(uint64_t taxid, uint64_t parent_taxid,
                                             TaxRank rank, std::string_view name,
                                             bool is_synthetic)>& cb) const;

    // Iterate all accession→taxid pairs
    void scan_accessions(const std::function<void(std::string_view accession,
                                                  uint64_t taxid)>& cb) const;

    // Internal: built-from-lineage path
    struct OwnedData {
        std::vector<TaxNodeRecord> nodes;        // sorted by taxid
        std::string                name_pool;
        std::vector<std::string>   acc_sorted;   // sorted accession strings
        std::vector<uint64_t>      acc_taxids;   // parallel to acc_sorted
        std::vector<uint32_t>      acc_buckets;  // FNV hash table (indices into acc_sorted)
        uint32_t                   n_buckets = 0;
    };
    explicit TaxonomyTree(OwnedData d)
        : owned_(std::make_shared<OwnedData>(std::move(d))) {}

    const OwnedData* owned_data() const { return owned_.get(); }

    // Internal: mmap'd path
    struct MmapView {
        const uint8_t*       base         = nullptr;
        uint64_t             size         = 0;
        const TxdbHeader*    hdr          = nullptr;
        const TaxNodeRecord* nodes        = nullptr;
        const char*          name_pool    = nullptr;
        const uint32_t*      acc_offsets  = nullptr;
        const uint64_t*      acc_taxids   = nullptr;  // uint64_t in v2
        const uint32_t*      acc_buckets  = nullptr;
        const char*          acc_strings  = nullptr;
        uint32_t n_nodes       = 0;
        uint32_t n_accessions  = 0;
        uint32_t n_buckets     = 0;
    };
    explicit TaxonomyTree(MmapView v)
        : mmap_(std::make_shared<MmapView>(v)) {}

private:
    std::shared_ptr<OwnedData> owned_;
    std::shared_ptr<MmapView>  mmap_;

    const TaxNodeRecord* nodes_ptr() const;
    uint32_t             nodes_count() const;
    const char*          name_pool_ptr() const;
    const uint32_t*      acc_offsets_ptr() const;
    const uint64_t*      acc_taxids_ptr() const;
    const uint32_t*      acc_buckets_ptr() const;
    const char*          acc_strings_ptr() const;
    uint32_t             acc_buckets_count() const;
    uint32_t             acc_count() const;

    const TaxNodeRecord* find_node(uint64_t taxid) const;
};

// ── TxdbWriter ────────────────────────────────────────────────────────────────
class TxdbWriter {
public:
    void add(std::string_view accession, std::string_view lineage_string);

    void set_kmer_profiles(
        const std::unordered_map<std::string, std::array<float, 136>>& profiles,
        float similarity_threshold = 0.97f);

    TaxonomyTree build_tree() const;

    SectionDesc finalize(AppendWriter& writer, uint64_t section_id);

private:
    std::vector<std::pair<std::string, std::string>> entries_;
    const std::unordered_map<std::string, std::array<float, 136>>* kmer_profiles_ = nullptr;
    float similarity_threshold_ = 0.97f;
};

// ── TxdbReader ───────────────────────────────────────────────────────────────
class TxdbReader {
public:
    void open(const uint8_t* base, uint64_t offset, uint64_t size);
    TaxonomyTree tree() const { return tree_; }
private:
    TaxonomyTree tree_;
};

} // namespace genopack
