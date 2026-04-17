#pragma once
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <cstdint>
#include <string>
#include <vector>

namespace genopack {

// ── GtaxAliasType ─────────────────────────────────────────────────────────────
enum class GtaxAliasType : uint8_t {
    MERGED         = 0,  // old_taxid merged into new_taxid; auto-follow
    RENAMED        = 1,  // old name → new name, same concept; auto-follow; old name queryable
    RETIRED_SPLIT  = 2,  // old_taxid split; new_taxid is ONE successor (may have multiple rows)
                         // NEVER auto-follow — caller must handle AMBIGUOUS
    TOMBSTONED     = 3,  // old_taxid removed, no successor; new_taxid=0
};

// ── GtaxHeader: 64 bytes ──────────────────────────────────────────────────────
struct GtaxHeader {
    uint32_t magic;      // SEC_GTAX
    uint16_t version;    // 1
    uint16_t flags;
    uint64_t n_entries;
    uint8_t  _pad[48];
};
static_assert(sizeof(GtaxHeader) == 64);

// ── GtaxEntry: 24 bytes ───────────────────────────────────────────────────────
// Entries are stored sorted by old_taxid for binary search.
// For RETIRED_SPLIT: multiple entries share the same old_taxid, each pointing
// to a different successor. Query: if multiple entries found → AMBIGUOUS.
struct GtaxEntry {
    uint64_t old_taxid;     // 8
    uint64_t new_taxid;     // 8 (0 for TOMBSTONED)
    uint32_t source_ver;    // 4 — YYYYMMDD of the release that introduced this alias
    uint8_t  alias_type;    // 1 — GtaxAliasType
    uint8_t  _pad[3];       // 3
};
static_assert(sizeof(GtaxEntry) == 24);

// ── GtaxWriter ────────────────────────────────────────────────────────────────
class GtaxWriter {
public:
    void add(uint64_t old_taxid, uint64_t new_taxid,
             GtaxAliasType type, uint32_t source_ver = 0);

    SectionDesc finalize(AppendWriter& writer, uint64_t section_id);

private:
    std::vector<GtaxEntry> entries_;
};

// ── GtaxReader ────────────────────────────────────────────────────────────────
// Result of resolving an old_taxid through the alias table.
struct GtaxResult {
    enum class Kind { NOT_FOUND, ACTIVE, MERGED, RENAMED, AMBIGUOUS, TOMBSTONED };
    Kind                   kind       = Kind::NOT_FOUND;
    uint64_t               taxid      = 0;   // resolved taxid (if unique)
    std::vector<uint64_t>  successors;        // populated for AMBIGUOUS (split)
};

class GtaxReader {
public:
    void open(const uint8_t* base, uint64_t offset, uint64_t size);

    // Resolve old_taxid through the alias table (one level only).
    // Callers must detect compound chains (merge→split) by checking Kind::AMBIGUOUS
    // on the result of following a MERGED chain.
    GtaxResult resolve(uint64_t old_taxid) const;

    size_t n_entries() const { return hdr_ ? static_cast<size_t>(hdr_->n_entries) : 0; }

private:
    const GtaxHeader* hdr_     = nullptr;
    const GtaxEntry*  entries_ = nullptr;

    // Returns index of first entry with old_taxid, or n_entries if not found.
    size_t lower_bound(uint64_t old_taxid) const;
};

} // namespace genopack
