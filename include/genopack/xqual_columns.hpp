#pragma once
// ── XQAL: external quality as a columnar store ────────────────────────────────
// SEC_XQAL holds quality measured by EXTERNAL tools (CheckM2, anvi'o, …) as
// provenance-carrying columns keyed by genome_id — the same colstore container as
// intrinsic QCOL, but with tool != "genopack". Ingest is ADDITIVE: a new tool's
// columns merge with whatever XQAL already exists (same-identity columns are
// overwritten), so CheckM2 and anvi'o measurements accumulate in one section.
// Fusing intrinsic + external into a single reported number is an explicit later
// step (CQAL / profiles), never silent here: a CheckM2 completeness and a genopack
// completeness live in distinct columns that only differ in their `tool` field.
#include "colstore.hpp"
#include "quality_schema.hpp"
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace genopack {

// One external-quality column to (re)write: full provenance identity + per-genome
// values. Values are held as float (interpreted per key.unit: Fraction01/Percent/…).
struct XqualColumn {
    ColumnKey                           key;
    std::unordered_map<uint64_t, float> values;   // genome_id -> value
};

// Write SEC_XQAL. The row set is the sorted union of every genome_id that appears
// in any column; each column carries a null bitmap (present iff that genome_id is
// in its map). Columns with no values are skipped.
SectionDesc xqual_write(AppendWriter& w, uint64_t section_id,
                        const std::vector<XqualColumn>& columns);

// Read every column of an existing SEC_XQAL back into XqualColumn form so a new
// ingest can merge into it. Reconstructs the full ColumnKey from the stored
// ColumnDesc, so identity is preserved byte-for-byte across a re-ingest.
std::vector<XqualColumn> xqual_read_all(const ColStoreReader& r);

// ── Generic colstore column helpers (shared by XQAL and CQAL) ──────────────────
// Write any SEC_* columnar section from a set of named float columns: the row set
// is the sorted union of genome_ids, each column null-bitmapped (present iff the
// genome_id is in its map). read_named_columns is the inverse (reconstructs the
// full ColumnKey from the stored ColumnDesc). xqual_write/xqual_read_all are thin
// wrappers fixing section_type = SEC_XQAL.
SectionDesc write_named_columns(uint32_t section_type, AppendWriter& w,
                                uint64_t section_id, const std::vector<XqualColumn>& columns);
std::vector<XqualColumn> read_named_columns(const ColStoreReader& r);

} // namespace genopack
