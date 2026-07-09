#pragma once
// ── QCOL: intrinsic quality as a columnar store ───────────────────────────────
// SEC_QCOL replaces the flat 80-byte SEC_QUAL POD. Each legacy QualRecord field
// becomes a named, provenance-carrying column (axis/tool=genopack/method/unit)
// stored byte-for-byte in its original dtype — so materializing a QualRecord back
// from QCOL is exact (no quantization drift), and check/calibrate/subset that
// consume QualRecord via ArchiveReader::scan_qual need no change. New columns
// (e.g. COMPLETENESS/genopack/AAMER_GENUS_CORE, referencing core_model_hash) are
// added alongside and simply bypass the QualRecord interchange.
#include "qual.hpp"
#include "colstore.hpp"
#include "quality_schema.hpp"
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <functional>
#include <unordered_map>
#include <vector>

namespace genopack {

// One legacy QualRecord field <-> one QCOL column. `offset` is the field's byte
// offset in QualRecord; the value is copied verbatim (sentinels intact), so the
// column dtype mirrors the field type and round-trip is byte-exact.
struct QcolField {
    const char* axis;
    const char* method;
    Unit        unit;
    ColDType    dtype;
    size_t      offset;
};

inline const std::vector<QcolField>& qcol_fields() {
    static const std::vector<QcolField> F = {
        {qual_axis::COMPLETENESS,  qual_method::CLUSTER_RELATIVE, Unit::Fraction01, ColDType::F32, offsetof(QualRecord, completeness_cluster_relative)},
        {qual_axis::COMPLETENESS,  qual_method::POST_DECONTAM,    Unit::Fraction01, ColDType::F32, offsetof(QualRecord, completeness_post_decontam)},
        {qual_axis::CONTAMINATION, qual_method::LEAKAGE,          Unit::Fraction01, ColDType::F32, offsetof(QualRecord, contamination_leakage)},
        {qual_axis::CONTAMINATION, qual_method::TNF_EXCESS,       Unit::Ratio,      ColDType::F32, offsetof(QualRecord, contamination_tnf_excess)},
        {qual_axis::SUPPORT,       qual_method::TIER,             Unit::Count,      ColDType::U8,  offsetof(QualRecord, support_tier)},
        {qual_axis::CONTAMINATION, qual_method::SPE_OUTLIER,      Unit::Fraction01, ColDType::U8,  offsetof(QualRecord, spe_outlier_u8)},
        {qual_axis::CONTAMINATION, qual_method::RHO_OUTLIER,      Unit::Fraction01, ColDType::U8,  offsetof(QualRecord, rho_outlier_u8)},
        {qual_axis::SUPPORT,       qual_method::FLAGS,            Unit::Count,      ColDType::U8,  offsetof(QualRecord, qual_flags)},
        {qual_axis::CONTAMINATION, qual_method::CONTIG_OUTLIER,   Unit::Fraction01, ColDType::U16, offsetof(QualRecord, contig_outlier_u16)},
        {qual_axis::CONTAMINATION, qual_method::FMH_MINORITY,     Unit::Fraction01, ColDType::U16, offsetof(QualRecord, fmh_minority_u16)},
        {qual_axis::COMPLETENESS,  qual_method::MARKER_SCG,       Unit::Fraction01, ColDType::U8,  offsetof(QualRecord, marker_completeness_u8)},
        {qual_axis::CONTAMINATION, qual_method::CROSS_GENUS,      Unit::Fraction01, ColDType::U8,  offsetof(QualRecord, cross_genus_u8)},
        {qual_axis::CONTAMINATION, qual_method::DUPLICATION,      Unit::Fraction01, ColDType::U16, offsetof(QualRecord, contamination_duplication_u16)},
        {qual_axis::SUPPORT,       qual_method::QUALITY_TIER,     Unit::Count,      ColDType::U8,  offsetof(QualRecord, quality_tier_u8)},
        // Phase-2 estimators (F32 byte-transpose; NAN sentinel preserved, so old QCOL sections
        // that lack these columns materialize back to make_empty's NAN default).
        {qual_axis::CONTAMINATION, qual_method::CORE_DUP_MASS,    Unit::Ratio,      ColDType::F32, offsetof(QualRecord, contamination_core_dup_mass)},
        {qual_axis::COMPLETENESS,  qual_method::ACCESSORY_RATIO,  Unit::Ratio,      ColDType::F32, offsetof(QualRecord, accessory_ratio)},
    };
    return F;
}

inline ColumnKey qcol_legacy_key(const QcolField& f) {
    return ColumnKey{f.axis, qual_tool::GENOPACK, f.method, f.unit, /*version=*/1};
}

// New (non-QualRecord) columns to add alongside the legacy transpose.
struct QcolExtraColumns {
    uint64_t                                        core_model_hash = 0;
    const std::unordered_map<uint64_t, float>*      aamer_genus_core = nullptr; // genome_id -> coverage
    uint64_t                                        family_core_model_hash = 0;
    const std::unordered_map<uint64_t, float>*      aamer_family_core = nullptr; // genome_id -> coverage (FCORE)
};

// Build SEC_QCOL from per-genome QualRecords (legacy-field transpose, byte-exact)
// plus any extra columns. `recs` is sorted by genome_id (the row key).
inline SectionDesc qcol_write(AppendWriter& w, uint64_t section_id,
                              std::vector<QualRecord> recs,
                              const QcolExtraColumns& extra = {}) {
    std::sort(recs.begin(), recs.end(),
              [](const QualRecord& a, const QualRecord& b) { return a.genome_id < b.genome_id; });
    // Deduplicate by genome_id: deduped archives can alias several accessions to one
    // genome_id, which yields multiple records here. The colstore row key must be
    // strictly ascending, so collapse to one record per genome_id (keep the first).
    recs.erase(std::unique(recs.begin(), recs.end(),
               [](const QualRecord& a, const QualRecord& b) { return a.genome_id == b.genome_id; }),
               recs.end());

    ColStoreWriter cw(SEC_QCOL);
    std::vector<uint64_t> rows;
    rows.reserve(recs.size());
    for (const auto& r : recs) rows.push_back(r.genome_id);
    cw.set_rows(rows);

    std::vector<uint8_t> buf;
    for (const QcolField& f : qcol_fields()) {
        const uint64_t sz = dtype_size(f.dtype);
        buf.assign(static_cast<size_t>(recs.size() * sz), 0);
        for (size_t i = 0; i < recs.size(); ++i)
            std::memcpy(buf.data() + i * sz,
                        reinterpret_cast<const uint8_t*>(&recs[i]) + f.offset, sz);
        cw.add_column_raw(qcol_legacy_key(f), f.dtype, buf.data(), recs.size());
    }

    if (extra.aamer_genus_core) {
        std::vector<float>   vals(recs.size(), 0.0f);
        std::vector<uint8_t> present(recs.size(), 0);
        for (size_t i = 0; i < recs.size(); ++i) {
            auto it = extra.aamer_genus_core->find(recs[i].genome_id);
            if (it != extra.aamer_genus_core->end()) { vals[i] = it->second; present[i] = 1; }
        }
        ColumnKey k{qual_axis::COMPLETENESS, qual_tool::GENOPACK, qual_method::AAMER_GENUS_CORE,
                    Unit::Fraction01, /*version=*/1, /*ref_db=*/0,
                    /*core_model=*/extra.core_model_hash, /*calib=*/0};
        cw.add_column<float>(k, vals, present);
    }

    if (extra.aamer_family_core) {
        std::vector<float>   vals(recs.size(), 0.0f);
        std::vector<uint8_t> present(recs.size(), 0);
        for (size_t i = 0; i < recs.size(); ++i) {
            auto it = extra.aamer_family_core->find(recs[i].genome_id);
            if (it != extra.aamer_family_core->end()) { vals[i] = it->second; present[i] = 1; }
        }
        ColumnKey k{qual_axis::COMPLETENESS, qual_tool::GENOPACK, qual_method::AAMER_FAMILY_CORE,
                    Unit::Fraction01, /*version=*/1, /*ref_db=*/0,
                    /*core_model=*/extra.family_core_model_hash, /*calib=*/0};
        cw.add_column<float>(k, vals, present);
    }

    return cw.finalize(w, section_id);
}

// Materialize each row of a SEC_QCOL store back into a QualRecord (legacy fields
// + aamer_core extra columns; absent columns keep make_empty defaults).
// Drop-in for QualReader::scan.
inline void qcol_scan(const ColStoreReader& r,
                      const std::function<void(const QualRecord&)>& cb) {
    const auto& fields = qcol_fields();
    std::vector<int> idx(fields.size());
    for (size_t fi = 0; fi < fields.size(); ++fi) {
        int ci = r.find_column(qcol_legacy_key(fields[fi]));
        // Skip if on-disk dtype doesn't match: stale U8 column for a field now stored as U16
        // (or any other format change). Lets the fallback path (make_empty defaults) stay consistent.
        if (ci >= 0 && r.column(static_cast<uint32_t>(ci)).dtype != fields[fi].dtype) ci = -1;
        idx[fi] = ci;
    }

    // Extra columns stored with core_model_hash in key — find by method name.
    int aamer_core_col = -1, aamer_family_col = -1;
    for (uint32_t c = 0; c < r.n_columns(); ++c) {
        auto cv = r.column(c);
        if (cv.method == qual_method::AAMER_GENUS_CORE)  aamer_core_col   = static_cast<int>(c);
        if (cv.method == qual_method::AAMER_FAMILY_CORE) aamer_family_col = static_cast<int>(c);
    }

    const uint32_t n = r.n_rows();
    for (uint32_t row = 0; row < n; ++row) {
        QualRecord q = QualRecord::make_empty(r.row_key(row));
        for (size_t fi = 0; fi < fields.size(); ++fi) {
            if (idx[fi] < 0) continue;
            auto cv = r.column(static_cast<uint32_t>(idx[fi]));
            const uint64_t sz = dtype_size(fields[fi].dtype);
            std::memcpy(reinterpret_cast<uint8_t*>(&q) + fields[fi].offset,
                        cv.values + static_cast<size_t>(row) * sz, sz);
        }
        if (aamer_core_col >= 0) {
            auto cv = r.column(static_cast<uint32_t>(aamer_core_col));
            if (cv.present(row)) q.completeness_aamer_core = cv.get<float>(row);
        }
        if (aamer_family_col >= 0) {
            auto cv = r.column(static_cast<uint32_t>(aamer_family_col));
            if (cv.present(row)) q.completeness_aamer_family_core = cv.get<float>(row);
        }
        cb(q);
    }
}

} // namespace genopack
