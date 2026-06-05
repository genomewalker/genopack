// Unit test for the columnar quality store (Phase 0 foundation):
//   1. Build a ColStore section (sorted genome_id rows + typed columns, one with
//      a null) via ColStoreWriter; write through AppendWriter to a temp file.
//   2. Re-read via ColStoreReader: row-key binary search, per-column values,
//      presence/null bitmap, get_f64 widening, find-by-identity, checksum verify.
//   3. Corrupt one value byte in a heap copy and assert verify_checksums() fails.
//   4. Column-identity hash: stable for equal keys, distinct for different method,
//      and the on-disk identity matches the recomputed key hash (merge key).

#include <genopack/colstore.hpp>
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/quality_schema.hpp>
#include <cmath>
#include <iostream>
#include <vector>
#include "test_utils.hpp"

using namespace genopack;
using namespace genopack_test;

static ColumnKey comp_key() {
    return ColumnKey{qual_axis::COMPLETENESS, qual_tool::GENOPACK, qual_method::AAMER_GENUS_CORE,
                     Unit::Fraction01, /*version=*/1, /*ref_db=*/0xABCDull,
                     /*core_model=*/0x1234ull, /*calib=*/0, "completeness (genus core)"};
}

int main() {
    const std::vector<uint64_t> rows = {10, 20, 30, 40, 50};
    const std::vector<float>    comp = {0.95f, 0.80f, 0.0f, 0.50f, 1.00f};
    const std::vector<uint8_t>  comp_present = {1, 1, 0, 1, 1};  // row 30 is null
    const std::vector<uint8_t>  tier = {0, 0, 0xFF, 0, 0};
    const std::vector<uint32_t> markers = {120, 118, 0, 60, 122};

    const ColumnKey ckey = comp_key();
    const ColumnKey tkey{qual_axis::SUPPORT, qual_tool::GENOPACK, qual_method::TIER, Unit::Count};
    const ColumnKey mkey{qual_axis::COMPLETENESS, qual_tool::GENOPACK, qual_method::MARKER_SCG,
                         Unit::Count};

    // ── identity-hash properties ──
    require(column_identity_hash(ckey) == column_identity_hash(comp_key()),
            "identity hash must be stable for equal keys");
    {
        ColumnKey other = ckey;
        other.method = qual_method::CLUSTER_RELATIVE;
        require(column_identity_hash(other) != column_identity_hash(ckey),
                "different method must change identity hash");
        ColumnKey label_only = ckey;
        label_only.label = "totally different label";
        require(column_identity_hash(label_only) == column_identity_hash(ckey),
                "label must NOT affect identity hash");
    }

    TempDir tmp("colstore");
    const auto path = tmp.path / "qcol.bin";

    SectionDesc sd;
    {
        ColStoreWriter w(SEC_QCOL);
        w.set_rows(rows);
        w.add_column<float>(ckey, comp, comp_present);
        w.add_column<uint8_t>(tkey, tier);
        w.add_column<uint32_t>(mkey, markers);

        AppendWriter aw;
        aw.create(path);
        sd = w.finalize(aw, /*section_id=*/7);
        aw.flush();
        aw.close();
    }
    require(sd.type == SEC_QCOL, "section type mismatch");
    require(sd.item_count == rows.size(), "item_count should equal row count");
    require(sd.aux0 == 3, "aux0 should equal column count");
    require(sd.header_size == sizeof(ColStoreHeader), "header_size must be set");

    // ── read back ──
    {
        MmapFileReader mm;
        mm.open(path);
        ColStoreReader r;
        r.open(mm.data(), sd.file_offset, sd.compressed_size);

        require(r.n_rows() == rows.size(), "n_rows mismatch");
        require(r.n_columns() == 3, "n_columns mismatch");
        require(r.verify_checksums(), "per-column checksums must verify on a clean read");

        // row-key binary search
        require(r.row_index(30) == 2, "row_index(30) should be 2");
        require(r.row_index(10) == 0, "row_index(10) should be 0");
        require(r.row_index(50) == 4, "row_index(50) should be 4");
        require(r.row_index(35) == -1, "row_index(35) absent should be -1");
        require(r.row_index(0) == -1,  "row_index(0) absent should be -1");

        // completeness column (by identity)
        int ci = r.find_column(ckey);
        require(ci >= 0, "completeness column must be found by key");
        auto cv = r.column(static_cast<uint32_t>(ci));
        require(cv.identity_hash == column_identity_hash(ckey), "on-disk identity must match key hash");
        require(cv.axis == "completeness" && cv.tool == "genopack"
                && cv.method == "aamer_genus_core", "resolved names mismatch");
        require(cv.label == "completeness (genus core)", "label must round-trip");
        require(cv.unit == Unit::Fraction01, "unit mismatch");
        require(cv.dtype == ColDType::F32, "dtype mismatch");
        require(cv.ref_db_hash == 0xABCDull, "ref_db_hash must round-trip");
        require(cv.core_model_hash == 0x1234ull, "core_model_hash must round-trip");

        require(cv.present(0) && !cv.present(2) && cv.present(4), "presence bitmap mismatch");
        require(std::abs(cv.get<float>(0) - 0.95f) < 1e-6f, "comp[0] value mismatch");
        require(std::isnan(cv.get_f64(2)), "null cell must read as NaN via get_f64");
        require(std::abs(cv.get_f64(4) - 1.0) < 1e-6, "comp[4] get_f64 mismatch");

        // tier column (all present, u8)
        int ti = r.find_column(tkey);
        require(ti >= 0, "tier column must be found");
        auto tv = r.column(static_cast<uint32_t>(ti));
        require(tv.dtype == ColDType::U8, "tier dtype mismatch");
        require(tv.present(2) && tv.get<uint8_t>(2) == 0xFF, "tier value/presence mismatch");

        // markers column (u32)
        int mi = r.find_column(mkey);
        require(mi >= 0, "markers column must be found");
        auto mv = r.column(static_cast<uint32_t>(mi));
        require(mv.dtype == ColDType::U32, "markers dtype mismatch");
        require(mv.get<uint32_t>(4) == 122, "markers[4] mismatch");
        require(std::abs(mv.get_f64(0) - 120.0) < 1e-9, "markers get_f64 widening mismatch");

        // a never-added column must not be found
        require(r.find_column(ColumnKey{qual_axis::CONTAMINATION, qual_tool::CHECKM2,
                                        qual_method::ML, Unit::Percent}) == -1,
                "absent column must not be found");
    }

    // ── corruption detection ──
    {
        MmapFileReader mm;
        mm.open(path);
        std::vector<uint8_t> buf(mm.data(), mm.data() + mm.size());

        ColStoreReader r;
        r.open(buf.data(), sd.file_offset, sd.compressed_size);
        require(r.verify_checksums(), "checksums should verify before corruption");

        auto cv = r.column(0);
        size_t off = static_cast<size_t>(cv.values - buf.data());
        buf[off] ^= 0xFFu;  // flip a byte of the first column's first value
        require(!r.verify_checksums(), "verify_checksums must detect a corrupted value byte");
    }

    std::cout << "genopack colstore test passed\n";
    return 0;
}
