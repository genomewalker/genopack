// Unit test for SEC_XQAL external-quality columns (Phase 3):
//   1. xqual_write emits the sorted union of genome_ids as rows; each column gets a
//      null bitmap (present iff that genome_id has a value for that column).
//   2. Round-trips through AppendWriter/ColStoreReader: values, presence, identity,
//      per-column checksums.
//   3. xqual_read_all reconstructs the full ColumnKey (so a re-ingest merges by
//      identity, not position).

#include <genopack/xqual_columns.hpp>
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <cmath>
#include <iostream>
#include "test_utils.hpp"

using namespace genopack;
using namespace genopack_test;

int main() {
    // CheckM2 completeness on genomes {10,20,30}; contamination only on {10,30}
    // (sparse — genome 20 must be null in the contamination column).
    XqualColumn comp;
    comp.key = ColumnKey{qual_axis::COMPLETENESS, qual_tool::CHECKM2, qual_method::ML,
                         Unit::Fraction01, 1, 0, 0, 0, "CheckM2 completeness"};
    comp.values = {{10, 0.95f}, {20, 0.50f}, {30, 0.10f}};

    XqualColumn cont;
    cont.key = ColumnKey{qual_axis::CONTAMINATION, qual_tool::CHECKM2, qual_method::ML,
                         Unit::Fraction01, 1, 0, 0, 0, "CheckM2 contamination"};
    cont.values = {{10, 0.02f}, {30, 0.07f}};

    TempDir tmp("xqual");
    const auto path = tmp.path / "xqual.bin";

    SectionDesc sd;
    {
        AppendWriter aw;
        aw.create(path);
        sd = xqual_write(aw, /*section_id=*/7, {comp, cont});
        aw.flush();
        aw.close();
    }
    require(sd.type == SEC_XQAL, "section type must be SEC_XQAL");

    MmapFileReader mm;
    mm.open(path);
    ColStoreReader r;
    r.open(mm.data(), sd.file_offset, sd.compressed_size);

    require(r.n_rows() == 3, "row set must be the union {10,20,30}");
    require(r.n_columns() == 2, "two columns expected");
    require(r.row_key(0) == 10 && r.row_key(1) == 20 && r.row_key(2) == 30,
            "row keys must be sorted ascending");
    require(r.verify_checksums(), "per-column checksums must verify");

    // Completeness column: all three present.
    const int ci = r.find_column(comp.key);
    require(ci >= 0, "completeness column must be found by identity");
    auto cv = r.column(static_cast<uint32_t>(ci));
    require(cv.tool == "checkm2" && cv.axis == "completeness", "reconstructed names");
    require(cv.unit == Unit::Fraction01, "unit preserved");
    require(std::abs(cv.get_f64(0) - 0.95) < 1e-6, "comp[10]");
    require(std::abs(cv.get_f64(1) - 0.50) < 1e-6, "comp[20]");
    require(std::abs(cv.get_f64(2) - 0.10) < 1e-6, "comp[30]");

    // Contamination column: genome 20 (row index 1) is absent -> NaN.
    const int ki = r.find_column(cont.key);
    require(ki >= 0, "contamination column must be found");
    auto kv = r.column(static_cast<uint32_t>(ki));
    require(kv.present(0) && !kv.present(1) && kv.present(2),
            "contamination present mask: {10,30} present, 20 null");
    require(std::isnan(kv.get_f64(1)), "absent cell must read NaN");
    require(std::abs(kv.get_f64(0) - 0.02) < 1e-6, "cont[10]");

    // Identity is content-derived, not positional.
    require(cv.identity_hash == column_identity_hash(comp.key),
            "stored identity must equal recomputed identity");
    require(cv.identity_hash != kv.identity_hash, "distinct columns distinct identity");

    // Read-all reconstructs keys + values for merge.
    auto back = xqual_read_all(r);
    require(back.size() == 2, "read_all returns both columns");
    bool saw_comp = false;
    for (const auto& x : back) {
        if (column_identity_hash(x.key) == column_identity_hash(comp.key)) {
            saw_comp = true;
            require(x.values.size() == 3, "completeness has 3 values");
            require(std::abs(x.values.at(30) - 0.10f) < 1e-6f, "round-trip comp[30]");
        }
    }
    require(saw_comp, "completeness column round-tripped");

    std::cout << "genopack xqual test passed\n";
    return 0;
}
