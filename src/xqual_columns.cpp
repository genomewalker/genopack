#include <genopack/xqual_columns.hpp>
#include <set>

namespace genopack {

SectionDesc write_named_columns(uint32_t section_type, AppendWriter& w,
                                uint64_t section_id, const std::vector<XqualColumn>& columns) {
    // Sorted union of all genome_ids -> row keys.
    std::set<uint64_t> row_set;
    for (const auto& c : columns)
        for (const auto& kv : c.values) row_set.insert(kv.first);
    std::vector<uint64_t> rows(row_set.begin(), row_set.end());

    std::unordered_map<uint64_t, uint32_t> idx;
    idx.reserve(rows.size() * 2 + 16);
    for (uint32_t i = 0; i < rows.size(); ++i) idx[rows[i]] = i;

    ColStoreWriter cw(section_type);
    cw.set_rows(rows);

    for (const auto& c : columns) {
        if (c.values.empty()) continue;
        std::vector<float>   vals(rows.size(), 0.0f);
        std::vector<uint8_t> present(rows.size(), 0);
        for (const auto& [gid, v] : c.values) {
            auto it = idx.find(gid);
            if (it == idx.end()) continue;
            vals[it->second]    = v;
            present[it->second] = 1;
        }
        cw.add_column<float>(c.key, vals, present);
    }
    return cw.finalize(w, section_id);
}

SectionDesc xqual_write(AppendWriter& w, uint64_t section_id,
                        const std::vector<XqualColumn>& columns) {
    return write_named_columns(SEC_XQAL, w, section_id, columns);
}

std::vector<XqualColumn> read_named_columns(const ColStoreReader& r) {
    std::vector<XqualColumn> out;
    const uint32_t nc = r.n_columns();
    out.reserve(nc);
    for (uint32_t c = 0; c < nc; ++c) {
        ColStoreReader::ColumnView cv = r.column(c);
        XqualColumn xc;
        xc.key.axis             = std::string(cv.axis);
        xc.key.tool             = std::string(cv.tool);
        xc.key.method           = std::string(cv.method);
        xc.key.unit             = cv.unit;
        xc.key.version          = cv.version;
        xc.key.ref_db_hash      = cv.ref_db_hash;
        xc.key.core_model_hash  = cv.core_model_hash;
        xc.key.calibration_hash = cv.calibration_hash;
        xc.key.label            = std::string(cv.label);
        for (uint32_t i = 0; i < cv.n_rows; ++i)
            if (cv.present(i))
                xc.values[r.row_key(i)] = static_cast<float>(cv.get_f64(i));
        out.push_back(std::move(xc));
    }
    return out;
}

std::vector<XqualColumn> xqual_read_all(const ColStoreReader& r) {
    return read_named_columns(r);
}

} // namespace genopack
