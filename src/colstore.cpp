#include <genopack/colstore.hpp>
#include <algorithm>
#include <cmath>
#include <cstring>

namespace genopack {

// ── checksum helper ───────────────────────────────────────────────────────────
static uint32_t checksum32(const uint8_t* a, size_t na, const uint8_t* b, size_t nb) {
    XXH3_state_t* st = XXH3_createState();
    XXH3_64bits_reset(st);
    if (na) XXH3_64bits_update(st, a, na);
    if (nb) XXH3_64bits_update(st, b, nb);
    uint64_t h = XXH3_64bits_digest(st);
    XXH3_freeState(st);
    return static_cast<uint32_t>(h);
}

// ── ColStoreWriter ────────────────────────────────────────────────────────────

void ColStoreWriter::set_rows(std::vector<uint64_t> row_keys) {
    for (size_t i = 1; i < row_keys.size(); ++i)
        if (row_keys[i] <= row_keys[i - 1])
            throw std::runtime_error("ColStoreWriter::set_rows: keys must be strictly ascending");
    rows_ = std::move(row_keys);
}

uint32_t ColStoreWriter::intern(const std::string& s) {
    if (s.empty()) return 0;  // offset 0 is the canonical empty string
    auto off = static_cast<uint32_t>(strtab_.size());
    strtab_.append(s);
    strtab_.push_back('\0');
    return off;
}

void ColStoreWriter::add_column_raw(const ColumnKey& key, ColDType dtype, const void* data,
                                    uint64_t n_values, const std::vector<uint8_t>& present) {
    if (n_values != rows_.size())
        throw std::runtime_error("ColStoreWriter::add_column: value count != row count");
    const uint64_t sz = dtype_size(dtype);
    if (sz == 0) throw std::runtime_error("ColStoreWriter::add_column: invalid dtype");

    PendingColumn col;
    col.key   = key;
    col.dtype = dtype;
    col.values.resize(static_cast<size_t>(n_values * sz));
    if (n_values) std::memcpy(col.values.data(), data, col.values.size());

    if (!present.empty()) {
        if (present.size() != n_values)
            throw std::runtime_error("ColStoreWriter::add_column: present mask size mismatch");
        bool any_absent = false;
        for (uint8_t p : present) if (!p) { any_absent = true; break; }
        if (any_absent) {
            col.has_nulls = true;
            col.nullbits.assign(static_cast<size_t>((n_values + 7) / 8), 0);
            for (uint64_t i = 0; i < n_values; ++i)
                if (present[i]) col.nullbits[i >> 3] |= static_cast<uint8_t>(1u << (i & 7));
        }
    }
    columns_.push_back(std::move(col));
}

SectionDesc ColStoreWriter::finalize(AppendWriter& w, uint64_t section_id) {
    const auto n_rows    = static_cast<uint32_t>(rows_.size());
    const auto n_columns = static_cast<uint32_t>(columns_.size());

    w.align(8);
    const uint64_t section_start = w.current_offset();

    // 1) placeholder header (patched at the end once offsets are known)
    ColStoreHeader hdr{};
    w.append(&hdr, sizeof(hdr));

    // 2) row keys
    w.align(8);
    const uint64_t rowkey_off = w.current_offset() - section_start;
    if (n_rows) w.append(rows_.data(), static_cast<uint64_t>(n_rows) * sizeof(uint64_t));

    // 3) column bodies (values then optional null bitmap), building the directory
    std::vector<ColumnDesc> descs(n_columns);
    for (uint32_t c = 0; c < n_columns; ++c) {
        PendingColumn& pc = columns_[c];
        w.align(8);
        const uint64_t data_off = w.current_offset() - section_start;
        if (!pc.values.empty()) w.append(pc.values.data(), pc.values.size());
        uint64_t nullbits_off = 0;
        if (pc.has_nulls) {
            nullbits_off = w.current_offset() - section_start;
            w.append(pc.nullbits.data(), pc.nullbits.size());
        }

        ColumnDesc& d = descs[c];
        d.col_identity_hash = column_identity_hash(pc.key);
        d.ref_db_hash       = pc.key.ref_db_hash;
        d.core_model_hash   = pc.key.core_model_hash;
        d.calibration_hash  = pc.key.calibration_hash;
        d.data_off          = data_off;
        d.data_len          = pc.values.size();
        d.nullbits_off      = nullbits_off;
        d.axis_name_off     = intern(pc.key.axis);
        d.tool_name_off     = intern(pc.key.tool);
        d.method_name_off   = intern(pc.key.method);
        d.label_name_off    = intern(pc.key.label);
        d.data_checksum32   = checksum32(pc.values.data(), pc.values.size(),
                                         pc.nullbits.data(), pc.nullbits.size());
        d.unit              = static_cast<uint16_t>(pc.key.unit);
        d.dtype             = static_cast<uint16_t>(pc.dtype);
        d.encoding          = static_cast<uint16_t>(ColEncoding::Raw);
        d.col_version       = pc.key.version;
        d.col_flags         = pc.has_nulls ? COL_HAS_NULLS : uint16_t{0};
    }

    // 4) column directory
    w.align(8);
    const uint64_t coldir_off = w.current_offset() - section_start;
    if (n_columns) w.append(descs.data(), static_cast<uint64_t>(n_columns) * sizeof(ColumnDesc));

    // 5) string table
    const uint64_t strtab_off = w.current_offset() - section_start;
    w.append(strtab_.data(), strtab_.size());

    const uint64_t section_end = w.current_offset();

    // 6) patch the header now that every offset is known
    hdr.magic          = CSTORE_MAGIC;
    hdr.format_version = CSTORE_FORMAT_VERSION;
    hdr.flags          = 0;
    hdr.n_rows         = n_rows;
    hdr.n_columns      = n_columns;
    hdr.rowkey_off     = rowkey_off;
    hdr.coldir_off     = coldir_off;
    hdr.strtab_off     = strtab_off;
    hdr.strtab_len     = strtab_.size();
    w.write_at(section_start, &hdr, sizeof(hdr));

    SectionDesc sd{};
    sd.type              = section_type_;
    sd.version           = CSTORE_FORMAT_VERSION;
    sd.flags             = 0;
    sd.section_id        = section_id;
    sd.file_offset       = section_start;
    sd.compressed_size   = section_end - section_start;
    sd.uncompressed_size = section_end - section_start;
    sd.item_count        = n_rows;
    sd.aux0              = n_columns;
    sd.header_size       = sizeof(ColStoreHeader);
    return sd;
}

// ── ColStoreReader ────────────────────────────────────────────────────────────

void ColStoreReader::open(const uint8_t* data, uint64_t offset, uint64_t size) {
    if (size < sizeof(ColStoreHeader))
        throw std::runtime_error("ColStore section too small");
    data_   = data + offset;
    header_ = reinterpret_cast<const ColStoreHeader*>(data_);
    if (header_->magic != CSTORE_MAGIC)
        throw std::runtime_error("ColStore: bad magic");
    if (header_->format_version > CSTORE_FORMAT_VERSION)
        throw std::runtime_error("ColStore: format version "
                                 + std::to_string(header_->format_version) + " too new");

    const uint64_t n_rows    = header_->n_rows;
    const uint64_t n_columns = header_->n_columns;
    const uint64_t rk_end    = header_->rowkey_off + n_rows * sizeof(uint64_t);
    const uint64_t cd_end    = header_->coldir_off + n_columns * sizeof(ColumnDesc);
    const uint64_t st_end    = header_->strtab_off + header_->strtab_len;
    if (rk_end > size || cd_end > size || st_end > size)
        throw std::runtime_error("ColStore section truncated");

    rowkeys_    = reinterpret_cast<const uint64_t*>(data_ + header_->rowkey_off);
    coldir_     = reinterpret_cast<const ColumnDesc*>(data_ + header_->coldir_off);
    strtab_     = reinterpret_cast<const char*>(data_ + header_->strtab_off);
    strtab_len_ = header_->strtab_len;

    // bounds-check every column body
    for (uint64_t c = 0; c < n_columns; ++c) {
        const ColumnDesc& d = coldir_[c];
        if (d.data_off + d.data_len > size)
            throw std::runtime_error("ColStore: column data out of bounds");
        if (d.col_flags & COL_HAS_NULLS) {
            const uint64_t nb = (n_rows + 7) / 8;
            if (d.nullbits_off + nb > size)
                throw std::runtime_error("ColStore: column null bitmap out of bounds");
        }
    }
}

std::string_view ColStoreReader::str_at(uint32_t off) const {
    if (off == 0 || off >= strtab_len_) return {};
    const char* p = strtab_ + off;
    size_t maxlen = static_cast<size_t>(strtab_len_ - off);
    size_t len = ::strnlen(p, maxlen);
    return {p, len};
}

int64_t ColStoreReader::row_index(uint64_t key) const {
    if (!header_ || header_->n_rows == 0) return -1;
    const uint64_t* end = rowkeys_ + header_->n_rows;
    const uint64_t* it = std::lower_bound(rowkeys_, end, key);
    if (it != end && *it == key) return it - rowkeys_;
    return -1;
}

ColStoreReader::ColumnView ColStoreReader::column(uint32_t c) const {
    const ColumnDesc& d = coldir_[c];
    ColumnView v;
    v.identity_hash    = d.col_identity_hash;
    v.axis             = str_at(d.axis_name_off);
    v.tool             = str_at(d.tool_name_off);
    v.method           = str_at(d.method_name_off);
    v.label            = str_at(d.label_name_off);
    v.unit             = static_cast<Unit>(d.unit);
    v.dtype            = static_cast<ColDType>(d.dtype);
    v.version          = d.col_version;
    v.ref_db_hash      = d.ref_db_hash;
    v.core_model_hash  = d.core_model_hash;
    v.calibration_hash = d.calibration_hash;
    v.n_rows           = header_->n_rows;
    v.values           = data_ + d.data_off;
    v.nullbits         = (d.col_flags & COL_HAS_NULLS) ? data_ + d.nullbits_off : nullptr;
    return v;
}

int ColStoreReader::find_column(uint64_t identity_hash) const {
    if (!header_) return -1;
    for (uint32_t c = 0; c < header_->n_columns; ++c)
        if (coldir_[c].col_identity_hash == identity_hash) return static_cast<int>(c);
    return -1;
}

bool ColStoreReader::verify_checksums() const {
    if (!header_) return false;
    const uint64_t nb = (static_cast<uint64_t>(header_->n_rows) + 7) / 8;
    for (uint32_t c = 0; c < header_->n_columns; ++c) {
        const ColumnDesc& d = coldir_[c];
        const uint8_t* nbp = (d.col_flags & COL_HAS_NULLS) ? data_ + d.nullbits_off : nullptr;
        uint32_t got = checksum32(data_ + d.data_off, d.data_len, nbp, nbp ? nb : 0);
        if (got != d.data_checksum32) return false;
    }
    return true;
}

double ColStoreReader::ColumnView::get_f64(uint32_t i) const {
    if (!present(i)) return std::nan("");
    switch (dtype) {
        case ColDType::F32: return get<float>(i);
        case ColDType::F64: return get<double>(i);
        case ColDType::U8:  return get<uint8_t>(i);
        case ColDType::U16: return get<uint16_t>(i);
        case ColDType::U32: return get<uint32_t>(i);
        case ColDType::U64: return static_cast<double>(get<uint64_t>(i));
        case ColDType::I8:  return get<int8_t>(i);
        case ColDType::I16: return get<int16_t>(i);
        case ColDType::I32: return get<int32_t>(i);
        case ColDType::I64: return static_cast<double>(get<int64_t>(i));
        default:            return std::nan("");
    }
}

} // namespace genopack
