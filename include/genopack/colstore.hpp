#pragma once
// ── ColumnStore: self-describing columnar section ─────────────────────────────
// A fixed-width-per-row columnar container used by the provenance-first quality
// sections (intrinsic QCOL, external XQAL, calibrated CQAL). Rows are keyed by a
// sorted uint64 row key (genome_id for quality). Each column is an independent
// typed array of n_rows values plus an optional null bitmap; a column carries
// its full provenance identity (axis/tool/method/unit/version + ref/core/calib
// hashes) so columns from different archives merge by content hash, not by
// position. The container is generic over the TOC section type (one writer class
// emits QCOL/XQAL/CQAL by passing the section magic), mirroring how GcovWriter
// emits both GCOV and FCOV.
//
// On-disk layout (within the section, all offsets from section start):
//   [ColStoreHeader 64B]
//   [row keys: n_rows * u64]                 (8-aligned)
//   [column 0 values][column 0 null bitmap]  (values 8-aligned; bitmap follows)
//   [column 1 values][column 1 null bitmap]
//   ...
//   [ColumnDesc[n_columns]]                  (8-aligned)
//   [string table]                           (NUL-delimited; offset 0 == "")
#include "format.hpp"
#include "mmap_file.hpp"
#include "quality_schema.hpp"
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace genopack {

static constexpr uint32_t CSTORE_MAGIC          = 0x52545343u; // "CSTR"
static constexpr uint16_t CSTORE_FORMAT_VERSION = 1;

// ── On-disk preamble ──────────────────────────────────────────────────────────
struct ColStoreHeader {              // 64 bytes
    uint32_t magic;                  //  0  CSTORE_MAGIC
    uint16_t format_version;         //  4  CSTORE_FORMAT_VERSION
    uint16_t flags;                  //  6
    uint32_t n_rows;                 //  8
    uint32_t n_columns;              // 12
    uint64_t rowkey_off;             // 16  -> n_rows * u64 (sorted ascending)
    uint64_t coldir_off;             // 24  -> ColumnDesc[n_columns]
    uint64_t strtab_off;             // 32  -> string table
    uint64_t strtab_len;             // 40
    uint64_t _reserved0;             // 48
    uint64_t _reserved1;             // 56
};
static_assert(sizeof(ColStoreHeader) == 64);

// Per-column directory entry. Self-describing: a reader can recompute and verify
// the identity hash from the resolved name strings + scalar fields.
struct ColumnDesc {                  // 96 bytes
    uint64_t col_identity_hash;      //  0  content hash of the provenance tuple
    uint64_t ref_db_hash;            //  8
    uint64_t core_model_hash;        // 16
    uint64_t calibration_hash;       // 24
    uint64_t data_off;              // 32  value array start (8-aligned)
    uint64_t data_len;              // 40  value bytes == n_rows * dtype_size
    uint64_t nullbits_off;          // 48  null bitmap start (0 if no nulls)
    uint32_t axis_name_off;         // 56  strtab offset (>=1)
    uint32_t tool_name_off;         // 60
    uint32_t method_name_off;       // 64
    uint32_t label_name_off;        // 68  0 == none
    uint32_t data_checksum32;       // 72  low32 XXH3_64 of [values || nullbits]
    uint16_t unit;                  // 76  Unit
    uint16_t dtype;                 // 78  ColDType
    uint16_t encoding;              // 80  ColEncoding
    uint16_t col_version;           // 82
    uint16_t col_flags;             // 84  bit0 = COL_HAS_NULLS
    uint16_t _pad0;                 // 86
    uint64_t _reserved;             // 88
};
static_assert(sizeof(ColumnDesc) == 96);
static_assert(offsetof(ColumnDesc, data_off) == 32);

static constexpr uint16_t COL_HAS_NULLS = 0x1;

// ── Writer ────────────────────────────────────────────────────────────────────
// dtype_of<T> maps a C++ scalar to its ColDType.
template <class T> struct col_dtype_of;
template <> struct col_dtype_of<float>    { static constexpr ColDType v = ColDType::F32; };
template <> struct col_dtype_of<double>   { static constexpr ColDType v = ColDType::F64; };
template <> struct col_dtype_of<uint8_t>  { static constexpr ColDType v = ColDType::U8;  };
template <> struct col_dtype_of<uint16_t> { static constexpr ColDType v = ColDType::U16; };
template <> struct col_dtype_of<uint32_t> { static constexpr ColDType v = ColDType::U32; };
template <> struct col_dtype_of<uint64_t> { static constexpr ColDType v = ColDType::U64; };
template <> struct col_dtype_of<int8_t>   { static constexpr ColDType v = ColDType::I8;  };
template <> struct col_dtype_of<int16_t>  { static constexpr ColDType v = ColDType::I16; };
template <> struct col_dtype_of<int32_t>  { static constexpr ColDType v = ColDType::I32; };
template <> struct col_dtype_of<int64_t>  { static constexpr ColDType v = ColDType::I64; };

class ColStoreWriter {
public:
    explicit ColStoreWriter(uint32_t section_type) : section_type_(section_type) {}

    // Set the row keys (genome_ids). Must be strictly ascending; every column's
    // values[i] corresponds to row key i. Call once before adding columns.
    void set_rows(std::vector<uint64_t> row_keys);

    // Add a typed column. `values.size()` must equal n_rows. `present` is an
    // optional per-row mask (1 = present, 0 = null); empty = all present. If all
    // rows are present no null bitmap is written.
    template <class T>
    void add_column(const ColumnKey& key, const std::vector<T>& values,
                    const std::vector<uint8_t>& present = {}) {
        add_column_raw(key, col_dtype_of<T>::v, values.data(),
                       static_cast<uint64_t>(values.size()), present);
    }

    // Type-erased column add (dtype given explicitly). `data` points at
    // n_values elements of `dtype`.
    void add_column_raw(const ColumnKey& key, ColDType dtype, const void* data,
                        uint64_t n_values, const std::vector<uint8_t>& present = {});

    // Emit the section through `w` and return its descriptor (checksum /
    // derivation_hash / semantic_schema_hash stamped later by the caller).
    SectionDesc finalize(AppendWriter& w, uint64_t section_id);

    uint32_t n_rows() const { return static_cast<uint32_t>(rows_.size()); }
    size_t   n_columns() const { return columns_.size(); }

private:
    struct PendingColumn {
        ColumnKey key;
        ColDType  dtype;
        std::vector<uint8_t> values;    // raw value bytes
        std::vector<uint8_t> nullbits;  // empty if all present
        bool has_nulls = false;
    };
    uint32_t intern(const std::string& s);  // -> strtab offset

    uint32_t section_type_;
    std::vector<uint64_t>      rows_;
    std::vector<PendingColumn> columns_;
    std::string                strtab_{std::string(1, '\0')};  // offset 0 == ""
};

// ── Reader ────────────────────────────────────────────────────────────────────
class ColStoreReader {
public:
    struct ColumnView {
        uint64_t         identity_hash = 0;
        std::string_view axis, tool, method, label;
        Unit             unit = Unit::Unknown;
        ColDType         dtype = ColDType::None;
        uint16_t         version = 0;
        uint64_t         ref_db_hash = 0, core_model_hash = 0, calibration_hash = 0;
        uint32_t         n_rows = 0;
        const uint8_t*   values = nullptr;    // n_rows * dtype_size bytes
        const uint8_t*   nullbits = nullptr;  // nullptr -> all present

        bool present(uint32_t i) const {
            return !nullbits || (nullbits[i >> 3] & (1u << (i & 7)));
        }
        // Numeric value at row i widened to double; NaN if absent.
        double get_f64(uint32_t i) const;
        template <class T> T get(uint32_t i) const {
            return reinterpret_cast<const T*>(values)[i];
        }
    };

    void open(const uint8_t* data, uint64_t offset, uint64_t size);
    bool is_open() const { return data_ != nullptr; }

    uint32_t n_rows() const { return header_ ? header_->n_rows : 0; }
    uint32_t n_columns() const { return header_ ? header_->n_columns : 0; }

    uint64_t row_key(uint32_t i) const { return rowkeys_[i]; }
    // Binary search the sorted row key; -1 if absent.
    int64_t  row_index(uint64_t key) const;

    ColumnView column(uint32_t c) const;
    // Find a column by identity hash / key; -1 if absent.
    int  find_column(uint64_t identity_hash) const;
    int  find_column(const ColumnKey& key) const { return find_column(column_identity_hash(key)); }

    // Verify every column's stored data_checksum32. Returns false on mismatch.
    bool verify_checksums() const;

private:
    std::string_view str_at(uint32_t off) const;

    const uint8_t*        data_     = nullptr;
    const ColStoreHeader* header_   = nullptr;
    const uint64_t*       rowkeys_  = nullptr;
    const ColumnDesc*     coldir_   = nullptr;
    const char*           strtab_   = nullptr;
    uint64_t              strtab_len_ = 0;
};

} // namespace genopack
