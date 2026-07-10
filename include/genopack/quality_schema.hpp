#pragma once
// ── Quality schema: provenance-first column identity ──────────────────────────
// Every quality measurement is a cell (genome, axis, source/tool, method). A
// column groups all cells that share one (axis, tool, method, unit, version,
// reference-db, core-model, calibration) identity. That identity is a content
// hash — numeric ids in the on-disk store are archive-local interning only, so a
// merge resolves columns by hash, never by raw id (raw ids collide across
// writers). This header defines the enums, the ColumnKey, and the hash; the
// physical columnar container lives in colstore.hpp.
#include <cstdint>
#include <string>
#include <string_view>
#include <xxhash.h>

namespace genopack {

// Semantic unit of a column's values. Replaces an ad-hoc per-column scale float:
// the unit tells a consumer how to interpret/format a value without a magic
// scale factor. Append-only — never renumber.
enum class Unit : uint16_t {
    Unknown   = 0,
    Fraction01 = 1,   // [0,1] proportion (completeness, contamination fraction)
    Percent    = 2,   // [0,100]
    Count      = 3,   // non-negative integer count
    Bits       = 4,   // information / entropy in bits
    ZScore     = 5,   // standardized score
    LogOdds    = 6,   // log-odds / logit
    Boolean    = 7,   // 0/1 flag
    Ratio      = 8,   // unbounded non-negative ratio
};

// Physical element type of a column. Append-only — never renumber.
enum class ColDType : uint16_t {
    None = 0,
    F32  = 1,
    F64  = 2,
    U8   = 3,
    U16  = 4,
    U32  = 5,
    U64  = 6,
    I8   = 7,
    I16  = 8,
    I32  = 9,
    I64  = 10,
};

// Column value encoding. Raw = values stored verbatim in their dtype. The
// columnar store is lossless: quantized u8/u16 encodings of the legacy POD are
// not reused; presence is an explicit null bitmap, not a value sentinel.
enum class ColEncoding : uint16_t {
    Raw = 0,
};

inline uint64_t dtype_size(ColDType t) {
    switch (t) {
        case ColDType::F32: case ColDType::U32: case ColDType::I32: return 4;
        case ColDType::F64: case ColDType::U64: case ColDType::I64: return 8;
        case ColDType::U8:  case ColDType::I8:                      return 1;
        case ColDType::U16: case ColDType::I16:                     return 2;
        default:                                                    return 0;
    }
}

// ── Canonical names ───────────────────────────────────────────────────────────
// Producers should use these constants so identity hashes match across the
// codebase. The on-disk store records whatever strings the producer passed; two
// columns collapse on merge iff their (axis,tool,method,unit,version,...) match
// byte-for-byte. New names may be added freely; do NOT rename an existing one
// (that changes the identity hash of every column that used it).

namespace qual_tool {
inline constexpr const char* GENOPACK = "genopack";  // intrinsic, computed in-archive
inline constexpr const char* CHECKM2  = "checkm2";
inline constexpr const char* ANVIO    = "anvio";
}

namespace qual_axis {
inline constexpr const char* COMPLETENESS   = "completeness";
inline constexpr const char* CONTAMINATION  = "contamination";
inline constexpr const char* COHERENCE      = "coherence";
inline constexpr const char* CHROM_SKEW     = "chromosome_skew";
inline constexpr const char* SUPPORT        = "support";
inline constexpr const char* INTERVAL       = "interval";
}

namespace qual_method {
// completeness methods
inline constexpr const char* AAMER_GENUS_CORE = "aamer_genus_core"; // genus prevalence core
inline constexpr const char* CLUSTER_RELATIVE = "cluster_relative"; // OPH containment vs genus
inline constexpr const char* SKETCH_FILL      = "sketch_fill";
inline constexpr const char* FRAGMENTATION    = "fragmentation";
inline constexpr const char* POST_DECONTAM    = "post_decontam";
inline constexpr const char* MARKER_SCG       = "marker_scg";       // GTDB bac120/ar53
inline constexpr const char* ML               = "ml";              // CheckM2 gradient-boost/NN
// contamination methods
inline constexpr const char* DUPLICATION      = "duplication";     // redundancy_fraction
inline constexpr const char* CORE_DUP_MASS    = "core_dup_mass";   // Phase-2: non-saturating SCC dup mass Σ(c-1)/Σc
inline constexpr const char* ACCESSORY_RATIO  = "accessory_ratio"; // Phase-2: c0_query / median(c0_all)
inline constexpr const char* LEAKAGE          = "leakage";
inline constexpr const char* TNF_EXCESS       = "tnf_excess";
inline constexpr const char* MIXTURE          = "mixture";
inline constexpr const char* MIXTURE_SOURCES  = "mixture_sources";
inline constexpr const char* MIXTURE_WINDOWS  = "mixture_windows";
inline constexpr const char* FLAGS            = "flags";
inline constexpr const char* SPECTRAL         = "spectral";        // fiedler
inline constexpr const char* CONTIG_OUTLIER   = "contig_outlier";
inline constexpr const char* SPE_OUTLIER      = "spe_outlier";
inline constexpr const char* SIBLING_OUTLIER  = "sibling_outlier";
inline constexpr const char* RHO_OUTLIER      = "rho_outlier";
inline constexpr const char* CROSS_GENUS      = "cross_genus";
inline constexpr const char* FMH_MINORITY     = "fmh_minority";
inline constexpr const char* MARKER_REDUNDANCY = "marker_redundancy";
// coherence / skew / support
inline constexpr const char* SELF_COHERENCE   = "self_coherence";
inline constexpr const char* CHARGAFF_PARITY  = "chargaff_parity";
inline constexpr const char* SPECTRAL_GAP     = "spectral_gap";
inline constexpr const char* SCALE_KINK       = "scale_kink";
inline constexpr const char* LEAKAGE_RESIDUAL = "leakage_residual";
inline constexpr const char* SKEW_CLOSURE     = "skew_closure";
inline constexpr const char* TIER             = "tier";
inline constexpr const char* QUALITY_TIER     = "quality_tier";
inline constexpr const char* WIDTH            = "width";
}

// ── Column identity ───────────────────────────────────────────────────────────
// The provenance tuple that uniquely names a column. `label` is a human-facing
// display string only and is NOT part of identity.
struct ColumnKey {
    std::string axis;
    std::string tool;
    std::string method;
    Unit        unit             = Unit::Unknown;
    uint16_t    version          = 1;
    uint64_t    ref_db_hash      = 0;  // reference DB / panel identity (0 = N/A)
    uint64_t    core_model_hash  = 0;  // CORE per-genus model identity (0 = N/A)
    uint64_t    calibration_hash = 0;  // calibration model identity (0 = uncalibrated)
    std::string label;
};

// Content hash of a column's identity. Canonical serialization: the three name
// strings NUL-delimited, then the scalar fields in fixed little-endian order.
// `label` is excluded. Stable across writers and runs.
inline uint64_t column_identity_hash(const ColumnKey& k) {
    std::string buf;
    buf.reserve(k.axis.size() + k.tool.size() + k.method.size() + 40);
    buf.append(k.axis);   buf.push_back('\0');
    buf.append(k.tool);   buf.push_back('\0');
    buf.append(k.method); buf.push_back('\0');
    auto put = [&](uint64_t v, int n) {
        for (int i = 0; i < n; ++i) buf.push_back(static_cast<char>((v >> (8 * i)) & 0xFF));
    };
    put(static_cast<uint64_t>(k.unit), 2);
    put(k.version, 2);
    put(k.ref_db_hash, 8);
    put(k.core_model_hash, 8);
    put(k.calibration_hash, 8);
    return XXH3_64bits(buf.data(), buf.size());
}

} // namespace genopack
