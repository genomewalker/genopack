#pragma once
#include "format.hpp"
#include "types.hpp"
#include "mmap_file.hpp"
#include <cstdint>
#include <cmath>
#include <vector>

namespace genopack {

// On-disk per-genome quality record (80 bytes, POD).
// Layout is fixed; add fields by shrinking _reserved only.
struct QualRecord {
    uint64_t genome_id;                     //  8
    float    completeness_cluster_relative; //  4
    float    completeness_fragmentation;    //  4
    float    completeness_post_decontam;    //  4
    float    contamination_leakage;         //  4
    float    contamination_tnf_excess;      //  4
    float    chromosome_skew_closure;       //  4
    float    leakage_residual;              //  4  per-DoF RMSE of slope-only log-linear fit
    uint8_t  support_tier;                  //  1
    uint8_t  spe_outlier_u8;               //  1  contamination_spe × 255; 0=clean
    uint8_t  sibling_outlier_u8;           //  1  contamination_sibling_outlier × 255; 0=clean
    uint8_t  rho_outlier_u8;              //  1  contamination_rho_outlier × 255; 0=clean
    float    interval_width;                //  4
    float    self_coherence;               //  4  (NAN = not computed)
    float    contamination_mixture;         //  4  (NAN = not computed)
    int16_t  mixture_sources;              //  2  (0 = not computed)
    uint16_t fiedler_u16;                  //  2  contamination_spectral × 65535; 0=pure
    uint16_t n_mix_windows;               //  2  windows used in mixture model
    uint8_t  qual_flags;                  //  1  bitfield — see QUAL_FLAG_* constants
    uint8_t  contig_outlier_u8;          //  1  contamination_contig_outlier × 255; 0=pure
    uint8_t  fmh_minority_u8;             //  1  LEGACY: fmh_minority_fraction × 255 (truncated); prefer fmh_minority_u16
    uint8_t  marker_completeness_u8;     //  1  marker_completeness encoded: 0=not_scored, 1-255=(v-1)/254.0
    uint16_t marker_redundancy_u16;      //  2  raw marker_redundancy × 65535 (0xFFFF = not scored)
    float    chargaff_parity;            //  4  Chargaff 2nd parity score [0,1]; NAN = not computed
    float    spectral_gap;               //  4  1 - |λ₂| from k3 Markov op; NAN = not computed
    float    scale_kink;                 //  4  log₂(W*) dyadic-window scale; NAN = not computed
    uint8_t  cross_genus_u8;             //  1  contamination_cross_genus × 255; 0=clean
    uint8_t  sketch_fill_u8;             //  1  completeness_sketch_fill × 200 (200=100%, >200 allowed up to 255)
    uint16_t contamination_duplication_u16; //  2  redundancy_fraction: 0=not scored, else round(f*65534)+1
    // total = 80 (v3 layout — kV3Stride)
    float    completeness_aamer_core;        //  4  CORE aamer coverage [0,1]; NAN = not scored
    float    completeness_aamer_family_core; //  4  FCORE aamer coverage [0,1]; NAN = not scored
    uint16_t fmh_minority_u16;              //  2  fmh_minority_fraction × 65535 (rounded); 0=not scored/clean
    uint16_t contig_outlier_u16;            //  2  contamination_contig_outlier × 65535 (rounded); 0=not scored/clean
    uint8_t  quality_tier_u8;               //  1  LEGACY discrete tier; 0=not_set,1=LQ,2=MQ,3=HQ (kept for back-compat)
    uint8_t  qscore_lo;                     //  1  low byte of continuous quality ∈ [0,1] × 65534 + 1; 0 = not scored
    uint8_t  qscore_hi;                     //  1  high byte (split into two u8 to keep quality_tier_u8 at a fixed offset)
    uint8_t  _reserved[1];                  //  1  reserved for future fields
    // total = 96 (kV4Stride — old readers stop here)
    float    contamination_core_dup_mass;   //  4  Phase-2: non-saturating SCC dup mass Σ(c-1)/Σc; NAN = not scored
    float    accessory_ratio;               //  4  Phase-2: c0_query / median(c0_all); NAN = not scored
    // total = 104

    // contamination_duplication encode/decode (0 reserved as not-scored sentinel,
    // so a genuinely-clean 0.0 stays distinguishable from an unscored genome).
    static uint16_t encode_dup(float f) {
        if (!(f >= 0.0f)) return 0;            // NaN / negative -> not scored
        if (f > 1.0f) f = 1.0f;
        return static_cast<uint16_t>(f * 65534.0f + 0.5f) + 1;
    }
    static float decode_dup(uint16_t v) {
        return v == 0 ? NAN : static_cast<float>(v - 1) / 65534.0f;
    }

    // qual_flags bits
    static constexpr uint8_t QUAL_FLAG_MIX_NO_DATA        = 0x01; // mixture model had < min_windows
    static constexpr uint8_t QUAL_FLAG_FRAG_GATED         = 0x02; // completeness_fragmentation < 0.65
    static constexpr uint8_t QUAL_FLAG_LEAKAGE_UNRELIABLE = 0x04; // leakage_residual above threshold
    static constexpr uint8_t QUAL_FLAG_COHERENCE_GATED    = 0x08; // self_coherence < 0.80
    static constexpr uint8_t QUAL_FLAG_MARKER_SCORED      = 0x10; // marker_redundancy_u16 is valid
    static constexpr uint8_t QUAL_FLAG_BUILD_SIGNALS      = 0x20; // AllSignals computed at build time (chargaff/spectral/scale_kink valid)
    static constexpr uint8_t QUAL_FLAG_GCOV_SCORED        = 0x40; // per-contig GCOV scoring cached (contig_outlier/spe/sibling/rho valid)
    static constexpr uint8_t QUAL_FLAG_FMH_AXIS_ABSENT    = 0x80; // FMH minority axis unavailable -> contamination tier used the sketch-leakage fallback

    // quality_tier_u8 values (0 = not written, backward-compatible with old packs)
    static constexpr uint8_t QTIER_NOT_SET = 0;
    static constexpr uint8_t QTIER_LQ      = 1;
    static constexpr uint8_t QTIER_MQ      = 2;
    static constexpr uint8_t QTIER_HQ      = 3;

    static uint8_t encode_qtier(const char* tier) {
        if (tier[0] == 'H') return QTIER_HQ;
        if (tier[0] == 'M') return QTIER_MQ;
        if (tier[0] == 'L') return QTIER_LQ;
        return QTIER_NOT_SET;
    }

    // Continuous quality score encode/decode (0 reserved as not-scored sentinel,
    // so a genuine 0.0 stays distinguishable from an unscored genome). Stored across
    // qscore_lo/qscore_hi so quality_tier_u8 keeps its fixed on-disk offset.
    static uint16_t encode_qscore(float s) {
        if (!(s >= 0.0f)) return 0;            // NaN / negative -> not scored
        if (s > 1.0f) s = 1.0f;
        return static_cast<uint16_t>(s * 65534.0f + 0.5f) + 1;
    }
    static float decode_qscore(uint16_t v) {
        return v == 0 ? NAN : static_cast<float>(v - 1) / 65534.0f;
    }
    void set_quality_score(float s) {
        const uint16_t v = encode_qscore(s);
        qscore_lo = static_cast<uint8_t>(v & 0xFF);
        qscore_hi = static_cast<uint8_t>(v >> 8);
    }
    float quality_score() const {
        return decode_qscore(static_cast<uint16_t>(qscore_lo) |
                             (static_cast<uint16_t>(qscore_hi) << 8));
    }

    // support_tier sentinel: no genus-level reference available (singleton genus)
    static constexpr uint8_t TIER_NO_GENUS_REF = 0xFF;

    static QualRecord make_empty(uint64_t gid) {
        QualRecord r{};
        r.genome_id                     = gid;
        r.completeness_cluster_relative = NAN;
        r.completeness_fragmentation    = NAN;
        r.completeness_post_decontam    = NAN;
        r.contamination_leakage         = NAN;
        r.contamination_tnf_excess      = NAN;
        r.chromosome_skew_closure       = NAN;
        r.leakage_residual              = NAN;
        r.interval_width                = NAN;
        r.self_coherence               = NAN;
        r.contamination_mixture        = NAN;
        r.mixture_sources              = 0;
        r.fiedler_u16                  = 0;
        r.n_mix_windows               = 0;
        r.spe_outlier_u8              = 0;
        r.sibling_outlier_u8          = 0;
        r.rho_outlier_u8              = 0;
        r.cross_genus_u8              = 0;
        r.sketch_fill_u8              = 0;
        r.marker_redundancy_u16       = 0xFFFF;
        r.chargaff_parity             = NAN;
        r.spectral_gap                = NAN;
        r.scale_kink                  = NAN;
        r.qual_flags                  = 0;
        r.fmh_minority_u8             = 0;
        r.fmh_minority_u16            = 0;
        r.contig_outlier_u16          = 0;
        r.marker_completeness_u8      = 0;
        r.contamination_duplication_u16 = 0; // 0 = not scored
        r.completeness_aamer_core        = NAN;
        r.completeness_aamer_family_core = NAN;
        r.quality_tier_u8                = QTIER_NOT_SET;
        r.qscore_lo                      = 0; // 0 = not scored
        r.qscore_hi                      = 0;
        r.contamination_core_dup_mass    = NAN; // Phase-2: default not-scored (also the old-pack upcast default)
        r.accessory_ratio                = NAN;
        return r;
    }
};
static_assert(sizeof(QualRecord) == 104);

// ── Writer ────────────────────────────────────────────────────────────────────

struct QualHeader {            // 32 bytes
    uint32_t magic;            // SEC_QUAL
    uint32_t n_records;
    uint64_t records_offset;   // bytes from section start
    uint64_t record_stride;    // sizeof(QualRecord)
    uint8_t  reserved[8];
};
static_assert(sizeof(QualHeader) == 32);

class QualWriter {
public:
    void add(const QualRecord& r) { records_.push_back(r); }

    SectionDesc finalize(AppendWriter& w, uint64_t section_id) {
        const uint32_t n = static_cast<uint32_t>(records_.size());
        w.align(8);
        const uint64_t section_start = w.current_offset();

        QualHeader hdr{};
        hdr.magic          = SEC_QUAL;
        hdr.n_records      = n;
        hdr.records_offset = sizeof(QualHeader);
        hdr.record_stride  = sizeof(QualRecord);
        w.append(&hdr, sizeof(hdr));
        if (n > 0)
            w.append(records_.data(), static_cast<uint64_t>(n) * sizeof(QualRecord));

        const uint64_t section_end = w.current_offset();
        SectionDesc sd{};
        sd.type              = SEC_QUAL;
        sd.version           = 1;
        sd.flags             = 0;
        sd.section_id        = section_id;
        sd.file_offset       = section_start;
        sd.compressed_size   = section_end - section_start;
        sd.uncompressed_size = section_end - section_start;
        sd.item_count        = n;
        return sd;
    }

    size_t size() const { return records_.size(); }

    const std::vector<QualRecord>& records() const { return records_; }
    std::vector<QualRecord> take() { return std::move(records_); }

private:
    std::vector<QualRecord> records_;
};

// ── Reader ────────────────────────────────────────────────────────────────────

class QualReader {
public:
    static constexpr uint64_t kOldStride    = 48; // pre-v2 layout without contamination_mixture
    static constexpr uint64_t kMediumStride = 64; // pre-v3 layout without chargaff/spectral/scale_kink
    static constexpr uint64_t kV3Stride     = 80; // v3 layout without completeness_aamer_core
    static constexpr uint64_t kV4Stride     = 96; // v4 layout without core_dup_mass/accessory_ratio

    void open(const uint8_t* data, uint64_t offset, uint64_t size) {
        if (size < sizeof(QualHeader))
            throw std::runtime_error("QUAL section too small");
        data_   = data + offset;
        header_ = reinterpret_cast<const QualHeader*>(data_);
        if (header_->magic != SEC_QUAL)
            throw std::runtime_error("QUAL: bad magic");

        const uint64_t stride = header_->record_stride;
        if (stride != sizeof(QualRecord) && stride != kOldStride && stride != kMediumStride
            && stride != kV3Stride && stride != kV4Stride)
            throw std::runtime_error("QUAL: unknown record_stride " + std::to_string(stride)
                                     + " — rebuild required");
        old_layout_    = (stride == kOldStride);
        medium_layout_ = (stride == kMediumStride);
        v3_layout_     = (stride == kV3Stride);
        v4_layout_     = (stride == kV4Stride);

        const uint64_t end = header_->records_offset
            + static_cast<uint64_t>(header_->n_records) * stride;
        if (end > size)
            throw std::runtime_error("QUAL section truncated");
        base_ = data_ + header_->records_offset;
    }

    bool is_open() const { return data_ != nullptr; }
    uint32_t n_records() const { return header_ ? header_->n_records : 0; }
    bool is_old_layout() const { return old_layout_; }

    // Linear scan. Old 48- and 64-byte records are upcast: new fields default to NAN/0.
    void scan(const std::function<void(const QualRecord&)>& cb) const {
        if (!data_) return;
        const uint64_t stride = header_->record_stride;
        for (uint32_t i = 0; i < header_->n_records; ++i) {
            if (!old_layout_ && !medium_layout_ && !v3_layout_ && !v4_layout_) {
                cb(*reinterpret_cast<const QualRecord*>(base_ + i * stride));
            } else {
                QualRecord r = QualRecord::make_empty(0);
                uint64_t copy_sz = old_layout_ ? kOldStride
                                 : medium_layout_ ? kMediumStride
                                 : v3_layout_ ? kV3Stride
                                 : kV4Stride;
                __builtin_memcpy(&r, base_ + i * stride, copy_sz);
                cb(r);
            }
        }
    }

private:
    const uint8_t*    data_          = nullptr;
    const QualHeader* header_        = nullptr;
    const uint8_t*    base_          = nullptr;
    bool              old_layout_    = false;
    bool              medium_layout_ = false;
    bool              v3_layout_     = false;
    bool              v4_layout_     = false;
};

} // namespace genopack
