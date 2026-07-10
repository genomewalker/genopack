#pragma once
#include "format.hpp"
#include "types.hpp"
#include "mmap_file.hpp"
#include <cstdint>
#include <cmath>
#include <vector>

namespace genopack {

// On-disk per-genome quality record (64 bytes, POD).
// Single fixed layout, size-sorted for tight packing.
struct QualRecord {
    uint64_t genome_id;                     //  8
    float    completeness_cluster_relative; //  4
    float    completeness_post_decontam;    //  4
    float    contamination_leakage;         //  4
    float    contamination_tnf_excess;      //  4
    float    completeness_aamer_core;        //  4  CORE aamer coverage [0,1]; NAN = not scored
    float    contamination_core_dup_mass;   //  4  non-saturating SCC dup mass Σ(c-1)/Σc; NAN = not scored
    float    accessory_ratio;               //  4  c0_query / median(c0_all); NAN = not scored
    uint16_t contamination_duplication_u16; //  2  redundancy_fraction: 0=not scored, else round(f*65534)+1
    uint16_t fmh_minority_u16;              //  2  fmh_minority_fraction × 65535 (rounded); 0=not scored/clean
    uint16_t contig_outlier_u16;            //  2  contamination_contig_outlier × 65535 (rounded); 0=not scored/clean
    uint8_t  support_tier;                  //  1
    uint8_t  spe_outlier_u8;                //  1  contamination_spe × 255; 0=clean
    uint8_t  rho_outlier_u8;                //  1  contamination_rho_outlier × 255; 0=clean
    uint8_t  qual_flags;                    //  1  bitfield — see QUAL_FLAG_* constants
    uint8_t  contig_outlier_u8;             //  1  contamination_contig_outlier × 255; 0=pure
    uint8_t  fmh_minority_u8;               //  1  LEGACY: fmh_minority_fraction × 255 (truncated); prefer fmh_minority_u16
    uint8_t  marker_completeness_u8;        //  1  marker_completeness encoded: 0=not_scored, 1-255=(v-1)/254.0
    uint8_t  cross_genus_u8;                //  1  contamination_cross_genus × 255; 0=clean
    uint8_t  quality_tier_u8;               //  1  LEGACY discrete tier; 0=not_set,1=LQ,2=MQ,3=HQ
    uint8_t  qscore_lo;                     //  1  low byte of continuous quality ∈ [0,1] × 65534 + 1; 0 = not scored
    uint8_t  qscore_hi;                     //  1  high byte (split so quality_tier_u8 keeps a fixed offset)
    // total = 56 (53 bytes used + 3 tail pad to 8-byte alignment)

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
    static constexpr uint8_t QUAL_FLAG_BUILD_SIGNALS      = 0x20; // AllSignals computed at build time (cache freshness)
    static constexpr uint8_t QUAL_FLAG_GCOV_SCORED        = 0x40; // per-contig GCOV scoring cached (contig_outlier/spe/rho valid)
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
        QualRecord r{};                          // zero-inits every field
        r.genome_id                     = gid;
        r.completeness_cluster_relative = NAN;
        r.completeness_post_decontam    = NAN;
        r.contamination_leakage         = NAN;
        r.contamination_tnf_excess      = NAN;
        r.completeness_aamer_core        = NAN;
        r.contamination_core_dup_mass    = NAN;  // default not-scored
        r.accessory_ratio                = NAN;
        r.quality_tier_u8                = QTIER_NOT_SET;
        return r;
    }
};
static_assert(sizeof(QualRecord) == 56);

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
    void open(const uint8_t* data, uint64_t offset, uint64_t size) {
        if (size < sizeof(QualHeader))
            throw std::runtime_error("QUAL section too small");
        data_   = data + offset;
        header_ = reinterpret_cast<const QualHeader*>(data_);
        if (header_->magic != SEC_QUAL)
            throw std::runtime_error("QUAL: bad magic");

        // Single fixed layout: the on-disk stride must equal the current record size.
        // No legacy-layout upcast — a mismatched pack is a hard error (rebuild).
        const uint64_t stride = header_->record_stride;
        if (stride != sizeof(QualRecord))
            throw std::runtime_error("QUAL: record_stride " + std::to_string(stride)
                                     + " != " + std::to_string(sizeof(QualRecord))
                                     + " — rebuild required (legacy QUAL layouts unsupported)");

        const uint64_t end = header_->records_offset
            + static_cast<uint64_t>(header_->n_records) * stride;
        if (end > size)
            throw std::runtime_error("QUAL section truncated");
        base_ = data_ + header_->records_offset;
    }

    bool is_open() const { return data_ != nullptr; }
    uint32_t n_records() const { return header_ ? header_->n_records : 0; }

    // Linear scan over the fixed 104-byte records.
    void scan(const std::function<void(const QualRecord&)>& cb) const {
        if (!data_) return;
        for (uint32_t i = 0; i < header_->n_records; ++i)
            cb(*reinterpret_cast<const QualRecord*>(base_ + i * sizeof(QualRecord)));
    }

private:
    const uint8_t*    data_   = nullptr;
    const QualHeader* header_ = nullptr;
    const uint8_t*    base_   = nullptr;
};

} // namespace genopack
