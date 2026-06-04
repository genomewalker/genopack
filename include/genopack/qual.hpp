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
    uint8_t  fmh_minority_u8;             //  1  fmh_minority_fraction × 255; 0=clean/not computed
    uint8_t  marker_completeness_u8;     //  1  marker_completeness encoded: 0=not_scored, 1-255=(v-1)/254.0
    uint16_t marker_redundancy_u16;      //  2  raw marker_redundancy × 65535 (0xFFFF = not scored)
    float    chargaff_parity;            //  4  Chargaff 2nd parity score [0,1]; NAN = not computed
    float    spectral_gap;               //  4  1 - |λ₂| from k3 Markov op; NAN = not computed
    float    scale_kink;                 //  4  log₂(W*) dyadic-window scale; NAN = not computed
    uint8_t  cross_genus_u8;             //  1  contamination_cross_genus × 255; 0=clean
    uint8_t  sketch_fill_u8;             //  1  completeness_sketch_fill × 200 (200=100%, >200 allowed up to 255)
    uint8_t  _pad[2];                    //  2  padding to reach 80 bytes
    // total = 80

    // qual_flags bits
    static constexpr uint8_t QUAL_FLAG_MIX_NO_DATA        = 0x01; // mixture model had < min_windows
    static constexpr uint8_t QUAL_FLAG_FRAG_GATED         = 0x02; // completeness_fragmentation < 0.65
    static constexpr uint8_t QUAL_FLAG_LEAKAGE_UNRELIABLE = 0x04; // leakage_residual above threshold
    static constexpr uint8_t QUAL_FLAG_COHERENCE_GATED    = 0x08; // self_coherence < 0.80
    static constexpr uint8_t QUAL_FLAG_MARKER_SCORED      = 0x10; // marker_redundancy_u16 is valid
    static constexpr uint8_t QUAL_FLAG_BUILD_SIGNALS      = 0x20; // AllSignals computed at build time (chargaff/spectral/scale_kink valid)
    static constexpr uint8_t QUAL_FLAG_GCOV_SCORED        = 0x40; // per-contig GCOV scoring cached (contig_outlier/spe/sibling/rho valid)
    static constexpr uint8_t QUAL_FLAG_FMH_AXIS_ABSENT    = 0x80; // FMH minority axis unavailable -> contamination tier used the sketch-leakage fallback

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
        r.marker_completeness_u8      = 0;
        return r;
    }
};
static_assert(sizeof(QualRecord) == 80);

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

private:
    std::vector<QualRecord> records_;
};

// ── Reader ────────────────────────────────────────────────────────────────────

class QualReader {
public:
    static constexpr uint64_t kOldStride    = 48; // pre-v2 layout without contamination_mixture
    static constexpr uint64_t kMediumStride = 64; // pre-v3 layout without chargaff/spectral/scale_kink

    void open(const uint8_t* data, uint64_t offset, uint64_t size) {
        if (size < sizeof(QualHeader))
            throw std::runtime_error("QUAL section too small");
        data_   = data + offset;
        header_ = reinterpret_cast<const QualHeader*>(data_);
        if (header_->magic != SEC_QUAL)
            throw std::runtime_error("QUAL: bad magic");

        const uint64_t stride = header_->record_stride;
        if (stride != sizeof(QualRecord) && stride != kOldStride && stride != kMediumStride)
            throw std::runtime_error("QUAL: unknown record_stride " + std::to_string(stride)
                                     + " — rebuild required");
        old_layout_    = (stride == kOldStride);
        medium_layout_ = (stride == kMediumStride);

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
            if (!old_layout_ && !medium_layout_) {
                cb(*reinterpret_cast<const QualRecord*>(base_ + i * stride));
            } else {
                QualRecord r = QualRecord::make_empty(0);
                __builtin_memcpy(&r, base_ + i * stride,
                                 old_layout_ ? kOldStride : kMediumStride);
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
};

} // namespace genopack
