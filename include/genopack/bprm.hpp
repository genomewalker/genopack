#pragma once
#include "format.hpp"
#include "mmap_file.hpp"
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <xxhash.h>

namespace genopack {

// SEC_BPRM — self-describing build parameters (one per archive).
// Records the concrete parameter VALUES every param-bearing section was built
// with (OPH sketch, FracMinHash, TNF, model versions) plus a params_hash for
// O(1) cross-section consistency checks. Replaces hardcoded literals (fmhr.cpp
// k=21,c=125) and anchors GPK3 content-addressed derivation (build_params_hash).

struct BprmHeader {                  // 128 bytes — explicit offsets, no natural-packing reliance
    uint32_t magic;                  //   0  SEC_BPRM
    uint16_t version;                //   4  schema version (1)
    uint16_t header_size;            //   6  sizeof(BprmHeader); readers tolerate growth
    uint32_t sketch_kmer_size;       //   8  single-k OPH path
    uint32_t sketch_size;            //  12  OPH bins
    uint32_t sketch_syncmer_s;       //  16  0 = disabled
    uint32_t n_kmer_sizes;           //  20  0 => single-k via sketch_kmer_size
    uint32_t kmer_sizes[8];          //  24  multi-k set (32 B)
    uint64_t sketch_seed;            //  56
    uint32_t fmh_k;                  //  64  FracMinHash k
    uint32_t fmh_c;                  //  68  FracMinHash density (1/c)
    uint32_t tnf_order;              //  72  TNF k
    uint32_t gstx_model_version;     //  76
    uint32_t gcov_model_version;     //  80
    uint8_t  taxonomy_rank;          //  84  'g' / 'f'
    uint8_t  _pad0[3];               //  85
    uint64_t params_hash;            //  88  XXH3_64 of this struct, params_hash field zeroed
    uint32_t build_flags;            //  96  section-set toggles (BprmBuildFlags)
    uint32_t micro_genus_threshold;  // 100  min genus members for a dedicated shard + consensus model
    float    core_theta;             // 104  CORE prevalence threshold (0 = CORE disabled)
    uint8_t  _reserved[20];          // 108  pad to 128
};
static_assert(sizeof(BprmHeader) == 128);

// Section-set toggles folded into build_flags (and thus params_hash): two
// archives with the same params_hash produce the same section SET, so a verbatim
// --from-gpk repack can reuse one for the other without losing/keeping a section
// the new flags would have dropped/added.
enum BprmBuildFlags : uint32_t {
    BPRM_F_CIDX     = 1u << 0,
    BPRM_F_SKETCH   = 1u << 1,
    BPRM_F_GSTX     = 1u << 2,
    BPRM_F_GCOV     = 1u << 3,
    BPRM_F_TAXGROUP = 1u << 4,
    BPRM_F_MARKERS  = 1u << 5,
    BPRM_F_CORE     = 1u << 6,
};
static_assert(offsetof(BprmHeader, sketch_seed) == 56);
static_assert(offsetof(BprmHeader, core_theta) == 104);
static_assert(offsetof(BprmHeader, params_hash) == 88);
static_assert(offsetof(BprmHeader, micro_genus_threshold) == 100);

// XXH3_64 over the full header with params_hash zeroed — the canonical value.
inline uint64_t bprm_compute_params_hash(const BprmHeader& h) {
    BprmHeader tmp = h;
    tmp.params_hash = 0;
    return XXH3_64bits(&tmp, sizeof(tmp));
}

// ── Writer (header-only; the section is a single fixed-size record) ─────────────
class BprmWriter {
public:
    explicit BprmWriter(BprmHeader params) : params_(params) {}

    SectionDesc finalize(AppendWriter& w, uint64_t section_id) {
        params_.magic       = SEC_BPRM;
        params_.version     = 1;
        params_.header_size = sizeof(BprmHeader);
        params_.params_hash = bprm_compute_params_hash(params_);

        w.align(8);
        const uint64_t start = w.current_offset();
        w.append(&params_, sizeof(params_));
        const uint64_t end = w.current_offset();

        SectionDesc sd{};
        sd.type              = SEC_BPRM;
        sd.version           = 1;
        sd.flags             = 0;
        sd.section_id        = section_id;
        sd.file_offset       = start;
        sd.compressed_size   = end - start;
        sd.uncompressed_size = end - start;
        sd.item_count        = 1;
        return sd;
    }

    uint64_t params_hash() const noexcept { return params_.params_hash; }

private:
    BprmHeader params_;
};

// ── Reader ──────────────────────────────────────────────────────────────────────
class BprmReader {
public:
    void open(const uint8_t* data, uint64_t offset, uint64_t size) {
        if (size < sizeof(BprmHeader))
            throw std::runtime_error("BPRM section too small");
        std::memcpy(&hdr_, data + offset, sizeof(BprmHeader));
        if (hdr_.magic != SEC_BPRM)
            throw std::runtime_error("BPRM: bad magic");
        if (hdr_.params_hash != bprm_compute_params_hash(hdr_))
            throw std::runtime_error("BPRM: params_hash mismatch (corrupt or tampered)");
        valid_ = true;
    }

    bool valid() const noexcept { return valid_; }
    const BprmHeader& header() const noexcept { return hdr_; }
    uint64_t params_hash() const noexcept { return hdr_.params_hash; }

private:
    BprmHeader hdr_{};
    bool valid_ = false;
};

} // namespace genopack
