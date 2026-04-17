#pragma once
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/types.hpp>
#include <cstdint>
#include <memory>
#include <mutex>
#include <functional>
#include <optional>
#include <vector>

namespace genopack {

// ── SKCH V4 (dual-seed) ──────────────────────────────────────────────────────
//
// V4 stores TWO planar signature blocks per k (seed=42 and seed=43), computed
// in a single k-mer scan. Masks are seed-independent (they record which OPH
// bins were observed by at least one real k-mer) and are stored once per k.
//
// No backwards compatibility: V1/V2/V3 archives are rejected.
//
// Section payload layout (all fields before frame data are UNCOMPRESSED):
//   [SkchSeekHdr         (96 bytes)]
//   [SkchFrameDesc       (24 bytes × n_frames)]
//   [uint64_t genome_ids[n_genomes]]        ← for binary search without decompression
//   [uint64_t genome_lengths[n_genomes]]    ← parallel to genome_ids
//   [Frame 0..n_frames-1: independent zstd frames]
//     Each frame covers rows [frame_idx*frame_size, min(n_genomes, (frame_idx+1)*frame_size))
//     Frame payload (compressed) — planar by k:
//       [uint32_t n_real_bins[n_k × frame_n]]
//       [uint16_t sigs1 [n_k × frame_n × sketch_size]]  (seed = seed1, default 42)
//       [uint16_t sigs2 [n_k × frame_n × sketch_size]]  (seed = seed2, default 43)
//       [uint64_t masks [n_k × frame_n × mask_words]]
//
// Offsets within a decompressed frame of raw size R, given
//     S    = sketch_size
//     W    = mask_words = ceil(S/64)
//     NK   = n_kmer_sizes
//     FN   = frame_n
//     SIG  = 2 * NK * FN * S   (u16 elements)
//     MSK  = NK * FN * W       (u64 elements)
//   nrb  @ byte 0
//   sigs1@ byte sizeof(u32) * NK * FN
//   sigs2@ byte sizeof(u32) * NK * FN + sizeof(u16) * NK * FN * S
//   masks@ byte sizeof(u32) * NK * FN + 2 * sizeof(u16) * NK * FN * S

static constexpr char     SKCH_V4_MAGIC_STR[4] = {'S','K','C','4'};
static constexpr uint32_t SKCH_V4_MAGIC       = 0x34434B53u; // "SKC4"
static constexpr uint32_t SKCH_V4_FRAME_SIZE  = 16384;        // genomes per frame (default)

// Retained only for rejection messages — V3 archives are no longer supported.
static constexpr uint32_t SKCH_V3_MAGIC       = 0x33434B53u; // "SKC3"

struct SkchSeekHdr {
    uint32_t magic;           // SKCH_V4_MAGIC
    uint32_t n_frames;
    uint32_t frame_size;      // genomes per frame (last frame may be smaller)
    uint32_t n_genomes;
    uint32_t sketch_size;
    uint32_t n_kmer_sizes;    // ≤ 8
    uint32_t kmer_sizes[8];   // sorted ascending
    uint32_t syncmer_s;
    uint32_t mask_words;
    uint64_t seed1;           // primary seed (default 42)
    uint64_t seed2;           // dual seed (default 43) — NOT the old "densification seed"
    uint8_t  reserved[16];
};
static_assert(sizeof(SkchSeekHdr) == 96);

struct SkchFrameDesc {
    uint64_t data_offset;     // byte offset from start of section payload to compressed frame
    uint32_t compressed_size;
    uint32_t n_genomes;       // genomes in this frame (≤ frame_size; last frame may be less)
};
static_assert(sizeof(SkchFrameDesc) == 16);

// ── SketchResult ─────────────────────────────────────────────────────────────
// V4 always yields BOTH sig1 and sig2 pointers (populated from the same
// decompressed frame). Pointers are valid only during the sketch_for_ids
// callback invocation, or until ensure_loaded() is released.

struct SketchResult {
    const uint16_t* sig;            // sketch_size elements, seed = seed1
    const uint16_t* sig2;           // sketch_size elements, seed = seed2 (V4 only)
    const uint64_t* mask;           // mask_words elements
    uint32_t        n_real_bins;
    uint32_t        mask_words;
    uint64_t        genome_length;
    uint32_t        sketch_size;
    uint32_t        kmer_size;      // k used to build these OPH signatures
};

// ── SkchWriter ───────────────────────────────────────────────────────────────
// V4 single-k writer. Sigs1/sigs2/masks are spilled to a tmpfile immediately
// on add() so only small per-genome fixed data lives in RAM during the build.

class SkchWriter {
public:
    SkchWriter(uint32_t sketch_size, uint32_t kmer_size,
               uint32_t syncmer_s = 0, uint64_t seed1 = 42, uint64_t seed2 = 43,
               std::string spill_dir = {});
    ~SkchWriter();
    SkchWriter(const SkchWriter&) = delete;
    SkchWriter& operator=(const SkchWriter&) = delete;

    // Dual-sketch add: sig1 uses seed1, sig2 uses seed2. Masks are seed-independent.
    void add(GenomeId genome_id,
             const std::vector<uint16_t>& oph_sig1,
             const std::vector<uint16_t>& oph_sig2,
             uint32_t n_real_bins,
             uint64_t genome_length,
             const std::vector<uint64_t>& mask);

    SectionDesc finalize(AppendWriter& writer, uint64_t section_id);

private:
    uint32_t sketch_size_;
    uint32_t kmer_size_;
    uint32_t syncmer_s_;
    uint64_t seed1_;
    uint64_t seed2_;
    uint32_t mask_words_;
    size_t   record_size_;  // 2*sig_bytes + mask_bytes (per genome in spill)

    std::vector<uint64_t> ids_;
    std::vector<uint32_t> n_real_bins_;
    std::vector<uint64_t> genome_lengths_;
    FILE*                 spill_fp_ = nullptr;
};

// ── SkchWriterMultiK ─────────────────────────────────────────────────────────
// V4 multi-k writer. Spills sigs1_k*, sigs2_k*, masks_k* per genome; finalize
// writes frames in the V4 planar layout.

class SkchWriterMultiK {
public:
    SkchWriterMultiK(std::vector<uint32_t> kmer_sizes, uint32_t sketch_size,
                     uint32_t syncmer_s = 0, uint64_t seed1 = 42, uint64_t seed2 = 43,
                     std::string spill_dir = {});
    ~SkchWriterMultiK();
    SkchWriterMultiK(const SkchWriterMultiK&) = delete;
    SkchWriterMultiK& operator=(const SkchWriterMultiK&) = delete;

    // sigs1_per_k[i] / sigs2_per_k[i] / masks_per_k[i] correspond to kmer_sizes[i].
    void add(GenomeId genome_id,
             uint64_t genome_length,
             const std::vector<std::vector<uint16_t>>& sigs1_per_k,
             const std::vector<std::vector<uint16_t>>& sigs2_per_k,
             const std::vector<uint32_t>&              n_real_bins_per_k,
             const std::vector<std::vector<uint64_t>>& masks_per_k);

    SectionDesc finalize(AppendWriter& writer, uint64_t section_id);

    const std::vector<uint32_t>& kmer_sizes() const { return kmer_sizes_; }

private:
    std::vector<uint32_t>              kmer_sizes_;
    uint32_t                           sketch_size_;
    uint32_t                           syncmer_s_;
    uint64_t                           seed1_, seed2_;
    uint32_t                           mask_words_;
    size_t                             spill_record_size_;  // bytes per genome in spill

    std::vector<uint64_t>              ids_;
    std::vector<uint64_t>              genome_lengths_;
    std::vector<std::vector<uint32_t>> n_real_bins_;  // [ki][genome_idx]
    FILE*                              spill_fp_ = nullptr;
};

// ── SkchReader ───────────────────────────────────────────────────────────────

class SkchReader {
public:
    // Opens a SKC4 section. Rejects V1/V2/V3 archives with a clear error.
    void open(const uint8_t* base, uint64_t offset, uint64_t compressed_size);

    // Check if a genome_id exists in this section (uses lightweight id index, no decompress).
    bool contains(GenomeId genome_id) const;

    // V4 stores both seeds natively.
    bool has_sig2() const { return true; }

    std::optional<SketchResult> sketch_for(GenomeId genome_id) const;

    // Param-aware: picks the right k slot, then hierarchically slices to
    // requested_sz bins if requested_sz < stored size.
    std::optional<SketchResult> sketch_for(GenomeId genome_id,
                                           uint32_t k, uint32_t requested_sz) const;

    uint32_t version()        const { return 4; }
    uint32_t n_genomes()      const { return n_genomes_; }
    uint32_t sketch_size()    const { return sketch_size_; }
    uint32_t kmer_size()      const { return kmer_size_; }   // first k
    uint32_t n_kmer_sizes()   const { return n_kmer_sizes_; }
    uint32_t kmer_size_at(uint32_t i) const { return (i < n_kmer_sizes_) ? kmer_sizes_[i] : 0; }
    uint32_t syncmer_s()  const { return syncmer_s_; }
    uint64_t seed1()      const { return seed1_; }
    uint64_t seed2()      const { return seed2_; }
    const std::vector<uint64_t>& genome_ids() const { return id_index_; }

    bool has_kmer_size(uint32_t k) const;

    // Batch lookup: decompresses only the frames that contain the requested
    // genome_ids. sorted_ids must be sorted ascending.
    //
    // Threading: cb may be invoked CONCURRENTLY from multiple OMP worker
    // threads, each with a UNIQUE row index. SketchResult pointers
    // (sig / sig2 / mask) are valid only during the cb invocation —
    // consumers must copy out before returning.
    using SketchCallback = std::function<void(size_t, const SketchResult&)>;
    void sketch_for_ids(const std::vector<GenomeId>& sorted_ids,
                        uint32_t k, uint32_t sz,
                        const SketchCallback& cb) const;

    // Peek version and stored kmer_sizes without building id_index. Returns
    // {4, kmer_sizes} for V4; throws for anything else.
    static std::pair<uint32_t, std::vector<uint32_t>>
        peek_params(const uint8_t* base, uint64_t offset, uint64_t compressed_sz);

    bool matches_params(uint32_t sketch_size, uint32_t kmer_size,
                        uint64_t seed1 = 42, uint64_t seed2 = 43) const {
        return sketch_size_ == sketch_size && kmer_size_ == kmer_size
            && seed1_ == seed1 && seed2_ == seed2;
    }

private:
    // Shared slice helper.
    static std::optional<SketchResult> apply_slice(SketchResult base, uint32_t requested_sz,
                                                    uint32_t mask_words);

    const uint8_t* cdata_    = nullptr;
    uint64_t       cdata_sz_ = 0;

    std::vector<uint64_t> id_index_;     // sorted genome_ids
    std::vector<uint64_t> genome_lengths_;

    uint32_t                    frame_sz_ = 0;
    std::vector<SkchFrameDesc>  frames_;
    const uint8_t*              section_base_ = nullptr;

    uint32_t n_kmer_sizes_ = 0;
    uint32_t kmer_sizes_[8] = {};
    uint32_t n_genomes_   = 0;
    uint32_t sketch_size_ = 0;
    uint32_t kmer_size_   = 0;   // kmer_sizes_[0]
    uint32_t syncmer_s_   = 0;
    uint32_t mask_words_  = 0;
    uint64_t seed1_       = 0;
    uint64_t seed2_       = 0;
};

} // namespace genopack
