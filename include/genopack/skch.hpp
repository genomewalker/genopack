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

// ── SkchHeader: 64 bytes ─────────────────────────────────────────────────────
// On-disk layout (zstd-compressed as single frame):
//   [SkchHeader (64 bytes)]
//   [uint64_t genome_ids[n_genomes]]     ← sorted ascending
//   [SkchEntryFixed entries[n_genomes]]   ← parallel to genome_ids
//   [uint16_t sigs[n_genomes * sketch_size]]  ← contiguous OPH signatures
//   [uint64_t masks[n_genomes * mask_words]]  ← contiguous real-bins masks

struct SkchHeader {
    uint32_t magic;         // SEC_SKCH
    uint32_t version;       // 1
    uint32_t n_genomes;
    uint32_t sketch_size;   // e.g. 10000
    uint32_t kmer_size;     // e.g. 21
    uint32_t syncmer_s;     // 0 = disabled
    uint64_t seed1;         // 42
    uint64_t seed2;         // 1337
    uint8_t  reserved[24];
};
static_assert(sizeof(SkchHeader) == 64);

// ── MultiKSkchHeader: 96 bytes ────────────────────────────────────────────────
// v2 section: single zstd frame storing sketches at multiple k values.
// First 16 bytes (magic, version, n_genomes, sketch_size) are identical to v1.
// On-disk layout after the 96-byte header:
//   [uint64_t genome_ids[n_genomes]]                   ← sorted ascending (same as v1)
//   [uint64_t genome_lengths[n_genomes]]
//   [uint32_t n_real_bins[n_kmer_sizes × n_genomes]]   ← planar by k
//   [uint16_t sigs[n_kmer_sizes × n_genomes × sketch_size]]  ← planar by k
//   [uint64_t masks[n_kmer_sizes × n_genomes × mask_words]]  ← planar by k

struct MultiKSkchHeader {
    uint32_t magic;             // SEC_SKCH
    uint32_t version;           // 2
    uint32_t n_genomes;
    uint32_t sketch_size;
    uint32_t n_kmer_sizes;      // number of k values stored (max 8)
    uint32_t kmer_sizes[8];     // sorted ascending, e.g. [16, 21, 31, 0, ...]
    uint32_t syncmer_s;
    uint32_t mask_words;        // ceil(sketch_size / 64)
    uint32_t pad_;              // explicit alignment pad for seed1
    uint64_t seed1;
    uint64_t seed2;
    uint8_t  reserved[16];
};
static_assert(sizeof(MultiKSkchHeader) == 96);

struct SkchEntryFixed {
    uint32_t n_real_bins;
    uint32_t mask_words;    // ceil(sketch_size / 64)
    uint64_t genome_length;
};
static_assert(sizeof(SkchEntryFixed) == 16);

// ── SketchResult ─────────────────────────────────────────────────────────────
// Pointers into the decompressed buffer owned by SkchReader.

struct SketchResult {
    const uint16_t* sig;            // sketch_size elements
    const uint64_t* mask;           // mask_words elements
    uint32_t        n_real_bins;
    uint32_t        mask_words;
    uint64_t        genome_length;
    uint32_t        sketch_size;
    uint32_t        kmer_size;      // k used to build these OPH signatures
};

// ── V3 seekable format structs ────────────────────────────────────────────────
//
// Section payload layout (all fields before frame data are UNCOMPRESSED):
//   [SkchSeekHdr         (96 bytes)]
//   [SkchFrameDesc       (24 bytes × n_frames)]
//   [uint64_t genome_ids[n_genomes]]        ← for binary search without decompression
//   [uint64_t genome_lengths[n_genomes]]    ← parallel to genome_ids
//   [Frame 0..n_frames-1: independent zstd frames]
//     Each frame covers rows [frame_idx*frame_size, min(n_genomes, (frame_idx+1)*frame_size))
//     Frame payload (compressed):
//       [uint32_t n_real_bins[n_k × frame_n]]
//       [uint16_t sigs      [n_k × frame_n × sketch_size]]
//       [uint64_t masks     [n_k × frame_n × mask_words]]

static constexpr uint32_t SKCH_V3_MAGIC       = 0x33434B53u; // "SKC3"
static constexpr uint32_t SKCH_V3_FRAME_SIZE  = 16384;        // genomes per frame (default)

struct SkchSeekHdr {
    uint32_t magic;           // SKCH_V3_MAGIC
    uint32_t n_frames;
    uint32_t frame_size;      // genomes per frame (last frame may be smaller)
    uint32_t n_genomes;
    uint32_t sketch_size;
    uint32_t n_kmer_sizes;    // ≤ 8
    uint32_t kmer_sizes[8];   // sorted ascending
    uint32_t syncmer_s;
    uint32_t mask_words;
    uint64_t seed1;
    uint64_t seed2;
    uint8_t  reserved[16];
};
static_assert(sizeof(SkchSeekHdr) == 96);

struct SkchFrameDesc {
    uint64_t data_offset;     // byte offset from start of section payload to compressed frame
    uint32_t compressed_size;
    uint32_t n_genomes;       // genomes in this frame (≤ frame_size; last frame may be less)
};
static_assert(sizeof(SkchFrameDesc) == 16);

// ── SkchWriter ───────────────────────────────────────────────────────────────
// Sigs and masks are spilled to a tmpfile immediately on add() so that only
// small per-genome fixed data (ids, SkchEntryFixed) lives in RAM during the
// build. finalize() reads the spill back in sorted order to build the section.

class SkchWriter {
public:
    SkchWriter(uint32_t sketch_size, uint32_t kmer_size,
               uint32_t syncmer_s = 0, uint64_t seed1 = 42, uint64_t seed2 = 1337,
               std::string spill_dir = {});
    ~SkchWriter();
    SkchWriter(const SkchWriter&) = delete;
    SkchWriter& operator=(const SkchWriter&) = delete;

    void add(GenomeId genome_id,
             const std::vector<uint16_t>& oph_sig,
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
    size_t   record_size_;

    std::vector<uint64_t>       ids_;
    std::vector<SkchEntryFixed> fixed_;
    FILE*                       spill_fp_ = nullptr;
};

// ── SkchWriterMultiK ─────────────────────────────────────────────────────────
// Writes a v2 SKCH section storing OPH sketches at multiple k-mer sizes in one
// compressed frame with planar layout for cache-efficient per-k access.
// Per-genome spill: sigs_k0[S] | sigs_k1[S] | ... | masks_k0[W] | masks_k1[W] | ...
// Finalize writes planar: all_sigs_k0 | all_sigs_k1 | ... | all_masks_k0 | ...

class SkchWriterMultiK {
public:
    SkchWriterMultiK(std::vector<uint32_t> kmer_sizes, uint32_t sketch_size,
                     uint32_t syncmer_s = 0, uint64_t seed1 = 42, uint64_t seed2 = 1337,
                     std::string spill_dir = {});
    ~SkchWriterMultiK();
    SkchWriterMultiK(const SkchWriterMultiK&) = delete;
    SkchWriterMultiK& operator=(const SkchWriterMultiK&) = delete;

    // sigs_per_k[i] and masks_per_k[i] correspond to kmer_sizes[i].
    void add(GenomeId genome_id,
             uint64_t genome_length,
             const std::vector<std::vector<uint16_t>>& sigs_per_k,
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

class SkchReader;  // forward declaration for SkchV2Repacker

// ── SkchV2Repacker ────────────────────────────────────────────────────────────
// Two-phase V2→V3 repack that avoids loading 50+ GB into RAM.
// Phase 1: prepare() — stream-decompress V2 frame to a temp file (call while source mmap is open).
// Phase 2: write_v3() — write seekable V3 section to AppendWriter from the temp file mmap.

class SkchV2Repacker {
public:
    SkchV2Repacker() = default;
    ~SkchV2Repacker();
    SkchV2Repacker(const SkchV2Repacker&) = delete;
    SkchV2Repacker& operator=(const SkchV2Repacker&) = delete;

    void        prepare(const SkchReader& src, const std::string& tmp_dir);
    bool        ready()  const { return fd_ >= 0; }
    SectionDesc write_v3(AppendWriter& writer, uint64_t section_id);

private:
    int      fd_           = -1;
    size_t   sz_           = 0;
    uint32_t n_genomes_    = 0;
    uint32_t sketch_size_  = 0;
    uint32_t n_kmer_sizes_ = 0;
    uint32_t mask_words_   = 0;
    uint32_t kmer_sizes_[8] = {};
    uint32_t syncmer_s_    = 0;
    uint64_t seed1_        = 0;
    uint64_t seed2_        = 0;
};

// ── SkchReader ───────────────────────────────────────────────────────────────

class SkchReader {
public:
    void open(const uint8_t* base, uint64_t offset, uint64_t compressed_size);

    // Check if a genome_id exists in this section (uses lightweight id index, no full decompress).
    bool contains(GenomeId genome_id) const;

    // Ensure full decompressed data is loaded (idempotent).
    void ensure_loaded() const;

    // Release the decompressed buffer to free memory. The reader can re-decompress on next access.
    void release() const;

    bool   is_loaded()    const { return !buf_.empty(); }
    size_t memory_bytes() const { return buf_.size(); }

    std::optional<SketchResult> sketch_for(GenomeId genome_id) const;

    // Param-aware: picks the right k slot (v2) or checks kmer_size (v1), then
    // hierarchically slices to requested_sz bins if requested_sz < stored size.
    std::optional<SketchResult> sketch_for(GenomeId genome_id,
                                           uint32_t k, uint32_t requested_sz) const;

    uint32_t version()        const { return version_; }
    uint32_t n_genomes()      const { return n_genomes_; }
    uint32_t sketch_size()    const { return sketch_size_; }
    uint32_t kmer_size()      const { return kmer_size_; }   // first k (compat)
    uint32_t n_kmer_sizes()   const { return n_kmer_sizes_; }
    uint32_t kmer_size_at(uint32_t i) const { return (i < n_kmer_sizes_) ? kmer_sizes_[i] : 0; }
    uint32_t syncmer_s()  const { return syncmer_s_; }
    uint64_t seed1()      const { return seed1_; }
    uint64_t seed2()      const { return seed2_; }
    // Sorted genome_id index (available without full decompression).
    const std::vector<uint64_t>& genome_ids() const { return id_index_; }

    // True if this section stores sketches for k (v1: kmer_size_==k; v2: k in kmer_sizes_[]).
    bool has_kmer_size(uint32_t k) const;

    // Batch lookup: for V3 archives, decompresses only the frames that contain
    // the requested genome_ids. sorted_ids must be sorted ascending.
    // Falls back to per-genome sketch_for() for V1/V2.
    // cb(row_in_sorted_ids, SketchResult) — called for each hit.
    //
    // Threading: cb may be invoked CONCURRENTLY from multiple OMP worker
    // threads. Each invocation carries a UNIQUE row index (no duplicates),
    // so pre-sized output indexed by row is race-free without locks. Shared
    // counters/maps inside cb require their own synchronisation.
    // SketchResult.sig/.mask point into a thread-local frame buffer that is
    // only valid while cb runs; consumers must copy out before returning.
    using SketchCallback = std::function<void(size_t, const SketchResult&)>;
    void sketch_for_ids(const std::vector<GenomeId>& sorted_ids,
                        uint32_t k, uint32_t sz,
                        const SketchCallback& cb) const;

    // Peek version and stored kmer_sizes from compressed data without building id_index.
    // Decompresses only the section header (~96 bytes). Thread-safe, stateless.
    static std::pair<uint32_t, std::vector<uint32_t>>
        peek_params(const uint8_t* base, uint64_t offset, uint64_t compressed_sz);

    bool matches_params(uint32_t sketch_size, uint32_t kmer_size,
                        uint64_t seed1 = 42, uint64_t seed2 = 1337) const {
        return sketch_size_ == sketch_size && kmer_size_ == kmer_size
            && seed1_ == seed1 && seed2_ == seed2;
    }

private:
    friend class SkchV2Repacker;

    void decompress_full() const;
    void parse_buf() const;

    // Shared slice helper used by both v1 and v2 sketch_for(id, k, sz).
    static std::optional<SketchResult> apply_slice(SketchResult base, uint32_t requested_sz,
                                                    uint32_t mask_words);

    mutable std::unique_ptr<std::mutex> load_mu_ = std::make_unique<std::mutex>();

    const uint8_t* cdata_    = nullptr;
    uint64_t       cdata_sz_ = 0;

    std::vector<uint64_t> id_index_;  // sorted genome_ids for O(log n) lookup

    mutable std::vector<uint8_t>       buf_;
    // v1 decompressed pointers
    mutable const uint64_t*            ids_     = nullptr;
    mutable const SkchEntryFixed*      entries_ = nullptr;
    mutable const uint16_t*            sigs_    = nullptr;
    mutable const uint64_t*            masks_   = nullptr;
    // v2 decompressed pointers (planar layout)
    mutable const uint64_t*            genome_lengths_v2_ = nullptr;
    mutable const uint32_t*            n_real_bins_v2_    = nullptr;  // [ki*n + pos]
    mutable const uint16_t*            sigs_v2_           = nullptr;  // [(ki*n+pos)*S]
    mutable const uint64_t*            masks_v2_          = nullptr;  // [(ki*n+pos)*W]

    // V3 seekable format state (populated by open() when magic == SKCH_V3_MAGIC).
    bool                        v3_          = false;
    uint32_t                    v3_frame_sz_ = 0;
    std::vector<SkchFrameDesc>  v3_frames_;
    // V3 uncompressed genome metadata (no decompression needed).
    std::vector<uint64_t>       v3_genome_lengths_;
    // base pointer into the mmap for per-frame decompression.
    const uint8_t*              v3_section_base_ = nullptr;

    uint32_t version_      = 1;
    uint32_t n_kmer_sizes_ = 0;
    uint32_t kmer_sizes_[8] = {};
    uint32_t n_genomes_   = 0;
    uint32_t sketch_size_ = 0;
    uint32_t kmer_size_   = 0;   // kmer_sizes_[0] for compat
    uint32_t syncmer_s_   = 0;
    uint32_t mask_words_  = 0;
    uint64_t seed1_       = 0;
    uint64_t seed2_       = 0;
};

} // namespace genopack
