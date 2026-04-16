#include <genopack/skch.hpp>
#include <zstd.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <unistd.h>

namespace genopack {

// ── SkchWriter ───────────────────────────────────────────────────────────────

SkchWriter::SkchWriter(uint32_t sketch_size, uint32_t kmer_size,
                       uint32_t syncmer_s, uint64_t seed1, uint64_t seed2,
                       std::string spill_dir)
    : sketch_size_(sketch_size)
    , kmer_size_(kmer_size)
    , syncmer_s_(syncmer_s)
    , seed1_(seed1)
    , seed2_(seed2)
    , mask_words_((sketch_size + 63) / 64)
    , record_size_(sizeof(uint16_t) * sketch_size + sizeof(uint64_t) * ((sketch_size + 63) / 64))
{
    // Resolve spill directory: explicit arg > env var > none (falls back to std::tmpfile)
    if (spill_dir.empty()) {
        const char* env = std::getenv("GENOPACK_SPILL_DIR");
        if (env && *env) spill_dir = env;
    }

    if (!spill_dir.empty()) {
        char path[4096];
        std::snprintf(path, sizeof(path), "%s/genopack_skch_XXXXXX", spill_dir.c_str());
        int fd = ::mkstemp(path);
        if (fd < 0)
            throw std::runtime_error(std::string("SkchWriter: mkstemp failed: ") + path);
        ::unlink(path);
        spill_fp_ = ::fdopen(fd, "w+b");
        if (!spill_fp_) { ::close(fd); throw std::runtime_error("SkchWriter: fdopen failed"); }
    } else {
        spill_fp_ = std::tmpfile();
        if (!spill_fp_)
            throw std::runtime_error("SkchWriter: cannot create spill tmpfile");
    }
}

SkchWriter::~SkchWriter() {
    if (spill_fp_) { std::fclose(spill_fp_); spill_fp_ = nullptr; }
}

void SkchWriter::add(GenomeId genome_id,
                     const std::vector<uint16_t>& oph_sig,
                     uint32_t n_real_bins,
                     uint64_t genome_length,
                     const std::vector<uint64_t>& mask) {
    if (oph_sig.size() != sketch_size_)
        throw std::runtime_error("SkchWriter::add: sig size mismatch");
    if (mask.size() != mask_words_)
        throw std::runtime_error("SkchWriter::add: mask size mismatch");

    ids_.push_back(genome_id);
    SkchEntryFixed ef{};
    ef.n_real_bins   = n_real_bins;
    ef.mask_words    = mask_words_;
    ef.genome_length = genome_length;
    fixed_.push_back(ef);

    // Spill sig + mask to tmpfile — no RAM retained for these
    std::fwrite(oph_sig.data(), sizeof(uint16_t), sketch_size_, spill_fp_);
    std::fwrite(mask.data(),    sizeof(uint64_t), mask_words_,  spill_fp_);
}

SectionDesc SkchWriter::finalize(AppendWriter& writer, uint64_t section_id) {
    uint32_t n = static_cast<uint32_t>(ids_.size());

    // Sort order by genome_id
    std::vector<uint32_t> order(n);
    std::iota(order.begin(), order.end(), 0u);
    std::sort(order.begin(), order.end(),
              [&](uint32_t a, uint32_t b) { return ids_[a] < ids_[b]; });

    size_t ids_bytes     = sizeof(uint64_t)      * n;
    size_t entries_bytes = sizeof(SkchEntryFixed) * n;
    size_t sigs_bytes    = sizeof(uint16_t) * static_cast<size_t>(n) * sketch_size_;
    size_t masks_bytes   = sizeof(uint64_t) * static_cast<size_t>(n) * mask_words_;
    size_t total         = sizeof(SkchHeader) + ids_bytes + entries_bytes + sigs_bytes + masks_bytes;

    // Stream-compress using ZSTD streaming API to avoid a ~27 GB flat buffer.
    ZSTD_CStream* cstream = ZSTD_createCStream();
    if (!cstream) throw std::runtime_error("SkchWriter: ZSTD_createCStream failed");
    ZSTD_initCStream(cstream, 3);
    ZSTD_CCtx_setPledgedSrcSize(cstream, total); // embeds content size in frame header for fast decompression

    // Write compressed output via a small staging buffer into AppendWriter.
    const size_t OUT_BUF = 4 << 20; // 4 MB
    std::vector<uint8_t> out_buf(OUT_BUF);
    uint64_t section_start = writer.current_offset();

    auto flush_stream = [&](bool final_flush) {
        ZSTD_outBuffer zout{out_buf.data(), out_buf.size(), 0};
        size_t ret = final_flush
            ? ZSTD_endStream(cstream, &zout)
            : ZSTD_flushStream(cstream, &zout);
        if (ZSTD_isError(ret))
            throw std::runtime_error(std::string("SkchWriter: zstd stream: ") + ZSTD_getErrorName(ret));
        if (zout.pos > 0)
            writer.append(out_buf.data(), zout.pos);
        return ret;
    };

    auto compress_chunk = [&](const void* data, size_t size) {
        ZSTD_inBuffer zin{data, size, 0};
        while (zin.pos < zin.size) {
            ZSTD_outBuffer zout{out_buf.data(), out_buf.size(), 0};
            size_t ret = ZSTD_compressStream(cstream, &zout, &zin);
            if (ZSTD_isError(ret))
                throw std::runtime_error(std::string("SkchWriter: zstd compress: ") + ZSTD_getErrorName(ret));
            if (zout.pos > 0)
                writer.append(out_buf.data(), zout.pos);
        }
    };

    // Header
    SkchHeader hdr{};
    hdr.magic       = SEC_SKCH;
    hdr.version     = 1;
    hdr.n_genomes   = n;
    hdr.sketch_size = sketch_size_;
    hdr.kmer_size   = kmer_size_;
    hdr.syncmer_s   = syncmer_s_;
    hdr.seed1       = seed1_;
    hdr.seed2       = seed2_;
    compress_chunk(&hdr, sizeof(hdr));

    // Sorted genome IDs
    for (uint32_t i : order) compress_chunk(&ids_[i], sizeof(uint64_t));

    // Sorted SkchEntryFixed
    for (uint32_t i : order) compress_chunk(&fixed_[i], sizeof(SkchEntryFixed));

    // Sigs (sorted): read each record from spill, emit sig bytes
    std::vector<uint8_t> rec(record_size_);
    size_t sig_bytes  = sizeof(uint16_t) * sketch_size_;
    size_t mask_bytes = sizeof(uint64_t) * mask_words_;

    for (uint32_t rank = 0; rank < n; ++rank) {
        uint32_t i = order[rank];
        long offset = static_cast<long>(static_cast<size_t>(i) * record_size_);
        if (std::fseek(spill_fp_, offset, SEEK_SET) != 0)
            throw std::runtime_error("SkchWriter::finalize: fseek failed");
        if (std::fread(rec.data(), 1, record_size_, spill_fp_) != record_size_)
            throw std::runtime_error("SkchWriter::finalize: spill read failed");
        compress_chunk(rec.data(), sig_bytes);
    }

    // Masks (sorted): second pass over spill
    for (uint32_t rank = 0; rank < n; ++rank) {
        uint32_t i = order[rank];
        long offset = static_cast<long>(static_cast<size_t>(i) * record_size_);
        if (std::fseek(spill_fp_, offset, SEEK_SET) != 0)
            throw std::runtime_error("SkchWriter::finalize: fseek failed (masks)");
        if (std::fread(rec.data(), 1, record_size_, spill_fp_) != record_size_)
            throw std::runtime_error("SkchWriter::finalize: spill read failed (masks)");
        compress_chunk(rec.data() + sig_bytes, mask_bytes);
    }

    // Flush and end stream
    size_t remaining = 1;
    while (remaining != 0) remaining = flush_stream(true);
    ZSTD_freeCStream(cstream);

    std::fclose(spill_fp_);
    spill_fp_ = nullptr;
    uint64_t section_end = writer.current_offset();

    SectionDesc desc{};
    desc.type              = SEC_SKCH;
    desc.version           = 1;
    desc.flags             = 0;
    desc.section_id        = section_id;
    desc.file_offset       = section_start;
    desc.compressed_size   = section_end - section_start;
    desc.uncompressed_size = total; // logical uncompressed size
    desc.item_count        = n;
    desc.aux0              = sketch_size_;
    desc.aux1              = kmer_size_;
    std::memset(desc.checksum, 0, sizeof(desc.checksum));
    return desc;
}

// ── SkchReader ───────────────────────────────────────────────────────────────

void SkchReader::open(const uint8_t* base, uint64_t offset, uint64_t compressed_size) {
    cdata_    = base + offset;
    cdata_sz_ = compressed_size;

    // V3: uncompressed header — detect by magic before attempting zstd.
    if (compressed_size >= sizeof(SkchSeekHdr)) {
        uint32_t maybe_magic;
        std::memcpy(&maybe_magic, cdata_, sizeof(uint32_t));
        if (maybe_magic == SKCH_V3_MAGIC) {
            const auto* h = reinterpret_cast<const SkchSeekHdr*>(cdata_);
            v3_           = true;
            version_      = 3;
            n_genomes_    = h->n_genomes;
            sketch_size_  = h->sketch_size;
            mask_words_   = h->mask_words;
            n_kmer_sizes_ = std::min(h->n_kmer_sizes, 8u);
            for (uint32_t i = 0; i < 8; ++i) kmer_sizes_[i] = h->kmer_sizes[i];
            kmer_size_    = (n_kmer_sizes_ > 0) ? kmer_sizes_[0] : 0;
            syncmer_s_    = h->syncmer_s;
            seed1_        = h->seed1;
            seed2_        = h->seed2;
            v3_frame_sz_  = h->frame_size;
            v3_section_base_ = cdata_;

            // Frame table immediately after header.
            const auto* ft = reinterpret_cast<const SkchFrameDesc*>(cdata_ + sizeof(SkchSeekHdr));
            v3_frames_.assign(ft, ft + h->n_frames);

            // Uncompressed genome_ids and genome_lengths follow the frame table.
            const uint8_t* ids_ptr = cdata_ + sizeof(SkchSeekHdr)
                                     + sizeof(SkchFrameDesc) * h->n_frames;
            id_index_.resize(n_genomes_);
            std::memcpy(id_index_.data(), ids_ptr, sizeof(uint64_t) * n_genomes_);

            const uint8_t* len_ptr = ids_ptr + sizeof(uint64_t) * n_genomes_;
            v3_genome_lengths_.resize(n_genomes_);
            std::memcpy(v3_genome_lengths_.data(), len_ptr, sizeof(uint64_t) * n_genomes_);
            return;
        }
    }

    ZSTD_DStream* ds = ZSTD_createDStream();
    ZSTD_initDStream(ds);
    ZSTD_inBuffer zin{cdata_, cdata_sz_, 0};

    // Helper: read exactly n bytes from the decompression stream.
    auto read_exact = [&](void* dst, size_t n) -> bool {
        ZSTD_outBuffer zo{dst, n, 0};
        while (zo.pos < n) {
            size_t r = ZSTD_decompressStream(ds, &zo, &zin);
            if (ZSTD_isError(r) || r == 0) return false;
        }
        return true;
    };

    // Step 1: Read the 16-byte common prefix (magic, version, n_genomes, sketch_size).
    // Both v1 (64 bytes) and v2 (96 bytes) share these four fields at offset 0.
    struct CommonPrefix { uint32_t magic, version, n_genomes, sketch_size; };
    static_assert(sizeof(CommonPrefix) == 16);
    CommonPrefix pfx{};
    if (!read_exact(&pfx, sizeof(pfx))) {
        ZSTD_freeDStream(ds);
        throw std::runtime_error("SkchReader: truncated header prefix");
    }
    if (pfx.magic != SEC_SKCH) {
        ZSTD_freeDStream(ds);
        throw std::runtime_error("SkchReader: bad magic");
    }

    version_     = pfx.version;
    n_genomes_   = pfx.n_genomes;
    sketch_size_ = pfx.sketch_size;
    mask_words_  = (sketch_size_ + 63) / 64;

    if (version_ == 1) {
        // v1 tail: kmer_size(4) syncmer_s(4) seed1(8) seed2(8) reserved(24) = 48 bytes
        struct V1Tail { uint32_t kmer_size, syncmer_s; uint64_t seed1, seed2; uint8_t reserved[24]; };
        static_assert(sizeof(V1Tail) == 48);
        V1Tail t{};
        if (!read_exact(&t, sizeof(t))) {
            ZSTD_freeDStream(ds);
            throw std::runtime_error("SkchReader: truncated v1 header");
        }
        kmer_size_    = t.kmer_size;
        syncmer_s_    = t.syncmer_s;
        seed1_        = t.seed1;
        seed2_        = t.seed2;
        n_kmer_sizes_ = 1;
        kmer_sizes_[0] = kmer_size_;
    } else if (version_ == 2) {
        // v2 tail: n_kmer_sizes(4) kmer_sizes[8](32) syncmer_s(4) mask_words(4)
        //          pad_(4) seed1(8) seed2(8) reserved(16) = 80 bytes
        struct V2Tail {
            uint32_t n_kmer_sizes, kmer_sizes[8], syncmer_s, mask_words_stored;
            uint32_t pad_;
            uint64_t seed1, seed2;
            uint8_t  reserved[16];
        };
        static_assert(sizeof(V2Tail) == 80);
        V2Tail t{};
        if (!read_exact(&t, sizeof(t))) {
            ZSTD_freeDStream(ds);
            throw std::runtime_error("SkchReader: truncated v2 header");
        }
        n_kmer_sizes_ = std::min(t.n_kmer_sizes, 8u);
        for (uint32_t i = 0; i < 8; ++i) kmer_sizes_[i] = t.kmer_sizes[i];
        syncmer_s_    = t.syncmer_s;
        seed1_        = t.seed1;
        seed2_        = t.seed2;
        kmer_size_    = (n_kmer_sizes_ > 0) ? kmer_sizes_[0] : 0;
    } else {
        ZSTD_freeDStream(ds);
        throw std::runtime_error("SkchReader: unknown SKCH version " + std::to_string(version_));
    }

    // Read genome_ids immediately after the header (same relative position in both versions).
    id_index_.resize(n_genomes_);
    if (!read_exact(id_index_.data(), sizeof(uint64_t) * n_genomes_)) {
        ZSTD_freeDStream(ds);
        throw std::runtime_error("SkchReader: truncated genome id index");
    }
    ZSTD_freeDStream(ds);
}

std::pair<uint32_t, std::vector<uint32_t>>
SkchReader::peek_params(const uint8_t* base, uint64_t offset, uint64_t compressed_sz) {
    const uint8_t* src = base + offset;
    ZSTD_DStream* ds = ZSTD_createDStream();
    ZSTD_initDStream(ds);
    ZSTD_inBuffer zin{src, compressed_sz, 0};

    auto read_exact = [&](void* dst, size_t n) -> bool {
        ZSTD_outBuffer zo{dst, n, 0};
        while (zo.pos < n) {
            size_t r = ZSTD_decompressStream(ds, &zo, &zin);
            if (ZSTD_isError(r) || r == 0) return false;
        }
        return true;
    };

    // Check for V3 (uncompressed header) before attempting zstd.
    if (compressed_sz >= sizeof(SkchSeekHdr)) {
        uint32_t maybe_magic;
        std::memcpy(&maybe_magic, src, sizeof(uint32_t));
        if (maybe_magic == SKCH_V3_MAGIC) {
            ZSTD_freeDStream(ds);
            const auto* h = reinterpret_cast<const SkchSeekHdr*>(src);
            uint32_t n = std::min(h->n_kmer_sizes, 8u);
            std::vector<uint32_t> ks(h->kmer_sizes, h->kmer_sizes + n);
            return {3, ks};
        }
    }

    struct CommonPrefix { uint32_t magic, version, n_genomes, sketch_size; };
    CommonPrefix pfx{};
    if (!read_exact(&pfx, sizeof(pfx)) || pfx.magic != SEC_SKCH) {
        ZSTD_freeDStream(ds);
        return {0, {}};
    }

    std::vector<uint32_t> ks;
    if (pfx.version == 1) {
        struct V1Tail { uint32_t kmer_size, syncmer_s; uint64_t seed1, seed2; uint8_t reserved[24]; };
        V1Tail t{};
        if (read_exact(&t, sizeof(t))) ks.push_back(t.kmer_size);
    } else if (pfx.version == 2) {
        struct V2Tail {
            uint32_t n_kmer_sizes, kmer_sizes[8], syncmer_s, mask_words;
            uint32_t pad_;
            uint64_t seed1, seed2; uint8_t reserved[16];
        };
        V2Tail t{};
        if (read_exact(&t, sizeof(t))) {
            uint32_t n = std::min(t.n_kmer_sizes, 8u);
            ks.assign(t.kmer_sizes, t.kmer_sizes + n);
        }
    }
    ZSTD_freeDStream(ds);
    return {pfx.version, ks};
}

bool SkchReader::has_kmer_size(uint32_t k) const {
    for (uint32_t i = 0; i < n_kmer_sizes_; ++i)
        if (kmer_sizes_[i] == k) return true;
    return false;
}

bool SkchReader::contains(GenomeId genome_id) const {
    if (id_index_.empty()) return false;
    uint32_t lo = 0, hi = n_genomes_;
    while (lo < hi) {
        uint32_t mid = lo + (hi - lo) / 2;
        if (id_index_[mid] == genome_id) return true;
        if (id_index_[mid] < genome_id) lo = mid + 1;
        else                             hi = mid;
    }
    return false;
}

void SkchReader::decompress_full() const {
    if (v3_) {
        // Assemble all frames into a v2-compatible flat buffer so that parse_buf()
        // and single-genome sketch_for() work correctly without code duplication.
        const uint32_t nk = n_kmer_sizes_;
        const size_t nrb_sz  = sizeof(uint32_t) * nk * n_genomes_;
        const size_t sigs_sz = sizeof(uint16_t) * nk * static_cast<size_t>(n_genomes_) * sketch_size_;
        const size_t msks_sz = sizeof(uint64_t) * nk * static_cast<size_t>(n_genomes_) * mask_words_;
        const size_t total   = sizeof(MultiKSkchHeader)
                             + sizeof(uint64_t) * n_genomes_ * 2  // ids + genome_lengths
                             + nrb_sz + sigs_sz + msks_sz;
        buf_.resize(total, 0);

        // Synthetic MultiKSkchHeader (parse_buf skips it via sizeof offset).
        auto* fhdr = reinterpret_cast<MultiKSkchHeader*>(buf_.data());
        fhdr->magic        = SEC_SKCH;
        fhdr->version      = 2;
        fhdr->n_genomes    = n_genomes_;
        fhdr->sketch_size  = sketch_size_;
        fhdr->n_kmer_sizes = nk;
        for (uint32_t i = 0; i < 8; ++i) fhdr->kmer_sizes[i] = kmer_sizes_[i];
        fhdr->syncmer_s    = syncmer_s_;
        fhdr->mask_words   = mask_words_;
        fhdr->seed1        = seed1_;
        fhdr->seed2        = seed2_;

        uint8_t* p = buf_.data() + sizeof(MultiKSkchHeader);
        std::memcpy(p, id_index_.data(), sizeof(uint64_t) * n_genomes_);
        p += sizeof(uint64_t) * n_genomes_;
        std::memcpy(p, v3_genome_lengths_.data(), sizeof(uint64_t) * n_genomes_);
        p += sizeof(uint64_t) * n_genomes_;

        uint32_t* dst_nrb = reinterpret_cast<uint32_t*>(p);  p += nrb_sz;
        uint16_t* dst_sig = reinterpret_cast<uint16_t*>(p);  p += sigs_sz;
        uint64_t* dst_msk = reinterpret_cast<uint64_t*>(p);

        uint32_t grow = 0;
        for (uint32_t fi = 0; fi < static_cast<uint32_t>(v3_frames_.size()); ++fi) {
            const SkchFrameDesc& fd = v3_frames_[fi];
            const uint32_t fn       = fd.n_genomes;
            const uint8_t* csrc     = v3_section_base_ + fd.data_offset;
            unsigned long long rsz  = ZSTD_getFrameContentSize(csrc, fd.compressed_size);
            if (rsz == ZSTD_CONTENTSIZE_ERROR || rsz == ZSTD_CONTENTSIZE_UNKNOWN)
                throw std::runtime_error("SkchReader V3 decompress_full: unknown frame size");
            std::vector<uint8_t> fb(static_cast<size_t>(rsz));
            if (ZSTD_isError(ZSTD_decompress(fb.data(), fb.size(), csrc, fd.compressed_size)))
                throw std::runtime_error("SkchReader V3 decompress_full: frame zstd error");

            const uint32_t* fnrb = reinterpret_cast<const uint32_t*>(fb.data());
            const uint16_t* fsig = reinterpret_cast<const uint16_t*>(
                fb.data() + sizeof(uint32_t) * nk * fn);
            const uint64_t* fmsk = reinterpret_cast<const uint64_t*>(
                fb.data() + sizeof(uint32_t) * nk * fn
                          + sizeof(uint16_t) * nk * fn * sketch_size_);

            for (uint32_t ki = 0; ki < nk; ++ki) {
                std::memcpy(dst_nrb + ki * n_genomes_ + grow,
                            fnrb   + ki * fn, sizeof(uint32_t) * fn);
                std::memcpy(dst_sig + (static_cast<size_t>(ki) * n_genomes_ + grow) * sketch_size_,
                            fsig   + ki * static_cast<size_t>(fn) * sketch_size_,
                            sizeof(uint16_t) * fn * sketch_size_);
                std::memcpy(dst_msk + (static_cast<size_t>(ki) * n_genomes_ + grow) * mask_words_,
                            fmsk   + ki * static_cast<size_t>(fn) * mask_words_,
                            sizeof(uint64_t) * fn * mask_words_);
            }
            grow += fn;
        }
        parse_buf();
        return;
    }

    unsigned long long raw_size = ZSTD_getFrameContentSize(cdata_, cdata_sz_);
    buf_.resize(static_cast<size_t>(raw_size));
    size_t dsize = ZSTD_decompress(buf_.data(), buf_.size(), cdata_, cdata_sz_);
    if (ZSTD_isError(dsize))
        throw std::runtime_error(std::string("SkchReader: zstd: ") + ZSTD_getErrorName(dsize));
    parse_buf();
}

void SkchReader::parse_buf() const {
    if (version_ == 1) {
        const uint8_t* p = buf_.data() + sizeof(SkchHeader);  // 64 bytes
        ids_     = reinterpret_cast<const uint64_t*>(p);
        p += sizeof(uint64_t) * n_genomes_;
        entries_ = reinterpret_cast<const SkchEntryFixed*>(p);
        p += sizeof(SkchEntryFixed) * n_genomes_;
        sigs_    = reinterpret_cast<const uint16_t*>(p);
        p += sizeof(uint16_t) * static_cast<size_t>(n_genomes_) * sketch_size_;
        masks_   = reinterpret_cast<const uint64_t*>(p);
    } else {  // v2 (or v3 assembled into v2 layout): MultiKSkchHeader = 96 bytes
        const uint8_t* p = buf_.data() + sizeof(MultiKSkchHeader);
        p += sizeof(uint64_t) * n_genomes_;               // genome_ids (same as id_index_)
        genome_lengths_v2_ = reinterpret_cast<const uint64_t*>(p);
        p += sizeof(uint64_t) * n_genomes_;
        n_real_bins_v2_ = reinterpret_cast<const uint32_t*>(p);
        p += sizeof(uint32_t) * n_kmer_sizes_ * n_genomes_;
        sigs_v2_  = reinterpret_cast<const uint16_t*>(p);
        p += sizeof(uint16_t) * n_kmer_sizes_ * static_cast<size_t>(n_genomes_) * sketch_size_;
        masks_v2_ = reinterpret_cast<const uint64_t*>(p);
    }
}

void SkchReader::ensure_loaded() const {
    std::lock_guard lock(*load_mu_);
    if (!buf_.empty()) return;
    decompress_full();
}

void SkchReader::release() const {
    std::lock_guard lock(*load_mu_);
    buf_.clear();
    buf_.shrink_to_fit();
    ids_     = nullptr; entries_ = nullptr; sigs_ = nullptr; masks_ = nullptr;
    genome_lengths_v2_ = nullptr; n_real_bins_v2_ = nullptr;
    sigs_v2_ = nullptr; masks_v2_ = nullptr;
}

// Binary search in id_index_; returns UINT32_MAX if not found.
static uint32_t find_genome_pos(const std::vector<uint64_t>& idx, uint64_t id) {
    uint32_t lo = 0, hi = static_cast<uint32_t>(idx.size());
    while (lo < hi) {
        uint32_t mid = lo + (hi - lo) / 2;
        if (idx[mid] == id) return mid;
        if (idx[mid] < id) lo = mid + 1;
        else               hi = mid;
    }
    return UINT32_MAX;
}

// Hierarchical prefix-slice: first requested_sz bins of an OPH sketch are a
// valid sub-sketch because OPH assigns bins by permutation index.
std::optional<SketchResult> SkchReader::apply_slice(SketchResult r,
                                                     uint32_t     requested_sz,
                                                     uint32_t     stored_mask_words) {
    if (requested_sz == r.sketch_size) return r;
    uint32_t new_mw = (requested_sz + 63) / 64;
    uint32_t cnt = 0;
    const uint64_t* m = r.mask;
    for (uint32_t w = 0; w + 1 < new_mw; ++w)
        cnt += static_cast<uint32_t>(__builtin_popcountll(m[w]));
    uint32_t rem = requested_sz % 64;
    uint64_t tail = (rem == 0) ? ~uint64_t{0} : ((uint64_t{1} << rem) - 1);
    cnt += static_cast<uint32_t>(__builtin_popcountll(m[new_mw - 1] & tail));
    r.sketch_size = requested_sz;
    r.mask_words  = new_mw;
    r.n_real_bins = cnt;
    (void)stored_mask_words;
    return r;
}

std::optional<SketchResult> SkchReader::sketch_for(GenomeId genome_id) const {
    if (n_genomes_ == 0) return std::nullopt;
    // For v2/v3: delegate to first k
    if ((version_ == 2 || v3_) && n_kmer_sizes_ > 0)
        return sketch_for(genome_id, kmer_sizes_[0], sketch_size_);

    uint32_t pos = find_genome_pos(id_index_, genome_id);
    if (pos == UINT32_MAX) return std::nullopt;
    ensure_loaded();

    SketchResult r{};
    r.sig           = sigs_ + static_cast<size_t>(pos) * sketch_size_;
    r.mask          = masks_ + static_cast<size_t>(pos) * mask_words_;
    r.n_real_bins   = entries_[pos].n_real_bins;
    r.mask_words    = entries_[pos].mask_words;
    r.genome_length = entries_[pos].genome_length;
    r.sketch_size   = sketch_size_;
    r.kmer_size     = kmer_size_;
    return r;
}

std::optional<SketchResult> SkchReader::sketch_for(GenomeId genome_id,
                                                    uint32_t k,
                                                    uint32_t requested_sz) const {
    if (requested_sz > sketch_size_) return std::nullopt;

    if (version_ == 2 || v3_) {
        // Find k index in kmer_sizes_[]
        uint32_t ki = UINT32_MAX;
        for (uint32_t i = 0; i < n_kmer_sizes_; ++i)
            if (kmer_sizes_[i] == k) { ki = i; break; }
        if (ki == UINT32_MAX) return std::nullopt;

        uint32_t pos = find_genome_pos(id_index_, genome_id);
        if (pos == UINT32_MAX) return std::nullopt;
        ensure_loaded();

        size_t gi = static_cast<size_t>(ki) * n_genomes_ + pos;
        SketchResult r{};
        r.sig           = sigs_v2_  + gi * sketch_size_;
        r.mask          = masks_v2_ + gi * mask_words_;
        r.n_real_bins   = n_real_bins_v2_[gi];
        r.mask_words    = mask_words_;
        r.genome_length = genome_lengths_v2_[pos];
        r.sketch_size   = sketch_size_;
        r.kmer_size     = k;
        return apply_slice(r, requested_sz, mask_words_);
    }

    // v1
    if (kmer_size_ != k) return std::nullopt;
    auto base = sketch_for(genome_id);
    if (!base) return std::nullopt;
    return apply_slice(*base, requested_sz, mask_words_);
}

void SkchReader::sketch_for_ids(const std::vector<GenomeId>& sorted_ids,
                                 uint32_t k, uint32_t sz,
                                 const SketchCallback& cb) const
{
    if (sorted_ids.empty()) return;

    // V1/V2: no seekable frames — fall back to per-genome lookup.
    if (!v3_) {
        for (size_t i = 0; i < sorted_ids.size(); ++i) {
            auto sk = (k > 0 && sz > 0)
                      ? sketch_for(sorted_ids[i], k, sz)
                      : sketch_for(sorted_ids[i]);
            if (sk) cb(i, *sk);
        }
        return;
    }

    // V3: group requests by frame, decompress each needed frame exactly once.
    uint32_t ki = UINT32_MAX;
    if (k > 0) {
        for (uint32_t i = 0; i < n_kmer_sizes_; ++i)
            if (kmer_sizes_[i] == k) { ki = i; break; }
        if (ki == UINT32_MAX) return;
    } else {
        ki = 0;
        k  = (n_kmer_sizes_ > 0) ? kmer_sizes_[0] : 0;
    }
    if (sz == 0 || sz > sketch_size_) sz = sketch_size_;

    // Map each sorted_id → (frame_idx, local_row, original index in sorted_ids).
    // sorted_ids is ascending, so id_index_ binary searches produce non-decreasing pos,
    // meaning frame_idx is also non-decreasing — the group scan below is a single pass.
    struct Hit { uint32_t frame_idx; uint32_t local_row; size_t orig_idx; };
    std::vector<Hit> hits;
    hits.reserve(sorted_ids.size());
    for (size_t i = 0; i < sorted_ids.size(); ++i) {
        uint32_t pos = find_genome_pos(id_index_, sorted_ids[i]);
        if (pos == UINT32_MAX) continue;
        hits.push_back({pos / v3_frame_sz_, pos % v3_frame_sz_, i});
    }

    size_t hi_start = 0;
    while (hi_start < hits.size()) {
        const uint32_t fi     = hits[hi_start].frame_idx;
        size_t         hi_end = hi_start;
        while (hi_end < hits.size() && hits[hi_end].frame_idx == fi) ++hi_end;

        const SkchFrameDesc& fd = v3_frames_[fi];
        const uint32_t frame_n  = fd.n_genomes;
        const uint8_t* csrc     = v3_section_base_ + fd.data_offset;

        unsigned long long raw_sz = ZSTD_getFrameContentSize(csrc, fd.compressed_size);
        if (raw_sz == ZSTD_CONTENTSIZE_ERROR || raw_sz == ZSTD_CONTENTSIZE_UNKNOWN) {
            hi_start = hi_end; continue;
        }
        std::vector<uint8_t> fbuf(static_cast<size_t>(raw_sz));
        if (ZSTD_isError(ZSTD_decompress(fbuf.data(), fbuf.size(), csrc, fd.compressed_size))) {
            hi_start = hi_end; continue;
        }

        // Frame layout (planar by k):
        //   n_real_bins: uint32_t [n_kmer_sizes_ * frame_n]
        //   sigs:        uint16_t [n_kmer_sizes_ * frame_n * sketch_size_]
        //   masks:       uint64_t [n_kmer_sizes_ * frame_n * mask_words_]
        const uint32_t* nrb = reinterpret_cast<const uint32_t*>(fbuf.data());
        const uint16_t* sigs = reinterpret_cast<const uint16_t*>(
            fbuf.data() + sizeof(uint32_t) * n_kmer_sizes_ * frame_n);
        const uint64_t* msks = reinterpret_cast<const uint64_t*>(
            fbuf.data() + sizeof(uint32_t) * n_kmer_sizes_ * frame_n
                        + sizeof(uint16_t) * n_kmer_sizes_ * frame_n * sketch_size_);

        for (size_t hi = hi_start; hi < hi_end; ++hi) {
            const uint32_t lr = hits[hi].local_row;
            const size_t   gi = static_cast<size_t>(ki) * frame_n + lr;
            const uint32_t global_pos = fi * v3_frame_sz_ + lr;

            SketchResult r{};
            r.sig           = sigs + gi * sketch_size_;
            r.mask          = msks + gi * mask_words_;
            r.n_real_bins   = nrb[gi];
            r.mask_words    = mask_words_;
            r.genome_length = v3_genome_lengths_[global_pos];
            r.sketch_size   = sketch_size_;
            r.kmer_size     = k;
            auto sliced = apply_slice(r, sz, mask_words_);
            if (sliced) cb(hits[hi].orig_idx, *sliced);
        }

        hi_start = hi_end;
    }
}

// ── SkchWriterMultiK ─────────────────────────────────────────────────────────

SkchWriterMultiK::SkchWriterMultiK(std::vector<uint32_t> kmer_sizes, uint32_t sketch_size,
                                    uint32_t syncmer_s, uint64_t seed1, uint64_t seed2,
                                    std::string spill_dir)
    : kmer_sizes_(std::move(kmer_sizes))
    , sketch_size_(sketch_size)
    , syncmer_s_(syncmer_s)
    , seed1_(seed1), seed2_(seed2)
    , mask_words_((sketch_size + 63) / 64)
    , spill_record_size_(kmer_sizes_.size() *
                         (sizeof(uint16_t) * sketch_size + sizeof(uint64_t) * ((sketch_size + 63) / 64)))
{
    std::sort(kmer_sizes_.begin(), kmer_sizes_.end());
    kmer_sizes_.erase(std::unique(kmer_sizes_.begin(), kmer_sizes_.end()), kmer_sizes_.end());
    if (kmer_sizes_.size() > 8)
        throw std::runtime_error("SkchWriterMultiK: at most 8 k values supported");

    n_real_bins_.resize(kmer_sizes_.size());

    if (spill_dir.empty()) {
        const char* env = std::getenv("GENOPACK_SPILL_DIR");
        if (env && *env) spill_dir = env;
    }
    if (!spill_dir.empty()) {
        char path[4096];
        std::snprintf(path, sizeof(path), "%s/genopack_skch_mk_XXXXXX", spill_dir.c_str());
        int fd = ::mkstemp(path);
        if (fd < 0) throw std::runtime_error(std::string("SkchWriterMultiK: mkstemp: ") + path);
        ::unlink(path);
        spill_fp_ = ::fdopen(fd, "w+b");
        if (!spill_fp_) { ::close(fd); throw std::runtime_error("SkchWriterMultiK: fdopen failed"); }
    } else {
        spill_fp_ = std::tmpfile();
        if (!spill_fp_) throw std::runtime_error("SkchWriterMultiK: cannot create spill tmpfile");
    }
}

SkchWriterMultiK::~SkchWriterMultiK() {
    if (spill_fp_) { std::fclose(spill_fp_); spill_fp_ = nullptr; }
}

void SkchWriterMultiK::add(GenomeId genome_id, uint64_t genome_length,
                            const std::vector<std::vector<uint16_t>>& sigs_per_k,
                            const std::vector<uint32_t>&              n_real_bins_per_k,
                            const std::vector<std::vector<uint64_t>>& masks_per_k)
{
    const size_t nk = kmer_sizes_.size();
    if (sigs_per_k.size() != nk || n_real_bins_per_k.size() != nk || masks_per_k.size() != nk)
        throw std::runtime_error("SkchWriterMultiK::add: k-count mismatch");

    ids_.push_back(genome_id);
    genome_lengths_.push_back(genome_length);
    for (size_t ki = 0; ki < nk; ++ki)
        n_real_bins_[ki].push_back(n_real_bins_per_k[ki]);

    // Spill record: sigs_k0 | sigs_k1 | ... | masks_k0 | masks_k1 | ...
    for (size_t ki = 0; ki < nk; ++ki)
        std::fwrite(sigs_per_k[ki].data(), sizeof(uint16_t), sketch_size_, spill_fp_);
    for (size_t ki = 0; ki < nk; ++ki)
        std::fwrite(masks_per_k[ki].data(), sizeof(uint64_t), mask_words_, spill_fp_);
}

SectionDesc SkchWriterMultiK::finalize(AppendWriter& writer, uint64_t section_id) {
    const uint32_t n  = static_cast<uint32_t>(ids_.size());
    const uint32_t nk = static_cast<uint32_t>(kmer_sizes_.size());

    std::vector<uint32_t> order(n);
    std::iota(order.begin(), order.end(), 0u);
    std::sort(order.begin(), order.end(),
              [&](uint32_t a, uint32_t b) { return ids_[a] < ids_[b]; });

    const size_t sig_bytes_k  = sizeof(uint16_t) * sketch_size_;
    const size_t mask_bytes_k = sizeof(uint64_t) * mask_words_;
    const uint32_t n_frames   = (n + SKCH_V3_FRAME_SIZE - 1) / SKCH_V3_FRAME_SIZE;
    const uint64_t section_start = writer.current_offset();

    // ── 1. Write uncompressed header ─────────────────────────────────────────
    SkchSeekHdr hdr{};
    hdr.magic        = SKCH_V3_MAGIC;
    hdr.n_frames     = n_frames;
    hdr.frame_size   = SKCH_V3_FRAME_SIZE;
    hdr.n_genomes    = n;
    hdr.sketch_size  = sketch_size_;
    hdr.n_kmer_sizes = nk;
    for (uint32_t i = 0; i < nk; ++i) hdr.kmer_sizes[i] = kmer_sizes_[i];
    hdr.syncmer_s    = syncmer_s_;
    hdr.mask_words   = mask_words_;
    hdr.seed1        = seed1_;
    hdr.seed2        = seed2_;
    writer.append(&hdr, sizeof(hdr));

    // ── 2. Write frame table placeholder (will be filled in below) ───────────
    const uint64_t frame_table_offset = writer.current_offset();
    std::vector<SkchFrameDesc> frame_descs(n_frames);
    writer.append(frame_descs.data(), sizeof(SkchFrameDesc) * n_frames);

    // ── 3. Write sorted genome_ids and genome_lengths (uncompressed) ─────────
    for (uint32_t i : order) writer.append(&ids_[i],            sizeof(uint64_t));
    for (uint32_t i : order) writer.append(&genome_lengths_[i], sizeof(uint64_t));

    // ── 4. Compress and write each frame, recording offset + size ────────────
    const size_t OUT_BUF = 4 << 20;
    std::vector<uint8_t> out_buf(OUT_BUF);
    std::vector<uint8_t> rec(spill_record_size_);

    for (uint32_t fi = 0; fi < n_frames; ++fi) {
        const uint32_t row_start = fi * SKCH_V3_FRAME_SIZE;
        const uint32_t row_end   = std::min(n, row_start + SKCH_V3_FRAME_SIZE);
        const uint32_t frame_n   = row_end - row_start;

        const size_t frame_raw_sz =
            sizeof(uint32_t) * nk * frame_n
            + sig_bytes_k  * nk * frame_n
            + mask_bytes_k * nk * frame_n;

        const uint64_t frame_start = writer.current_offset();
        frame_descs[fi].data_offset  = frame_start - section_start;
        frame_descs[fi].n_genomes    = frame_n;

        ZSTD_CStream* cs = ZSTD_createCStream();
        ZSTD_initCStream(cs, 3);
        ZSTD_CCtx_setPledgedSrcSize(cs, frame_raw_sz);

        auto compress = [&](const void* data, size_t size) {
            ZSTD_inBuffer zi{data, size, 0};
            while (zi.pos < zi.size) {
                ZSTD_outBuffer zo{out_buf.data(), out_buf.size(), 0};
                size_t r = ZSTD_compressStream(cs, &zo, &zi);
                if (ZSTD_isError(r))
                    throw std::runtime_error(std::string("SkchWriterMultiK V3: ") + ZSTD_getErrorName(r));
                if (zo.pos) writer.append(out_buf.data(), zo.pos);
            }
        };

        // n_real_bins planar
        for (uint32_t ki = 0; ki < nk; ++ki)
            for (uint32_t rank = row_start; rank < row_end; ++rank)
                compress(&n_real_bins_[ki][order[rank]], sizeof(uint32_t));

        // sigs planar
        for (uint32_t ki = 0; ki < nk; ++ki) {
            size_t sig_off = ki * sig_bytes_k;
            for (uint32_t rank = row_start; rank < row_end; ++rank) {
                long off = static_cast<long>(static_cast<size_t>(order[rank]) * spill_record_size_);
                if (std::fseek(spill_fp_, off, SEEK_SET) != 0)
                    throw std::runtime_error("SkchWriterMultiK V3: fseek sigs");
                if (std::fread(rec.data(), 1, spill_record_size_, spill_fp_) != spill_record_size_)
                    throw std::runtime_error("SkchWriterMultiK V3: fread sigs");
                compress(rec.data() + sig_off, sig_bytes_k);
            }
        }

        // masks planar
        size_t masks_base = nk * sig_bytes_k;
        for (uint32_t ki = 0; ki < nk; ++ki) {
            size_t mask_off = masks_base + ki * mask_bytes_k;
            for (uint32_t rank = row_start; rank < row_end; ++rank) {
                long off = static_cast<long>(static_cast<size_t>(order[rank]) * spill_record_size_);
                if (std::fseek(spill_fp_, off, SEEK_SET) != 0)
                    throw std::runtime_error("SkchWriterMultiK V3: fseek masks");
                if (std::fread(rec.data(), 1, spill_record_size_, spill_fp_) != spill_record_size_)
                    throw std::runtime_error("SkchWriterMultiK V3: fread masks");
                compress(rec.data() + mask_off, mask_bytes_k);
            }
        }

        size_t remaining = 1;
        while (remaining) {
            ZSTD_outBuffer zo{out_buf.data(), out_buf.size(), 0};
            remaining = ZSTD_endStream(cs, &zo);
            if (ZSTD_isError(remaining))
                throw std::runtime_error(std::string("SkchWriterMultiK V3 end: ") + ZSTD_getErrorName(remaining));
            if (zo.pos) writer.append(out_buf.data(), zo.pos);
        }
        ZSTD_freeCStream(cs);

        frame_descs[fi].compressed_size = static_cast<uint32_t>(writer.current_offset() - frame_start);
    }

    std::fclose(spill_fp_); spill_fp_ = nullptr;

    // ── 5. Seek back and write the completed frame table ─────────────────────
    const uint64_t section_end = writer.current_offset();
    writer.seek_to(frame_table_offset);
    writer.append(frame_descs.data(), sizeof(SkchFrameDesc) * n_frames);
    writer.seek_to(section_end);

    SectionDesc desc{};
    desc.type              = SEC_SKCH;
    desc.version           = 3;
    desc.section_id        = section_id;
    desc.file_offset       = section_start;
    desc.compressed_size   = section_end - section_start;
    desc.uncompressed_size = 0;
    desc.item_count        = n;
    desc.aux0              = sketch_size_;
    desc.aux1              = (nk > 0) ? kmer_sizes_[0] : 0;
    std::memset(desc.checksum, 0, sizeof(desc.checksum));
    return desc;
}

} // namespace genopack
