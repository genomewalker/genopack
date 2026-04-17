#include <genopack/skch.hpp>
#include <zstd.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <spdlog/spdlog.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace genopack {

// ── Spill dir resolution ─────────────────────────────────────────────────────
// Priority: explicit arg > GENOPACK_SKETCH_SPILL_DIR > GENOPACK_SPILL_DIR > std::tmpfile
static std::string resolve_spill_dir(std::string spill_dir) {
    if (!spill_dir.empty()) return spill_dir;
    if (const char* env = std::getenv("GENOPACK_SKETCH_SPILL_DIR"); env && *env) return env;
    if (const char* env = std::getenv("GENOPACK_SPILL_DIR");        env && *env) return env;
    return {};
}

static FILE* open_spill(const std::string& spill_dir, const char* tag) {
    if (!spill_dir.empty()) {
        char path[4096];
        std::snprintf(path, sizeof(path), "%s/genopack_%s_XXXXXX", spill_dir.c_str(), tag);
        int fd = ::mkstemp(path);
        if (fd < 0) throw std::runtime_error(std::string("SkchWriter: mkstemp failed: ") + path);
        ::unlink(path);
        FILE* fp = ::fdopen(fd, "w+b");
        if (!fp) { ::close(fd); throw std::runtime_error("SkchWriter: fdopen failed"); }
        return fp;
    }
    FILE* fp = std::tmpfile();
    if (!fp) throw std::runtime_error("SkchWriter: cannot create spill tmpfile");
    return fp;
}

// ── SkchWriter (V4, single-k, dual-seed) ─────────────────────────────────────

SkchWriter::SkchWriter(uint32_t sketch_size, uint32_t kmer_size,
                       uint32_t syncmer_s, uint64_t seed1, uint64_t seed2,
                       std::string spill_dir)
    : sketch_size_(sketch_size)
    , kmer_size_(kmer_size)
    , syncmer_s_(syncmer_s)
    , seed1_(seed1)
    , seed2_(seed2)
    , mask_words_((sketch_size + 63) / 64)
    , record_size_(2 * sizeof(uint16_t) * sketch_size
                   + sizeof(uint64_t) * ((sketch_size + 63) / 64))
{
    spill_fp_ = open_spill(resolve_spill_dir(std::move(spill_dir)), "skch_v4");
}

SkchWriter::~SkchWriter() {
    if (spill_fp_) { std::fclose(spill_fp_); spill_fp_ = nullptr; }
}

void SkchWriter::add(GenomeId genome_id,
                     const std::vector<uint16_t>& oph_sig1,
                     const std::vector<uint16_t>& oph_sig2,
                     uint32_t n_real_bins,
                     uint64_t genome_length,
                     const std::vector<uint64_t>& mask) {
    if (oph_sig1.size() != sketch_size_ || oph_sig2.size() != sketch_size_)
        throw std::runtime_error("SkchWriter::add: sig size mismatch");
    if (mask.size() != mask_words_)
        throw std::runtime_error("SkchWriter::add: mask size mismatch");

    ids_.push_back(genome_id);
    n_real_bins_.push_back(n_real_bins);
    genome_lengths_.push_back(genome_length);

    // Spill record: sig1 | sig2 | mask
    std::fwrite(oph_sig1.data(), sizeof(uint16_t), sketch_size_, spill_fp_);
    std::fwrite(oph_sig2.data(), sizeof(uint16_t), sketch_size_, spill_fp_);
    std::fwrite(mask.data(),     sizeof(uint64_t), mask_words_,  spill_fp_);
}

SectionDesc SkchWriter::finalize(AppendWriter& writer, uint64_t section_id) {
    const uint32_t n = static_cast<uint32_t>(ids_.size());

    std::vector<uint32_t> order(n);
    std::iota(order.begin(), order.end(), 0u);
    std::sort(order.begin(), order.end(),
              [&](uint32_t a, uint32_t b) { return ids_[a] < ids_[b]; });

    const size_t sig_bytes   = sizeof(uint16_t) * sketch_size_;
    const size_t mask_bytes  = sizeof(uint64_t) * mask_words_;
    const uint32_t n_frames  = (n + SKCH_V4_FRAME_SIZE - 1) / SKCH_V4_FRAME_SIZE;
    const uint64_t section_start = writer.current_offset();

    // Header (SKCH_V4_MAGIC, single k).
    SkchSeekHdr hdr{};
    hdr.magic        = SKCH_V4_MAGIC;
    hdr.n_frames     = n_frames;
    hdr.frame_size   = SKCH_V4_FRAME_SIZE;
    hdr.n_genomes    = n;
    hdr.sketch_size  = sketch_size_;
    hdr.n_kmer_sizes = 1;
    hdr.kmer_sizes[0] = kmer_size_;
    hdr.syncmer_s    = syncmer_s_;
    hdr.mask_words   = mask_words_;
    hdr.seed1        = seed1_;
    hdr.seed2        = seed2_;
    writer.append(&hdr, sizeof(hdr));

    const uint64_t frame_table_offset = writer.current_offset();
    std::vector<SkchFrameDesc> frame_descs(n_frames);
    writer.append(frame_descs.data(), sizeof(SkchFrameDesc) * n_frames);

    // Sorted genome_ids and genome_lengths (uncompressed).
    for (uint32_t i : order) writer.append(&ids_[i],            sizeof(uint64_t));
    for (uint32_t i : order) writer.append(&genome_lengths_[i], sizeof(uint64_t));

    const size_t OUT_BUF = 4 << 20;
    std::vector<uint8_t> out_buf(OUT_BUF);
    std::vector<uint8_t> rec(record_size_);

    for (uint32_t fi = 0; fi < n_frames; ++fi) {
        const uint32_t row_start = fi * SKCH_V4_FRAME_SIZE;
        const uint32_t row_end   = std::min(n, row_start + SKCH_V4_FRAME_SIZE);
        const uint32_t frame_n   = row_end - row_start;

        const size_t frame_raw_sz =
            sizeof(uint32_t) * frame_n
            + 2 * sig_bytes * frame_n
            + mask_bytes * frame_n;

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
                    throw std::runtime_error(std::string("SkchWriter V4: ") + ZSTD_getErrorName(r));
                if (zo.pos) writer.append(out_buf.data(), zo.pos);
            }
        };

        // n_real_bins: one plane (NK=1).
        for (uint32_t rank = row_start; rank < row_end; ++rank)
            compress(&n_real_bins_[order[rank]], sizeof(uint32_t));

        // sigs1 plane
        for (uint32_t rank = row_start; rank < row_end; ++rank) {
            long off = static_cast<long>(static_cast<size_t>(order[rank]) * record_size_);
            if (std::fseek(spill_fp_, off, SEEK_SET) != 0)
                throw std::runtime_error("SkchWriter V4: fseek sig1");
            if (std::fread(rec.data(), 1, record_size_, spill_fp_) != record_size_)
                throw std::runtime_error("SkchWriter V4: fread sig1");
            compress(rec.data(), sig_bytes);
        }
        // sigs2 plane
        for (uint32_t rank = row_start; rank < row_end; ++rank) {
            long off = static_cast<long>(static_cast<size_t>(order[rank]) * record_size_);
            if (std::fseek(spill_fp_, off, SEEK_SET) != 0)
                throw std::runtime_error("SkchWriter V4: fseek sig2");
            if (std::fread(rec.data(), 1, record_size_, spill_fp_) != record_size_)
                throw std::runtime_error("SkchWriter V4: fread sig2");
            compress(rec.data() + sig_bytes, sig_bytes);
        }
        // masks plane
        for (uint32_t rank = row_start; rank < row_end; ++rank) {
            long off = static_cast<long>(static_cast<size_t>(order[rank]) * record_size_);
            if (std::fseek(spill_fp_, off, SEEK_SET) != 0)
                throw std::runtime_error("SkchWriter V4: fseek mask");
            if (std::fread(rec.data(), 1, record_size_, spill_fp_) != record_size_)
                throw std::runtime_error("SkchWriter V4: fread mask");
            compress(rec.data() + 2 * sig_bytes, mask_bytes);
        }

        size_t remaining = 1;
        while (remaining) {
            ZSTD_outBuffer zo{out_buf.data(), out_buf.size(), 0};
            remaining = ZSTD_endStream(cs, &zo);
            if (ZSTD_isError(remaining))
                throw std::runtime_error(std::string("SkchWriter V4 end: ") + ZSTD_getErrorName(remaining));
            if (zo.pos) writer.append(out_buf.data(), zo.pos);
        }
        ZSTD_freeCStream(cs);

        frame_descs[fi].compressed_size = static_cast<uint32_t>(writer.current_offset() - frame_start);
    }

    std::fclose(spill_fp_); spill_fp_ = nullptr;

    const uint64_t section_end = writer.current_offset();
    writer.seek_to(frame_table_offset);
    writer.append(frame_descs.data(), sizeof(SkchFrameDesc) * n_frames);
    writer.seek_to(section_end);

    SectionDesc desc{};
    desc.type              = SEC_SKCH;
    desc.version           = 4;
    desc.section_id        = section_id;
    desc.file_offset       = section_start;
    desc.compressed_size   = section_end - section_start;
    desc.uncompressed_size = 0;
    desc.item_count        = n;
    desc.aux0              = sketch_size_;
    desc.aux1              = kmer_size_;
    std::memset(desc.checksum, 0, sizeof(desc.checksum));
    return desc;
}

// ── SkchWriterMultiK (V4, multi-k, dual-seed) ────────────────────────────────

SkchWriterMultiK::SkchWriterMultiK(std::vector<uint32_t> kmer_sizes, uint32_t sketch_size,
                                    uint32_t syncmer_s, uint64_t seed1, uint64_t seed2,
                                    std::string spill_dir)
    : kmer_sizes_(std::move(kmer_sizes))
    , sketch_size_(sketch_size)
    , syncmer_s_(syncmer_s)
    , seed1_(seed1), seed2_(seed2)
    , mask_words_((sketch_size + 63) / 64)
    , spill_record_size_(kmer_sizes_.size() *
                         (2 * sizeof(uint16_t) * sketch_size
                          + sizeof(uint64_t) * ((sketch_size + 63) / 64)))
{
    std::sort(kmer_sizes_.begin(), kmer_sizes_.end());
    kmer_sizes_.erase(std::unique(kmer_sizes_.begin(), kmer_sizes_.end()), kmer_sizes_.end());
    if (kmer_sizes_.size() > 8)
        throw std::runtime_error("SkchWriterMultiK: at most 8 k values supported");
    n_real_bins_.resize(kmer_sizes_.size());

    // Rough tmpfs budget warning: 2 * sig_bytes * n_k per genome.
    const size_t per_genome = 2 * sizeof(uint16_t) * sketch_size * kmer_sizes_.size();
    if (per_genome > (256u << 20))
        spdlog::warn("SkchWriterMultiK: per-genome spill {} MB is large; "
                     "set GENOPACK_SKETCH_SPILL_DIR to a big-enough volume",
                     per_genome >> 20);

    spill_fp_ = open_spill(resolve_spill_dir(std::move(spill_dir)), "skch_v4_mk");
}

SkchWriterMultiK::~SkchWriterMultiK() {
    if (spill_fp_) { std::fclose(spill_fp_); spill_fp_ = nullptr; }
}

void SkchWriterMultiK::add(GenomeId genome_id, uint64_t genome_length,
                            const std::vector<std::vector<uint16_t>>& sigs1_per_k,
                            const std::vector<std::vector<uint16_t>>& sigs2_per_k,
                            const std::vector<uint32_t>&              n_real_bins_per_k,
                            const std::vector<std::vector<uint64_t>>& masks_per_k)
{
    const size_t nk = kmer_sizes_.size();
    if (sigs1_per_k.size() != nk || sigs2_per_k.size() != nk
        || n_real_bins_per_k.size() != nk || masks_per_k.size() != nk)
        throw std::runtime_error("SkchWriterMultiK::add: k-count mismatch");

    ids_.push_back(genome_id);
    genome_lengths_.push_back(genome_length);
    for (size_t ki = 0; ki < nk; ++ki)
        n_real_bins_[ki].push_back(n_real_bins_per_k[ki]);

    // Spill record: sigs1_k0 | sigs2_k0 | masks_k0 | sigs1_k1 | sigs2_k1 | masks_k1 | ...
    // (contiguous per-k triples — finalize can read one k at a time per genome).
    for (size_t ki = 0; ki < nk; ++ki) {
        if (sigs1_per_k[ki].size() != sketch_size_ || sigs2_per_k[ki].size() != sketch_size_
            || masks_per_k[ki].size() != mask_words_)
            throw std::runtime_error("SkchWriterMultiK::add: per-k size mismatch");
        std::fwrite(sigs1_per_k[ki].data(), sizeof(uint16_t), sketch_size_, spill_fp_);
        std::fwrite(sigs2_per_k[ki].data(), sizeof(uint16_t), sketch_size_, spill_fp_);
        std::fwrite(masks_per_k[ki].data(), sizeof(uint64_t), mask_words_,  spill_fp_);
    }
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
    const size_t k_triple_sz  = 2 * sig_bytes_k + mask_bytes_k;  // per-k bytes in spill
    const uint32_t n_frames   = (n + SKCH_V4_FRAME_SIZE - 1) / SKCH_V4_FRAME_SIZE;
    const uint64_t section_start = writer.current_offset();

    SkchSeekHdr hdr{};
    hdr.magic        = SKCH_V4_MAGIC;
    hdr.n_frames     = n_frames;
    hdr.frame_size   = SKCH_V4_FRAME_SIZE;
    hdr.n_genomes    = n;
    hdr.sketch_size  = sketch_size_;
    hdr.n_kmer_sizes = nk;
    for (uint32_t i = 0; i < nk; ++i) hdr.kmer_sizes[i] = kmer_sizes_[i];
    hdr.syncmer_s    = syncmer_s_;
    hdr.mask_words   = mask_words_;
    hdr.seed1        = seed1_;
    hdr.seed2        = seed2_;
    writer.append(&hdr, sizeof(hdr));

    const uint64_t frame_table_offset = writer.current_offset();
    std::vector<SkchFrameDesc> frame_descs(n_frames);
    writer.append(frame_descs.data(), sizeof(SkchFrameDesc) * n_frames);

    for (uint32_t i : order) writer.append(&ids_[i],            sizeof(uint64_t));
    for (uint32_t i : order) writer.append(&genome_lengths_[i], sizeof(uint64_t));

    const size_t OUT_BUF = 4 << 20;
    std::vector<uint8_t> out_buf(OUT_BUF);
    std::vector<uint8_t> rec(spill_record_size_);

    for (uint32_t fi = 0; fi < n_frames; ++fi) {
        const uint32_t row_start = fi * SKCH_V4_FRAME_SIZE;
        const uint32_t row_end   = std::min(n, row_start + SKCH_V4_FRAME_SIZE);
        const uint32_t frame_n   = row_end - row_start;

        const size_t frame_raw_sz =
            sizeof(uint32_t) * nk * frame_n
            + 2 * sig_bytes_k  * nk * frame_n
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
                    throw std::runtime_error(std::string("SkchWriterMultiK V4: ") + ZSTD_getErrorName(r));
                if (zo.pos) writer.append(out_buf.data(), zo.pos);
            }
        };

        // n_real_bins planar by k.
        for (uint32_t ki = 0; ki < nk; ++ki)
            for (uint32_t rank = row_start; rank < row_end; ++rank)
                compress(&n_real_bins_[ki][order[rank]], sizeof(uint32_t));

        // sigs1 planar by k.
        for (uint32_t ki = 0; ki < nk; ++ki) {
            const size_t k_off = ki * k_triple_sz;  // offset of this k's triple in the record
            for (uint32_t rank = row_start; rank < row_end; ++rank) {
                long off = static_cast<long>(static_cast<size_t>(order[rank]) * spill_record_size_);
                if (std::fseek(spill_fp_, off, SEEK_SET) != 0)
                    throw std::runtime_error("SkchWriterMultiK V4: fseek sig1");
                if (std::fread(rec.data(), 1, spill_record_size_, spill_fp_) != spill_record_size_)
                    throw std::runtime_error("SkchWriterMultiK V4: fread sig1");
                compress(rec.data() + k_off, sig_bytes_k);
            }
        }
        // sigs2 planar by k.
        for (uint32_t ki = 0; ki < nk; ++ki) {
            const size_t k_off = ki * k_triple_sz + sig_bytes_k;
            for (uint32_t rank = row_start; rank < row_end; ++rank) {
                long off = static_cast<long>(static_cast<size_t>(order[rank]) * spill_record_size_);
                if (std::fseek(spill_fp_, off, SEEK_SET) != 0)
                    throw std::runtime_error("SkchWriterMultiK V4: fseek sig2");
                if (std::fread(rec.data(), 1, spill_record_size_, spill_fp_) != spill_record_size_)
                    throw std::runtime_error("SkchWriterMultiK V4: fread sig2");
                compress(rec.data() + k_off, sig_bytes_k);
            }
        }
        // masks planar by k.
        for (uint32_t ki = 0; ki < nk; ++ki) {
            const size_t k_off = ki * k_triple_sz + 2 * sig_bytes_k;
            for (uint32_t rank = row_start; rank < row_end; ++rank) {
                long off = static_cast<long>(static_cast<size_t>(order[rank]) * spill_record_size_);
                if (std::fseek(spill_fp_, off, SEEK_SET) != 0)
                    throw std::runtime_error("SkchWriterMultiK V4: fseek mask");
                if (std::fread(rec.data(), 1, spill_record_size_, spill_fp_) != spill_record_size_)
                    throw std::runtime_error("SkchWriterMultiK V4: fread mask");
                compress(rec.data() + k_off, mask_bytes_k);
            }
        }

        size_t remaining = 1;
        while (remaining) {
            ZSTD_outBuffer zo{out_buf.data(), out_buf.size(), 0};
            remaining = ZSTD_endStream(cs, &zo);
            if (ZSTD_isError(remaining))
                throw std::runtime_error(std::string("SkchWriterMultiK V4 end: ") + ZSTD_getErrorName(remaining));
            if (zo.pos) writer.append(out_buf.data(), zo.pos);
        }
        ZSTD_freeCStream(cs);

        frame_descs[fi].compressed_size = static_cast<uint32_t>(writer.current_offset() - frame_start);
    }

    std::fclose(spill_fp_); spill_fp_ = nullptr;

    const uint64_t section_end = writer.current_offset();
    writer.seek_to(frame_table_offset);
    writer.append(frame_descs.data(), sizeof(SkchFrameDesc) * n_frames);
    writer.seek_to(section_end);

    SectionDesc desc{};
    desc.type              = SEC_SKCH;
    desc.version           = 4;
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

// ── SkchReader (V4) ──────────────────────────────────────────────────────────

// Rejection helper: called whenever we encounter a non-SKC4 section.
[[noreturn]] static void reject_old_version(uint32_t magic_seen) {
    std::string hint;
    if (magic_seen == SKCH_V3_MAGIC) hint = " (found SKC3/V3)";
    throw std::runtime_error("genopack SKCH version unsupported (V4 required; "
                             "rebuild archive with genopack build)" + hint);
}

void SkchReader::open(const uint8_t* base, uint64_t offset, uint64_t compressed_size) {
    cdata_    = base + offset;
    cdata_sz_ = compressed_size;

    if (compressed_size < sizeof(SkchSeekHdr))
        throw std::runtime_error("SkchReader: section too small for V4 header");

    uint32_t maybe_magic;
    std::memcpy(&maybe_magic, cdata_, sizeof(uint32_t));
    if (maybe_magic != SKCH_V4_MAGIC) reject_old_version(maybe_magic);

    const auto* h = reinterpret_cast<const SkchSeekHdr*>(cdata_);
    n_genomes_    = h->n_genomes;
    sketch_size_  = h->sketch_size;
    mask_words_   = h->mask_words;
    n_kmer_sizes_ = std::min(h->n_kmer_sizes, 8u);
    for (uint32_t i = 0; i < 8; ++i) kmer_sizes_[i] = h->kmer_sizes[i];
    kmer_size_    = (n_kmer_sizes_ > 0) ? kmer_sizes_[0] : 0;
    syncmer_s_    = h->syncmer_s;
    seed1_        = h->seed1;
    seed2_        = h->seed2;
    frame_sz_     = h->frame_size;
    section_base_ = cdata_;

    const auto* ft = reinterpret_cast<const SkchFrameDesc*>(cdata_ + sizeof(SkchSeekHdr));
    frames_.assign(ft, ft + h->n_frames);

    const uint8_t* ids_ptr = cdata_ + sizeof(SkchSeekHdr)
                             + sizeof(SkchFrameDesc) * h->n_frames;
    id_index_.resize(n_genomes_);
    std::memcpy(id_index_.data(), ids_ptr, sizeof(uint64_t) * n_genomes_);

    const uint8_t* len_ptr = ids_ptr + sizeof(uint64_t) * n_genomes_;
    genome_lengths_.resize(n_genomes_);
    std::memcpy(genome_lengths_.data(), len_ptr, sizeof(uint64_t) * n_genomes_);
}

std::pair<uint32_t, std::vector<uint32_t>>
SkchReader::peek_params(const uint8_t* base, uint64_t offset, uint64_t compressed_sz) {
    const uint8_t* src = base + offset;
    if (compressed_sz < sizeof(SkchSeekHdr))
        throw std::runtime_error("SkchReader::peek_params: section too small for V4 header");

    uint32_t maybe_magic;
    std::memcpy(&maybe_magic, src, sizeof(uint32_t));
    if (maybe_magic != SKCH_V4_MAGIC) reject_old_version(maybe_magic);

    const auto* h = reinterpret_cast<const SkchSeekHdr*>(src);
    uint32_t n = std::min(h->n_kmer_sizes, 8u);
    std::vector<uint32_t> ks(h->kmer_sizes, h->kmer_sizes + n);
    return {4, ks};
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
                                                     uint32_t     /*stored_mask_words*/) {
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
    return r;
}

// sketch_for(): single-genome access. Decompresses exactly the one frame
// that contains the genome, copies the relevant slice out to a per-thread
// static buffer, and returns a SketchResult referring to it.
//
// NOTE: this is for convenience only; batch lookups via sketch_for_ids()
// are far more efficient and the main code path.
namespace {
thread_local std::vector<uint8_t> tl_skch_frame_buf;
thread_local std::vector<uint16_t> tl_skch_sig1_buf;
thread_local std::vector<uint16_t> tl_skch_sig2_buf;
thread_local std::vector<uint64_t> tl_skch_mask_buf;
thread_local uint32_t tl_skch_n_real_bins = 0;
thread_local uint64_t tl_skch_genome_len  = 0;
}

std::optional<SketchResult> SkchReader::sketch_for(GenomeId genome_id) const {
    if (n_genomes_ == 0) return std::nullopt;
    return sketch_for(genome_id, kmer_size_, sketch_size_);
}

std::optional<SketchResult> SkchReader::sketch_for(GenomeId genome_id,
                                                    uint32_t k,
                                                    uint32_t requested_sz) const {
    if (requested_sz > sketch_size_) return std::nullopt;

    uint32_t ki = UINT32_MAX;
    for (uint32_t i = 0; i < n_kmer_sizes_; ++i)
        if (kmer_sizes_[i] == k) { ki = i; break; }
    if (ki == UINT32_MAX) return std::nullopt;

    const uint32_t pos = find_genome_pos(id_index_, genome_id);
    if (pos == UINT32_MAX) return std::nullopt;

    const uint32_t fi    = pos / frame_sz_;
    const uint32_t lr    = pos % frame_sz_;
    const SkchFrameDesc& fd = frames_[fi];
    const uint8_t* csrc  = section_base_ + fd.data_offset;

    unsigned long long raw_sz = ZSTD_getFrameContentSize(csrc, fd.compressed_size);
    if (raw_sz == ZSTD_CONTENTSIZE_ERROR || raw_sz == ZSTD_CONTENTSIZE_UNKNOWN)
        return std::nullopt;
    tl_skch_frame_buf.resize(static_cast<size_t>(raw_sz));
    if (ZSTD_isError(ZSTD_decompress(tl_skch_frame_buf.data(), tl_skch_frame_buf.size(),
                                     csrc, fd.compressed_size)))
        return std::nullopt;

    const uint32_t nk = n_kmer_sizes_;
    const uint32_t fn = fd.n_genomes;
    const uint32_t* nrb = reinterpret_cast<const uint32_t*>(tl_skch_frame_buf.data());
    const uint16_t* sigs1 = reinterpret_cast<const uint16_t*>(
        tl_skch_frame_buf.data() + sizeof(uint32_t) * nk * fn);
    const uint16_t* sigs2 = reinterpret_cast<const uint16_t*>(
        tl_skch_frame_buf.data() + sizeof(uint32_t) * nk * fn
                                 + sizeof(uint16_t) * nk * fn * sketch_size_);
    const uint64_t* msks = reinterpret_cast<const uint64_t*>(
        tl_skch_frame_buf.data() + sizeof(uint32_t) * nk * fn
                                 + 2 * sizeof(uint16_t) * nk * fn * sketch_size_);

    const size_t gj = static_cast<size_t>(ki) * fn + lr;
    tl_skch_sig1_buf.assign(sigs1 + gj * sketch_size_,
                            sigs1 + (gj + 1) * sketch_size_);
    tl_skch_sig2_buf.assign(sigs2 + gj * sketch_size_,
                            sigs2 + (gj + 1) * sketch_size_);
    tl_skch_mask_buf.assign(msks + gj * mask_words_,
                            msks + (gj + 1) * mask_words_);
    tl_skch_n_real_bins = nrb[gj];
    tl_skch_genome_len  = genome_lengths_[pos];

    SketchResult r{};
    r.sig           = tl_skch_sig1_buf.data();
    r.sig2          = tl_skch_sig2_buf.data();
    r.mask          = tl_skch_mask_buf.data();
    r.n_real_bins   = tl_skch_n_real_bins;
    r.mask_words    = mask_words_;
    r.genome_length = tl_skch_genome_len;
    r.sketch_size   = sketch_size_;
    r.kmer_size     = k;
    return apply_slice(r, requested_sz, mask_words_);
}

void SkchReader::sketch_for_ids(const std::vector<GenomeId>& sorted_ids,
                                 uint32_t k, uint32_t sz,
                                 const SketchCallback& cb) const
{
    if (sorted_ids.empty()) return;

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

    struct Hit { uint32_t frame_idx; uint32_t local_row; size_t orig_idx; };
    std::vector<Hit> hits;
    hits.reserve(sorted_ids.size());
    for (size_t i = 0; i < sorted_ids.size(); ++i) {
        uint32_t pos = find_genome_pos(id_index_, sorted_ids[i]);
        if (pos == UINT32_MAX) continue;
        hits.push_back({pos / frame_sz_, pos % frame_sz_, i});
    }

    struct FrameGroup { uint32_t fi; size_t hi_start; size_t hi_end; };
    std::vector<FrameGroup> groups;
    groups.reserve(frames_.size());
    {
        size_t s = 0;
        while (s < hits.size()) {
            const uint32_t fi = hits[s].frame_idx;
            size_t e = s;
            while (e < hits.size() && hits[e].frame_idx == fi) ++e;
            groups.push_back({fi, s, e});
            s = e;
        }
    }

    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t gi_idx = 0; gi_idx < groups.size(); ++gi_idx) {
        const uint32_t fi     = groups[gi_idx].fi;
        const size_t   hi_s   = groups[gi_idx].hi_start;
        const size_t   hi_e   = groups[gi_idx].hi_end;

        const SkchFrameDesc& fd = frames_[fi];
        const uint32_t frame_n  = fd.n_genomes;
        const uint8_t* csrc     = section_base_ + fd.data_offset;

        unsigned long long raw_sz = ZSTD_getFrameContentSize(csrc, fd.compressed_size);
        if (raw_sz == ZSTD_CONTENTSIZE_ERROR || raw_sz == ZSTD_CONTENTSIZE_UNKNOWN)
            continue;
        std::vector<uint8_t> fbuf(static_cast<size_t>(raw_sz));
        if (ZSTD_isError(ZSTD_decompress(fbuf.data(), fbuf.size(), csrc, fd.compressed_size)))
            continue;

        // V4 frame layout (planar by k):
        //   n_real_bins: uint32_t [NK * frame_n]
        //   sigs1:       uint16_t [NK * frame_n * sketch_size_]
        //   sigs2:       uint16_t [NK * frame_n * sketch_size_]
        //   masks:       uint64_t [NK * frame_n * mask_words_]
        const uint32_t nk = n_kmer_sizes_;
        const uint32_t* nrb = reinterpret_cast<const uint32_t*>(fbuf.data());
        const uint16_t* sigs1 = reinterpret_cast<const uint16_t*>(
            fbuf.data() + sizeof(uint32_t) * nk * frame_n);
        const uint16_t* sigs2 = reinterpret_cast<const uint16_t*>(
            fbuf.data() + sizeof(uint32_t) * nk * frame_n
                        + sizeof(uint16_t) * nk * frame_n * sketch_size_);
        const uint64_t* msks = reinterpret_cast<const uint64_t*>(
            fbuf.data() + sizeof(uint32_t) * nk * frame_n
                        + 2 * sizeof(uint16_t) * nk * frame_n * sketch_size_);

        for (size_t hi = hi_s; hi < hi_e; ++hi) {
            const uint32_t lr = hits[hi].local_row;
            const size_t   gj = static_cast<size_t>(ki) * frame_n + lr;
            const uint32_t global_pos = fi * frame_sz_ + lr;

            SketchResult r{};
            r.sig           = sigs1 + gj * sketch_size_;
            r.sig2          = sigs2 + gj * sketch_size_;
            r.mask          = msks  + gj * mask_words_;
            r.n_real_bins   = nrb[gj];
            r.mask_words    = mask_words_;
            r.genome_length = genome_lengths_[global_pos];
            r.sketch_size   = sketch_size_;
            r.kmer_size     = k;
            auto sliced = apply_slice(r, sz, mask_words_);
            if (sliced) cb(hits[hi].orig_idx, *sliced);
        }
    }
}

} // namespace genopack
