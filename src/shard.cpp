#include <genopack/shard.hpp>
#include <zstd.h>
#include <zdict.h>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace genopack {

// ── ShardWriterV1::Impl ─────────────────────────────────────────────────────────

struct ShardWriterV1::Impl {
    std::filesystem::path   path;
    ShardId                 shard_id;
    uint32_t                cluster_id;
    Config                  cfg;

    // Genome records accumulated before finalize
    struct GenomeRecord {
        GenomeId genome_id;
        uint32_t flags;
        uint64_t oph_fingerprint;
        std::string compressed_blob;  // zstd-compressed FASTA
    };
    std::vector<GenomeRecord> genomes;

    // Training samples for shared dictionary
    std::vector<std::string> dict_samples;
    std::string shared_dict;

    // Streaming context for dict-trained compression
    ZSTD_CCtx* cctx = nullptr;

    Impl(const std::filesystem::path& p, ShardId sid, uint32_t cid, Config c)
        : path(p), shard_id(sid), cluster_id(cid), cfg(c)
    {
        cctx = ZSTD_createCCtx();
        if (!cctx) throw std::runtime_error("ZSTD_createCCtx failed");
        apply_zstd_params();
    }

    ~Impl() {
        ZSTD_freeCCtx(cctx);
    }

    void apply_zstd_params() {
        ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, cfg.zstd_level);
        if (cfg.use_long_match) {
            ZSTD_CCtx_setParameter(cctx, ZSTD_c_enableLongDistanceMatching, 1);
            ZSTD_CCtx_setParameter(cctx, ZSTD_c_windowLog, cfg.zstd_wlog);
        }
    }

    void train_dict_if_ready() {
        if (!cfg.train_dict) return;
        if (dict_samples.size() < std::min(cfg.dict_samples, genomes.size())) return;
        if (!shared_dict.empty()) return;  // already trained

        // Collect sample buffers + sizes
        std::vector<size_t> sample_sizes;
        std::string concat;
        for (const auto& s : dict_samples) {
            sample_sizes.push_back(s.size());
            concat += s;
        }

        shared_dict.resize(cfg.dict_size);
        size_t dict_size = ZDICT_trainFromBuffer(
            shared_dict.data(), cfg.dict_size,
            concat.data(), sample_sizes.data(), sample_sizes.size());

        if (ZDICT_isError(dict_size)) {
            spdlog::warn("ZDICT training failed: {} — proceeding without dict",
                         ZDICT_getErrorName(dict_size));
            shared_dict.clear();
        } else {
            shared_dict.resize(dict_size);
            spdlog::debug("Shard {}: trained {}B dictionary from {} samples",
                          shard_id, dict_size, dict_samples.size());
            // Load dictionary into context
            ZSTD_CCtx_loadDictionary(cctx, shared_dict.data(), shared_dict.size());
        }
    }

    std::string compress(const char* data, size_t len) {
        size_t bound = ZSTD_compressBound(len);
        std::string out(bound, '\0');
        size_t csize = ZSTD_compress2(cctx, out.data(), bound, data, len);
        if (ZSTD_isError(csize))
            throw std::runtime_error(std::string("ZSTD_compress2: ") + ZSTD_getErrorName(csize));
        out.resize(csize);
        return out;
    }
};

ShardWriterV1::ShardWriterV1(const std::filesystem::path& path, ShardId shard_id,
                          uint32_t cluster_id, Config cfg)
    : impl_(std::make_unique<Impl>(path, shard_id, cluster_id, cfg))
{}

ShardWriterV1::~ShardWriterV1() = default;

void ShardWriterV1::add_genome(GenomeId id, uint64_t oph_fingerprint,
                              const char* fasta_data, size_t fasta_len, uint32_t flags)
{
    // Collect training samples for dict
    if (impl_->cfg.train_dict &&
        impl_->dict_samples.size() < impl_->cfg.dict_samples) {
        impl_->dict_samples.emplace_back(fasta_data, fasta_len);
    }

    // Train dict once we have enough samples
    impl_->train_dict_if_ready();

    std::string blob = impl_->compress(fasta_data, fasta_len);

    impl_->genomes.push_back({
        id, flags, oph_fingerprint, std::move(blob)
    });
}

void ShardWriterV1::finalize() {
    std::ofstream f(impl_->path, std::ios::binary | std::ios::trunc);
    if (!f) throw std::runtime_error("Cannot open shard for writing: " + impl_->path.string());

    const size_t n = impl_->genomes.size();

    // Compute layout offsets
    const size_t header_size   = sizeof(ShardHeader);
    const size_t dir_size      = n * sizeof(GenomeDirEntry);
    const size_t dict_size     = impl_->shared_dict.size();
    const size_t genome_dir_offset     = header_size;
    const size_t dict_offset           = genome_dir_offset + dir_size;
    const size_t blob_area_offset      = dict_offset + dict_size;

    // Compute blob offsets
    std::vector<GenomeDirEntry> dir(n);
    size_t blob_cursor = blob_area_offset;
    for (size_t i = 0; i < n; ++i) {
        const auto& g = impl_->genomes[i];
        dir[i].genome_id       = g.genome_id;
        dir[i]._reserved0      = 0;
        dir[i].flags           = g.flags;
        dir[i].oph_fingerprint = g.oph_fingerprint;
        dir[i].blob_offset     = blob_cursor;
        dir[i].blob_len_cmp    = static_cast<uint32_t>(g.compressed_blob.size());
        dir[i].blob_len_raw    = 0;  // placeholder (we store compressed only for now)
        dir[i].n_checkpoints   = 0;  // v1: no intra-genome checkpoints
        dir[i].checkpoint_idx  = 0;
        blob_cursor += g.compressed_blob.size();
    }

    const size_t checkpoint_area_offset = blob_cursor;
    const size_t footer_offset          = checkpoint_area_offset;  // v1: no checkpoints

    // Write header (placeholder, overwritten at end)
    ShardHeader hdr{};
    hdr.magic                  = GPKS_MAGIC;
    hdr.version                = FORMAT_VERSION;
    hdr.shard_id               = impl_->shard_id;
    hdr.n_genomes              = static_cast<uint32_t>(n);
    hdr.n_deleted              = static_cast<uint32_t>(
        std::count_if(impl_->genomes.begin(), impl_->genomes.end(),
                      [](const auto& g){ return g.flags & GenomeMeta::FLAG_DELETED; }));
    hdr.cluster_id             = impl_->cluster_id;
    hdr.dict_size              = static_cast<uint32_t>(dict_size);
    hdr.genome_dir_offset      = genome_dir_offset;
    hdr.dict_offset            = dict_offset;
    hdr.blob_area_offset       = blob_area_offset;
    hdr.checkpoint_area_offset = checkpoint_area_offset;
    hdr.footer_offset          = footer_offset;

    f.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));

    // Write genome directory
    f.write(reinterpret_cast<const char*>(dir.data()), dir_size);

    // Write shared dictionary
    if (!impl_->shared_dict.empty())
        f.write(impl_->shared_dict.data(), dict_size);

    // Write genome blobs
    for (const auto& g : impl_->genomes)
        f.write(g.compressed_blob.data(), g.compressed_blob.size());

    spdlog::info("Shard {}: wrote {} genomes, {} bytes",
                 impl_->shard_id, n, static_cast<size_t>(f.tellp()));
}

size_t ShardWriterV1::n_genomes() const { return impl_->genomes.size(); }
size_t ShardWriterV1::compressed_size() const {
    size_t total = 0;
    for (const auto& g : impl_->genomes) total += g.compressed_blob.size();
    return total;
}

// ── ShardReaderV1::Impl ─────────────────────────────────────────────────────────

struct ShardReaderV1::Impl {
    std::filesystem::path      path;
    std::vector<char>          file_data;  // entire shard file (mmap-style read)
    const ShardHeader*         header = nullptr;
    const GenomeDirEntry*      dir    = nullptr;
    std::string                dict;       // shared dictionary

    ZSTD_DDict* ddict = nullptr;
    ZSTD_DCtx*  dctx  = nullptr;

    void open(const std::filesystem::path& p) {
        path = p;
        std::ifstream f(p, std::ios::binary | std::ios::ate);
        if (!f) throw std::runtime_error("Cannot open shard: " + p.string());
        size_t size = static_cast<size_t>(f.tellg());
        f.seekg(0);
        file_data.resize(size);
        f.read(file_data.data(), size);
        if (!f) throw std::runtime_error("Read error: " + p.string());

        header = reinterpret_cast<const ShardHeader*>(file_data.data());
        if (header->magic != GPKS_MAGIC)
            throw std::runtime_error("Invalid shard magic: " + p.string());
        if (header->version != FORMAT_VERSION)
            throw std::runtime_error("Shard version mismatch: " + p.string());

        dir = reinterpret_cast<const GenomeDirEntry*>(
            file_data.data() + header->genome_dir_offset);

        dctx = ZSTD_createDCtx();
        if (!dctx) throw std::runtime_error("ZSTD_createDCtx failed");

        if (header->dict_size > 0) {
            const char* dict_ptr = file_data.data() + header->dict_offset;
            ddict = ZSTD_createDDict(dict_ptr, header->dict_size);
            if (!ddict) throw std::runtime_error("ZSTD_createDDict failed");
        }
    }

    void close() {
        if (ddict) { ZSTD_freeDDict(ddict); ddict = nullptr; }
        if (dctx)  { ZSTD_freeDCtx(dctx);  dctx  = nullptr; }
        file_data.clear();
        header = nullptr;
        dir    = nullptr;
    }

    std::string decompress_blob(const GenomeDirEntry& e) const {
        const char* src = file_data.data() + e.blob_offset;
        size_t src_size = e.blob_len_cmp;
        size_t raw_size = ZSTD_getFrameContentSize(src, src_size);
        if (raw_size == ZSTD_CONTENTSIZE_UNKNOWN || raw_size == ZSTD_CONTENTSIZE_ERROR) {
            // Fallback: decompress into a growing buffer
            std::string out(4 << 20, '\0');
            size_t written = ddict
                ? ZSTD_decompress_usingDDict(dctx, out.data(), out.size(), src, src_size, ddict)
                : ZSTD_decompressDCtx(dctx, out.data(), out.size(), src, src_size);
            while (ZSTD_isError(written) && out.size() < (256 << 20)) {
                out.resize(out.size() * 2);
                written = ddict
                    ? ZSTD_decompress_usingDDict(dctx, out.data(), out.size(), src, src_size, ddict)
                    : ZSTD_decompressDCtx(dctx, out.data(), out.size(), src, src_size);
            }
            if (ZSTD_isError(written))
                throw std::runtime_error(std::string("ZSTD decompress: ") + ZSTD_getErrorName(written));
            out.resize(written);
            return out;
        }
        std::string out(raw_size, '\0');
        size_t written = ddict
            ? ZSTD_decompress_usingDDict(dctx, out.data(), raw_size, src, src_size, ddict)
            : ZSTD_decompressDCtx(dctx, out.data(), raw_size, src, src_size);
        if (ZSTD_isError(written))
            throw std::runtime_error(std::string("ZSTD decompress: ") + ZSTD_getErrorName(written));
        return out;
    }

    const GenomeDirEntry* find_genome(GenomeId id) const {
        // Linear scan for now; O(n) fine for typical shard sizes.
        // Binary search on genome_id column is a future optimization.
        for (size_t i = 0; i < header->n_genomes; ++i) {
            if (dir[i].genome_id == id) return &dir[i];
        }
        return nullptr;
    }
};

ShardReaderV1::ShardReaderV1()              : impl_(std::make_unique<Impl>()) {}
ShardReaderV1::~ShardReaderV1()             = default;
ShardReaderV1::ShardReaderV1(ShardReaderV1&&) noexcept = default;
ShardReaderV1& ShardReaderV1::operator=(ShardReaderV1&&) noexcept = default;

void ShardReaderV1::open(const std::filesystem::path& path) { impl_->open(path); }
void ShardReaderV1::close() { impl_->close(); }
ShardId  ShardReaderV1::shard_id()  const { return impl_->header->shard_id; }
uint32_t ShardReaderV1::cluster_id() const { return impl_->header->cluster_id; }
size_t   ShardReaderV1::n_genomes() const { return impl_->header->n_genomes; }

std::string ShardReaderV1::fetch_genome(GenomeId id) const {
    const GenomeDirEntry* e = impl_->find_genome(id);
    if (!e) throw std::runtime_error("genome_id not found in shard");
    if (e->flags & GenomeMeta::FLAG_DELETED) return {};
    return impl_->decompress_blob(*e);
}

const GenomeDirEntry* ShardReaderV1::dir_begin() const { return impl_->dir; }
const GenomeDirEntry* ShardReaderV1::dir_end()   const {
    return impl_->dir + impl_->header->n_genomes;
}

// ── ShardWriter::Impl ───────────────────────────────────────────────────────

struct ShardWriter::Impl {
    uint32_t shard_id;
    uint32_t cluster_id;
    Config   cfg;

    struct PendingGenome {
        GenomeId genome_id;
        uint64_t oph_fingerprint;
        uint32_t flags;
        uint32_t meta_row_id = 0;
        uint64_t raw_offset;   // byte offset into raw_buffer
        uint32_t raw_len;
    };
    std::vector<PendingGenome> pending;

    std::vector<char> raw_buffer;     // append-only store of raw FASTAs
    uint64_t          total_raw_bytes = 0;

    std::vector<std::string> dict_samples;
    std::string              shared_dict;

    Impl(uint32_t sid, uint32_t cid, Config c)
        : shard_id(sid), cluster_id(cid), cfg(c)
    {}
};

ShardWriter::ShardWriter(uint32_t shard_id, uint32_t cluster_id, Config cfg)
    : impl_(std::make_unique<Impl>(shard_id, cluster_id, cfg))
{}

ShardWriter::~ShardWriter() = default;

void ShardWriter::add_genome(GenomeId id, uint64_t oph_fingerprint,
                                const char* fasta_data, size_t fasta_len,
                                uint32_t flags)
{
    static constexpr size_t MAX_SAMPLE = 65536;

    if (impl_->cfg.train_dict &&
        impl_->dict_samples.size() < impl_->cfg.dict_samples) {
        impl_->dict_samples.emplace_back(fasta_data, std::min(fasta_len, MAX_SAMPLE));
    }

    Impl::PendingGenome pg;
    pg.genome_id       = id;
    pg.oph_fingerprint = oph_fingerprint;
    pg.flags           = flags;
    pg.meta_row_id     = 0;
    pg.raw_offset      = static_cast<uint64_t>(impl_->raw_buffer.size());
    pg.raw_len         = static_cast<uint32_t>(fasta_len);
    impl_->raw_buffer.insert(impl_->raw_buffer.end(), fasta_data, fasta_data + fasta_len);
    impl_->total_raw_bytes += fasta_len;
    impl_->pending.push_back(std::move(pg));
}

uint64_t ShardWriter::finalize(AppendWriter& writer) {
    const size_t n = impl_->pending.size();

    // ── 1. Train dictionary ───────────────────────────────────────────────────
    bool use_dict = false;
    if (impl_->cfg.train_dict && !impl_->dict_samples.empty()) {
        std::string concat;
        std::vector<size_t> sizes;
        concat.reserve(impl_->dict_samples.size() * 32768);
        for (const auto& s : impl_->dict_samples) {
            sizes.push_back(s.size());
            concat += s;
        }
        impl_->shared_dict.resize(impl_->cfg.dict_size);
        size_t trained = ZDICT_trainFromBuffer(
            impl_->shared_dict.data(), impl_->cfg.dict_size,
            concat.data(), sizes.data(), sizes.size());
        if (ZDICT_isError(trained)) {
            spdlog::warn("ShardWriter shard {}: ZDICT training failed: {} — no dict",
                         impl_->shard_id, ZDICT_getErrorName(trained));
            impl_->shared_dict.clear();
        } else {
            impl_->shared_dict.resize(trained);
            spdlog::debug("ShardWriter shard {}: trained {}B dict from {} samples",
                          impl_->shard_id, trained, impl_->dict_samples.size());
            use_dict = true;
        }
    }

    // ── 2. Set up ZSTD compression context ───────────────────────────────────
    ZSTD_CCtx* cctx = ZSTD_createCCtx();
    if (!cctx) throw std::runtime_error("ZSTD_createCCtx failed");

    ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, impl_->cfg.zstd_level);
    if (impl_->cfg.use_long_match) {
        ZSTD_CCtx_setParameter(cctx, ZSTD_c_enableLongDistanceMatching, 1);
        ZSTD_CCtx_setParameter(cctx, ZSTD_c_windowLog, impl_->cfg.zstd_wlog);
    }
    if (use_dict)
        ZSTD_CCtx_loadDictionary(cctx, impl_->shared_dict.data(), impl_->shared_dict.size());

    // ── 3. Compress all pending genomes ──────────────────────────────────────
    struct CompressedBlob {
        std::vector<char> data;
        uint32_t          raw_len;
    };
    std::vector<CompressedBlob> blobs;
    blobs.reserve(n);

    for (const auto& pg : impl_->pending) {
        const char* src     = impl_->raw_buffer.data() + pg.raw_offset;
        size_t      src_len = pg.raw_len;
        size_t      bound   = ZSTD_compressBound(src_len);
        CompressedBlob cb;
        cb.data.resize(bound);
        cb.raw_len = static_cast<uint32_t>(src_len);
        size_t csize = ZSTD_compress2(cctx, cb.data.data(), bound, src, src_len);
        if (ZSTD_isError(csize)) {
            ZSTD_freeCCtx(cctx);
            throw std::runtime_error(std::string("ZSTD_compress2: ") + ZSTD_getErrorName(csize));
        }
        cb.data.resize(csize);
        blobs.push_back(std::move(cb));
    }
    ZSTD_freeCCtx(cctx);

    // ── 4. Compute layout (all offsets relative to section start) ─────────────
    const uint64_t header_size       = sizeof(ShardHeaderV2);
    const uint64_t dir_size          = n * sizeof(GenomeDirEntryV2);
    const uint64_t dict_bytes        = use_dict ? impl_->shared_dict.size() : 0;
    const uint64_t genome_dir_offset = header_size;
    const uint64_t dict_offset       = genome_dir_offset + dir_size;
    const uint64_t blob_area_offset  = dict_offset + dict_bytes;

    // Compute per-genome blob offsets (relative to blob_area_offset)
    std::vector<GenomeDirEntryV2> dir(n);
    uint64_t blob_cursor      = 0;
    uint64_t total_compressed = 0;
    for (size_t i = 0; i < n; ++i) {
        const auto& pg = impl_->pending[i];
        const auto& cb = blobs[i];
        dir[i].genome_id       = pg.genome_id;
        dir[i].oph_fingerprint = pg.oph_fingerprint;
        dir[i].blob_offset     = blob_cursor;
        dir[i].blob_len_cmp    = static_cast<uint32_t>(cb.data.size());
        dir[i].blob_len_raw    = cb.raw_len;
        dir[i].checkpoint_idx  = 0;
        dir[i].n_checkpoints   = 0;
        dir[i].flags           = pg.flags;
        dir[i].meta_row_id     = pg.meta_row_id;
        std::memset(dir[i].reserved, 0, sizeof(dir[i].reserved));
        blob_cursor      += cb.data.size();
        total_compressed += cb.data.size();
    }

    // ── 5. Write section ──────────────────────────────────────────────────────
    const uint64_t section_start = writer.current_offset();

    // 5a. Build and write header
    ShardHeaderV2 hdr{};
    hdr.magic                  = GPKS_MAGIC;
    hdr.version                = 2;
    hdr.flags                  = 0;
    hdr.shard_id               = impl_->shard_id;
    hdr.cluster_id             = impl_->cluster_id;
    hdr.n_genomes              = static_cast<uint32_t>(n);
    hdr.n_deleted              = static_cast<uint32_t>(
        std::count_if(impl_->pending.begin(), impl_->pending.end(),
                      [](const Impl::PendingGenome& pg){
                          return pg.flags & GenomeMeta::FLAG_DELETED;
                      }));
    hdr.codec                  = use_dict ? 1u : 0u;
    hdr.dict_size              = static_cast<uint32_t>(dict_bytes);
    hdr.genome_dir_offset      = genome_dir_offset;
    hdr.dict_offset            = dict_offset;
    hdr.blob_area_offset       = blob_area_offset;
    hdr.shard_raw_bp           = impl_->total_raw_bytes;
    hdr.shard_compressed_bytes = total_compressed;
    std::memset(hdr.checksum, 0, sizeof(hdr.checksum));
    std::memset(hdr.reserved, 0, sizeof(hdr.reserved));
    writer.append(&hdr, sizeof(hdr));

    // 5b. Directory
    writer.append(dir.data(), dir_size);

    // 5c. Dictionary
    if (use_dict)
        writer.append(impl_->shared_dict.data(), dict_bytes);

    // 5d. Blob area
    for (const auto& cb : blobs)
        writer.append(cb.data.data(), cb.data.size());

    spdlog::info("ShardWriter shard {}: wrote {} genomes, raw {}B, compressed {}B",
                 impl_->shard_id, n, impl_->total_raw_bytes, total_compressed);

    return section_start;
}

size_t ShardWriter::n_genomes()   const { return impl_->pending.size(); }
size_t ShardWriter::n_bytes_raw() const { return impl_->total_raw_bytes; }

// ── ShardReader::Impl ───────────────────────────────────────────────────────

struct ShardReader::Impl {
    const uint8_t*          base_         = nullptr;  // start of shard section
    uint64_t                section_size_ = 0;
    const ShardHeaderV2*    header_       = nullptr;
    const GenomeDirEntryV2* dir_          = nullptr;

    std::vector<uint8_t>    owned_data_;  // non-empty when file-based open

    ZSTD_DDict* ddict_ = nullptr;
    ZSTD_DCtx*  dctx_  = nullptr;

    void setup(const uint8_t* section_base, uint64_t section_size) {
        base_         = section_base;
        section_size_ = section_size;

        header_ = reinterpret_cast<const ShardHeaderV2*>(base_);
        if (header_->magic != GPKS_MAGIC)
            throw std::runtime_error("ShardReader: invalid shard magic");
        if (header_->version != 2)
            throw std::runtime_error("ShardReader: expected shard version 2");

        dir_ = reinterpret_cast<const GenomeDirEntryV2*>(
            base_ + header_->genome_dir_offset);

        dctx_ = ZSTD_createDCtx();
        if (!dctx_) throw std::runtime_error("ZSTD_createDCtx failed");

        if (header_->dict_size > 0) {
            const uint8_t* dict_ptr = base_ + header_->dict_offset;
            ddict_ = ZSTD_createDDict(dict_ptr, header_->dict_size);
            if (!ddict_) throw std::runtime_error("ZSTD_createDDict failed");
        }
    }

    void reset() {
        if (ddict_) { ZSTD_freeDDict(ddict_); ddict_ = nullptr; }
        if (dctx_)  { ZSTD_freeDCtx(dctx_);  dctx_  = nullptr; }
        base_         = nullptr;
        section_size_ = 0;
        header_       = nullptr;
        dir_          = nullptr;
        owned_data_.clear();
    }

    const GenomeDirEntryV2* find_genome(GenomeId id) const {
        for (uint32_t i = 0; i < header_->n_genomes; ++i) {
            if (dir_[i].genome_id == id) return &dir_[i];
        }
        return nullptr;
    }

    std::string decompress_blob(const GenomeDirEntryV2& e) const {
        const uint8_t* src      = base_ + header_->blob_area_offset + e.blob_offset;
        const size_t   src_size = e.blob_len_cmp;
        const size_t   raw_size = e.blob_len_raw;

        std::string out(raw_size, '\0');
        size_t written = ddict_
            ? ZSTD_decompress_usingDDict(dctx_, out.data(), raw_size, src, src_size, ddict_)
            : ZSTD_decompressDCtx(dctx_, out.data(), raw_size, src, src_size);

        if (ZSTD_isError(written)) {
            // Fall back to frame-size detection if stored raw_len was wrong
            size_t frame_size = ZSTD_getFrameContentSize(src, src_size);
            if (frame_size != ZSTD_CONTENTSIZE_UNKNOWN &&
                frame_size != ZSTD_CONTENTSIZE_ERROR &&
                frame_size != raw_size) {
                out.resize(frame_size);
                written = ddict_
                    ? ZSTD_decompress_usingDDict(dctx_, out.data(), frame_size,
                                                 src, src_size, ddict_)
                    : ZSTD_decompressDCtx(dctx_, out.data(), frame_size, src, src_size);
            }
            if (ZSTD_isError(written))
                throw std::runtime_error(std::string("ShardReader decompress: ") +
                                         ZSTD_getErrorName(written));
            out.resize(written);
        }
        return out;
    }
};

ShardReader::ShardReader() = default;

ShardReader::~ShardReader() {
    if (impl_) impl_->reset();
}

ShardReader::ShardReader(ShardReader&&) noexcept = default;
ShardReader& ShardReader::operator=(ShardReader&&) noexcept = default;

void ShardReader::open(const uint8_t* base, uint64_t section_offset, uint64_t section_size) {
    if (!impl_) impl_ = std::make_unique<Impl>();
    else        impl_->reset();
    impl_->setup(base + section_offset, section_size);
}

void ShardReader::open_file(const std::filesystem::path& path) {
    if (!impl_) impl_ = std::make_unique<Impl>();
    else        impl_->reset();

    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) throw std::runtime_error("ShardReader: cannot open: " + path.string());
    size_t size = static_cast<size_t>(f.tellg());
    f.seekg(0);
    impl_->owned_data_.resize(size);
    f.read(reinterpret_cast<char*>(impl_->owned_data_.data()), size);
    if (!f) throw std::runtime_error("ShardReader: read error: " + path.string());

    impl_->setup(impl_->owned_data_.data(), size);
}

bool     ShardReader::is_open()   const { return impl_ && impl_->header_ != nullptr; }
uint32_t ShardReader::shard_id()  const { return impl_->header_->shard_id; }
uint32_t ShardReader::n_genomes() const { return impl_->header_->n_genomes; }

std::string ShardReader::fetch_genome(GenomeId id) const {
    const GenomeDirEntryV2* e = impl_->find_genome(id);
    if (!e) throw std::runtime_error("ShardReader: genome_id not found");
    if (e->flags & GenomeMeta::FLAG_DELETED) return {};
    return impl_->decompress_blob(*e);
}

const GenomeDirEntryV2* ShardReader::dir_begin() const { return impl_->dir_; }
const GenomeDirEntryV2* ShardReader::dir_end()   const {
    return impl_->dir_ + impl_->header_->n_genomes;
}

} // namespace genopack
