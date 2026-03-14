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

// ── ShardWriter::Impl ─────────────────────────────────────────────────────────

struct ShardWriter::Impl {
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

ShardWriter::ShardWriter(const std::filesystem::path& path, ShardId shard_id,
                          uint32_t cluster_id, Config cfg)
    : impl_(std::make_unique<Impl>(path, shard_id, cluster_id, cfg))
{}

ShardWriter::~ShardWriter() = default;

void ShardWriter::add_genome(GenomeId id, uint64_t oph_fingerprint,
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

void ShardWriter::finalize() {
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

size_t ShardWriter::n_genomes() const { return impl_->genomes.size(); }
size_t ShardWriter::compressed_size() const {
    size_t total = 0;
    for (const auto& g : impl_->genomes) total += g.compressed_blob.size();
    return total;
}

// ── ShardReader::Impl ─────────────────────────────────────────────────────────

struct ShardReader::Impl {
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

ShardReader::ShardReader()              : impl_(std::make_unique<Impl>()) {}
ShardReader::~ShardReader()             = default;
ShardReader::ShardReader(ShardReader&&) noexcept = default;
ShardReader& ShardReader::operator=(ShardReader&&) noexcept = default;

void ShardReader::open(const std::filesystem::path& path) { impl_->open(path); }
void ShardReader::close() { impl_->close(); }
ShardId  ShardReader::shard_id()  const { return impl_->header->shard_id; }
uint32_t ShardReader::cluster_id() const { return impl_->header->cluster_id; }
size_t   ShardReader::n_genomes() const { return impl_->header->n_genomes; }

std::string ShardReader::fetch_genome(GenomeId id) const {
    const GenomeDirEntry* e = impl_->find_genome(id);
    if (!e) throw std::runtime_error("genome_id not found in shard");
    if (e->flags & GenomeMeta::FLAG_DELETED) return {};
    return impl_->decompress_blob(*e);
}

const GenomeDirEntry* ShardReader::dir_begin() const { return impl_->dir; }
const GenomeDirEntry* ShardReader::dir_end()   const {
    return impl_->dir + impl_->header->n_genomes;
}

} // namespace genopack
