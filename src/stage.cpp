#include "genopack/stage.hpp"
#include "genopack/mmap_file.hpp"

#include <zstd.h>

#include <cstring>
#include <stdexcept>
#include <vector>

namespace genopack {

namespace {

constexpr int kStageZstdLevel = 3;

inline void put_u16(std::string& b, uint16_t v) { b.append(reinterpret_cast<const char*>(&v), 2); }
inline void put_u32(std::string& b, uint32_t v) { b.append(reinterpret_cast<const char*>(&v), 4); }
inline void put_u64(std::string& b, uint64_t v) { b.append(reinterpret_cast<const char*>(&v), 8); }

// Compress `payload`, prefixed with its u64 uncompressed size. The returned blob is
// [u64 usize][zstd frame]; block-table comp_size covers the whole blob.
std::string frame_compress(const std::string& payload) {
    size_t bound = ZSTD_compressBound(payload.size());
    std::string out;
    out.resize(8 + bound);
    uint64_t usize = payload.size();
    std::memcpy(out.data(), &usize, 8);
    size_t c = ZSTD_compress(out.data() + 8, bound,
                             payload.data(), payload.size(), kStageZstdLevel);
    if (ZSTD_isError(c))
        throw std::runtime_error(std::string("stage: ZSTD_compress: ") + ZSTD_getErrorName(c));
    out.resize(8 + c);
    return out;
}

std::string frame_decompress(const uint8_t* blob, uint64_t blob_len) {
    if (blob_len < 8) throw std::runtime_error("stage: truncated block");
    uint64_t usize;
    std::memcpy(&usize, blob, 8);
    std::string out;
    out.resize(usize);
    size_t d = ZSTD_decompress(out.data(), usize, blob + 8, blob_len - 8);
    if (ZSTD_isError(d))
        throw std::runtime_error(std::string("stage: ZSTD_decompress: ") + ZSTD_getErrorName(d));
    if (d != usize) throw std::runtime_error("stage: block size mismatch");
    return out;
}

struct Header {
    uint64_t magic;
    uint64_t n_genomes;
    uint64_t n_blocks;
    uint64_t block_table_off;
    uint64_t meta_off;
    uint64_t meta_csize;
};
static_assert(sizeof(Header) == 48, "stage header must be 48 bytes");

struct BlockEntry { uint64_t off, csize, n_recs; };
static_assert(sizeof(BlockEntry) == 24, "block table entry must be 24 bytes");

} // namespace

// ── StageWriter ───────────────────────────────────────────────────────────────

struct StageWriter::Impl {
    AppendWriter            w;
    uint64_t                block_bytes = 64ull << 20;
    std::string             block_buf;     // uncompressed records of the open block
    uint64_t                block_recs = 0;
    std::vector<BlockEntry> table;
    std::string             meta_buf;      // (acc, tax) directory, accumulated
    uint64_t                n_genomes = 0;

    void flush_block() {
        if (block_buf.empty()) return;
        std::string framed = frame_compress(block_buf);
        uint64_t off = w.append(framed.data(), framed.size());
        table.push_back({off, framed.size(), block_recs});
        block_buf.clear();
        block_recs = 0;
    }
};

StageWriter::StageWriter() : impl_(std::make_unique<Impl>()) {}
StageWriter::~StageWriter() = default;
uint64_t StageWriter::n_genomes() const { return impl_->n_genomes; }

void StageWriter::open(const std::filesystem::path& path, uint64_t block_bytes) {
    impl_->block_bytes = block_bytes ? block_bytes : (64ull << 20);
    impl_->w.create(path);
    Header h{};
    std::memset(&h, 0, sizeof(h));   // patched in finalize()
    impl_->w.append(&h, sizeof(h));
}

void StageWriter::add(const std::string& accession, const std::string& taxonomy,
                      const std::string& sequence) {
    auto& I = *impl_;
    // The .gstage record headers store lengths as u16 (accession, taxonomy) and u32
    // (sequence); reject inputs that would silently truncate and corrupt the stream.
    if (accession.size() > UINT16_MAX)
        throw std::runtime_error("stage: accession exceeds 65535 bytes: " + accession.substr(0, 64));
    if (taxonomy.size() > UINT16_MAX)
        throw std::runtime_error("stage: taxonomy exceeds 65535 bytes for " + accession);
    if (sequence.size() > UINT32_MAX)
        throw std::runtime_error("stage: sequence exceeds 4 GiB for " + accession);
    // Data record: u16 acc_len, u32 seq_len, acc, seq
    put_u16(I.block_buf, static_cast<uint16_t>(accession.size()));
    put_u32(I.block_buf, static_cast<uint32_t>(sequence.size()));
    I.block_buf.append(accession);
    I.block_buf.append(sequence);
    ++I.block_recs;
    // Meta record: u16 acc_len, u16 tax_len, acc, tax
    put_u16(I.meta_buf, static_cast<uint16_t>(accession.size()));
    put_u16(I.meta_buf, static_cast<uint16_t>(taxonomy.size()));
    I.meta_buf.append(accession);
    I.meta_buf.append(taxonomy);
    ++I.n_genomes;
    if (I.block_buf.size() >= I.block_bytes) I.flush_block();
}

void StageWriter::finalize() {
    auto& I = *impl_;
    I.flush_block();

    uint64_t block_table_off = I.w.current_offset();
    std::string tb;
    tb.reserve(I.table.size() * sizeof(BlockEntry));
    for (const auto& e : I.table) { put_u64(tb, e.off); put_u64(tb, e.csize); put_u64(tb, e.n_recs); }
    I.w.append(tb.data(), tb.size());

    uint64_t meta_off = I.w.current_offset();
    std::string meta_framed = frame_compress(I.meta_buf);
    I.w.append(meta_framed.data(), meta_framed.size());

    Header h;
    h.magic           = kStageMagic;
    h.n_genomes       = I.n_genomes;
    h.n_blocks        = I.table.size();
    h.block_table_off = block_table_off;
    h.meta_off        = meta_off;
    h.meta_csize      = meta_framed.size();
    I.w.write_at(0, &h, sizeof(h));
    I.w.flush();
    I.w.close();
}

// ── StageReader ───────────────────────────────────────────────────────────────

struct StageReader::Impl {
    MmapFileReader          m;
    Header                  h{};
    std::vector<BlockEntry> table;
};

StageReader::StageReader() : impl_(std::make_unique<Impl>()) {}
StageReader::~StageReader() = default;
uint64_t StageReader::n_genomes() const { return impl_->h.n_genomes; }
uint64_t StageReader::n_blocks() const { return impl_->h.n_blocks; }

void StageReader::open(const std::filesystem::path& path) {
    auto& I = *impl_;
    I.m.open(path);
    if (I.m.size() < sizeof(Header)) throw std::runtime_error("stage: file too small");
    std::memcpy(&I.h, I.m.data(), sizeof(Header));
    if (I.h.magic != kStageMagic)
        throw std::runtime_error("stage: bad magic (not a .gstage file)");
    auto sp = I.m.view(I.h.block_table_off, I.h.n_blocks * sizeof(BlockEntry));
    I.table.resize(I.h.n_blocks);
    if (I.h.n_blocks) std::memcpy(I.table.data(), sp.data(), I.h.n_blocks * sizeof(BlockEntry));
}

void StageReader::scan_meta(
    const std::function<void(const std::string&, const std::string&)>& cb) const {
    const auto& I = *impl_;
    auto sp = I.m.view(I.h.meta_off, I.h.meta_csize);
    std::string meta = frame_decompress(sp.data(), sp.size());
    const char* p   = meta.data();
    const char* end = p + meta.size();
    while (p + 4 <= end) {
        uint16_t al, tl;
        std::memcpy(&al, p, 2);
        std::memcpy(&tl, p + 2, 2);
        p += 4;
        if (p + al + tl > end) break;
        std::string acc(p, al); p += al;
        std::string tax(p, tl); p += tl;
        cb(acc, tax);
    }
}

void StageReader::scan(const std::function<void(StageRecord&&)>& cb) const {
    const auto& I = *impl_;
    for (const auto& e : I.table) {
        auto sp = I.m.view(e.off, e.csize);
        std::string blk = frame_decompress(sp.data(), sp.size());
        const char* p   = blk.data();
        const char* end = p + blk.size();
        for (uint64_t r = 0; r < e.n_recs; ++r) {
            if (p + 6 > end) throw std::runtime_error("stage: truncated record header");
            uint16_t al;
            uint32_t sl;
            std::memcpy(&al, p, 2);
            std::memcpy(&sl, p + 2, 4);
            p += 6;
            if (p + static_cast<size_t>(al) + sl > end)
                throw std::runtime_error("stage: truncated record body");
            StageRecord rec;
            rec.accession.assign(p, al); p += al;
            rec.sequence.assign(p, sl);  p += sl;
            cb(std::move(rec));
        }
    }
}

} // namespace genopack
