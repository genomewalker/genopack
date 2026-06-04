#include <genopack/stage.hpp>
#include <zstd.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>
#include <stdexcept>
#include <string>

namespace genopack {

static constexpr uint64_t STAGE_MAGIC   = 0x4753544147450001ULL; // "GSTAGE\x00\x01"
static constexpr uint64_t HEADER_BYTES  = 32;
static constexpr uint64_t BLOCKDESC_BYTES = 24;

// ── helpers ───────────────────────────────────────────────────────────────────

static void checked_write(int fd, const void* buf, size_t len) {
    const char* p = static_cast<const char*>(buf);
    while (len > 0) {
        ssize_t n = ::write(fd, p, len);
        if (n <= 0) throw std::runtime_error("stage write error");
        p += n; len -= static_cast<size_t>(n);
    }
}

static void checked_pread(int fd, void* buf, size_t len, off_t off) {
    char* p = static_cast<char*>(buf);
    while (len > 0) {
        ssize_t n = ::pread(fd, p, len, off);
        if (n <= 0) throw std::runtime_error("stage read error");
        p += n; off += n; len -= static_cast<size_t>(n);
    }
}

static void append_u16(std::vector<char>& buf, uint16_t v) {
    buf.push_back(static_cast<char>(v & 0xFF));
    buf.push_back(static_cast<char>(v >> 8));
}
static void append_u32(std::vector<char>& buf, uint32_t v) {
    for (int i = 0; i < 4; ++i) buf.push_back(static_cast<char>((v >> (8*i)) & 0xFF));
}
static void append_bytes(std::vector<char>& buf, const std::string& s) {
    buf.insert(buf.end(), s.begin(), s.end());
}

static uint16_t read_u16(const char* p) {
    return static_cast<uint16_t>(static_cast<uint8_t>(p[0]))
         | (static_cast<uint16_t>(static_cast<uint8_t>(p[1])) << 8);
}
static uint32_t read_u32(const char* p) {
    uint32_t v = 0;
    for (int i = 0; i < 4; ++i) v |= static_cast<uint32_t>(static_cast<uint8_t>(p[i])) << (8*i);
    return v;
}
static uint64_t read_u64(const char* p) {
    uint64_t v = 0;
    for (int i = 0; i < 8; ++i) v |= static_cast<uint64_t>(static_cast<uint8_t>(p[i])) << (8*i);
    return v;
}
static void write_u64(char* p, uint64_t v) {
    for (int i = 0; i < 8; ++i) p[i] = static_cast<char>((v >> (8*i)) & 0xFF);
}

// ── StageWriter ───────────────────────────────────────────────────────────────

void StageWriter::open(const std::filesystem::path& path, uint64_t block_bytes) {
    block_bytes_ = block_bytes;
    fd_ = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd_ < 0) throw std::runtime_error("stage: cannot open for write: " + path.string());

    // Reserve header + block table space: write placeholder zeros.
    // Block table size is unknown upfront; we'll patch it at finalize.
    // Strategy: write header now with n_blocks=0, seek past it as data grows,
    // then at finalize we rewrite the whole header+block-table prefix.
    char hdr[HEADER_BYTES] = {};
    checked_write(fd_, hdr, HEADER_BYTES);
    write_offset_ = HEADER_BYTES;
    data_start_   = HEADER_BYTES; // will be updated at finalize
}

void StageWriter::add(StageRecord r) {
    // Uppercase in-place
    for (char& c : r.sequence) if (c >= 'a' && c <= 'z') c -= 32;

    append_u16(pending_, static_cast<uint16_t>(r.accession.size()));
    append_u16(pending_, static_cast<uint16_t>(r.taxonomy.size()));
    append_u32(pending_, static_cast<uint32_t>(r.sequence.size()));
    append_bytes(pending_, r.accession);
    append_bytes(pending_, r.taxonomy);
    append_bytes(pending_, r.sequence);
    ++pending_n_;
    ++n_genomes_;

    if (pending_.size() >= block_bytes_)
        flush_block();
}

void StageWriter::flush_block() {
    if (pending_.empty()) return;

    size_t bound = ZSTD_compressBound(pending_.size());
    std::vector<char> cbuf(bound);
    size_t csize = ZSTD_compress(cbuf.data(), bound,
                                  pending_.data(), pending_.size(), 3);
    if (ZSTD_isError(csize))
        throw std::runtime_error(std::string("stage zstd: ") + ZSTD_getErrorName(csize));

    blocks_.push_back({write_offset_, csize, pending_n_});
    checked_write(fd_, cbuf.data(), csize);
    write_offset_ += csize;
    pending_.clear();
    pending_n_ = 0;
}

void StageWriter::finalize() {
    flush_block();

    const uint64_t n_blocks = blocks_.size();
    // Block table size
    const uint64_t btable_bytes = n_blocks * BLOCKDESC_BYTES;
    // data_start = HEADER_BYTES + btable_bytes (blocks are stored after block table)
    // But we already wrote blocks starting at HEADER_BYTES.
    // We need to shift: prepend the block table. Do this by rewriting the whole file header.
    // Simpler: write block table at current end, then write header at start via pwrite.

    // Write block table at current end of file
    const uint64_t btable_offset = write_offset_;
    for (const auto& b : blocks_) {
        char buf[BLOCKDESC_BYTES];
        write_u64(buf,     b.offset);
        write_u64(buf + 8, b.csize);
        write_u64(buf +16, b.n_genomes);
        checked_write(fd_, buf, BLOCKDESC_BYTES);
    }
    write_offset_ += btable_bytes;

    // Write final header at byte 0 via pwrite
    char hdr[HEADER_BYTES];
    write_u64(hdr,      STAGE_MAGIC);
    write_u64(hdr + 8,  static_cast<uint64_t>(n_genomes_));
    write_u64(hdr + 16, n_blocks);
    write_u64(hdr + 24, btable_offset); // store btable_offset in flags field
    if (::pwrite(fd_, hdr, HEADER_BYTES, 0) != HEADER_BYTES)
        throw std::runtime_error("stage: failed to write header");

    ::close(fd_);
    fd_ = -1;
}

// ── StageReader ───────────────────────────────────────────────────────────────

void StageReader::open(const std::filesystem::path& path) {
    path_ = path;
    int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) throw std::runtime_error("stage: cannot open: " + path.string());

    char hdr[HEADER_BYTES];
    if (::pread(fd, hdr, HEADER_BYTES, 0) != HEADER_BYTES) {
        ::close(fd);
        throw std::runtime_error("stage: truncated header");
    }
    if (read_u64(hdr) != STAGE_MAGIC) {
        ::close(fd);
        throw std::runtime_error("stage: bad magic — not a .gstage file");
    }
    n_genomes_         = read_u64(hdr + 8);
    n_blocks_          = read_u64(hdr + 16);
    uint64_t bt_offset = read_u64(hdr + 24);

    blocks_.resize(n_blocks_);
    std::vector<char> btbuf(n_blocks_ * BLOCKDESC_BYTES);
    if (n_blocks_ > 0) {
        checked_pread(fd, btbuf.data(), btbuf.size(), static_cast<off_t>(bt_offset));
        for (uint64_t i = 0; i < n_blocks_; ++i) {
            const char* p = btbuf.data() + i * BLOCKDESC_BYTES;
            blocks_[i] = {read_u64(p), read_u64(p + 8), read_u64(p + 16)};
        }
    }
    ::close(fd);
    data_start_ = HEADER_BYTES;
}

void StageReader::scan_blocks(uint64_t block_start, uint64_t block_end,
                               const std::function<void(StageRecord&&)>& cb) const {
    if (block_start >= n_blocks_) return;
    block_end = std::min(block_end, n_blocks_);

    int fd = ::open(path_.c_str(), O_RDONLY);
    if (fd < 0) throw std::runtime_error("stage: cannot open for read: " + path_.string());

    for (uint64_t bi = block_start; bi < block_end; ++bi) {
        const auto& bd = blocks_[bi];
        std::vector<char> cbuf(bd.csize);
        checked_pread(fd, cbuf.data(), bd.csize, static_cast<off_t>(bd.offset));

        size_t dbound = ZSTD_getFrameContentSize(cbuf.data(), bd.csize);
        if (dbound == ZSTD_CONTENTSIZE_ERROR || dbound == ZSTD_CONTENTSIZE_UNKNOWN)
            dbound = bd.n_genomes * (512 + 512 + 5000000); // generous fallback
        std::vector<char> dbuf(dbound);
        size_t dsize = ZSTD_decompress(dbuf.data(), dbound, cbuf.data(), bd.csize);
        if (ZSTD_isError(dsize))
            throw std::runtime_error(std::string("stage decompress: ") + ZSTD_getErrorName(dsize));

        const char* p   = dbuf.data();
        const char* end = p + dsize;
        for (uint64_t gi = 0; gi < bd.n_genomes; ++gi) {
            if (p + 8 > end) throw std::runtime_error("stage: block truncated");
            uint16_t acc_len = read_u16(p);     p += 2;
            uint16_t tax_len = read_u16(p);     p += 2;
            uint32_t seq_len = read_u32(p);     p += 4;
            if (p + acc_len + tax_len + seq_len > end)
                throw std::runtime_error("stage: record overrun");
            StageRecord r;
            r.accession.assign(p, acc_len); p += acc_len;
            r.taxonomy .assign(p, tax_len); p += tax_len;
            r.sequence .assign(p, seq_len); p += seq_len;
            cb(std::move(r));
        }
    }
    ::close(fd);
}

} // namespace genopack
