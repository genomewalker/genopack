#include <genopack/ntdb.hpp>
#include <zstd.h>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace genopack {

static std::string read_file(const std::filesystem::path& p) {
    std::ifstream f(p, std::ios::binary);
    if (!f) throw std::runtime_error("ntdb: cannot open " + p.string());
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

static std::string zstd_compress(const std::string& src, int level = 3) {
    size_t bound = ZSTD_compressBound(src.size());
    std::string dst(bound, '\0');
    size_t n = ZSTD_compress(dst.data(), bound, src.data(), src.size(), level);
    if (ZSTD_isError(n))
        throw std::runtime_error(std::string("ntdb: zstd compress error: ") + ZSTD_getErrorName(n));
    dst.resize(n);
    return dst;
}

static std::string zstd_decompress(const uint8_t* src, size_t src_size, size_t raw_size) {
    std::string dst(raw_size, '\0');
    size_t n = ZSTD_decompress(dst.data(), raw_size, src, src_size);
    if (ZSTD_isError(n))
        throw std::runtime_error(std::string("ntdb: zstd decompress error: ") + ZSTD_getErrorName(n));
    return dst;
}

// ── NtdbWriter ────────────────────────────────────────────────────────────────

void NtdbWriter::load(const std::filesystem::path& taxdump_dir, uint64_t taxdump_date) {
    nodes_raw_     = read_file(taxdump_dir / "nodes.dmp");
    names_raw_     = read_file(taxdump_dir / "names.dmp");
    taxdump_date_  = taxdump_date;
}

SectionDesc NtdbWriter::finalize(AppendWriter& writer, uint64_t section_id) {
    if (nodes_raw_.empty() || names_raw_.empty())
        throw std::runtime_error("NtdbWriter: call load() before finalize()");

    std::string nodes_z = zstd_compress(nodes_raw_);
    std::string names_z = zstd_compress(names_raw_);

    NtdbHeader hdr{};
    hdr.magic           = SEC_NTDB;
    hdr.version         = 1;
    hdr.taxdump_date    = taxdump_date_;
    hdr.nodes_raw_size  = nodes_raw_.size();
    hdr.nodes_zstd_size = nodes_z.size();
    hdr.names_raw_size  = names_raw_.size();
    hdr.names_zstd_size = names_z.size();

    uint64_t section_start = writer.current_offset();

    writer.append(&hdr,          sizeof(hdr));
    writer.append(nodes_z.data(),nodes_z.size());
    writer.append(names_z.data(),names_z.size());

    uint64_t payload_size = writer.current_offset() - section_start;

    SectionDesc desc{};
    desc.type              = SEC_NTDB;
    desc.version           = 1;
    desc.section_id        = section_id;
    desc.file_offset       = section_start;
    desc.compressed_size   = payload_size;
    desc.uncompressed_size = sizeof(hdr) + nodes_raw_.size() + names_raw_.size();
    return desc;
}

// ── NtdbReader ────────────────────────────────────────────────────────────────

void NtdbReader::open(const uint8_t* base, uint64_t offset, uint64_t size) {
    if (size < sizeof(NtdbHeader))
        throw std::runtime_error("NtdbReader: section too small");

    const uint8_t* payload = base + offset;
    hdr_ = reinterpret_cast<const NtdbHeader*>(payload);

    if (hdr_->magic != SEC_NTDB)
        throw std::runtime_error("NtdbReader: bad magic");

    nodes_blob_ = payload + sizeof(NtdbHeader);
    names_blob_ = nodes_blob_ + hdr_->nodes_zstd_size;
}

std::string NtdbReader::nodes_dmp() const {
    if (!hdr_) return {};
    return zstd_decompress(nodes_blob_, hdr_->nodes_zstd_size,
                           static_cast<size_t>(hdr_->nodes_raw_size));
}

std::string NtdbReader::names_dmp() const {
    if (!hdr_) return {};
    return zstd_decompress(names_blob_, hdr_->names_zstd_size,
                           static_cast<size_t>(hdr_->names_raw_size));
}

} // namespace genopack
