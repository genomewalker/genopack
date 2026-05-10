#include "genopack/derep_view.hpp"
#include "genopack/archive.hpp"

#define XXH_STATIC_LINKING_ONLY
#define XXH_IMPLEMENTATION
#include <xxhash.h>

#include <zstd.h>
#include <zlib.h>

#include <algorithm>
#include <cstring>
#include <fcntl.h>
#include <stdexcept>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

namespace genopack {

namespace {

constexpr uint32_t FILE_MAGIC = 0x46445047u;
constexpr uint32_t TAIL_MAGIC = 0x54445047u;
constexpr uint32_t HDR_MAGIC  = 0x48445047u;
constexpr uint32_t ASOF_MAGIC = 0x4F534147u;
constexpr uint32_t ARMP_MAGIC = 0x4D524147u;
constexpr uint32_t RTBL_MAGIC = 0x42545247u;
constexpr uint32_t G2RM_MAGIC = 0x4D523247u;
constexpr uint32_t EMBD_MAGIC = 0x424D4547u;
constexpr uint32_t TOC_MAGIC  = 0x434F5447u;

constexpr uint32_t GPD_SEC_HDR = 0x52444847u;
constexpr uint32_t GPD_SEC_AST = 0x54534147u;
constexpr uint32_t GPD_SEC_ASO = 0x4F534147u;
constexpr uint32_t GPD_SEC_ARM = 0x4D524147u;
constexpr uint32_t GPD_SEC_RTB = 0x42545247u;
constexpr uint32_t GPD_SEC_G2R = 0x4D523247u;
constexpr uint32_t GPD_SEC_EMB = 0x424D4547u;

constexpr uint32_t SENTINEL_UNCLUSTERED   = 0xFFFFFFFEu;
constexpr uint32_t SENTINEL_TOMBSTONE_OR_EMPTY = 0xFFFFFFFFu;

#pragma pack(push, 1)
struct GpdFileHeader {
    uint32_t magic;
    uint16_t format_major;
    uint16_t format_minor;
    uint64_t toc_offset;
    uint64_t toc_size;
    uint64_t reserved[5];
};
static_assert(sizeof(GpdFileHeader) == 64);

struct GpdTailLocator {
    uint64_t toc_offset;
    uint32_t magic;
    uint32_t crc32;
};
static_assert(sizeof(GpdTailLocator) == 16);

struct GpdSectionDesc {
    uint32_t type;
    uint32_t flags;
    uint64_t file_offset;
    uint64_t compressed_size;
    uint64_t uncompressed_size;
    uint64_t section_id;
    uint64_t reserved[2];
};
static_assert(sizeof(GpdSectionDesc) == 56);

struct GpdHeaderRaw {
    uint32_t magic;
    uint16_t format_major;
    uint16_t format_minor;
    uint64_t created_at_unix;
    uint8_t  run_id[16];
    uint16_t n_parts;
    uint16_t embedding_dim;
    uint8_t  embedding_dtype;
    uint8_t  has_cstats;
    uint8_t  pad0[2];
    uint64_t n_genomes;
    uint64_t n_reps;
    uint64_t n_unclustered;
};
static_assert(sizeof(GpdHeaderRaw) == 64);

struct GpdSourcePart {
    uint8_t  archive_uuid[16];
    uint64_t generation;
    uint64_t n_genomes_total;
    uint64_t n_genomes_live;
    uint64_t accession_set_hash;
};
static_assert(sizeof(GpdSourcePart) == 48);

struct GpdRepEntry {
    uint32_t rep_acc_ord;
    uint32_t cluster_size;
    uint64_t source_locator;
    uint16_t sketch_kmer;
    uint8_t  flags;
    uint8_t  pad;
    uint32_t cstat_offset;
};
static_assert(sizeof(GpdRepEntry) == 24);

struct GpdArmpEntry {
    uint64_t hash;
    uint32_t ordinal;
    uint32_t pad;
};
static_assert(sizeof(GpdArmpEntry) == 16);
#pragma pack(pop)

struct SectionDescView {
    uint32_t type;
    uint32_t flags;
    uint64_t file_offset;
    uint64_t compressed_size;
    uint64_t uncompressed_size;
};

float f16_to_f32(uint16_t h) {
    uint32_t sign = (h & 0x8000u) << 16;
    uint32_t exp  = (h & 0x7C00u) >> 10;
    uint32_t mant = (h & 0x03FFu);
    uint32_t bits;
    if (exp == 0) {
        if (mant == 0) {
            bits = sign;
        } else {
            // subnormal — normalise
            int e = -1;
            do { ++e; mant <<= 1; } while ((mant & 0x0400u) == 0);
            mant &= 0x03FFu;
            bits = sign | ((127u - 15u - static_cast<uint32_t>(e)) << 23) | (mant << 13);
        }
    } else if (exp == 31) {
        bits = sign | 0x7F800000u | (mant << 13);
    } else {
        bits = sign | ((exp + (127u - 15u)) << 23) | (mant << 13);
    }
    float f;
    std::memcpy(&f, &bits, 4);
    return f;
}

uint64_t accession_set_hash_for_part(const ArchiveReader& reader) {
    std::vector<std::string> accs;
    reader.scan_genome_accessions([&](std::string_view a, GenomeId) {
        accs.emplace_back(a);
    });
    std::sort(accs.begin(), accs.end());
    XXH3_state_t* st = XXH3_createState();
    XXH3_64bits_reset(st);
    for (size_t i = 0; i < accs.size(); ++i) {
        if (i > 0) { char nl = '\n'; XXH3_64bits_update(st, &nl, 1); }
        XXH3_64bits_update(st, accs[i].data(), accs[i].size());
    }
    uint64_t h = XXH3_64bits_digest(st);
    XXH3_freeState(st);
    return h;
}

bool read_archive_uuid(const std::filesystem::path& p, uint8_t out[16]) {
    auto try_read = [&](const std::filesystem::path& f) {
        int fd = ::open(f.c_str(), O_RDONLY);
        if (fd < 0) return false;
        uint8_t buf[24] = {};
        bool ok = (::read(fd, buf, sizeof(buf)) == static_cast<ssize_t>(sizeof(buf)));
        ::close(fd);
        if (!ok) return false;
        std::memcpy(out,     buf + 8,  8);
        std::memcpy(out + 8, buf + 16, 8);
        return true;
    };
    if (try_read(p)) return true;
    return try_read(p / "toc.bin");
}

} // namespace

struct DerepView::Impl {
    int                  fd      = -1;
    const uint8_t*       map     = nullptr;
    size_t               map_len = 0;

    GpdHeaderRaw                       hdr_raw{};
    std::vector<GpdSourcePart>         parts;
    std::string                        geodesic_version;
    uint8_t                            kmer_sizes[7]{};
    uint8_t                            n_kmer_sizes = 0;
    uint32_t                           sketch_size  = 0;
    uint64_t                           sig1_seed = 0, sig2_seed = 0;
    float                              jaccard_thresh = 0.f;

    std::vector<uint8_t>               astr;
    std::vector<uint32_t>              asof;
    std::vector<GpdArmpEntry>          armp_entries;
    uint32_t                           armp_n_buckets = 0;
    bool                               has_armp = false;

    std::vector<GpdRepEntry>           rtbl;
    std::unordered_map<uint32_t,uint32_t> ord_to_repid;

    std::vector<uint32_t>              g2rm;

    bool                               has_embd = false;
    uint16_t                           embd_dim = 0;
    uint8_t                            embd_dtype = 0;
    std::vector<uint8_t>               embd;

    void unmap() {
        if (map) { ::munmap(const_cast<uint8_t*>(map), map_len); map = nullptr; map_len = 0; }
        if (fd >= 0) { ::close(fd); fd = -1; }
    }

    std::vector<uint8_t> load_section(const SectionDescView& s) const {
        if (s.file_offset + s.compressed_size > map_len)
            throw std::runtime_error("gpd: section overruns file");
        std::vector<uint8_t> out;
        if (s.flags & 1u) {
            out.resize(s.uncompressed_size);
            size_t n = ZSTD_decompress(out.data(), out.size(),
                                       map + s.file_offset, s.compressed_size);
            if (ZSTD_isError(n))
                throw std::runtime_error(std::string("gpd: zstd decompress: ") + ZSTD_getErrorName(n));
            if (n != s.uncompressed_size)
                throw std::runtime_error("gpd: zstd: size mismatch");
        } else {
            out.assign(map + s.file_offset,
                       map + s.file_offset + s.compressed_size);
        }
        return out;
    }

    std::string_view acc_at(uint32_t ord) const {
        uint32_t a = asof[ord];
        uint32_t b = asof[ord + 1];
        return std::string_view(reinterpret_cast<const char*>(astr.data() + a), b - a);
    }

    int32_t lookup_ord(std::string_view acc) const {
        if (has_armp && armp_n_buckets > 0) {
            XXH128_hash_t h128 = XXH3_128bits(acc.data(), acc.size());
            uint64_t hkey  = h128.high64;
            uint32_t mask  = armp_n_buckets - 1;
            uint32_t b     = static_cast<uint32_t>(hkey & mask);
            for (uint32_t probes = 0; probes < armp_n_buckets; ++probes) {
                const auto& e = armp_entries[b];
                if (e.ordinal == SENTINEL_TOMBSTONE_OR_EMPTY) return -1;
                if (e.hash == hkey && acc_at(e.ordinal) == acc)
                    return static_cast<int32_t>(e.ordinal);
                b = (b + 1) & mask;
            }
            return -1;
        }
        // ASOF binary search fallback
        uint32_t lo = 0, hi = static_cast<uint32_t>(asof.size() ? asof.size() - 1 : 0);
        while (lo < hi) {
            uint32_t mid = lo + (hi - lo) / 2;
            std::string_view m = acc_at(mid);
            int c = m.compare(acc);
            if (c == 0)      return static_cast<int32_t>(mid);
            else if (c < 0)  lo = mid + 1;
            else             hi = mid;
        }
        return -1;
    }
};

DerepView::DerepView() : impl_(std::make_unique<Impl>()) {}
DerepView::~DerepView() { if (impl_) impl_->unmap(); }
DerepView::DerepView(DerepView&&) noexcept = default;
DerepView& DerepView::operator=(DerepView&&) noexcept = default;

void DerepView::close() { if (impl_) impl_->unmap(); impl_ = std::make_unique<Impl>(); }
bool DerepView::is_open() const { return impl_ && impl_->map != nullptr; }

void DerepView::open(const std::filesystem::path& gpd_path) {
    close();
    auto& I = *impl_;

    I.fd = ::open(gpd_path.c_str(), O_RDONLY);
    if (I.fd < 0)
        throw std::runtime_error("gpd: cannot open: " + gpd_path.string());

    struct stat st{};
    if (::fstat(I.fd, &st) != 0)
        throw std::runtime_error("gpd: fstat failed");
    I.map_len = static_cast<size_t>(st.st_size);
    if (I.map_len < sizeof(GpdFileHeader) + sizeof(GpdTailLocator))
        throw std::runtime_error("gpd: file too short");

    void* m = ::mmap(nullptr, I.map_len, PROT_READ, MAP_PRIVATE, I.fd, 0);
    if (m == MAP_FAILED)
        throw std::runtime_error("gpd: mmap failed");
    I.map = static_cast<const uint8_t*>(m);

    GpdFileHeader fh{};
    std::memcpy(&fh, I.map, sizeof(fh));
    if (fh.magic != FILE_MAGIC)
        throw std::runtime_error("gpd: bad file magic");
    if (fh.format_major != 1)
        throw std::runtime_error("gpd: unsupported format_major");

    GpdTailLocator tail{};
    std::memcpy(&tail, I.map + I.map_len - sizeof(tail), sizeof(tail));
    if (tail.magic != TAIL_MAGIC)
        throw std::runtime_error("gpd: bad tail magic");
    if (tail.toc_offset != fh.toc_offset)
        throw std::runtime_error("gpd: toc_offset mismatch");

    if (fh.toc_offset + fh.toc_size > I.map_len)
        throw std::runtime_error("gpd: toc overruns file");

    const uint8_t* tocp = I.map + fh.toc_offset;
    uint32_t toc_magic, n_sections, crc_stored, pad;
    std::memcpy(&toc_magic,  tocp + 0,  4);
    std::memcpy(&n_sections, tocp + 4,  4);
    std::memcpy(&crc_stored, tocp + 8,  4);
    std::memcpy(&pad,        tocp + 12, 4);
    if (toc_magic != TOC_MAGIC)
        throw std::runtime_error("gpd: bad TOC magic");

    const uint8_t* descp = tocp + 16;
    size_t descs_bytes  = size_t{n_sections} * sizeof(GpdSectionDesc);
    if (16 + descs_bytes > fh.toc_size)
        throw std::runtime_error("gpd: TOC undersized");

    uLong crc = ::crc32(0L, reinterpret_cast<const Bytef*>(descp),
                        static_cast<uInt>(descs_bytes));
    if (static_cast<uint32_t>(crc) != crc_stored)
        throw std::runtime_error("gpd: TOC crc mismatch");
    if (tail.crc32 != crc_stored)
        throw std::runtime_error("gpd: tail crc mismatch");

    std::vector<SectionDescView> secs(n_sections);
    for (uint32_t i = 0; i < n_sections; ++i) {
        GpdSectionDesc d{};
        std::memcpy(&d, descp + i * sizeof(GpdSectionDesc), sizeof(d));
        secs[i] = SectionDescView{d.type, d.flags, d.file_offset,
                                  d.compressed_size, d.uncompressed_size};
    }

    auto find_sec = [&](uint32_t type) -> const SectionDescView* {
        for (auto& s : secs) if (s.type == type) return &s;
        return nullptr;
    };

    const SectionDescView* hdr_s = find_sec(GPD_SEC_HDR);
    if (!hdr_s) throw std::runtime_error("gpd: missing HDR");
    auto hdr_buf = I.load_section(*hdr_s);
    if (hdr_buf.size() < sizeof(GpdHeaderRaw))
        throw std::runtime_error("gpd: HDR too small");
    std::memcpy(&I.hdr_raw, hdr_buf.data(), sizeof(GpdHeaderRaw));
    if (I.hdr_raw.magic != HDR_MAGIC)
        throw std::runtime_error("gpd: bad HDR magic");

    size_t off = sizeof(GpdHeaderRaw);
    I.parts.resize(I.hdr_raw.n_parts);
    if (off + I.hdr_raw.n_parts * sizeof(GpdSourcePart) > hdr_buf.size())
        throw std::runtime_error("gpd: HDR truncated (parts)");
    std::memcpy(I.parts.data(), hdr_buf.data() + off,
                I.hdr_raw.n_parts * sizeof(GpdSourcePart));
    off += I.hdr_raw.n_parts * sizeof(GpdSourcePart);

    if (off + 36 <= hdr_buf.size()) {
        I.n_kmer_sizes = hdr_buf[off];
        std::memcpy(I.kmer_sizes, hdr_buf.data() + off + 1, 7);
        std::memcpy(&I.sketch_size,    hdr_buf.data() + off + 8,  4);
        std::memcpy(&I.sig1_seed,      hdr_buf.data() + off + 12, 8);
        std::memcpy(&I.sig2_seed,      hdr_buf.data() + off + 20, 8);
        std::memcpy(&I.jaccard_thresh, hdr_buf.data() + off + 28, 4);
        uint16_t ver_len;
        std::memcpy(&ver_len, hdr_buf.data() + off + 32, 2);
        off += 36;
        if (off + ver_len <= hdr_buf.size())
            I.geodesic_version.assign(reinterpret_cast<const char*>(hdr_buf.data() + off), ver_len);
    }

    const SectionDescView* astr_s = find_sec(GPD_SEC_AST);
    const SectionDescView* asof_s = find_sec(GPD_SEC_ASO);
    const SectionDescView* g2rm_s = find_sec(GPD_SEC_G2R);
    const SectionDescView* rtbl_s = find_sec(GPD_SEC_RTB);
    if (!astr_s || !asof_s || !g2rm_s || !rtbl_s)
        throw std::runtime_error("gpd: required section missing");

    I.astr = I.load_section(*astr_s);
    auto asof_buf = I.load_section(*asof_s);
    if (asof_buf.size() < 16)
        throw std::runtime_error("gpd: ASOF too small");
    uint32_t asof_magic, asof_n;
    std::memcpy(&asof_magic, asof_buf.data() + 0, 4);
    std::memcpy(&asof_n,     asof_buf.data() + 4, 4);
    if (asof_magic != ASOF_MAGIC)
        throw std::runtime_error("gpd: bad ASOF magic");
    I.asof.resize(size_t{asof_n} + 1);
    if (asof_buf.size() < 16 + (size_t{asof_n} + 1) * sizeof(uint32_t))
        throw std::runtime_error("gpd: ASOF truncated");
    std::memcpy(I.asof.data(), asof_buf.data() + 16,
                (size_t{asof_n} + 1) * sizeof(uint32_t));

    auto rtbl_buf = I.load_section(*rtbl_s);
    if (rtbl_buf.size() < 16)
        throw std::runtime_error("gpd: RTBL too small");
    uint32_t rtbl_magic, rtbl_n;
    std::memcpy(&rtbl_magic, rtbl_buf.data() + 0, 4);
    std::memcpy(&rtbl_n,     rtbl_buf.data() + 4, 4);
    if (rtbl_magic != RTBL_MAGIC)
        throw std::runtime_error("gpd: bad RTBL magic");
    I.rtbl.resize(rtbl_n);
    if (rtbl_buf.size() < 16 + size_t{rtbl_n} * sizeof(GpdRepEntry))
        throw std::runtime_error("gpd: RTBL truncated");
    std::memcpy(I.rtbl.data(), rtbl_buf.data() + 16,
                size_t{rtbl_n} * sizeof(GpdRepEntry));
    I.ord_to_repid.reserve(rtbl_n * 2);
    for (uint32_t i = 0; i < rtbl_n; ++i)
        I.ord_to_repid[I.rtbl[i].rep_acc_ord] = i;

    auto g2rm_buf = I.load_section(*g2rm_s);
    if (g2rm_buf.size() < 16)
        throw std::runtime_error("gpd: G2RM too small");
    uint32_t g2rm_magic, g2rm_n;
    std::memcpy(&g2rm_magic, g2rm_buf.data() + 0, 4);
    std::memcpy(&g2rm_n,     g2rm_buf.data() + 4, 4);
    if (g2rm_magic != G2RM_MAGIC)
        throw std::runtime_error("gpd: bad G2RM magic");
    I.g2rm.resize(g2rm_n);
    if (g2rm_buf.size() < 16 + size_t{g2rm_n} * sizeof(uint32_t))
        throw std::runtime_error("gpd: G2RM truncated");
    std::memcpy(I.g2rm.data(), g2rm_buf.data() + 16,
                size_t{g2rm_n} * sizeof(uint32_t));

    if (const SectionDescView* armp_s = find_sec(GPD_SEC_ARM)) {
        auto buf = I.load_section(*armp_s);
        if (buf.size() < 16)
            throw std::runtime_error("gpd: ARMP too small");
        uint32_t magic, nb;
        std::memcpy(&magic, buf.data() + 0, 4);
        std::memcpy(&nb,    buf.data() + 4, 4);
        if (magic != ARMP_MAGIC)
            throw std::runtime_error("gpd: bad ARMP magic");
        I.armp_n_buckets = nb;
        I.armp_entries.resize(nb);
        if (buf.size() < 16 + size_t{nb} * sizeof(GpdArmpEntry))
            throw std::runtime_error("gpd: ARMP truncated");
        std::memcpy(I.armp_entries.data(), buf.data() + 16,
                    size_t{nb} * sizeof(GpdArmpEntry));
        I.has_armp = true;
    }

    if (const SectionDescView* embd_s = find_sec(GPD_SEC_EMB)) {
        auto buf = I.load_section(*embd_s);
        // EMBD header per DEREP_FORMAT.md is 16 bytes:
        //   uint32 magic; uint16 dim; uint8 dtype; uint8 pad0;
        //   uint32 n_reps; uint32 pad1;
        constexpr size_t EMBD_HDR_BYTES = 16;
        if (buf.size() < EMBD_HDR_BYTES)
            throw std::runtime_error("gpd: EMBD too small");
        uint32_t magic;
        uint16_t dim;
        uint8_t  dtype;
        uint32_t n_reps;
        std::memcpy(&magic,  buf.data() + 0, 4);
        std::memcpy(&dim,    buf.data() + 4, 2);
        dtype = buf[6];
        std::memcpy(&n_reps, buf.data() + 8, 4);
        if (magic != EMBD_MAGIC)
            throw std::runtime_error("gpd: bad EMBD magic");
        size_t bytes_per = size_t{dim} * (dtype == 1 ? 2u : 4u);
        size_t need = EMBD_HDR_BYTES + size_t{n_reps} * bytes_per;
        if (buf.size() < need)
            throw std::runtime_error("gpd: EMBD truncated");
        I.embd_dim   = dim;
        I.embd_dtype = dtype;
        I.embd.assign(buf.begin() + EMBD_HDR_BYTES, buf.begin() + need);
        I.has_embd = true;
    }
}

DerepStaleness DerepView::check(const ArchiveSetReader& archive) const {
    const auto& I = *impl_;
    const auto& paths = archive.part_paths();

    if (paths.size() != I.parts.size())
        return DerepStaleness::Mismatch;

    DerepStaleness level = DerepStaleness::Valid;
    for (size_t i = 0; i < paths.size(); ++i) {
        ArchiveReader rd;
        rd.open(paths[i]);
        auto stats = rd.archive_stats();
        uint8_t live_uuid[16] = {};
        read_archive_uuid(paths[i], live_uuid);
        uint64_t live_hash = accession_set_hash_for_part(rd);

        const auto& fp = I.parts[i];
        if (live_hash == fp.accession_set_hash) {
            bool same_layout = (std::memcmp(live_uuid, fp.archive_uuid, 16) == 0)
                            && (stats.generation == fp.generation);
            if (!same_layout && level == DerepStaleness::Valid)
                level = DerepStaleness::LayoutChangedSameLiveSet;
        } else {
            if (stats.n_genomes_live > fp.n_genomes_live) {
                if (level < DerepStaleness::StaleNewGenomes)
                    level = DerepStaleness::StaleNewGenomes;
            } else if (stats.n_genomes_live < fp.n_genomes_live) {
                if (level < DerepStaleness::StaleTombstones)
                    level = DerepStaleness::StaleTombstones;
            } else {
                level = DerepStaleness::Mismatch;
            }
        }
    }
    return level;
}

RepStatus DerepView::status_for_accession(std::string_view accession) const {
    const auto& I = *impl_;
    int32_t ord = I.lookup_ord(accession);
    RepStatus rs;
    if (ord < 0) { rs.kind = RepStatus::Kind::Absent; return rs; }
    uint32_t rep = I.g2rm[ord];
    if (rep == SENTINEL_UNCLUSTERED) {
        rs.kind = RepStatus::Kind::Unclustered;
        return rs;
    }
    if (rep == SENTINEL_TOMBSTONE_OR_EMPTY) {
        rs.kind = RepStatus::Kind::Tombstoned;
        return rs;
    }
    if (rep >= I.rtbl.size()) {
        rs.kind = RepStatus::Kind::Absent;
        return rs;
    }
    const auto& re = I.rtbl[rep];
    rs.rep_id        = rep;
    rs.rep_accession = I.acc_at(re.rep_acc_ord);
    rs.kind = (re.rep_acc_ord == static_cast<uint32_t>(ord))
              ? RepStatus::Kind::Representative
              : RepStatus::Kind::Member;
    return rs;
}

void DerepView::scan_representatives(
    const std::function<void(uint32_t, std::string_view, uint32_t)>& cb) const
{
    const auto& I = *impl_;
    for (uint32_t r = 0; r < I.rtbl.size(); ++r) {
        const auto& e = I.rtbl[r];
        cb(r, I.acc_at(e.rep_acc_ord), e.cluster_size);
    }
}

bool DerepView::embedding_for_rep(uint32_t rep_id, std::span<float> out) const {
    const auto& I = *impl_;
    if (!I.has_embd) return false;
    if (rep_id >= I.rtbl.size()) return false;
    if (out.size() < I.embd_dim) return false;
    size_t bytes_per = size_t{I.embd_dim} * (I.embd_dtype == 1 ? 2u : 4u);
    const uint8_t* row = I.embd.data() + size_t{rep_id} * bytes_per;
    if (I.embd_dtype == 0) {
        std::memcpy(out.data(), row, bytes_per);
    } else {
        for (uint16_t i = 0; i < I.embd_dim; ++i) {
            uint16_t h;
            std::memcpy(&h, row + i * 2, 2);
            out[i] = f16_to_f32(h);
        }
    }
    return true;
}

DerepStats DerepView::stats() const {
    const auto& I = *impl_;
    DerepStats s;
    s.n_genomes_indexed = I.hdr_raw.n_genomes;
    s.n_reps            = static_cast<uint32_t>(I.rtbl.size());
    s.n_unclustered     = static_cast<uint32_t>(I.hdr_raw.n_unclustered);
    s.n_singletons      = 0;
    uint64_t members = 0;
    for (const auto& e : I.rtbl) {
        if (e.cluster_size == 1) ++s.n_singletons;
        if (e.cluster_size > 1)  members += (e.cluster_size - 1);
    }
    s.n_members = members;
    return s;
}

uint32_t DerepView::embedding_dim() const {
    return impl_->has_embd ? impl_->embd_dim : 0u;
}

uint64_t DerepView::accession_set_hash() const {
    XXH3_state_t* st = XXH3_createState();
    XXH3_64bits_reset(st);
    for (const auto& p : impl_->parts)
        XXH3_64bits_update(st, &p.accession_set_hash, sizeof(p.accession_set_hash));
    uint64_t h = XXH3_64bits_digest(st);
    XXH3_freeState(st);
    return h;
}

uint32_t DerepView::format_major() const { return impl_->hdr_raw.format_major; }
uint32_t DerepView::format_minor() const { return impl_->hdr_raw.format_minor; }

std::array<uint8_t,16> DerepView::run_id() const {
    std::array<uint8_t,16> out;
    std::memcpy(out.data(), impl_->hdr_raw.run_id, 16);
    return out;
}
uint64_t DerepView::created_at_unix() const { return impl_->hdr_raw.created_at_unix; }
uint16_t DerepView::source_n_parts() const { return impl_->hdr_raw.n_parts; }

} // namespace genopack
