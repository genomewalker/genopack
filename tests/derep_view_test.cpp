#include "genopack/derep_view.hpp"

#define XXH_STATIC_LINKING_ONLY
#include <xxhash.h>

#include <zstd.h>
#include <zlib.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

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

constexpr uint32_t SENT_UNCL  = 0xFFFFFFFEu;
constexpr uint32_t SENT_TOMB  = 0xFFFFFFFFu;

#pragma pack(push, 1)
struct FileHeader {
    uint32_t magic;
    uint16_t fmaj;
    uint16_t fmin;
    uint64_t toc_offset;
    uint64_t toc_size;
    uint64_t reserved[5];
};
struct Tail {
    uint64_t toc_offset;
    uint32_t magic;
    uint32_t crc32;
};
struct SecDesc {
    uint32_t type;
    uint32_t flags;
    uint64_t file_offset;
    uint64_t comp;
    uint64_t uncomp;
    uint64_t section_id;
    uint64_t reserved[2];
};
struct HdrRaw {
    uint32_t magic;
    uint16_t fmaj;
    uint16_t fmin;
    uint64_t created;
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
struct DerepParams {
    uint8_t  n_kmer_sizes;
    uint8_t  kmer_sizes[7];
    uint32_t sketch_size;
    uint64_t sig1_seed;
    uint64_t sig2_seed;
    float    jaccard_thresh;
    uint16_t ver_len;
    uint8_t  pad1[2];
};
struct RepEntry {
    uint32_t rep_acc_ord;
    uint32_t cluster_size;
    uint64_t source_locator;
    uint16_t sketch_kmer;
    uint8_t  flags;
    uint8_t  pad;
    uint32_t cstat_offset;
};
struct ArmpEntry {
    uint64_t hash;
    uint32_t ordinal;
    uint32_t pad;
};
#pragma pack(pop)

uint64_t align8(uint64_t x) { return (x + 7u) & ~uint64_t{7}; }

std::vector<uint8_t> zstd_pack(const void* d, size_t n) {
    std::vector<uint8_t> out(ZSTD_compressBound(n));
    size_t c = ZSTD_compress(out.data(), out.size(), d, n, 3);
    if (ZSTD_isError(c)) throw std::runtime_error("zstd");
    out.resize(c);
    return out;
}

uint16_t f32_to_f16(float f) {
    uint32_t b;
    std::memcpy(&b, &f, 4);
    uint32_t sign = (b >> 16) & 0x8000u;
    int32_t  e    = static_cast<int32_t>((b >> 23) & 0xFF) - 127 + 15;
    uint32_t m    = b & 0x7FFFFFu;
    if (e <= 0)  return static_cast<uint16_t>(sign);
    if (e >= 31) return static_cast<uint16_t>(sign | 0x7C00u);
    return static_cast<uint16_t>(sign | (uint32_t(e) << 10) | (m >> 13));
}

struct GpdBuffer { std::vector<uint8_t> bytes; };

struct BuildOpts {
    uint32_t n_genomes   = 32;
    uint32_t n_reps      = 8;
    uint32_t n_unclus    = 2;
    uint32_t n_tombstone = 2;
    uint16_t embd_dim    = 0;     // 0 → no EMBD
    bool     emit_armp   = true;
    bool     corrupt_section_crc = false;
};

struct GenomeRec {
    std::string acc;
    enum { Rep, Member, Unclus, Tomb } kind;
    uint32_t rep_idx;
};

GpdBuffer build_gpd(const BuildOpts& o) {
    std::vector<GenomeRec> recs;
    recs.reserve(o.n_genomes);
    auto pad_acc = [](uint32_t i) {
        char buf[16];
        std::snprintf(buf, sizeof(buf), "ACC_%05u", i);
        return std::string(buf);
    };
    for (uint32_t i = 0; i < o.n_reps; ++i)
        recs.push_back({pad_acc(i), GenomeRec::Rep, i});
    uint32_t next = o.n_reps;
    uint32_t n_members = o.n_genomes - o.n_reps - o.n_unclus - o.n_tombstone;
    for (uint32_t m = 0; m < n_members; ++m, ++next)
        recs.push_back({pad_acc(next), GenomeRec::Member, m % o.n_reps});
    for (uint32_t u = 0; u < o.n_unclus; ++u, ++next)
        recs.push_back({pad_acc(next), GenomeRec::Unclus, 0});
    for (uint32_t t = 0; t < o.n_tombstone; ++t, ++next)
        recs.push_back({pad_acc(next), GenomeRec::Tomb, 0});

    std::vector<uint32_t> order(recs.size());
    for (uint32_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](uint32_t a, uint32_t b) { return recs[a].acc < recs[b].acc; });
    std::vector<uint32_t> ord_of(recs.size());
    for (uint32_t o2 = 0; o2 < order.size(); ++o2) ord_of[order[o2]] = o2;

    std::string astr;
    std::vector<uint32_t> offs(recs.size() + 1);
    for (uint32_t i = 0; i < order.size(); ++i) {
        offs[i] = static_cast<uint32_t>(astr.size());
        astr += recs[order[i]].acc;
    }
    offs.back() = static_cast<uint32_t>(astr.size());

    std::vector<uint8_t> asof_raw;
    asof_raw.resize(16 + (recs.size() + 1) * sizeof(uint32_t));
    uint32_t am = ASOF_MAGIC, an = static_cast<uint32_t>(recs.size());
    std::memcpy(asof_raw.data() + 0, &am, 4);
    std::memcpy(asof_raw.data() + 4, &an, 4);
    std::memcpy(asof_raw.data() + 16, offs.data(), (recs.size() + 1) * sizeof(uint32_t));

    struct RepRow { uint32_t ord; uint32_t orig; };
    std::vector<RepRow> reps;
    for (uint32_t i = 0; i < recs.size(); ++i)
        if (recs[i].kind == GenomeRec::Rep)
            reps.push_back({ord_of[i], i});
    std::sort(reps.begin(), reps.end(),
              [](const RepRow& a, const RepRow& b) { return a.ord < b.ord; });

    std::vector<uint32_t> orig_to_repid(recs.size(), UINT32_MAX);
    for (uint32_t r = 0; r < reps.size(); ++r)
        orig_to_repid[reps[r].orig] = r;

    std::vector<uint32_t> g2rm(recs.size(), 0);
    for (uint32_t i = 0; i < recs.size(); ++i) {
        uint32_t ord = ord_of[i];
        switch (recs[i].kind) {
        case GenomeRec::Rep:    g2rm[ord] = orig_to_repid[i]; break;
        case GenomeRec::Member: g2rm[ord] = orig_to_repid[recs[i].rep_idx]; break;
        case GenomeRec::Unclus: g2rm[ord] = SENT_UNCL; break;
        case GenomeRec::Tomb:   g2rm[ord] = SENT_TOMB; break;
        }
    }

    std::vector<uint8_t> g2rm_raw;
    g2rm_raw.resize(16 + recs.size() * 4);
    uint32_t gm = G2RM_MAGIC, gn = static_cast<uint32_t>(recs.size());
    std::memcpy(g2rm_raw.data() + 0, &gm, 4);
    std::memcpy(g2rm_raw.data() + 4, &gn, 4);
    std::memcpy(g2rm_raw.data() + 16, g2rm.data(), recs.size() * 4);

    std::vector<RepEntry> rtbl_entries(reps.size());
    for (uint32_t r = 0; r < reps.size(); ++r) {
        rtbl_entries[r].rep_acc_ord   = reps[r].ord;
        rtbl_entries[r].cluster_size  = (r == 0) ? 1u : 2u;
        rtbl_entries[r].source_locator= 0;
        rtbl_entries[r].sketch_kmer   = 16;
        rtbl_entries[r].flags         = 0x01;
        rtbl_entries[r].pad           = 0;
        rtbl_entries[r].cstat_offset  = UINT32_MAX;
    }
    std::vector<uint8_t> rtbl_raw;
    rtbl_raw.resize(16 + reps.size() * sizeof(RepEntry));
    uint32_t rm = RTBL_MAGIC, rn = static_cast<uint32_t>(reps.size());
    std::memcpy(rtbl_raw.data() + 0, &rm, 4);
    std::memcpy(rtbl_raw.data() + 4, &rn, 4);
    std::memcpy(rtbl_raw.data() + 16, rtbl_entries.data(), reps.size() * sizeof(RepEntry));

    std::vector<uint8_t> armp_raw;
    if (o.emit_armp) {
        uint32_t nb = 1;
        while (nb < recs.size() * 2) nb <<= 1;
        std::vector<ArmpEntry> entries(nb);
        for (auto& e : entries) { e.hash = 0; e.ordinal = SENT_TOMB; e.pad = 0; }
        for (uint32_t o2 = 0; o2 < recs.size(); ++o2) {
            const auto& acc = recs[order[o2]].acc;
            XXH128_hash_t h = XXH3_128bits(acc.data(), acc.size());
            uint32_t b = static_cast<uint32_t>(h.high64 & (nb - 1));
            while (entries[b].ordinal != SENT_TOMB) b = (b + 1) & (nb - 1);
            entries[b].hash    = h.high64;
            entries[b].ordinal = o2;
        }
        armp_raw.resize(16 + nb * sizeof(ArmpEntry));
        uint32_t mg = ARMP_MAGIC, nb_w = nb, seed = 0, pad = 0;
        std::memcpy(armp_raw.data() + 0,  &mg,  4);
        std::memcpy(armp_raw.data() + 4,  &nb_w,4);
        std::memcpy(armp_raw.data() + 8,  &seed,4);
        std::memcpy(armp_raw.data() + 12, &pad, 4);
        std::memcpy(armp_raw.data() + 16, entries.data(), nb * sizeof(ArmpEntry));
    }

    HdrRaw hr{};
    hr.magic = HDR_MAGIC;
    hr.fmaj = 1; hr.fmin = 0;
    hr.created = 1700000000;
    for (int i = 0; i < 16; ++i) hr.run_id[i] = static_cast<uint8_t>(i);
    hr.n_parts        = 1;
    hr.embedding_dim  = o.embd_dim;
    hr.embedding_dtype = 1;
    hr.has_cstats     = 0;
    hr.n_genomes      = recs.size();
    hr.n_reps         = reps.size();
    hr.n_unclustered  = o.n_unclus;

    DerepParams dp{};
    dp.n_kmer_sizes = 1;
    dp.kmer_sizes[0] = 16;
    dp.sketch_size = 1024;
    dp.sig1_seed = 1; dp.sig2_seed = 2;
    dp.jaccard_thresh = 0.95f;
    dp.ver_len = 0;

    std::vector<uint8_t> hdr_raw;
    hdr_raw.resize(sizeof(HdrRaw) + sizeof(uint8_t) * 48 + sizeof(DerepParams));
    size_t pos = 0;
    std::memcpy(hdr_raw.data() + pos, &hr, sizeof(HdrRaw)); pos += sizeof(HdrRaw);
    {
        // n_parts=1 GpdSourcePart record (zeroed; check() not exercised by basic test)
        std::memset(hdr_raw.data() + pos, 0, 48); pos += 48;
    }
    std::memcpy(hdr_raw.data() + pos, &dp, sizeof(DerepParams)); pos += sizeof(DerepParams);
    hdr_raw.resize(pos);

    std::vector<uint8_t> embd_raw;
    if (o.embd_dim > 0) {
        // EMBD header per spec is 16 bytes (magic4 + dim2 + dtype1 + pad0_1
        // + n_reps4 + pad1_4).
        constexpr size_t EMBD_HDR = 16;
        size_t bytes_per = size_t{o.embd_dim} * 2u;
        embd_raw.resize(EMBD_HDR + reps.size() * bytes_per);
        uint32_t mg = EMBD_MAGIC;
        uint16_t dim = o.embd_dim;
        uint8_t  dt = 1, p0 = 0;
        uint32_t nrp = static_cast<uint32_t>(reps.size()), pad1 = 0;
        std::memcpy(embd_raw.data() + 0, &mg,  4);
        std::memcpy(embd_raw.data() + 4, &dim, 2);
        embd_raw[6] = dt;
        embd_raw[7] = p0;
        std::memcpy(embd_raw.data() + 8,  &nrp,  4);
        std::memcpy(embd_raw.data() + 12, &pad1, 4);
        for (uint32_t r = 0; r < reps.size(); ++r) {
            for (uint16_t k = 0; k < o.embd_dim; ++k) {
                float v = static_cast<float>(r) + 0.25f * static_cast<float>(k);
                uint16_t h = f32_to_f16(v);
                std::memcpy(embd_raw.data() + EMBD_HDR + r * bytes_per + k * 2,
                            &h, 2);
            }
        }
    }

    std::vector<uint8_t> astr_c = zstd_pack(astr.data(), astr.size());
    std::vector<uint8_t> asof_c = zstd_pack(asof_raw.data(), asof_raw.size());
    std::vector<uint8_t> g2rm_c = zstd_pack(g2rm_raw.data(), g2rm_raw.size());
    std::vector<uint8_t> armp_c;
    if (!armp_raw.empty()) armp_c = zstd_pack(armp_raw.data(), armp_raw.size());

    GpdBuffer gb;
    auto& out = gb.bytes;
    out.resize(sizeof(FileHeader));

    struct Sec { uint32_t type, flags; uint64_t off, csz, usz; };
    std::vector<Sec> secs;
    auto place = [&](uint32_t type, uint32_t flags,
                     const std::vector<uint8_t>& data, uint64_t usz) {
        uint64_t off = out.size();
        out.insert(out.end(), data.begin(), data.end());
        while (out.size() % 8) out.push_back(0);
        secs.push_back({type, flags, off, data.size(), usz});
    };
    place(GPD_SEC_HDR, 0, hdr_raw, hdr_raw.size());
    place(GPD_SEC_AST, 1, astr_c,  astr.size());
    place(GPD_SEC_ASO, 1, asof_c,  asof_raw.size());
    if (!armp_c.empty()) place(GPD_SEC_ARM, 1, armp_c, armp_raw.size());
    place(GPD_SEC_RTB, 0, rtbl_raw, rtbl_raw.size());
    place(GPD_SEC_G2R, 1, g2rm_c,   g2rm_raw.size());
    if (!embd_raw.empty()) place(GPD_SEC_EMB, 0, embd_raw, embd_raw.size());

    uint64_t toc_off = out.size();
    std::vector<SecDesc> descs(secs.size());
    for (size_t i = 0; i < secs.size(); ++i) {
        descs[i].type = secs[i].type;
        descs[i].flags = secs[i].flags;
        descs[i].file_offset = secs[i].off;
        descs[i].comp = secs[i].csz;
        descs[i].uncomp = secs[i].usz;
        descs[i].section_id = i + 1;
        descs[i].reserved[0] = descs[i].reserved[1] = 0;
    }
    uint32_t n_sec = static_cast<uint32_t>(descs.size());
    uLong crc = ::crc32(0L, reinterpret_cast<const Bytef*>(descs.data()),
                        static_cast<uInt>(descs.size() * sizeof(SecDesc)));
    if (o.corrupt_section_crc) crc ^= 0xDEADBEEFu;

    uint32_t toc_magic = TOC_MAGIC;
    uint32_t pad = 0;
    uint32_t crc32_field = static_cast<uint32_t>(crc);
    out.insert(out.end(), reinterpret_cast<uint8_t*>(&toc_magic),
                          reinterpret_cast<uint8_t*>(&toc_magic) + 4);
    out.insert(out.end(), reinterpret_cast<uint8_t*>(&n_sec),
                          reinterpret_cast<uint8_t*>(&n_sec) + 4);
    out.insert(out.end(), reinterpret_cast<uint8_t*>(&crc32_field),
                          reinterpret_cast<uint8_t*>(&crc32_field) + 4);
    out.insert(out.end(), reinterpret_cast<uint8_t*>(&pad),
                          reinterpret_cast<uint8_t*>(&pad) + 4);
    out.insert(out.end(),
               reinterpret_cast<uint8_t*>(descs.data()),
               reinterpret_cast<uint8_t*>(descs.data()) + descs.size() * sizeof(SecDesc));
    uint64_t toc_size = 16 + descs.size() * sizeof(SecDesc);
    while (out.size() % 8) out.push_back(0);

    Tail tail{};
    tail.toc_offset = toc_off;
    tail.magic      = TAIL_MAGIC;
    tail.crc32      = static_cast<uint32_t>(crc);
    out.insert(out.end(), reinterpret_cast<uint8_t*>(&tail),
                          reinterpret_cast<uint8_t*>(&tail) + sizeof(Tail));

    FileHeader fh{};
    fh.magic = FILE_MAGIC;
    fh.fmaj = 1; fh.fmin = 0;
    fh.toc_offset = toc_off;
    fh.toc_size = toc_size;
    std::memcpy(out.data(), &fh, sizeof(fh));
    return gb;
}

std::string write_temp(const GpdBuffer& g) {
    char tmpl[] = "/tmp/gpd_view_test_XXXXXX";
    int fd = ::mkstemp(tmpl);
    if (fd < 0) throw std::runtime_error("mkstemp");
    ssize_t w = ::write(fd, g.bytes.data(), g.bytes.size());
    (void)w;
    ::close(fd);
    return tmpl;
}

} // namespace

int main() {
    using namespace genopack;

    {
        BuildOpts o;
        auto buf = build_gpd(o);
        auto path = write_temp(buf);

        DerepView v;
        v.open(path);
        auto s = v.stats();
        if (s.n_genomes_indexed != 32) { std::cerr << "n_genomes\n"; return 1; }
        if (s.n_reps != 8)              { std::cerr << "n_reps\n"; return 1; }
        if (s.n_unclustered != 2)       { std::cerr << "n_unclus\n"; return 1; }
        if (v.embedding_dim() != 0)     { std::cerr << "embd_dim\n"; return 1; }

        uint32_t count = 0;
        v.scan_representatives([&](uint32_t rid, std::string_view, uint32_t) {
            if (rid != count) { std::cerr << "rid order\n"; std::exit(1); }
            ++count;
        });
        if (count != 8) { std::cerr << "scan_count\n"; return 1; }

        // Rep itself
        char acc[16];
        std::snprintf(acc, sizeof(acc), "ACC_00000");
        auto st = v.status_for_accession(acc);
        if (st.kind != RepStatus::Kind::Representative) { std::cerr << "rep kind\n"; return 1; }

        // A member (index 8 = first member, rep is ACC_00000)
        std::snprintf(acc, sizeof(acc), "ACC_00008");
        st = v.status_for_accession(acc);
        if (st.kind != RepStatus::Kind::Member) { std::cerr << "mem kind\n"; return 1; }
        if (st.rep_accession != "ACC_00000")     { std::cerr << "mem rep_acc\n"; return 1; }

        // Unclus: index n_reps + n_members = 8 + 20 = 28
        std::snprintf(acc, sizeof(acc), "ACC_00028");
        st = v.status_for_accession(acc);
        if (st.kind != RepStatus::Kind::Unclustered) { std::cerr << "unc kind\n"; return 1; }

        // Tombstone: index 30
        std::snprintf(acc, sizeof(acc), "ACC_00030");
        st = v.status_for_accession(acc);
        if (st.kind != RepStatus::Kind::Tombstoned) { std::cerr << "tomb kind\n"; return 1; }

        // Absent
        st = v.status_for_accession("DOES_NOT_EXIST");
        if (st.kind != RepStatus::Kind::Absent) { std::cerr << "absent kind\n"; return 1; }

        std::remove(path.c_str());
    }

    {
        BuildOpts o;
        o.embd_dim = 4;
        auto buf = build_gpd(o);
        auto path = write_temp(buf);

        DerepView v;
        v.open(path);
        if (v.embedding_dim() != 4) { std::cerr << "embd_dim=4\n"; return 1; }
        for (uint32_t r = 0; r < 8; ++r) {
            float out[4] = {0,0,0,0};
            if (!v.embedding_for_rep(r, std::span<float>(out, 4))) {
                std::cerr << "embd_for_rep\n"; return 1;
            }
            for (uint16_t k = 0; k < 4; ++k) {
                float expected = static_cast<float>(r) + 0.25f * static_cast<float>(k);
                if (std::fabs(out[k] - expected) > 1e-2f) {
                    std::cerr << "embd value " << r << "," << k
                              << " got=" << out[k] << " exp=" << expected << "\n";
                    return 1;
                }
            }
        }
        std::remove(path.c_str());
    }

    {
        BuildOpts o;
        o.corrupt_section_crc = true;
        auto buf = build_gpd(o);
        auto path = write_temp(buf);
        DerepView v;
        bool threw = false;
        try { v.open(path); }
        catch (const std::exception&) { threw = true; }
        if (!threw) { std::cerr << "expected throw on corrupt CRC\n"; return 1; }
        std::remove(path.c_str());
    }

    std::cout << "derep_view ok\n";
    return 0;
}
