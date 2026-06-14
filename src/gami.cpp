#include <genopack/gami.hpp>
#include <genopack/mmap_file.hpp>   // AppendWriter
#include <genopack/format.hpp>      // SectionDesc, SEC_GAMI
#include <zstd.h>
#include <sys/mman.h>
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace genopack {

SectionDesc GamiV2Writer::finalize(AppendWriter& w, uint64_t section_id,
                                   const GlobalMultiplicityIndex& gmi,
                                   uint64_t frac_max) {
    // Collect all (hash, count) pairs, skip empty slots.
    struct Pair { uint64_t h; uint8_t c; };
    std::vector<Pair> pairs;
    pairs.reserve(gmi.count());
    const auto& keys = gmi.keys_;
    const auto& vals = gmi.vals_;
    for (size_t i = 0; i < keys.size(); ++i)
        if (keys[i] != GlobalMultiplicityIndex::kEmpty)
            pairs.push_back({keys[i], vals[i]});

    std::sort(pairs.begin(), pairs.end(), [](const Pair& a, const Pair& b){ return a.h < b.h; });

    // Serialize: split into parallel arrays for better compression.
    // Layout: [N × uint64_t hashes][N × uint8_t counts]
    const size_t N = pairs.size();
    std::vector<uint8_t> raw(N * sizeof(uint64_t) + N);
    auto* hptr = reinterpret_cast<uint64_t*>(raw.data());
    uint8_t* cptr = raw.data() + N * sizeof(uint64_t);
    for (size_t i = 0; i < N; ++i) { hptr[i] = pairs[i].h; cptr[i] = pairs[i].c; }

    const size_t cbound = ZSTD_compressBound(raw.size());
    std::vector<uint8_t> compressed(cbound);
    const size_t csz = ZSTD_compress(compressed.data(), cbound, raw.data(), raw.size(), 3);
    if (ZSTD_isError(csz))
        throw std::runtime_error(std::string("gami_v2: zstd compress failed: ") + ZSTD_getErrorName(csz));
    compressed.resize(csz);

    GamiHeaderV2 hdr;
    hdr.n_entries     = static_cast<uint64_t>(N);
    hdr.frac_max      = frac_max;
    hdr.payload_bytes = static_cast<uint64_t>(csz);

    w.align(8);
    const uint64_t sec_start = w.current_offset();
    w.append(&hdr, sizeof(hdr));
    w.append(compressed.data(), compressed.size());
    const uint64_t sec_end = w.current_offset();

    SectionDesc sd{};
    sd.type              = SEC_GAMI;
    sd.version           = 2;
    sd.section_id        = section_id;
    sd.file_offset       = sec_start;
    sd.compressed_size   = sec_end - sec_start;
    sd.uncompressed_size = sec_end - sec_start;
    sd.item_count        = static_cast<uint64_t>(N);
    sd.header_size       = sizeof(GamiHeaderV2);
    return sd;
}

void gami_v2_load(const uint8_t* section_data, uint64_t section_size,
                  GlobalMultiplicityIndex& out) {
    if (section_size < sizeof(GamiHeaderV2))
        throw std::runtime_error("gami_v2_load: section too small");
    const auto* hdr = reinterpret_cast<const GamiHeaderV2*>(section_data);
    if (hdr->magic != GAMI_MAGIC_V2)
        throw std::runtime_error("gami_v2_load: bad magic");

    const uint64_t N   = hdr->n_entries;
    const uint64_t psz = hdr->payload_bytes;
    if (sizeof(GamiHeaderV2) + psz > section_size)
        throw std::runtime_error("gami_v2_load: truncated section");

    // Pre-warm NFS mmap pages sequentially so ZSTD doesn't page-fault 4 KB at a time.
    ::madvise(const_cast<void*>(static_cast<const void*>(section_data + sizeof(GamiHeaderV2))),
              static_cast<size_t>(psz), MADV_SEQUENTIAL);

    const size_t raw_size = N * sizeof(uint64_t) + N;
    std::vector<uint8_t> raw(raw_size);
    const size_t dsz = ZSTD_decompress(raw.data(), raw_size,
                                       section_data + sizeof(GamiHeaderV2), psz);
    if (ZSTD_isError(dsz))
        throw std::runtime_error(std::string("gami_v2_load: zstd decompress failed: ") + ZSTD_getErrorName(dsz));
    if (dsz != raw_size)
        throw std::runtime_error("gami_v2_load: size mismatch after decompress");

    const size_t n = static_cast<size_t>(N);
    std::vector<uint64_t> hashes(n);
    std::memcpy(hashes.data(), raw.data(), N * sizeof(uint64_t));
    out.ef_.build(hashes.data(), n, UINT64_MAX);

    const uint8_t* src = raw.data() + N * sizeof(uint64_t);
    out.packed_vals_.assign((n + 31) / 32, 0);
    for (size_t i = 0; i < n; ++i) {
        const uint8_t gc = src[i];
        write_pv(out.packed_vals_.data(), i, gc >= 3 ? uint8_t(3) : gc);
    }
}

} // namespace genopack
