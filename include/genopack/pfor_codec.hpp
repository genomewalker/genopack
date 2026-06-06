#pragma once
// ── PFOR-delta codec for sorted uint64 sequences ──────────────────────────────
// Frame-of-reference, block-adaptive delta coding for the PCORE v1 aamer runs.
// Each genus run is a strictly-ascending, unique sequence of FMH-reduced uint64
// hashes; gaps are ~uniform, so per-block bit-packed deltas approach the gap
// entropy (~log2(U/n) bits/elem) at ~2x the naive 64-bit cost, with zero query
// penalty because PCORE is consumed by sequential scan only.
//
// Stream layout (single continuous LSB-first bitstream per run):
//   put(v[0], 64)                       // first element, absolute
//   for each group of up to BLK deltas:
//     put(nbits, 7)                     // bit-width for this group (1..64)
//     put(delta, nbits) x group_size    // delta[k] = v[k]-v[k-1]  (>=1)
// The decoder reconstructs v[] given the element count n (stored in PcoreEntry).
#include <cstdint>
#include <cstddef>
#include <vector>

namespace genopack::pfor {

inline constexpr size_t BLK = 128;

inline int bit_len(uint64_t x) noexcept {
    int b = 0;
    while (x) { ++b; x >>= 1; }
    return b;
}

// LSB-first bit writer over a byte vector. A 128-bit accumulator lets a single
// put() safely append up to 64 bits without losing carry across the 64-bit line.
struct BitWriter {
    std::vector<uint8_t>& out;
    unsigned __int128 acc = 0;
    int nbits = 0;
    explicit BitWriter(std::vector<uint8_t>& o) : out(o) {}
    void put(uint64_t v, int bits) {              // caller guarantees v < 2^bits
        if (bits == 0) return;
        acc |= (static_cast<unsigned __int128>(v)) << nbits;
        nbits += bits;
        while (nbits >= 8) { out.push_back(static_cast<uint8_t>(acc & 0xff)); acc >>= 8; nbits -= 8; }
    }
    void flush() { if (nbits > 0) { out.push_back(static_cast<uint8_t>(acc & 0xff)); acc = 0; nbits = 0; } }
};

struct BitReader {
    const uint8_t* p;
    size_t nbytes;
    size_t bytepos = 0;
    unsigned __int128 acc = 0;
    int nbits = 0;
    BitReader(const uint8_t* d, size_t n) : p(d), nbytes(n) {}
    uint64_t get(int bits) {
        if (bits == 0) return 0;
        while (nbits < bits) {
            const uint64_t b = (bytepos < nbytes) ? p[bytepos++] : 0;
            acc |= (static_cast<unsigned __int128>(b)) << nbits;
            nbits += 8;
        }
        const unsigned __int128 mask = ((static_cast<unsigned __int128>(1) << bits) - 1);
        const uint64_t v = static_cast<uint64_t>(acc & mask);
        acc >>= bits;
        nbits -= bits;
        return v;
    }
};

// Encode a strictly-ascending unique sequence vals[0..n). Returns encoded bytes.
inline std::vector<uint8_t> encode_sorted(const uint64_t* vals, size_t n) {
    std::vector<uint8_t> out;
    if (n == 0) return out;
    out.reserve(8 + n);                            // generous; gaps usually <8 bytes/elem
    BitWriter bw(out);
    bw.put(vals[0], 64);
    size_t k = 1;
    while (k < n) {
        const size_t m = (n - k < BLK) ? (n - k) : BLK;   // deltas in this group
        uint64_t maxd = 0;
        for (size_t j = 0; j < m; ++j) {
            const uint64_t d = vals[k + j] - vals[k + j - 1];
            if (d > maxd) maxd = d;
        }
        const int nb = bit_len(maxd);              // >=1 (unique ascending => d>=1)
        bw.put(static_cast<uint64_t>(nb), 7);
        for (size_t j = 0; j < m; ++j)
            bw.put(vals[k + j] - vals[k + j - 1], nb);
        k += m;
    }
    bw.flush();
    return out;
}

inline std::vector<uint8_t> encode_sorted(const std::vector<uint64_t>& v) {
    return encode_sorted(v.data(), v.size());
}

// Decode n elements from the bitstream into out (cleared first).
inline void decode_into(const uint8_t* data, size_t nbytes, size_t n, std::vector<uint64_t>& out) {
    out.clear();
    if (n == 0) return;
    out.resize(n);
    BitReader br(data, nbytes);
    out[0] = br.get(64);
    size_t k = 1;
    while (k < n) {
        const size_t m = (n - k < BLK) ? (n - k) : BLK;
        const int nb = static_cast<int>(br.get(7));
        for (size_t j = 0; j < m; ++j) out[k + j] = out[k + j - 1] + br.get(nb);
        k += m;
    }
}

inline std::vector<uint64_t> decode(const uint8_t* data, size_t nbytes, size_t n) {
    std::vector<uint64_t> out;
    decode_into(data, nbytes, n, out);
    return out;
}

} // namespace genopack::pfor
