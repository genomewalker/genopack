#include <genopack/oph_sketch.hpp>
#include <algorithm>
#include <array>
#include <cstring>
#include <limits>
#ifdef __AVX2__
#include <immintrin.h>
#endif

namespace genopack {

namespace {

// Lookup table for base encoding: A=0, C=1, G=2, T=3, invalid=255
alignas(64) constexpr uint8_t BASE_ENCODE[256] = {
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,  0,255,  1,255,255,255,  2,255,255,255,255,255,255,255,255,
    255,255,255,255,  3,255,255,255,255,255,255,255,255,255,255,255,
    255,  0,255,  1,255,255,255,  2,255,255,255,255,255,255,255,255,
    255,255,255,255,  3,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
};

[[gnu::always_inline]] inline uint64_t encode_base(char c) {
    return BASE_ENCODE[static_cast<uint8_t>(c)];
}

template <uint32_t Bins>
[[gnu::always_inline]] inline uint32_t fast_range_u64(uint64_t h) {
    static_assert(Bins > 0);
    return static_cast<uint32_t>((static_cast<__uint128_t>(h) * Bins) >> 64);
}

[[gnu::always_inline]] inline uint32_t fast_range_u64(uint64_t h, uint32_t bins) {
    return static_cast<uint32_t>((static_cast<__uint128_t>(h) * bins) >> 64);
}

[[gnu::always_inline]] inline uint64_t wymix(uint64_t a, uint64_t b) {
    const __uint128_t product = static_cast<__uint128_t>(a) * b;
    return static_cast<uint64_t>(product) ^ static_cast<uint64_t>(product >> 64);
}

[[gnu::always_inline]] inline uint64_t oph_hash_wymix(uint64_t canonical, uint64_t seed) {
    constexpr uint64_t P0 = 0xa0761d6478bd642fULL;
    constexpr uint64_t P1 = 0xe7037ed1a0b428dbULL;
    return wymix(canonical ^ (seed + P0), canonical ^ P1);
}

[[gnu::always_inline]] inline uint32_t oph_sig32(uint64_t canonical_hash) {
    return static_cast<uint32_t>(canonical_hash >> 32);
}

struct OPHKmerState {
    uint64_t fwd = 0;
    uint64_t rev = 0;
    int valid = 0;
    void reset() { fwd = 0; rev = 0; valid = 0; }
};

inline void finalize_oph_sketch(OPHSketchResult& result, int m) {
    const uint32_t EMPTY = std::numeric_limits<uint32_t>::max();

    for (auto w : result.real_bins_bitmask)
        result.n_real_bins += static_cast<uint32_t>(__builtin_popcountll(w));

    auto mix = [](uint64_t x) -> uint32_t {
        x ^= x >> 30;
        x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27;
        x *= 0x94d049bb133111ebULL;
        x ^= x >> 31;
        return static_cast<uint32_t>(x >> 32);
    };
    for (int t = 1; t < m; ++t)
        if (result.signature[t] == EMPTY && result.signature[t - 1] != EMPTY)
            result.signature[t] = mix(
                static_cast<uint64_t>(result.signature[t - 1]) ^ static_cast<uint64_t>(t));
    for (int t = m - 2; t >= 0; --t)
        if (result.signature[t] == EMPTY && result.signature[t + 1] != EMPTY)
            result.signature[t] = mix(
                static_cast<uint64_t>(result.signature[t + 1]) ^ static_cast<uint64_t>(t));
}

template <uint32_t FixedBins>
[[gnu::always_inline]] inline void update_oph_bin(
        uint32_t* __restrict signature,
        uint64_t* __restrict real_bins,
        uint32_t runtime_bins,
        uint64_t canonical_hash) {
    const uint32_t bin = [&]() -> uint32_t {
        if constexpr (FixedBins != 0) return fast_range_u64<FixedBins>(canonical_hash);
        return fast_range_u64(canonical_hash, runtime_bins);
    }();

    const uint32_t h = oph_sig32(canonical_hash);
    uint32_t& cur = signature[bin];
    if (__builtin_expect(h < cur, 0)) {
        cur = h;
        real_bins[bin >> 6] |= (1ULL << (bin & 63));
    }
}

// AVX2 ACGT run scanner: returns length of valid ACGT run starting at data[0].
[[gnu::always_inline]] static inline size_t scan_valid_run(const char* data, size_t len) noexcept {
#ifdef __AVX2__
    if (__builtin_expect(len >= 32, 1)) {
        const __m256i v_mask = _mm256_set1_epi8(static_cast<char>(0xDF));
        const __m256i v_A    = _mm256_set1_epi8('A');
        const __m256i v_C    = _mm256_set1_epi8('C');
        const __m256i v_G    = _mm256_set1_epi8('G');
        const __m256i v_T    = _mm256_set1_epi8('T');
        size_t j = 0;
        for (; j + 32 <= len; j += 32) {
            __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(data + j));
            __m256i upper = _mm256_and_si256(chunk, v_mask);
            __m256i ok    = _mm256_or_si256(
                _mm256_or_si256(_mm256_cmpeq_epi8(upper, v_A), _mm256_cmpeq_epi8(upper, v_C)),
                _mm256_or_si256(_mm256_cmpeq_epi8(upper, v_G), _mm256_cmpeq_epi8(upper, v_T)));
            uint32_t bad = static_cast<uint32_t>(~_mm256_movemask_epi8(ok));
            if (bad != 0)
                return j + static_cast<size_t>(__builtin_ctz(bad));
        }
        for (; j < len; ++j)
            if (BASE_ENCODE[static_cast<uint8_t>(data[j])] == 255) return j;
        return len;
    }
#endif
    size_t j = 0;
    while (j < len && BASE_ENCODE[static_cast<uint8_t>(data[j])] != 255) ++j;
    return j;
}

// Process a pre-validated ACGT run [i, run_end).
template <uint32_t FixedBins>
[[gnu::always_inline]] inline void process_valid_run_bounded(
        const char* data,
        size_t i,
        size_t run_end,
        int k,
        uint64_t seed,
        uint64_t k_mask,
        int rev_shift,
        uint32_t runtime_bins,
        OPHKmerState& state,
        uint32_t* __restrict signature,
        uint64_t* __restrict real_bins) {
    uint64_t fwd   = state.fwd;
    uint64_t rev   = state.rev;
    int      valid = state.valid;
    size_t   j     = i;

    for (; j < run_end && valid < k; ++j) {
        const uint8_t base = BASE_ENCODE[static_cast<uint8_t>(data[j])];
        fwd = ((fwd << 2) | base) & k_mask;
        rev = (rev >> 2) | ((3ULL - base) << rev_shift);
        ++valid;
        if (valid >= k) {
            const uint64_t canonical = fwd ^ ((fwd ^ rev) & -(uint64_t)(fwd > rev));
            const uint64_t h = oph_hash_wymix(canonical, seed);
            update_oph_bin<FixedBins>(signature, real_bins, runtime_bins, h);
        }
    }

    for (; j < run_end; ++j) {
        const uint8_t base = BASE_ENCODE[static_cast<uint8_t>(data[j])];
        fwd = ((fwd << 2) | base) & k_mask;
        rev = (rev >> 2) | ((3ULL - base) << rev_shift);
        const uint64_t canonical = fwd ^ ((fwd ^ rev) & -(uint64_t)(fwd > rev));
        const uint64_t h = oph_hash_wymix(canonical, seed);
        update_oph_bin<FixedBins>(signature, real_bins, runtime_bins, h);
    }

    state.fwd   = fwd;
    state.rev   = rev;
    state.valid = valid;
}

constexpr uint32_t DEFAULT_OPH_BINS = 10000;

template <uint32_t FixedBins>
OPHSketchResult sketch_oph_impl(const char* data, size_t size,
                                 int m, const OPHSketchConfig& cfg) {
    const uint32_t EMPTY = std::numeric_limits<uint32_t>::max();

    OPHSketchResult result;
    result.genome_length = 0;
    result.signature.assign(m, EMPTY);
    result.real_bins_bitmask.assign((m + 63) / 64, 0ULL);

    const int k = cfg.kmer_size;
    if (k <= 0 || k > 32) return result;

    const uint64_t seed = cfg.seed;
    const uint64_t k_mask = (k == 32) ? UINT64_MAX : ((1ULL << (2 * k)) - 1);
    const int rev_shift = 2 * (k - 1);
    const uint32_t runtime_bins = static_cast<uint32_t>(m);

    OPHKmerState state;
    bool in_header = false;

    uint32_t* __restrict signature = result.signature.data();
    uint64_t* __restrict real_bins = result.real_bins_bitmask.data();

    size_t i = 0;
    while (i < size) {
        const char c = data[i];

        if (c == '>') {
            ++result.n_contigs;
            state.reset();
            in_header = true;
            const char* nl = static_cast<const char*>(std::memchr(data + i, '\n', size - i));
            if (nl) {
                i = static_cast<size_t>(nl - data) + 1;
                in_header = false;
            } else {
                break;
            }
            continue;
        }
        if (in_header) {
            if (c == '\n') in_header = false;
            ++i;
            continue;
        }
        if (c == '\n' || c == '\r') {
            ++i;
            continue;
        }

        if (encode_base(c) == 255) {
            state.reset();
            ++i;
            continue;
        }

        const size_t run_start = i;
        const size_t run_end = i + scan_valid_run(data + i, size - i);
        process_valid_run_bounded<FixedBins>(
            data, i, run_end, k, seed, k_mask, rev_shift, runtime_bins,
            state, signature, real_bins);
        i = run_end;

        result.genome_length += static_cast<uint64_t>(i - run_start);

        if (i < size) {
            const char stop = data[i];
            if (stop != '\n' && stop != '\r' && stop != '>') {
                state.reset();
                ++i;
            }
        }
    }

    finalize_oph_sketch(result, m);
    return result;
}

} // anonymous namespace

OPHSketchResult sketch_oph_from_buffer(const char* data, size_t len,
                                        const OPHSketchConfig& cfg) {
    const int m = cfg.sketch_size;
    if (m == static_cast<int>(DEFAULT_OPH_BINS))
        return sketch_oph_impl<DEFAULT_OPH_BINS>(data, len, m, cfg);
    return sketch_oph_impl<0>(data, len, m, cfg);
}

// ── Dual-seed OPH: one k-mer scan, two signatures ───────────────────────────
// The single scanner walks canonical k-mers; each k-mer is hashed twice
// (once per seed) and placed into its per-seed (bin, value) slot. Bins are
// chosen from the seed=seed1 hash to keep sig1 byte-identical to what a
// single-seed run with seed=seed1 would produce; sig2 uses the seed=seed2
// hash for BOTH the bin assignment and the value — i.e. sig2 is byte-equal
// to an independent OPH run with seed=seed2. real_bins_bitmask tracks any
// observed bin (either seed).

namespace {

inline void finalize_signature(std::vector<uint32_t>& signature, int m) {
    const uint32_t EMPTY = std::numeric_limits<uint32_t>::max();
    auto mix = [](uint64_t x) -> uint32_t {
        x ^= x >> 30;
        x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27;
        x *= 0x94d049bb133111ebULL;
        x ^= x >> 31;
        return static_cast<uint32_t>(x >> 32);
    };
    for (int t = 1; t < m; ++t)
        if (signature[t] == EMPTY && signature[t - 1] != EMPTY)
            signature[t] = mix(static_cast<uint64_t>(signature[t - 1]) ^ static_cast<uint64_t>(t));
    for (int t = m - 2; t >= 0; --t)
        if (signature[t] == EMPTY && signature[t + 1] != EMPTY)
            signature[t] = mix(static_cast<uint64_t>(signature[t + 1]) ^ static_cast<uint64_t>(t));
}

// Shared dual-seed bin update: bin index from its own seed hash; value from
// the same hash. This matches what two independent sketch_oph_impl runs
// would produce because each seed sees the same canonical k-mer stream.
[[gnu::always_inline]] inline void update_dual_bin(
        uint32_t* __restrict sig1, uint64_t* __restrict real_bins1,
        uint32_t* __restrict sig2, uint64_t* __restrict real_bins2,
        uint64_t h1, uint64_t h2, uint32_t bins) {
    const uint32_t bin1 = fast_range_u64(h1, bins);
    const uint32_t bin2 = fast_range_u64(h2, bins);
    const uint32_t v1 = oph_sig32(h1);
    const uint32_t v2 = oph_sig32(h2);
    if (__builtin_expect(v1 < sig1[bin1], 0)) {
        sig1[bin1] = v1;
        real_bins1[bin1 >> 6] |= (1ULL << (bin1 & 63));
    }
    if (__builtin_expect(v2 < sig2[bin2], 0)) {
        sig2[bin2] = v2;
        real_bins2[bin2 >> 6] |= (1ULL << (bin2 & 63));
    }
}

template <uint32_t FixedBins>
void process_valid_run_dual(
        const char* data,
        size_t i, size_t run_end,
        int k,
        uint64_t seed1, uint64_t seed2,
        uint64_t k_mask, int rev_shift,
        uint32_t runtime_bins,
        OPHKmerState& state,
        uint32_t* __restrict sig1, uint64_t* __restrict real_bins1,
        uint32_t* __restrict sig2, uint64_t* __restrict real_bins2) {
    uint64_t fwd   = state.fwd;
    uint64_t rev   = state.rev;
    int      valid = state.valid;
    size_t   j     = i;

    for (; j < run_end && valid < k; ++j) {
        const uint8_t base = BASE_ENCODE[static_cast<uint8_t>(data[j])];
        fwd = ((fwd << 2) | base) & k_mask;
        rev = (rev >> 2) | ((3ULL - base) << rev_shift);
        ++valid;
        if (valid >= k) {
            const uint64_t canonical = fwd ^ ((fwd ^ rev) & -(uint64_t)(fwd > rev));
            const uint64_t h1 = oph_hash_wymix(canonical, seed1);
            const uint64_t h2 = oph_hash_wymix(canonical, seed2);
            if constexpr (FixedBins != 0)
                update_dual_bin(sig1, real_bins1, sig2, real_bins2, h1, h2, FixedBins);
            else
                update_dual_bin(sig1, real_bins1, sig2, real_bins2, h1, h2, runtime_bins);
        }
    }
    for (; j < run_end; ++j) {
        const uint8_t base = BASE_ENCODE[static_cast<uint8_t>(data[j])];
        fwd = ((fwd << 2) | base) & k_mask;
        rev = (rev >> 2) | ((3ULL - base) << rev_shift);
        const uint64_t canonical = fwd ^ ((fwd ^ rev) & -(uint64_t)(fwd > rev));
        const uint64_t h1 = oph_hash_wymix(canonical, seed1);
        const uint64_t h2 = oph_hash_wymix(canonical, seed2);
        if constexpr (FixedBins != 0)
            update_dual_bin(sig1, real_bins1, sig2, real_bins2, h1, h2, FixedBins);
        else
            update_dual_bin(sig1, real_bins1, sig2, real_bins2, h1, h2, runtime_bins);
    }
    state.fwd   = fwd;
    state.rev   = rev;
    state.valid = valid;
}

template <uint32_t FixedBins>
OPHDualSketchResult sketch_oph_dual_impl(const char* data, size_t size,
                                          int m, int k,
                                          uint64_t seed1, uint64_t seed2) {
    const uint32_t EMPTY = std::numeric_limits<uint32_t>::max();

    OPHDualSketchResult result;
    result.signature1.assign(m, EMPTY);
    result.signature2.assign(m, EMPTY);
    result.real_bins_bitmask.assign((m + 63) / 64, 0ULL);

    if (k <= 0 || k > 32) return result;

    // Work with per-seed bitmasks during the scan; union them at the end
    // so real_bins_bitmask reflects "any bin hit by any seed" (equivalent
    // to the single-seed mask in practice, since all bins are covered when
    // a genome is long enough to saturate OPH).
    std::vector<uint64_t> rbins1(result.real_bins_bitmask.size(), 0ULL);
    std::vector<uint64_t> rbins2(result.real_bins_bitmask.size(), 0ULL);

    const uint64_t k_mask = (k == 32) ? UINT64_MAX : ((1ULL << (2 * k)) - 1);
    const int rev_shift = 2 * (k - 1);
    const uint32_t runtime_bins = static_cast<uint32_t>(m);

    OPHKmerState state;
    bool in_header = false;

    uint32_t* __restrict s1 = result.signature1.data();
    uint32_t* __restrict s2 = result.signature2.data();
    uint64_t* __restrict rb1 = rbins1.data();
    uint64_t* __restrict rb2 = rbins2.data();

    size_t i = 0;
    while (i < size) {
        const char c = data[i];
        if (c == '>') {
            ++result.n_contigs;
            state.reset();
            in_header = true;
            const char* nl = static_cast<const char*>(std::memchr(data + i, '\n', size - i));
            if (nl) {
                i = static_cast<size_t>(nl - data) + 1;
                in_header = false;
            } else {
                break;
            }
            continue;
        }
        if (in_header) { if (c == '\n') in_header = false; ++i; continue; }
        if (c == '\n' || c == '\r') { ++i; continue; }
        if (encode_base(c) == 255) { state.reset(); ++i; continue; }

        const size_t run_start = i;
        const size_t run_end = i + scan_valid_run(data + i, size - i);
        process_valid_run_dual<FixedBins>(
            data, i, run_end, k, seed1, seed2, k_mask, rev_shift, runtime_bins,
            state, s1, rb1, s2, rb2);
        i = run_end;
        result.genome_length += static_cast<uint64_t>(i - run_start);
        if (i < size) {
            const char stop = data[i];
            if (stop != '\n' && stop != '\r' && stop != '>') { state.reset(); ++i; }
        }
    }

    // Union the per-seed bitmasks. In practice every bin ends up covered in
    // both when there are enough k-mers, but the union is the safest
    // seed-independent representation of "real" bins.
    for (size_t w = 0; w < result.real_bins_bitmask.size(); ++w)
        result.real_bins_bitmask[w] = rbins1[w] | rbins2[w];
    for (auto w : result.real_bins_bitmask)
        result.n_real_bins += static_cast<uint32_t>(__builtin_popcountll(w));

    finalize_signature(result.signature1, m);
    finalize_signature(result.signature2, m);
    return result;
}

} // anonymous namespace

OPHDualSketchResult sketch_oph_dual_from_buffer(const char* data, size_t len,
                                                int kmer_size, int sketch_size,
                                                int /*syncmer_s*/,
                                                uint64_t seed1, uint64_t seed2) {
    if (sketch_size == static_cast<int>(DEFAULT_OPH_BINS))
        return sketch_oph_dual_impl<DEFAULT_OPH_BINS>(data, len, sketch_size, kmer_size, seed1, seed2);
    return sketch_oph_dual_impl<0>(data, len, sketch_size, kmer_size, seed1, seed2);
}

} // namespace genopack
