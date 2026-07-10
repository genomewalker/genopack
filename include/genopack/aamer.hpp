#pragma once
#include <algorithm>
#include <cstdint>
#include <functional>
#include <string_view>
#include <vector>
#ifdef __AVX2__
#  include <immintrin.h>
#endif

namespace genopack {

// Amino acid k=8 aamer primitives for SCG completeness/contamination estimation.
// Build side: protein input from GTDB-Tk MSA (extract_aamers_aa).
// Check side: DNA input via 6-frame translation (extract_aamers_dna / translate_6frame).
// Hash encoding is identical on both sides — required for membership lookup to work.

static constexpr int AAMER_K = 8;

// ── Encoding tables ───────────────────────────────────────────────────────────

// AA index: alphabetical single-letter order (A=0,C=1,D=2,E=3,F=4,G=5,H=6,I=7,
//   K=8,L=9,M=10,N=11,P=12,Q=13,R=14,S=15,T=16,V=17,W=18,Y=19). 0xFF = stop/invalid.
static constexpr uint8_t AA_STOP = 0xFF;
static constexpr uint8_t AA_COUNT = 20;

// Char → AA index (handles upper and lower case, '-', '*', 'X').
// clang-format off
alignas(64) static constexpr uint8_t AA_ENC[256] = {
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF, // 0x00
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF, // 0x10
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF, // 0x20  '-'=0xFF
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF, // 0x30
    0xFF,   0,0xFF,   1,   2,   3,   4,   5,   6,   7,0xFF,   8,   9,  10,  11,0xFF, // A..O
       12,  13,  14,  15,  16,   1,  17,  18,0xFF,  19,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF, // P..Z  U→C
    0xFF,   0,0xFF,   1,   2,   3,   4,   5,   6,   7,0xFF,   8,   9,  10,  11,0xFF, // a..o
       12,  13,  14,  15,  16,   1,  17,  18,0xFF,  19,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF, // p..z
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
};
// clang-format on

// DNA base encoding: A=0, C=1, G=2, T/U=3, ambiguous/other=0xFF.
// clang-format off
alignas(64) static constexpr uint8_t DNA_ENC[256] = {
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,   0,0xFF,   1,0xFF,0xFF,0xFF,   2,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF, // A C G
    0xFF,0xFF,0xFF,0xFF,   3,   3,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF, // T U
    0xFF,   0,0xFF,   1,0xFF,0xFF,0xFF,   2,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF, // a c g
    0xFF,0xFF,0xFF,0xFF,   3,   3,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF, // t u
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
};
// clang-format on

// Complement of encoded base: A↔T (0↔3), C↔G (1↔2).
static constexpr uint8_t BASE_COMP[4] = { 3, 2, 1, 0 };

// ── Codon table ───────────────────────────────────────────────────────────────

// Genetic code 11 (bacteria/archaea). Index = first_base*16 + second_base*4 + third_base.
// A=0,C=1,G=2,T=3. AA values per earlier encoding.
// clang-format off
static constexpr uint8_t CODON11[64] = {
    // b0=A (indices 0-15): AAA=K,AAC=N,AAG=K,AAT=N, ACA=T,ACC=T,ACG=T,ACT=T,
    //                       AGA=R,AGC=S,AGG=R,AGT=S, ATA=I,ATC=I,ATG=M,ATT=I
     8, 11,  8, 11,   16, 16, 16, 16,   14, 15, 14, 15,    7,  7, 10,  7,
    // b0=C (indices 16-31): CAA=Q,CAC=H,CAG=Q,CAT=H, CCA=P,CCC=P,CCG=P,CCT=P,
    //                        CGA=R,CGC=R,CGG=R,CGT=R, CTA=L,CTC=L,CTG=L,CTT=L
    13,  6, 13,  6,   12, 12, 12, 12,   14, 14, 14, 14,    9,  9,  9,  9,
    // b0=G (indices 32-47): GAA=E,GAC=D,GAG=E,GAT=D, GCA=A,GCC=A,GCG=A,GCT=A,
    //                        GGA=G,GGC=G,GGG=G,GGT=G, GTA=V,GTC=V,GTG=V,GTT=V
     3,  2,  3,  2,    0,  0,  0,  0,    5,  5,  5,  5,   17, 17, 17, 17,
    // b0=T (indices 48-63): TAA=*,TAC=Y,TAG=*,TAT=Y, TCA=S,TCC=S,TCG=S,TCT=S,
    //                        TGA=*,TGC=C,TGG=W,TGT=C, TTA=L,TTC=F,TTG=L,TTT=F
    0xFF,19,0xFF,19,   15, 15, 15, 15, 0xFF,  1, 18,  1,    9,  4,  9,  4,
};
// clang-format on

// ── Reduced alphabets ────────────────────────────────────────────────────────

// Murphy 10-class reduction (index = AA_ENC value, result = class 0..9).
// Groups: {LVIM}=0  {C}=1  {A}=2  {G}=3  {ST}=4  {P}=5  {FYW}=6  {EDNQ}=7  {KR}=8  {H}=9
// 10^8 = 100M k-mer space (vs 6^8=1.7M for Dayhoff6) — survives cross-family filter.
static constexpr uint8_t AA_MURPHY10[20] = {
    2, 1, 7, 7, 6, 3, 9, 0, 8, 0, 0, 7, 5, 7, 8, 4, 4, 0, 6, 6
};

// ── Hash ─────────────────────────────────────────────────────────────────────

// Hash one k=8 amino acid k-mer (aa[0..7] in 0-19 encoding).
// Packs 8 AAs at 5 bits each (40 bits total) then applies splitmix64 mixing.
// Must be called with identical encoding on build side and check side.
[[nodiscard]] inline uint64_t aamer_hash(const uint8_t* aa) noexcept {
    uint64_t v =  (uint64_t)aa[0]
               | ((uint64_t)aa[1] <<  5)
               | ((uint64_t)aa[2] << 10)
               | ((uint64_t)aa[3] << 15)
               | ((uint64_t)aa[4] << 20)
               | ((uint64_t)aa[5] << 25)
               | ((uint64_t)aa[6] << 30)
               | ((uint64_t)aa[7] << 35);
    v ^= v >> 30;
    v *= 0xbf58476d1ce4e5b9ULL;
    v ^= v >> 27;
    v *= 0x94d049bb133111ebULL;
    v ^= v >> 31;
    return v;
}

// Murphy 10-class variant: 10 classes × 4 bits = 32-bit packed input, then splitmix64.
// 10^8 = 100M k-mer space — large enough to survive without cross-family filtering.
[[nodiscard]] inline uint64_t aamer_hash_murphy10(const uint8_t* aa) noexcept {
    uint64_t v =  (uint64_t)AA_MURPHY10[aa[0]]
               | ((uint64_t)AA_MURPHY10[aa[1]] <<  4)
               | ((uint64_t)AA_MURPHY10[aa[2]] <<  8)
               | ((uint64_t)AA_MURPHY10[aa[3]] << 12)
               | ((uint64_t)AA_MURPHY10[aa[4]] << 16)
               | ((uint64_t)AA_MURPHY10[aa[5]] << 20)
               | ((uint64_t)AA_MURPHY10[aa[6]] << 24)
               | ((uint64_t)AA_MURPHY10[aa[7]] << 28);
    v ^= v >> 30;
    v *= 0xbf58476d1ce4e5b9ULL;
    v ^= v >> 27;
    v *= 0x94d049bb133111ebULL;
    v ^= v >> 31;
    return v;
}

// Returns true if the k-mer contains a run of >= (k+1)/2 identical amino acids.
[[nodiscard]] inline bool aamer_is_low_complexity(const uint8_t* aa, int k) noexcept {
    int run = 1, max_run = 1;
    for (int i = 1; i < k; ++i) {
        run = (aa[i] == aa[i-1]) ? run + 1 : 1;
        if (run > max_run) max_run = run;
    }
    return max_run >= (k + 1) / 2;
}

// Append all k-mers from one inter-stop AA segment to `out`, skipping low-complexity.
inline void emit_aamers(const uint8_t* seg, int len, int k,
                          std::vector<uint64_t>& out) {
    for (int i = 0; i + k <= len; ++i) {
        if (!aamer_is_low_complexity(seg + i, k))
            out.push_back(aamer_hash(seg + i));
    }
}

// FracMinHash variant with rolling packed-value: avoids re-packing 8 bytes per slide.
// v_{i+1} = (v_i - seg[i]) >> 5 | (seg[i+k] << (k-1)*5)  — 3 ops vs 22 for full repack.
// Hash output is bit-identical to aamer_hash(seg+i); GMI index is unaffected.
// seg[i] values are 0..19 < 32, so low 5 bits of v always equal seg[i] (no borrow on subtract).

inline void emit_aamers_frac(const uint8_t* seg, int len, int k,
                               uint64_t max_hash, std::vector<uint64_t>& out) {
    if (len < k) return;
    uint64_t v = 0;
    for (int j = 0; j < k; ++j) v |= (uint64_t)seg[j] << (j * 5);
    const int tail_shift = (k - 1) * 5;
    for (int i = 0; i + k <= len; ++i) {
        // Hash-gate FIRST: the FMH subsample (h <= max_hash) rejects ~1-1/N of
        // k-mers, so the low-complexity scan runs only on the ~1/N survivors.
        // Output is identical (both predicates must hold; neither has side effects).
        uint64_t h = v;
        h ^= h >> 30; h *= 0xbf58476d1ce4e5b9ULL;
        h ^= h >> 27; h *= 0x94d049bb133111ebULL;
        h ^= h >> 31;
        if (h <= max_hash && !aamer_is_low_complexity(seg + i, k))
            out.push_back(h);
        if (i + k < len)
            v = ((v - (uint64_t)seg[i]) >> 5) | ((uint64_t)seg[i + k] << tail_shift);
    }
}

// ── SIMD helpers for 6-frame translation (AVX2/SSSE3/SSE4.1) ─────────────────
#ifdef __AVX2__
namespace aamer_detail {

// CODON11 split into four 16-byte SSSE3 shuffle tables indexed by (codon & 0x0F).
// Quarter selected by (codon >> 4). Stops (0xFF in scalar) stored as 0x7F to keep
// them distinct from the 0xFF ambiguous-base sentinel used in SIMD logic.
alignas(16) static constexpr int8_t C11Q0[16] = {  // b0=A, indices 0..15
     8, 11,  8, 11, 16, 16, 16, 16, 14, 15, 14, 15,  7,  7, 10,  7};
alignas(16) static constexpr int8_t C11Q1[16] = {  // b0=C, indices 16..31
    13,  6, 13,  6, 12, 12, 12, 12, 14, 14, 14, 14,  9,  9,  9,  9};
alignas(16) static constexpr int8_t C11Q2[16] = {  // b0=G, indices 32..47
     3,  2,  3,  2,  0,  0,  0,  0,  5,  5,  5,  5, 17, 17, 17, 17};
alignas(16) static constexpr int8_t C11Q3[16] = {  // b0=T, indices 48..63 (TAA/TAG/TGA→0x7F)
    0x7F, 19, 0x7F, 19, 15, 15, 15, 15, 0x7F, 1, 18, 1, 9, 4, 9, 4};

// Deinterleave 24 input bytes (lo=p[0..15], hi=p[8..23]) into three 8-wide base streams.
// Byte positions: b0 at {0,3,6,9,12,15,18,21}, b1 at {1,4,7,10,13,16,19,22},
//                b2 at {2,5,8,11,14,17,20,23}.  Lanes 8..15 zeroed (0x80 mask).
alignas(16) static constexpr int8_t DIB0LO[16] = {
     0,  3,  6,  9, 12, 15, -128,-128,-128,-128,-128,-128,-128,-128,-128,-128};
alignas(16) static constexpr int8_t DIB0HI[16] = {
    -128,-128,-128,-128,-128,-128, 10, 13,-128,-128,-128,-128,-128,-128,-128,-128};
alignas(16) static constexpr int8_t DIB1LO[16] = {
     1,  4,  7, 10, 13,-128,-128,-128,-128,-128,-128,-128,-128,-128,-128,-128};
alignas(16) static constexpr int8_t DIB1HI[16] = {
    -128,-128,-128,-128,-128,  8, 11, 14,-128,-128,-128,-128,-128,-128,-128,-128};
alignas(16) static constexpr int8_t DIB2LO[16] = {
     2,  5,  8, 11, 14,-128,-128,-128,-128,-128,-128,-128,-128,-128,-128,-128};
alignas(16) static constexpr int8_t DIB2HI[16] = {
    -128,-128,-128,-128,-128,  9, 12, 15,-128,-128,-128,-128,-128,-128,-128,-128};

// Encode 16 ASCII DNA bytes → 0(A),1(C),2(G),3(T/U) or 0xFF (invalid/ambiguous).
// Clears bit 5 only to normalise case (a→A, t→T, u→U) then does exact comparison.
static inline __m128i enc16(const __m128i v) noexcept {
    // 0xDF = 1101 1111: clears bit 5 only → maps acgtu → ACGTU, high-bit bytes
    // (e.g. 0xC1) keep bit 7 and stay rejected — matches the scalar DNA_ENC table.
    const __m128i upr  = _mm_and_si128(v, _mm_set1_epi8((char)0xDF));
    const __m128i isA  = _mm_cmpeq_epi8(upr, _mm_set1_epi8('A'));
    const __m128i isC  = _mm_cmpeq_epi8(upr, _mm_set1_epi8('C'));
    const __m128i isG  = _mm_cmpeq_epi8(upr, _mm_set1_epi8('G'));
    const __m128i isT  = _mm_cmpeq_epi8(upr, _mm_set1_epi8('T'));
    const __m128i isU  = _mm_cmpeq_epi8(upr, _mm_set1_epi8('U'));
    const __m128i isTU = _mm_or_si128(isT, isU);
    const __m128i base = _mm_or_si128(
        _mm_or_si128(_mm_and_si128(isC, _mm_set1_epi8(1)),
                     _mm_and_si128(isG, _mm_set1_epi8(2))),
        _mm_and_si128(isTU, _mm_set1_epi8(3)));
    const __m128i valid = _mm_or_si128(_mm_or_si128(isA, isC), _mm_or_si128(isG, isTU));
    return _mm_or_si128(base, _mm_andnot_si128(valid, _mm_set1_epi8(-1)));
}

// 4-shuffle CODON11 lookup: 8 6-bit indices → 8 AA bytes (0..19) or 0x7F (stop).
static inline __m128i codon8(const __m128i idx6) noexcept {
    const __m128i lo4 = _mm_and_si128(idx6, _mm_set1_epi8(0x0F));
    const __m128i q0  = _mm_shuffle_epi8(_mm_load_si128(reinterpret_cast<const __m128i*>(C11Q0)), lo4);
    const __m128i q1  = _mm_shuffle_epi8(_mm_load_si128(reinterpret_cast<const __m128i*>(C11Q1)), lo4);
    const __m128i q2  = _mm_shuffle_epi8(_mm_load_si128(reinterpret_cast<const __m128i*>(C11Q2)), lo4);
    const __m128i q3  = _mm_shuffle_epi8(_mm_load_si128(reinterpret_cast<const __m128i*>(C11Q3)), lo4);
    // Per-byte >>4 without lane bleed: mask to high 3 bits (idx<64 so bit6=0 always),
    // then 16-bit right-shift 4 — the zero lower nibble prevents cross-byte contamination.
    const __m128i hi2 = _mm_srli_epi16(_mm_and_si128(idx6, _mm_set1_epi8(0x70)), 4);
    const __m128i is1 = _mm_cmpeq_epi8(hi2, _mm_set1_epi8(1));
    const __m128i is2 = _mm_cmpeq_epi8(hi2, _mm_set1_epi8(2));
    const __m128i is3 = _mm_cmpeq_epi8(hi2, _mm_set1_epi8(3));
    __m128i r = _mm_blendv_epi8(q0, q1, is1);
    r = _mm_blendv_epi8(r, q2, is2);
    r = _mm_blendv_epi8(r, q3, is3);
    return r;
}

// Translate 8 codons at p[0..23]. Writes 8 AA bytes to out[0..7].
// Returns 8-bit mask: bit i set → out[i] is a valid AA (not stop/ambiguous).
// Caller must ensure p[0..23] are all readable.
static inline unsigned xlate8(const char* p, uint8_t* out) noexcept {
    // Encode ONCE per input vector, then deinterleave the encoded bytes.
    // enc16 is per-byte so it commutes with the shuffle permutation: this turns
    // 3 enc16 (15 cmpeq) into 2 (10 cmpeq). The -128 shuffle indices zero inactive
    // lanes 8..15 (to 0x00 here vs 0xFF before) — those lanes are unused (storel
    // writes only 0..7, brk mask is &0xFF), so output on lanes 0..7 is identical.
    const __m128i elo = enc16(_mm_loadu_si128(reinterpret_cast<const __m128i*>(p)));
    const __m128i ehi = enc16(_mm_loadu_si128(reinterpret_cast<const __m128i*>(p + 8)));
    const __m128i b0 = _mm_or_si128(
        _mm_shuffle_epi8(elo, _mm_load_si128(reinterpret_cast<const __m128i*>(DIB0LO))),
        _mm_shuffle_epi8(ehi, _mm_load_si128(reinterpret_cast<const __m128i*>(DIB0HI))));
    const __m128i b1 = _mm_or_si128(
        _mm_shuffle_epi8(elo, _mm_load_si128(reinterpret_cast<const __m128i*>(DIB1LO))),
        _mm_shuffle_epi8(ehi, _mm_load_si128(reinterpret_cast<const __m128i*>(DIB1HI))));
    const __m128i b2 = _mm_or_si128(
        _mm_shuffle_epi8(elo, _mm_load_si128(reinterpret_cast<const __m128i*>(DIB2LO))),
        _mm_shuffle_epi8(ehi, _mm_load_si128(reinterpret_cast<const __m128i*>(DIB2HI))));
    // Ambiguous: any base has bit7 set (0xFF in signed = negative)
    const __m128i amb  = _mm_cmpgt_epi8(_mm_setzero_si128(),
                                         _mm_or_si128(_mm_or_si128(b0, b1), b2));
    // Codon index = b0*16 + b1*4 + b2. Mask to 0..3 first to neutralise 0xFF bases.
    // slli_epi16 on values 0..3 stays within each byte (max 3*16=48, 3*4=12 < 256).
    const __m128i b0s = _mm_and_si128(b0, _mm_set1_epi8(3));
    const __m128i b1s = _mm_and_si128(b1, _mm_set1_epi8(3));
    const __m128i b2s = _mm_and_si128(b2, _mm_set1_epi8(3));
    const __m128i idx = _mm_add_epi8(
        _mm_add_epi8(_mm_slli_epi16(b0s, 4), _mm_slli_epi16(b1s, 2)), b2s);
    const __m128i aa   = codon8(idx);
    const __m128i stop = _mm_cmpeq_epi8(aa, _mm_set1_epi8(0x7F));
    const __m128i brk  = _mm_or_si128(stop, amb);
    _mm_storel_epi64(reinterpret_cast<__m128i*>(out), aa);
    return (~static_cast<unsigned>(_mm_movemask_epi8(brk))) & 0xFFu;
}

// Dual forward+reverse-complement translate of 8 codons at p[0..23].
// Shares ONE enc16 pass (and the b0/b1/b2 deinterleave) between both strands:
//   forward codon index  = 16*b0 + 4*b1 + b2
//   RC codon index       = 63 - (16*b2 + 4*b1 + b0)   (BASE_COMP is x->3-x, i.e. ^3)
// out_fwd[0..7] = raw forward AA bytes (0..19 / 0x7F stop) — mask decides validity,
//   byte-identical to xlate8's output; forward caller path is unchanged.
// out_rc[0..7]  = RC AA bytes with break lanes (RC-stop OR ambiguous) set to 0xFF
//   (== AA_STOP sentinel) so the RC segment scan can read the array directly.
// Returns the forward valid mask (same contract as xlate8); RC breaks are in out_rc.
static inline unsigned xlate8_dual(const char* p, uint8_t* out_fwd, uint8_t* out_rc) noexcept {
    const __m128i elo = enc16(_mm_loadu_si128(reinterpret_cast<const __m128i*>(p)));
    const __m128i ehi = enc16(_mm_loadu_si128(reinterpret_cast<const __m128i*>(p + 8)));
    const __m128i b0 = _mm_or_si128(
        _mm_shuffle_epi8(elo, _mm_load_si128(reinterpret_cast<const __m128i*>(DIB0LO))),
        _mm_shuffle_epi8(ehi, _mm_load_si128(reinterpret_cast<const __m128i*>(DIB0HI))));
    const __m128i b1 = _mm_or_si128(
        _mm_shuffle_epi8(elo, _mm_load_si128(reinterpret_cast<const __m128i*>(DIB1LO))),
        _mm_shuffle_epi8(ehi, _mm_load_si128(reinterpret_cast<const __m128i*>(DIB1HI))));
    const __m128i b2 = _mm_or_si128(
        _mm_shuffle_epi8(elo, _mm_load_si128(reinterpret_cast<const __m128i*>(DIB2LO))),
        _mm_shuffle_epi8(ehi, _mm_load_si128(reinterpret_cast<const __m128i*>(DIB2HI))));
    const __m128i amb = _mm_cmpgt_epi8(_mm_setzero_si128(),
                                        _mm_or_si128(_mm_or_si128(b0, b1), b2));
    const __m128i b0s = _mm_and_si128(b0, _mm_set1_epi8(3));
    const __m128i b1s = _mm_and_si128(b1, _mm_set1_epi8(3));
    const __m128i b2s = _mm_and_si128(b2, _mm_set1_epi8(3));
    // Forward: idx = 16*b0 + 4*b1 + b2
    const __m128i idx = _mm_add_epi8(
        _mm_add_epi8(_mm_slli_epi16(b0s, 4), _mm_slli_epi16(b1s, 2)), b2s);
    const __m128i aa   = codon8(idx);
    const __m128i stop = _mm_cmpeq_epi8(aa, _mm_set1_epi8(0x7F));
    const __m128i brk  = _mm_or_si128(stop, amb);
    _mm_storel_epi64(reinterpret_cast<__m128i*>(out_fwd), aa);
    // RC: rev = 16*b2 + 4*b1 + b0 ; rc_idx = 63 - rev = 63 ^ rev (rev in 0..63)
    const __m128i rev = _mm_add_epi8(
        _mm_add_epi8(_mm_slli_epi16(b2s, 4), _mm_slli_epi16(b1s, 2)), b0s);
    const __m128i rcidx = _mm_xor_si128(rev, _mm_set1_epi8(63));
    const __m128i rcaa  = codon8(rcidx);
    const __m128i rcstop = _mm_cmpeq_epi8(rcaa, _mm_set1_epi8(0x7F));
    const __m128i rcbrk  = _mm_or_si128(rcstop, amb);
    // break lanes -> 0xFF sentinel; valid lanes -> AA (0..19)
    const __m128i rcout = _mm_blendv_epi8(rcaa, _mm_set1_epi8(-1), rcbrk);
    _mm_storel_epi64(reinterpret_cast<__m128i*>(out_rc), rcout);
    return (~static_cast<unsigned>(_mm_movemask_epi8(brk))) & 0xFFu;
}

} // namespace aamer_detail
#endif // __AVX2__

// ── 6-frame translation ───────────────────────────────────────────────────────

// For each inter-stop AA segment of length >= min_aa_len in all 6 reading frames,
// calls cb(frame, aa_ptr, aa_len, nt_start, nt_end) where [nt_start, nt_end) is the
// half-open interval on the FORWARD strand (reverse-strand ORFs are complemented back).
// Frame 0-2 = forward, 3-5 = reverse complement.
// aa_ptr is valid only for the duration of the call (points into a thread-local buffer).
// Ambiguous bases (N etc.) break the reading frame like a stop codon.
// Genetic code 11 (standard bacterial/archaeal).
template <typename Cb>
inline void translate_6frame(std::string_view seq, int min_aa_len, Cb&& cb) {
    const int n = static_cast<int>(seq.size());
    const char* s = seq.data();
    // Thread-local AA segment buffer: eliminates 6 per-call heap allocations.
    thread_local std::vector<uint8_t> seg;

#ifdef __AVX2__
    // Folded 6-frame: the forward SIMD pass ALSO derives the reverse-complement
    // codons from the SAME enc16/deinterleave (xlate8_dual), so 6 translation
    // passes collapse to 3. RC frame rf shares the triplet boundaries of forward
    // frame ff=(n-rf)%3 (length-dependent — never hardcoded), so RC frame rf's AA
    // array is exactly forward frame ff's codons walked BACKWARD. We materialise
    // rc_arr[rf] (0xFF = stop/ambiguous break sentinel) during the forward pass,
    // then emit RC segments after all forward frames — preserving the original
    // callback order (fwd 0,1,2 then rc 3,4,5).
    thread_local std::vector<uint8_t> rc_arr[3];

    // One segment buffer sized to the max codons in any frame (<= n/3), written
    // through a raw cursor: no per-AA capacity checks, no per-segment clear cost.
    const int cap = (n / 3) + 1;
    seg.resize(cap);
    uint8_t* const segp = seg.data();

    for (int frame = 0; frame < 3; ++frame) {
        uint8_t* cur = segp;                          // segment write cursor
        int seg_nt_start = frame;
        const int ncod = (n - frame) / 3;
        const int rf = ((n - frame) % 3 + 3) % 3;   // paired RC frame
        rc_arr[rf].resize(ncod > 0 ? ncod : 0);
        uint8_t* rp = ncod > 0 ? rc_arr[rf].data() : nullptr;
        const char* base = s + frame;
        int c = 0;

        uint8_t tmp[8], rtmp[8];
        // SIMD loop: 8 codons per iteration.
        // Loop guard: c+8 <= ncod ensures bytes base[c*3..c*3+23] are all in-bounds
        // (since (c+8)*3 <= ncod*3 <= n-frame, so frame+(c+8)*3 <= n).
        for (; c + 8 <= ncod; c += 8) {
            const unsigned vmask = aamer_detail::xlate8_dual(base + c * 3, tmp, rtmp);
            std::copy(rtmp, rtmp + 8, rp + c);        // stash RC codons for this frame
            if (vmask == 0xFFu) {
                // all 8 valid → single 8-byte store into the cursor (cur+8 <= cap)
                _mm_storel_epi64(reinterpret_cast<__m128i*>(cur),
                                 _mm_loadl_epi64(reinterpret_cast<const __m128i*>(tmp)));
                cur += 8;
                continue;
            }
            for (int j = 0; j < 8; ++j) {
                if (vmask & (1u << j)) {
                    *cur++ = tmp[j];
                } else {
                    const int nt = frame + (c + j) * 3;
                    const int len = (int)(cur - segp);
                    if (len >= min_aa_len)
                        cb(frame, segp, len, seg_nt_start, nt);
                    cur = segp;
                    seg_nt_start = nt + 3;
                }
            }
        }

        // Scalar tail (remainder after SIMD): emit forward AND fill rc_arr[rf].
        for (int i = frame + c * 3, cc = c; i + 2 < n; i += 3, ++cc) {
            const uint8_t b0 = DNA_ENC[(uint8_t)s[i]];
            const uint8_t b1 = DNA_ENC[(uint8_t)s[i+1]];
            const uint8_t b2 = DNA_ENC[(uint8_t)s[i+2]];
            if (b0 == 0xFF || b1 == 0xFF || b2 == 0xFF) {
                rp[cc] = 0xFF;                          // ambiguous breaks both strands
                const int len = (int)(cur - segp);
                if (len >= min_aa_len)
                    cb(frame, segp, len, seg_nt_start, i);
                cur = segp;
                seg_nt_start = i + 3;
                continue;
            }
            // RC codon = complement+reverse of forward bases: (3^b2,3^b1,3^b0).
            // CODON11 stores AA_STOP (0xFF) for stops, which is our break sentinel.
            rp[cc] = CODON11[(3 ^ b2) * 16 + (3 ^ b1) * 4 + (3 ^ b0)];
            const uint8_t aa = CODON11[b0 * 16 + b1 * 4 + b2];
            if (aa == AA_STOP) {
                const int len = (int)(cur - segp);
                if (len >= min_aa_len)
                    cb(frame, segp, len, seg_nt_start, i);
                cur = segp;
                seg_nt_start = i + 3;
            } else {
                *cur++ = aa;
            }
        }
        const int flen = (int)(cur - segp);
        if (flen >= min_aa_len)
            cb(frame, segp, flen, seg_nt_start, n);
    }

    // Reverse-complement frames 3-5, emitted from the stashed arrays. RC codons
    // emerge in reverse of forward order, so walk each frame's array BACKWARD.
    for (int rf = 0; rf < 3; ++rf) {
        const int ncod = (int)rc_arr[rf].size();
        if (ncod == 0) continue;
        const int ff = ((n - rf) % 3 + 3) % 3;         // producing forward frame
        const uint8_t* rp = rc_arr[rf].data();
        uint8_t* cur = segp;
        int seg_rc_end = n - rf;                        // = (n-1-rf) + 1
        for (int cc = ncod - 1; cc >= 0; --cc) {
            const int pos_start = ff + 3 * cc;          // fwd-strand start of triplet
            const uint8_t a = rp[cc];
            if (a == 0xFF) {                            // stop or ambiguous → break
                const int len = (int)(cur - segp);
                if (len >= min_aa_len)
                    cb(3 + rf, segp, len, n - seg_rc_end, n - pos_start);
                cur = segp;
                seg_rc_end = pos_start;
            } else {
                *cur++ = a;
            }
        }
        const int rlen = (int)(cur - segp);
        if (rlen >= min_aa_len)
            cb(3 + rf, segp, rlen, n - seg_rc_end, n);
    }
#else
    // Scalar fallback (no AVX2): original independent forward + RC passes.
    for (int frame = 0; frame < 3; ++frame) {
        seg.clear();
        int seg_nt_start = frame;
        for (int i = frame; i + 2 < n; i += 3) {
            const uint8_t b0 = DNA_ENC[(uint8_t)s[i]];
            const uint8_t b1 = DNA_ENC[(uint8_t)s[i+1]];
            const uint8_t b2 = DNA_ENC[(uint8_t)s[i+2]];
            if (b0 == 0xFF || b1 == 0xFF || b2 == 0xFF) {
                if ((int)seg.size() >= min_aa_len)
                    cb(frame, seg.data(), (int)seg.size(), seg_nt_start, i);
                seg.clear();
                seg_nt_start = i + 3;
                continue;
            }
            const uint8_t aa = CODON11[b0 * 16 + b1 * 4 + b2];
            if (aa == AA_STOP) {
                if ((int)seg.size() >= min_aa_len)
                    cb(frame, seg.data(), (int)seg.size(), seg_nt_start, i);
                seg.clear();
                seg_nt_start = i + 3;
            } else {
                seg.push_back(aa);
            }
        }
        if ((int)seg.size() >= min_aa_len)
            cb(frame, seg.data(), (int)seg.size(), seg_nt_start, n);
    }
    for (int frame = 0; frame < 3; ++frame) {
        seg.clear();
        const int rc_start = n - 1 - frame;
        int seg_rc_end = rc_start + 1;
        for (int pos = rc_start; pos - 2 >= 0; pos -= 3) {
            const uint8_t e0 = DNA_ENC[(uint8_t)s[pos]];
            const uint8_t e1 = DNA_ENC[(uint8_t)s[pos-1]];
            const uint8_t e2 = DNA_ENC[(uint8_t)s[pos-2]];
            if (e0 == 0xFF || e1 == 0xFF || e2 == 0xFF) {
                if ((int)seg.size() >= min_aa_len)
                    cb(3+frame, seg.data(), (int)seg.size(), n - seg_rc_end, n - (pos - 2));
                seg.clear();
                seg_rc_end = pos - 2;
                continue;
            }
            const uint8_t b0 = BASE_COMP[e0];
            const uint8_t b1 = BASE_COMP[e1];
            const uint8_t b2 = BASE_COMP[e2];
            const uint8_t aa = CODON11[b0 * 16 + b1 * 4 + b2];
            if (aa == AA_STOP) {
                if ((int)seg.size() >= min_aa_len)
                    cb(3+frame, seg.data(), (int)seg.size(), n - seg_rc_end, n - (pos - 2));
                seg.clear();
                seg_rc_end = pos - 2;
            } else {
                seg.push_back(aa);
            }
        }
        if ((int)seg.size() >= min_aa_len)
            cb(3+frame, seg.data(), (int)seg.size(), n - seg_rc_end, n);
    }
#endif
}

// Forward-only 3-frame translation (frames 0-2). Identical to the forward half of
// translate_6frame, but skips the reverse-complement frames 3-5 entirely.
// Used by marker scoring, which only needs forward ORFs and pays no RC overhead.
template <typename Cb>
inline void translate_3frame_fwd(std::string_view seq, int min_aa_len, Cb&& cb) {
    const int n = static_cast<int>(seq.size());
    const char* s = seq.data();
    thread_local std::vector<uint8_t> seg;

    // Cursor-based segment buffer (see translate_6frame): one resize, no per-AA
    // capacity checks, single 8-byte store on the all-valid fast path.
    const int cap = (n / 3) + 1;
    seg.resize(cap);
    uint8_t* const segp = seg.data();

    for (int frame = 0; frame < 3; ++frame) {
        uint8_t* cur = segp;
        int seg_nt_start = frame;
        const int ncod = (n - frame) / 3;
        const char* base = s + frame;
        int c = 0;

#ifdef __AVX2__
        uint8_t tmp[8];
        for (; c + 8 <= ncod; c += 8) {
            const unsigned vmask = aamer_detail::xlate8(base + c * 3, tmp);
            if (vmask == 0xFFu) {
                _mm_storel_epi64(reinterpret_cast<__m128i*>(cur),
                                 _mm_loadl_epi64(reinterpret_cast<const __m128i*>(tmp)));
                cur += 8;
                continue;
            }
            for (int j = 0; j < 8; ++j) {
                if (vmask & (1u << j)) {
                    *cur++ = tmp[j];
                } else {
                    const int nt = frame + (c + j) * 3;
                    const int len = (int)(cur - segp);
                    if (len >= min_aa_len)
                        cb(frame, segp, len, seg_nt_start, nt);
                    cur = segp;
                    seg_nt_start = nt + 3;
                }
            }
        }
#endif

        for (int i = frame + c * 3; i + 2 < n; i += 3) {
            const uint8_t b0 = DNA_ENC[(uint8_t)s[i]];
            const uint8_t b1 = DNA_ENC[(uint8_t)s[i+1]];
            const uint8_t b2 = DNA_ENC[(uint8_t)s[i+2]];
            if (b0 == 0xFF || b1 == 0xFF || b2 == 0xFF) {
                const int len = (int)(cur - segp);
                if (len >= min_aa_len)
                    cb(frame, segp, len, seg_nt_start, i);
                cur = segp;
                seg_nt_start = i + 3;
                continue;
            }
            const uint8_t aa = CODON11[b0 * 16 + b1 * 4 + b2];
            if (aa == AA_STOP) {
                const int len = (int)(cur - segp);
                if (len >= min_aa_len)
                    cb(frame, segp, len, seg_nt_start, i);
                cur = segp;
                seg_nt_start = i + 3;
            } else {
                *cur++ = aa;
            }
        }
        const int flen = (int)(cur - segp);
        if (flen >= min_aa_len)
            cb(frame, segp, flen, seg_nt_start, n);
    }
}

// ── Convenience extractors ────────────────────────────────────────────────────

// Extract sorted, deduplicated aamer hashes from a DNA sequence (6-frame translation).
// min_seg_aa: minimum inter-stop segment length to emit k-mers from (default = k).
[[nodiscard]] inline std::vector<uint64_t>
extract_aamers_dna(std::string_view seq, int k = AAMER_K, int min_seg_aa = AAMER_K) {
    std::vector<uint64_t> hashes;
    hashes.reserve(seq.size() / 4);
    translate_6frame(seq, min_seg_aa, [&](int, const uint8_t* seg, int len, int, int) {
        emit_aamers(seg, len, k, hashes);
    });
    std::sort(hashes.begin(), hashes.end());
    hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());
    return hashes;
}

// FracMinHash variant: keep only hashes ≤ max_hash during extraction.
[[nodiscard]] inline std::vector<uint64_t>
extract_aamers_dna(std::string_view seq, int k, int min_seg_aa, uint64_t max_hash) {
    std::vector<uint64_t> hashes;
    hashes.reserve(seq.size() / (4 * 32));
    translate_6frame(seq, min_seg_aa, [&](int, const uint8_t* seg, int len, int, int) {
        emit_aamers_frac(seg, len, k, max_hash, hashes);
    });
    std::sort(hashes.begin(), hashes.end());
    hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());
    return hashes;
}

// In-place append variant: appends unsorted hashes to an existing vector (caller sorts/deduplicates).
// Use with a thread_local buffer to avoid per-call allocation.
inline void
extract_aamers_dna_into(std::string_view seq, int k, int min_seg_aa, uint64_t max_hash,
                          std::vector<uint64_t>& out) {
    translate_6frame(seq, min_seg_aa, [&](int, const uint8_t* seg, int len, int, int) {
        emit_aamers_frac(seg, len, k, max_hash, out);
    });
}

// Murphy10 variant of extract_aamers_aa: same extraction but uses aamer_hash_murphy10.
// Used for building and querying reduced-alphabet marker pools.
[[nodiscard]] inline std::vector<uint64_t>
extract_aamers_aa_murphy10(std::string_view protein, int k = AAMER_K) {
    std::vector<uint64_t> hashes;
    hashes.reserve(protein.size());
    std::vector<uint8_t> seg;
    seg.reserve(protein.size());
    auto flush = [&] {
        if ((int)seg.size() >= k) {
            for (int i = 0; i + k <= (int)seg.size(); ++i) {
                if (!aamer_is_low_complexity(seg.data() + i, k))
                    hashes.push_back(aamer_hash_murphy10(seg.data() + i));
            }
        }
        seg.clear();
    };
    for (unsigned char c : protein) {
        if (c == '-') continue;
        const uint8_t aa = AA_ENC[c];
        if (aa == AA_STOP) flush();
        else seg.push_back(aa);
    }
    flush();
    std::sort(hashes.begin(), hashes.end());
    hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());
    return hashes;
}

// Extract sorted, deduplicated aamer hashes from a protein sequence.
// Handles uppercase/lowercase, '-' gaps (skip), '*' stops, 'X'/'B'/'Z' (break).
// Used for building the marker panel from GTDB-Tk MSA sequences.
[[nodiscard]] inline std::vector<uint64_t>
extract_aamers_aa(std::string_view protein, int k = AAMER_K) {
    std::vector<uint64_t> hashes;
    hashes.reserve(protein.size());

    std::vector<uint8_t> seg;
    seg.reserve(protein.size());

    auto flush = [&] {
        if ((int)seg.size() >= k) {
            emit_aamers(seg.data(), (int)seg.size(), k, hashes);
        }
        seg.clear();
    };

    for (unsigned char c : protein) {
        if (c == '-') continue;          // MSA gap — skip without breaking segment
        const uint8_t aa = AA_ENC[c];
        if (aa == AA_STOP) {
            flush();
        } else {
            seg.push_back(aa);
        }
    }
    flush();

    std::sort(hashes.begin(), hashes.end());
    hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());
    return hashes;
}

// ── Dayhoff-6 reduced alphabet ────────────────────────────────────────────────
// Groups (index into AA_ENC order ACDEFGHIKLMNPQRSTVWY → 0..19):
//   0=C  1=AGPST  2=DENQ  3=HKR  4=ILMV  5=FWY
static constexpr uint8_t AA_DAYHOFF6[20] = {
    1, 0, 2, 2, 5, 1, 3, 4, 3, 4, 4, 2, 1, 2, 3, 1, 1, 4, 5, 5
};

static constexpr int AAMER_K_D6 = 12; // Dayhoff-6 k-mer length
static constexpr int AAMER_S_D6 =  4; // syncmer inner s-mer length

// Dayhoff-6 k=12 hash: pack 12 groups × 3 bits = 36 bits → splitmix64.
// d6[0..11] must be AA_DAYHOFF6-encoded (values 0..5).
[[nodiscard]] inline uint64_t aamer_hash_d6(const uint8_t* d6) noexcept {
    uint64_t v = 0;
    for (int i = 0; i < AAMER_K_D6; ++i)
        v |= static_cast<uint64_t>(d6[i]) << (3 * i);
    v ^= v >> 30;
    v *= 0xbf58476d1ce4e5b9ULL;
    v ^= v >> 27;
    v *= 0x94d049bb133111ebULL;
    v ^= v >> 31;
    return v;
}

// Closed syncmer (t=0): select k-mer if the minimum s-mer (packed 3×s bits, no mixing)
// is at position 0. Gives ~1/(k-s+1) = ~11% selection at k=12, s=4.
// Identical criterion at pool-build and query time → consistent per-ORF sketch density.
[[nodiscard]] inline bool aamer_is_syncmer_d6(const uint8_t* d6) noexcept {
    // Pack each s-mer as s×3-bit integer (12 bits for s=4) — no mixing, just ordering.
    auto smer = [&](int pos) -> uint16_t {
        uint16_t v = 0;
        for (int i = 0; i < AAMER_S_D6; ++i)
            v |= static_cast<uint16_t>(d6[pos + i]) << (3 * i);
        return v;
    };
    const uint16_t s0 = smer(0);
    for (int i = 1; i <= AAMER_K_D6 - AAMER_S_D6; ++i)
        if (smer(i) < s0) return false;
    return true;
}

// Append Dayhoff-6 k=12 syncmer hashes from a DNA sequence (6-frame translation).
// Only k-mers passing the closed-syncmer criterion (t=0) are emitted.
// Low-complexity k-mers (run ≥ (k+1)/2 of same Dayhoff group) are skipped.
// No FracMinHash subsampling: syncmers provide consistent ~11% density already.
inline void extract_d6_dna_into(std::string_view seq, int min_seg_aa,
                                std::vector<uint64_t>& out) {
    constexpr int k = AAMER_K_D6;
    translate_6frame(seq, min_seg_aa, [&](int, const uint8_t* seg, int len, int, int) {
        for (int i = 0; i + k <= len; ++i) {
            // Recode to Dayhoff-6 and check for invalid AA (0xFF from AA_ENC).
            uint8_t d6[AAMER_K_D6];
            bool valid = true;
            for (int j = 0; j < k; ++j) {
                if (seg[i + j] >= 20) { valid = false; break; }
                d6[j] = AA_DAYHOFF6[seg[i + j]];
            }
            if (!valid) continue;
            if (!aamer_is_syncmer_d6(d6)) continue;
            if (aamer_is_low_complexity(d6, k)) continue;
            out.push_back(aamer_hash_d6(d6));
        }
    });
}

// Extract Dayhoff-6 k=12 syncmer hashes from one pre-translated ORF segment using a
// sliding-window minimum deque — O(len) vs the per-position O(WIN) smer recomputation.
// seg: AA_ENC-encoded values 0..19 (translate_6frame guarantees this).
// Appends to `out` without clearing.
inline void extract_d6_orf_syncmers(const uint8_t* seg, int len,
                                     std::vector<uint64_t>& out) {
    constexpr int K = AAMER_K_D6;
    constexpr int S = AAMER_S_D6;
    constexpr int W = K - S + 1;  // 9: number of s-mers per k-mer window
    if (len < K) return;

    thread_local std::vector<uint8_t>  d6;
    thread_local std::vector<uint16_t> sv;
    d6.resize(len);
    for (int i = 0; i < len; ++i) {
        if (seg[i] >= 20) return;  // ambiguous AA — skip ORF
        d6[i] = AA_DAYHOFF6[seg[i]];
    }
    const int nsm = len - S + 1;
    sv.resize(nsm);
    for (int i = 0; i < nsm; ++i)
        sv[i] = static_cast<uint16_t>(d6[i] | (d6[i+1]<<3) | (d6[i+2]<<6) | (d6[i+3]<<9));

    // Monotonic deque for sliding-window minimum (window = [i, i+W-1]).
    // Fixed-capacity 16 (> W=9). Stores s-mer indices; back >= sv comparison keeps
    // invariant that sv[dq[head]] is the minimum in the current window.
    struct Deque {
        int buf[16]; int h = 0, t = 0;
        bool empty()     const { return h == t; }
        int  front()     const { return buf[h & 15]; }
        int  back()      const { return buf[(t-1) & 15]; }
        void push(int v)       { buf[t++ & 15] = v; }
        void pop_front()       { ++h; }
        void pop_back()        { --t; }
    } dq;

    // Seed deque with first window [0, W-1].
    for (int j = 0; j < W && j < nsm; ++j) {
        while (!dq.empty() && sv[dq.back()] >= sv[j]) dq.pop_back();
        dq.push(j);
    }

    const int nk = len - K + 1;
    for (int i = 0; i < nk; ++i) {
        // Expire front if it left the window.
        while (!dq.empty() && dq.front() < i) dq.pop_front();

        // Closed syncmer t=0: emit if window minimum is at left edge i.
        if (!dq.empty() && dq.front() == i)
            if (!aamer_is_low_complexity(d6.data() + i, K))
                out.push_back(aamer_hash_d6(d6.data() + i));

        // Enqueue next s-mer (right edge of next window).
        const int nx = i + W;
        if (nx < nsm) {
            while (!dq.empty() && sv[dq.back()] >= sv[nx]) dq.pop_back();
            dq.push(nx);
        }
    }
}

// Dayhoff-6 variant of extract_aamers_aa: for pool building from MSA protein strings.
// Gaps ('-') skip without breaking the segment (same as extract_aamers_aa).
// Only emits syncmer k-mers to keep pool sparse and density-consistent.
[[nodiscard]] inline std::vector<uint64_t>
extract_aamers_d6_aa(std::string_view protein) {
    std::vector<uint64_t> hashes;
    hashes.reserve(protein.size() / 8);
    constexpr int k = AAMER_K_D6;

    std::vector<uint8_t> seg; // Dayhoff-6 encoded segment
    seg.reserve(protein.size());

    auto flush = [&] {
        for (int i = 0; i + k <= (int)seg.size(); ++i) {
            if (!aamer_is_syncmer_d6(seg.data() + i)) continue;
            if (aamer_is_low_complexity(seg.data() + i, k)) continue;
            hashes.push_back(aamer_hash_d6(seg.data() + i));
        }
        seg.clear();
    };

    for (unsigned char c : protein) {
        if (c == '-') continue;
        const uint8_t aa = AA_ENC[c];
        if (aa == AA_STOP || aa == 0xFF) { flush(); continue; }
        seg.push_back(AA_DAYHOFF6[aa]);
    }
    flush();

    std::sort(hashes.begin(), hashes.end());
    hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());
    return hashes;
}

} // namespace genopack
