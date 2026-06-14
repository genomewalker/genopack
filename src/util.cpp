#include <genopack/util.hpp>
#include <Eigen/Eigenvalues>
#include <zlib.h>
#ifdef GENOPACK_HAVE_LIBDEFLATE
#  include <libdeflate.h>
#endif
#include <thread>
#include <algorithm>
#include <array>
#include <chrono>
#include <numeric>
#include <cmath>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <set>
#include <unistd.h>
#include <sstream>
#include <stdexcept>

namespace genopack {

#ifdef GENOPACK_HAVE_LIBDEFLATE
// Thread-local libdeflate decompressor — reused across calls, avoids alloc/free
// per genome. libdeflate_decompressor is stateless between calls and thread-safe
// when each thread owns its own instance.
namespace {
struct TlDecompressor {
    libdeflate_decompressor* d;
    TlDecompressor() : d(libdeflate_alloc_decompressor()) {}
    ~TlDecompressor() { if (d) libdeflate_free_decompressor(d); }
};
static thread_local TlDecompressor tl_decomp;
} // namespace
#endif

static std::string read_file_bytes(const std::filesystem::path& path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) throw std::runtime_error("Cannot open: " + path.string());
    size_t sz = static_cast<size_t>(f.tellg());
    f.seekg(0);
    std::string buf(sz, '\0');
    f.read(buf.data(), static_cast<std::streamsize>(sz));
    if (!f) throw std::runtime_error("Read error: " + path.string());
    return buf;
}

// Read from an already-open fd (closes fd before returning).
static std::string read_fd_bytes(int fd, const std::filesystem::path& path) {
    off_t sz = ::lseek(fd, 0, SEEK_END);
    if (sz < 0) { ::close(fd); throw std::runtime_error("lseek failed: " + path.string()); }
    ::lseek(fd, 0, SEEK_SET);
    std::string buf(static_cast<size_t>(sz), '\0');
    size_t total = 0;
    while (total < static_cast<size_t>(sz)) {
        ssize_t n = ::read(fd, buf.data() + total, static_cast<size_t>(sz) - total);
        if (n < 0) {
            if (errno == EINTR) continue; // retry on signal interrupt
            ::close(fd); throw std::runtime_error("Read error: " + path.string());
        }
        if (n == 0) break; // NFS may report stale inode size; treat early EOF as end of data
        total += static_cast<size_t>(n);
    }
    ::close(fd);
    buf.resize(total);
    return buf;
}

static bool is_gzip(const std::string& buf) {
    return buf.size() >= 2 &&
           static_cast<uint8_t>(buf[0]) == 0x1f &&
           static_cast<uint8_t>(buf[1]) == 0x8b;
}

// Shared decompression core — operates on already-loaded compressed bytes.
static std::string decompress_gz_buf(const std::string& compressed,
                                     const std::filesystem::path& path) {
    if (!is_gzip(compressed))
        return compressed; // plain FASTA

#ifdef GENOPACK_HAVE_LIBDEFLATE
    // Fast path: libdeflate single-shot decompression.
    // Get uncompressed size from gzip ISIZE trailer (last 4 bytes, little-endian).
    // ISIZE is mod 2^32; safe for files < 4 GB.
    if (compressed.size() < 18)
        throw std::runtime_error("decompress_gz: truncated gzip: " + path.string());
    uint32_t isize = 0;
    std::memcpy(&isize, compressed.data() + compressed.size() - 4, 4);

    // Allocate output (add 1 MB slack for edge cases near 2^32 boundary)
    std::string out(static_cast<size_t>(isize) + (1u << 20), '\0');

    struct libdeflate_decompressor* d = tl_decomp.d;
    if (!d) throw std::runtime_error("libdeflate_alloc_decompressor failed");
    size_t actual_out = 0;
    auto result = libdeflate_gzip_decompress(
        d,
        compressed.data(), compressed.size(),
        out.data(), out.size(),
        &actual_out);

    if (result == LIBDEFLATE_SUCCESS) {
        out.resize(actual_out);
        return out;
    }
    // Fall through to zlib on failure (multi-member gzip, etc.)
#endif

    // zlib fallback
    {
        z_stream zs{};
        if (inflateInit2(&zs, 15 + 32) != Z_OK) // +32 = auto-detect gzip/zlib
            throw std::runtime_error("inflateInit2 failed");
        zs.next_in  = reinterpret_cast<Bytef*>(const_cast<char*>(compressed.data()));
        zs.avail_in = static_cast<uInt>(compressed.size());
        std::string zlib_out;
        char buf[1 << 20];
        int ret;
        do {
            zs.next_out  = reinterpret_cast<Bytef*>(buf);
            zs.avail_out = sizeof(buf);
            ret = inflate(&zs, Z_NO_FLUSH);
            if (ret != Z_OK && ret != Z_STREAM_END)
                throw std::runtime_error("inflate error: " + path.string());
            zlib_out.append(buf, sizeof(buf) - zs.avail_out);
        } while (ret != Z_STREAM_END);
        inflateEnd(&zs);
        return zlib_out;
    }
}

std::string decompress_gz(const std::filesystem::path& path) {
    return decompress_gz_buf(read_file_bytes(path), path);
}

std::string decompress_gz_fd(int fd, const std::filesystem::path& path) {
    return decompress_gz_buf(read_fd_bytes(fd, path), path);
}

uint64_t genome_minhash(const std::string& fasta, int k) {
    uint64_t min_hash = UINT64_MAX;
    uint64_t kmer = 0;
    int kmer_len = 0;
    for (char c : fasta) {
        if (c == '>' || c == '\n' || c == '\r') { kmer_len = 0; kmer = 0; continue; }
        char u = static_cast<char>(c & ~0x20);
        uint8_t base;
        if      (u == 'A') base = 0;
        else if (u == 'C') base = 1;
        else if (u == 'G') base = 2;
        else if (u == 'T') base = 3;
        else { kmer_len = 0; kmer = 0; continue; }
        kmer = (kmer << 2) | base;
        if (++kmer_len >= k) {
            uint64_t v = kmer & ((1ULL << (2 * k)) - 1);
            // MurmurHash3 finalizer mix
            v ^= v >> 33;
            v *= 0xff51afd7ed558ccdULL;
            v ^= v >> 33;
            v *= 0xc4ceb9fe1a85ec53ULL;
            v ^= v >> 33;
            if (v < min_hash) min_hash = v;
        }
    }
    return min_hash;
}

// Build canonical k=4 index at static init time.
// Base encoding: A=0, C=1, G=2, T=3
// 4-mer encoding: b0*64 + b1*16 + b2*4 + b3  (b0 = leftmost)
// revcomp4(x) using complement table comp[] = {3,2,1,0}
static const std::array<uint8_t, 136>& canonical_kmer4_ids() {
    static std::array<uint8_t, 136> ids = []() {
        const uint8_t comp[4] = {3, 2, 1, 0};
        std::array<uint8_t, 136> result{};
        int n = 0;
        for (int i = 0; i < 256; ++i) {
            // revcomp of encoded 4-mer i
            uint8_t rc = static_cast<uint8_t>(
                (comp[i & 3] << 6) | (comp[(i >> 2) & 3] << 4) |
                (comp[(i >> 4) & 3] << 2) | comp[(i >> 6) & 3]);
            if (i <= rc) {
                result[n++] = static_cast<uint8_t>(i);
            }
        }
        // n must equal 136
        return result;
    }();
    return ids;
}

// Map each of 256 4-mer encodings to its index in the 136-dim canonical array.
// Returns 255 for any encoding not at the canonical (min) end.
static const std::array<uint8_t, 256>& kmer4_canon_index() {
    static std::array<uint8_t, 256> idx = []() {
        const uint8_t comp[4] = {3, 2, 1, 0};
        const auto& ids = canonical_kmer4_ids();
        std::array<uint8_t, 256> result;
        result.fill(255);
        for (int j = 0; j < 136; ++j) {
            uint8_t i  = ids[j];
            uint8_t rc = static_cast<uint8_t>(
                (comp[i & 3] << 6) | (comp[(i >> 2) & 3] << 4) |
                (comp[(i >> 4) & 3] << 2) | comp[(i >> 6) & 3]);
            result[i]  = static_cast<uint8_t>(j);
            result[rc] = static_cast<uint8_t>(j); // palindromes map to same index
        }
        return result;
    }();
    return idx;
}

// Single 256-entry LUT encoding all per-character info in one byte:
//   bits 1:0 = base (0=A,1=C,2=G,3=T)
//   bit  2   = valid ACGT base
//   bit  3   = is GC
//   bit  4   = is newline (\n or \r)
//   bit  5   = is header start (>)
static const std::array<uint8_t, 256>& char_lut() {
    static std::array<uint8_t, 256> lut = []() {
        std::array<uint8_t, 256> t{};
        t['\n'] = 0x10; t['\r'] = 0x10;
        t['>']  = 0x20;
        // lowercase and uppercase ACGT
        for (unsigned char c : {'A','a'}) t[c] = 0x04 | 0; // valid, AT, base=0
        for (unsigned char c : {'C','c'}) t[c] = 0x04 | 0x08 | 1; // valid, GC, base=1
        for (unsigned char c : {'G','g'}) t[c] = 0x04 | 0x08 | 2; // valid, GC, base=2
        for (unsigned char c : {'T','t'}) t[c] = 0x04 | 3; // valid, AT, base=3
        return t;
    }();
    return lut;
}

FastaStats compute_fasta_stats(const std::string& fasta, int k) {
    FastaStats s;

    // Single pass: GC/length/n_contigs + k=4 canonical profile + OPH minhash
    const auto& lut       = char_lut();
    const auto& canon_idx = kmer4_canon_index();

    std::array<uint32_t, 136> counts{};
    uint64_t gc  = 0;
    uint64_t at  = 0;
    uint8_t  enc = 0;
    int      len = 0;
    bool     in_header = false;
    int64_t  gc_skew_cum = 0, gc_skew_max = 0;

    // OPH minhash inline (k passed in, default 21)
    uint64_t min_hash  = UINT64_MAX;
    uint64_t oph_kmer  = 0;
    int      oph_len   = 0;
    const uint64_t oph_mask = (1ULL << (2 * k)) - 1;

    for (unsigned char c : fasta) {
        uint8_t info = lut[c];

        if (info & 0x20) {          // '>'
            ++s.n_contigs;
            in_header = true;
            enc = 0; len = 0;
            oph_kmer = 0; oph_len = 0;
            continue;
        }
        if (in_header) {
            if (info & 0x10) in_header = false;
            continue;
        }
        if (info & 0x10) continue;  // newline

        if (!(info & 0x04)) {       // non-ACGT (N, ambiguous)
            enc = 0; len = 0;
            oph_kmer = 0; oph_len = 0;
            continue;
        }

        ++s.genome_length;
        uint8_t base = info & 0x03;
        if (info & 0x08) ++gc; else ++at;
        if (base == 2) ++gc_skew_cum;
        else if (base == 1) --gc_skew_cum;
        { int64_t a = gc_skew_cum < 0 ? -gc_skew_cum : gc_skew_cum; if (a > gc_skew_max) gc_skew_max = a; }

        // k=4 canonical profile
        enc = static_cast<uint8_t>((enc << 2) | base);
        if (++len >= 4)
            ++counts[canon_idx[enc]];

        // OPH minhash
        oph_kmer = ((oph_kmer << 2) | base) & oph_mask;
        if (++oph_len >= k) {
            uint64_t v = oph_kmer;
            v ^= v >> 33; v *= 0xff51afd7ed558ccdULL;
            v ^= v >> 33; v *= 0xc4ceb9fe1a85ec53ULL;
            v ^= v >> 33;
            if (v < min_hash) min_hash = v;
        }
    }

    s.oph_fingerprint = min_hash;
    uint64_t total = gc + at;
    s.gc_pct_x100 = total > 0 ? static_cast<uint16_t>(gc * 10000 / total) : 0;
    if (s.genome_length >= 10000 && gc_skew_max > 0) {
        int64_t abs_final = gc_skew_cum < 0 ? -gc_skew_cum : gc_skew_cum;
        s.gc_skew_closure = 1.0f - static_cast<float>(abs_final) / static_cast<float>(gc_skew_max);
    }

    // L2-normalise k=4 profile
    uint64_t total_kmers = 0;
    for (uint32_t ct : counts) total_kmers += ct;

    if (total_kmers > 0) {
        float inv = 1.0f / static_cast<float>(total_kmers);
        float sumsq = 0.0f;
        for (int i = 0; i < 136; ++i) {
            float v = static_cast<float>(counts[i]) * inv;
            s.kmer4_profile[i] = v;
            sumsq += v * v;
        }
        if (sumsq > 0.0f) {
            float norm = 1.0f / std::sqrt(sumsq);
            for (int i = 0; i < 136; ++i)
                s.kmer4_profile[i] *= norm;
        }
    }

    return s;
}

float compute_completeness_post_decontam(
    const std::string& fasta,
    const float*       tnf_centroid,
    float              threshold)
{
    if (!tnf_centroid) return 1.0f;

    float centroid_norm2 = 0.0f;
    for (int i = 0; i < 136; ++i) centroid_norm2 += tnf_centroid[i] * tnf_centroid[i];
    const float centroid_norm = std::sqrt(centroid_norm2);
    if (centroid_norm < 1e-8f) return 1.0f;

    const auto& lut       = char_lut();
    const auto& canon_idx = kmer4_canon_index();

    uint64_t clean_len = 0, total_len = 0;
    const char* p   = fasta.data();
    const char* end = p + fasta.size();

    while (p < end) {
        while (p < end && *p != '>') ++p;
        if (p >= end) break;
        ++p;
        while (p < end && *p != '\n') ++p;
        if (p < end) ++p;

        std::array<uint32_t, 136> counts{};
        uint64_t seq_len = 0;
        uint8_t enc = 0;
        int klen = 0;

        while (p < end && *p != '>') {
            uint8_t info = lut[static_cast<uint8_t>(*p++)];
            if (info & 0x10) continue;
            if (!(info & 0x04)) { enc = 0; klen = 0; continue; }
            ++seq_len;
            enc = static_cast<uint8_t>((enc << 2) | (info & 0x03));
            if (++klen >= 4) ++counts[canon_idx[enc]];
        }

        total_len += seq_len;

        if (seq_len < 5000) { clean_len += seq_len; continue; }

        uint64_t total_kmers = 0;
        for (uint32_t ct : counts) total_kmers += ct;
        if (total_kmers == 0) { clean_len += seq_len; continue; }

        float profile[136];
        float inv = 1.0f / static_cast<float>(total_kmers);
        float sumsq = 0.0f;
        for (int i = 0; i < 136; ++i) { profile[i] = static_cast<float>(counts[i]) * inv; sumsq += profile[i] * profile[i]; }
        if (sumsq > 0.0f) {
            float norm = 1.0f / std::sqrt(sumsq);
            for (int i = 0; i < 136; ++i) profile[i] *= norm;
        }

        float d2 = 0.0f;
        for (int i = 0; i < 136; ++i) { float diff = profile[i] - tnf_centroid[i]; d2 += diff * diff; }
        if (std::sqrt(d2) / centroid_norm < threshold)
            clean_len += seq_len;
    }

    return total_len > 0 ? static_cast<float>(clean_len) / static_cast<float>(total_len) : NAN;
}

float compute_self_coherence(
    const std::string& fasta,
    int   min_contig_len,
    int   min_contigs,
    float mad_multiplier)
{
    const auto& lut       = char_lut();
    const auto& canon_idx = kmer4_canon_index();

    struct ContigProf { float p[136]; uint64_t bp; };
    std::vector<ContigProf> profs;

    const char* ptr = fasta.data();
    const char* end = ptr + fasta.size();
    while (ptr < end) {
        while (ptr < end && *ptr != '>') ++ptr;
        if (ptr >= end) break;
        ++ptr;
        while (ptr < end && *ptr != '\n') ++ptr;
        if (ptr < end) ++ptr;

        std::array<uint32_t, 136> counts{};
        uint64_t seq_len = 0;
        uint8_t enc = 0; int klen = 0;
        while (ptr < end && *ptr != '>') {
            uint8_t info = lut[static_cast<uint8_t>(*ptr++)];
            if (info & 0x10) continue;
            if (!(info & 0x04)) { enc = 0; klen = 0; continue; }
            ++seq_len;
            enc = static_cast<uint8_t>((enc << 2) | (info & 0x03));
            if (++klen >= 4) ++counts[canon_idx[enc]];
        }
        if (static_cast<int>(seq_len) < min_contig_len) continue;

        uint64_t total_kmers = 0;
        for (uint32_t ct : counts) total_kmers += ct;
        if (total_kmers == 0) continue;

        ContigProf cp; cp.bp = seq_len;
        float inv = 1.0f / static_cast<float>(total_kmers), sumsq = 0.0f;
        for (int i = 0; i < 136; ++i) { cp.p[i] = static_cast<float>(counts[i]) * inv; sumsq += cp.p[i] * cp.p[i]; }
        if (sumsq > 0.0f) { float n = 1.0f / std::sqrt(sumsq); for (int i = 0; i < 136; ++i) cp.p[i] *= n; }
        profs.push_back(cp);
    }

    if (static_cast<int>(profs.size()) < min_contigs) return NAN;

    // Centroid = mean of L2-normalised profiles, then L2-normalise again.
    float centroid[136] = {};
    for (const auto& cp : profs)
        for (int i = 0; i < 136; ++i) centroid[i] += cp.p[i];
    {
        float n = static_cast<float>(profs.size()), sumsq = 0.0f;
        for (int i = 0; i < 136; ++i) { centroid[i] /= n; sumsq += centroid[i] * centroid[i]; }
        if (sumsq > 0.0f) { float nrm = 1.0f / std::sqrt(sumsq); for (int i = 0; i < 136; ++i) centroid[i] *= nrm; }
    }

    // Per-contig distances to self-centroid.
    std::vector<float> dists(profs.size());
    for (size_t i = 0; i < profs.size(); ++i) {
        float d2 = 0.0f;
        for (int j = 0; j < 136; ++j) { float df = profs[i].p[j] - centroid[j]; d2 += df * df; }
        dists[i] = std::sqrt(d2);
    }

    // Robust threshold: median + mad_multiplier * 1.4826 * MAD.
    std::vector<float> sd = dists;
    std::sort(sd.begin(), sd.end());
    const size_t h = sd.size() / 2;
    const float  med = (sd.size() % 2 == 0) ? 0.5f * (sd[h-1] + sd[h]) : sd[h];

    std::vector<float> ad(profs.size());
    for (size_t i = 0; i < profs.size(); ++i) ad[i] = std::fabs(dists[i] - med);
    std::sort(ad.begin(), ad.end());
    const float mad = (ad.size() % 2 == 0) ? 0.5f * (ad[h-1] + ad[h]) : ad[h];
    const float threshold = med + mad_multiplier * 1.4826f * mad;

    uint64_t clean_bp = 0, total_bp = 0;
    for (size_t i = 0; i < profs.size(); ++i) {
        total_bp += profs[i].bp;
        if (dists[i] <= threshold) clean_bp += profs[i].bp;
    }
    return total_bp > 0 ? static_cast<float>(clean_bp) / static_cast<float>(total_bp) : NAN;
}

// ── Chargaff parity + transition-operator signals ─────────────────────────────

float compute_chargaff_parity(const std::string& fasta, int min_contig_len) {
    const auto& lut = char_lut();
    uint64_t counts[256] = {};

    const char* p = fasta.data(), *end = p + fasta.size();
    while (p < end) {
        while (p < end && *p != '>') ++p;
        if (p >= end) break;
        ++p; while (p < end && *p != '\n') ++p; if (p < end) ++p;

        uint32_t local[256] = {};
        uint8_t enc = 0; int klen = 0; uint64_t seq_len = 0;
        while (p < end && *p != '>') {
            uint8_t info = lut[static_cast<uint8_t>(*p++)];
            if (info & 0x10) continue;
            if (!(info & 0x04)) { enc = 0; klen = 0; continue; }
            ++seq_len;
            enc = static_cast<uint8_t>((enc << 2) | (info & 0x03));
            if (++klen >= 4) ++local[enc];
        }
        if (static_cast<int>(seq_len) >= min_contig_len)
            for (int i = 0; i < 256; ++i) counts[i] += local[i];
    }

    // For k=4 encoding enc=(b3<<6)|(b2<<4)|(b1<<2)|b0 (b3=first base added)
    // RC: complement+reverse = ((3-b0)<<6)|((3-b1)<<4)|((3-b2)<<2)|(3-b3)
    double sum = 0.0; int n = 0;
    for (int i = 0; i < 256; ++i) {
        uint8_t b0=i&3, b1=(i>>2)&3, b2=(i>>4)&3, b3=(i>>6)&3;
        int rc = ((3-b0)<<6)|((3-b1)<<4)|((3-b2)<<2)|(3-b3);
        if (rc <= i) continue;  // skip palindromes (rc==i) and already-seen pairs (rc<i)
        uint64_t ci = counts[i], crc = counts[static_cast<uint8_t>(rc)];
        uint64_t tot = ci + crc;
        if (tot == 0) continue;
        sum += 1.0 - static_cast<double>(ci > crc ? ci - crc : crc - ci) / static_cast<double>(tot);
        ++n;
    }
    return n > 0 ? static_cast<float>(sum / n) : NAN;
}

namespace {

// Accumulate k=3 transition counts from FASTA-region [p, end).
// T[64][64] is updated in-place (not normalised here).
static void accum_k3_transitions(const char* p, const char* end,
                                  double T[64][64],
                                  const std::array<uint8_t,256>& lut)
{
    constexpr int K = 3, N = 64;
    uint32_t state = 0; int klen = 0;
    while (p < end && *p != '>') {
        uint8_t info = lut[static_cast<uint8_t>(*p++)];
        if (info & 0x10) continue;
        if (!(info & 0x04)) { klen = 0; state = 0; continue; }
        uint8_t b = info & 0x03;
        uint32_t ns = ((state << 2) | b) & (N-1);
        if (klen >= K) T[state][ns] += 1.0;
        state = ns; ++klen;
    }
}

// Row-normalise a 64×64 matrix; rows with zero sum get a self-loop.
static void normalise_rows(double T[64][64]) {
    for (int i = 0; i < 64; ++i) {
        double s = 0; for (int j = 0; j < 64; ++j) s += T[i][j];
        if (s > 0) for (int j = 0; j < 64; ++j) T[i][j] /= s;
        else T[i][i] = 1.0;
    }
}

// Eigenvalue moduli (sorted descending) of a 64×64 row-stochastic matrix.
static std::vector<double> eigen_mods(const double T[64][64]) {
    Eigen::MatrixXd M(64, 64);
    for (int i = 0; i < 64; ++i) for (int j = 0; j < 64; ++j) M(i,j) = T[i][j];
    Eigen::EigenSolver<Eigen::MatrixXd> es(M, /*computeEigenvectors=*/false);
    auto ev = es.eigenvalues();
    std::vector<double> mods(64);
    for (int i = 0; i < 64; ++i) mods[i] = std::abs(ev[i]);
    std::sort(mods.rbegin(), mods.rend());
    return mods;
}

} // namespace

TransitionResult compute_transition_result(const std::string& fasta, int min_contig_len) {
    TransitionResult res;
    const auto& lut = char_lut();

    // Collect spans of long contigs (pointer pairs into fasta).
    struct Span { const char* begin; const char* end_ptr; uint64_t bp; };
    std::vector<Span> spans;
    {
        const char* p = fasta.data(), *end = p + fasta.size();
        while (p < end) {
            while (p < end && *p != '>') ++p; if (p >= end) break;
            ++p; while (p < end && *p != '\n') ++p; if (p < end) ++p;
            const char* cb = p; uint64_t bp = 0;
            while (p < end && *p != '>') { if (lut[static_cast<uint8_t>(*p++)] & 0x04) ++bp; }
            if (static_cast<int>(bp) >= min_contig_len) spans.push_back({cb, p, bp});
        }
    }
    if (spans.size() < 2) return res;

    // ── Full-genome transition matrix ─────────────────────────────────────────
    double T[64][64] = {};
    for (const auto& s : spans) accum_k3_transitions(s.begin, s.end_ptr, T, lut);
    normalise_rows(T);
    auto mods = eigen_mods(T);

    res.spectral_gap = static_cast<float>(1.0 - mods[1]);

    // ── Dyadic scale ladder → kink detection ──────────────────────────────────
    constexpr int scales[]  = {8192, 16384, 32768, 65536};
    constexpr int n_scales  = 4;
    float  gaps[n_scales]   = {};
    bool   valid[n_scales]  = {};

    for (int si = 0; si < n_scales; ++si) {
        const uint64_t W = static_cast<uint64_t>(scales[si]);
        double TW[64][64] = {};
        int n_win = 0;
        for (const auto& sp : spans) {
            if (sp.bp < W) continue;
            const char* wp = sp.begin;
            while (wp < sp.end_ptr) {
                const char* wend = wp; uint64_t bp = 0;
                while (wend < sp.end_ptr && bp < W) {
                    if (lut[static_cast<uint8_t>(*wend++)] & 0x04) ++bp;
                }
                if (bp >= W / 2) { accum_k3_transitions(wp, wend, TW, lut); ++n_win; }
                wp = wend;
            }
        }
        if (n_win == 0) continue;
        normalise_rows(TW);
        auto mw = eigen_mods(TW);
        gaps[si]  = static_cast<float>(1.0 - mw[1]);
        valid[si] = true;
    }

    float max_diff = -1.0f; int kink = -1;
    for (int i = 0; i + 1 < n_scales; ++i) {
        if (!valid[i] || !valid[i+1]) continue;
        float d = std::fabs(gaps[i+1] - gaps[i]);
        if (d > max_diff) { max_diff = d; kink = i; }
    }
    if (kink >= 0 && max_diff > 0.02f)
        res.scale_kink = static_cast<float>(std::log2(static_cast<double>(scales[kink])));

    return res;
}


MixtureResult compute_mixture_model(const std::string& fasta, int min_window, int min_windows) {
    MixtureResult res;
    const auto& lut       = char_lut();
    const auto& canon_idx = kmer4_canon_index();
    constexpr int D = 136; // canonical k=4 space; eliminates strand-bias noise

    // Thread-local storage: avoids per-genome heap allocation across 9.3M calls.
    struct Window { std::array<float, D> prof; uint64_t bp; double norm_sq; };
    thread_local std::vector<Window> windows;
    thread_local std::vector<int>    labels;
    thread_local std::vector<float>  bpi_U;   // n×3 column-major, spectral scratch
    thread_local std::vector<int>    sp_order; // window sort order for spectral cut
    windows.clear();

    struct Span { const char* begin; const char* end_ptr; };
    std::vector<Span> spans;
    {
        const char* p = fasta.data(), *end = p + fasta.size();
        while (p < end) {
            while (p < end && *p != '>') ++p; if (p >= end) break;
            ++p; while (p < end && *p != '\n') ++p; if (p < end) ++p;
            const char* cb = p;
            while (p < end && *p != '>') ++p;
            spans.push_back({cb, p});
        }
    }

    auto flush_window = [&](std::array<uint32_t, D>& buf, uint64_t bp) {
        float total = 0;
        for (int i = 0; i < D; ++i) total += buf[i];
        Window w; w.bp = bp; w.norm_sq = 0;
        if (total > 0)
            for (int i = 0; i < D; ++i) {
                w.prof[i] = buf[i] / total;
                w.norm_sq += static_cast<double>(w.prof[i]) * w.prof[i];
            }
        windows.push_back(w);
        buf = {};
    };

    for (const auto& sp : spans) {
        std::array<uint32_t, D> buf{};
        uint64_t bp = 0;
        uint8_t state = 0; int klen = 0; // state stays 8-bit (0-255); canon_idx maps to [0,135]
        for (const char* p = sp.begin; p < sp.end_ptr; ) {
            uint8_t info = lut[static_cast<uint8_t>(*p++)];
            if (info & 0x10) continue;
            if (!(info & 0x04)) { klen = 0; state = 0; continue; }
            state = static_cast<uint8_t>((state << 2) | (info & 0x03));
            if (++klen >= 4 && ++bp == static_cast<uint64_t>(min_window)) {
                buf[canon_idx[state]]++;
                flush_window(buf, bp); bp = 0; klen = 0;
            } else if (klen >= 4) {
                buf[canon_idx[state]]++;
            }
        }
        if (bp >= static_cast<uint64_t>(min_window) / 2 && bp > 0)
            flush_window(buf, bp);
    }

    const int n = static_cast<int>(windows.size());
    res.n_mix_windows = static_cast<uint16_t>(std::min(n, 65535));
    if (n < min_windows) { res.mix_no_data = true; return res; }

    // K=1: centroid + RSS via dot identity: ||x-c||² = ||x||² - 2·x·c + ||c||²
    std::array<double, D> c1{};
    for (const auto& w : windows)
        for (int d = 0; d < D; ++d) c1[d] += w.prof[d];
    for (int d = 0; d < D; ++d) c1[d] /= n;

    // ── Spectral split: Pearson-residual PCA + 1D GMM BIC ───────────────────────
    // Transform: X_r[i][d] = (prof[i][d] - mean[d]) / sqrt(mean[d] + ε)
    // This is the chi-squared / Correspondence Analysis metric — correct for compositional data.
    // Block power iteration finds top-3 left singular vectors without external libraries.

    // Pearson residual scale: pr_scale[d] = 1/sqrt(c1[d] + ε)
    std::array<float, D> pr_scale{};
    for (int d = 0; d < D; ++d)
        pr_scale[d] = 1.0f / std::sqrt(static_cast<float>(c1[d]) + 1e-6f);

    // V in R^{D×3}, U in R^{n×3} (column-major)
    constexpr int K3 = 3;
    std::array<float, D * K3> bpi_V{};
    bpi_U.resize(n * K3);

    // Deterministic seed: unit vectors at evenly spaced indices
    bpi_V[0]             = 1.0f;
    bpi_V[D + D/2]       = 1.0f;
    bpi_V[2*D + D - 1]   = 1.0f;

    auto xr_times_v = [&](const float* vk, float* uk) {
        for (int i = 0; i < n; ++i) {
            double s = 0;
            for (int d = 0; d < D; ++d)
                s += (windows[i].prof[d] - static_cast<float>(c1[d])) * pr_scale[d] * vk[d];
            uk[i] = static_cast<float>(s);
        }
    };
    auto xr_T_times_u = [&](const float* uk, float* vk) {
        std::fill(vk, vk + D, 0.0f);
        for (int i = 0; i < n; ++i) {
            const float ui = uk[i];
            for (int d = 0; d < D; ++d)
                vk[d] += (windows[i].prof[d] - static_cast<float>(c1[d])) * pr_scale[d] * ui;
        }
    };
    auto gs3 = [&]() { // Gram-Schmidt on 3 columns of bpi_V
        float* v0 = bpi_V.data(); float* v1 = v0+D; float* v2 = v1+D;
        auto nrm = [](float* v, int sz) {
            double s=0; for(int d=0;d<sz;++d) s+=v[d]*v[d];
            if(s>1e-30){float r=1.0f/std::sqrt(static_cast<float>(s)); for(int d=0;d<sz;++d) v[d]*=r;}
        };
        auto proj = [](const float* a, const float* b, int sz) {
            double s=0; for(int d=0;d<sz;++d) s+=a[d]*b[d]; return static_cast<float>(s);
        };
        nrm(v0, D);
        float p01 = proj(v1,v0,D); for(int d=0;d<D;++d) v1[d]-=p01*v0[d]; nrm(v1, D);
        float p02 = proj(v2,v0,D), p12 = proj(v2,v1,D);
        for(int d=0;d<D;++d) v2[d]-=p02*v0[d]+p12*v1[d]; nrm(v2, D);
    };

    for (int iter = 0; iter < 15; ++iter) {
        for (int k = 0; k < K3; ++k)
            xr_times_v(bpi_V.data()+k*D, bpi_U.data()+k*n);
        for (int k = 0; k < K3; ++k)
            xr_T_times_u(bpi_U.data()+k*n, bpi_V.data()+k*D);
        gs3();
    }
    // Final U pass to get PC1 scores
    for (int k = 0; k < K3; ++k)
        xr_times_v(bpi_V.data()+k*D, bpi_U.data()+k*n);

    // ── 1D Gaussian mixture BIC on PC1 projection ────────────────────────────
    const float* u0 = bpi_U.data();

    // K=1 log-likelihood
    double u_mean = 0;
    for (int i = 0; i < n; ++i) u_mean += u0[i];
    u_mean /= n;
    double u_var = 0;
    for (int i = 0; i < n; ++i) { double d = u0[i]-u_mean; u_var += d*d; }
    u_var /= n;
    const double ll1_1d = -0.5*n*std::log(2*M_PI*(u_var+1e-15)) - 0.5*n;
    const double bic1_1d = -2*ll1_1d + 2*std::log(static_cast<double>(n));

    // K=2 EM: initialise by median split of u0
    sp_order.resize(n);
    std::iota(sp_order.begin(), sp_order.end(), 0);
    std::sort(sp_order.begin(), sp_order.end(), [&](int a, int b){ return u0[a] < u0[b]; });
    const double u_med = u0[sp_order[n/2]];

    double mu1 = 0, mu2 = 0; int cnt1 = 0, cnt2 = 0;
    for (int i = 0; i < n; ++i) {
        if (u0[i] <= u_med) { mu1 += u0[i]; ++cnt1; }
        else                 { mu2 += u0[i]; ++cnt2; }
    }
    const double u_std = std::sqrt(u_var + 1e-15);
    mu1 = cnt1 > 0 ? mu1/cnt1 : u_mean - u_std;
    mu2 = cnt2 > 0 ? mu2/cnt2 : u_mean + u_std;
    double s1 = u_std, s2 = u_std, pi1 = static_cast<double>(cnt1)/n;

    std::vector<double> gmm_r1(n), gmm_r2(n);
    constexpr double kSqrt2Pi = 2.5066282746310002;
    auto gauss1d = [&](double x, double mu, double s) -> double {
        double z = (x-mu)/s; return std::exp(-0.5*z*z) / (s * kSqrt2Pi);
    };
    for (int iter = 0; iter < 15; ++iter) {
        for (int i = 0; i < n; ++i) {
            double p1 = pi1*(gauss1d)(u0[i],mu1,s1);
            double p2 = (1-pi1)*(gauss1d)(u0[i],mu2,s2);
            double tot = p1+p2+1e-300;
            gmm_r1[i] = p1/tot; gmm_r2[i] = p2/tot;
        }
        double n1=0,n2=0,m1=0,m2=0;
        for (int i=0;i<n;++i){ n1+=gmm_r1[i]; n2+=gmm_r2[i]; m1+=gmm_r1[i]*u0[i]; m2+=gmm_r2[i]*u0[i]; }
        n1=std::max(n1,1e-6); n2=std::max(n2,1e-6);
        mu1=m1/n1; mu2=m2/n2;
        double v1=0,v2=0;
        for (int i=0;i<n;++i){
            v1+=gmm_r1[i]*(u0[i]-mu1)*(u0[i]-mu1);
            v2+=gmm_r2[i]*(u0[i]-mu2)*(u0[i]-mu2);
        }
        s1=std::sqrt(v1/n1+1e-12); s2=std::sqrt(v2/n2+1e-12);
        pi1=n1/n;
    }
    double ll2_1d = 0;
    for (int i = 0; i < n; ++i)
        ll2_1d += std::log(pi1*(gauss1d)(u0[i],mu1,s1) + (1-pi1)*(gauss1d)(u0[i],mu2,s2) + 1e-300);
    const double bic2_1d = -2*ll2_1d + 5*std::log(static_cast<double>(n));

    // Pooled-σ separation: key continuous contamination signal
    // Cohen's d separation between the two fitted Gaussians.
    // Strand asymmetry (universal in bacteria) gives sep ≈ 0.5-1.5; cross-species contamination ≥ 2.
    // fiedler maps sep to [0,1): 0=unimodal, 0.49 at sep=2, 0.86 at sep=6.
    const double pooled_s = std::sqrt(0.5*(s1*s1+s2*s2) + 1e-12);
    const double sep      = std::fabs(mu2-mu1) / pooled_s;
    res.fiedler_value = static_cast<float>(1.0 - std::exp(-sep / 3.0));

    // bp-weighted cluster assignment from soft GMM responsibilities
    uint64_t len_a = 0, len_b = 0;
    for (int i = 0; i < n; ++i) {
        if (gmm_r1[i] >= 0.5) len_a += windows[i].bp;
        else                   len_b += windows[i].bp;
    }
    const uint64_t total = len_a + len_b;
    const double minor_frac = total > 0
        ? static_cast<double>(std::min(len_a,len_b)) / static_cast<double>(total) : 0.0;

    // Accept K=2 by pure BIC model selection — no secondary threshold guards
    if (bic1_1d > bic2_1d) {
        res.mixture_sources       = 2;
        res.contamination_mixture = total > 0 ? static_cast<float>(minor_frac) : NAN;
    } else {
        res.mixture_sources       = 1;
        res.contamination_mixture = 0.0f;
    }
    return res;
}

AllSignals compute_all_signals(const std::string& fasta,
                               int  min_contig_len,
                               int  mix_window,
                               int  min_mix_windows,
                               bool skip_mixture)
{
    AllSignals res;
    const auto& lut       = char_lut();
    const auto& canon_idx = kmer4_canon_index();

    constexpr int D4 = 256, D4c = 136, D3 = 64;

    // Collect contig spans (one scan of fasta bytes, no base decoding).
    struct Span { const char* begin; const char* end_ptr; };
    std::vector<Span> spans;
    {
        const char* p = fasta.data(), *end = p + fasta.size();
        while (p < end) {
            while (p < end && *p != '>') ++p; if (p >= end) break;
            ++p; while (p < end && *p != '\n') ++p; if (p < end) ++p;
            const char* cb = p;
            while (p < end && *p != '>') ++p;
            spans.push_back({cb, p});
        }
    }

    // Global accumulators (committed only from contigs >= min_contig_len).
    struct ContigProf { float p[D4c]; uint64_t bp; };
    thread_local std::vector<ContigProf> contig_profs;
    contig_profs.clear();
    uint64_t chargaff_counts[D4] = {};
    double   T_global[D3][D3]    = {};

    constexpr int n_sc = 4;
    constexpr int sc_sizes[n_sc] = {8192, 16384, 32768, 65536};
    double T_sc[n_sc][D3][D3] = {};

    // Mix-window accumulators (thread-local to avoid per-call alloc).
    // D4c=136 canonical space — matches compute_mixture_model, eliminates strand-bias noise.
    struct MixWin { std::array<float, D4c> prof; uint64_t bp; double norm_sq; };
    thread_local std::vector<MixWin> mix_wins;
    thread_local std::vector<int>    mix_labels;
    thread_local std::vector<float>  bpi_U;
    thread_local std::vector<int>    sp_order;
    mix_wins.clear();

    auto flush_mix_win = [&](uint32_t buf[D4c], uint64_t bp) {
        float total = 0;
        for (int i = 0; i < D4c; ++i) total += buf[i];
        MixWin w; w.bp = bp; w.norm_sq = 0;
        if (total > 0)
            for (int i = 0; i < D4c; ++i) {
                w.prof[i] = buf[i] / total;
                w.norm_sq += static_cast<double>(w.prof[i]) * w.prof[i];
            }
        mix_wins.push_back(w);
    };

    // Single walk per contig — accumulate into local buffers, commit at end.
    for (const auto& sp : spans) {
        uint32_t k4_raw[D4] = {};
        uint32_t k4_can[D4c] = {};
        double   T_loc[D3][D3] = {};
        uint32_t win_buf[D4c] = {}; // canonical 136-dim mix windows
        uint64_t contig_bp = 0, win_bp = 0;
        uint8_t  k4st = 0; int klen = 0;

        for (const char* p = sp.begin; p < sp.end_ptr; ) {
            uint8_t info = lut[static_cast<uint8_t>(*p++)];
            if (info & 0x10) continue;
            if (!(info & 0x04)) { k4st = 0; klen = 0; continue; }
            uint8_t b = info & 0x03;
            uint8_t old_k4 = k4st;
            k4st = static_cast<uint8_t>((k4st << 2) | b);

            if (klen >= 2) {
                // k3 transition: old k3 = old_k4 & 0x3F, new k3 = k4st & 0x3F
                T_loc[old_k4 & 0x3F][k4st & 0x3F] += 1.0;
            }
            if (klen >= 3) {
                ++k4_raw[k4st];
                ++k4_can[canon_idx[k4st]];
                ++win_buf[canon_idx[k4st]]; // canonical index for mix window
                ++contig_bp;
                if (++win_bp == static_cast<uint64_t>(mix_window)) {
                    flush_mix_win(win_buf, win_bp);
                    std::fill(win_buf, win_buf + D4c, 0u);
                    win_bp = 0;
                }
            }
            ++klen;
        }

        if (static_cast<int>(contig_bp) < min_contig_len) continue;

        // Partial tail mix window (>= half).
        if (win_bp >= static_cast<uint64_t>(mix_window) / 2 && win_bp > 0)
            flush_mix_win(win_buf, win_bp);

        // Commit k3 to global and per-scale T matrices.
        for (int i = 0; i < D3; ++i)
            for (int j = 0; j < D3; ++j) {
                T_global[i][j] += T_loc[i][j];
                for (int s = 0; s < n_sc; ++s)
                    if (contig_bp >= static_cast<uint64_t>(sc_sizes[s]) / 2)
                        T_sc[s][i][j] += T_loc[i][j];
            }

        // Commit k4 raw counts for chargaff.
        for (int i = 0; i < D4; ++i) chargaff_counts[i] += k4_raw[i];

        // Build L2-normalised canonical profile for self_coherence.
        uint64_t total_k = 0;
        for (int i = 0; i < D4c; ++i) total_k += k4_can[i];
        if (total_k > 0) {
            ContigProf cp; cp.bp = contig_bp;
            float inv = 1.0f / static_cast<float>(total_k), sumsq = 0;
            for (int i = 0; i < D4c; ++i) { cp.p[i] = k4_can[i] * inv; sumsq += cp.p[i]*cp.p[i]; }
            if (sumsq > 0) { float nrm = 1.0f / std::sqrt(sumsq); for (int i = 0; i < D4c; ++i) cp.p[i] *= nrm; }
            contig_profs.push_back(cp);
        }
    }

    // ── self_coherence ────────────────────────────────────────────────────────
    if (static_cast<int>(contig_profs.size()) >= 3) {
        float centroid[D4c] = {};
        for (const auto& cp : contig_profs)
            for (int i = 0; i < D4c; ++i) centroid[i] += cp.p[i];
        {
            float n = static_cast<float>(contig_profs.size()), sumsq = 0;
            for (int i = 0; i < D4c; ++i) { centroid[i] /= n; sumsq += centroid[i]*centroid[i]; }
            if (sumsq > 0) { float nrm = 1.0f/std::sqrt(sumsq); for (int i=0;i<D4c;++i) centroid[i]*=nrm; }
        }
        std::vector<float> dists(contig_profs.size());
        for (size_t i = 0; i < contig_profs.size(); ++i) {
            float d2 = 0;
            for (int j = 0; j < D4c; ++j) { float df = contig_profs[i].p[j]-centroid[j]; d2+=df*df; }
            dists[i] = std::sqrt(d2);
        }
        std::vector<float> sd = dists;
        std::sort(sd.begin(), sd.end());
        const size_t h = sd.size()/2;
        const float med = (sd.size()%2==0) ? 0.5f*(sd[h-1]+sd[h]) : sd[h];
        std::vector<float> ad(contig_profs.size());
        for (size_t i = 0; i < contig_profs.size(); ++i) ad[i] = std::fabs(dists[i]-med);
        std::sort(ad.begin(), ad.end());
        const float mad = (ad.size()%2==0) ? 0.5f*(ad[h-1]+ad[h]) : ad[h];
        const float thr = med + 3.0f * 1.4826f * mad;
        uint64_t clean_bp = 0, total_bp = 0;
        for (size_t i = 0; i < contig_profs.size(); ++i) {
            total_bp += contig_profs[i].bp;
            if (dists[i] <= thr) clean_bp += contig_profs[i].bp;
        }
        res.self_coherence = total_bp > 0 ? static_cast<float>(clean_bp)/static_cast<float>(total_bp) : NAN;
    }


    // ── chargaff_parity ───────────────────────────────────────────────────────
    {
        double sum = 0; int n = 0;
        for (int i = 0; i < D4; ++i) {
            uint8_t b0=i&3,b1=(i>>2)&3,b2=(i>>4)&3,b3=(i>>6)&3;
            int rc=((3-b0)<<6)|((3-b1)<<4)|((3-b2)<<2)|(3-b3);
            if (rc <= i) continue;
            uint64_t ci=chargaff_counts[i], crc=chargaff_counts[static_cast<uint8_t>(rc)];
            uint64_t tot=ci+crc; if (tot==0) continue;
            sum += 1.0 - static_cast<double>(ci>crc?ci-crc:crc-ci)/static_cast<double>(tot);
            ++n;
        }
        if (n > 0) res.chargaff_parity = static_cast<float>(sum/n);
    }

    // ── spectral_gap + scale_kink (shared eigendecomposition) ────────────────
    {
        normalise_rows(T_global);
        auto mods = eigen_mods(T_global);
        if (mods.size() >= 2) res.spectral_gap = static_cast<float>(1.0 - mods[1]);

        float gaps[n_sc] = {}; bool valid[n_sc] = {};
        for (int s = 0; s < n_sc; ++s) {
            normalise_rows(T_sc[s]);
            auto mw = eigen_mods(T_sc[s]);
            if (mw.size() >= 2) { gaps[s] = static_cast<float>(1.0-mw[1]); valid[s] = true; }
        }
        float max_diff = -1; int kink = -1;
        for (int i = 0; i+1 < n_sc; ++i) {
            if (!valid[i]||!valid[i+1]) continue;
            float d = std::fabs(gaps[i+1]-gaps[i]);
            if (d > max_diff) { max_diff = d; kink = i; }
        }
        if (kink >= 0 && max_diff > 0.02f)
            res.scale_kink = static_cast<float>(std::log2(static_cast<double>(sc_sizes[kink])));
    }

    // ── contamination_mixture (reuse optimised K-means from compute_mixture_model) ──
    // Skipped for tnf-only flagged genomes (skip_mixture=true) to save PCA+GMM cost.
    {
        const int n = static_cast<int>(mix_wins.size());
        res.n_mix_windows = static_cast<uint16_t>(std::min(n, 65535));
        if (!skip_mixture && n >= min_mix_windows) {
            constexpr int D = D4c; // 136-dim canonical, matches compute_mixture_model
            std::array<double, D> c1{};
            for (const auto& w : mix_wins)
                for (int d = 0; d < D; ++d) c1[d] += w.prof[d];
            for (int d = 0; d < D; ++d) c1[d] /= n;

            // ── Spectral split: Pearson-residual PCA + 1D GMM BIC ────────────────
            std::array<float, D> pr_scale{};
            for (int d = 0; d < D; ++d)
                pr_scale[d] = 1.0f / std::sqrt(static_cast<float>(c1[d]) + 1e-6f);

            constexpr int K3 = 3;
            std::array<float, D * K3> bpi_V{};
            bpi_U.resize(n * K3);

            bpi_V[0]           = 1.0f;
            bpi_V[D + D/2]     = 1.0f;
            bpi_V[2*D + D - 1] = 1.0f;

            auto xr_times_v_m = [&](const float* vk, float* uk) {
                for (int i = 0; i < n; ++i) {
                    double s = 0;
                    for (int d = 0; d < D; ++d)
                        s += (mix_wins[i].prof[d] - static_cast<float>(c1[d])) * pr_scale[d] * vk[d];
                    uk[i] = static_cast<float>(s);
                }
            };
            auto xr_T_times_u_m = [&](const float* uk, float* vk) {
                std::fill(vk, vk + D, 0.0f);
                for (int i = 0; i < n; ++i) {
                    const float ui = uk[i];
                    for (int d = 0; d < D; ++d)
                        vk[d] += (mix_wins[i].prof[d] - static_cast<float>(c1[d])) * pr_scale[d] * ui;
                }
            };
            auto gs3_m = [&]() {
                float* v0 = bpi_V.data(); float* v1 = v0+D; float* v2 = v1+D;
                auto nrm = [](float* v, int sz) {
                    double s=0; for(int d=0;d<sz;++d) s+=v[d]*v[d];
                    if(s>1e-30){float r=1.0f/std::sqrt(static_cast<float>(s)); for(int d=0;d<sz;++d) v[d]*=r;}
                };
                auto proj = [](const float* a, const float* b, int sz) {
                    double s=0; for(int d=0;d<sz;++d) s+=a[d]*b[d]; return static_cast<float>(s);
                };
                nrm(v0, D);
                float p01 = proj(v1,v0,D); for(int d=0;d<D;++d) v1[d]-=p01*v0[d]; nrm(v1, D);
                float p02 = proj(v2,v0,D), p12 = proj(v2,v1,D);
                for(int d=0;d<D;++d) v2[d]-=p02*v0[d]+p12*v1[d]; nrm(v2, D);
            };

            for (int iter = 0; iter < 15; ++iter) {
                for (int k = 0; k < K3; ++k)
                    xr_times_v_m(bpi_V.data()+k*D, bpi_U.data()+k*n);
                for (int k = 0; k < K3; ++k)
                    xr_T_times_u_m(bpi_U.data()+k*n, bpi_V.data()+k*D);
                gs3_m();
            }
            for (int k = 0; k < K3; ++k)
                xr_times_v_m(bpi_V.data()+k*D, bpi_U.data()+k*n);

            // ── 1D Gaussian mixture BIC on PC1 projection ────────────────────────
            const float* u0 = bpi_U.data();

            double u_mean = 0;
            for (int i = 0; i < n; ++i) u_mean += u0[i];
            u_mean /= n;
            double u_var = 0;
            for (int i = 0; i < n; ++i) { double d = u0[i]-u_mean; u_var += d*d; }
            u_var /= n;
            const double ll1_1d = -0.5*n*std::log(2*M_PI*(u_var+1e-15)) - 0.5*n;
            const double bic1_1d = -2*ll1_1d + 2*std::log(static_cast<double>(n));

            sp_order.resize(n);
            std::iota(sp_order.begin(), sp_order.end(), 0);
            std::sort(sp_order.begin(), sp_order.end(), [&](int a, int b){ return u0[a] < u0[b]; });
            const double u_med = u0[sp_order[n/2]];

            double mu1 = 0, mu2 = 0; int cnt1 = 0, cnt2 = 0;
            for (int i = 0; i < n; ++i) {
                if (u0[i] <= u_med) { mu1 += u0[i]; ++cnt1; }
                else                 { mu2 += u0[i]; ++cnt2; }
            }
            const double u_std = std::sqrt(u_var + 1e-15);
            mu1 = cnt1 > 0 ? mu1/cnt1 : u_mean - u_std;
            mu2 = cnt2 > 0 ? mu2/cnt2 : u_mean + u_std;
            double s1 = u_std, s2 = u_std, pi1 = static_cast<double>(cnt1)/n;

            constexpr double kSqrt2Pi_as = 2.5066282746310002;
            auto gauss1d_as = [&](double x, double mu, double s) -> double {
                double z = (x-mu)/s; return std::exp(-0.5*z*z) / (s * kSqrt2Pi_as);
            };
            std::vector<double> gmm_r1(n), gmm_r2(n);
            for (int iter = 0; iter < 15; ++iter) {
                for (int i = 0; i < n; ++i) {
                    double p1 = pi1*(gauss1d_as)(u0[i],mu1,s1);
                    double p2 = (1-pi1)*(gauss1d_as)(u0[i],mu2,s2);
                    double tot = p1+p2+1e-300;
                    gmm_r1[i] = p1/tot; gmm_r2[i] = p2/tot;
                }
                double n1=0,n2=0,m1=0,m2=0;
                for (int i=0;i<n;++i){ n1+=gmm_r1[i]; n2+=gmm_r2[i]; m1+=gmm_r1[i]*u0[i]; m2+=gmm_r2[i]*u0[i]; }
                n1=std::max(n1,1e-6); n2=std::max(n2,1e-6);
                mu1=m1/n1; mu2=m2/n2;
                double v1=0,v2=0;
                for (int i=0;i<n;++i){
                    v1+=gmm_r1[i]*(u0[i]-mu1)*(u0[i]-mu1);
                    v2+=gmm_r2[i]*(u0[i]-mu2)*(u0[i]-mu2);
                }
                s1=std::sqrt(v1/n1+1e-12); s2=std::sqrt(v2/n2+1e-12);
                pi1=n1/n;
            }
            double ll2_1d = 0;
            for (int i = 0; i < n; ++i)
                ll2_1d += std::log(pi1*(gauss1d_as)(u0[i],mu1,s1) + (1-pi1)*(gauss1d_as)(u0[i],mu2,s2) + 1e-300);
            const double bic2_1d = -2*ll2_1d + 5*std::log(static_cast<double>(n));

            const double pooled_s = std::sqrt(0.5*(s1*s1+s2*s2) + 1e-12);
            const double sep      = std::fabs(mu2-mu1) / pooled_s;
            res.fiedler_value = static_cast<float>(1.0 - std::exp(-sep / 3.0));

            uint64_t len_a = 0, len_b = 0;
            for (int i = 0; i < n; ++i) {
                if (gmm_r1[i] >= 0.5) len_a += mix_wins[i].bp;
                else                   len_b += mix_wins[i].bp;
            }
            const uint64_t total = len_a + len_b;
            const double minor_frac = total > 0
                ? static_cast<double>(std::min(len_a,len_b)) / static_cast<double>(total) : 0.0;

            if (bic1_1d > bic2_1d) {
                res.mixture_sources       = 2;
                res.contamination_mixture = total > 0 ? static_cast<float>(minor_frac) : NAN;
            } else {
                res.mixture_sources       = 1;
                res.contamination_mixture = 0.0f;
            }
        } else {
            res.mix_no_data = true;
        }
    }

    return res;
}

// Pre-parsed contig overload — eliminates the raw FASTA re-scan and duplicate TNF construction.
// Walks each contig's seq once for: k4_raw (chargaff), k3 transitions, mix windows.
// Uses caller-provided tnf136 (already computed by compute_tnf) for self_coherence contig_profs.
static std::atomic<int64_t> s_cas_scan{0}, s_cas_coh{0}, s_cas_chargaff{0}, s_cas_spectral{0}, s_cas_mix{0};
AllSignals compute_all_signals(std::span<const ContigAccum> contigs,
                               int  min_contig_len,
                               int  mix_window,
                               int  min_mix_windows,
                               bool skip_mixture)
{
    AllSignals res;
    using Clock = std::chrono::steady_clock;
    const auto t_cas0 = Clock::now();
    const auto& lut       = char_lut();
    const auto& canon_idx = kmer4_canon_index();

    constexpr int D4 = 256, D4c = 136, D3 = 64;

    struct ContigProf { float p[D4c]; uint64_t bp; };
    thread_local std::vector<ContigProf> contig_profs;
    contig_profs.clear();

    uint64_t chargaff_counts[D4] = {};
    double   T_global[D3][D3]    = {};

    constexpr int n_sc = 4;
    constexpr int sc_sizes[n_sc] = {8192, 16384, 32768, 65536};
    double T_sc[n_sc][D3][D3] = {};

    struct MixWin { std::array<float, D4c> prof; uint64_t bp; double norm_sq; };
    thread_local std::vector<MixWin> mix_wins;
    thread_local std::vector<int>    mix_labels;
    thread_local std::vector<float>  bpi_U;
    thread_local std::vector<int>    sp_order;
    mix_wins.clear();

    auto flush_mix_win = [&](uint32_t buf[D4c], uint64_t bp_w) {
        float total = 0;
        for (int i = 0; i < D4c; ++i) total += buf[i];
        MixWin w; w.bp = bp_w; w.norm_sq = 0;
        if (total > 0)
            for (int i = 0; i < D4c; ++i) {
                w.prof[i] = buf[i] / total;
                w.norm_sq += static_cast<double>(w.prof[i]) * w.prof[i];
            }
        mix_wins.push_back(w);
    };

    for (const auto& ca : contigs) {
        if (static_cast<int>(ca.len) < min_contig_len) continue;

        uint32_t k4_raw[D4]  = {};
        double   T_loc[D3][D3] = {};
        uint32_t win_buf[D4c] = {};
        uint64_t contig_bp = 0, win_bp = 0;
        uint8_t  k4st = 0; int klen = 0;

        for (const char* p = ca.seq.data(), *pe = p + ca.seq.size(); p < pe; ) {
            uint8_t info = lut[static_cast<uint8_t>(*p++)];
            if (info & 0x10) continue;
            if (!(info & 0x04)) { k4st = 0; klen = 0; continue; }
            uint8_t b     = info & 0x03;
            uint8_t old_k4 = k4st;
            k4st = static_cast<uint8_t>((k4st << 2) | b);
            if (klen >= 2) T_loc[old_k4 & 0x3F][k4st & 0x3F] += 1.0;
            if (klen >= 3) {
                ++k4_raw[k4st];
                ++win_buf[canon_idx[k4st]];
                ++contig_bp;
                if (++win_bp == static_cast<uint64_t>(mix_window)) {
                    flush_mix_win(win_buf, win_bp);
                    std::fill(win_buf, win_buf + D4c, 0u);
                    win_bp = 0;
                }
            }
            ++klen;
        }

        if (win_bp >= static_cast<uint64_t>(mix_window) / 2 && win_bp > 0)
            flush_mix_win(win_buf, win_bp);

        for (int i = 0; i < D3; ++i)
            for (int j = 0; j < D3; ++j) {
                T_global[i][j] += T_loc[i][j];
                for (int s = 0; s < n_sc; ++s)
                    if (contig_bp >= static_cast<uint64_t>(sc_sizes[s]) / 2)
                        T_sc[s][i][j] += T_loc[i][j];
            }
        for (int i = 0; i < D4; ++i) chargaff_counts[i] += k4_raw[i];

        // Use pre-computed TNF for contig_profs (self_coherence); skip rebuilding from k4_can.
        if (ca.tnf136) {
            ContigProf cp; cp.bp = ca.len;
            float sumsq = 0;
            for (int i = 0; i < D4c; ++i) { cp.p[i] = ca.tnf136[i]; sumsq += cp.p[i]*cp.p[i]; }
            if (sumsq > 0) {
                float nrm = 1.0f / std::sqrt(sumsq);
                for (int i = 0; i < D4c; ++i) cp.p[i] *= nrm;
            }
            contig_profs.push_back(cp);
        }
    }
    const auto t_cas1 = Clock::now();

    // ── self_coherence, chargaff, spectral, mixture ───────────────────────────
    // (identical post-processing as the original compute_all_signals)
    if (static_cast<int>(contig_profs.size()) >= 3) {
        float centroid[D4c] = {};
        for (const auto& cp : contig_profs)
            for (int i = 0; i < D4c; ++i) centroid[i] += cp.p[i];
        {
            float n = static_cast<float>(contig_profs.size()), sumsq = 0;
            for (int i = 0; i < D4c; ++i) { centroid[i] /= n; sumsq += centroid[i]*centroid[i]; }
            if (sumsq > 0) { float nrm = 1.0f/std::sqrt(sumsq); for (int i=0;i<D4c;++i) centroid[i]*=nrm; }
        }
        std::vector<float> dists(contig_profs.size());
        for (size_t i = 0; i < contig_profs.size(); ++i) {
            float d2 = 0;
            for (int j = 0; j < D4c; ++j) { float df = contig_profs[i].p[j]-centroid[j]; d2+=df*df; }
            dists[i] = std::sqrt(d2);
        }
        std::vector<float> sd = dists;
        std::sort(sd.begin(), sd.end());
        const size_t h = sd.size()/2;
        const float med = (sd.size()%2==0) ? 0.5f*(sd[h-1]+sd[h]) : sd[h];
        std::vector<float> ad(contig_profs.size());
        for (size_t i = 0; i < contig_profs.size(); ++i) ad[i] = std::fabs(dists[i]-med);
        std::sort(ad.begin(), ad.end());
        const float mad = (ad.size()%2==0) ? 0.5f*(ad[h-1]+ad[h]) : ad[h];
        const float thr = med + 3.0f * 1.4826f * mad;
        uint64_t clean_bp = 0, total_bp = 0;
        for (size_t i = 0; i < contig_profs.size(); ++i) {
            total_bp += contig_profs[i].bp;
            if (dists[i] <= thr) clean_bp += contig_profs[i].bp;
        }
        res.self_coherence = total_bp > 0 ? static_cast<float>(clean_bp)/static_cast<float>(total_bp) : NAN;
    }
    const auto t_cas2 = Clock::now();
    {
        double sum = 0; int n = 0;
        for (int i = 0; i < D4; ++i) {
            uint8_t b0=i&3,b1=(i>>2)&3,b2=(i>>4)&3,b3=(i>>6)&3;
            int rc=((3-b0)<<6)|((3-b1)<<4)|((3-b2)<<2)|(3-b3);
            if (rc <= i) continue;
            uint64_t ci=chargaff_counts[i], crc=chargaff_counts[static_cast<uint8_t>(rc)];
            uint64_t tot=ci+crc; if (tot==0) continue;
            sum += 1.0 - static_cast<double>(ci>crc?ci-crc:crc-ci)/static_cast<double>(tot);
            ++n;
        }
        if (n > 0) res.chargaff_parity = static_cast<float>(sum/n);
    }
    const auto t_cas3 = Clock::now();
    {
        normalise_rows(T_global);
        auto mods = eigen_mods(T_global);
        if (mods.size() >= 2) res.spectral_gap = static_cast<float>(1.0 - mods[1]);
        float gaps[n_sc] = {}; bool valid[n_sc] = {};
        for (int s = 0; s < n_sc; ++s) {
            normalise_rows(T_sc[s]);
            auto mw = eigen_mods(T_sc[s]);
            if (mw.size() >= 2) { gaps[s] = static_cast<float>(1.0-mw[1]); valid[s] = true; }
        }
        float max_diff = -1; int kink = -1;
        for (int i = 0; i+1 < n_sc; ++i) {
            if (!valid[i]||!valid[i+1]) continue;
            float d = std::fabs(gaps[i+1]-gaps[i]);
            if (d > max_diff) { max_diff = d; kink = i; }
        }
        if (kink >= 0 && max_diff > 0.02f)
            res.scale_kink = static_cast<float>(std::log2(static_cast<double>(sc_sizes[kink])));
    }
    const auto t_cas4 = Clock::now();
    {
        const int n = static_cast<int>(mix_wins.size());
        res.n_mix_windows = static_cast<uint16_t>(std::min(n, 65535));
        if (!skip_mixture && n >= min_mix_windows) {
            // Delegate to the original function's mixture logic by reconstructing a minimal
            // fasta and calling through — the mixture model is a small fraction of total time
            // and this avoids duplicating the 100-line GMM code.
            // Build a synthetic fasta string from contig sequences so we can call original.
            // Mix wins are already in mix_wins (thread_local, shared) — use them directly.

            // ── Spectral split: Pearson-residual PCA + 1D GMM BIC (duplicated from original) ──
            constexpr int D = D4c;
            std::array<double, D> c1{};
            for (const auto& w : mix_wins)
                for (int d = 0; d < D; ++d) c1[d] += w.prof[d];
            for (int d = 0; d < D; ++d) c1[d] /= n;

            std::array<float, D> pr_scale{};
            for (int d = 0; d < D; ++d)
                pr_scale[d] = 1.0f / std::sqrt(static_cast<float>(c1[d]) + 1e-6f);

            constexpr int K3 = 3;
            std::array<float, D * K3> bpi_V{};
            bpi_U.resize(n * K3);
            bpi_V[0]           = 1.0f;
            bpi_V[D + D/2]     = 1.0f;
            bpi_V[2*D + D - 1] = 1.0f;

            auto xr_times_v_m = [&](const float* vk, float* uk) {
                for (int i = 0; i < n; ++i) {
                    double s = 0;
                    for (int d = 0; d < D; ++d)
                        s += (mix_wins[i].prof[d] - static_cast<float>(c1[d])) * pr_scale[d] * vk[d];
                    uk[i] = static_cast<float>(s);
                }
            };
            auto xr_T_times_u_m = [&](const float* uk, float* vk) {
                std::fill(vk, vk + D, 0.0f);
                for (int i = 0; i < n; ++i) {
                    const float ui = uk[i];
                    for (int d = 0; d < D; ++d)
                        vk[d] += (mix_wins[i].prof[d] - static_cast<float>(c1[d])) * pr_scale[d] * ui;
                }
            };
            auto gs3_m = [&]() {
                float* v0 = bpi_V.data(); float* v1 = v0+D; float* v2 = v1+D;
                auto nrm = [](float* v, int sz) {
                    double s=0; for(int d=0;d<sz;++d) s+=v[d]*v[d];
                    if(s>1e-30){float r=1.0f/std::sqrt(static_cast<float>(s)); for(int d=0;d<sz;++d) v[d]*=r;}
                };
                auto proj = [](const float* a, const float* b, int sz) {
                    double s=0; for(int d=0;d<sz;++d) s+=a[d]*b[d]; return static_cast<float>(s);
                };
                nrm(v0, D);
                float p01 = proj(v1,v0,D); for(int d=0;d<D;++d) v1[d]-=p01*v0[d]; nrm(v1, D);
                float p02 = proj(v2,v0,D), p12 = proj(v2,v1,D);
                for(int d=0;d<D;++d) v2[d]-=p02*v0[d]+p12*v1[d]; nrm(v2, D);
            };
            for (int iter = 0; iter < 15; ++iter) {
                for (int k = 0; k < K3; ++k) xr_times_v_m(bpi_V.data()+k*D, bpi_U.data()+k*n);
                for (int k = 0; k < K3; ++k) xr_T_times_u_m(bpi_U.data()+k*n, bpi_V.data()+k*D);
                gs3_m();
            }
            for (int k = 0; k < K3; ++k) xr_times_v_m(bpi_V.data()+k*D, bpi_U.data()+k*n);

            const float* u0 = bpi_U.data();
            double u_mean = 0;
            for (int i = 0; i < n; ++i) u_mean += u0[i];
            u_mean /= n;
            double u_var = 0;
            for (int i = 0; i < n; ++i) { double d = u0[i]-u_mean; u_var += d*d; }
            u_var /= n;
            const double ll1_1d = -0.5*n*std::log(2*M_PI*(u_var+1e-15)) - 0.5*n;
            const double bic1_1d = -2*ll1_1d + 2*std::log(static_cast<double>(n));

            sp_order.resize(n);
            std::iota(sp_order.begin(), sp_order.end(), 0);
            std::sort(sp_order.begin(), sp_order.end(), [&](int a, int b){ return u0[a] < u0[b]; });
            const double u_med = u0[sp_order[n/2]];

            double mu1 = 0, mu2 = 0; int cnt1 = 0, cnt2 = 0;
            for (int i = 0; i < n; ++i) {
                if (u0[i] <= u_med) { mu1 += u0[i]; ++cnt1; }
                else                 { mu2 += u0[i]; ++cnt2; }
            }
            const double u_std = std::sqrt(u_var + 1e-15);
            mu1 = cnt1 > 0 ? mu1/cnt1 : u_mean - u_std;
            mu2 = cnt2 > 0 ? mu2/cnt2 : u_mean + u_std;
            double s1 = u_std, s2 = u_std, pi1 = static_cast<double>(cnt1)/n;

            constexpr double kSqrt2Pi_as = 2.5066282746310002;
            auto gauss1d_as = [&](double x, double mu, double s) -> double {
                double z = (x-mu)/s; return std::exp(-0.5*z*z) / (s * kSqrt2Pi_as);
            };
            std::vector<double> gmm_r1(n), gmm_r2(n);
            for (int iter = 0; iter < 15; ++iter) {
                for (int i = 0; i < n; ++i) {
                    double p1 = pi1*(gauss1d_as)(u0[i],mu1,s1);
                    double p2 = (1-pi1)*(gauss1d_as)(u0[i],mu2,s2);
                    double tot = p1+p2+1e-300; gmm_r1[i] = p1/tot; gmm_r2[i] = p2/tot;
                }
                double n1=0,n2=0,m1=0,m2=0;
                for (int i=0;i<n;++i){ n1+=gmm_r1[i]; n2+=gmm_r2[i]; m1+=gmm_r1[i]*u0[i]; m2+=gmm_r2[i]*u0[i]; }
                n1=std::max(n1,1e-6); n2=std::max(n2,1e-6);
                mu1=m1/n1; mu2=m2/n2;
                double v1=0,v2=0;
                for (int i=0;i<n;++i){
                    v1+=gmm_r1[i]*(u0[i]-mu1)*(u0[i]-mu1);
                    v2+=gmm_r2[i]*(u0[i]-mu2)*(u0[i]-mu2);
                }
                s1=std::sqrt(v1/n1+1e-12); s2=std::sqrt(v2/n2+1e-12); pi1=n1/n;
            }
            double ll2_1d = 0;
            for (int i = 0; i < n; ++i)
                ll2_1d += std::log(pi1*(gauss1d_as)(u0[i],mu1,s1) + (1-pi1)*(gauss1d_as)(u0[i],mu2,s2) + 1e-300);
            const double bic2_1d = -2*ll2_1d + 5*std::log(static_cast<double>(n));
            const double pooled_s = std::sqrt(0.5*(s1*s1+s2*s2) + 1e-12);
            const double sep      = std::fabs(mu2-mu1) / pooled_s;
            res.fiedler_value = static_cast<float>(1.0 - std::exp(-sep / 3.0));

            uint64_t len_a = 0, len_b = 0;
            for (int i = 0; i < n; ++i) {
                if (gmm_r1[i] >= 0.5) len_a += mix_wins[i].bp;
                else                   len_b += mix_wins[i].bp;
            }
            const uint64_t total = len_a + len_b;
            const double minor_frac = total > 0
                ? static_cast<double>(std::min(len_a,len_b)) / static_cast<double>(total) : 0.0;
            if (bic1_1d > bic2_1d) {
                res.mixture_sources       = 2;
                res.contamination_mixture = total > 0 ? static_cast<float>(minor_frac) : NAN;
            } else {
                res.mixture_sources       = 1;
                res.contamination_mixture = 0.0f;
            }
        } else {
            res.mix_no_data = true;
        }
    }
    const auto t_cas5 = Clock::now();
    s_cas_scan.fetch_add((t_cas1-t_cas0).count(), std::memory_order_relaxed);
    s_cas_coh.fetch_add((t_cas2-t_cas1).count(), std::memory_order_relaxed);
    s_cas_chargaff.fetch_add((t_cas3-t_cas2).count(), std::memory_order_relaxed);
    s_cas_spectral.fetch_add((t_cas4-t_cas3).count(), std::memory_order_relaxed);
    s_cas_mix.fetch_add((t_cas5-t_cas4).count(), std::memory_order_relaxed);
    return res;
}

CasTimings get_cas_timings() {
    return {s_cas_scan.load(), s_cas_coh.load(), s_cas_chargaff.load(),
            s_cas_spectral.load(), s_cas_mix.load()};
}

uint32_t days_since_epoch() {
    using namespace std::chrono;
    auto now = system_clock::now();
    auto tp  = floor<days>(now);
    auto ref = sys_days{year{2024}/January/1};
    return static_cast<uint32_t>((tp - ref).count());
}

std::vector<BuildRecord> parse_tsv_records(const std::filesystem::path& tsv_path) {
    // Read entire file in one syscall — avoids thousands of small NFS round-trips
    // that std::getline/ifstream (8 KB buffer) would generate on slow NFS mounts.
    std::string raw = read_file_bytes(tsv_path);
    const char* p   = raw.data();
    const char* end = raw.data() + raw.size();

    // Helper: consume one tab-delimited field from [p,end), advance p.
    auto next_field = [&](char delim) -> std::string_view {
        const char* start = p;
        while (p < end && *p != delim && *p != '\n' && *p != '\r') ++p;
        std::string_view sv(start, static_cast<size_t>(p - start));
        if (p < end && *p == delim) ++p;
        return sv;
    };
    auto skip_line = [&]() {
        while (p < end && *p != '\n') ++p;
        if (p < end) ++p;
    };

    // Parse header line into column names.
    std::vector<std::string> cols;
    {
        const char* line_end = p;
        while (line_end < end && *line_end != '\n') ++line_end;
        std::string header_line(p, static_cast<size_t>(line_end - p));
        if (!header_line.empty() && header_line.back() == '\r')
            header_line.pop_back();
        p = (line_end < end) ? line_end + 1 : end;
        std::istringstream ss(header_line);
        std::string tok;
        while (std::getline(ss, tok, '\t')) cols.push_back(tok);
    }
    auto find_col = [&](std::initializer_list<const char*> names) -> int {
        for (const char* name : names)
            for (int i = 0; i < (int)cols.size(); ++i)
                if (cols[i] == name) return i;
        return -1;
    };

    int idx_acc   = find_col({"accession", "acc"});
    int idx_path  = find_col({"file_path", "path", "fasta_path", "fasta", "file"});
    int idx_comp  = find_col({"completeness"});
    int idx_cont  = find_col({"contamination"});
    int idx_glen  = find_col({"genome_length"});
    int idx_ctg   = find_col({"n_contigs"});
    int idx_gtype = find_col({"genome_type"});

    if (idx_acc  < 0) throw std::runtime_error("TSV missing 'accession' column");
    if (idx_path < 0) throw std::runtime_error("TSV missing file path column"
        " (expected: file_path, path, fasta_path, fasta, or file)");

    std::set<int> known;
    if (idx_acc   >= 0) known.insert(idx_acc);
    if (idx_path  >= 0) known.insert(idx_path);
    if (idx_comp  >= 0) known.insert(idx_comp);
    if (idx_cont  >= 0) known.insert(idx_cont);
    if (idx_glen  >= 0) known.insert(idx_glen);
    if (idx_ctg   >= 0) known.insert(idx_ctg);
    if (idx_gtype >= 0) known.insert(idx_gtype);

    std::vector<int>         extra_indices;
    std::vector<std::string> extra_names;
    for (int i = 0; i < (int)cols.size(); ++i) {
        if (!known.count(i)) {
            extra_indices.push_back(i);
            extra_names.push_back(cols[i]);
        }
    }

    std::vector<BuildRecord> records;
    int max_col = std::max(idx_acc, idx_path);
    while (p < end) {
        if (*p == '#') { skip_line(); continue; }
        std::vector<std::string_view> fields;
        while (p < end && *p != '\n' && *p != '\r')
            fields.push_back(next_field('\t'));
        while (p < end && (*p == '\r' || *p == '\n')) ++p;
        if (fields.empty()) continue;
        if ((int)fields.size() <= max_col) continue;

        BuildRecord r;
        r.accession = std::string(fields[idx_acc]);
        r.file_path = std::string(fields[idx_path]);
        auto sv2f = [](std::string_view sv) { return std::stof(std::string(sv)); };
        auto sv2u = [](std::string_view sv) { return std::stoull(std::string(sv)); };
        if (idx_comp >= 0 && idx_comp < (int)fields.size())
            r.completeness  = sv2f(fields[idx_comp]);
        if (idx_cont >= 0 && idx_cont < (int)fields.size())
            r.contamination = sv2f(fields[idx_cont]);
        if (idx_glen >= 0 && idx_glen < (int)fields.size())
            r.genome_length  = sv2u(fields[idx_glen]);
        if (idx_ctg   >= 0 && idx_ctg  < (int)fields.size())
            r.n_contigs      = static_cast<uint32_t>(sv2u(fields[idx_ctg]));
        if (idx_gtype >= 0 && idx_gtype < (int)fields.size())
            r.genome_type    = parse_genome_type(fields[idx_gtype]);
        for (int j = 0; j < (int)extra_indices.size(); ++j) {
            int ci = extra_indices[j];
            r.extra_fields.emplace_back(extra_names[j],
                ci < (int)fields.size() ? std::string(fields[ci]) : "");
        }
        records.push_back(std::move(r));
    }
    return records;
}

std::vector<ContigProfile> compute_long_contig_profiles(const std::string& fasta,
                                                         uint32_t min_bp)
{
    const auto& lut  = char_lut();
    const auto& cidx = kmer4_canon_index();
    std::vector<ContigProfile> out;

    const char* ptr = fasta.data();
    const char* end = ptr + fasta.size();
    while (ptr < end) {
        while (ptr < end && *ptr != '>') ++ptr;
        if (ptr >= end) break;
        ++ptr;
        while (ptr < end && *ptr != '\n') ++ptr;
        if (ptr < end) ++ptr;

        std::array<uint32_t, 136> counts{};
        uint32_t mono[4] = {}, di[16] = {};
        uint64_t seq_len = 0;
        uint8_t enc = 0; int klen = 0; int prev_base = -1;
        while (ptr < end && *ptr != '>') {
            uint8_t info = lut[static_cast<uint8_t>(*ptr++)];
            if (info & 0x10) continue;
            if (!(info & 0x04)) { enc = 0; klen = 0; prev_base = -1; continue; }
            ++seq_len;
            int b = static_cast<int>(info & 0x03);
            mono[b]++;
            if (prev_base >= 0) di[prev_base*4 + b]++;
            prev_base = b;
            enc = static_cast<uint8_t>((enc << 2) | b);
            if (++klen >= 4) ++counts[cidx[enc]];
        }
        if (seq_len < min_bp) continue;

        uint64_t total = 0;
        for (uint32_t ct : counts) total += ct;
        if (total == 0) continue;

        ContigProfile cp;
        cp.bp = static_cast<uint32_t>(std::min<uint64_t>(seq_len, UINT32_MAX));
        float inv = 1.0f / static_cast<float>(total), sumsq = 0.0f;
        for (int i = 0; i < 136; ++i) { cp.p[i] = static_cast<float>(counts[i]) * inv; sumsq += cp.p[i] * cp.p[i]; }
        cp.freq_norm = std::sqrt(sumsq);
        if (sumsq > 0.0f) { float nrm = 1.0f / cp.freq_norm; for (int i = 0; i < 136; ++i) cp.p[i] *= nrm; }

        // Compute ρ*: dinucleotide relative abundance
        uint32_t n_mono = mono[0]+mono[1]+mono[2]+mono[3];
        uint32_t n_di = 0; for (int i=0;i<16;i++) n_di += di[i];
        if (n_mono > 0 && n_di > 0) {
            float fm[4], fd[16];
            for (int i=0;i<4;i++) fm[i] = static_cast<float>(mono[i]) / n_mono;
            for (int i=0;i<16;i++) fd[i] = static_cast<float>(di[i]) / n_di;
            for (int i=0;i<4;i++)
                for (int j=0;j<4;j++) {
                    float den = fm[i]*fm[j];
                    cp.rho[i*4+j] = (den > 1e-10f) ? fd[i*4+j] / den : 1.0f;
                }
        } else {
            for (int i=0;i<16;i++) cp.rho[i] = 1.0f;
        }
        out.push_back(cp);
    }
    return out;
}

} // namespace genopack
