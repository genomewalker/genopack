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
    for (int attempt = 0; attempt < 5; ++attempt) {
        std::ifstream f(path, std::ios::binary | std::ios::ate);
        if (f) {
            size_t sz = static_cast<size_t>(f.tellg());
            f.seekg(0);
            std::string buf(sz, '\0');
            f.read(buf.data(), static_cast<std::streamsize>(sz));
            if (!f) throw std::runtime_error("Read error: " + path.string());
            return buf;
        }
        if (attempt < 4)
            std::this_thread::sleep_for(std::chrono::seconds(1 << attempt)); // 1, 2, 4, 8 s
    }
    throw std::runtime_error("Cannot open: " + path.string());
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



namespace {

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

// Pre-parsed contig overload — eliminates the raw FASTA re-scan and duplicate TNF construction.
// Walks each contig's seq once for: k4_raw (chargaff), k3 transitions, mix windows.
// Walks each contig's sequence once to accumulate canonical-k4 mix windows.
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

    constexpr int D4c = 136;

    struct MixWin { std::array<float, D4c> prof; uint64_t bp; double norm_sq; };
    thread_local std::vector<MixWin> mix_wins;
    thread_local std::vector<int>    mix_win_contig;   // originating contig index per window (near-clade guard)
    thread_local std::vector<float>  bpi_U;
    thread_local std::vector<int>    sp_order;
    mix_wins.clear();
    mix_win_contig.clear();
    int mix_cur_contig = -1;

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
        mix_win_contig.push_back(mix_cur_contig);
    };

    // Single walk per contig: accumulate canonical-k4 windows for the near-clade TNF-GMM.
    for (const auto& ca : contigs) {
        if (static_cast<int>(ca.len) < min_contig_len) continue;
        ++mix_cur_contig;

        uint32_t win_buf[D4c] = {};
        uint64_t win_bp = 0;
        uint8_t  k4st = 0; int klen = 0;

        for (const char* p = ca.seq.data(), *pe = p + ca.seq.size(); p < pe; ) {
            uint8_t info = lut[static_cast<uint8_t>(*p++)];
            if (info & 0x10) continue;
            if (!(info & 0x04)) { k4st = 0; klen = 0; continue; }
            uint8_t b     = info & 0x03;
            k4st = static_cast<uint8_t>((k4st << 2) | b);
            if (klen >= 3) {
                ++win_buf[canon_idx[k4st]];
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
    }
    const auto t_cas1 = Clock::now();

    // ── near-clade TNF-GMM minority mass (contamination_tnf_minor) ─────────────
    {
        const int n = static_cast<int>(mix_wins.size());
        if (!skip_mixture && n >= min_mix_windows) {
            // Pearson-residual PCA + 1D GMM BIC over the k4 windows.
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

            uint64_t len_a = 0, len_b = 0;
            for (int i = 0; i < n; ++i) {
                if (gmm_r1[i] >= 0.5) len_a += mix_wins[i].bp;
                else                   len_b += mix_wins[i].bp;
            }
            const uint64_t total = len_a + len_b;
            const double minor_frac = total > 0
                ? static_cast<double>(std::min(len_a,len_b)) / static_cast<double>(total) : 0.0;

            // Near-clade contamination: surface the TNF-GMM minority mass even when BIC
            // prefers K=1. Near-clade (same-species/genus) sub-populations are composition-
            // ally close, so the BIC penalty (∝ log n, harsh on few-window genomes) rejects
            // a coherent minority that a duplicated-marker method (CheckM2) still detects.
            // Guarded against legitimate strain accessory genome: (a) Cohen's-d separation
            // above the strand-asymmetry band, and (b) the minority must span >=2 distinct
            // contigs — a plasmid/accessory island lives on one contig.
            // ceiling: sep floor 1.6 sits above normal strand asymmetry (0.5-1.5) and below
            //   cross-species (>=2); an identical-strain contaminant (sep~0) stays invisible
            //   to composition alone and needs read coverage. Thresholds are calibrated by
            //   check/validate/chimera_panel.sh.
            float tnf_minor = 0.0f;
            if (total > 0 && sep >= 1.6 && minor_frac >= 0.03) {
                const bool minority_in_a = (len_a <= len_b);
                int prev_contig = -1, n_minor_contigs = 0;   // recorded contig ids are always >=0
                for (int i = 0; i < n; ++i) {
                    const bool win_in_a = (gmm_r1[i] >= 0.5);
                    if (win_in_a == minority_in_a) {            // window belongs to minority component
                        if (mix_win_contig[i] != prev_contig) { // windows are emitted in contig order
                            ++n_minor_contigs;
                            prev_contig = mix_win_contig[i];
                        }
                    }
                }
                if (n_minor_contigs >= 2) tnf_minor = static_cast<float>(minor_frac);
            }
            res.contamination_tnf_minor = tnf_minor;
        }
    }
    const auto t_cas5 = Clock::now();
    s_cas_scan.fetch_add((t_cas1-t_cas0).count(), std::memory_order_relaxed);
    s_cas_mix.fetch_add((t_cas5-t_cas1).count(), std::memory_order_relaxed);
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
