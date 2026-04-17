#include <genopack/util.hpp>
#include <zlib.h>
#ifdef GENOPACK_HAVE_LIBDEFLATE
#  include <libdeflate.h>
#endif
#include <thread>
#include <algorithm>
#include <array>
#include <chrono>
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
        if (n <= 0) { ::close(fd); throw std::runtime_error("Read error: " + path.string()); }
        total += static_cast<size_t>(n);
    }
    ::close(fd);
    return buf;
}

static bool is_gzip(const std::string& buf) {
    return buf.size() >= 2 &&
           static_cast<uint8_t>(buf[0]) == 0x1f &&
           static_cast<uint8_t>(buf[1]) == 0x8b;
}

// Shared decompression core — operates on already-loaded compressed bytes.
static std::string decompress_gz_buf(std::string compressed,
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

} // namespace genopack
