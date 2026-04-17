#include <genopack/mem_delta.hpp>
#include <zstd.h>
#include <algorithm>
#include <cstring>
#include <stdexcept>
#ifdef __AVX2__
#  include <immintrin.h>
#endif

namespace genopack {

static constexpr int      K          = 31;
static constexpr uint32_t INDEX_STEP = 16;
static constexpr uint32_t PROBE_STEP = 15;
static constexpr int      BLOB_ZSTD_LEVEL = 3; // fast + effective for sorted integers
static constexpr uint32_t MEMD_MAGIC = 0x32444D47u; // "GMD2"

// ── Base encoding ─────────────────────────────────────────────────────────────

static inline uint8_t base2bit(char c) {
    static const uint8_t lut[256] = {
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
    return lut[static_cast<uint8_t>(c)];
}

static uint64_t kmer62(const char* s) {
    uint64_t h = 0;
    for (int i = 0; i < K; ++i) {
        uint8_t b = base2bit(s[i]);
        if (b > 3) return UINT64_MAX;
        h = (h << 2) | b;
    }
    return h;
}

// ── Flat open-addressing hash table ───────────────────────────────────────────

static size_t next_pow2(size_t n) {
    size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

AnchorIndex::AnchorIndex(size_t n_entries) {
    size_t cap = next_pow2(n_entries * 3 + 1);  // ~33% load factor
    keys.assign(cap, UINT64_MAX);
    vals.assign(cap, 0);
    mask = cap - 1;
}

void AnchorIndex::insert(uint64_t key, uint32_t val) {
    uint64_t slot = key & mask;
    while (keys[slot] != UINT64_MAX) {
        if (keys[slot] == key) return;  // duplicate: keep first position
        slot = (slot + 1) & mask;
    }
    keys[slot] = key;
    vals[slot] = val;
}

uint32_t AnchorIndex::find(uint64_t key) const {
    uint64_t slot = key & mask;
    while (keys[slot] != UINT64_MAX && keys[slot] != key)
        slot = (slot + 1) & mask;
    return keys[slot] == key ? vals[slot] : UINT32_MAX;
}

// ── FASTA parsing ─────────────────────────────────────────────────────────────

FastaComponents extract_fasta_components(const char* fasta, size_t len) {
    FastaComponents fc;
    const char* end = fasta + len;
    const char* p   = fasta;

    fc.seq.reserve(len);
    std::string cur_seq;
    cur_seq.reserve(1 << 20);

    while (p < end) {
        if (*p == '>') {
            if (!cur_seq.empty()) {
                fc.contig_ends.push_back(
                    static_cast<uint32_t>(fc.seq.size() + cur_seq.size()));
                fc.seq += cur_seq;
                cur_seq.clear();
            }
            ++p;
            const char* hstart = p;
            while (p < end && *p != '\n' && *p != '\r') ++p;
            fc.headers.emplace_back(hstart, p);
            while (p < end && (*p == '\n' || *p == '\r')) ++p;
        } else if (*p == '\n' || *p == '\r') {
            ++p;
        } else {
            const char* line_start = p;
            while (p < end && *p != '\n' && *p != '\r') ++p;
            size_t line_len = static_cast<size_t>(p - line_start);
            size_t old_size = cur_seq.size();
            cur_seq.resize(old_size + line_len);
            char* dst = &cur_seq[old_size];
            const char* src = line_start;
#ifdef __AVX2__
            {
                const __m256i vlo   = _mm256_set1_epi8('a' - 1);
                const __m256i vhi   = _mm256_set1_epi8('z');
                const __m256i vdiff = _mm256_set1_epi8(32);
                size_t i = 0;
                for (; i + 32 <= line_len; i += 32) {
                    __m256i v      = _mm256_loadu_si256((const __m256i*)(src + i));
                    __m256i gt_lo  = _mm256_cmpgt_epi8(v, vlo);
                    __m256i gt_hi  = _mm256_cmpgt_epi8(v, vhi);
                    __m256i is_lc  = _mm256_andnot_si256(gt_hi, gt_lo);
                    __m256i result = _mm256_blendv_epi8(v, _mm256_sub_epi8(v, vdiff), is_lc);
                    _mm256_storeu_si256((__m256i*)(dst + i), result);
                }
                for (; i < line_len; ++i) {
                    char c = src[i];
                    dst[i] = (c >= 'a' && c <= 'z') ? (char)(c - 32) : c;
                }
            }
#else
            for (size_t i = 0; i < line_len; ++i) {
                char c = src[i];
                dst[i] = (c >= 'a' && c <= 'z') ? (char)(c - 32) : c;
            }
#endif
        }
    }
    if (!cur_seq.empty()) {
        fc.contig_ends.push_back(
            static_cast<uint32_t>(fc.seq.size() + cur_seq.size()));
        fc.seq += cur_seq;
    }
    return fc;
}

// ── Index building ────────────────────────────────────────────────────────────

AnchorIndex build_anchor_index(const std::string& seq) {
    const size_t n = seq.size();
    if (n < static_cast<size_t>(K)) return AnchorIndex(0);

    size_t n_entries = (n / INDEX_STEP) + 1;
    AnchorIndex idx(n_entries);

    for (size_t i = 0; i + K <= n; i += INDEX_STEP) {
        uint64_t h = kmer62(seq.data() + i);
        if (h != UINT64_MAX)
            idx.insert(h, static_cast<uint32_t>(i));
    }
    return idx;
}

// ── MEM finding ───────────────────────────────────────────────────────────────

static std::vector<MemEntry> find_mems(const std::string& anchor,
                                        const std::string& query,
                                        const AnchorIndex& idx)
{
    const uint32_t alen = static_cast<uint32_t>(anchor.size());
    const uint32_t qlen = static_cast<uint32_t>(query.size());

    uint32_t covered_until = 0;
    std::vector<MemEntry> mems;

    static constexpr uint32_t PF_AHEAD = 12;

    for (uint32_t qi = 0; qi + static_cast<uint32_t>(K) <= qlen; qi += PROBE_STEP) {
        const uint32_t pf_qi = qi + PF_AHEAD * PROBE_STEP;
        if (pf_qi + static_cast<uint32_t>(K) <= qlen) {
            uint64_t ph = kmer62(query.data() + pf_qi);
            if (ph != UINT64_MAX) {
                uint64_t slot = ph & idx.mask;
                __builtin_prefetch(idx.keys.data() + slot, 0, 1);
                __builtin_prefetch(idx.vals.data() + slot, 0, 1);
            }
        }

        if (qi < covered_until) continue;

        uint64_t h = kmer62(query.data() + qi);
        if (h == UINT64_MAX) continue;

        uint32_t ri = idx.find(h);
        if (ri == UINT32_MAX) continue;

        // 62-bit hash → no false positives for ACGT k=31; memcmp unnecessary
        uint32_t qpos = qi, rpos = ri;
        while (qpos > covered_until && rpos > 0 &&
               query[qpos - 1] == anchor[rpos - 1]) {
            --qpos; --rpos;
        }

        uint32_t qend = qi + K, rend = ri + K;
        while (qend < qlen && rend < alen && query[qend] == anchor[rend]) {
            ++qend; ++rend;
        }

        uint32_t length = qend - qpos;
        covered_until = std::max(covered_until, qpos + length);
        mems.push_back({qpos, rpos, length});
    }

    // MEMs emitted in non-decreasing qpos order — no sort needed
    return mems;
}

// ── Serialisation helpers ─────────────────────────────────────────────────────

static void push_bytes(std::vector<uint8_t>& out, const void* d, size_t n) {
    const auto* p = static_cast<const uint8_t*>(d);
    out.insert(out.end(), p, p + n);
}

static void push_varint(std::vector<uint8_t>& out, uint64_t v) {
    while (v >= 0x80) {
        out.push_back(static_cast<uint8_t>(v) | 0x80u);
        v >>= 7;
    }
    out.push_back(static_cast<uint8_t>(v));
}

static uint32_t read32(const uint8_t*& p, const uint8_t* end) {
    if (p + 4 > end) throw std::runtime_error("decode_mem_delta: truncated");
    uint32_t v;
    std::memcpy(&v, p, 4);
    p += 4;
    return v;
}

static uint64_t read_varint(const uint8_t*& p, const uint8_t* end) {
    uint64_t value = 0;
    unsigned shift = 0;
    while (p < end && shift <= 63) {
        uint8_t byte = *p++;
        value |= static_cast<uint64_t>(byte & 0x7Fu) << shift;
        if ((byte & 0x80u) == 0)
            return value;
        shift += 7;
    }
    throw std::runtime_error("decode_mem_delta: malformed varint");
}

static double sample_anchor_hit_rate(const std::string& query,
                                     const AnchorIndex& idx,
                                     uint32_t start,
                                     uint32_t len) {
    if (query.size() < static_cast<size_t>(K) || len < static_cast<uint32_t>(K) || idx.keys.empty())
        return 0.0;
    uint32_t end = std::min<uint32_t>(static_cast<uint32_t>(query.size()), start + len);
    if (end <= start || end - start < static_cast<uint32_t>(K))
        return 0.0;

    uint32_t step = std::max<uint32_t>(256u, len / 8u);
    uint32_t hits = 0;
    uint32_t total = 0;
    for (uint32_t pos = start; pos + static_cast<uint32_t>(K) <= end; pos += step) {
        uint64_t h = kmer62(query.data() + pos);
        if (h != UINT64_MAX && idx.find(h) != UINT32_MAX)
            ++hits;
        ++total;
    }
    return total > 0 ? static_cast<double>(hits) / total : 0.0;
}

static std::vector<std::pair<uint32_t, uint32_t>>
build_query_chunks(const std::string& query,
                   const AnchorIndex& idx,
                   size_t chunk_bases) {
    std::vector<std::pair<uint32_t, uint32_t>> chunks;
    uint32_t query_len = static_cast<uint32_t>(query.size());
    if (query_len == 0) {
        chunks.emplace_back(0, 0);
        return chunks;
    }

    uint32_t step = static_cast<uint32_t>(std::max<size_t>(1, chunk_bases));
    uint32_t min_chunk = std::max<uint32_t>(8192u, step / 2u);
    uint32_t max_chunk = std::max<uint32_t>(step, step * 4u);
    for (uint32_t start = 0; start < query_len; ) {
        double hit_rate = sample_anchor_hit_rate(query, idx, start,
                                                 std::min<uint32_t>(step, query_len - start));
        uint32_t target = step;
        if (hit_rate >= 0.85) target = max_chunk;
        else if (hit_rate >= 0.65) target = std::min<uint32_t>(max_chunk, step * 2u);
        else if (hit_rate < 0.15) target = std::max<uint32_t>(8192u, step / 4u);
        else if (hit_rate < 0.35) target = min_chunk;

        uint32_t len = std::min<uint32_t>(target, query_len - start);
        chunks.emplace_back(start, len);
        start += len;
    }
    return chunks;
}

static std::vector<uint8_t> compress_payload(const std::vector<uint8_t>& payload) {
    size_t bound = ZSTD_compressBound(payload.size());
    std::vector<uint8_t> out(bound);
    size_t csize = ZSTD_compress(out.data(), bound,
                                 payload.data(), payload.size(),
                                 BLOB_ZSTD_LEVEL);
    if (ZSTD_isError(csize))
        throw std::runtime_error(std::string("encode_mem_delta: compress: ") +
                                 ZSTD_getErrorName(csize));
    out.resize(csize);
    return out;
}

static std::vector<uint8_t> decompress_payload(const uint8_t* blob,
                                               size_t blob_len,
                                               const char* err_prefix) {
    unsigned long long raw_size = ZSTD_getFrameContentSize(blob, blob_len);
    if (raw_size == ZSTD_CONTENTSIZE_ERROR || raw_size == ZSTD_CONTENTSIZE_UNKNOWN)
        throw std::runtime_error(std::string(err_prefix) + ": cannot determine decompressed size");

    std::vector<uint8_t> payload(static_cast<size_t>(raw_size));
    size_t dsize = ZSTD_decompress(payload.data(), payload.size(), blob, blob_len);
    if (ZSTD_isError(dsize))
        throw std::runtime_error(std::string(err_prefix) + ": " + ZSTD_getErrorName(dsize));
    payload.resize(dsize);
    return payload;
}

struct ParsedChunkPayload {
    const uint8_t* compressed_key = nullptr;
    size_t compressed_len = 0;
    std::vector<uint8_t> storage;
    std::vector<MemEntry> mems;
    size_t verbatim_offset = 0;
    size_t verbatim_len = 0;
};

static const ParsedChunkPayload& parse_chunk_payload(const uint8_t* blob,
                                                     size_t blob_len) {
    static thread_local ParsedChunkPayload cache;
    if (cache.compressed_key == blob && cache.compressed_len == blob_len)
        return cache;

    cache.compressed_key = blob;
    cache.compressed_len = blob_len;
    cache.storage = decompress_payload(blob, blob_len, "decode_mem_delta chunk decompress");
    cache.mems.clear();
    cache.verbatim_offset = 0;
    cache.verbatim_len = 0;

    const uint8_t* p = cache.storage.data();
    const uint8_t* end = cache.storage.data() + cache.storage.size();

    uint64_t n_mems = read_varint(p, end);
    cache.mems.reserve(static_cast<size_t>(n_mems));
    uint32_t prev_qpos = 0;
    for (uint64_t i = 0; i < n_mems; ++i) {
        MemEntry m;
        m.qpos = prev_qpos + static_cast<uint32_t>(read_varint(p, end));
        m.rpos = static_cast<uint32_t>(read_varint(p, end));
        m.length = static_cast<uint32_t>(read_varint(p, end));
        cache.mems.push_back(m);
        prev_qpos = m.qpos;
    }

    cache.verbatim_len = static_cast<size_t>(read_varint(p, end));
    if (p + cache.verbatim_len > end)
        throw std::runtime_error("decode_mem_delta: chunk verbatim overflow");
    cache.verbatim_offset = static_cast<size_t>(p - cache.storage.data());
    return cache;
}

static std::string render_chunk_slice(const std::string& anchor_seq,
                                      const ParsedChunkPayload& parsed,
                                      uint32_t query_len,
                                      uint32_t start,
                                      uint32_t length) {
    if (length == 0 || start >= query_len) return {};
    uint32_t end = std::min<uint32_t>(query_len, start + length);
    std::string out;
    out.reserve(end - start);

    const uint8_t* verbatim = parsed.storage.data() + parsed.verbatim_offset;
    size_t verbatim_pos = 0;
    uint32_t cursor = 0;

    for (const auto& m : parsed.mems) {
        if (m.rpos + m.length > anchor_seq.size() || m.qpos + m.length > query_len)
            throw std::runtime_error("decode_mem_delta: chunk MEM out of bounds");

        if (cursor < m.qpos) {
            uint32_t literal_start = cursor;
            uint32_t literal_end = m.qpos;
            if (end > literal_start && start < literal_end) {
                uint32_t copy_start = std::max(start, literal_start);
                uint32_t copy_end = std::min(end, literal_end);
                out.append(reinterpret_cast<const char*>(verbatim + verbatim_pos + (copy_start - literal_start)),
                           copy_end - copy_start);
            }
            verbatim_pos += literal_end - literal_start;
            cursor = literal_end;
        }

        uint32_t mem_end = m.qpos + m.length;
        if (end > m.qpos && start < mem_end) {
            uint32_t copy_start = std::max(start, m.qpos);
            uint32_t copy_end = std::min(end, mem_end);
            out.append(anchor_seq.data() + m.rpos + (copy_start - m.qpos),
                       copy_end - copy_start);
        }
        cursor = std::max(cursor, mem_end);
    }

    if (cursor < query_len) {
        if (end > cursor && start < query_len) {
            uint32_t copy_start = std::max(start, cursor);
            uint32_t copy_end = std::min(end, query_len);
            out.append(reinterpret_cast<const char*>(verbatim + verbatim_pos + (copy_start - cursor)),
                       copy_end - copy_start);
        }
        verbatim_pos += query_len - cursor;
    }

    if (verbatim_pos != parsed.verbatim_len)
        throw std::runtime_error("decode_mem_delta: chunk verbatim length mismatch");
    return out;
}

struct ParsedMemDeltaV2 {
    MemDeltaHeader hdr{};
    const uint32_t* contig_ends = nullptr;
    const char* headers = nullptr;
    const MemDeltaChunkDesc* chunks = nullptr;
    const uint8_t* base = nullptr;
    size_t size = 0;
};

static ParsedMemDeltaV2 parse_mem_delta_v2(const uint8_t* blob, size_t blob_len) {
    if (blob_len < sizeof(MemDeltaHeader))
        throw std::runtime_error("decode_mem_delta: blob too small");
    ParsedMemDeltaV2 parsed;
    parsed.base = blob;
    parsed.size = blob_len;
    std::memcpy(&parsed.hdr, blob, sizeof(MemDeltaHeader));
    if (parsed.hdr.magic != MEMD_MAGIC || parsed.hdr.version != 2)
        throw std::runtime_error("decode_mem_delta: bad v2 header");
    if (parsed.hdr.contig_ends_offset +
            static_cast<uint64_t>(parsed.hdr.n_contigs) * sizeof(uint32_t) > blob_len ||
        parsed.hdr.headers_offset + parsed.hdr.headers_len > blob_len ||
        parsed.hdr.chunks_offset +
            static_cast<uint64_t>(parsed.hdr.n_chunks) * sizeof(MemDeltaChunkDesc) > blob_len)
        throw std::runtime_error("decode_mem_delta: v2 header overflow");

    parsed.contig_ends = reinterpret_cast<const uint32_t*>(blob + parsed.hdr.contig_ends_offset);
    parsed.headers = reinterpret_cast<const char*>(blob + parsed.hdr.headers_offset);
    parsed.chunks = reinterpret_cast<const MemDeltaChunkDesc*>(blob + parsed.hdr.chunks_offset);
    for (uint32_t i = 0; i < parsed.hdr.n_chunks; ++i) {
        const auto& ch = parsed.chunks[i];
        if (static_cast<uint64_t>(ch.payload_offset) + ch.payload_size > blob_len)
            throw std::runtime_error("decode_mem_delta: chunk payload overflow");
    }
    return parsed;
}

static std::string reconstruct_fasta_from_seq(const std::string& seq,
                                              const uint32_t* contig_ends,
                                              uint32_t n_contigs,
                                              const char* headers,
                                              uint32_t headers_len) {
    std::vector<std::string> header_list;
    {
        size_t start = 0;
        for (size_t i = 0; i <= headers_len; ++i) {
            if (i == headers_len || headers[i] == '\0') {
                if (i > start) header_list.emplace_back(headers + start, i - start);
                start = i + 1;
            }
        }
    }

    std::string fasta;
    fasta.reserve(seq.size() + n_contigs * 64);
    uint32_t prev_end = 0;
    for (uint32_t ci = 0; ci < n_contigs; ++ci) {
        uint32_t this_end = contig_ends[ci];
        fasta += '>';
        if (ci < header_list.size()) fasta += header_list[ci];
        fasta += '\n';
        fasta.append(seq.data() + prev_end, this_end - prev_end);
        fasta += '\n';
        prev_end = this_end;
    }
    return fasta;
}

static std::string decode_mem_delta_legacy_full(const std::string& anchor_seq,
                                                const uint8_t* blob,
                                                size_t blob_len) {
    auto payload = decompress_payload(blob, blob_len, "decode_mem_delta decompress");
    const uint8_t* p   = payload.data();
    const uint8_t* end = payload.data() + payload.size();

    uint32_t anchor_seq_len = read32(p, end);
    uint32_t query_seq_len  = read32(p, end);
    if (anchor_seq_len > 0 && anchor_seq.size() != anchor_seq_len)
        throw std::runtime_error("decode_mem_delta: anchor sequence length mismatch");

    uint32_t n_contigs = read32(p, end);
    std::vector<uint32_t> contig_ends(n_contigs);
    for (auto& ce : contig_ends) ce = read32(p, end);

    uint32_t headers_len = read32(p, end);
    if (p + headers_len > end) throw std::runtime_error("decode_mem_delta: headers overflow");
    const char* headers = reinterpret_cast<const char*>(p);
    p += headers_len;

    uint32_t n_mems = read32(p, end);
    if (p + static_cast<size_t>(n_mems) * 12 > end)
        throw std::runtime_error("decode_mem_delta: mems overflow");
    std::vector<MemEntry> mems(n_mems);
    for (auto& m : mems) {
        m.qpos   = read32(p, end);
        m.rpos   = read32(p, end);
        m.length = read32(p, end);
    }

    uint32_t verbatim_len = read32(p, end);
    if (p + verbatim_len > end)
        throw std::runtime_error("decode_mem_delta: verbatim overflow");
    const uint8_t* verbatim = p;

    std::string seq(query_seq_len, 'N');
    std::vector<uint8_t> cov(query_seq_len, 0);
    for (const auto& m : mems) {
        if (m.rpos + m.length > anchor_seq.size() || m.qpos + m.length > query_seq_len)
            throw std::runtime_error("decode_mem_delta: MEM out of bounds");
        std::memcpy(&seq[m.qpos], anchor_seq.data() + m.rpos, m.length);
        std::fill(cov.begin() + m.qpos, cov.begin() + m.qpos + m.length, 1);
    }

    size_t vi = 0;
    for (uint32_t i = 0; i < query_seq_len; ++i) {
        if (!cov[i]) {
            if (vi >= verbatim_len)
                throw std::runtime_error("decode_mem_delta: verbatim underflow");
            seq[i] = static_cast<char>(verbatim[vi++]);
        }
    }

    return reconstruct_fasta_from_seq(seq, contig_ends.data(), n_contigs, headers, headers_len);
}

// ── Encode ────────────────────────────────────────────────────────────────────

std::vector<uint8_t> encode_mem_delta(const FastaComponents& anchor,
                                       const AnchorIndex&     anchor_idx,
                                       const FastaComponents& query,
                                       size_t                 chunk_bases)
{
    const uint32_t qlen = static_cast<uint32_t>(query.seq.size());

    // Flatten headers
    std::string headers_flat;
    for (const auto& h : query.headers) {
        headers_flat += h;
        headers_flat += '\0';
    }

    auto chunks = build_query_chunks(query.seq, anchor_idx, chunk_bases);
    std::vector<MemDeltaChunkDesc> descs;
    descs.reserve(chunks.size());
    std::vector<std::vector<uint8_t>> payloads;
    payloads.reserve(chunks.size());

    uint32_t payload_offset = static_cast<uint32_t>(
        sizeof(MemDeltaHeader) +
        query.contig_ends.size() * sizeof(uint32_t) +
        headers_flat.size() +
        chunks.size() * sizeof(MemDeltaChunkDesc));

    for (const auto& [chunk_start, chunk_len] : chunks) {
        std::string chunk_query = query.seq.substr(chunk_start, chunk_len);
        auto mems = find_mems(anchor.seq, chunk_query, anchor_idx);

        std::vector<uint8_t> verbatim;
        verbatim.reserve(chunk_len / 20 + 16);
        uint32_t next_verb = 0;
        for (const auto& m : mems) {
            for (uint32_t i = next_verb; i < m.qpos; ++i)
                verbatim.push_back(static_cast<uint8_t>(chunk_query[i]));
            if (m.qpos + m.length > next_verb)
                next_verb = m.qpos + m.length;
        }
        for (uint32_t i = next_verb; i < chunk_len; ++i)
            verbatim.push_back(static_cast<uint8_t>(chunk_query[i]));

        std::vector<uint8_t> payload;
        payload.reserve(8 + mems.size() * 12 + verbatim.size());
        push_varint(payload, static_cast<uint32_t>(mems.size()));
        uint32_t prev_qpos = 0;
        for (const auto& m : mems) {
            push_varint(payload, m.qpos - prev_qpos);
            push_varint(payload, m.rpos);
            push_varint(payload, m.length);
            prev_qpos = m.qpos;
        }
        push_varint(payload, static_cast<uint32_t>(verbatim.size()));
        push_bytes(payload, verbatim.data(), verbatim.size());

        auto compressed = compress_payload(payload);
        descs.push_back(MemDeltaChunkDesc{
            chunk_start,
            chunk_len,
            payload_offset,
            static_cast<uint32_t>(compressed.size()),
        });
        payload_offset += static_cast<uint32_t>(compressed.size());
        payloads.push_back(std::move(compressed));
    }

    MemDeltaHeader hdr{};
    hdr.magic = MEMD_MAGIC;
    hdr.version = 2;
    hdr.flags = 0;
    hdr.anchor_seq_len = static_cast<uint32_t>(anchor.seq.size());
    hdr.query_seq_len = qlen;
    hdr.n_contigs = static_cast<uint32_t>(query.contig_ends.size());
    hdr.n_chunks = static_cast<uint32_t>(descs.size());
    hdr.headers_len = static_cast<uint32_t>(headers_flat.size());
    hdr.reserved = 0;
    hdr.contig_ends_offset = sizeof(MemDeltaHeader);
    hdr.headers_offset = hdr.contig_ends_offset + query.contig_ends.size() * sizeof(uint32_t);
    hdr.chunks_offset = hdr.headers_offset + headers_flat.size();

    std::vector<uint8_t> out;
    out.reserve(payload_offset);
    push_bytes(out, &hdr, sizeof(hdr));
    if (!query.contig_ends.empty())
        push_bytes(out, query.contig_ends.data(), query.contig_ends.size() * sizeof(uint32_t));
    if (!headers_flat.empty())
        push_bytes(out, headers_flat.data(), headers_flat.size());
    if (!descs.empty())
        push_bytes(out, descs.data(), descs.size() * sizeof(MemDeltaChunkDesc));
    for (const auto& payload : payloads)
        push_bytes(out, payload.data(), payload.size());
    return out;
}

// ── Decode ────────────────────────────────────────────────────────────────────

std::string decode_mem_delta(const std::string& anchor_seq,
                              const uint8_t*     blob,
                              size_t             blob_len)
{
    if (blob_len >= sizeof(uint32_t) &&
        *reinterpret_cast<const uint32_t*>(blob) == MEMD_MAGIC) {
        auto parsed = parse_mem_delta_v2(blob, blob_len);
        if (parsed.hdr.anchor_seq_len > 0 &&
            parsed.hdr.anchor_seq_len != anchor_seq.size())
            throw std::runtime_error("decode_mem_delta: anchor sequence length mismatch");

        std::string seq;
        seq.reserve(parsed.hdr.query_seq_len);
        for (uint32_t i = 0; i < parsed.hdr.n_chunks; ++i) {
            const auto& ch = parsed.chunks[i];
            const auto& payload = parse_chunk_payload(parsed.base + ch.payload_offset,
                                                      ch.payload_size);
            std::string chunk = render_chunk_slice(anchor_seq, payload, ch.query_len, 0, ch.query_len);
            seq += chunk;
        }
        return reconstruct_fasta_from_seq(seq,
                                          parsed.contig_ends,
                                          parsed.hdr.n_contigs,
                                          parsed.headers,
                                          parsed.hdr.headers_len);
    }
    return decode_mem_delta_legacy_full(anchor_seq, blob, blob_len);
}

std::string decode_mem_delta_slice(const std::string& anchor_seq,
                                   const uint8_t*     blob,
                                   size_t             blob_len,
                                   uint64_t           start,
                                   uint64_t           length)
{
    if (length == 0) return {};
    if (blob_len >= sizeof(uint32_t) &&
        *reinterpret_cast<const uint32_t*>(blob) == MEMD_MAGIC) {
        auto parsed = parse_mem_delta_v2(blob, blob_len);
        if (start >= parsed.hdr.query_seq_len) return {};
        uint64_t end = std::min<uint64_t>(parsed.hdr.query_seq_len, start + length);

        std::string out;
        out.reserve(end - start);
        for (uint32_t i = 0; i < parsed.hdr.n_chunks && out.size() < end - start; ++i) {
            const auto& ch = parsed.chunks[i];
            uint64_t chunk_end = static_cast<uint64_t>(ch.query_start) + ch.query_len;
            if (chunk_end <= start || ch.query_start >= end) continue;
            uint64_t local_start = (start > ch.query_start) ? (start - ch.query_start) : 0;
            uint64_t local_end = std::min<uint64_t>(ch.query_len, end - ch.query_start);
            const auto& payload = parse_chunk_payload(parsed.base + ch.payload_offset,
                                                      ch.payload_size);
            out += render_chunk_slice(anchor_seq,
                                      payload,
                                      ch.query_len,
                                      static_cast<uint32_t>(local_start),
                                      static_cast<uint32_t>(local_end - local_start));
        }
        return out;
    }

    std::string fasta = decode_mem_delta_legacy_full(anchor_seq, blob, blob_len);
    auto fc = extract_fasta_components(fasta.data(), fasta.size());
    if (start >= fc.seq.size()) return {};
    uint64_t end = std::min<uint64_t>(fc.seq.size(), start + length);
    return fc.seq.substr(start, end - start);
}

} // namespace genopack
