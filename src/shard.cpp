#include <genopack/shard.hpp>
#include <genopack/mem_delta.hpp>
#include <zstd.h>
#include <zdict.h>
#include <spdlog/spdlog.h>
#include <sys/mman.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <fstream>
#include <future>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string_view>
#include <thread>
#include <vector>

namespace genopack {

namespace {

struct ChunkBoundary {
    uint32_t byte_start = 0;
    uint32_t byte_len = 0;
    uint64_t symbol_offset = 0;
};

struct PackedCheckpointLayout {
    uint16_t seq_byte_offset = 0;
    uint16_t line_bases = 0;
};

static constexpr uint32_t SEQB_MAGIC             = 0x42514553u; // "SEQB"
static constexpr uint16_t SEQB_FLAG_PACKED_2BIT  = 0x0001;     // nucleotides packed 2 bits/base
static constexpr float    TWOBIT_AMBIG_THRESHOLD = 0.05f;      // skip 2-bit if >5% ambiguous

// Ambiguity sidecar entry stored in blob prefix when SEQB_FLAG_PACKED_2BIT is set.
#pragma pack(push, 1)
struct AmbigEntry {
    uint32_t pos;       // base position (0-based)
    uint8_t  base;      // original ASCII character
    uint8_t  _pad[3];
};
#pragma pack(pop)
static_assert(sizeof(AmbigEntry) == 8, "AmbigEntry layout changed");

static inline uint8_t base_to_2bit(char c) noexcept {
    switch (c | 0x20) {
        case 'a': return 0;
        case 'c': return 1;
        case 'g': return 2;
        case 't': case 'u': return 3;
        default:  return 0xFF; // ambiguous
    }
}

// Pack ASCII nucleotide sequence to 2 bits/base.
// Returns false if ambiguity fraction exceeds threshold (caller should skip packing).
// hdr.reserved is set to n_ambiguous; ambig_out holds sidecar entries for the prefix.
static bool pack_sequence_2bit(const char* seq, size_t seq_len,
                               std::vector<uint8_t>& packed_out,
                               std::vector<AmbigEntry>& ambig_out) noexcept
{
    if (seq_len == 0) { packed_out.clear(); ambig_out.clear(); return true; }

    // Count ambiguous bases first
    size_t n_ambig = 0;
    for (size_t i = 0; i < seq_len; ++i)
        if (base_to_2bit(seq[i]) == 0xFF) ++n_ambig;

    if (n_ambig > static_cast<size_t>(seq_len * TWOBIT_AMBIG_THRESHOLD))
        return false;

    ambig_out.clear();
    ambig_out.reserve(n_ambig);
    size_t packed_len = (seq_len + 3) / 4;
    packed_out.assign(packed_len, 0);

    for (size_t i = 0; i < seq_len; ++i) {
        uint8_t bits = base_to_2bit(seq[i]);
        if (bits == 0xFF) {
            ambig_out.push_back(AmbigEntry{static_cast<uint32_t>(i),
                                           static_cast<uint8_t>(seq[i]), {0,0,0}});
            bits = 0; // encode as A in packed stream
        }
        // MSB-first: position 0 → bits 7:6, position 1 → bits 5:4, etc.
        packed_out[i / 4] |= static_cast<uint8_t>(bits << (6 - (i % 4) * 2));
    }
    return true;
}

// Unpack 2-bit packed sequence back to ASCII, applying ambiguity sidecar patches.
static std::string unpack_sequence_2bit(const uint8_t* packed, size_t /*packed_len*/,
                                        const AmbigEntry* ambig, uint32_t n_ambig,
                                        uint32_t seq_len)
{
    static constexpr char BITS_TO_BASE[4] = {'A', 'C', 'G', 'T'};
    std::string seq(seq_len, 'A');
    for (uint32_t i = 0; i < seq_len; ++i) {
        uint8_t bits = (packed[i / 4] >> (6 - (i % 4) * 2)) & 0x03;
        seq[i] = BITS_TO_BASE[bits];
    }
    for (uint32_t j = 0; j < n_ambig; ++j)
        if (ambig[j].pos < seq_len)
            seq[ambig[j].pos] = static_cast<char>(ambig[j].base);
    return seq;
}

#pragma pack(push, 1)
struct SequenceBlobHeader {
    uint32_t magic;
    uint16_t version;
    uint16_t flags;
    uint32_t seq_len;
    uint32_t n_contigs;
    uint32_t headers_len;
    uint32_t reserved;
    uint64_t contig_ends_offset;
    uint64_t headers_offset;
    uint64_t blocks_offset;
};
#pragma pack(pop)
static_assert(sizeof(SequenceBlobHeader) == 48, "SequenceBlobHeader layout changed");

struct ParsedSequenceBlob {
    SequenceBlobHeader hdr{};
    const uint32_t* contig_ends = nullptr;
    const char* headers = nullptr;
    const uint8_t* blocks = nullptr;
};

struct SequenceBlobLayout {
    std::string seq;
    std::vector<uint32_t> contig_ends;
    std::string headers_flat;
    std::vector<uint8_t> prefix;
    uint32_t blocks_offset = 0;
};

static std::vector<ChunkBoundary> build_sequence_chunks(size_t seq_len,
                                                        size_t target_bases)
{
    std::vector<ChunkBoundary> chunks;
    if (seq_len == 0) return chunks;

    const size_t interval = std::max<size_t>(1, target_bases);
    for (size_t start = 0; start < seq_len; start += interval) {
        size_t end = std::min(seq_len, start + interval);
        chunks.push_back(ChunkBoundary{
            static_cast<uint32_t>(start),
            static_cast<uint32_t>(end - start),
            static_cast<uint64_t>(start),
        });
    }
    return chunks;
}

static std::string reconstruct_fasta_from_seq(const std::string& seq,
                                              const uint32_t* contig_ends,
                                              uint32_t n_contigs,
                                              const char* headers,
                                              uint32_t headers_len)
{
    std::vector<std::string_view> header_list;
    {
        size_t start = 0;
        for (size_t i = 0; i <= headers_len; ++i) {
            if (i == headers_len || headers[i] == '\0') {
                if (i > start)
                    header_list.emplace_back(headers + start, i - start);
                start = i + 1;
            }
        }
    }

    std::string fasta;
    fasta.reserve(seq.size() + static_cast<size_t>(n_contigs) * 64);
    uint32_t prev_end = 0;
    for (uint32_t ci = 0; ci < n_contigs; ++ci) {
        uint32_t this_end = contig_ends[ci];
        fasta += '>';
        if (ci < header_list.size())
            fasta.append(header_list[ci].data(), header_list[ci].size());
        fasta += '\n';
        fasta.append(seq.data() + prev_end, this_end - prev_end);
        fasta += '\n';
        prev_end = this_end;
    }
    return fasta;
}

static SequenceBlobLayout build_sequence_blob_layout(const char* fasta, size_t len,
                                                     bool use_2bit = false)
{
    auto fc = extract_fasta_components(fasta, len);

    SequenceBlobLayout layout;
    layout.seq = std::move(fc.seq);
    layout.contig_ends = std::move(fc.contig_ends);
    for (const auto& h : fc.headers) {
        layout.headers_flat += h;
        layout.headers_flat.push_back('\0');
    }

    // Attempt 2-bit packing if requested
    bool packed = false;
    std::vector<uint8_t> packed_bytes;
    std::vector<AmbigEntry> ambig_entries;
    if (use_2bit) {
        packed = pack_sequence_2bit(layout.seq.data(), layout.seq.size(),
                                    packed_bytes, ambig_entries);
    }

    SequenceBlobHeader hdr{};
    hdr.magic = SEQB_MAGIC;
    hdr.version = 1;
    hdr.flags = packed ? SEQB_FLAG_PACKED_2BIT : 0;
    hdr.seq_len = static_cast<uint32_t>(layout.seq.size()); // always original base count
    hdr.n_contigs = static_cast<uint32_t>(layout.contig_ends.size());
    hdr.headers_len = static_cast<uint32_t>(layout.headers_flat.size());
    hdr.reserved = packed ? static_cast<uint32_t>(ambig_entries.size()) : 0;
    hdr.contig_ends_offset = sizeof(SequenceBlobHeader);
    hdr.headers_offset = hdr.contig_ends_offset +
        static_cast<uint64_t>(layout.contig_ends.size()) * sizeof(uint32_t);
    // Ambiguity sidecar (when packed) is placed after headers in the uncompressed prefix
    uint64_t sidecar_size = packed
        ? static_cast<uint64_t>(ambig_entries.size()) * sizeof(AmbigEntry)
        : 0;
    hdr.blocks_offset = hdr.headers_offset + layout.headers_flat.size() + sidecar_size;

    layout.prefix.reserve(static_cast<size_t>(hdr.blocks_offset));
    auto push_prefix = [&](const void* data, size_t size) {
        const auto* p = static_cast<const uint8_t*>(data);
        layout.prefix.insert(layout.prefix.end(), p, p + size);
    };
    push_prefix(&hdr, sizeof(hdr));
    if (!layout.contig_ends.empty())
        push_prefix(layout.contig_ends.data(),
                    layout.contig_ends.size() * sizeof(uint32_t));
    if (!layout.headers_flat.empty())
        push_prefix(layout.headers_flat.data(), layout.headers_flat.size());
    if (packed && !ambig_entries.empty())
        push_prefix(ambig_entries.data(), ambig_entries.size() * sizeof(AmbigEntry));

    layout.blocks_offset = static_cast<uint32_t>(hdr.blocks_offset);

    // Replace seq with packed bytes so compress_range operates on packed data
    if (packed)
        layout.seq.assign(reinterpret_cast<const char*>(packed_bytes.data()), packed_bytes.size());

    return layout;
}

static ParsedSequenceBlob parse_sequence_blob(const uint8_t* blob, size_t blob_len)
{
    if (blob_len < sizeof(SequenceBlobHeader))
        throw std::runtime_error("ShardReader: sequence blob too small");

    ParsedSequenceBlob parsed;
    std::memcpy(&parsed.hdr, blob, sizeof(SequenceBlobHeader));
    if (parsed.hdr.magic != SEQB_MAGIC || parsed.hdr.version != 1)
        throw std::runtime_error("ShardReader: bad sequence blob header");
    if (parsed.hdr.contig_ends_offset +
            static_cast<uint64_t>(parsed.hdr.n_contigs) * sizeof(uint32_t) > blob_len ||
        parsed.hdr.headers_offset + parsed.hdr.headers_len > blob_len ||
        parsed.hdr.blocks_offset > blob_len)
        throw std::runtime_error("ShardReader: sequence blob overflow");

    parsed.contig_ends = parsed.hdr.n_contigs == 0
        ? nullptr
        : reinterpret_cast<const uint32_t*>(blob + parsed.hdr.contig_ends_offset);
    parsed.headers = parsed.hdr.headers_len == 0
        ? nullptr
        : reinterpret_cast<const char*>(blob + parsed.hdr.headers_offset);
    parsed.blocks = blob + parsed.hdr.blocks_offset;
    return parsed;
}

static std::string extract_sequence_slice_from_fasta_block(std::string_view fasta,
                                                           uint64_t start,
                                                           uint64_t length)
{
    if (length == 0) return {};

    std::string seq;
    seq.reserve(static_cast<size_t>(length));

    uint64_t seq_pos = 0;
    size_t i = 0;
    while (i < fasta.size() && seq.size() < length) {
        if (fasta[i] == '>') {
            while (i < fasta.size() && fasta[i] != '\n' && fasta[i] != '\r') ++i;
            while (i < fasta.size() && (fasta[i] == '\n' || fasta[i] == '\r')) ++i;
            continue;
        }
        char c = fasta[i++];
        if (c == '\n' || c == '\r')
            continue;
        if (seq_pos >= start)
            seq.push_back(c);
        ++seq_pos;
    }
    return seq;
}

static PackedCheckpointLayout unpack_checkpoint_layout(const CheckpointEntry& cp)
{
    PackedCheckpointLayout layout;
    layout.seq_byte_offset = static_cast<uint16_t>(cp._pad & 0xffffu);
    layout.line_bases = static_cast<uint16_t>(cp._pad >> 16);
    return layout;
}

static bool append_packed_sequence_slice(std::string& out,
                                         std::string_view fasta,
                                         const CheckpointEntry& cp,
                                         uint64_t local_start,
                                         uint64_t length)
{
    std::string tmp;
    tmp.reserve(static_cast<size_t>(length));

    PackedCheckpointLayout layout = unpack_checkpoint_layout(cp);
    if (layout.line_bases == 0 || layout.seq_byte_offset > fasta.size())
        return false;

    uint64_t byte_pos = layout.seq_byte_offset + local_start +
                        (local_start / layout.line_bases);
    size_t probe = static_cast<size_t>(layout.seq_byte_offset) + layout.line_bases;
    uint64_t newline_bytes = 0;
    if (probe < fasta.size()) {
        if (fasta[probe] == '\r' && probe + 1 < fasta.size() && fasta[probe + 1] == '\n')
            newline_bytes = 2;
        else if (fasta[probe] == '\n' || fasta[probe] == '\r')
            newline_bytes = 1;
    }

    if (byte_pos > fasta.size())
        return false;

    uint64_t remaining = length;
    uint64_t in_line = local_start % layout.line_bases;
    while (remaining > 0 && byte_pos <= fasta.size()) {
        uint64_t take = std::min<uint64_t>(layout.line_bases - in_line, remaining);
        if (byte_pos + take > fasta.size())
            take = fasta.size() - byte_pos;
        tmp.append(fasta.data() + byte_pos, static_cast<size_t>(take));
        remaining -= take;
        byte_pos += take;
        if (remaining == 0 || byte_pos >= fasta.size())
            break;
        if (newline_bytes > 0) {
            if (byte_pos + newline_bytes > fasta.size())
                break;
            byte_pos += newline_bytes;
        }
        in_line = 0;
    }
    if (remaining != 0)
        return false;
    out += tmp;
    return true;
}

static std::vector<size_t> choose_sample_indices(size_t population,
                                                 size_t sample_count,
                                                 uint64_t seed)
{
    std::vector<size_t> indices(population);
    std::iota(indices.begin(), indices.end(), 0);
    std::mt19937_64 rng(seed);
    std::shuffle(indices.begin(), indices.end(), rng);
    if (indices.size() > sample_count)
        indices.resize(sample_count);
    std::sort(indices.begin(), indices.end());
    return indices;
}

static std::vector<size_t> choose_spaced_indices(size_t population,
                                                 size_t sample_count)
{
    std::vector<size_t> out;
    if (population == 0 || sample_count == 0)
        return out;
    if (population <= sample_count) {
        out.resize(population);
        std::iota(out.begin(), out.end(), 0);
        return out;
    }
    if (sample_count == 1) {
        out.push_back(population / 2);
        return out;
    }
    out.reserve(sample_count);
    for (size_t i = 0; i < sample_count; ++i)
        out.push_back((i * (population - 1)) / (sample_count - 1));
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

static size_t effective_checkpoint_bases(size_t base_bases, bool delta_mode)
{
    if (!delta_mode) return base_bases;
    size_t boosted = (std::max<size_t>(base_bases, 65536) * 3) / 2;
    return std::min<size_t>(boosted, 196608);
}

} // namespace

// ── ShardWriter::Impl ───────────────────────────────────────────────────────

struct ShardWriter::Impl {
    uint32_t shard_id;
    uint32_t cluster_id;
    Config   cfg;

    struct PendingGenome {
        GenomeId genome_id;
        uint64_t oph_fingerprint;
        uint32_t flags;
        uint32_t meta_row_id = 0;
        uint64_t raw_offset;   // byte offset into raw_buffer
        uint32_t raw_len;
    };
    std::vector<PendingGenome> pending;

    std::vector<char> raw_buffer;     // append-only store of raw FASTAs
    uint64_t          total_raw_bytes = 0;

    std::vector<std::string> dict_samples;
    std::string              shared_dict;

    Impl(uint32_t sid, uint32_t cid, Config c)
        : shard_id(sid), cluster_id(cid), cfg(c)
    {
        raw_buffer.reserve(c.max_shard_size_bytes);
    }
};

ShardWriter::ShardWriter(uint32_t shard_id, uint32_t cluster_id, Config cfg)
    : impl_(std::make_unique<Impl>(shard_id, cluster_id, cfg))
{}

ShardWriter::~ShardWriter() = default;

void ShardWriter::add_genome(GenomeId id, uint64_t oph_fingerprint,
                                const char* fasta_data, size_t fasta_len,
                                uint32_t flags)
{
    static constexpr size_t MAX_SAMPLE = 65536;

    if (impl_->cfg.train_dict &&
        impl_->dict_samples.size() < impl_->cfg.dict_samples) {
        impl_->dict_samples.emplace_back(fasta_data, std::min(fasta_len, MAX_SAMPLE));
    }

    Impl::PendingGenome pg;
    pg.genome_id       = id;
    pg.oph_fingerprint = oph_fingerprint;
    pg.flags           = flags;
    pg.meta_row_id     = 0;
    pg.raw_offset      = static_cast<uint64_t>(impl_->raw_buffer.size());
    pg.raw_len         = static_cast<uint32_t>(fasta_len);
    impl_->raw_buffer.insert(impl_->raw_buffer.end(), fasta_data, fasta_data + fasta_len);
    impl_->total_raw_bytes += fasta_len;
    impl_->pending.push_back(std::move(pg));
}

FrozenShard ShardWriter::freeze() {
    const size_t n = impl_->pending.size();
    const bool allow_auto_codec = impl_->cfg.auto_codec;
    const bool try_mem_delta = impl_->cfg.use_mem_delta;

    // ── 0. MEM-delta path (codec=4): position-independent k-mer matching ──────
    if (try_mem_delta && n >= 1) {
        struct CompressedBlob {
            std::vector<char> data;
            uint32_t raw_len = 0;
            std::vector<CheckpointEntry> checkpoints;
            bool use_mem_delta = false;
        };
        std::vector<CompressedBlob> blobs(n);

        const size_t n_ref = std::min(impl_->cfg.mem_delta_ref_panel, n);
        std::vector<size_t> panel_indices;
        panel_indices.reserve(n_ref);
        if (n_ref == 1) {
            panel_indices.push_back(n / 2);
        } else {
            for (size_t i = 0; i < n_ref; ++i)
                panel_indices.push_back((i * (n - 1)) / (n_ref - 1));
        }

        std::vector<uint8_t> is_panel(n, 0);
        for (size_t idx : panel_indices)
            is_panel[idx] = 1;

        std::vector<size_t> order;
        order.reserve(n);
        for (size_t idx : panel_indices)
            order.push_back(idx);
        for (size_t i = 0; i < n; ++i)
            if (!is_panel[i]) order.push_back(i);

        auto compress_plain_blob = [&](const Impl::PendingGenome& pg, CompressedBlob& out) {
            ZSTD_CCtx* cctx = ZSTD_createCCtx();
            if (!cctx) throw std::runtime_error("ZSTD_createCCtx failed");
            ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, impl_->cfg.zstd_level);
            if (impl_->cfg.use_long_match) {
                ZSTD_CCtx_setParameter(cctx, ZSTD_c_enableLongDistanceMatching, 1);
                ZSTD_CCtx_setParameter(cctx, ZSTD_c_windowLog, impl_->cfg.zstd_wlog);
            }

            const char* src = impl_->raw_buffer.data() + pg.raw_offset;
            auto blob_layout = build_sequence_blob_layout(src, pg.raw_len,
                                                          impl_->cfg.use_2bit_pack);
            out.raw_len = pg.raw_len;
            out.data.clear();
            out.checkpoints.clear();
            out.data.insert(out.data.end(), blob_layout.prefix.begin(), blob_layout.prefix.end());

            auto ranges = build_sequence_chunks(blob_layout.seq.size(), impl_->cfg.checkpoint_bases);
            out.checkpoints.reserve(ranges.size());

            uint32_t block_offset = blob_layout.blocks_offset;
            for (const auto& range : ranges) {
                size_t bound = ZSTD_compressBound(range.byte_len);
                size_t old_size = out.data.size();
                out.data.resize(old_size + bound);
                size_t csize = ZSTD_compress2(cctx,
                                              out.data.data() + old_size,
                                              bound,
                                              blob_layout.seq.data() + range.byte_start,
                                              range.byte_len);
                if (ZSTD_isError(csize)) {
                    ZSTD_freeCCtx(cctx);
                    throw std::runtime_error(std::string("ZSTD_compress2: ") +
                                             ZSTD_getErrorName(csize));
                }
                out.data.resize(old_size + csize);
                out.checkpoints.push_back(CheckpointEntry{
                    range.symbol_offset,
                    block_offset,
                    0,
                });
                block_offset += static_cast<uint32_t>(csize);
            }

            ZSTD_freeCCtx(cctx);
        };

        for (size_t i = 0; i < n_ref; ++i) {
            const auto& pg = impl_->pending[order[i]];
            compress_plain_blob(pg, blobs[i]);
        }

        FastaComponents panel_fc;
        for (size_t i = 0; i < n_ref; ++i) {
            const auto& pg = impl_->pending[order[i]];
            auto fc = extract_fasta_components(impl_->raw_buffer.data() + pg.raw_offset, pg.raw_len);
            panel_fc.seq += fc.seq;
        }
        AnchorIndex anchor_idx = build_anchor_index(panel_fc.seq);

        const size_t n_delta   = n - n_ref;
        const size_t n_threads = impl_->cfg.compress_threads > 0
            ? impl_->cfg.compress_threads
            : std::min(size_t(16), std::max(size_t(1),
                static_cast<size_t>(std::thread::hardware_concurrency())));
        const size_t actual_threads = std::min(n_threads, n_delta > 0 ? n_delta : size_t(1));
        const size_t chunk = n_delta > 0
            ? (n_delta + actual_threads - 1) / actual_threads : 0;

        auto encode_range = [&](size_t start, size_t end_idx) {
            for (size_t i = start; i < end_idx; ++i) {
                const size_t out_idx = n_ref + i;
                const auto& pg = impl_->pending[order[out_idx]];
                const char* src = impl_->raw_buffer.data() + pg.raw_offset;
                FastaComponents qfc = extract_fasta_components(src, pg.raw_len);

                auto delta = encode_mem_delta(panel_fc, anchor_idx, qfc, impl_->cfg.checkpoint_bases);

                CompressedBlob plain;
                compress_plain_blob(pg, plain);
                size_t plain_cost = plain.data.size() +
                    plain.checkpoints.size() * sizeof(CheckpointEntry);

                if (!delta.empty() && delta.size() * 100 <= plain_cost * 70) {
                    blobs[out_idx].raw_len = pg.raw_len;
                    blobs[out_idx].data.assign(reinterpret_cast<const char*>(delta.data()),
                                               reinterpret_cast<const char*>(delta.data()) + delta.size());
                    blobs[out_idx].checkpoints.clear();
                    blobs[out_idx].use_mem_delta = true;
                } else {
                    blobs[out_idx] = std::move(plain);
                    blobs[out_idx].use_mem_delta = false;
                }
            }
        };

        if (n_delta == 0) {
            // all genomes are in the reference panel; nothing to delta-encode
        } else if (actual_threads == 1 || n_delta == 1) {
            encode_range(0, n_delta);
        } else {
            std::vector<std::future<void>> futs;
            futs.reserve(actual_threads);
            for (size_t t = 0; t < actual_threads; ++t) {
                size_t s = t * chunk;
                size_t e = std::min(s + chunk, n_delta);
                if (s >= n_delta) break;
                futs.push_back(std::async(std::launch::async, encode_range, s, e));
            }
            for (auto& f : futs) f.get();
        }

        const uint64_t header_size       = sizeof(ShardHeader);
        const uint64_t dir_size          = n * sizeof(GenomeDirEntry);
        const uint64_t genome_dir_offset = header_size;
        const uint64_t blob_area_offset  = genome_dir_offset + dir_size;

        std::vector<GenomeDirEntry> dir(n);
        std::vector<CheckpointEntry> checkpoint_entries;
        uint64_t blob_cursor = 0;
        uint64_t total_compressed = 0;
        size_t mem_delta_count = 0;
        size_t plain_count = 0;
        for (size_t i = 0; i < n; ++i) {
            const auto& pg = impl_->pending[order[i]];
            const auto& cb = blobs[i];
            dir[i].genome_id       = pg.genome_id;
            dir[i].oph_fingerprint = pg.oph_fingerprint;
            dir[i].blob_offset     = blob_cursor;
            dir[i].blob_len_cmp    = static_cast<uint32_t>(cb.data.size());
            dir[i].blob_len_raw    = cb.raw_len;
            dir[i].checkpoint_idx  = static_cast<uint32_t>(checkpoint_entries.size());
            dir[i].n_checkpoints   = static_cast<uint32_t>(cb.checkpoints.size());
            dir[i].flags           = pg.flags | (cb.use_mem_delta ? GenomeMeta::FLAG_DELTA : 0u);
            dir[i].meta_row_id     = pg.meta_row_id;
            std::memset(dir[i].reserved, 0, sizeof(dir[i].reserved));
            blob_cursor      += cb.data.size();
            total_compressed += cb.data.size();
            checkpoint_entries.insert(checkpoint_entries.end(),
                                      cb.checkpoints.begin(),
                                      cb.checkpoints.end());
            if (cb.use_mem_delta) ++mem_delta_count;
            else ++plain_count;
        }
        const uint64_t checkpoint_area_offset = blob_area_offset + blob_cursor;

        FrozenShard frozen;
        frozen.shard_id  = impl_->shard_id;
        frozen.n_genomes = static_cast<uint32_t>(n);
        frozen.raw_bytes = impl_->total_raw_bytes;

        const uint64_t total_section_size = checkpoint_area_offset +
            checkpoint_entries.size() * sizeof(CheckpointEntry);
        frozen.bytes.reserve(total_section_size);
        auto push_bytes = [&](const void* d, size_t sz) {
            const auto* p = static_cast<const uint8_t*>(d);
            frozen.bytes.insert(frozen.bytes.end(), p, p + sz);
        };

        ShardHeader hdr{};
        hdr.magic                  = GPKS_MAGIC;
        hdr.version                = 4;
        hdr.flags                  = 0;
        hdr.shard_id               = impl_->shard_id;
        hdr.cluster_id             = impl_->cluster_id;
        hdr.n_genomes              = static_cast<uint32_t>(n);
        hdr.n_deleted              = static_cast<uint32_t>(
            std::count_if(impl_->pending.begin(), impl_->pending.end(),
                          [](const Impl::PendingGenome& pg){ return pg.flags & GenomeMeta::FLAG_DELETED; }));
        hdr.codec                  = 4u;
        hdr.dict_size              = static_cast<uint32_t>(n_ref);
        hdr.genome_dir_offset      = genome_dir_offset;
        hdr.dict_offset            = blob_area_offset;
        hdr.blob_area_offset       = blob_area_offset;
        hdr.checkpoint_area_offset = checkpoint_entries.empty() ? 0 : checkpoint_area_offset;
        hdr.checkpoint_count       = checkpoint_entries.size();
        hdr.shard_raw_bp           = impl_->total_raw_bytes;
        hdr.shard_compressed_bytes = total_compressed;
        std::memset(hdr.checksum, 0, sizeof(hdr.checksum));
        std::memset(hdr.reserved, 0, sizeof(hdr.reserved));
        push_bytes(&hdr, sizeof(hdr));
        push_bytes(dir.data(), dir_size);
        for (const auto& cb : blobs)
            push_bytes(cb.data.data(), cb.data.size());
        if (!checkpoint_entries.empty())
            push_bytes(checkpoint_entries.data(),
                       checkpoint_entries.size() * sizeof(CheckpointEntry));

        spdlog::info("ShardWriter shard {} (MEM-delta): {} genomes ({} ref panel), raw {}B, mem-delta {}, plain {}, bytes {}B",
                     impl_->shard_id, n, n_ref, impl_->total_raw_bytes,
                     mem_delta_count, plain_count, total_compressed);
        return frozen;
    }

    // ── 1. Build dictionary ───────────────────────────────────────────────────
    // 2-bit packing is incompatible with delta compression: the delta prefix is the
    // packed byte stream, but the decompressor unpacks the reference to ASCII before
    // using it as a prefix, causing a byte-stream mismatch. Disable delta when 2-bit.
    bool use_delta_mode = impl_->cfg.use_delta && !impl_->cfg.use_2bit_pack;
    bool use_dict     = false;
    bool use_ref_dict = false;
    uint32_t delta_ref_index = 0;
    SequenceBlobLayout delta_ref_blob;

    // Single CCtx shared across all estimate_chunked_cost calls — avoids
    // malloc/free overhead for up to 4×4 candidate × target combinations.
    ZSTD_CCtx* estimate_cctx = ZSTD_createCCtx();
    if (!estimate_cctx) throw std::runtime_error("ZSTD_createCCtx failed (estimate)");
    struct EstimateCCtxGuard {
        ZSTD_CCtx*& p;
        ~EstimateCCtxGuard() { ZSTD_freeCCtx(p); }
    } estimate_cctx_guard{estimate_cctx};

    // Reusable scratch buffer for codec sampling — grows to max chunk bound and stays.
    std::vector<char> estimate_scratch;

    auto estimate_chunked_cost = [&](const Impl::PendingGenome& pg,
                                     bool use_prefix,
                                     std::string_view ref_seq) {
            // Reset session + parameters so each call starts clean.
            ZSTD_CCtx_reset(estimate_cctx, ZSTD_reset_session_and_parameters);
            ZSTD_CCtx_setParameter(estimate_cctx, ZSTD_c_compressionLevel, impl_->cfg.zstd_level);
            if (impl_->cfg.use_long_match) {
                ZSTD_CCtx_setParameter(estimate_cctx, ZSTD_c_enableLongDistanceMatching, 1);
                ZSTD_CCtx_setParameter(estimate_cctx, ZSTD_c_windowLog, impl_->cfg.zstd_wlog);
            }

            const char* src = impl_->raw_buffer.data() + pg.raw_offset;
            auto blob_layout = build_sequence_blob_layout(src, pg.raw_len,
                                                          impl_->cfg.use_2bit_pack);
            const char* ref_ptr = ref_seq.data();
            size_t ref_len = ref_seq.size();
            size_t target_bases = effective_checkpoint_bases(impl_->cfg.checkpoint_bases, use_prefix);
            if (use_prefix) {
                int wlog = 23;
                for (size_t sz = 1u << wlog; sz < target_bases * 2; sz <<= 1) ++wlog;
                ZSTD_CCtx_setParameter(estimate_cctx, ZSTD_c_windowLog, wlog);
            }

            size_t total = 0;
            auto ranges = build_sequence_chunks(blob_layout.seq.size(), target_bases);
            auto sampled_ranges = choose_spaced_indices(ranges.size(), std::min<size_t>(3, ranges.size()));
            size_t max_bound = 0;
            for (size_t ridx : sampled_ranges)
                max_bound = std::max(max_bound, ZSTD_compressBound(ranges[ridx].byte_len));
            if (estimate_scratch.size() < max_bound)
                estimate_scratch.resize(max_bound);
            for (size_t ridx : sampled_ranges) {
                const auto& range = ranges[ridx];
                if (use_prefix) {
                    size_t ref_chunk_start = std::min<uint64_t>(range.symbol_offset, ref_len);
                    size_t ref_chunk_len = ref_chunk_start < ref_len
                        ? std::min<size_t>(range.byte_len, ref_len - ref_chunk_start)
                        : 0;
                    if (ref_chunk_len > 0)
                        ZSTD_CCtx_refPrefix(estimate_cctx, ref_ptr + ref_chunk_start, ref_chunk_len);
                }
                size_t bound = ZSTD_compressBound(range.byte_len);
                size_t csize = ZSTD_compress2(estimate_cctx, estimate_scratch.data(), bound,
                                              blob_layout.seq.data() + range.byte_start,
                                              range.byte_len);
                if (ZSTD_isError(csize))
                    throw std::runtime_error(std::string("ZSTD_compress2: ") + ZSTD_getErrorName(csize));
                total += csize;
            }
            if (!sampled_ranges.empty())
                total = (total * ranges.size() + sampled_ranges.size() - 1) / sampled_ranges.size();
            return total + ranges.size() * sizeof(CheckpointEntry);
    };

    if (use_delta_mode && n >= 1) {
        uint64_t seed = (static_cast<uint64_t>(impl_->shard_id) << 32) ^
            static_cast<uint64_t>(impl_->total_raw_bytes) ^ 0x44454c5441415554ULL;
        if (n >= 2) {
            auto candidate_indices = choose_sample_indices(n, std::min<size_t>(4, n), seed ^ 0xA11CE55EEDULL);
            auto target_indices = choose_sample_indices(n, std::min<size_t>(4, n), seed ^ 0xD371A5EULL);

            size_t best_delta_cost = std::numeric_limits<size_t>::max();
            for (size_t cand_idx : candidate_indices) {
                auto candidate_blob = build_sequence_blob_layout(
                    impl_->raw_buffer.data() + impl_->pending[cand_idx].raw_offset,
                    impl_->pending[cand_idx].raw_len);
                size_t total_cost = 0;
                size_t compared = 0;
                for (size_t target_idx : target_indices) {
                    if (target_idx == cand_idx) continue;
                    total_cost += estimate_chunked_cost(impl_->pending[target_idx], true, candidate_blob.seq);
                    ++compared;
                }
                if (compared == 0) continue;
                if (total_cost < best_delta_cost) {
                    best_delta_cost = total_cost;
                    delta_ref_index = static_cast<uint32_t>(cand_idx);
                    delta_ref_blob = std::move(candidate_blob);
                }
            }
            if (delta_ref_blob.seq.empty()) {
                delta_ref_index = 0;
                delta_ref_blob = build_sequence_blob_layout(
                    impl_->raw_buffer.data() + impl_->pending[0].raw_offset,
                    impl_->pending[0].raw_len);
            }

            if (allow_auto_codec) {
                size_t delta_wins = 0;
                size_t plain_cost_sum = 0;
                size_t delta_cost_sum = 0;
                size_t compared = 0;
                for (size_t target_idx : target_indices) {
                    if (target_idx == delta_ref_index) continue;
                    const auto& pg = impl_->pending[target_idx];
                    size_t plain_cost = estimate_chunked_cost(pg, false, {});
                    size_t delta_cost = estimate_chunked_cost(pg, true, delta_ref_blob.seq);
                    plain_cost_sum += plain_cost;
                    delta_cost_sum += delta_cost;
                    if (delta_cost < plain_cost) ++delta_wins;
                    ++compared;
                }

                if (compared == 0 ||
                    delta_wins * 2 < compared ||
                    delta_cost_sum * 100 > plain_cost_sum * 97) {
                    use_delta_mode = false;
                    spdlog::debug("ShardWriter shard {}: auto codec chose plain zstd over delta (wins {}/{})",
                                  impl_->shard_id, delta_wins, compared);
                }
            }
        } else {
            delta_ref_blob = build_sequence_blob_layout(
                impl_->raw_buffer.data() + impl_->pending[0].raw_offset,
                impl_->pending[0].raw_len);
        }
    }

    if (impl_->cfg.use_reference_dict && n >= 1) {
        // Reference content dictionary: use the first genome as raw dictionary
        // content via ZDICT_finalizeDictionary. Genomes with high sequence
        // similarity to the reference compress much smaller; decompression is
        // faster because blobs are smaller and the DDict is loaded once per
        // shard open.
        const auto& ref     = impl_->pending[0];
        const char* ref_ptr = impl_->raw_buffer.data() + ref.raw_offset;
        const size_t ref_len = ref.raw_len;

        // Target dict size: embed as much of the reference genome as possible.
        const size_t target = std::min(ref_len + 4096, impl_->cfg.ref_dict_max_size);
        impl_->shared_dict.resize(target);

        // Use the reference genome as both forced content and the training sample.
        // ZDICT_finalizeDictionary embeds up to `target` bytes from the END of
        // dictContent (the recent sequence context zstd expects) plus hash tables.
        const size_t sample_size = std::min(ref_len, size_t(65536));
        ZDICT_params_t params{};
        params.compressionLevel = impl_->cfg.zstd_level;

        size_t dict_bytes = ZDICT_finalizeDictionary(
            impl_->shared_dict.data(), target,
            ref_ptr, ref_len,           // forced content = full reference genome
            ref_ptr, &sample_size, 1,   // single training sample
            params);

        if (ZDICT_isError(dict_bytes)) {
            spdlog::warn("ShardWriter shard {}: reference dict build failed: {} — no dict",
                         impl_->shard_id, ZDICT_getErrorName(dict_bytes));
            impl_->shared_dict.clear();
        } else {
            impl_->shared_dict.resize(dict_bytes);
            spdlog::debug("ShardWriter shard {}: reference dict {}B from {}B reference genome",
                          impl_->shard_id, dict_bytes, ref_len);
            use_dict     = true;
            use_ref_dict = true;
        }
    } else if (impl_->cfg.train_dict && impl_->dict_samples.size() >= 8) {
        std::string concat;
        std::vector<size_t> sizes;
        concat.reserve(impl_->dict_samples.size() * 32768);
        for (const auto& s : impl_->dict_samples) {
            sizes.push_back(s.size());
            concat += s;
        }
        impl_->shared_dict.resize(impl_->cfg.dict_size);
        size_t trained = ZDICT_trainFromBuffer(
            impl_->shared_dict.data(), impl_->cfg.dict_size,
            concat.data(), sizes.data(), sizes.size());
        if (ZDICT_isError(trained)) {
            spdlog::warn("ShardWriter shard {}: ZDICT training failed: {} — no dict",
                         impl_->shard_id, ZDICT_getErrorName(trained));
            impl_->shared_dict.clear();
        } else {
            impl_->shared_dict.resize(trained);
            spdlog::debug("ShardWriter shard {}: trained {}B dict from {} samples",
                          impl_->shard_id, trained, impl_->dict_samples.size());
            use_dict = true;
        }
    }

    // ── 2. Build shared CDict (immutable, shared across threads) ─────────────
    ZSTD_CDict* cdict = nullptr;
    if (use_dict) {
        cdict = ZSTD_createCDict(impl_->shared_dict.data(), impl_->shared_dict.size(),
                                 impl_->cfg.zstd_level);
        if (!cdict) throw std::runtime_error("ZSTD_createCDict failed");
    }

    // ── 3. Compress all pending genomes in parallel ───────────────────────────
    struct CompressedBlob {
        std::vector<char> data;
        uint32_t          raw_len;
        std::vector<CheckpointEntry> checkpoints;
    };
    std::vector<CompressedBlob> blobs(n);

    const size_t n_threads = impl_->cfg.compress_threads > 0
        ? impl_->cfg.compress_threads
        : std::min(size_t(16), std::max(size_t(1), static_cast<size_t>(std::thread::hardware_concurrency())));
    const size_t actual_threads = std::min(n_threads, n);
    const size_t chunk = (n + actual_threads - 1) / actual_threads;

    // For delta mode: reference sequence bytes are the prefix for all non-reference blobs.
    const char* delta_ref_ptr = use_delta_mode && n >= 1
        ? delta_ref_blob.seq.data() : nullptr;
    const size_t delta_ref_len = use_delta_mode && n >= 1
        ? delta_ref_blob.seq.size() : 0;
    const size_t delta_target_bases = effective_checkpoint_bases(impl_->cfg.checkpoint_bases, use_delta_mode);
    int delta_window_log = 23;
    if (delta_ref_ptr) {
        for (size_t sz = 1u << delta_window_log; sz < delta_target_bases * 2; sz <<= 1)
            ++delta_window_log;
    }

    auto compress_range = [&](size_t start, size_t end) {
        ZSTD_CCtx* cctx = ZSTD_createCCtx();
        if (!cctx) throw std::runtime_error("ZSTD_createCCtx failed");
        ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, impl_->cfg.zstd_level);
        if (impl_->cfg.use_long_match) {
            ZSTD_CCtx_setParameter(cctx, ZSTD_c_enableLongDistanceMatching, 1);
            ZSTD_CCtx_setParameter(cctx, ZSTD_c_windowLog, impl_->cfg.zstd_wlog);
        }
        if (delta_ref_ptr) {
            ZSTD_CCtx_setParameter(cctx, ZSTD_c_windowLog, delta_window_log);
        }
        if (cdict)
            ZSTD_CCtx_refCDict(cctx, cdict);
        for (size_t i = start; i < end; ++i) {
            const auto& pg  = impl_->pending[i];
            const char* src = impl_->raw_buffer.data() + pg.raw_offset;
            auto blob_layout = build_sequence_blob_layout(src, pg.raw_len,
                                                          impl_->cfg.use_2bit_pack);
            blobs[i].raw_len = pg.raw_len;

            auto ranges = build_sequence_chunks(blob_layout.seq.size(), delta_target_bases);
            blobs[i].checkpoints.clear();
            blobs[i].data.clear();
            blobs[i].checkpoints.reserve(ranges.size());
            blobs[i].data.reserve(blob_layout.prefix.size() + blob_layout.seq.size());
            blobs[i].data.insert(blobs[i].data.end(),
                                 blob_layout.prefix.begin(),
                                 blob_layout.prefix.end());

            uint32_t block_offset = blob_layout.blocks_offset;
            for (const auto& range : ranges) {
                // Delta: non-reference chunks set the reference genome as prefix.
                // ZSTD_CCtx_refPrefix is consumed by the next compress2 call.
                if (delta_ref_ptr && i != delta_ref_index) {
                    size_t ref_chunk_start = std::min<uint64_t>(range.symbol_offset, delta_ref_len);
                    size_t ref_chunk_len = ref_chunk_start < delta_ref_len
                        ? std::min<size_t>(range.byte_len, delta_ref_len - ref_chunk_start)
                        : 0;
                    if (ref_chunk_len > 0)
                        ZSTD_CCtx_refPrefix(cctx,
                                            delta_ref_ptr + ref_chunk_start,
                                            ref_chunk_len);
                }

                // Write compressed output directly into blobs[i].data — no scratch
                // buffer or memcpy needed.
                size_t bound = ZSTD_compressBound(range.byte_len);
                size_t old_size = blobs[i].data.size();
                blobs[i].data.resize(old_size + bound);
                size_t csize = ZSTD_compress2(cctx,
                                              blobs[i].data.data() + old_size,
                                              bound,
                                              blob_layout.seq.data() + range.byte_start,
                                              range.byte_len);
                if (ZSTD_isError(csize)) {
                    ZSTD_freeCCtx(cctx);
                    throw std::runtime_error(std::string("ZSTD_compress2: ") + ZSTD_getErrorName(csize));
                }
                blobs[i].data.resize(old_size + csize);
                blobs[i].checkpoints.push_back(CheckpointEntry{
                    range.symbol_offset,
                    block_offset,
                    0,
                });
                block_offset += static_cast<uint32_t>(csize);
            }
        }
        ZSTD_freeCCtx(cctx);
    };

    if (actual_threads == 1) {
        compress_range(0, n);
    } else {
        std::vector<std::future<void>> compress_futs;
        compress_futs.reserve(actual_threads);
        for (size_t t = 0; t < actual_threads; ++t) {
            size_t s = t * chunk;
            size_t e = std::min(s + chunk, n);
            if (s >= n) break;
            compress_futs.push_back(std::async(std::launch::async, compress_range, s, e));
        }
        for (auto& f : compress_futs) f.get();
    }

    if (cdict) ZSTD_freeCDict(cdict);

    // ── 4. Compute layout ─────────────────────────────────────────────────────
    const uint64_t header_size       = sizeof(ShardHeader);
    const uint64_t dir_size          = n * sizeof(GenomeDirEntry);
    const uint64_t dict_bytes        = use_dict ? impl_->shared_dict.size() : 0;
    const uint64_t genome_dir_offset = header_size;
    const uint64_t dict_offset       = genome_dir_offset + dir_size;
    const uint64_t blob_area_offset  = dict_offset + dict_bytes;

    std::vector<GenomeDirEntry> dir(n);
    std::vector<CheckpointEntry> checkpoint_entries;
    uint64_t blob_cursor      = 0;
    uint64_t total_compressed = 0;
    for (size_t i = 0; i < n; ++i) {
        const auto& pg = impl_->pending[i];
        const auto& cb = blobs[i];
        dir[i].genome_id       = pg.genome_id;
        dir[i].oph_fingerprint = pg.oph_fingerprint;
        dir[i].blob_offset     = blob_cursor;
        dir[i].blob_len_cmp    = static_cast<uint32_t>(cb.data.size());
        dir[i].blob_len_raw    = cb.raw_len;
        dir[i].checkpoint_idx  = static_cast<uint32_t>(checkpoint_entries.size());
        dir[i].n_checkpoints   = static_cast<uint32_t>(cb.checkpoints.size());
        dir[i].flags           = pg.flags | (use_delta_mode && i != delta_ref_index ? GenomeMeta::FLAG_DELTA : 0u);
        dir[i].meta_row_id     = pg.meta_row_id;
        std::memset(dir[i].reserved, 0, sizeof(dir[i].reserved));
        blob_cursor      += cb.data.size();
        total_compressed += cb.data.size();
        checkpoint_entries.insert(checkpoint_entries.end(),
                                  cb.checkpoints.begin(),
                                  cb.checkpoints.end());
    }
    const uint64_t checkpoint_area_offset = blob_area_offset + blob_cursor;

    // ── 5. Serialize to bytes ─────────────────────────────────────────────────
    FrozenShard frozen;
    frozen.shard_id  = impl_->shard_id;
    frozen.n_genomes = static_cast<uint32_t>(n);
    frozen.raw_bytes = impl_->total_raw_bytes;

    const uint64_t total_section_size = checkpoint_area_offset +
        checkpoint_entries.size() * sizeof(CheckpointEntry);
    frozen.bytes.reserve(total_section_size);

    auto push_bytes = [&](const void* d, size_t sz) {
        const auto* p = static_cast<const uint8_t*>(d);
        frozen.bytes.insert(frozen.bytes.end(), p, p + sz);
    };

    ShardHeader hdr{};
    hdr.magic                  = GPKS_MAGIC;
    hdr.version                = 4;
    hdr.flags                  = use_delta_mode ? delta_ref_index : 0;
    hdr.shard_id               = impl_->shard_id;
    hdr.cluster_id             = impl_->cluster_id;
    hdr.n_genomes              = static_cast<uint32_t>(n);
    hdr.n_deleted              = static_cast<uint32_t>(
        std::count_if(impl_->pending.begin(), impl_->pending.end(),
                      [](const Impl::PendingGenome& pg){ return pg.flags & GenomeMeta::FLAG_DELETED; }));
    hdr.codec                  = use_delta_mode ? 3u : (use_ref_dict ? 2u : (use_dict ? 1u : 0u));
    hdr.dict_size              = static_cast<uint32_t>(dict_bytes);
    hdr.genome_dir_offset      = genome_dir_offset;
    hdr.dict_offset            = dict_offset;
    hdr.blob_area_offset       = blob_area_offset;
    hdr.checkpoint_area_offset = checkpoint_entries.empty() ? 0 : checkpoint_area_offset;
    hdr.checkpoint_count       = checkpoint_entries.size();
    hdr.shard_raw_bp           = impl_->total_raw_bytes;
    hdr.shard_compressed_bytes = total_compressed;
    std::memset(hdr.checksum, 0, sizeof(hdr.checksum));
    std::memset(hdr.reserved, 0, sizeof(hdr.reserved));
    push_bytes(&hdr, sizeof(hdr));
    push_bytes(dir.data(), dir_size);
    if (use_dict)
        push_bytes(impl_->shared_dict.data(), dict_bytes);
    for (const auto& cb : blobs)
        push_bytes(cb.data.data(), cb.data.size());
    if (!checkpoint_entries.empty())
        push_bytes(checkpoint_entries.data(),
                   checkpoint_entries.size() * sizeof(CheckpointEntry));

    spdlog::info("ShardWriter shard {}: wrote {} genomes, raw {}B, compressed {}B",
                 impl_->shard_id, n, impl_->total_raw_bytes, total_compressed);
    return frozen;
}

uint64_t ShardWriter::finalize(AppendWriter& writer) {
    FrozenShard frozen = freeze();
    uint64_t section_start = writer.current_offset();
    writer.append(frozen.bytes.data(), frozen.bytes.size());
    return section_start;
}

size_t ShardWriter::n_genomes()   const { return impl_->pending.size(); }
size_t ShardWriter::n_bytes_raw() const { return impl_->total_raw_bytes; }

// ── ShardReader::Impl ───────────────────────────────────────────────────────

struct ShardReader::Impl {
    const uint8_t*          base_         = nullptr;  // start of shard section
    uint64_t                section_size_ = 0;
    const ShardHeader*    header_       = nullptr;
    const GenomeDirEntry* dir_          = nullptr;
    const CheckpointEntry* checkpoints_ = nullptr;
    uint64_t               checkpoint_count_ = 0;

    std::vector<uint8_t>    owned_data_;  // non-empty when file-based open

    ZSTD_DDict*  ddict_    = nullptr;   // immutable after setup, safe to share across threads
    mutable std::string  ref_prefix_;   // cached reference prefix for legacy delta shards
    mutable std::mutex   ref_prefix_mx_; // protects ref_prefix_ from concurrent initialization
    std::string  mem_anchor_seq_;      // pure sequence of anchor genome for MEM-delta shards (codec==4)

    // Thread-local ZSTD decompression context with RAII cleanup on thread exit
    struct TlDctx {
        ZSTD_DCtx* ctx;
        TlDctx() : ctx(ZSTD_createDCtx()) {}
        ~TlDctx() { if (ctx) ZSTD_freeDCtx(ctx); }
    };

    static ZSTD_DCtx* get_thread_dctx() {
        static thread_local TlDctx tl;
        return tl.ctx;
    }

    struct CachedFrame {
        const uint8_t* src = nullptr;
        size_t src_size = 0;
        const char* ref_ptr = nullptr;
        size_t ref_size = 0;
        bool use_prefix = false;
        uint64_t age = 0;
        std::string data;
    };

    struct TlFrameCache {
        std::array<CachedFrame, 4> entries{};
        uint64_t next_age = 1;
    };

    static TlFrameCache& get_thread_frame_cache() {
        static thread_local TlFrameCache cache;
        return cache;
    }

    uint32_t delta_ref_index() const {
        if (!header_ || header_->codec != 3 || header_->n_genomes == 0)
            return 0;
        if (header_->version >= 4 && header_->flags < header_->n_genomes)
            return header_->flags;
        return 0;
    }

    void setup(const uint8_t* section_base, uint64_t section_size) {
        base_         = section_base;
        section_size_ = section_size;

        header_ = reinterpret_cast<const ShardHeader*>(base_);
        if (header_->magic != GPKS_MAGIC)
            throw std::runtime_error("ShardReader: invalid shard magic");
        if (header_->version != 2 && header_->version != 3 && header_->version != 4)
            throw std::runtime_error("ShardReader: expected shard version 2, 3, or 4");

        dir_ = reinterpret_cast<const GenomeDirEntry*>(
            base_ + header_->genome_dir_offset);

        if (header_->checkpoint_area_offset > 0 &&
            header_->checkpoint_area_offset <= section_size_) {
            checkpoints_ = reinterpret_cast<const CheckpointEntry*>(
                base_ + header_->checkpoint_area_offset);
            checkpoint_count_ = header_->checkpoint_count;
        } else {
            checkpoints_ = nullptr;
            checkpoint_count_ = 0;
        }

        if (header_->dict_size > 0 && (header_->codec == 1 || header_->codec == 2)) {
            const uint8_t* dict_ptr = base_ + header_->dict_offset;
            ddict_ = ZSTD_createDDict(dict_ptr, header_->dict_size);
            if (!ddict_) throw std::runtime_error("ZSTD_createDDict failed");
        }

        // Delta codec: decompress reference genome (first entry) once and cache it.
        // All delta blobs (FLAG_DELTA) require this prefix at decompression time.
        if (header_->codec == 3 && header_->n_genomes >= 1) {
            if (header_->version < 4)
                ref_prefix_ = decompress_blob(dir_[0]);
        }

        // MEM-delta codec: decompress all reference panel genomes (dict_size = n_ref)
        // and concatenate their pure sequences into mem_anchor_seq_.
        if (header_->codec == 4 && header_->n_genomes >= 1) {
            const uint32_t n_ref = header_->dict_size > 0
                ? std::min(header_->dict_size, header_->n_genomes) : 1u;
            for (uint32_t i = 0; i < n_ref; ++i) {
                std::string ref_fasta = decompress_blob(dir_[i]);
                auto fc = extract_fasta_components(ref_fasta.data(), ref_fasta.size());
                mem_anchor_seq_ += fc.seq;
            }
        }
    }

    void reset() {
        if (ddict_) { ZSTD_freeDDict(ddict_); ddict_ = nullptr; }
        ref_prefix_.clear();
        mem_anchor_seq_.clear();
        base_         = nullptr;
        section_size_ = 0;
        header_       = nullptr;
        dir_          = nullptr;
        checkpoints_  = nullptr;
        checkpoint_count_ = 0;
        owned_data_.clear();
    }

    const GenomeDirEntry* find_genome(GenomeId id) const {
        for (uint32_t i = 0; i < header_->n_genomes; ++i) {
            if (dir_[i].genome_id == id) return &dir_[i];
        }
        return nullptr;
    }

    const GenomeDirEntry* entry_at(uint32_t dir_index) const {
        if (!header_ || dir_index >= header_->n_genomes) return nullptr;
        return &dir_[dir_index];
    }

    const CheckpointEntry* checkpoints_for(const GenomeDirEntry& e) const {
        if (!checkpoints_ || e.n_checkpoints == 0) return nullptr;
        uint64_t end = static_cast<uint64_t>(e.checkpoint_idx) + e.n_checkpoints;
        if (end > checkpoint_count_) return nullptr;
        return checkpoints_ + e.checkpoint_idx;
    }

    const std::string& decompress_zstd_frame_cached_with_ref(const uint8_t* src,
                                                             size_t src_size,
                                                             const char* ref_ptr,
                                                             size_t ref_size) const
    {
        bool use_prefix = ref_ptr != nullptr && ref_size > 0;
        auto& frame_cache = get_thread_frame_cache();
        for (auto& entry : frame_cache.entries) {
            if (entry.src == src &&
                entry.src_size == src_size &&
                entry.ref_ptr == ref_ptr &&
                entry.ref_size == ref_size &&
                entry.use_prefix == use_prefix) {
                entry.age = frame_cache.next_age++;
                return entry.data;
            }
        }

        ZSTD_DCtx* dctx = get_thread_dctx();
        if (ref_ptr)
            ZSTD_DCtx_refPrefix(dctx, ref_ptr, ref_size);

        unsigned long long frame_size = ZSTD_getFrameContentSize(src, src_size);
        if (frame_size == ZSTD_CONTENTSIZE_ERROR)
            throw std::runtime_error("ShardReader decompress: invalid zstd frame");
        if (frame_size == ZSTD_CONTENTSIZE_UNKNOWN)
            throw std::runtime_error("ShardReader decompress: unknown zstd frame size");

        std::string out(static_cast<size_t>(frame_size), '\0');
        size_t written = ddict_
            ? ZSTD_decompress_usingDDict(dctx, out.data(), out.size(), src, src_size, ddict_)
            : ZSTD_decompressDCtx(dctx, out.data(), out.size(), src, src_size);
        if (ZSTD_isError(written))
            throw std::runtime_error(std::string("ShardReader decompress: ") +
                                     ZSTD_getErrorName(written));
        out.resize(written);

        auto* victim = &frame_cache.entries[0];
        for (auto& entry : frame_cache.entries) {
            if (entry.src == nullptr) {
                victim = &entry;
                break;
            }
            if (entry.age < victim->age)
                victim = &entry;
        }
        victim->src = src;
        victim->src_size = src_size;
        victim->ref_ptr = ref_ptr;
        victim->ref_size = ref_size;
        victim->use_prefix = use_prefix;
        victim->age = frame_cache.next_age++;
        victim->data = std::move(out);
        return victim->data;
    }

    const std::string& decompress_zstd_frame_cached(const uint8_t* src,
                                                    size_t src_size,
                                                    bool use_prefix) const
    {
        const char* ref_ptr = use_prefix && !ref_prefix_.empty() ? ref_prefix_.data() : nullptr;
        size_t ref_size = use_prefix ? ref_prefix_.size() : 0;
        return decompress_zstd_frame_cached_with_ref(src, src_size, ref_ptr, ref_size);
    }

    const std::string& decompress_v4_chunk(const GenomeDirEntry& e,
                                           const CheckpointEntry* cps,
                                           uint32_t chunk_idx) const {
        const uint8_t* blob = base_ + header_->blob_area_offset + e.blob_offset;
        uint32_t block_offset = cps[chunk_idx].block_offset;
        uint32_t next_offset = (chunk_idx + 1 < e.n_checkpoints)
            ? cps[chunk_idx + 1].block_offset
            : e.blob_len_cmp;
        return decompress_zstd_frame_cached(blob + block_offset,
                                            next_offset - block_offset,
                                            false);
    }

    const std::string& decompress_v4_delta_chunk(const GenomeDirEntry& e,
                                                 const ParsedSequenceBlob& parsed,
                                                 const CheckpointEntry* cps,
                                                 uint32_t chunk_idx) const {
        const uint8_t* blob = base_ + header_->blob_area_offset + e.blob_offset;
        uint32_t block_offset = cps[chunk_idx].block_offset;
        uint32_t next_offset = (chunk_idx + 1 < e.n_checkpoints)
            ? cps[chunk_idx + 1].block_offset
            : e.blob_len_cmp;

        if ((e.flags & GenomeMeta::FLAG_DELTA) == 0 || header_->codec != 3)
            return decompress_zstd_frame_cached(blob + block_offset,
                                                next_offset - block_offset,
                                                false);

        uint32_t ref_idx = delta_ref_index();
        const GenomeDirEntry& ref_e = dir_[ref_idx];
        const CheckpointEntry* ref_cps = checkpoints_for(ref_e);
        if (!ref_cps || ref_e.n_checkpoints == 0) {
            {
                std::unique_lock<std::mutex> lk(ref_prefix_mx_);
                if (ref_prefix_.empty())
                    ref_prefix_ = decompress_sequence_blob(ref_e);
            }
            return decompress_zstd_frame_cached(blob + block_offset,
                                                next_offset - block_offset,
                                                true);
        }

        uint64_t block_start = cps[chunk_idx].symbol_offset;
        uint64_t block_end = (chunk_idx + 1 < e.n_checkpoints)
            ? cps[chunk_idx + 1].symbol_offset
            : parsed.hdr.seq_len;
        uint64_t block_len = block_end - block_start;

        auto ref_begin = ref_cps;
        auto ref_end = ref_cps + ref_e.n_checkpoints;
        auto ref_it = std::upper_bound(ref_begin, ref_end, block_start,
                                       [](uint64_t value, const CheckpointEntry& cp) {
                                           return value < cp.symbol_offset;
                                       });
        uint32_t ref_chunk_idx = (ref_it == ref_begin)
            ? 0u
            : static_cast<uint32_t>((ref_it - ref_begin) - 1);
        const std::string& ref_chunk = decompress_v4_chunk(ref_e, ref_cps, ref_chunk_idx);
        uint64_t ref_chunk_start = ref_cps[ref_chunk_idx].symbol_offset;
        uint64_t ref_offset = block_start > ref_chunk_start ? block_start - ref_chunk_start : 0;
        size_t ref_size = ref_offset < ref_chunk.size()
            ? static_cast<size_t>(std::min<uint64_t>(block_len, ref_chunk.size() - ref_offset))
            : 0;
        const char* ref_ptr = ref_size > 0 ? ref_chunk.data() + ref_offset : nullptr;
        return decompress_zstd_frame_cached_with_ref(blob + block_offset,
                                                     next_offset - block_offset,
                                                     ref_ptr,
                                                     ref_size);
    }

    const std::string& decompress_v4_delta_chunk_full_ref(const GenomeDirEntry& e,
                                                          const ParsedSequenceBlob& parsed,
                                                          const CheckpointEntry* cps,
                                                          uint32_t chunk_idx) const {
        const uint8_t* blob = base_ + header_->blob_area_offset + e.blob_offset;
        uint32_t block_offset = cps[chunk_idx].block_offset;
        uint32_t next_offset = (chunk_idx + 1 < e.n_checkpoints)
            ? cps[chunk_idx + 1].block_offset
            : e.blob_len_cmp;

        if ((e.flags & GenomeMeta::FLAG_DELTA) == 0 || header_->codec != 3)
            return decompress_zstd_frame_cached(blob + block_offset,
                                                next_offset - block_offset,
                                                false);

        uint32_t ref_idx = delta_ref_index();
        const GenomeDirEntry& ref_e = dir_[ref_idx];
        {
            std::unique_lock<std::mutex> lk(ref_prefix_mx_);
            if (ref_prefix_.empty())
                ref_prefix_ = decompress_sequence_blob(ref_e);
        }

        uint64_t block_start = cps[chunk_idx].symbol_offset;
        uint64_t block_end = (chunk_idx + 1 < e.n_checkpoints)
            ? cps[chunk_idx + 1].symbol_offset
            : parsed.hdr.seq_len;
        uint64_t block_len = block_end - block_start;
        size_t ref_size = block_start < ref_prefix_.size()
            ? static_cast<size_t>(std::min<uint64_t>(block_len, ref_prefix_.size() - block_start))
            : 0;
        const char* ref_ptr = ref_size > 0 ? ref_prefix_.data() + block_start : nullptr;
        return decompress_zstd_frame_cached_with_ref(blob + block_offset,
                                                     next_offset - block_offset,
                                                     ref_ptr,
                                                     ref_size);
    }

    // Decompress sequence directly into dst without touching the LRU frame cache.
    // Uses symbol_offset deltas (known from v4 checkpoint metadata) to pre-allocate
    // the exact output size upfront, avoiding ZSTD_getFrameContentSize and per-chunk
    // resize calls. For sequential linear scans where frames are never revisited.
    void decompress_sequence_blob_into(const GenomeDirEntry& e, std::string& dst) const {
        const uint8_t* blob = base_ + header_->blob_area_offset + e.blob_offset;
        const size_t blob_len = e.blob_len_cmp;
        auto parsed = parse_sequence_blob(blob, blob_len);

        ZSTD_DCtx* dctx = get_thread_dctx();

        const CheckpointEntry* cps = checkpoints_for(e);
        if (!cps || e.n_checkpoints == 0) {
            // Single frame
            if (parsed.hdr.flags & SEQB_FLAG_PACKED_2BIT) {
                uint32_t packed_size = (parsed.hdr.seq_len + 3) / 4;
                dst.resize(packed_size);
                size_t written = ddict_
                    ? ZSTD_decompress_usingDDict(dctx, dst.data(), packed_size,
                                                 blob + parsed.hdr.blocks_offset,
                                                 blob_len - parsed.hdr.blocks_offset, ddict_)
                    : ZSTD_decompressDCtx(dctx, dst.data(), packed_size,
                                          blob + parsed.hdr.blocks_offset,
                                          blob_len - parsed.hdr.blocks_offset);
                if (ZSTD_isError(written))
                    throw std::runtime_error(std::string("ShardReader decompress: ") +
                                             ZSTD_getErrorName(written));
                uint32_t n_ambig = parsed.hdr.reserved;
                const AmbigEntry* ambig = n_ambig > 0
                    ? reinterpret_cast<const AmbigEntry*>(
                        blob + parsed.hdr.headers_offset + parsed.hdr.headers_len)
                    : nullptr;
                dst = unpack_sequence_2bit(
                    reinterpret_cast<const uint8_t*>(dst.data()), written,
                    ambig, n_ambig, parsed.hdr.seq_len);
                return;
            }
            dst.resize(parsed.hdr.seq_len);
            size_t written = ddict_
                ? ZSTD_decompress_usingDDict(dctx, dst.data(), parsed.hdr.seq_len,
                                             blob + parsed.hdr.blocks_offset,
                                             blob_len - parsed.hdr.blocks_offset, ddict_)
                : ZSTD_decompressDCtx(dctx, dst.data(), parsed.hdr.seq_len,
                                      blob + parsed.hdr.blocks_offset,
                                      blob_len - parsed.hdr.blocks_offset);
            if (ZSTD_isError(written))
                throw std::runtime_error(std::string("ShardReader decompress: ") +
                                         ZSTD_getErrorName(written));
            dst.resize(written);
            return;
        }

        // Multi-chunk: pre-allocate full output; write each chunk at its known offset.
        // Chunk i decompresses to [cps[i].symbol_offset, cps[i+1].symbol_offset).
        dst.resize(parsed.hdr.seq_len);
        size_t out_written = 0;
        for (uint32_t i = 0; i < e.n_checkpoints; ++i) {
            uint64_t chunk_start = cps[i].symbol_offset;
            uint64_t chunk_end   = (i + 1 < e.n_checkpoints)
                ? cps[i + 1].symbol_offset
                : static_cast<uint64_t>(parsed.hdr.seq_len);
            size_t chunk_out = static_cast<size_t>(chunk_end - chunk_start);

            uint32_t block_offset = cps[i].block_offset;
            uint32_t next_block   = (i + 1 < e.n_checkpoints)
                ? cps[i + 1].block_offset : static_cast<uint32_t>(blob_len);

            size_t written = ddict_
                ? ZSTD_decompress_usingDDict(dctx, dst.data() + chunk_start, chunk_out,
                                             blob + block_offset,
                                             next_block - block_offset, ddict_)
                : ZSTD_decompressDCtx(dctx, dst.data() + chunk_start, chunk_out,
                                      blob + block_offset,
                                      next_block - block_offset);
            if (ZSTD_isError(written))
                throw std::runtime_error(std::string("ShardReader decompress: ") +
                                         ZSTD_getErrorName(written));
            out_written += written;
        }
        dst.resize(out_written);
    }

    std::string decompress_sequence_blob(const GenomeDirEntry& e) const {
        const uint8_t* blob = base_ + header_->blob_area_offset + e.blob_offset;
        const size_t blob_len = e.blob_len_cmp;
        auto parsed = parse_sequence_blob(blob, blob_len);

        const CheckpointEntry* cps = checkpoints_for(e);
        if (!cps || e.n_checkpoints == 0) {
            if (parsed.hdr.flags & SEQB_FLAG_PACKED_2BIT) {
                // Bypass cache: decompress packed bytes and unpack to ASCII
                ZSTD_DCtx* dctx = get_thread_dctx();
                uint32_t packed_size = (parsed.hdr.seq_len + 3) / 4;
                std::string packed(packed_size, '\0');
                size_t written = ddict_
                    ? ZSTD_decompress_usingDDict(dctx, packed.data(), packed_size,
                                                 blob + parsed.hdr.blocks_offset,
                                                 blob_len - parsed.hdr.blocks_offset, ddict_)
                    : ZSTD_decompressDCtx(dctx, packed.data(), packed_size,
                                          blob + parsed.hdr.blocks_offset,
                                          blob_len - parsed.hdr.blocks_offset);
                if (ZSTD_isError(written))
                    throw std::runtime_error(std::string("ShardReader decompress 2bit: ") +
                                             ZSTD_getErrorName(written));
                uint32_t n_ambig = parsed.hdr.reserved;
                const AmbigEntry* ambig = n_ambig > 0
                    ? reinterpret_cast<const AmbigEntry*>(
                        blob + parsed.hdr.headers_offset + parsed.hdr.headers_len)
                    : nullptr;
                return unpack_sequence_2bit(
                    reinterpret_cast<const uint8_t*>(packed.data()), written,
                    ambig, n_ambig, parsed.hdr.seq_len);
            }
            const std::string& seq = decompress_zstd_frame_cached(
                blob + parsed.hdr.blocks_offset,
                blob_len - parsed.hdr.blocks_offset,
                (e.flags & GenomeMeta::FLAG_DELTA) != 0);
            return seq;
        }

        std::string out;
        out.reserve(parsed.hdr.seq_len);
        for (uint32_t i = 0; i < e.n_checkpoints; ++i) {
            const std::string& chunk = decompress_v4_delta_chunk_full_ref(e, parsed, cps, i);
            out.append(chunk);
        }
        // Checkpointed 2-bit blobs: each chunk decompressed to packed bytes;
        // unpack the accumulated packed bytes to ASCII now.
        if (parsed.hdr.flags & SEQB_FLAG_PACKED_2BIT) {
            uint32_t n_ambig = parsed.hdr.reserved;
            const AmbigEntry* ambig = n_ambig > 0
                ? reinterpret_cast<const AmbigEntry*>(
                    blob + parsed.hdr.headers_offset + parsed.hdr.headers_len)
                : nullptr;
            return unpack_sequence_2bit(
                reinterpret_cast<const uint8_t*>(out.data()), out.size(),
                ambig, n_ambig, parsed.hdr.seq_len);
        }
        return out;
    }

    // Hint to the OS that the compressed blob pages are no longer needed after
    // decompression. The decompressed FASTA is returned as a string; the kernel
    // can reclaim the mmap'd page-cache pages to reduce RSS under memory pressure.
    static void release_blob_pages(const uint8_t* src, size_t src_size) noexcept {
#ifdef MADV_DONTNEED
        static const size_t PAGE = 4096;
        uintptr_t addr    = reinterpret_cast<uintptr_t>(src);
        uintptr_t aligned = addr & ~(PAGE - 1);
        size_t    len     = (addr + src_size) - aligned;
        ::madvise(reinterpret_cast<void*>(aligned), len, MADV_DONTNEED);
#endif
    }

    std::string decompress_blob(const GenomeDirEntry& e) const {
        const uint8_t* src      = base_ + header_->blob_area_offset + e.blob_offset;
        const size_t   src_size = e.blob_len_cmp;

        // MEM-delta blob: decode DeltaBlob using cached anchor sequence.
        if ((e.flags & GenomeMeta::FLAG_DELTA) && header_->codec == 4)
            return decode_mem_delta(mem_anchor_seq_, src, src_size);

        if (header_->version >= 4) {
            auto parsed = parse_sequence_blob(src, src_size);
            std::string seq = decompress_sequence_blob(e);
            // Do NOT release blob pages here — caller handles MADV_DONTNEED
            // at the shard level after all genomes in a batch are fetched.
            return reconstruct_fasta_from_seq(seq,
                                              parsed.contig_ends,
                                              parsed.hdr.n_contigs,
                                              parsed.headers ? parsed.headers : "",
                                              parsed.hdr.headers_len);
        }

        const CheckpointEntry* cps = checkpoints_for(e);
        if (!cps || e.n_checkpoints == 0) {
            return decompress_zstd_frame_cached(src, src_size, (e.flags & GenomeMeta::FLAG_DELTA) != 0);
        }

        std::string out;
        out.reserve(e.blob_len_raw);
        for (uint32_t i = 0; i < e.n_checkpoints; ++i) {
            uint32_t block_offset = cps[i].block_offset;
            uint32_t next_offset = (i + 1 < e.n_checkpoints)
                ? cps[i + 1].block_offset
                : e.blob_len_cmp;
            const std::string& chunk = decompress_zstd_frame_cached(src + block_offset,
                                                                    next_offset - block_offset,
                                                                    (e.flags & GenomeMeta::FLAG_DELTA) != 0);
            out.append(chunk);
        }
        return out;
    }

    std::string decompress_sequence_slice(const GenomeDirEntry& e,
                                          uint64_t start,
                                          uint64_t length) const
    {
        if (length == 0) return {};
        uint64_t end = start + length;

        const uint8_t* src = base_ + header_->blob_area_offset + e.blob_offset;
        const size_t src_size = e.blob_len_cmp;

        if ((e.flags & GenomeMeta::FLAG_DELTA) && header_->codec == 4)
            return decode_mem_delta_slice(mem_anchor_seq_, src, src_size, start, length);

        if (header_->version >= 4) {
            auto parsed = parse_sequence_blob(src, src_size);
            if (start >= parsed.hdr.seq_len) return {};
            uint64_t clamped_end = std::min<uint64_t>(parsed.hdr.seq_len, end);

            const CheckpointEntry* cps = checkpoints_for(e);
            if (!cps || e.n_checkpoints == 0) {
                const std::string& seq = decompress_zstd_frame_cached(
                    src + parsed.hdr.blocks_offset,
                    src_size - parsed.hdr.blocks_offset,
                    (e.flags & GenomeMeta::FLAG_DELTA) != 0);
                return seq.substr(static_cast<size_t>(start),
                                  static_cast<size_t>(clamped_end - start));
            }

            auto begin = cps;
            auto end_it = cps + e.n_checkpoints;
            auto it = std::upper_bound(begin, end_it, start,
                                       [](uint64_t value, const CheckpointEntry& cp) {
                                           return value < cp.symbol_offset;
                                       });
            uint32_t chunk_idx = (it == begin)
                ? 0u
                : static_cast<uint32_t>((it - begin) - 1);

            std::string out;
            out.reserve(clamped_end - start);
            for (uint32_t i = chunk_idx; i < e.n_checkpoints && out.size() < clamped_end - start; ++i) {
                uint64_t block_start = cps[i].symbol_offset;
                if (block_start >= clamped_end)
                    break;
                const std::string& seq = decompress_v4_delta_chunk(e, parsed, cps, i);
                uint64_t block_end = block_start + seq.size();
                if (block_end <= start) continue;
                uint64_t copy_start = std::max(start, block_start);
                uint64_t copy_end = std::min(clamped_end, block_end);
                out.append(seq.data() + (copy_start - block_start), copy_end - copy_start);
            }
            return out;
        }

        const CheckpointEntry* cps = checkpoints_for(e);
        if (!cps || e.n_checkpoints == 0) {
            const std::string& fasta = decompress_zstd_frame_cached(src, src_size, (e.flags & GenomeMeta::FLAG_DELTA) != 0);
            return extract_sequence_slice_from_fasta_block(fasta, start, length);
        }

        auto begin = cps;
        auto end_it = cps + e.n_checkpoints;
        auto it = std::upper_bound(begin, end_it, start,
                                   [](uint64_t value, const CheckpointEntry& cp) {
                                       return value < cp.symbol_offset;
                                   });
        uint32_t chunk_idx = (it == begin)
            ? 0u
            : static_cast<uint32_t>((it - begin) - 1);

        std::string out;
        out.reserve(length);
        for (uint32_t i = chunk_idx; i < e.n_checkpoints && out.size() < length; ++i) {
            uint64_t block_start = cps[i].symbol_offset;
            if (block_start >= end)
                break;
            uint32_t block_offset = cps[i].block_offset;
            uint32_t next_offset = (i + 1 < e.n_checkpoints)
                ? cps[i + 1].block_offset
                : e.blob_len_cmp;
            uint64_t local_start = start > block_start ? start - block_start : 0;
            uint64_t want = length - out.size();
            const std::string& fasta = decompress_zstd_frame_cached(src + block_offset,
                                                                    next_offset - block_offset,
                                                                    (e.flags & GenomeMeta::FLAG_DELTA) != 0);
            if (!append_packed_sequence_slice(out, fasta, cps[i], local_start, want))
                out += extract_sequence_slice_from_fasta_block(fasta, local_start, want);
        }
        return out;
    }
};

ShardReader::ShardReader() = default;

ShardReader::~ShardReader() {
    if (impl_) impl_->reset();
}

ShardReader::ShardReader(ShardReader&&) noexcept = default;
ShardReader& ShardReader::operator=(ShardReader&&) noexcept = default;

void ShardReader::open(const uint8_t* base, uint64_t section_offset, uint64_t section_size) {
    if (!impl_) impl_ = std::make_unique<Impl>();
    else        impl_->reset();
    impl_->setup(base + section_offset, section_size);
}

void ShardReader::open_file(const std::filesystem::path& path) {
    if (!impl_) impl_ = std::make_unique<Impl>();
    else        impl_->reset();

    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) throw std::runtime_error("ShardReader: cannot open: " + path.string());
    size_t size = static_cast<size_t>(f.tellg());
    f.seekg(0);
    impl_->owned_data_.resize(size);
    f.read(reinterpret_cast<char*>(impl_->owned_data_.data()), size);
    if (!f) throw std::runtime_error("ShardReader: read error: " + path.string());

    impl_->setup(impl_->owned_data_.data(), size);
}

bool     ShardReader::is_open()   const { return impl_ && impl_->header_ != nullptr; }
uint32_t ShardReader::shard_id()  const { return impl_->header_->shard_id; }
uint32_t ShardReader::n_genomes() const { return impl_->header_->n_genomes; }

void ShardReader::release_pages() const noexcept {
#ifdef MADV_DONTNEED
    if (!impl_ || !impl_->base_) return;
    // Release the entire shard section (header + dir + blobs).
    // base_ points to section start; section_size_ covers everything.
    // We advise DONTNEED on the whole region so subsequent accesses to
    // OTHER shards are not penalised — the kernel re-reads this shard
    // from NFS only if it is accessed again.
    static const size_t PAGE = 4096;
    uintptr_t addr    = reinterpret_cast<uintptr_t>(impl_->base_);
    uintptr_t aligned = addr & ~(PAGE - 1);
    size_t    len     = impl_->section_size_
                        ? ((addr + impl_->section_size_) - aligned)
                        : 0;
    if (len > 0)
        ::madvise(reinterpret_cast<void*>(aligned), len, MADV_DONTNEED);
#endif
}

std::string ShardReader::fetch_genome(GenomeId id) const {
    const GenomeDirEntry* e = impl_->find_genome(id);
    if (!e) throw std::runtime_error("ShardReader: genome_id not found");
    if (e->flags & GenomeMeta::FLAG_DELETED) return {};
    return impl_->decompress_blob(*e);
}

std::string ShardReader::fetch_genome_at(uint32_t dir_index) const {
    const GenomeDirEntry* e = impl_->entry_at(dir_index);
    if (!e) throw std::runtime_error("ShardReader: dir_index out of range");
    if (e->flags & GenomeMeta::FLAG_DELETED) return {};
    return impl_->decompress_blob(*e);
}

void ShardReader::fetch_sequence_at_into(uint32_t dir_index,
                                         std::string& seq_buf,
                                         const uint32_t*& contig_ends,
                                         uint32_t& n_contigs) const
{
    contig_ends = nullptr;
    n_contigs   = 0;

    const GenomeDirEntry* e = impl_->entry_at(dir_index);
    if (!e) throw std::runtime_error("ShardReader: dir_index out of range");
    if (e->flags & GenomeMeta::FLAG_DELETED) { seq_buf.clear(); return; }

    // v4 non-MEM-delta: uncompressed metadata prefix lives in mmap, sequence in zstd blocks
    if (impl_->header_->version >= 4 && impl_->header_->codec != 4) {
        const uint8_t* blob = impl_->base_ + impl_->header_->blob_area_offset + e->blob_offset;
        auto parsed = parse_sequence_blob(blob, e->blob_len_cmp);
        impl_->decompress_sequence_blob_into(*e, seq_buf);  // uncached, writes directly into seq_buf
        contig_ends = parsed.contig_ends;
        n_contigs   = parsed.hdr.n_contigs;
        return;
    }

    // Fallback (v2/v3 or codec=4): decompress full FASTA and strip headers
    std::string fasta = impl_->decompress_blob(*e);
    seq_buf.clear();
    seq_buf.reserve(fasta.size());
    for (size_t i = 0; i < fasta.size(); ) {
        if (fasta[i] == '>') {
            while (i < fasta.size() && fasta[i] != '\n') ++i;
            if (i < fasta.size()) ++i;
            continue;
        }
        if (fasta[i] != '\n' && fasta[i] != '\r') seq_buf.push_back(fasta[i]);
        ++i;
    }
}

std::string ShardReader::fetch_sequence_slice(GenomeId id, uint64_t start, uint64_t length) const {
    const GenomeDirEntry* e = impl_->find_genome(id);
    if (!e) throw std::runtime_error("ShardReader: genome_id not found");
    if (e->flags & GenomeMeta::FLAG_DELETED) return {};
    return impl_->decompress_sequence_slice(*e, start, length);
}

std::string ShardReader::fetch_sequence_slice_at(uint32_t dir_index,
                                                 uint64_t start,
                                                 uint64_t length) const {
    const GenomeDirEntry* e = impl_->entry_at(dir_index);
    if (!e) throw std::runtime_error("ShardReader: dir_index out of range");
    if (e->flags & GenomeMeta::FLAG_DELETED) return {};
    return impl_->decompress_sequence_slice(*e, start, length);
}

const GenomeDirEntry* ShardReader::dir_begin() const { return impl_->dir_; }
const GenomeDirEntry* ShardReader::dir_end()   const {
    return impl_->dir_ + impl_->header_->n_genomes;
}

const GenomeDirEntry* ShardReader::dir_entry(uint32_t dir_index) const {
    return impl_->entry_at(dir_index);
}

} // namespace genopack
