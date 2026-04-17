#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace genopack {

// ── MEM-delta genomic compression ────────────────────────────────────────────
//
// Encodes a query FASTA relative to a reference panel using Maximal Exact
// Matches (MEMs) seeded with k=31 k-mers (step=16). MEMs are found via a
// flat open-addressing hash table (3-5× faster than unordered_map), then the
// full DeltaBlob (MEM list + verbatim residue) is zstd-compressed as a single
// frame — sorted integer arrays in the MEM list compress 5-8×.
//
// DeltaBlob v2 wire format (little-endian):
//   [MemDeltaHeader]
//   [uint32_t contig_ends[n_contigs]]
//   [headers bytes ('\0'-separated)]
//   [MemDeltaChunkDesc[n_chunks]]
//   [zstd-compressed chunk payloads]
//
// Each chunk payload stores only its local MEM triples and verbatim residue, so
// readers can decode overlapping chunks for coordinate slices without expanding
// the full genome delta blob.

struct MemEntry {
    uint32_t qpos;
    uint32_t rpos;
    uint32_t length;
};

struct FastaComponents {
    std::string              seq;         // concatenated pure sequence (ACGT/N, uppercase)
    std::vector<uint32_t>    contig_ends; // exclusive end positions in seq per contig
    std::vector<std::string> headers;    // header text (no '>')
};

#pragma pack(push, 1)
struct MemDeltaHeader {
    uint32_t magic;
    uint16_t version;
    uint16_t flags;
    uint32_t anchor_seq_len;
    uint32_t query_seq_len;
    uint32_t n_contigs;
    uint32_t n_chunks;
    uint32_t headers_len;
    uint32_t reserved;
    uint64_t contig_ends_offset;
    uint64_t headers_offset;
    uint64_t chunks_offset;
};
static_assert(sizeof(MemDeltaHeader) == 56);

struct MemDeltaChunkDesc {
    uint32_t query_start;
    uint32_t query_len;
    uint32_t payload_offset;
    uint32_t payload_size;
};
static_assert(sizeof(MemDeltaChunkDesc) == 16);
#pragma pack(pop)

// ── Flat open-addressing k-mer index ─────────────────────────────────────────
// Power-of-2 capacity, linear probing, ~50% load factor.
// 3-5× faster than unordered_map for the lookup-heavy MEM finding phase
// because the key/value arrays are dense and cache-prefetch-friendly.

struct AnchorIndex {
    std::vector<uint64_t> keys; // UINT64_MAX = empty slot
    std::vector<uint32_t> vals;
    uint64_t              mask; // capacity - 1

    AnchorIndex() : mask(0) {}
    explicit AnchorIndex(size_t n_entries);

    void     insert(uint64_t key, uint32_t val);
    uint32_t find(uint64_t key) const; // returns UINT32_MAX if not found
};

// Strip headers and newlines from a FASTA string, returning components.
FastaComponents extract_fasta_components(const char* fasta, size_t len);

// Build a step-16 k=31 anchor k-mer index over the reference sequence.
AnchorIndex build_anchor_index(const std::string& seq);

// Encode query relative to anchor. Returns chunked MEM-delta bytes.
std::vector<uint8_t> encode_mem_delta(const FastaComponents& anchor,
                                       const AnchorIndex&     anchor_idx,
                                       const FastaComponents& query,
                                       size_t                 chunk_bases = 65536);

// Decode a DeltaBlob back to the original FASTA string.
// anchor_seq: concatenated pure sequence of the reference panel.
std::string decode_mem_delta(const std::string& anchor_seq,
                              const uint8_t*     blob,
                              size_t             blob_len);

// Decode only the requested pure-sequence slice [start, start + length).
std::string decode_mem_delta_slice(const std::string& anchor_seq,
                                   const uint8_t*     blob,
                                   size_t             blob_len,
                                   uint64_t           start,
                                   uint64_t           length);

} // namespace genopack
