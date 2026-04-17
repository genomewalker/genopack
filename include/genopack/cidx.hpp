#pragma once
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <cstdint>
#include <functional>
#include <string_view>
#include <vector>

namespace genopack {

// ── CIDX — Contig Accession Index ─────────────────────────────────────────────
//
// Maps contig accession strings → genome_id using FNV-1a-64 sorted array.
// One CIDX section is written per incremental build batch; multiple sections
// coexist in the archive and are searched newest→oldest by MergedCidxReader.
// genopack merge consolidates them into a single sorted section.
//
// Section layout:
//   CidxHeader (32 bytes)
//   CidxEntry[] sorted ascending by acc_hash (16 bytes each)
//
// Batch lookup (MergedCidxReader::batch_find):
//   1. Hash all N accessions in parallel across T threads           O(N/T)
//   2. Sort queries by hash (pdq/introsort)                         O(N log N)
//   3. Merge-join sorted queries against each sorted CIDX section   O(N + M) per section
//      with hardware prefetch; only unfound queries enter next pass

// 32-byte section header
struct CidxHeader {
    uint32_t magic;        // SEC_CIDX
    uint32_t version;      // 1
    uint64_t n_entries;    // number of CidxEntry records
    uint64_t batch_id;     // monotonically increasing; newest section = highest batch_id
    uint8_t  reserved[8];
};
static_assert(sizeof(CidxHeader) == 32);

// 16-byte entry, sorted by acc_hash ascending
struct CidxEntry {
    uint64_t acc_hash;   // FNV-1a-64(contig_accession_string)
    uint32_t genome_id;  // owning genome's genome_id
    uint32_t _pad;       // alignment
};
static_assert(sizeof(CidxEntry) == 16);

// Internal query record used by batch_find
struct CidxQuery {
    uint64_t hash;
    uint32_t orig_idx;   // position in the caller's input/output arrays
    uint32_t _pad;
};
static_assert(sizeof(CidxQuery) == 16);

// FNV-1a-64 hash of a contig accession string
uint64_t cidx_hash(std::string_view acc);

// Parse FASTA headers and call cb(contig_acc_token) for each '>' line.
// acc_token is the first whitespace-delimited token after '>'.
void parse_fasta_contig_accessions(std::string_view fasta,
                                   const std::function<void(std::string_view)>& cb);

// ── CidxWriter ────────────────────────────────────────────────────────────────

class CidxWriter {
public:
    void add(std::string_view contig_acc, uint32_t genome_id);
    void add_hash(uint64_t acc_hash, uint32_t genome_id);

    // Sort by acc_hash, write header + entries, return filled SectionDesc.
    SectionDesc finalize(AppendWriter& writer, uint64_t section_id, uint64_t batch_id = 0);

    size_t size() const { return entries_.size(); }

private:
    std::vector<CidxEntry> entries_;
};

// ── CidxReader ────────────────────────────────────────────────────────────────

class CidxReader {
public:
    void open(const uint8_t* base, uint64_t offset, uint64_t size);

    // Single lookup — binary search. Returns UINT32_MAX if not found.
    uint32_t find(uint64_t acc_hash) const;
    uint32_t find(std::string_view contig_acc) const;

    uint64_t         batch_id()  const;
    size_t           n_entries() const { return n_; }
    const CidxEntry* data()      const { return entries_; }

    void scan(const std::function<void(uint64_t acc_hash, uint32_t genome_id)>& cb) const;

private:
    const CidxHeader* hdr_     = nullptr;
    const CidxEntry*  entries_ = nullptr;
    uint64_t          n_       = 0;
};

// ── MergedCidxReader ──────────────────────────────────────────────────────────
// Searches multiple CIDX sections newest→oldest (by batch_id).

class MergedCidxReader {
public:
    void add_section(const uint8_t* base, uint64_t offset, uint64_t size);

    // Single lookup. Returns UINT32_MAX if not found in any section.
    uint32_t find(std::string_view contig_acc) const;
    uint32_t find_hash(uint64_t acc_hash) const;

    // ── High-performance batch lookup ─────────────────────────────────────────
    //
    // Maps N contig accession strings → genome_ids in-place.
    // out_genome_ids[i] = genome_id for accs[i], or UINT32_MAX if not found.
    //
    // Algorithm:
    //   Phase 1: hash all accs in parallel (n_threads workers)
    //   Phase 2: sort queries by hash
    //   Phase 3: merge-join sorted queries against each CIDX section
    //            (newest first; unfound queries recycled into next pass)
    //
    // n_threads = 0 → use hardware_concurrency()
    void batch_find(const std::string_view* accs,
                    uint32_t*               out_genome_ids,
                    size_t                  n,
                    size_t                  n_threads = 0) const;

    bool   empty()         const { return readers_.empty(); }
    size_t total_entries() const;

    // Merge-sort all sections into one sorted vector (used by genopack merge).
    std::vector<CidxEntry> merge_all() const;

private:
    std::vector<CidxReader> readers_; // sorted newest→oldest by batch_id

    // Merge-join one sorted CIDX section against sorted queries.
    // Fills out_genome_ids for matches; appends unfound queries to unfound_out.
    static void merge_join_section(const CidxReader&     reader,
                                   const CidxQuery*      queries,
                                   size_t                n_queries,
                                   uint32_t*             out_genome_ids,
                                   std::vector<CidxQuery>& unfound_out);
};

} // namespace genopack
