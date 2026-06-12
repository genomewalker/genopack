#pragma once
// PcoreSpillWriter — external-sort PCORE accumulator for large-scale builds.
//
// Key insight: each genome's aamer array is sorted-unique, so the count of
// (genus_hash, aamer_hash) tuple occurrences across all genomes directly equals
// aamer prevalence.  genome_id is NOT stored in the spill — n_mem is tracked via
// a separate per-genus member counter (add_member, one call per genome).
//
// Partitions by aamer_hash top-11 bits so tuples for one genus spread evenly
// across all 2048 partitions regardless of genus size.  A 600k-member genus ×
// 100 aamers contributes ~29k tuples per partition (not 60M into one partition).
//
// finalize() two-phase algorithm:
//   Phase 1 (parallel): for each partition, load spill file (or in-memory buffer
//     if small enough to avoid file I/O), sort by (genus_hash, aamer_hash), scan,
//     collect {aamers[], counts[]} per genus.
//   Phase 2 (serial): for each genus, concatenate sorted aamer runs from all
//     partitions (disjoint hash ranges → already globally sorted), look up n_mem
//     from member_count_, call cb(genus_hash, aamers, counts, n_mem).
//
// Memory at finalize = sum of unique (genus, aamer) pairs × 12 bytes + phase1
// overhead.  genome_ids are never stored; n_mem is O(n_genera).
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace genopack {

class PcoreSpillWriter {
public:
    static constexpr uint32_t N_PARTS   = 2048;
    static constexpr int      PART_BITS = 11;
    static constexpr int      BUF_TUPLES = 65536;

    struct Tuple {
        uint64_t genus_hash;
        uint64_t aamer_hash;
    };

    using Callback = std::function<void(uint64_t genus_hash,
                                        std::vector<uint64_t>& aamers,
                                        std::vector<uint32_t>& counts,
                                        uint32_t n_mem)>;

    explicit PcoreSpillWriter(std::filesystem::path spill_dir,
                              int buf_tuples = BUF_TUPLES);
    ~PcoreSpillWriter();

    PcoreSpillWriter(const PcoreSpillWriter&)            = delete;
    PcoreSpillWriter& operator=(const PcoreSpillWriter&) = delete;

    // Thread-safe: add all aamers for one genome in a single batch.
    // sorted_aamers must be sorted ascending (standard aamer order); since the
    // partition index is the top PART_BITS of each aamer_hash, a sorted aamer
    // array is also partition-ordered, so each partition's mutex is acquired
    // exactly once per genome instead of once per aamer.
    void add_genome(uint64_t genus_hash, const uint64_t* sorted_aamers, size_t n);
    void add_genome(uint64_t genus_hash, const std::vector<uint64_t>& sorted_aamers) {
        add_genome(genus_hash, sorted_aamers.data(), sorted_aamers.size());
    }

    // Thread-safe: record one genome contributing to this genus.  Call once per
    // genome (not per aamer).
    void add_member(uint64_t genus_hash);

    // Flush → sort partitions → assemble per-genus (aamers, counts, n_mem) → cb.
    void finalize(int threads, Callback cb);

    bool empty() const { return total_added_.load(std::memory_order_relaxed) == 0; }

private:
    struct Part {
        std::mutex         mtx;
        std::vector<Tuple> buf;
        FILE*              fp          = nullptr;
        bool               file_exists = false;
        std::string        path;
    };

    std::filesystem::path          spill_dir_;
    int                            buf_cap_;
    std::unique_ptr<Part[]>        parts_;
    std::atomic<uint64_t>          total_added_{0};
    std::mutex                     n_mem_mtx_;
    std::unordered_map<uint64_t, uint32_t> member_count_;
    // Only partitions that overflowed their in-memory buffer (and opened a spill
    // file) appear here.  Avoids touching all 2048 Part structs in finalize().
    std::mutex                     dirty_mtx_;
    std::vector<uint32_t>          dirty_parts_;

    void flush_part(Part& p);
    void cleanup();
};

} // namespace genopack
