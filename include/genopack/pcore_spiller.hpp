#pragma once
// PcoreSpillWriter — external-sort PCORE accumulator for large-scale builds.
//
// Key insight: each genome's aamer array is sorted-unique, so the count of
// (genus_hash, aamer_hash) tuple occurrences across all genomes directly equals
// aamer prevalence.  genome_id is NOT stored in the spill — n_mem is tracked via
// a separate per-genus member counter (add_member, one call per genome).
//
// Partitions by genus_hash top-11 bits: all aamers for genus G land in exactly
// one partition (pi = genus_hash >> 53).  finalize() processes each partition
// independently and emits its genera immediately — peak RAM = largest single
// partition, not the sum of all unique (genus, aamer) pairs.
//
// finalize() algorithm (requires quiescence — no concurrent add_genome/add_member):
//   Flush remaining in-memory buffers.
//   For each dirty partition (parallel, schedule dynamic):
//     Close write handle, load spill file, sort by (genus_hash, aamer_hash),
//     scan into local per_genus map, emit each genus via cb immediately.
//   Peak memory ≈ O(largest partition size) instead of O(all unique pairs).
//
// Memory at finalize ≈ max(partition_file_size) × 2 + O(n_genera_in_partition).
// genome_ids are never stored; n_mem is O(n_genera).
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
    // All aamers route to the single partition for genus_hash — one lock-acquire
    // per genome regardless of aamer count.
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
