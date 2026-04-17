#pragma once
#include "archive.hpp"
#include "types.hpp"
#include "catalog.hpp"
#include <functional>
#include <cstddef>

namespace genopack {

// ── ScanEngine ───────────────────────────────────────────────────────────────
// High-throughput streaming scan over archive genomes.
//
// Physical layout:
//   - 1-2 I/O threads prefetch chunk-sized slabs from disk
//   - N-1 worker threads decompress + process from RAM
//   - Shard-granular scheduling: thread i owns shard stripe i
//   - Reusable slab buffers: no per-genome heap allocation
//
// NFS: 1-2 sequential I/O streams, each reading contiguous shard stripes.
// NVMe: all workers issue reads directly (high queue depth, no I/O threads needed).
//
// Usage:
//   ScanEngine engine;
//   engine.scan_all(archive, [](GenomeId id, std::string_view fasta,
//                               const GenomeMeta& meta) {
//       // process genome -- called from worker thread, may be concurrent
//   });

class ScanEngine {
public:
    struct Config {
        size_t io_threads      = 1;   // dedicated I/O threads (0 = worker threads do I/O)
        size_t worker_threads  = 15;  // genome-processing threads
        size_t slab_size_mb    = 64;  // per-slab buffer size (should match chunk size)
        size_t prefetch_depth  = 3;   // chunks prefetched ahead of current
        bool   nfs_mode        = true;// if true: sequential I/O, coarse slabs
                                      // if false: NVMe mode, fine-grained parallel reads
    };

    using GenomeCallback = std::function<void(GenomeId id,
                                              std::string_view fasta,
                                              const GenomeMeta& meta)>;

    ScanEngine() = default;
    explicit ScanEngine(Config cfg) : cfg_(cfg) {}

    // Scan all live genomes. Callback is called from worker threads (may be concurrent).
    // Order is chunk-sequential within a stripe; no global ordering guarantee.
    void scan_all(const ArchiveReader& archive, GenomeCallback cb);

    // Scan only genomes matching query. Applies CATL + taxonomy pre-filter;
    // decompresses only survivor chunks.
    void scan_filtered(const ArchiveReader& archive,
                       const ExtractQuery& q,
                       GenomeCallback cb);

private:
    Config cfg_;
    // Implementation in scan_engine.cpp (Phase 4)
};

} // namespace genopack
