#pragma once
#include "shard.hpp"
#include <filesystem>

namespace genopack {

struct RepackConfig {
    char              taxonomy_rank        = 'g';              // 'g'=genus, 'f'=family
    ShardWriterConfig shard_cfg;
    size_t            threads              = 1;                // OMP threads for parallel decompression
    uint64_t          max_bucket_bytes     = 32ULL << 30;      // flush all buckets when total exceeds this (32 GB)
    bool              verbose              = false;
};

// Read `input_gpk`, re-shard genomes by taxonomy (genus by default),
// write a new taxonomy-grouped archive to `output_gpk`.
//
// Genome IDs are preserved. ACCX/TAXN/KMRX/CIDX/HNSW/TXDB sections are
// rebuilt or copied verbatim (they are indexed by genome_id, not shard).
// Only CATL and GIDX sections are rebuilt (they embed shard_id).
//
// Effect: genomes from the same genus cluster into contiguous shards,
// so visit_shard_batches() reads only the shards for the requested taxon
// instead of the entire archive.
void repack_archive(const std::filesystem::path& input_gpk,
                    const std::filesystem::path& output_gpk,
                    const RepackConfig& cfg = {});

} // namespace genopack
