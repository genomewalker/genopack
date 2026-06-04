#pragma once
#include <cstdint>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

namespace genopack {

struct StageRecord {
    std::string accession;
    std::string taxonomy;
    std::string sequence;   // gzip-decoded, uppercase, N/IUPAC preserved
};

// ── Writer ────────────────────────────────────────────────────────────────────

class StageWriter {
public:
    // block_bytes: target uncompressed size per zstd block (~64MB default)
    void open(const std::filesystem::path& path, uint64_t block_bytes = 64u << 20);
    void add(StageRecord r);
    void finalize();
    size_t n_genomes() const { return n_genomes_; }

private:
    void flush_block();

    struct BlockDesc { uint64_t offset; uint64_t csize; uint64_t n_genomes; };

    int                   fd_          = -1;
    uint64_t              block_bytes_ = 64u << 20;
    uint64_t              write_offset_= 0;
    size_t                n_genomes_   = 0;
    std::vector<char>     pending_;        // uncompressed accumulation buffer
    size_t                pending_n_      = 0;
    std::vector<BlockDesc>blocks_;
    uint64_t              data_start_     = 0; // byte offset where blocks begin
};

// ── Reader ────────────────────────────────────────────────────────────────────

class StageReader {
public:
    void open(const std::filesystem::path& path);
    uint64_t n_genomes() const { return n_genomes_; }
    uint64_t n_blocks()  const { return n_blocks_; }

    void scan(const std::function<void(StageRecord&&)>& cb) const {
        scan_blocks(0, n_blocks_, cb);
    }
    void scan_blocks(uint64_t block_start, uint64_t block_end,
                     const std::function<void(StageRecord&&)>& cb) const;

private:
    struct BlockDesc { uint64_t offset; uint64_t csize; uint64_t n_genomes; };
    std::filesystem::path path_;
    uint64_t              n_genomes_ = 0;
    uint64_t              n_blocks_  = 0;
    uint64_t              data_start_= 0;
    std::vector<BlockDesc>blocks_;
};

// ── CLI entry point ───────────────────────────────────────────────────────────
int cmd_stage(const std::filesystem::path& input_tsv,
              const std::filesystem::path& output_path,
              int n_threads = 48,
              uint64_t block_mb = 64);

} // namespace genopack
