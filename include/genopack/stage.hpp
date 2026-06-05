#pragma once
#include <cstdint>
#include <filesystem>
#include <functional>
#include <memory>
#include <string>

namespace genopack {

// ── .gstage: durable local sequence cache ─────────────────────────────────────
// A flat, block-compressed store of gzip-decoded genome FASTA bytes, produced once
// by a high-concurrency NFS reader (`genopack stage`) and replayed by sequential
// local reads for fast iterative rebuilds (`genopack build --from-stage`).
//
// Lossless: stores the decoded FASTA verbatim (no case folding), recompressed with
// zstd. A rebuild from a .gstage yields a byte-identical archive to the NFS build.
//
// On-disk layout (see stage.cpp):
//   [Header 48 B]
//   [Data block 0] .. [Data block n-1]   zstd, ~block_bytes uncompressed each
//   [Block table: n_blocks x 24 B]
//   [Meta block]                         zstd: (accession, taxonomy) directory
//
// Byte order is native (this is a local-node cache, not a portable archive).

constexpr uint64_t kStageMagic = 0x4753544147450001ULL; // "GSTAGE\0\1"

struct StageRecord {
    std::string accession;
    std::string taxonomy;   // populated by scan_meta(); empty in scan()
    std::string sequence;   // decoded FASTA bytes, verbatim
};

class StageWriter {
public:
    StageWriter();
    ~StageWriter();
    StageWriter(const StageWriter&) = delete;
    StageWriter& operator=(const StageWriter&) = delete;

    void open(const std::filesystem::path& path, uint64_t block_bytes = 64ull << 20);
    // accession + decoded FASTA bytes (verbatim); taxonomy may be empty.
    void add(const std::string& accession, const std::string& taxonomy,
             const std::string& sequence);
    void finalize();
    uint64_t n_genomes() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class StageReader {
public:
    StageReader();
    ~StageReader();
    StageReader(const StageReader&) = delete;
    StageReader& operator=(const StageReader&) = delete;

    void open(const std::filesystem::path& path);
    uint64_t n_genomes() const;
    uint64_t n_blocks() const;

    // (accession, taxonomy) directory — cheap, no sequence decompression.
    void scan_meta(const std::function<void(const std::string& acc,
                                            const std::string& tax)>& cb) const;
    // Full records in storage order (sequential block decompression). taxonomy empty.
    void scan(const std::function<void(StageRecord&&)>& cb) const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace genopack
