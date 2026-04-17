#include <genopack/scan_engine.hpp>
#include <genopack/shard.hpp>
#include <genopack/catalog.hpp>
#include <genopack/toc.hpp>
#include <genopack/format.hpp>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

namespace genopack {

// ── Work item: one shard section to decompress ──────────────────────────────

struct ShardWork {
    uint64_t offset;
    uint64_t size;
    uint32_t shard_id;
    // For filtered scan: genome_ids to extract (empty = all)
    std::vector<GenomeId> survivors;
};

// ── Bounded work queue ──────────────────────────────────────────────────────

class ShardQueue {
public:
    explicit ShardQueue(size_t capacity) : capacity_(capacity) {}

    void push(ShardWork item) {
        std::unique_lock<std::mutex> lk(mu_);
        not_full_.wait(lk, [&] { return q_.size() < capacity_ || done_; });
        if (done_) return;
        q_.push(std::move(item));
        not_empty_.notify_one();
    }

    bool pop(ShardWork& out) {
        std::unique_lock<std::mutex> lk(mu_);
        not_empty_.wait(lk, [&] { return !q_.empty() || done_; });
        if (q_.empty()) return false;
        out = std::move(q_.front());
        q_.pop();
        not_full_.notify_one();
        return true;
    }

    void finish() {
        std::lock_guard<std::mutex> lk(mu_);
        done_ = true;
        not_empty_.notify_all();
        not_full_.notify_all();
    }

private:
    std::mutex mu_;
    std::condition_variable not_empty_, not_full_;
    std::queue<ShardWork> q_;
    size_t capacity_;
    bool done_ = false;
};

// ── Internal helpers ────────────────────────────────────────────────────────

static void process_shard_all(const uint8_t* mmap_base,
                              const ShardWork& work,
                              const ArchiveReader& archive,
                              const ScanEngine::GenomeCallback& cb)
{
    ShardReader reader;
    reader.open(mmap_base, work.offset, work.size);

    uint32_t dir_idx = 0;
    for (auto* de = reader.dir_begin(); de != reader.dir_end(); ++de, ++dir_idx) {
        if (de->flags & GenomeMeta::FLAG_DELETED) continue;

        GenomeId gid = de->genome_id;
        std::string fasta = reader.fetch_genome_at(dir_idx);

        auto meta_opt = archive.genome_meta(gid);
        if (!meta_opt) continue;

        cb(gid, fasta, *meta_opt);
    }
}

static void process_shard_filtered(const uint8_t* mmap_base,
                                   const ShardWork& work,
                                   const ArchiveReader& archive,
                                   const ScanEngine::GenomeCallback& cb)
{
    ShardReader reader;
    reader.open(mmap_base, work.offset, work.size);

    std::unordered_set<GenomeId> survivor_set(work.survivors.begin(),
                                               work.survivors.end());

    uint32_t dir_idx = 0;
    for (auto* de = reader.dir_begin(); de != reader.dir_end(); ++de, ++dir_idx) {
        if (de->flags & GenomeMeta::FLAG_DELETED) continue;
        if (survivor_set.find(de->genome_id) == survivor_set.end()) continue;

        GenomeId gid = de->genome_id;
        std::string fasta = reader.fetch_genome_at(dir_idx);

        auto meta_opt = archive.genome_meta(gid);
        if (!meta_opt) continue;

        cb(gid, fasta, *meta_opt);
    }
}

// ── ScanEngine::scan_all ────────────────────────────────────────────────────

void ScanEngine::scan_all(const ArchiveReader& archive, GenomeCallback cb) {
    // Collect shard layout from archive
    struct ShardInfo {
        uint64_t offset;
        uint64_t size;
        uint32_t shard_id;
    };
    const uint8_t* mmap_base = nullptr;
    std::vector<ShardInfo> shards;
    archive.scan_shards([&](const uint8_t* data, uint64_t offset,
                            uint64_t size, uint32_t shard_id) {
        if (!mmap_base) mmap_base = data;
        shards.push_back({offset, size, shard_id});
    });

    if (shards.empty()) {
        spdlog::warn("ScanEngine::scan_all: no v2 shards found");
        return;
    }

    spdlog::info("ScanEngine: scanning {} shards, {} I/O + {} worker threads",
                 shards.size(), cfg_.io_threads, cfg_.worker_threads);

    int archive_fd = archive.fd();
    ShardQueue queue(cfg_.prefetch_depth);

    std::thread io_thread([&] {
        size_t lookahead = cfg_.prefetch_depth;
        for (size_t i = 0; i < shards.size(); ++i) {
            for (size_t j = i; j < std::min(i + lookahead, shards.size()); ++j) {
                if (cfg_.nfs_mode && archive_fd >= 0) {
                    posix_fadvise(archive_fd, shards[j].offset, shards[j].size,
                                  POSIX_FADV_WILLNEED);
                }
            }

            ShardWork work;
            work.offset   = shards[i].offset;
            work.size     = shards[i].size;
            work.shard_id = shards[i].shard_id;
            queue.push(std::move(work));
        }
        queue.finish();
    });

    // Worker threads: pop shard, decompress, call callback
    size_t n_workers = std::max(cfg_.worker_threads, size_t{1});
    std::vector<std::thread> workers;
    workers.reserve(n_workers);

    for (size_t i = 0; i < n_workers; ++i) {
        workers.emplace_back([&] {
            ShardWork work;
            while (queue.pop(work)) {
                process_shard_all(mmap_base, work, archive, cb);
            }
        });
    }

    io_thread.join();
    for (auto& w : workers) w.join();
}

// ── ScanEngine::scan_filtered ───────────────────────────────────────────────

void ScanEngine::scan_filtered(const ArchiveReader& archive,
                               const ExtractQuery& q,
                               GenomeCallback cb)
{
    // Pre-filter: get surviving genomes grouped by shard
    auto metas = archive.filter_meta(q);

    std::unordered_map<uint32_t, std::vector<GenomeId>> by_shard;
    for (const auto& m : metas) {
        if (!m.is_deleted())
            by_shard[m.shard_id].push_back(m.genome_id);
    }

    // Collect shard layout, only for shards with survivors
    struct ShardInfo {
        uint64_t offset;
        uint64_t size;
        uint32_t shard_id;
    };
    const uint8_t* mmap_base = nullptr;
    std::vector<ShardInfo> shards;
    archive.scan_shards([&](const uint8_t* data, uint64_t offset,
                            uint64_t size, uint32_t shard_id) {
        if (!mmap_base) mmap_base = data;
        if (by_shard.count(shard_id))
            shards.push_back({offset, size, shard_id});
    });

    if (shards.empty()) return;

    spdlog::info("ScanEngine: filtered scan: {} shards ({} genomes), {} I/O + {} workers",
                 shards.size(), metas.size(), cfg_.io_threads, cfg_.worker_threads);

    int archive_fd = archive.fd();
    ShardQueue queue(cfg_.prefetch_depth);

    std::thread io_thread([&] {
        size_t lookahead = cfg_.prefetch_depth;
        for (size_t i = 0; i < shards.size(); ++i) {
            for (size_t j = i; j < std::min(i + lookahead, shards.size()); ++j) {
                if (cfg_.nfs_mode && archive_fd >= 0) {
                    posix_fadvise(archive_fd, shards[j].offset, shards[j].size,
                                  POSIX_FADV_WILLNEED);
                }
            }

            ShardWork work;
            work.offset    = shards[i].offset;
            work.size      = shards[i].size;
            work.shard_id  = shards[i].shard_id;
            work.survivors = std::move(by_shard[shards[i].shard_id]);
            queue.push(std::move(work));
        }
        queue.finish();
    });

    size_t n_workers = std::max(cfg_.worker_threads, size_t{1});
    std::vector<std::thread> workers;
    workers.reserve(n_workers);

    for (size_t i = 0; i < n_workers; ++i) {
        workers.emplace_back([&] {
            ShardWork work;
            while (queue.pop(work)) {
                process_shard_filtered(mmap_base, work, archive, cb);
            }
        });
    }

    io_thread.join();
    for (auto& w : workers) w.join();
}

} // namespace genopack
