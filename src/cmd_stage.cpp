#include "genopack/stage.hpp"
#include "genopack/util.hpp"
#include "genopack/types.hpp"

#include <spdlog/spdlog.h>

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <filesystem>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <vector>

namespace genopack {

namespace {

std::string record_taxonomy(const BuildRecord& r) {
    for (const auto& [k, v] : r.extra_fields)
        if (k == "taxonomy") return v;
    return {};
}

// Genus bucket for ordering only — approximate is fine (does not affect the final
// archive; within-shard order is recomputed at build). Keeps same-genus genomes
// contiguous in the cache so build's per-genus shard buckets stay dense.
std::string genus_key(const std::string& tax) {
    auto pos = tax.find("g__");
    if (pos == std::string::npos) return "~unclassified";
    size_t end = tax.find(';', pos);
    if (end == std::string::npos) end = tax.size();
    return tax.substr(pos, end - pos);
}

std::string read_one(const BuildRecord& rec) {
    int fd = ::open(rec.file_path.c_str(), O_RDONLY);
    if (fd < 0) throw std::runtime_error("open failed: " + rec.file_path.string());
#ifdef POSIX_FADV_SEQUENTIAL
    ::posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
    std::string seq;
    try {
        seq = decompress_gz_fd(fd, rec.file_path);
    } catch (...) {
        ::close(fd);
        throw;
    }
    ::close(fd);
    return seq;
}

} // namespace

int cmd_stage(const std::string& input_tsv, const std::string& output,
              int threads, int block_mb) {
    auto t0 = std::chrono::steady_clock::now();

    std::vector<BuildRecord> records = parse_tsv_records(input_tsv);
    const size_t N = records.size();
    spdlog::info("stage: {} genomes from {}", N, input_tsv);

    std::stable_sort(records.begin(), records.end(),
        [](const BuildRecord& a, const BuildRecord& b) {
            return genus_key(record_taxonomy(a)) < genus_key(record_taxonomy(b));
        });

    StageWriter w;
    w.open(output, static_cast<uint64_t>(std::max(1, block_mb)) << 20);

    const int n_threads = std::max(1, threads);

    struct Item { std::string accession, taxonomy, sequence; };
    std::queue<Item>        q;
    std::mutex              mx;
    std::condition_variable cv_full, cv_empty;
    const size_t            qmax = static_cast<size_t>(n_threads) * 4 + 8;

    std::atomic<size_t> next{0};
    std::atomic<size_t> failed{0};
    std::atomic<bool>   readers_done{false};

    auto reader = [&]() {
        for (;;) {
            size_t i = next.fetch_add(1, std::memory_order_relaxed);
            if (i >= N) break;
            const BuildRecord& rec = records[i];
            std::string seq;
            try {
                seq = read_one(rec);
            } catch (const std::exception& ex) {
                spdlog::warn("stage: skip {}: {}", rec.accession, ex.what());
                failed.fetch_add(1, std::memory_order_relaxed);
                continue;
            }
            Item it{rec.accession, record_taxonomy(rec), std::move(seq)};
            std::unique_lock lk(mx);
            cv_full.wait(lk, [&] { return q.size() < qmax; });
            q.push(std::move(it));
            lk.unlock();
            cv_empty.notify_one();
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(n_threads);
    for (int t = 0; t < n_threads; ++t) pool.emplace_back(reader);

    std::thread joiner([&] {
        for (auto& th : pool) th.join();
        { std::lock_guard lk(mx); readers_done.store(true); }
        cv_empty.notify_all();
    });

    size_t written = 0;
    for (;;) {
        std::unique_lock lk(mx);
        cv_empty.wait(lk, [&] { return !q.empty() || readers_done.load(); });
        if (q.empty()) {
            if (readers_done.load()) break;
            continue;
        }
        Item it = std::move(q.front());
        q.pop();
        lk.unlock();
        cv_full.notify_one();
        w.add(it.accession, it.taxonomy, it.sequence);
        if (++written % 50000 == 0) spdlog::info("stage: {} / {} written", written, N);
    }
    joiner.join();

    w.finalize();

    double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    uint64_t fsz = 0;
    std::error_code ec;
    fsz = std::filesystem::file_size(output, ec);
    spdlog::info("stage: wrote {} genomes ({} failed) to {} — {:.1f} MB, {:.1f}s",
                 written, failed.load(), output,
                 static_cast<double>(fsz) / 1e6, dt);
    return 0;
}

} // namespace genopack
