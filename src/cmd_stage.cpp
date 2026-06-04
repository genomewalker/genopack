#include <genopack/stage.hpp>
#include <genopack/util.hpp>
#include <spdlog/spdlog.h>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <filesystem>
#include <map>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include <dirent.h>
#include <fcntl.h>
#include <unistd.h>

namespace genopack {

int cmd_stage(const std::filesystem::path& input_tsv,
              const std::filesystem::path& output_path,
              int n_threads,
              uint64_t block_mb)
{
    namespace fs = std::filesystem;

    // Parse input TSV
    auto records = parse_tsv_records(input_tsv);
    if (records.empty()) throw std::runtime_error("stage: no records in TSV");
    spdlog::info("stage: {} genomes from {}", records.size(), input_tsv.string());

    // READDIRPLUS-style pre-warm: sweep the common parent once to warm NFS
    // dentry/attr cache. Skip if directory appears too large (>500k entries)
    // since readdir over millions of entries is slower than the files themselves.
    {
        std::string common_dir;
        for (const auto& r : records) {
            auto d = r.file_path.parent_path().string();
            if (common_dir.empty()) { common_dir = d; continue; }
            if (d != common_dir) { common_dir.clear(); break; }
        }
        if (!common_dir.empty()) {
            if (DIR* dp = ::opendir(common_dir.c_str())) {
                struct dirent* de;
                size_t n = 0;
                while ((de = ::readdir(dp)) && n < 500000) ++n;
                ::closedir(dp);
                if (n < 500000)
                    spdlog::info("stage: pre-warmed NFS cache for {} ({} entries)", common_dir, n);
                else
                    spdlog::info("stage: skipping pre-warm — directory too large (>500k entries)");
            }
        }
    }

    // Bounded work queue: index into records → output queue of StageRecords
    struct Task { size_t idx; };
    struct Result { size_t idx; StageRecord rec; bool ok; std::string err; };

    const size_t q_max = static_cast<size_t>(n_threads) * 4;
    std::queue<Task>   task_q;
    std::queue<Result> done_q;
    std::mutex         task_mx, done_mx;
    std::condition_variable task_cv, done_cv;
    bool               tasks_done = false;

    std::atomic<size_t> n_failed{0};

    // Producer: push all tasks
    std::thread producer([&]{
        for (size_t i = 0; i < records.size(); ++i) {
            std::unique_lock lk(task_mx);
            task_cv.wait(lk, [&]{ return task_q.size() < q_max; });
            task_q.push({i});
            lk.unlock();
            task_cv.notify_one();
        }
        // poison pills
        {
            std::unique_lock lk(task_mx);
            task_cv.wait(lk, [&]{ return task_q.size() < q_max; });
            tasks_done = true;
        }
        task_cv.notify_all();
    });

    // Worker threads: open FASTA, decompress, uppercase
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(n_threads));
    for (int wi = 0; wi < n_threads; ++wi) {
        workers.emplace_back([&]{
            while (true) {
                Task t;
                {
                    std::unique_lock lk(task_mx);
                    task_cv.wait(lk, [&]{ return !task_q.empty() || tasks_done; });
                    if (task_q.empty()) return;
                    t = task_q.front(); task_q.pop();
                }
                task_cv.notify_one();

                Result res;
                res.idx = t.idx;
                res.ok  = false;
                try {
                    const auto& rec = records[t.idx];
                    int fd = ::open(rec.file_path.c_str(), O_RDONLY);
                    std::string seq = (fd >= 0)
                        ? decompress_gz_fd(fd, rec.file_path)  // fd closed inside
                        : decompress_gz(rec.file_path);

                    // Get taxonomy from extra_fields
                    std::string tax;
                    for (const auto& [k, v] : rec.extra_fields)
                        if (k == "taxonomy") { tax = v; break; }

                    res.rec = StageRecord{rec.accession, std::move(tax), std::move(seq)};
                    res.ok  = true;
                } catch (const std::exception& ex) {
                    res.err = ex.what();
                    ++n_failed;
                }

                {
                    std::unique_lock lk(done_mx);
                    done_cv.wait(lk, [&]{ return done_q.size() < q_max; });
                    done_q.push(std::move(res));
                }
                done_cv.notify_one();
            }
        });
    }

    // Main thread: drain done_q → StageWriter in input order via sliding reorder
    // window. Max out-of-orderness ≤ q_max (bounded by concurrent task count),
    // so the map never grows beyond ~4×n_threads entries. O(n_threads) memory.
    StageWriter writer;
    writer.open(output_path, block_mb << 20);

    auto t0 = std::chrono::steady_clock::now();

    std::map<size_t, Result> reorder;
    size_t next_out  = 0;
    size_t n_received = 0;

    while (n_received < records.size()) {
        // Drain available results into reorder map
        {
            std::unique_lock lk(done_mx);
            done_cv.wait(lk, [&]{ return !done_q.empty(); });
            while (!done_q.empty()) {
                auto r = std::move(done_q.front()); done_q.pop();
                reorder.emplace(r.idx, std::move(r));
                ++n_received;
            }
        }
        done_cv.notify_all();

        // Emit all consecutive results starting from next_out
        while (true) {
            auto it = reorder.find(next_out);
            if (it == reorder.end()) break;
            auto& res = it->second;
            if (!res.ok)
                spdlog::warn("stage: skipping {}: {}", records[res.idx].accession, res.err);
            else
                writer.add(std::move(res.rec));
            reorder.erase(it);
            ++next_out;
        }

        if (n_received % 10000 == 0)
            spdlog::info("stage: {}/{} read ({} failed)", n_received, records.size(),
                         n_failed.load());
    }

    producer.join();
    for (auto& w : workers) w.join();
    writer.finalize();

    auto elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    const auto fsize = fs::file_size(output_path);
    spdlog::info("stage: {} genomes → {} ({:.1f} GB) in {:.0f}s ({:.0f} genomes/s)",
                 writer.n_genomes(),
                 output_path.string(),
                 fsize / 1e9,
                 elapsed,
                 writer.n_genomes() / elapsed);

    if (n_failed > 0)
        spdlog::warn("stage: {} genomes failed (see warnings above)", n_failed.load());
    return 0;
}

} // namespace genopack
