#include <genopack/pcore_spiller.hpp>
#include <algorithm>
#include <stdexcept>
#include <sys/mman.h>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace genopack {

PcoreSpillWriter::PcoreSpillWriter(std::filesystem::path spill_dir, int buf_tuples)
    : spill_dir_(std::move(spill_dir)), buf_cap_(buf_tuples),
      parts_(new Part[N_PARTS])
{
    std::filesystem::create_directories(spill_dir_);
    // Pin the Part headers (229 KB) so they survive the high-memory shard-write
    // phase without being paged out.  finalize() would otherwise take 1-2 s on
    // page-faults when first touching each Part struct.  Ignore failure gracefully.
    ::mlock(parts_.get(), N_PARTS * sizeof(Part));
    for (uint32_t i = 0; i < N_PARTS; ++i) {
        parts_[i].path = (spill_dir_ / ("pspill_" + std::to_string(i) + ".bin")).string();
        parts_[i].buf.reserve(static_cast<size_t>(buf_tuples));
    }
}

PcoreSpillWriter::~PcoreSpillWriter() { cleanup(); }

void PcoreSpillWriter::flush_part(Part& p) {
    if (p.buf.empty()) return;
    if (!p.fp) {
        p.fp = std::fopen(p.path.c_str(), "wb");
        if (!p.fp) throw std::runtime_error("pcore_spiller: cannot open " + p.path);
        p.file_exists = true;
        const uint32_t pi = static_cast<uint32_t>(&p - parts_.get());
        std::lock_guard<std::mutex> dlk(dirty_mtx_);
        dirty_parts_.push_back(pi);
    }
    if (std::fwrite(p.buf.data(), sizeof(Tuple), p.buf.size(), p.fp) != p.buf.size())
        throw std::runtime_error("pcore_spiller: write error on " + p.path);
    p.buf.clear();
}

void PcoreSpillWriter::add_genome(uint64_t genus_hash,
                                   const uint64_t* sorted_aamers, size_t n) {
    if (n == 0) return;
    size_t i = 0;
    while (i < n) {
        const uint32_t pi = static_cast<uint32_t>(sorted_aamers[i] >> (64 - PART_BITS));
        Part& p = parts_[pi];
        std::lock_guard<std::mutex> lk(p.mtx);
        // Push all aamers that belong to this partition (contiguous run since
        // the input is sorted and partition index is the top PART_BITS bits).
        while (i < n && static_cast<uint32_t>(sorted_aamers[i] >> (64 - PART_BITS)) == pi) {
            p.buf.push_back({genus_hash, sorted_aamers[i]});
            ++i;
        }
        if (static_cast<int>(p.buf.size()) >= buf_cap_) flush_part(p);
    }
    total_added_.fetch_add(static_cast<uint64_t>(n), std::memory_order_relaxed);
}

void PcoreSpillWriter::add_member(uint64_t genus_hash) {
    std::lock_guard<std::mutex> lk(n_mem_mtx_);
    ++member_count_[genus_hash];
}

struct PartialGenus {
    std::vector<uint64_t> aamers;
    std::vector<uint32_t> counts;
};

static bool tuple_lt(const PcoreSpillWriter::Tuple& a, const PcoreSpillWriter::Tuple& b) {
    if (a.genus_hash != b.genus_hash) return a.genus_hash < b.genus_hash;
    return a.aamer_hash < b.aamer_hash;
}

static void scan_partition(const PcoreSpillWriter::Tuple* begin,
                           const PcoreSpillWriter::Tuple* end,
                           std::unordered_map<uint64_t, PartialGenus>& per_genus) {
    using T = PcoreSpillWriter::Tuple;
    const T* it = begin;
    while (it != end) {
        const uint64_t gh = it->genus_hash;
        PartialGenus& pg = per_genus[gh];
        while (it != end && it->genus_hash == gh) {
            const uint64_t ah = it->aamer_hash;
            uint32_t cnt = 0;
            while (it != end && it->genus_hash == gh && it->aamer_hash == ah) {
                ++cnt;
                ++it;
            }
            pg.aamers.push_back(ah);
            pg.counts.push_back(cnt);
        }
    }
}

void PcoreSpillWriter::finalize(int threads, Callback cb) {
    // Streaming batch design: process BATCH_SIZE partitions at a time so that
    // phase1 RAM = O(BATCH × n_genera × aamers_per_partition) not
    // O(N_PARTS × n_genera × aamers_per_partition).
    //
    // Partitions are processed in ascending pi order within each batch so that
    // appending aamer runs from batch 0, 1, … produces a globally-sorted
    // aamer list per genus (disjoint aamer_hash ranges across partitions).
    const uint32_t BATCH = static_cast<uint32_t>(std::max(1, threads) * 8);

    // Running accumulator: genus_hash → {aamers, counts} built incrementally.
    struct AccumGenus { std::vector<uint64_t> aamers; std::vector<uint32_t> counts; };
    std::unordered_map<uint64_t, AccumGenus> accum;

    for (uint32_t batch_start = 0; batch_start < N_PARTS; batch_start += BATCH) {
        const uint32_t batch_end = std::min(batch_start + BATCH, N_PARTS);
        const uint32_t batch_sz  = batch_end - batch_start;
        std::vector<std::unordered_map<uint64_t, PartialGenus>> batch_phase1(batch_sz);

        #pragma omp parallel for schedule(dynamic,1) num_threads(threads)
        for (uint32_t bi = 0; bi < batch_sz; ++bi) {
            const uint32_t pi = batch_start + bi;
            Part& p = parts_[pi];
            std::vector<Tuple> tuples;

            // Close write handle in parallel so GPFS sync overlaps across threads.
            if (p.fp) { std::fclose(p.fp); p.fp = nullptr; }

            FILE* f = p.file_exists ? std::fopen(p.path.c_str(), "rb") : nullptr;
            if (f) {
                std::fseek(f, 0, SEEK_END);
                const size_t n = static_cast<size_t>(std::ftell(f)) / sizeof(Tuple);
                std::fseek(f, 0, SEEK_SET);
                tuples.resize(n);
                if (n > 0 && std::fread(tuples.data(), sizeof(Tuple), n, f) != n)
                    throw std::runtime_error("pcore_spiller: read error on " + p.path);
                std::fclose(f);
                std::remove(p.path.c_str());
            }

            tuples.insert(tuples.end(), p.buf.begin(), p.buf.end());
            p.buf.clear();

            if (tuples.empty()) continue;
            std::sort(tuples.begin(), tuples.end(), tuple_lt);
            scan_partition(tuples.data(), tuples.data() + tuples.size(), batch_phase1[bi]);
        }

        // Serial merge: append this batch's aamer runs into the accumulator.
        // bi is iterated in ascending order so aamer_hash ranges concatenate
        // in the correct globally-sorted order.
        for (uint32_t bi = 0; bi < batch_sz; ++bi) {
            for (auto& [gh, pg] : batch_phase1[bi]) {
                AccumGenus& ag = accum[gh];
                ag.aamers.insert(ag.aamers.end(), pg.aamers.begin(), pg.aamers.end());
                ag.counts.insert(ag.counts.end(), pg.counts.begin(), pg.counts.end());
            }
            // Free this partition's partial results immediately.
            batch_phase1[bi].clear();
        }
    }

    // Emit one callback per genus.  n_mem from the member counter.
    for (auto& [gh, ag] : accum) {
        uint32_t n_mem = 0;
        {
            std::lock_guard<std::mutex> lk(n_mem_mtx_);
            auto it = member_count_.find(gh);
            if (it != member_count_.end()) n_mem = it->second;
        }
        cb(gh, ag.aamers, ag.counts, n_mem);
    }
}

void PcoreSpillWriter::cleanup() {
    for (uint32_t pi : dirty_parts_)
        std::remove(parts_[pi].path.c_str());
}

} // namespace genopack
