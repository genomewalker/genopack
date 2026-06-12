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
    // Route by genus_hash: all aamers for this genome go to the single partition
    // for its genus.  finalize() can emit each genus as soon as its partition is
    // processed — no global accumulator.  One lock-acquire per genome.
    const uint32_t pi = static_cast<uint32_t>(genus_hash >> (64 - PART_BITS));
    Part& p = parts_[pi];
    std::lock_guard<std::mutex> lk(p.mtx);
    for (size_t i = 0; i < n; ++i)
        p.buf.push_back({genus_hash, sorted_aamers[i]});
    if (static_cast<int>(p.buf.size()) >= buf_cap_) flush_part(p);
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
    // Requires quiescence: no concurrent add_genome() or add_member() calls.
    //
    // Because add_genome() routes by genus_hash, every tuple for genus G is in
    // exactly one partition.  We can load each partition, sort, scan, and emit
    // its genera immediately — no global accumulator, peak RAM = largest partition.
    const int nt = std::max(1, threads);

    // Flush any in-memory buffers that never hit buf_cap_.
    for (uint32_t pi = 0; pi < N_PARTS; ++pi) {
        Part& p = parts_[pi];
        if (!p.buf.empty()) flush_part(p);
    }

    const auto n_dirty = static_cast<uint32_t>(dirty_parts_.size());
    std::exception_ptr caught;
    std::mutex caught_mtx;

    #pragma omp parallel for schedule(dynamic,1) num_threads(nt)
    for (uint32_t di = 0; di < n_dirty; ++di) {
        try {
            {
                std::lock_guard<std::mutex> elk(caught_mtx);
                if (caught) continue;
            }
            const uint32_t pi = dirty_parts_[di];
            Part& p = parts_[pi];

            if (p.fp) {
                if (std::fclose(p.fp) != 0)
                    throw std::runtime_error("pcore_spiller: fclose failed: " + p.path);
                p.fp = nullptr;
            }

            FILE* f = std::fopen(p.path.c_str(), "rb");
            if (!f) throw std::runtime_error("pcore_spiller: cannot open " + p.path);
            std::fseek(f, 0, SEEK_END);
            const long fsz = std::ftell(f);
            if (fsz < 0) { std::fclose(f); throw std::runtime_error("pcore_spiller: ftell: " + p.path); }
            std::fseek(f, 0, SEEK_SET);
            const size_t n = static_cast<size_t>(fsz) / sizeof(Tuple);
            std::vector<Tuple> tuples(n);
            if (n > 0 && std::fread(tuples.data(), sizeof(Tuple), n, f) != n) {
                std::fclose(f);
                throw std::runtime_error("pcore_spiller: read error: " + p.path);
            }
            std::fclose(f);
            std::remove(p.path.c_str());

            if (tuples.empty()) continue;
            std::sort(tuples.begin(), tuples.end(), tuple_lt);

            std::unordered_map<uint64_t, PartialGenus> per_genus;
            scan_partition(tuples.data(), tuples.data() + tuples.size(), per_genus);
            { std::vector<Tuple>().swap(tuples); }

            for (auto& [gh, pg] : per_genus) {
                uint32_t n_mem = 0;
                auto it = member_count_.find(gh); // read-only under quiescence
                if (it != member_count_.end()) n_mem = it->second;
                cb(gh, pg.aamers, pg.counts, n_mem);
            }
        } catch (...) {
            std::lock_guard<std::mutex> elk(caught_mtx);
            if (!caught) caught = std::current_exception();
        }
    }

    if (caught) {
        cleanup();
        std::rethrow_exception(caught);
    }
}

void PcoreSpillWriter::cleanup() {
    for (uint32_t pi : dirty_parts_) {
        Part& p = parts_[pi];
        if (p.fp) { std::fclose(p.fp); p.fp = nullptr; }
        if (p.file_exists) std::remove(p.path.c_str());
    }
}

} // namespace genopack
