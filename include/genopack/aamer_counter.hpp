#pragma once
// ── GenusAamerCounter — hash-partitioned cache-resident exact prevalence counter ─
// Replaces std::unordered_map<uint64,uint32> for per-genus aamer member-counting.
// The map is DRAM-latency bound (~1e7 ++ops/s: every operator[] is a ~100ns cache
// miss + pointer chase + periodic rehash), which makes keep-all counting (~6.5M
// aamers/genome × members ≈ 1e13 ops/part) infeasible.
//
// Insight (room room-63dcac93): the input is already sorted-unique per member and
// the genus aamer union SATURATES, so partition by the hash HIGH bits into P flat
// open-addressing tables — each holds ~union/P distinct keys and stays L2-resident,
// turning ++count into a ~5ns hit (~20-50× single-threaded; partitions are
// independent so it also parallelizes lock-free across worker threads). Output is
// globally hash-sorted by concatenating the per-partition (high-bit-ordered) runs.
// Counts are EXACT (no fingerprint approximation) — singletons and the ceil(θ·n)
// CORE boundary must be exact, which rules out CQF/BQF.
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace genopack {

class GenusAamerCounter {
public:
    static constexpr int      PBITS = 11;            // 2048 partitions by top bits
    static constexpr uint32_t NP    = 1u << PBITS;

    GenusAamerCounter() : parts_(NP) {}

    // Fold one member's sorted-unique aamers (each counted once for this member).
    void add_sorted(const uint64_t* a, size_t n) {
        for (size_t i = 0; i < n; ++i)
            parts_[static_cast<uint32_t>(a[i] >> (64 - PBITS))].add(a[i]);
    }
    void add_sorted(const std::vector<uint64_t>& a) { add_sorted(a.data(), a.size()); }

    // Cache-resident PARALLEL batch fold of a chunk's member arrays. Processes
    // PARTITION-by-partition (outer) — each partition's flat table stays hot in L2
    // while every member's slice for that partition streams through it — and the
    // partitions are independent so OpenMP parallelises across them with zero locks
    // (thread t writes only the partitions it owns). This is the ~1000× path; doing
    // it per-genome (touching all 2048 tables per genome) defeats cache residency.
    void add_chunk(const std::vector<const std::vector<uint64_t>*>& members) {
        const size_t M = members.size();
        if (M == 0) return;
        // Per-member partition boundaries: off[m][p] = first index of partition p in
        // members[m] (sorted on hash, so each partition is a contiguous slice).
        std::vector<std::vector<uint32_t>> off(M);
        for (size_t m = 0; m < M; ++m) {
            const auto& a = *members[m];
            auto& o = off[m]; o.resize(NP + 1);
            size_t i = 0;
            for (uint32_t p = 0; p < NP; ++p) {
                o[p] = static_cast<uint32_t>(i);
                while (i < a.size() && static_cast<uint32_t>(a[i] >> (64 - PBITS)) == p) ++i;
            }
            o[NP] = static_cast<uint32_t>(a.size());
        }
        #pragma omp parallel for schedule(dynamic, 8)
        for (int p = 0; p < static_cast<int>(NP); ++p) {
            Part& part = parts_[static_cast<uint32_t>(p)];
            for (size_t m = 0; m < M; ++m) {
                const uint64_t* a = members[m]->data();
                const uint32_t e = off[m][p + 1];
                for (uint32_t i = off[m][p]; i < e; ++i) part.add(a[i]);
            }
        }
    }

    size_t distinct() const { size_t s = 0; for (const auto& p : parts_) s += p.occ; return s; }
    bool   empty()    const { for (const auto& p : parts_) if (p.occ) return false; return true; }

    // Emit globally-sorted (aamer, count). Partition index == hash-high-bits, so
    // concatenating per-partition-sorted runs is globally sorted ascending.
    void finalize_sorted(std::vector<uint64_t>& aamers, std::vector<uint32_t>& counts) const {
        aamers.clear(); counts.clear();
        // Sort each partition independently (parallel); concatenate in partition order
        // (== hash-high-bit order) for a globally-sorted result.
        std::vector<std::vector<std::pair<uint64_t, uint32_t>>> per(NP);
        #pragma omp parallel for schedule(dynamic, 8)
        for (int p = 0; p < static_cast<int>(NP); ++p) {
            const auto& part = parts_[static_cast<uint32_t>(p)];
            if (part.occ == 0) continue;
            auto& v = per[static_cast<uint32_t>(p)]; v.reserve(part.occ);
            for (uint32_t s = 0; s < part.size; ++s)
                if (part.counts[s]) v.emplace_back(part.keys[s], part.counts[s]);
            std::sort(v.begin(), v.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
        }
        size_t nd = 0; for (const auto& v : per) nd += v.size();
        aamers.reserve(nd); counts.reserve(nd);
        for (const auto& v : per) for (const auto& [k, c] : v) { aamers.push_back(k); counts.push_back(c); }
    }

private:
    struct Part {
        std::vector<uint64_t> keys;
        std::vector<uint32_t> counts;     // 0 = empty slot; present key has count ≥ 1
        uint32_t size = 0, mask = 0, occ = 0;
        void reset(uint32_t cap) { size = cap; mask = cap - 1; keys.assign(cap, 0); counts.assign(cap, 0); occ = 0; }
        void grow() {
            std::vector<uint64_t> ok = std::move(keys);
            std::vector<uint32_t> oc = std::move(counts);
            const uint32_t osz = size;
            reset(size << 1);
            for (uint32_t i = 0; i < osz; ++i) if (oc[i]) put_(ok[i], oc[i]);
        }
        inline void put_(uint64_t h, uint32_t c) {       // raw insert during rehash
            uint32_t s = static_cast<uint32_t>(h) & mask;
            while (counts[s]) { if (keys[s] == h) { counts[s] += c; return; } s = (s + 1) & mask; }
            keys[s] = h; counts[s] = c; ++occ;
        }
        inline void add(uint64_t h) {
            if (size == 0) reset(256);
            if ((occ + 1) * 10 >= size * 7) grow();      // keep load < 0.7
            uint32_t s = static_cast<uint32_t>(h) & mask;
            while (counts[s]) { if (keys[s] == h) { ++counts[s]; return; } s = (s + 1) & mask; }
            keys[s] = h; counts[s] = 1; ++occ;
        }
    };
    std::vector<Part> parts_;
};

} // namespace genopack
