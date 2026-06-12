#pragma once
// ── GenusAamerCounter — hash-partitioned cache-resident exact prevalence counter ─
// 2048 partitions by top-11 aamer-hash bits; each partition is a flat
// open-addressing table {uint64_t keys[], uint32_t counts[]} at ≤70% load.
// 12 bytes/slot × (1/0.7) ≈ 17 bytes/entry — less than ankerl::unordered_dense
// (22 bytes/entry) for this write-heavy, small-partition workload where DRAM
// bandwidth per entry dominates over load-factor differences.
// Counts are EXACT — CQF/BQF ruled out; singleton and ceil(θ·n) boundaries must
// be exact for CORE correctness.
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace genopack {

class GenusAamerCounter {
public:
    static constexpr int      PBITS = 11;
    static constexpr uint32_t NP    = 1u << PBITS;

    GenusAamerCounter() : parts_(NP) {}

    void add_sorted(const uint64_t* a, size_t n) {
        for (size_t i = 0; i < n; ++i)
            parts_[static_cast<uint32_t>(a[i] >> (64 - PBITS))].add(a[i]);
    }
    void add_sorted(const std::vector<uint64_t>& a) { add_sorted(a.data(), a.size()); }

    // Cache-resident PARALLEL batch fold: outer loop over partitions keeps each
    // partition's table L2-resident while every member's slice streams through it.
    void add_chunk(const std::vector<const std::vector<uint64_t>*>& members) {
        const size_t M = members.size();
        if (M == 0) return;
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

    void finalize_sorted(std::vector<uint64_t>& aamers, std::vector<uint32_t>& counts) const {
        // Prefix-sum occupancies → per-partition write offsets.
        // Writes directly into pre-allocated output; avoids the per[NP] Pair-buffer
        // that previously peaked at nd×16 B simultaneously with the output nd×12 B.
        size_t nd = 0;
        std::vector<uint32_t> off(NP + 1);
        for (uint32_t p = 0; p < NP; ++p) { off[p] = static_cast<uint32_t>(nd); nd += parts_[p].occ; }
        off[NP] = static_cast<uint32_t>(nd);

        aamers.resize(nd); counts.resize(nd);

        #pragma omp parallel for schedule(dynamic, 8)
        for (int p = 0; p < static_cast<int>(NP); ++p) {
            const auto& part = parts_[static_cast<uint32_t>(p)];
            if (part.occ == 0) continue;
            const uint32_t base = off[static_cast<uint32_t>(p)];
            const uint32_t n    = part.occ;
            uint32_t idx = base;
            for (uint32_t s = 0; s < part.size; ++s)
                if (part.counts[s]) { aamers[idx] = part.keys[s]; counts[idx++] = part.counts[s]; }
            // Sort this partition's slice by key. Thread-local scratch avoids
            // simultaneous heap allocation across all 2048 partitions; each thread
            // reuses its buffer across iterations (capacity grows to the per-thread max).
            // Partitions are disjoint by top-11 bits so the global output is already
            // sorted at partition granularity — only intra-partition sort is needed.
            thread_local std::vector<std::pair<uint64_t, uint32_t>> scratch;
            scratch.resize(n);
            for (uint32_t i = 0; i < n; ++i) scratch[i] = {aamers[base + i], counts[base + i]};
            std::sort(scratch.begin(), scratch.end(),
                [](const auto& a, const auto& b) { return a.first < b.first; });
            for (uint32_t i = 0; i < n; ++i) { aamers[base + i] = scratch[i].first; counts[base + i] = scratch[i].second; }
        }
    }

private:
    struct Part {
        std::vector<uint64_t> keys;
        std::vector<uint16_t> counts;     // 0 = empty slot; present key ≥ 1; safe because
                                          // per-genome aamers are sort+unique'd before insertion,
                                          // so max count = genus member count (≤65535, guarded at build start)
        uint32_t size = 0, mask = 0, occ = 0;
        void reset(uint32_t cap) { size=cap; mask=cap-1; keys.assign(cap,0); counts.assign(cap,0); occ=0; }
        void grow() {
            std::vector<uint64_t> ok = std::move(keys);
            std::vector<uint16_t> oc = std::move(counts);
            const uint32_t osz = size;
            reset(size << 1);
            for (uint32_t i = 0; i < osz; ++i) if (oc[i]) put_(ok[i], oc[i]);
        }
        inline void put_(uint64_t h, uint16_t c) {
            uint32_t s = static_cast<uint32_t>(h) & mask;
            while (counts[s]) { if (keys[s]==h) { counts[s]+=c; return; } s=(s+1)&mask; }
            keys[s]=h; counts[s]=c; ++occ;
        }
        inline void add(uint64_t h) {
            if (size == 0) reset(256);
            if ((occ+1)*10 >= size*7) grow();
            uint32_t s = static_cast<uint32_t>(h) & mask;
            while (counts[s]) { if (keys[s]==h) { ++counts[s]; return; } s=(s+1)&mask; }
            keys[s]=h; counts[s]=1; ++occ;
        }
    };
    std::vector<Part> parts_;
};

} // namespace genopack
