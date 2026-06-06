// Unit test for SEC_PCORE v1 — 3-run stratified PFOR + quantized prevalence.
//   1. PcoreWriter stratifies the dense union into singleton/multi/core runs with
//      TRUE-count-derived core membership; PcoreView::materialize round-trips the
//      union + dequantized prevalence; core_into yields the exact conserved core.
//   2. pcore_core_coverage reproduces the conserved-core (>=ceil(theta*n)) coverage.
//   3. REGRESSION: a >255-member genus keeps correct prevalence (~1.0, not 255/n)
//      and a non-empty core (the old u8-count-cap produced 0.255 + NaN coverage).

#include <genopack/pcore.hpp>
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "test_utils.hpp"

using namespace genopack;
using namespace genopack_test;

// prevalence of a specific aamer from a materialized (aamers,prevf) pair, or -1.
static float prev_of(const std::vector<uint64_t>& a, const std::vector<float>& p, uint64_t h) {
    for (size_t i = 0; i < a.size(); ++i) if (a[i] == h) return p[i];
    return -1.0f;
}

int main() {
    TempDir tmp("pcore");

    // ── Case 1: small genus — 3 members, counts 1->3, 2->3, 3->2, 4->1 ──────────
    // theta=0.90 → C=ceil(2.7)=3 → core={1,2}, multi={3}, singleton={4}.
    {
        std::vector<std::vector<uint64_t>> members = {{1,2,3,4}, {1,2,3}, {1,2}};
        const auto path = tmp.path / "small.bin";
        SectionDesc sd;
        { AppendWriter aw; aw.create(path);
          PcoreWriter pw(8, 8, 0.90f);
          pw.add_from_members(42, members);
          sd = pw.finalize(aw, 5, 123); aw.flush(); aw.close(); }
        require(sd.type == SEC_PCORE && sd.version == 2, "SEC_PCORE v1 (version 2)");

        MmapFileReader mm; mm.open(path);
        PcoreReader r; r.open(mm.data(), sd.file_offset, sd.compressed_size);
        require(r.n_entries() == 1 && r.codec() == PCORE_CODEC_V1, "one entry, v1 codec");
        require(std::abs(r.theta()-0.90f) < 1e-6f && r.frac_max_hash()==123, "theta/frac round-trip");

        PcoreView v = r.lookup(42);
        require(v.valid() && v.total() == 4 && v.n_members == 3, "union of 4 aamers, 3 members");
        require(v.n_singleton==1 && v.n_multi==1 && v.n_core==2, "runs 1/1/2 (single/multi/core)");
        require(!r.lookup(99).valid(), "absent key → invalid view");

        std::vector<uint64_t> a; std::vector<float> p; v.materialize(a, p);
        require(a.size()==4, "materialized union size 4");
        require(std::abs(prev_of(a,p,1)-1.0f)<1e-3f && std::abs(prev_of(a,p,2)-1.0f)<1e-3f, "core prev=1.0");
        require(prev_of(a,p,3) > 0.5f && prev_of(a,p,3) < 0.8f, "multi prev ~0.67");
        require(std::abs(prev_of(a,p,4) - 1.0f/3.0f) < 1e-3f, "singleton prev=1/3");

        require(pcore_core_threshold(3, 0.90f) == 3, "ceil(0.9*3)=3");
        std::vector<uint64_t> c; v.core_into(c);
        require(c.size()==2 && c[0]==1 && c[1]==2, "exact core run {1,2} sorted");
        require(std::abs(pcore_core_coverage(v, {1,2,9}, 0.90f) - 1.0f) < 1e-6f, "coverage 2/2 = 1.0");
        require(std::abs(pcore_core_coverage(v, {1,9},   0.90f) - 0.5f) < 1e-6f, "coverage 1/2 = 0.5");
    }

    // ── Case 2: REGRESSION — 1000-member genus (>255) ───────────────────────────
    // aamer 100 in all 1000 members (count=1000=n), 200 in 500, plus per-member
    // singletons. theta=0.90 → C=900. OLD u8-count cap: prev(100)=255/1000=0.255 and
    // C=900>255 → core empty → coverage NaN. v1 must give prev~1.0 and core={100}.
    {
        const uint32_t N = 1000;
        std::vector<std::vector<uint64_t>> members(N);
        for (uint32_t i = 0; i < N; ++i) {
            members[i].push_back(100);                 // in all
            if (i < 500) members[i].push_back(200);    // in half
            members[i].push_back(1000 + i);            // unique singleton
            std::sort(members[i].begin(), members[i].end());
        }
        const auto path = tmp.path / "big.bin";
        SectionDesc sd;
        { AppendWriter aw; aw.create(path);
          PcoreWriter pw(8, 8, 0.90f);
          pw.add_from_members(7, members);
          sd = pw.finalize(aw, 1, 0); aw.flush(); aw.close(); }

        MmapFileReader mm; mm.open(path);
        PcoreReader r; r.open(mm.data(), sd.file_offset, sd.compressed_size);
        PcoreView v = r.lookup(7);
        require(v.valid() && v.n_members == 1000, "1000 members");
        require(v.total() == 1002, "union: 100, 200, + 1000 singletons");

        std::vector<uint64_t> a; std::vector<float> p; v.materialize(a, p);
        const float p100 = prev_of(a,p,100), p200 = prev_of(a,p,200), p1500 = prev_of(a,p,1500);
        require(p100 > 0.95f, "BUGFIX: aamer-in-all prevalence ~1.0 (was 0.255 capped)");
        require(p200 > 0.4f && p200 < 0.6f, "half-prevalence aamer ~0.5");
        require(p1500 > 0.0f && p1500 < 0.01f, "singleton prevalence ~1/1000");

        std::vector<uint64_t> c; v.core_into(c);
        require(c.size()==1 && c[0]==100, "BUGFIX: core non-empty = {100} (was empty→NaN)");
        require(std::abs(pcore_core_coverage(v, {100}, 0.90f) - 1.0f) < 1e-6f,
                "BUGFIX: large-genus coverage = 1.0 (was NaN)");
    }

    std::cout << "genopack pcore test passed\n";
    return 0;
}
