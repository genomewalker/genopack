// Unit test for SEC_PCORE — the unified dense prevalence-annotated reference.
//   1. PcoreWriter stores every aamer (prev>=1) with a member COUNT; PcoreReader
//      round-trips aamers + counts + n_members.
//   2. pcore_core_coverage reproduces the conserved-core (>=ceil(theta*n)) coverage
//      that CORE computed — the parity guarantee for retiring CORE.

#include <genopack/pcore.hpp>
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <cmath>
#include <iostream>
#include "test_utils.hpp"

using namespace genopack;
using namespace genopack_test;

int main() {
    // 3 members. aamer prevalence counts: 1->3, 2->3, 3->2, 4->1.
    std::vector<std::vector<uint64_t>> members = {{1,2,3,4}, {1,2,3}, {1,2}};

    TempDir tmp("pcore");
    const auto path = tmp.path / "pcore.bin";
    SectionDesc sd;
    {
        AppendWriter aw; aw.create(path);
        PcoreWriter pw(8, 8, /*theta=*/0.90f);
        pw.add_from_members(42, members);
        sd = pw.finalize(aw, /*section_id=*/5, /*frac_max=*/123);
        aw.flush(); aw.close();
    }
    require(sd.type == SEC_PCORE, "section type must be SEC_PCORE");

    MmapFileReader mm; mm.open(path);
    PcoreReader r; r.open(mm.data(), sd.file_offset, sd.compressed_size);
    require(r.n_entries() == 1, "one entry");
    require(std::abs(r.theta() - 0.90f) < 1e-6f, "theta round-trips");
    require(r.frac_max_hash() == 123, "frac_max round-trips");

    PcoreView v = r.lookup(42);
    require(v.valid() && v.n_aamers == 4 && v.n_members == 3, "dense union of 4 aamers, 3 members");
    require(v.aamers[0]==1 && v.aamers[1]==2 && v.aamers[2]==3 && v.aamers[3]==4, "aamers sorted");
    require(v.prev[0]==3 && v.prev[1]==3 && v.prev[2]==2 && v.prev[3]==1, "member counts (not percent)");
    require(r.lookup(99).valid() == false, "absent key → invalid view");

    // theta=0.90 → threshold = ceil(2.7) = 3 → conserved core = {1,2}.
    require(pcore_core_threshold(3, 0.90f) == 3, "ceil(0.9*3)=3");
    {
        std::vector<uint64_t> q = {1, 2, 9};     // both core aamers present
        require(std::abs(pcore_core_coverage(v, q, 0.90f) - 1.0f) < 1e-6f, "coverage 2/2 = 1.0");
    }
    {
        std::vector<uint64_t> q = {1, 9};        // only one core aamer
        require(std::abs(pcore_core_coverage(v, q, 0.90f) - 0.5f) < 1e-6f, "coverage 1/2 = 0.5");
    }
    {
        // theta=0.50 → threshold = ceil(1.5) = 2 → core = {1,2,3}.
        require(pcore_core_threshold(3, 0.50f) == 2, "ceil(0.5*3)=2");
        std::vector<uint64_t> q = {1, 2, 9};     // 2 of the 3 core aamers
        require(std::abs(pcore_core_coverage(v, q, 0.50f) - (2.0f/3.0f)) < 1e-6f, "coverage 2/3");
    }

    std::cout << "genopack pcore test passed\n";
    return 0;
}
