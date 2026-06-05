// Unit test for SEC_CORE per-genus prevalence cores (Phase 1):
//   1. add_from_members builds the prevalence core (aamers in >= ceil(theta*N)
//      members); direct add() stores a precomputed core.
//   2. core_model_hash is deterministic (insertion-order independent), folds the
//      params (theta), and changes when a core changes.
//   3. Round-trip through AppendWriter/CoreReader: genus lookup returns the sorted
//      core aamers + sorted member ids; header carries theta/k/frac/model hash.

#include <genopack/core_section.hpp>
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <iostream>
#include <vector>
#include "test_utils.hpp"

using namespace genopack;
using namespace genopack_test;

static constexpr uint64_t GA = 0xAAAA0001ull;  // genus A hash
static constexpr uint64_t GB = 0xBBBB0002ull;  // genus B hash

int main() {
    // genus A: 3 members; theta 0.90 -> ceil(2.7)=3 -> core = aamers in ALL members
    const std::vector<std::vector<uint64_t>> A_members = {
        {11, 100, 200, 300},
        {22, 100, 200, 300},
        {33, 100, 200, 300},
    };
    const std::vector<uint64_t> A_ids = {30, 10, 20};  // stored sorted

    // ── core_model_hash determinism + sensitivity ──
    uint64_t h_order1, h_order2;
    {
        CoreWriter w1(8, 8, 0.90f);
        w1.add_from_members(GA, A_members, A_ids);
        w1.add(GB, {500, 600}, {8, 7});
        h_order1 = w1.core_model_hash();

        CoreWriter w2(8, 8, 0.90f);            // same content, reversed insertion order
        w2.add(GB, {500, 600}, {7, 8});
        w2.add_from_members(GA, A_members, A_ids);
        h_order2 = w2.core_model_hash();
    }
    require(h_order1 == h_order2, "core_model_hash must be insertion-order independent");
    {
        CoreWriter wt(8, 8, 0.50f);            // different theta, same direct core
        wt.add(GA, {100, 200, 300}, {1});
        CoreWriter wt2(8, 8, 0.90f);
        wt2.add(GA, {100, 200, 300}, {1});
        require(wt.core_model_hash() != wt2.core_model_hash(),
                "core_model_hash must fold theta");
        CoreWriter wc(8, 8, 0.90f);
        wc.add(GA, {100, 200, 301}, {1});      // one aamer differs
        require(wc.core_model_hash() != wt2.core_model_hash(),
                "core_model_hash must change when a core changes");
    }

    // ── write + read back ──
    TempDir tmp("core_section");
    const auto path = tmp.path / "core.bin";
    const uint64_t FRAC = 0x00000000FFFFFFFFull;

    SectionDesc sd;
    uint64_t model_hash;
    {
        CoreWriter w(8, 8, 0.90f);
        w.add_from_members(GA, A_members, A_ids);
        w.add(GB, {500, 600}, {8, 7});
        model_hash = w.core_model_hash();
        require(w.n_genera() == 2, "two genera expected");

        AppendWriter aw;
        aw.create(path);
        sd = w.finalize(aw, /*section_id=*/3, FRAC);
        aw.flush();
        aw.close();
    }
    require(sd.type == SEC_CORE, "section type mismatch");
    require(sd.item_count == 2, "item_count should equal genus count");
    require(sd.aux0 == model_hash, "aux0 should carry core_model_hash");

    {
        MmapFileReader mm;
        mm.open(path);
        CoreReader r;
        r.open(mm.data(), sd.file_offset, sd.compressed_size);

        require(r.n_genera() == 2, "n_genera mismatch");
        require(r.k() == 8, "k mismatch");
        require(r.theta() == 0.90f, "theta mismatch");
        require(r.frac_max_hash() == FRAC, "frac_max_hash mismatch");
        require(r.core_model_hash() == model_hash, "core_model_hash mismatch");

        CoreView a = r.lookup(GA);
        require(a.valid(), "genus A must be found");
        require(a.n_aamers == 3, "genus A core size should be 3");
        require(a.aamers[0] == 100 && a.aamers[1] == 200 && a.aamers[2] == 300,
                "genus A core aamers mismatch (sorted prevalence core)");
        require(a.n_members == 3, "genus A member count mismatch");
        require(a.members[0] == 10 && a.members[1] == 20 && a.members[2] == 30,
                "genus A member ids must be sorted");

        CoreView b = r.lookup(GB);
        require(b.valid(), "genus B must be found");
        require(b.n_aamers == 2 && b.aamers[0] == 500 && b.aamers[1] == 600,
                "genus B core mismatch");
        require(b.n_members == 2 && b.members[0] == 7 && b.members[1] == 8,
                "genus B member ids mismatch");

        require(!r.lookup(0xDEADBEEFull).valid(), "absent genus must be invalid");
    }

    std::cout << "genopack core-section test passed\n";
    return 0;
}
