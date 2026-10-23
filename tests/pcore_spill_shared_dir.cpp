// Two PcoreSpillWriters sharing a spill dir (as concurrent builds share /tmp) must
// not clash. A 1-tuple buffer flushes every genome to disk, and both genus hashes
// route to the same partition, so shared file names would truncate each other.
#include <genopack/pcore_spiller.hpp>
#include <iostream>
#include <map>
#include <set>
#include "test_utils.hpp"

using namespace genopack;

using Seen = std::map<uint64_t, std::set<uint64_t>>;

static Seen drain(PcoreSpillWriter& w) {
    Seen seen;
    w.finalize(1, [&](uint64_t genus, std::vector<uint64_t>& aamers, std::vector<uint32_t>&, uint32_t) {
        seen[genus].insert(aamers.begin(), aamers.end());
    });
    return seen;
}

int main() {
    using namespace genopack_test;
    TempDir tmp("genopack_pcore_spill_shared");
    const uint64_t ga = 0x1234000000000001ULL, gb = 0x1234000000000002ULL;  // same top 13 bits
    Seen a_seen, b_seen;
    {
        PcoreSpillWriter a(tmp.path, 1), b(tmp.path, 1);
        a.add_genome(ga, std::vector<uint64_t>{10, 20, 30}); a.add_member(ga);
        b.add_genome(gb, std::vector<uint64_t>{40, 50});     b.add_member(gb);
        a.add_genome(ga, std::vector<uint64_t>{11, 21});     a.add_member(ga);
        b.add_genome(gb, std::vector<uint64_t>{41});         b.add_member(gb);
        a_seen = drain(a);
        b_seen = drain(b);
    }
    require(a_seen == Seen{{ga, {10, 11, 20, 21, 30}}}, "writer A lost or mixed its spill data");
    require(b_seen == Seen{{gb, {40, 41, 50}}}, "writer B lost or mixed its spill data");
    for (const auto& e : std::filesystem::directory_iterator(tmp.path))
        require(e.path().filename().string().rfind("pspill", 0) != 0, "spill leftovers: " + e.path().string());
    std::cout << "PASS\n";
}
