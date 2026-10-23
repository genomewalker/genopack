// compute_ani's af is containment scaled so skani's GTDB cut AF >= 0.5 falls at
// containment 0.405, capped at 1; ANI and the directional containments are unscaled.
#include <genopack/skani.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "test_utils.hpp"

using namespace genopack;

// 1000-hash sketch a and sketch b sharing `shared` of them; b has `b_size` hashes.
static SkaniResult pair_with(size_t shared, size_t b_size) {
    SkaniSketch a, b;
    for (uint64_t h = 0; h < 1000; ++h) a.hashes.push_back(2 * h);
    for (uint64_t h = 0; h < shared; ++h) b.hashes.push_back(2 * h);
    for (uint64_t h = 0; b.hashes.size() < b_size; ++h) b.hashes.push_back(1000000 + 2 * h + 1);
    std::sort(b.hashes.begin(), b.hashes.end());
    return compute_ani(a, b, 21);
}

int main() {
    using genopack_test::require;
    const auto above = pair_with(406, 5000), below = pair_with(404, 5000);
    std::cerr << "c=0.406 af " << above.af << ", c=0.404 af " << below.af << "\n";
    require(above.af >= 0.5, "containment 0.406 must pass the GTDB AF cut");
    require(below.af < 0.5, "containment 0.404 must fail the GTDB AF cut");
    require(std::abs(above.c_ab - 0.406) < 1e-12 && std::abs(above.c_ba - 0.406 / 5) < 1e-12, "containments are unscaled");
    require(std::abs(above.ani - 100 * std::pow(0.406, 1.0 / 21)) < 1e-9, "ANI is unscaled");
    const auto high = pair_with(900, 1000);
    require(high.af == 1.0, "af is capped at 1");
    require(pair_with(0, 1000).af == 0.0, "no shared hashes means af 0");
    std::cout << "PASS\n";
}
