// Unit test for SEC_PROF named reporting/fusion profiles (Phase 5):
//   1. ProfWriter emits profiles + per-axis selections; ProfReader round-trips
//      names, descriptions, axis→column-identity (+ fallback) selections.
//   2. policy_hash is deterministic and order-independent (sorted by axis), and
//      excludes the cosmetic description.
//   3. find(name) locates a stored profile; absent names return -1.

#include <genopack/prof.hpp>
#include <genopack/quality_schema.hpp>
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <iostream>
#include "test_utils.hpp"

using namespace genopack;
using namespace genopack_test;

int main() {
    // "best" profile: completeness from a calibrated genopack column with an
    // intrinsic fallback; contamination from a pinned external column.
    Profile best;
    best.name = "site_best";
    best.description = "calibrated completeness, external contamination";
    best.selections = {
        {std::string(qual_axis::COMPLETENESS),  0xAAAAAAAAAAAAAAAAull, 0xBBBBBBBBBBBBBBBBull},
        {std::string(qual_axis::CONTAMINATION), 0xCCCCCCCCCCCCCCCCull, 0},
    };

    // Same policy, selections listed in the OTHER order, different description.
    Profile reordered;
    reordered.name = "site_best";
    reordered.description = "DIFFERENT cosmetic text";
    reordered.selections = {
        {std::string(qual_axis::CONTAMINATION), 0xCCCCCCCCCCCCCCCCull, 0},
        {std::string(qual_axis::COMPLETENESS),  0xAAAAAAAAAAAAAAAAull, 0xBBBBBBBBBBBBBBBBull},
    };
    require(profile_policy_hash(best) == profile_policy_hash(reordered),
            "policy hash must be order-independent and ignore description");

    Profile intrinsic;
    intrinsic.name = "site_intrinsic";
    intrinsic.description = "";
    intrinsic.selections = {
        {std::string(qual_axis::COMPLETENESS), 0x1111111111111111ull, 0},
    };
    require(profile_policy_hash(intrinsic) != profile_policy_hash(best),
            "distinct policies have distinct hashes");

    TempDir tmp("prof");
    const auto path = tmp.path / "prof.bin";

    SectionDesc sd;
    {
        AppendWriter aw;
        aw.create(path);
        ProfWriter pw;
        pw.add(best);
        pw.add(intrinsic);
        sd = pw.finalize(aw, /*section_id=*/9);
        aw.flush();
        aw.close();
    }
    require(sd.type == SEC_PROF, "section type must be SEC_PROF");
    require(sd.item_count == 2, "two profiles written");

    MmapFileReader mm;
    mm.open(path);
    ProfReader r;
    r.open(mm.data(), sd.file_offset, sd.compressed_size);

    require(r.n_profiles() == 2, "two profiles read back");

    const int bi = r.find("site_best");
    require(bi >= 0, "find(site_best) must locate the profile");
    require(r.find("nope") < 0, "absent name returns -1");

    Profile got = r.profile(static_cast<uint32_t>(bi));
    require(got.name == "site_best", "name round-trips");
    require(got.description == "calibrated completeness, external contamination",
            "description round-trips");
    require(got.selections.size() == 2, "two selections round-trip");

    bool saw_comp = false, saw_cont = false;
    for (const auto& s : got.selections) {
        if (s.axis == qual_axis::COMPLETENESS) {
            saw_comp = true;
            require(s.column_identity == 0xAAAAAAAAAAAAAAAAull, "completeness identity");
            require(s.fallback_identity == 0xBBBBBBBBBBBBBBBBull, "completeness fallback identity");
        }
        if (s.axis == qual_axis::CONTAMINATION) {
            saw_cont = true;
            require(s.column_identity == 0xCCCCCCCCCCCCCCCCull, "contamination identity");
            require(s.fallback_identity == 0, "contamination has no fallback");
        }
    }
    require(saw_comp && saw_cont, "both axes present");

    require(r.policy_hash(static_cast<uint32_t>(bi)) == profile_policy_hash(best),
            "stored policy hash equals recomputed");

    std::cout << "genopack prof test passed\n";
    return 0;
}
