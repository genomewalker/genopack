// Unit test for SEC_QCONTIG per-contig overlay (Phase 6b):
//   1. QcontigWriter sorts records by (genome_id, contig_index); finalize emits
//      a compact record section.
//   2. QcontigReader::contigs_for binary-searches a genome's contiguous range.
//   3. Per-contig fields (offset/length/tnf/leakage) round-trip, NaN preserved.

#include <genopack/qcontig.hpp>
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <cmath>
#include <iostream>
#include "test_utils.hpp"

using namespace genopack;
using namespace genopack_test;

struct TestFlag {
    uint32_t contig_offset;
    uint32_t contig_length;
    float    tnf_score;
    float    leakage_score;
    float    gcov_t2_pct;
    float    gcov_spe_pct;
    float    prot_foreign_frac;
    float    prot_loglr;
    uint16_t prot_classifiable;
    uint16_t prot_host_specific;
    uint16_t prot_foreign_specific;
    uint16_t prot_family_hits;
    uint32_t prot_best_genus;
    uint8_t  prot_flags;
};

int main() {
    // Add genomes out of order to exercise the sort. Last field group = protein channel.
    std::vector<TestFlag> g20 = {
        {0, 5000, 0.9f, NAN, 0.4f, 0.99f, 0.80f, 4.0f, 120, 24, 96, 5, 777, 0},   // foreign-heavy
        {5000, 3000, 0.1f, 0.2f, 0.1f, 0.05f, 0.02f, -1.0f, 90, 88, 2, 60, 0, 0}};
    std::vector<TestFlag> g10 = {
        {0, 8000, NAN, NAN, NAN, NAN, NAN, NAN, 0, 0, 0, 0, 0, 0x1},               // abstain
        {8000, 2000, 0.5f, 0.7f, 0.3f, 0.8f, 0.5f, 2.0f, 40, 20, 20, 3, 42, 0},
        {10000, 1500, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 35, 35, 0, 1, 0, 0}};

    TempDir tmp("qcontig");
    const auto path = tmp.path / "qcontig.bin";

    SectionDesc sd;
    {
        AppendWriter aw;
        aw.create(path);
        QcontigWriter qw;
        qw.add_genome(20, g20);
        qw.add_genome(10, g10);
        sd = qw.finalize(aw, /*section_id=*/3);
        aw.flush();
        aw.close();
    }
    require(sd.type == SEC_QCONTIG, "section type must be SEC_QCONTIG");
    require(sd.item_count == 5, "5 contig records total");

    MmapFileReader mm;
    mm.open(path);
    QcontigReader r;
    r.open(mm.data(), sd.file_offset, sd.compressed_size);
    require(r.n_rows() == 5, "5 rows");

    auto c10 = r.contigs_for(10);
    require(c10.size() == 3, "genome 10 has 3 contigs");
    require(c10[0].contig_index == 0 && c10[1].contig_index == 1 && c10[2].contig_index == 2,
            "contig indices ascending within genome");
    require(c10[0].contig_length == 8000, "g10 contig0 length");
    require(std::isnan(c10[0].tnf_score), "NaN tnf preserved");
    require(std::isnan(c10[0].gcov_spe_pct), "NaN spe preserved");
    require(std::abs(c10[1].leakage_score - 0.7f) < 1e-6f, "g10 contig1 leakage");
    require(std::abs(c10[1].gcov_spe_pct - 0.8f) < 1e-6f, "g10 contig1 SPE percentile");

    auto c20 = r.contigs_for(20);
    require(c20.size() == 2, "genome 20 has 2 contigs");
    require(c20[0].contig_length == 5000 && c20[1].contig_offset == 5000, "g20 fields");
    require(std::abs(c20[0].tnf_score - 0.9f) < 1e-6f, "g20 contig0 tnf");
    require(std::abs(c20[0].gcov_spe_pct - 0.99f) < 1e-6f, "g20 contig0 SPE percentile (foreign-like)");
    // Protein foreign-containment channel round-trips.
    require(std::abs(c20[0].prot_foreign_frac - 0.80f) < 1e-6f, "g20 contig0 foreign fraction");
    require(c20[0].prot_foreign_specific == 96 && c20[0].prot_classifiable == 120, "g20 contig0 prot counts");
    require(c20[0].prot_best_genus == 777, "g20 contig0 best foreign genus tag");
    require(c10[0].prot_flags == 0x1, "g10 contig0 ABSTAIN_LOW_N flag preserved");

    require(r.contigs_for(99).empty(), "absent genome yields empty range");

    std::cout << "genopack qcontig test passed\n";
    return 0;
}
