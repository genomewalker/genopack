#include <genopack/archive.hpp>
#include <genopack/qual.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include "test_utils.hpp"

using namespace genopack;

// Round-trip coverage for `genopack check`: the check subsystem had zero tests,
// and check writes QUAL via an append that rewrites the TOC footer. This guards
// that (a) check runs and writes a QUAL record per genome, and (b) the footer
// survives repeated QUAL-append rewrites so the archive keeps re-opening and
// verifying.
int main() {
    using namespace genopack_test;
    TempDir tmp("genopack_check_rt");

    // Four genomes, two genera (2 members each) so check has a genus reference.
    const auto g1a = make_random_sequence(1, 200000);
    const auto g1b = mutate_every(g1a, 131);
    const auto g2a = make_random_sequence(2, 180000);
    const auto g2b = mutate_every(g2a, 137);

    write_fasta(tmp.path / "g1a.fa", "g1a", g1a);
    write_fasta(tmp.path / "g1b.fa", "g1b", g1b);
    write_fasta(tmp.path / "g2a.fa", "g2a", g2a);
    write_fasta(tmp.path / "g2b.fa", "g2b", g2b);

    // 3-column TSV (accession, taxonomy, file_path) so check can group by genus.
    {
        std::ofstream out(tmp.path / "in.tsv");
        out << "accession\ttaxonomy\tfile_path\n";
        auto row = [&](const char* acc, const char* g, const std::filesystem::path& f) {
            out << acc << "\td__Bacteria;p__P;c__C;o__O;f__F;g__" << g
                << ";s__" << g << " sp\t" << f.string() << "\n";
        };
        row("ACC_G1A", "Genus1", tmp.path / "g1a.fa");
        row("ACC_G1B", "Genus1", tmp.path / "g1b.fa");
        row("ACC_G2A", "Genus2", tmp.path / "g2a.fa");
        row("ACC_G2B", "Genus2", tmp.path / "g2b.fa");
    }

    const std::string bin = GENOPACK_BIN;
    const auto gpk = tmp.path / "a.gpk";
    run_checked(shell_quote(bin) + " build -i " + shell_quote((tmp.path / "in.tsv").string()) +
                " -o " + shell_quote(gpk.string()));

    auto count_qual = [&]() {
        ArchiveReader ar;
        ar.open(gpk);
        require(ar.has_qual(), "no QUAL section present");
        size_t n = 0;
        ar.scan_qual([&](const QualRecord&) { ++n; });
        return n;
    };

    // First check: writes QUAL via an append + footer rewrite.
    run_checked(shell_quote(bin) + " check " + shell_quote(gpk.string()) +
                " -o " + shell_quote((tmp.path / "q.tsv").string()));
    require(count_qual() == 4, "QUAL record count != 4 after first check");

    // Second check (--recompute): repeated append + footer rewrite must still
    // re-open cleanly — this is the footer-consistency regression guard.
    run_checked(shell_quote(bin) + " check " + shell_quote(gpk.string()) +
                " -o " + shell_quote((tmp.path / "q2.tsv").string()) + " --recompute");
    require(count_qual() == 4, "QUAL record count != 4 after second check");

    // The rewritten footer + checksums must still verify.
    run_checked(shell_quote(bin) + " verify " + shell_quote(gpk.string()));

    std::cout << "check round-trip OK\n";
    return 0;
}
