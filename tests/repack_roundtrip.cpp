// Regression coverage for `genopack repack`.
//
// repack re-shards an archive by taxonomy. No test covered it, yet it is the
// target of a deferred raw-blob-copy optimization. This pins the invariants a
// correct repack must preserve: genome count, every genome's full sequence
// (byte-for-byte), and per-accession taxonomy.
#include <genopack/archive.hpp>
#include <genopack/mem_delta.hpp>
#include <iostream>
#include <map>
#include <vector>
#include "test_utils.hpp"

using namespace genopack;
using namespace genopack_test;

namespace {

std::string seq_of(ArchiveReader& ar, const std::string& acc) {
    auto g = ar.fetch_by_accession(acc);
    require(g.has_value(), "missing genome: " + acc);
    auto fc = extract_fasta_components(g->fasta.data(), g->fasta.size());
    return fc.seq;
}

} // namespace

int main() {
    TempDir tmp("genopack_repack");
    const std::string bin = GENOPACK_BIN;

    // Genomes across two genera so a by-genus repack actually regroups shards.
    struct Rec { std::string acc, file, lin, seq; };
    std::vector<Rec> recs;
    auto add = [&](const std::string& acc, const std::string& genus,
                   const std::string& fam, uint64_t seed, size_t n) {
        std::string s = make_random_sequence(seed, n);
        std::string f = (tmp.path / (acc + ".fa")).string();
        write_fasta(f, acc + "_ctg", s);
        std::string lin = "d__Bacteria;p__P;c__C;o__O;f__" + fam + ";g__" + genus +
                          ";s__" + genus + " " + acc;
        recs.push_back({acc, f, lin, s});
    };
    add("ACC_E1", "Escherichia", "Enterobacteriaceae", 1, 150000);
    add("ACC_E2", "Escherichia", "Enterobacteriaceae", 2, 130000);
    add("ACC_S1", "Salmonella",  "Enterobacteriaceae", 3, 140000);
    add("ACC_B1", "Bacillus",    "Bacillaceae",        4, 120000);

    {
        std::ofstream out(tmp.path / "all.tsv");
        out << "accession\tfile_path\ttaxonomy\n";
        for (const auto& r : recs) out << r.acc << "\t" << r.file << "\t" << r.lin << "\n";
    }

    run_checked(shell_quote(bin) + " build -i " + shell_quote((tmp.path / "all.tsv").string()) +
                " -o " + shell_quote((tmp.path / "src.gpk").string()));

    // Capture source truth.
    std::map<std::string, std::string> src_seq, src_tax;
    {
        ArchiveReader ar; ar.open(tmp.path / "src.gpk");
        for (const auto& r : recs) {
            src_seq[r.acc] = seq_of(ar, r.acc);
            auto t = ar.taxonomy_for_accession(r.acc);
            require(t.has_value(), "missing taxonomy in source: " + r.acc);
            src_tax[r.acc] = *t;
        }
    }

    run_checked(shell_quote(bin) + " repack " + shell_quote((tmp.path / "src.gpk").string()) +
                " " + shell_quote((tmp.path / "repacked.gpk").string()) + " --taxonomy-rank g");

    ArchiveReader rp; rp.open(tmp.path / "repacked.gpk");
    auto stats = rp.archive_stats();
    require(stats.n_genomes_total == recs.size(), "repack changed genome count");
    require(stats.n_genomes_live  == recs.size(), "repack changed live count");

    for (const auto& r : recs) {
        require(seq_of(rp, r.acc) == src_seq[r.acc],
                "sequence changed across repack: " + r.acc);
        require(src_seq[r.acc] == r.seq, "source sequence mismatch (sanity): " + r.acc);
        auto t = rp.taxonomy_for_accession(r.acc);
        require(t.has_value() && *t == src_tax[r.acc],
                "taxonomy changed across repack: " + r.acc);
    }

    std::cout << "genopack repack roundtrip test passed\n";
    return 0;
}
