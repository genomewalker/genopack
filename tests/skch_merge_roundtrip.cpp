// Regression coverage for single-k SKCH merge.
//
// The existing smoke test only exercises *multi-k* SKCH merge. The single-k
// merge path (merger.cpp, !multi_k branch) had no test, which is exactly the
// path a future streaming-SKCH-merge rewrite would touch. This test pins the
// invariant that merge preserves every genome's OPH signature byte-for-byte
// (sketches are content-derived; merge only copies/re-emits them, so a correct
// merge must not alter sig/sig2/mask).
#include <genopack/archive.hpp>
#include <genopack/skch.hpp>
#include <iostream>
#include <map>
#include <vector>
#include "test_utils.hpp"

using namespace genopack;
using namespace genopack_test;

namespace {

struct Sketch {
    std::vector<uint16_t> sig, sig2;
    std::vector<uint64_t> mask;
    uint32_t              n_real_bins = 0;
    uint64_t             genome_length = 0;
};

Sketch sketch_of(ArchiveReader& ar, const std::string& acc) {
    auto meta = ar.genome_meta_by_accession(acc);
    require(meta.has_value(), "accession not found: " + acc);
    auto ks = ar.available_sketch_kmer_sizes();
    require(ks.size() == 1, "expected exactly one sketch k-mer size (single-k)");
    auto sk = ar.sketch_for(meta->genome_id, ks[0], ar.sketch_sketch_size());
    require(sk.has_value(), "no sketch for accession: " + acc);
    Sketch out;
    out.sig.assign(sk->sig, sk->sig + sk->sketch_size);
    out.sig2.assign(sk->sig2, sk->sig2 + sk->sketch_size);
    out.mask.assign(sk->mask, sk->mask + sk->mask_words);
    out.n_real_bins   = sk->n_real_bins;
    out.genome_length = sk->genome_length;
    return out;
}

} // namespace

int main() {
    TempDir tmp("genopack_skch_merge");
    const std::string bin = GENOPACK_BIN;

    // Two near-identical genomes + one distinct, split across two archives so the
    // merge actually concatenates two single-k SKCH sections.
    const std::string g1 = make_random_sequence(11, 120000);
    const std::string g2 = mutate_every(g1, 800);
    const std::string g3 = make_random_sequence(97, 90000);
    write_fasta(tmp.path / "g1.fa", "g1", g1);
    write_fasta(tmp.path / "g2.fa", "g2", g2);
    write_fasta(tmp.path / "g3.fa", "g3", g3);

    write_tsv(tmp.path / "a.tsv", {{"ACC1", tmp.path / "g1.fa"},
                                   {"ACC2", tmp.path / "g2.fa"}});
    write_tsv(tmp.path / "b.tsv", {{"ACC3", tmp.path / "g3.fa"}});

    auto build = [&](const char* tsv, const char* out) {
        run_checked(shell_quote(bin) + " build -i " +
                    shell_quote((tmp.path / tsv).string()) + " -o " +
                    shell_quote((tmp.path / out).string()) + " --sketch-kmers 16");
    };
    build("a.tsv", "a.gpk");
    build("b.tsv", "b.gpk");

    std::map<std::string, Sketch> src;
    { ArchiveReader a; a.open(tmp.path / "a.gpk");
      src["ACC1"] = sketch_of(a, "ACC1");
      src["ACC2"] = sketch_of(a, "ACC2"); }
    { ArchiveReader b; b.open(tmp.path / "b.gpk");
      src["ACC3"] = sketch_of(b, "ACC3"); }

    run_checked(shell_quote(bin) + " merge " +
                shell_quote((tmp.path / "a.gpk").string()) + " " +
                shell_quote((tmp.path / "b.gpk").string()) + " -o " +
                shell_quote((tmp.path / "merged.gpk").string()));

    ArchiveReader m;
    m.open(tmp.path / "merged.gpk");
    require(m.available_sketch_kmer_sizes().size() == 1,
            "merged archive lost single-k SKCH section");

    for (const auto& acc : {"ACC1", "ACC2", "ACC3"}) {
        Sketch got = sketch_of(m, acc);
        const Sketch& want = src[acc];
        require(got.sig == want.sig,   std::string("sig changed across merge: ") + acc);
        require(got.sig2 == want.sig2, std::string("sig2 changed across merge: ") + acc);
        require(got.mask == want.mask, std::string("mask changed across merge: ") + acc);
        require(got.n_real_bins == want.n_real_bins,
                std::string("n_real_bins changed across merge: ") + acc);
        require(got.genome_length == want.genome_length,
                std::string("genome_length changed across merge: ") + acc);
    }

    std::cout << "genopack skch merge roundtrip test passed\n";
    return 0;
}
