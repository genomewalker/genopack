#include <genopack/archive.hpp>
#include <genopack/mem_delta.hpp>
#include <genopack/shard.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include "test_utils.hpp"

using namespace genopack;

int main() {
    using namespace genopack_test;
    TempDir tmp("genopack_smoke");

    const std::string a1 = make_sequence("ACGT", 140000);
    const std::string a2 = mutate_every(a1, 97);
    const std::string b1 = make_sequence("TTGCAACG", 90000);
    const std::string mc1 = make_sequence("AACCGGTT", 30000);
    const std::string mc2 = mutate_every(make_sequence("TTAACCGG", 25000), 53);
    const std::string mc_seq = mc1 + mc2;

    write_fasta(tmp.path / "a1.fa", "a1", a1);
    write_fasta(tmp.path / "a2.fa", "a2", a2);
    write_fasta(tmp.path / "b1.fa", "b1", b1);
    {
        std::ofstream out(tmp.path / "mc.fa");
        out << ">mc1\n" << mc1 << "\n>mc2\n" << mc2 << "\n";
    }

    write_tsv(tmp.path / "a.tsv", {
        {"ACC_A1", tmp.path / "a1.fa"},
        {"ACC_A2", tmp.path / "a2.fa"},
    });
    write_tsv(tmp.path / "b.tsv", {
        {"ACC_B1", tmp.path / "b1.fa"},
    });
    write_tsv(tmp.path / "mc.tsv", {
        {"ACC_MC", tmp.path / "mc.fa"},
    });

    const std::string bin = GENOPACK_BIN;
    run_checked(shell_quote(bin) + " build -i " + shell_quote((tmp.path / "a.tsv").string()) +
                " -o " + shell_quote((tmp.path / "a.gpk").string()) + " --delta");
    run_checked(shell_quote(bin) + " build -i " + shell_quote((tmp.path / "b.tsv").string()) +
                " -o " + shell_quote((tmp.path / "b.gpk").string()));
    run_checked(shell_quote(bin) + " build -i " + shell_quote((tmp.path / "mc.tsv").string()) +
                " -o " + shell_quote((tmp.path / "mc.gpk").string()) + " --delta");

    {
        ArchiveReader ar;
        ar.open(tmp.path / "a.gpk");
        auto slice = ar.fetch_sequence_slice_by_accession("ACC_A2", 65000, 200);
        require(slice.has_value(), "missing slice for ACC_A2");
        require(*slice == a2.substr(65000, 200), "delta slice mismatch");

        bool saw_checkpointed = false;
        ar.scan_shards([&](const uint8_t* data, uint64_t offset, uint64_t size, uint32_t) {
            ShardReader sr;
            sr.open(data, offset, size);
            for (const auto* de = sr.dir_begin(); de != sr.dir_end(); ++de) {
                if (de->n_checkpoints > 1) {
                    saw_checkpointed = true;
                    break;
                }
            }
        });
        require(saw_checkpointed, "expected chunk checkpoints in built shard");
    }

    {
        ArchiveReader ar;
        ar.open(tmp.path / "mc.gpk");
        auto genome = ar.fetch_by_accession("ACC_MC");
        require(genome.has_value(), "missing multi-contig genome");
        auto fc = extract_fasta_components(genome->fasta.data(), genome->fasta.size());
        require(fc.seq == mc_seq, "multi-contig sequence reconstruction mismatch");
        require(fc.contig_ends.size() == 2, "multi-contig contig count mismatch");
        require(fc.contig_ends[0] == mc1.size(), "multi-contig first boundary mismatch");
        require(fc.contig_ends[1] == mc_seq.size(), "multi-contig final boundary mismatch");

        auto slice = ar.fetch_sequence_slice_by_accession("ACC_MC", mc1.size() - 50, 120);
        require(slice.has_value(), "missing multi-contig slice");
        require(*slice == mc_seq.substr(mc1.size() - 50, 120), "multi-contig slice mismatch");
    }

    run_checked(shell_quote(bin) + " rm " + shell_quote((tmp.path / "a.gpk").string()) + " ACC_A2");
    run_checked(shell_quote(bin) + " merge " + shell_quote((tmp.path / "a.gpk").string()) + " " +
                shell_quote((tmp.path / "b.gpk").string()) + " -o " +
                shell_quote((tmp.path / "merged.gpk").string()));

    {
        ArchiveReader ar;
        ar.open(tmp.path / "merged.gpk");
        auto stats = ar.archive_stats();
        require(stats.n_genomes_total == 3, "merged total count mismatch");
        require(stats.n_genomes_live == 2, "merged live count mismatch");
        require(!ar.fetch_by_accession("ACC_A2").has_value(), "deleted accession resurrected");

        auto slice = ar.fetch_sequence_slice_by_accession("ACC_B1", 12345, 180);
        require(slice.has_value(), "missing slice for ACC_B1");
        require(*slice == b1.substr(12345, 180), "merged slice mismatch");

    }

    run_checked(shell_quote(bin) + " reindex " + shell_quote((tmp.path / "merged.gpk").string()) + " --force");

    {
        ArchiveReader ar;
        ar.open(tmp.path / "merged.gpk");
        auto stats = ar.archive_stats();
        require(stats.generation == 2, "reindex did not advance generation");
        auto slice = ar.fetch_sequence_slice_by_accession("ACC_A1", 70000, 128);
        require(slice.has_value(), "missing post-reindex slice");
        require(*slice == a1.substr(70000, 128), "post-reindex slice mismatch");
    }

    run_checked(shell_quote(bin) + " slice " + shell_quote((tmp.path / "merged.gpk").string()) +
                " ACC_A1 --start 70000 --length 128 --fasta > " +
                shell_quote((tmp.path / "slice.out").string()));
    std::string slice_out = read_text(tmp.path / "slice.out");
    require(slice_out.find(a1.substr(70000, 128)) != std::string::npos,
            "slice CLI output mismatch");

    std::cout << "genopack smoke test passed\n";
    return 0;
}
