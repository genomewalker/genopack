#include <genopack/archive.hpp>
#include <genopack/mem_delta.hpp>
#include <genopack/shard.hpp>
#include <iostream>
#include <string>
#include "test_utils.hpp"

using namespace genopack;

int main() {
    using namespace genopack_test;

    const std::string anchor_seq = make_sequence("ACGTTGCA", 220000);
    const std::string query_seq  = mutate_every(anchor_seq, 101);

    FastaComponents anchor_fc;
    anchor_fc.seq = anchor_seq;
    anchor_fc.contig_ends = {static_cast<uint32_t>(anchor_seq.size())};
    anchor_fc.headers = {"anchor"};

    FastaComponents query_fc;
    query_fc.seq = query_seq;
    query_fc.contig_ends = {static_cast<uint32_t>(query_seq.size())};
    query_fc.headers = {"query"};

    AnchorIndex anchor_idx = build_anchor_index(anchor_fc.seq);
    auto blob = encode_mem_delta(anchor_fc, anchor_idx, query_fc, 32768);
    std::string decoded = decode_mem_delta(anchor_fc.seq, blob.data(), blob.size());
    auto decoded_fc = extract_fasta_components(decoded.data(), decoded.size());
    require(decoded_fc.seq == query_seq, "full mem-delta decode mismatch");

    std::string slice = decode_mem_delta_slice(anchor_fc.seq, blob.data(), blob.size(), 70000, 4096);
    require(slice == query_seq.substr(70000, 4096), "mem-delta slice mismatch");

    TempDir tmp("genopack_mem_delta");
    const std::string related_0 = mutate_every(anchor_seq, 131);
    const std::string related_1 = mutate_every(anchor_seq, 149);
    const std::string related_2 = mutate_every(anchor_seq, 167);
    const std::string related_3 = mutate_every(anchor_seq, 181);
    const std::string related_4 = mutate_every(anchor_seq, 193);
    const std::string related_tail = mutate_every(anchor_seq, 211);
    const std::string mem_delta_candidate = mutate_every(anchor_seq, 97);
    const std::string fallback_candidate = mutate_every(anchor_seq, 37);

    write_fasta(tmp.path / "r0.fa", "r0", related_0);
    write_fasta(tmp.path / "r1.fa", "r1", related_1);
    write_fasta(tmp.path / "r2.fa", "r2", related_2);
    write_fasta(tmp.path / "r3.fa", "r3", related_3);
    write_fasta(tmp.path / "r4.fa", "r4", related_4);
    write_fasta(tmp.path / "r5.fa", "r5", related_tail);
    write_fasta(tmp.path / "mem.fa", "mem", mem_delta_candidate);
    write_fasta(tmp.path / "plain.fa", "plain", fallback_candidate);
    write_tsv(tmp.path / "input.tsv", {
        {"ACC_R0", tmp.path / "r0.fa"},
        {"ACC_R1", tmp.path / "r1.fa"},
        {"ACC_R2", tmp.path / "r2.fa"},
        {"ACC_R3", tmp.path / "r3.fa"},
        {"ACC_R4", tmp.path / "r4.fa"},
        {"ACC_MEM", tmp.path / "mem.fa"},
        {"ACC_PLAIN", tmp.path / "plain.fa"},
        {"ACC_R5", tmp.path / "r5.fa"},
    });

    const std::string bin = GENOPACK_BIN;
    run_checked(shell_quote(bin) + " build -i " + shell_quote((tmp.path / "input.tsv").string()) +
                " -o " + shell_quote((tmp.path / "db.gpk").string()) + " --mem-delta --no-hnsw");

    ArchiveReader ar;
    ar.open(tmp.path / "db.gpk");
    auto mem_slice = ar.fetch_sequence_slice_by_accession("ACC_MEM", 65536, 8192);
    require(mem_slice.has_value(), "archive mem-delta slice missing");
    require(*mem_slice == mem_delta_candidate.substr(65536, 8192), "archive mem-delta slice mismatch");

    auto plain_slice = ar.fetch_sequence_slice_by_accession("ACC_PLAIN", 65536, 8192);
    require(plain_slice.has_value(), "archive fallback slice missing");
    require(*plain_slice == fallback_candidate.substr(65536, 8192), "archive fallback slice mismatch");

    auto anchor_slice = ar.fetch_sequence_slice_by_accession("ACC_R0", 65536, 8192);
    require(anchor_slice.has_value(), "archive anchor slice missing");
    require(*anchor_slice == related_0.substr(65536, 8192), "archive anchor slice mismatch");

    size_t mem_delta_nonref = 0;
    size_t plain_total = 0;
    ar.scan_shards([&](const uint8_t* data, uint64_t offset, uint64_t size, uint32_t) {
        const auto* hdr = reinterpret_cast<const ShardHeader*>(data + offset);
        require(hdr->codec == 4, "expected MEM-delta shard");
        ShardReader sr;
        sr.open(data, offset, size);
        for (uint32_t i = 0; i < hdr->dict_size; ++i) {
            const auto* de = sr.dir_entry(i);
            require(de != nullptr, "missing anchor dir entry");
            if ((de->flags & GenomeMeta::FLAG_DELTA) == 0) ++plain_total;
        }
        for (uint32_t i = hdr->dict_size; i < sr.n_genomes(); ++i) {
            const auto* de = sr.dir_entry(i);
            require(de != nullptr, "missing shard dir entry");
            if (de->flags & GenomeMeta::FLAG_DELTA) ++mem_delta_nonref;
            else ++plain_total;
        }
    });
    require(mem_delta_nonref + plain_total == 8, "unexpected codec accounting");
    require(plain_total >= 1, "expected at least one plain genome in codec-4 shard");

    std::cout << "genopack mem-delta test passed\n";
    return 0;
}
