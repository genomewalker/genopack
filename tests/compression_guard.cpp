// Compression-ratio guard (relative).
//
// Codec selection (MEM-delta vs plain, kmer-NN ordering) is a build-time quality
// decision the roundtrip tests cannot observe: a regression that falls back to
// plain zstd still decodes correctly, it just bloats the archive. This builds a
// set of near-identical genomes twice — once forced plain, once with MEM-delta —
// and asserts MEM-delta beats plain by a wide margin. The check is relative, so
// it does not drift with zstd versions; it fails only if MEM-delta stops winning.
#include <genopack/archive.hpp>
#include <iostream>
#include <vector>
#include "test_utils.hpp"

using namespace genopack;
using namespace genopack_test;

namespace {

uint64_t build_and_measure(const std::filesystem::path& tsv,
                           const std::filesystem::path& out,
                           const std::string& bin,
                           const std::string& codec_flags) {
    run_checked(shell_quote(bin) + " build -i " + shell_quote(tsv.string()) +
                " -o " + shell_quote(out.string()) + " " + codec_flags + " --no-sketch");
    ArchiveReader ar; ar.open(out);
    return ar.archive_stats().total_compressed_bytes;
}

} // namespace

int main() {
    TempDir tmp("genopack_compguard");
    const std::string bin = GENOPACK_BIN;

    // 30 near-identical genomes: one base + copies differing at ~25 positions.
    // MEM-delta should store one panel + tiny diffs; plain stores each at ~4x.
    const std::string base = make_random_sequence(7, 200000);
    std::vector<std::pair<std::string, std::filesystem::path>> rows;
    for (int i = 0; i < 30; ++i) {
        std::string s = (i == 0) ? base : mutate_every(base, 8000 + i * 13);
        auto fa = tmp.path / ("g" + std::to_string(i) + ".fa");
        write_fasta(fa, "g" + std::to_string(i), s);
        rows.push_back({"ACC" + std::to_string(i), fa});
    }
    write_tsv(tmp.path / "in.tsv", rows);

    uint64_t plain = build_and_measure(tmp.path / "in.tsv", tmp.path / "plain.gpk",
                                       bin, "--no-mem-delta --no-dict");
    uint64_t mdelta = build_and_measure(tmp.path / "in.tsv", tmp.path / "mdelta.gpk",
                                        bin, "--mem-delta");

    double speedup = static_cast<double>(plain) / static_cast<double>(mdelta);
    std::cout << "plain=" << plain << " B, mem-delta=" << mdelta
              << " B, MEM-delta is " << speedup << "x smaller\n";

    // On near-identical genomes MEM-delta is ~5x smaller than plain in current
    // code; require at least 2x so the guard tolerates tuning noise but fails
    // loudly if codec selection stops choosing MEM-delta (speedup -> ~1x).
    require(speedup >= 2.0,
            "MEM-delta no longer beats plain by 2x on near-identical genomes "
            "(got " + std::to_string(speedup) + "x) — codec selection regressed");

    std::cout << "genopack compression guard test passed\n";
    return 0;
}
