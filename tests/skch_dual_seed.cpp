// Round-trip test for SKCH V4 (dual-seed OPH):
//
//   1. For 10 synthetic genomes, compute BOTH a reference sig1 (seed=42) and
//      a reference sig2 (seed=43) by calling the single-seed
//      sketch_oph_from_buffer twice (ground truth).
//   2. Build an in-memory SKC4 archive via SkchWriter.
//   3. Re-read via SkchReader::sketch_for_ids and assert that sig1/sig2 are
//      byte-exact against the reference sketches and that mask matches too.
//   4. Synthesise an SKC3 magic payload and assert SkchReader::peek_params
//      throws the V4-required error message.

#include <genopack/format.hpp>
#include <genopack/oph_sketch.hpp>
#include <genopack/skch.hpp>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include "test_utils.hpp"

using namespace genopack;
using namespace genopack_test;

static std::string random_sequence(std::mt19937& rng, size_t len) {
    static const char ABC[4] = {'A','C','G','T'};
    std::string out;
    out.reserve(len);
    std::uniform_int_distribution<int> d(0, 3);
    for (size_t i = 0; i < len; ++i) out.push_back(ABC[d(rng)]);
    return out;
}

int main() {
    constexpr uint32_t kK           = 16;
    constexpr uint32_t kSketchSize  = 1024;
    constexpr int      kN           = 10;
    constexpr uint64_t kSeed1       = 42;
    constexpr uint64_t kSeed2       = 43;

    // 1. Generate genomes + compute reference dual-seed sketches.
    std::mt19937 rng(123456);
    std::vector<std::string>        genomes(kN);
    std::vector<OPHSketchResult>    ref_s1(kN);
    std::vector<OPHSketchResult>    ref_s2(kN);
    for (int i = 0; i < kN; ++i) {
        genomes[i] = random_sequence(rng, 150000 + 1000 * i);
        // Wrap in a FASTA header so the sketch buffer parser accepts it.
        std::string fasta = ">g" + std::to_string(i) + "\n" + genomes[i] + "\n";
        genomes[i] = std::move(fasta);

        OPHSketchConfig sc;
        sc.kmer_size   = kK;
        sc.sketch_size = kSketchSize;
        sc.syncmer_s   = 0;
        sc.seed        = kSeed1;
        ref_s1[i] = sketch_oph_from_buffer(genomes[i].data(), genomes[i].size(), sc);
        sc.seed   = kSeed2;
        ref_s2[i] = sketch_oph_from_buffer(genomes[i].data(), genomes[i].size(), sc);
    }

    // 2. Write a SkchWriter archive to a temp file.
    TempDir tmp("skch_dual_seed");
    const auto gpk_path = tmp.path / "skch.bin";

    SectionDesc sd;
    {
        SkchWriter w(kSketchSize, kK, /*syncmer_s=*/0, kSeed1, kSeed2, tmp.path.string());
        for (int i = 0; i < kN; ++i) {
            const auto& r1 = ref_s1[i];
            const auto& r2 = ref_s2[i];
            std::vector<uint16_t> sig1(r1.signature.size());
            std::vector<uint16_t> sig2(r2.signature.size());
            for (size_t j = 0; j < r1.signature.size(); ++j) {
                sig1[j] = static_cast<uint16_t>(r1.signature[j] >> 16);
                sig2[j] = static_cast<uint16_t>(r2.signature[j] >> 16);
            }
            // Masks are seed-independent; both refs should agree bit-for-bit.
            require(r1.real_bins_bitmask == r2.real_bins_bitmask,
                    "reference masks differ between seeds — unexpected");

            w.add(static_cast<GenomeId>(i + 1),
                  sig1, sig2,
                  r1.n_real_bins,
                  r1.genome_length,
                  r1.real_bins_bitmask);
        }

        AppendWriter aw;
        aw.create(gpk_path);
        sd = w.finalize(aw, /*section_id=*/1);
        aw.flush();
        aw.close();
    }

    // 3. Read back and compare.
    {
        MmapFileReader mm;
        mm.open(gpk_path);

        SkchReader r;
        r.open(mm.data(), sd.file_offset, sd.compressed_size);
        require(r.version()    == 4, "SkchReader should report V4");
        require(r.n_genomes()  == static_cast<uint32_t>(kN), "n_genomes mismatch");
        require(r.sketch_size()== kSketchSize, "sketch_size mismatch");
        require(r.kmer_size()  == kK, "kmer_size mismatch");
        require(r.seed1()      == kSeed1, "seed1 mismatch");
        require(r.seed2()      == kSeed2, "seed2 mismatch");
        require(r.has_sig2(),             "V4 must report has_sig2()");

        std::vector<GenomeId> ids;
        for (int i = 0; i < kN; ++i) ids.push_back(static_cast<GenomeId>(i + 1));
        std::sort(ids.begin(), ids.end());

        std::vector<int> seen(kN, 0);
        r.sketch_for_ids(ids, kK, kSketchSize,
            [&](size_t row, const SketchResult& sr) {
                GenomeId gid = ids[row];
                int i = static_cast<int>(gid) - 1;
                seen[i]++;
                require(sr.sketch_size == kSketchSize, "cb sketch_size mismatch");
                require(sr.kmer_size   == kK,           "cb kmer_size mismatch");

                const auto& r1 = ref_s1[i];
                const auto& r2 = ref_s2[i];
                for (uint32_t j = 0; j < kSketchSize; ++j) {
                    uint16_t ref1 = static_cast<uint16_t>(r1.signature[j] >> 16);
                    uint16_t ref2 = static_cast<uint16_t>(r2.signature[j] >> 16);
                    require(sr.sig[j]  == ref1, "sig1 mismatch");
                    require(sr.sig2[j] == ref2, "sig2 mismatch");
                }
                require(sr.n_real_bins  == r1.n_real_bins,  "n_real_bins mismatch");
                require(sr.genome_length== r1.genome_length,"genome_length mismatch");
                require(sr.mask_words   == static_cast<uint32_t>(r1.real_bins_bitmask.size()),
                        "mask_words mismatch");
                for (uint32_t w = 0; w < sr.mask_words; ++w)
                    require(sr.mask[w] == r1.real_bins_bitmask[w], "mask word mismatch");
            });

        for (int i = 0; i < kN; ++i)
            require(seen[i] == 1, "genome row not visited exactly once");
    }

    // 4. Synthesise an SKC3 magic payload and assert rejection.
    {
        SkchSeekHdr fake{};
        fake.magic       = SKCH_V3_MAGIC;
        fake.n_frames    = 0;
        fake.frame_size  = 16384;
        fake.n_genomes   = 0;
        fake.sketch_size = kSketchSize;
        fake.n_kmer_sizes= 1;
        fake.kmer_sizes[0]= kK;
        fake.syncmer_s   = 0;
        fake.mask_words  = (kSketchSize + 63) / 64;
        fake.seed1       = kSeed1;
        fake.seed2       = kSeed2;

        auto fake_path = tmp.path / "fake_v3.bin";
        FILE* fp = std::fopen(fake_path.c_str(), "wb");
        require(fp != nullptr, "cannot open fake_v3.bin");
        std::fwrite(&fake, sizeof(fake), 1, fp);
        std::fclose(fp);

        MmapFileReader mm;
        mm.open(fake_path);

        bool threw = false;
        std::string msg;
        try {
            SkchReader::peek_params(mm.data(), 0, mm.size());
        } catch (const std::exception& ex) {
            threw = true;
            msg   = ex.what();
        }
        require(threw, "peek_params did not throw on SKC3 payload");
        require(msg.find("V4 required") != std::string::npos,
                std::string("unexpected rejection message: ") + msg);
    }

    std::cout << "genopack skch dual-seed test passed\n";
    return 0;
}
