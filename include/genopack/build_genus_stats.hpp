#pragma once
#include "gcov.hpp"
#include "fmhr.hpp"
#include "util.hpp"
#include <array>
#include <cstring>
#include <string>
#include <vector>

namespace genopack {

class ArchiveReader;

struct ProfileEntry {
    std::array<float, 136> p;
    float    freq_norm;
    uint32_t bp;
    std::array<float, 16> rho;
};

struct GenusAccum {
    std::vector<ProfileEntry> flat;            // all contigs, all genomes, contiguous
    std::vector<uint32_t>     genome_offsets;  // G+1 entries (CSR); empty until first add_genome

    void add_genome(const std::vector<ContigProfile>& cps) {
        if (cps.empty()) return;
        if (genome_offsets.empty()) genome_offsets.push_back(0);
        for (const auto& cp : cps) {
            std::array<float, 136> a; std::memcpy(a.data(), cp.p,  136 * sizeof(float));
            std::array<float,  16> r; std::memcpy(r.data(), cp.rho,  16 * sizeof(float));
            flat.push_back({a, cp.freq_norm, cp.bp, r});
        }
        genome_offsets.push_back(static_cast<uint32_t>(flat.size()));
    }
    uint32_t n_genomes_with_long() const {
        return genome_offsets.empty() ? 0u : static_cast<uint32_t>(genome_offsets.size() - 1);
    }
    uint32_t n_long_contigs() const { return static_cast<uint32_t>(flat.size()); }
};

// Finalize accumulator → GcovEntry; adds to writer and returns the entry for immediate scoring.
GcovEntry finalize_and_add_genus(const std::string& key, uint32_t n_members,
                                  GenusAccum& acc, GcovWriter& writer);

GcovWriter build_genus_covariance(const ArchiveReader& ar, int threads = 1);
GcovWriter build_family_covariance(const ArchiveReader& ar, int threads = 1);
FmhrWriter build_fmh_genus_section(const ArchiveReader& ar, int threads = 1);

// Fused single-pass build: GCOV + FCOV + FMHR in one shard scan.
// Decode each genome FASTA once; dispatch TNF profiles to genus+family accumulators
// and FMH hashes to genus hash sets simultaneously.
struct GcovFcovFmhrResult {
    GcovWriter gcov;
    GcovWriter fcov;
    FmhrWriter fmhr;
};
GcovFcovFmhrResult build_gcov_fcov_fmhr(const ArchiveReader& ar, int threads = 1);

} // namespace genopack
