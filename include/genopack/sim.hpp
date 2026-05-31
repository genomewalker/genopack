#pragma once
#include <cstdint>
#include <filesystem>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace genopack {

struct SimRef {
    std::filesystem::path fasta;
    std::string           taxonomy;
};

struct SimConfig {
    std::vector<SimRef>       refs;
    std::filesystem::path     contam;
    std::vector<double>       completeness {0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
    std::vector<double>       contamination {0.0};
    int                       reps        = 3;
    uint64_t                  seed        = 42;
    int                       chunk_size  = 20000;
    int                       min_chunk   = 1000;
    std::filesystem::path     output_dir;
    std::filesystem::path     output_tsv;
    std::filesystem::path     manifest_tsv;
    int                       threads     = 4;
};

std::vector<std::pair<std::string, std::string>>
parse_fasta_named(const std::string& fasta_buf);

std::vector<std::pair<std::string, std::string>>
make_fragment(const std::vector<std::pair<std::string, std::string>>& contigs,
              double comp_frac, int chunk_size, int min_chunk, std::mt19937_64& rng);

double
mix_contamination(std::vector<std::pair<std::string, std::string>>& fragments,
                  const std::vector<std::pair<std::string, std::string>>& contam_chunks,
                  double cont_frac, std::mt19937_64& rng);

int run_sim(const SimConfig& cfg);

} // namespace genopack
