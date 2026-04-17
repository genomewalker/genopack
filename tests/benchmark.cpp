#include <genopack/archive.hpp>
#include <genopack/cidx.hpp>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <vector>
#include <unistd.h>
#include "test_utils.hpp"

using namespace genopack;

struct BenchDir {
    std::filesystem::path path;

    explicit BenchDir(const std::filesystem::path& base_dir) {
        path = base_dir / (".genopack_bench_" + std::to_string(::getpid()));
        std::error_code ec;
        std::filesystem::remove_all(path, ec);
        std::filesystem::create_directories(path);
    }

    ~BenchDir() {
        std::error_code ec;
        std::filesystem::remove_all(path, ec);
    }
};

static std::string first_accession_from_tsv(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open benchmark tsv");
    std::string line;
    std::getline(in, line); // header
    if (!std::getline(in, line))
        throw std::runtime_error("benchmark tsv is empty");
    size_t tab = line.find('\t');
    if (tab == std::string::npos)
        throw std::runtime_error("benchmark tsv missing accession column");
    return line.substr(0, tab);
}

static std::vector<std::string> sample_accessions_from_tsv(const std::filesystem::path& path,
                                                           size_t max_count) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open benchmark tsv");
    std::string line;
    std::getline(in, line); // header
    std::vector<std::string> out;
    while (out.size() < max_count && std::getline(in, line)) {
        size_t tab = line.find('\t');
        if (tab == std::string::npos) continue;
        out.push_back(line.substr(0, tab));
    }
    return out;
}

int main(int argc, char** argv) {
    using namespace genopack_test;
    std::filesystem::path bench_base =
        std::filesystem::current_path() / "build";
    if (const char* env_dir = std::getenv("GENOPACK_BENCH_DIR"); env_dir && *env_dir)
        bench_base = env_dir;
    std::filesystem::create_directories(bench_base);
    BenchDir tmp(bench_base);

    std::filesystem::path input_tsv;
    std::string probe_accession;
    std::vector<std::string> probe_accessions;
    std::vector<std::string> selected_modes;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--tsv" && i + 1 < argc) {
            input_tsv = argv[++i];
        } else if (arg == "--accession" && i + 1 < argc) {
            probe_accession = argv[++i];
        } else if (arg == "--modes" && i + 1 < argc) {
            std::string csv = argv[++i];
            size_t start = 0;
            while (start <= csv.size()) {
                size_t comma = csv.find(',', start);
                std::string mode = csv.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
                if (!mode.empty())
                    selected_modes.push_back(mode);
                if (comma == std::string::npos)
                    break;
                start = comma + 1;
            }
        } else {
            throw std::runtime_error("usage: genopack_benchmark [--tsv path] [--accession acc] [--modes plain,delta]");
        }
    }

    if (input_tsv.empty()) {
        const size_t n_genomes = 24;
        const size_t genome_len = 160000;
        const std::string seed = make_sequence("ACGTTGCAAGTCCTGA", genome_len);

        std::vector<std::pair<std::string, std::filesystem::path>> rows;
        rows.reserve(n_genomes);
        for (size_t i = 0; i < n_genomes; ++i) {
            std::string seq = seed;
            size_t step = 89 + i * 7;
            for (size_t j = step; j < seq.size(); j += step)
                seq[j] = "ACGT"[(j / step + i) & 3];
            auto fasta = tmp.path / ("g" + std::to_string(i) + ".fa");
            write_fasta(fasta, "g" + std::to_string(i), seq);
            rows.emplace_back("ACC_" + std::to_string(i), fasta);
        }
        input_tsv = tmp.path / "input.tsv";
        write_tsv(input_tsv, rows);
        for (size_t i = 0; i < std::min<size_t>(8, rows.size()); ++i)
            probe_accessions.push_back(rows[i].first);
        if (probe_accession.empty())
            probe_accession = "ACC_7";
    } else if (probe_accession.empty()) {
        probe_accession = first_accession_from_tsv(input_tsv);
    }
    if (probe_accessions.empty())
        probe_accessions = sample_accessions_from_tsv(input_tsv, 8);
    if (probe_accessions.empty())
        probe_accessions.push_back(probe_accession);

    struct Mode {
        std::string name;
        std::string flags;
    };
    const std::vector<Mode> modes = {
        {"auto", "--no-hnsw"},
        {"plain", "--no-dict --no-hnsw"},
        {"delta", "--delta --no-hnsw"},
        {"mem_delta", "--mem-delta --no-hnsw"},
    };
    std::vector<Mode> active_modes;
    if (selected_modes.empty()) {
        active_modes = modes;
    } else {
        for (const auto& wanted : selected_modes) {
            auto it = std::find_if(modes.begin(), modes.end(), [&](const Mode& mode) {
                return mode.name == wanted;
            });
            if (it == modes.end())
                throw std::runtime_error("unknown benchmark mode: " + wanted);
            active_modes.push_back(*it);
        }
    }

    const std::string bin = GENOPACK_BIN;

    std::cout << "mode,build_sec,open_usec,file_bytes,ratio,cold_slice_usec,warm_slice_usec,fetch_usec\n";
    for (const auto& mode : active_modes) {
        auto out = tmp.path / (mode.name + ".gpk");

        auto t0 = std::chrono::steady_clock::now();
        run_checked(shell_quote(bin) + " build -i " + shell_quote(input_tsv.string()) +
                    " -o " + shell_quote(out.string()) + " " + mode.flags);
        auto t1 = std::chrono::steady_clock::now();

        auto o0 = std::chrono::steady_clock::now();
        ArchiveReader ar;
        ar.open(out);
        auto o1 = std::chrono::steady_clock::now();
        auto stats = ar.archive_stats();

        constexpr int cold_reps = 8;
        auto c0 = std::chrono::steady_clock::now();
        for (int i = 0; i < cold_reps; ++i) {
            ArchiveReader cold;
            cold.open(out);
            const std::string& acc = probe_accessions[static_cast<size_t>(i) % probe_accessions.size()];
            uint64_t start = 8000 + static_cast<uint64_t>((i * 7919) % 50000);
            auto seq = cold.fetch_sequence_slice_by_accession(acc, start, 4096);
            require(seq.has_value() && seq->size() == 4096, "benchmark cold slice failed");
        }
        auto c1 = std::chrono::steady_clock::now();

        constexpr int warm_reps = 32;
        auto s0 = std::chrono::steady_clock::now();
        for (int i = 0; i < warm_reps; ++i) {
            const std::string& acc = probe_accessions[static_cast<size_t>(i) % probe_accessions.size()];
            uint64_t start = 10000 + static_cast<uint64_t>((i * 7919) % 75000);
            auto seq = ar.fetch_sequence_slice_by_accession(acc, start, 4096);
            require(seq.has_value() && seq->size() == 4096, "benchmark warm slice failed");
        }
        auto s1 = std::chrono::steady_clock::now();

        constexpr int fetch_reps = 8;
        auto f0 = std::chrono::steady_clock::now();
        for (int i = 0; i < fetch_reps; ++i) {
            const std::string& acc = probe_accessions[static_cast<size_t>(i) % probe_accessions.size()];
            auto genome = ar.fetch_by_accession(acc);
            require(genome.has_value() && !genome->fasta.empty(), "benchmark fetch failed");
        }
        auto f1 = std::chrono::steady_clock::now();

        double build_sec =
            std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
        auto open_us =
            std::chrono::duration_cast<std::chrono::microseconds>(o1 - o0).count();
        auto cold_slice_us =
            std::chrono::duration_cast<std::chrono::microseconds>(c1 - c0).count() / cold_reps;
        auto warm_slice_us =
            std::chrono::duration_cast<std::chrono::microseconds>(s1 - s0).count() / warm_reps;
        auto fetch_us =
            std::chrono::duration_cast<std::chrono::microseconds>(f1 - f0).count() / fetch_reps;

        std::cout << mode.name << ","
                  << std::fixed << std::setprecision(4) << build_sec << ","
                  << open_us << ","
                  << std::filesystem::file_size(out) << ","
                  << std::setprecision(4) << stats.compression_ratio << ","
                  << cold_slice_us << ","
                  << warm_slice_us << ","
                  << fetch_us << "\n";
    }

    // ── CIDX benchmark ────────────────────────────────────────────────────────
    // Build one plain archive, then benchmark contig lookup at various scales.
    {
        auto cidx_out = tmp.path / "cidx_bench.gpk";
        run_checked(shell_quote(bin) + " build -i " + shell_quote(input_tsv.string()) +
                    " -o " + shell_quote(cidx_out.string()) + " --no-hnsw");

        ArchiveReader ar;
        ar.open(cidx_out);

        // Collect contig accession strings by fetching all genomes and parsing headers
        std::vector<std::string> contig_accs;
        for (const auto& genome_acc : probe_accessions) {
            auto genome = ar.fetch_by_accession(genome_acc);
            if (!genome) continue;
            parse_fasta_contig_accessions(genome->fasta,
                [&](std::string_view acc) { contig_accs.emplace_back(acc); });
        }
        // Also scan all genomes in archive for a larger pool
        ar.scan_shards([&](const uint8_t*, uint64_t, uint64_t, uint32_t shard_id) {
            if (contig_accs.size() >= 200000) return;
            (void)shard_id;
        });
        // Fallback: if FASTA parsing gave nothing, synthesise fake accession strings
        if (contig_accs.empty()) {
            for (size_t i = 0; i < 10000; ++i)
                contig_accs.push_back("ACC_" + std::to_string(i));
        }

        // Build query sets at several sizes, both scattered (random) and clustered (first N)
        const std::vector<size_t> batch_sizes = {100, 1000, 10000,
                                                  std::min<size_t>(100000, contig_accs.size())};
        const std::vector<size_t> thread_counts = {1, 4,
                                                    std::min<size_t>(8u, std::thread::hardware_concurrency())};

        std::mt19937_64 rng(42);

        // Header
        std::cout << "\ncidx_mode,n_queries,n_threads,pattern,total_us,per_query_ns,found\n";

        // Single lookup baseline
        {
            const int single_reps = 1000;
            std::vector<std::string_view> sv_pool;
            for (const auto& s : contig_accs) sv_pool.push_back(s);

            auto s0 = std::chrono::steady_clock::now();
            size_t found = 0;
            for (int i = 0; i < single_reps; ++i) {
                const auto& acc = contig_accs[static_cast<size_t>(i) % contig_accs.size()];
                if (ar.find_contig_genome_id(acc) != UINT32_MAX) ++found;
            }
            auto s1 = std::chrono::steady_clock::now();
            auto total_us = std::chrono::duration_cast<std::chrono::microseconds>(s1 - s0).count();
            std::cout << "cidx_single," << single_reps << ",1,scattered,"
                      << total_us << ","
                      << (total_us * 1000 / single_reps) << ","
                      << found << "\n";
        }

        // batch_find at various sizes × thread counts × patterns
        for (size_t n : batch_sizes) {
            // Build scattered query (random sample from pool)
            std::vector<std::string> scattered(n);
            std::uniform_int_distribution<size_t> dist(0, contig_accs.size() - 1);
            for (size_t i = 0; i < n; ++i)
                scattered[i] = contig_accs[dist(rng)];

            // Clustered query: first N contigs (all from a narrow hash range)
            std::vector<std::string> clustered(n);
            for (size_t i = 0; i < n; ++i)
                clustered[i] = contig_accs[i % contig_accs.size()];

            for (size_t threads : thread_counts) {
                for (const char* pattern : {"scattered", "clustered"}) {
                    const auto& queries = (pattern[0] == 's') ? scattered : clustered;
                    std::vector<std::string_view> sv(n);
                    for (size_t i = 0; i < n; ++i) sv[i] = queries[i];

                    std::vector<uint32_t> out(n);
                    auto b0 = std::chrono::steady_clock::now();
                    ar.batch_find_contig_genome_ids(sv.data(), out.data(), n, threads);
                    auto b1 = std::chrono::steady_clock::now();

                    size_t found = 0;
                    for (uint32_t gid : out) if (gid != UINT32_MAX) ++found;

                    auto total_us = std::chrono::duration_cast<std::chrono::microseconds>(b1 - b0).count();
                    std::cout << "cidx_batch," << n << "," << threads << "," << pattern << ","
                              << total_us << ","
                              << (n > 0 ? total_us * 1000 / n : 0) << ","
                              << found << "\n";
                }
            }
        }
    }

    return 0;
}
