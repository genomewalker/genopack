#include <genopack/archive.hpp>
#include <CLI/CLI.hpp>
#include <spdlog/spdlog.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

using namespace genopack;

// ── genopack build ─────────────────────────────────────────────────────────────
static int cmd_build(const std::string& input_tsv, const std::string& output_dir,
                     int threads, int zstd_level, bool no_dict, bool verbose) {
    ArchiveBuilder::Config cfg;
    cfg.io_threads           = static_cast<size_t>(threads);
    cfg.verbose              = verbose;
    cfg.shard_cfg.zstd_level = zstd_level;
    cfg.shard_cfg.train_dict = !no_dict;

    ArchiveBuilder builder(output_dir, cfg);
    builder.add_from_tsv(input_tsv);
    builder.finalize();
    return 0;
}

// ── genopack extract ───────────────────────────────────────────────────────────
static int cmd_extract(const std::string& archive_dir,
                       const std::vector<std::string>& accessions,
                       const std::string& accessions_file,
                       float min_completeness, float max_contamination,
                       const std::string& out_fasta) {
    ArchiveReader ar;
    ar.open(archive_dir);

    ExtractQuery q;
    q.min_completeness  = min_completeness;
    q.max_contamination = max_contamination;
    q.accessions = accessions;

    if (!accessions_file.empty()) {
        std::ifstream af(accessions_file);
        if (!af) { spdlog::error("Cannot open accessions file: {}", accessions_file); return 1; }
        std::string line;
        while (std::getline(af, line))
            if (!line.empty()) q.accessions.push_back(line);
    }

    std::vector<ExtractedGenome> results;
    if (!q.accessions.empty()) {
        for (const auto& acc : q.accessions) {
            auto eg = ar.fetch_by_accession(acc);
            if (eg) results.push_back(std::move(*eg));
            else spdlog::warn("Accession not found: {}", acc);
        }
    } else {
        results = ar.extract(q);
    }
    spdlog::info("Extracted {} genomes", results.size());

    std::ostream* out_stream = &std::cout;
    std::ofstream out_file;
    if (!out_fasta.empty() && out_fasta != "-") {
        out_file.open(out_fasta);
        if (!out_file) { spdlog::error("Cannot open output: {}", out_fasta); return 1; }
        out_stream = &out_file;
    }

    for (const auto& eg : results)
        *out_stream << eg.fasta;

    return 0;
}

// ── genopack stat ──────────────────────────────────────────────────────────────
static int cmd_stat(const std::string& archive_dir, bool json) {
    ArchiveReader ar;
    ar.open(archive_dir);

    auto s = ar.archive_stats();
    if (json) {
        std::cout << "{\n"
                  << "  \"generation\": "       << s.generation     << ",\n"
                  << "  \"n_shards\": "         << s.n_shards       << ",\n"
                  << "  \"n_genomes_total\": "  << s.n_genomes_total<< ",\n"
                  << "  \"n_genomes_live\": "   << s.n_genomes_live << ",\n"
                  << "  \"total_raw_bp\": "     << s.total_raw_bp   << ",\n"
                  << "  \"total_compressed_bytes\": " << s.total_compressed_bytes << ",\n"
                  << "  \"compression_ratio\": " << s.compression_ratio << "\n"
                  << "}\n";
    } else {
        std::cout << "Archive: " << archive_dir << "\n"
                  << "  Generation:       " << s.generation     << "\n"
                  << "  Shards:           " << s.n_shards       << "\n"
                  << "  Genomes (total):  " << s.n_genomes_total<< "\n"
                  << "  Genomes (live):   " << s.n_genomes_live << "\n"
                  << "  Total raw bp:     " << s.total_raw_bp   << "\n"
                  << "  Compressed bytes: " << s.total_compressed_bytes << "\n"
                  << "  Compression ratio: " << s.compression_ratio << "x\n";
    }
    return 0;
}

// ── genopack rm ────────────────────────────────────────────────────────────────
static int cmd_rm(const std::string& archive_dir,
                  const std::vector<std::string>& genome_ids) {
    ArchiveAppender app(archive_dir);
    for (const auto& id_str : genome_ids) {
        try {
            app.remove_by_accession(id_str);
        } catch (const std::exception& e) {
            spdlog::warn("rm {}: {}", id_str, e.what());
        }
    }
    app.commit();
    return 0;
}

// ── main ──────────────────────────────────────────────────────────────────────
int main(int argc, char** argv) {
    CLI::App app{"genopack — genome archive"};
    app.require_subcommand(1);

    // genopack build
    auto* build = app.add_subcommand("build", "Build a new archive from a genome TSV");
    std::string build_input, build_output;
    int build_threads = 4, build_level = 6;
    bool build_no_dict = false, build_verbose = false;
    build->add_option("-i,--input",  build_input,  "Input TSV (accession, file_path, ...)")->required();
    build->add_option("-o,--output", build_output, "Output archive directory (.gpk)")->required();
    build->add_option("-t,--threads", build_threads, "I/O threads");
    build->add_option("-z,--zstd-level", build_level, "zstd compression level (1-22)");
    build->add_flag("--no-dict", build_no_dict, "Disable shared dictionary training");
    build->add_flag("-v,--verbose", build_verbose, "Verbose progress");
    build->callback([&]() {
        std::exit(cmd_build(build_input, build_output,
                             build_threads, build_level, build_no_dict, build_verbose));
    });

    // genopack extract
    auto* extract = app.add_subcommand("extract", "Extract genomes by quality or accession");
    std::string ext_archive, ext_out, ext_acc_file;
    std::vector<std::string> ext_accessions;
    float ext_min_comp = 0, ext_max_contam = 100;
    extract->add_option("archive", ext_archive, "Archive directory")->required();
    extract->add_option("--accession", ext_accessions, "Accession to extract (repeatable)");
    extract->add_option("--accessions-file", ext_acc_file, "File with one accession per line");
    extract->add_option("--min-completeness",  ext_min_comp,   "Minimum completeness %");
    extract->add_option("--max-contamination", ext_max_contam, "Maximum contamination %");
    extract->add_option("-o,--out", ext_out, "Output FASTA (default: stdout)");
    extract->callback([&]() {
        std::exit(cmd_extract(ext_archive, ext_accessions, ext_acc_file,
                              ext_min_comp, ext_max_contam, ext_out));
    });

    // genopack stat
    auto* stat = app.add_subcommand("stat", "Show archive statistics");
    std::string stat_archive;
    bool stat_json = false;
    stat->add_option("archive", stat_archive, "Archive directory")->required();
    stat->add_flag("--json", stat_json, "Output JSON");
    stat->callback([&]() {
        std::exit(cmd_stat(stat_archive, stat_json));
    });

    // genopack add
    auto* add_cmd = app.add_subcommand("add", "Append genomes to an existing archive");
    std::string add_archive, add_input;
    add_cmd->add_option("archive", add_archive, "Archive directory")->required();
    add_cmd->add_option("-i,--input", add_input, "Input TSV")->required();
    add_cmd->callback([&]() {
        ArchiveAppender app2(add_archive);
        app2.add_from_tsv(add_input);
        app2.commit();
        std::exit(0);
    });

    // genopack rm
    auto* rm_cmd = app.add_subcommand("rm", "Tombstone (delete) genomes by accession");
    std::string rm_archive;
    std::vector<std::string> rm_ids;
    rm_cmd->add_option("archive", rm_archive, "Archive directory")->required();
    rm_cmd->add_option("genome-ids", rm_ids, "Accessions to remove")->required();
    rm_cmd->callback([&]() {
        std::exit(cmd_rm(rm_archive, rm_ids));
    });

    CLI11_PARSE(app, argc, argv);
    return 0;
}
