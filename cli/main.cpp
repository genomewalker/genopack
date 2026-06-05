#include "check/run_check.hpp"
#include "check/pass_b.hpp"
#include "run_ingest.hpp"
#include "run_report.hpp"
#include "run_fcore.hpp"
#include "run_pcore.hpp"
#include <genopack/score_bin.hpp>
#include <genopack/fmhr.hpp>
#include <set>
#include "bench/bench_grid.hpp"
#include <genopack/sim.hpp>
#include <genopack/markers_build.hpp>
#include <genopack/markers.hpp>
#include <genopack/aamer.hpp>
#include <genopack/gcov.hpp>
#include <genopack/cladesplit.hpp>
#include <fstream>
#include <limits>
#include "calibrate/calibrate.hpp"
#include <genopack/build_genus_stats.hpp>
#include <genopack/accx.hpp>
#include <genopack/archive.hpp>
#include <genopack/archive_set_reader.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <genopack/repack.hpp>
#include <genopack/catalog.hpp>
#include <genopack/cidx.hpp>
#include <genopack/checksum.hpp>
#include <genopack/section_checksum.hpp>
#include <genopack/format.hpp>
#include <genopack/gidx.hpp>
#include <genopack/qual.hpp>
#include <genopack/kmrx.hpp>
#include <genopack/merger.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/shard.hpp>
#include <genopack/toc.hpp>
#include <genopack/taxn.hpp>
#include <genopack/txdb.hpp>
#include <genopack/util.hpp>
#include <genopack/ncbi_taxdb.hpp>
#include <genopack/skch.hpp>
#include <genopack/oph_sketch.hpp>
#include <genopack/coordinator.hpp>
#include <CLI/CLI.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <atomic>
#include <fcntl.h>
#include <unistd.h>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <optional>
#include <thread>
#include <unordered_map>
#include <unordered_set>

using namespace genopack;

// Re-derive genus/quality sections on a merged archive (defined below; called by
// the parallel build after the internal merge).
static int cmd_reindex(const std::string& archive_path, bool force, bool build_txdb,
                       const std::string& cidx_tsv, int cidx_threads,
                       bool build_skch, int skch_threads,
                       int skch_kmer, int skch_size, int skch_syncmer,
                       std::vector<int> skch_kmers,
                       bool skip_gidx,
                       int repack_threads,
                       bool build_gstx,
                       bool build_qual,
                       bool build_gcov);

namespace genopack {
int cmd_stage(const std::string& input_tsv, const std::string& output,
              int threads, int block_mb);
}

// Resolve a reference data file that ships separately (markers panel, contamination
// panel). Precedence: explicit path (user override) > $GENOPACK_DATA/<name> >
// <exe_dir>/../share/genopack/<name>. Returns "" if nothing is found, so the caller
// silently skips that axis — the reference is never a required build argument.
static std::string resolve_data_file(const std::string& explicit_path, const char* name) {
    namespace fs = std::filesystem;
    if (!explicit_path.empty()) return explicit_path;   // user override wins verbatim
    std::error_code ec;
    if (const char* d = std::getenv("GENOPACK_DATA")) {
        fs::path p = fs::path(d) / name;
        if (fs::exists(p, ec)) return p.string();
    }
    fs::path exe = fs::read_symlink("/proc/self/exe", ec);
    if (!ec) {
        fs::path p = exe.parent_path() / ".." / "share" / "genopack" / name;
        if (fs::exists(p, ec)) return fs::weakly_canonical(p, ec).string();
    }
    return {};
}

// ── genopack build ─────────────────────────────────────────────────────────────
static int cmd_build(const std::string& input_tsv, const std::string& output_dir,
                     int threads, int zstd_level, bool no_dict, bool ref_dict,
                     bool delta, bool mem_delta, bool verbose, int n_parallel,
                     bool no_cidx, bool use_2bit, bool kmer_nn_sort,
                     bool taxonomy_group, const std::string& taxonomy_rank,
                     bool sketch, int sketch_kmer, int sketch_size, int sketch_syncmer,
                     std::vector<int> sketch_kmers = {},
                     bool no_gstx = false,
                     std::string markers_path = {},
                     std::string from_gpk = {},
                     uint32_t micro_genus_threshold = 0,
                     std::string from_stage = {},
                     std::string contam_panel = {}) {
    ArchiveBuilder::Config cfg;
    const bool explicit_codec = no_dict || ref_dict || delta || mem_delta;
    cfg.io_threads                        = static_cast<size_t>(std::max(1, threads));
    cfg.shard_cfg.compress_threads        = static_cast<size_t>(std::max(1, threads));
    cfg.verbose                           = verbose;
    cfg.build_cidx                        = !no_cidx;
    cfg.build_gstx                        = !no_gstx;
    cfg.micro_genus_threshold             = micro_genus_threshold;
    cfg.shard_cfg.zstd_level              = zstd_level;
    cfg.shard_cfg.auto_codec              = !explicit_codec;
    cfg.shard_cfg.train_dict              = false;
    cfg.shard_cfg.use_reference_dict      = ref_dict && !delta && !mem_delta;
    cfg.shard_cfg.use_delta               = (!explicit_codec || delta) && !mem_delta;
    cfg.shard_cfg.use_mem_delta           = mem_delta;
    cfg.shard_cfg.use_2bit_pack           = use_2bit;
    cfg.kmer_nn_sort                      = kmer_nn_sort;
    cfg.taxonomy_group                    = taxonomy_group;
    cfg.taxonomy_rank                     = taxonomy_rank;
    cfg.build_sketch                      = sketch;
    cfg.sketch_kmer_size                  = sketch_kmer;
    cfg.sketch_size                       = sketch_size;
    cfg.sketch_syncmer_s                  = sketch_syncmer;
    if (sketch_kmers.size() > 1) {
        std::sort(sketch_kmers.begin(), sketch_kmers.end());
        sketch_kmers.erase(std::unique(sketch_kmers.begin(), sketch_kmers.end()), sketch_kmers.end());
        cfg.sketch_kmer_sizes = sketch_kmers;
        cfg.sketch_kmer_size  = sketch_kmers[0];
    }
    cfg.markers_path      = resolve_data_file(markers_path, "markers.mrk");
    cfg.contam_panel_path = resolve_data_file(contam_panel, "contamination.csp");
    if (cfg.markers_path != markers_path && !cfg.markers_path.empty())
        spdlog::info("build: auto-resolved markers panel → {}", cfg.markers_path);
    if (!cfg.contam_panel_path.empty())
        spdlog::info("build: contamination panel → {}", cfg.contam_panel_path);

    if (!from_stage.empty()) {
        // Rebuild from a local .gstage cache: stream decoded sequence sequentially
        // (no per-genome NFS opens). Single process — parallelism is in the worker
        // pool; GSTX/GCOV/FMHR are computed inline over the full corpus, no merge.
        cfg.from_stage_source = from_stage;
        ArchiveBuilder builder(output_dir, cfg);
        builder.add_from_stage(from_stage);
        builder.finalize();
        return 0;
    }

    if (!from_gpk.empty()) {
        // Rebuild from an existing .gpk: stream decoded sequence from its shards
        // (sequential reads) instead of opening per-genome FASTA files on NFS.
        // Single process — parallelism is within the worker pool, no merge needed.
        cfg.from_gpk_source = from_gpk;
        // Reuse oracle: if params are unchanged, repack verbatim (no recompute).
        if (genopack::from_gpk_try_verbatim_reuse(from_gpk, output_dir, cfg))
            return 0;
        // Params changed — recompute by streaming decoded sequence from the shards.
        ArchiveBuilder builder(output_dir, cfg);
        builder.add_from_gpk(from_gpk);
        builder.finalize();
        return 0;
    }

    if (n_parallel <= 1) {
        // Single-process build
        ArchiveBuilder builder(output_dir, cfg);
        builder.add_from_tsv(input_tsv);
        builder.finalize();
        return 0;
    }

    // ── Parallel build: split TSV → N temp archives → merge ──────────────────
    spdlog::info("Parallel build: {} workers, {} io_threads each", n_parallel, cfg.io_threads);

    auto records = parse_tsv_records(input_tsv);
    if (records.empty()) { spdlog::warn("No records to build"); return 0; }

    const size_t total   = records.size();
    const size_t n_parts = static_cast<size_t>(n_parallel);

    // Determine temp dir next to the output file
    std::filesystem::path out_path = output_dir;
    if (out_path.extension() != ".gpk")
        out_path = std::filesystem::path(out_path.string() + ".gpk");
    std::filesystem::path tmp_dir = out_path.parent_path() / (out_path.stem().string() + "_tmp_parts");
    std::filesystem::create_directories(tmp_dir);

    // Split records into N contiguous parts
    std::vector<std::filesystem::path> part_paths;
    part_paths.reserve(n_parts);

    std::vector<std::future<std::exception_ptr>> futs;
    futs.reserve(n_parts);

    GenomeId gid_cursor = 1;
    for (size_t p = 0; p < n_parts; ++p) {
        size_t start = (p * total) / n_parts;
        size_t end   = ((p + 1) * total) / n_parts;
        if (start >= end) continue;

        std::filesystem::path part_path = tmp_dir / ("part_" + std::to_string(p) + ".gpk");
        part_paths.push_back(part_path);

        // Slice of records for this part
        std::vector<BuildRecord> slice(
            std::make_move_iterator(records.begin() + static_cast<ptrdiff_t>(start)),
            std::make_move_iterator(records.begin() + static_cast<ptrdiff_t>(end)));

        // Each part gets a unique genome_id range: [gid_cursor, gid_cursor + slice.size())
        ArchiveBuilder::Config part_cfg = cfg;
        part_cfg.starting_genome_id = gid_cursor;
        gid_cursor += static_cast<GenomeId>(slice.size());

        futs.push_back(std::async(std::launch::async,
            [part_cfg, part_path, slice = std::move(slice), p]() mutable -> std::exception_ptr {
                try {
                    // A sealed part (valid TailLocator at EOF) is already done.
                    if (std::filesystem::exists(part_path)) {
                        bool sealed = false;
                        try {
                            MmapFileReader mm;
                            mm.open(part_path);
                            if (mm.size() >= sizeof(TailLocator)) {
                                const auto* tail = mm.ptr_at<TailLocator>(
                                    mm.size() - sizeof(TailLocator));
                                sealed = (tail->magic == GPKT_MAGIC);
                            }
                        } catch (...) {}
                        if (sealed) {
                            spdlog::info("Part {}: skipping — sealed archive exists ({})",
                                         p, part_path.string());
                            return nullptr;
                        }
                        // A stale .gpk without a valid footer is removed. Resumable
                        // work lives in <part>.partial + .ckpt, which the builder
                        // picks up — we NEVER delete those (P2/P3).
                        spdlog::warn("Part {}: stale .gpk without footer, removing; "
                                     "resume (if any) uses .partial", p);
                        std::filesystem::remove(part_path);
                    }
                    spdlog::info("Part {}: {} genomes (gid_start={}) → {}", p, slice.size(),
                                 part_cfg.starting_genome_id, part_path.string());
                    ArchiveBuilder builder(part_path, part_cfg);
                    for (auto& r : slice) builder.add(r);
                    builder.finalize();
                    spdlog::info("Part {} done", p);
                    return nullptr;
                } catch (...) {
                    return std::current_exception();
                }
            }));
    }

    // Fault-tolerant join: collect per-part failures instead of aborting on the
    // first. Completed parts and resumable .partial/.ckpt state stay in tmp_dir so
    // a re-run resumes only the parts that need it (P28).
    size_t n_failed_parts = 0;
    std::exception_ptr first_error;
    for (auto& f : futs) {
        std::exception_ptr e = f.get();
        if (e) {
            ++n_failed_parts;
            if (!first_error) first_error = e;
        }
    }
    if (n_failed_parts > 0) {
        spdlog::error("{} of {} parts failed; keeping {} for resume on re-run",
                      n_failed_parts, part_paths.size(), tmp_dir.string());
        std::rethrow_exception(first_error);
    }
    spdlog::info("All {} parts built, merging…", part_paths.size());

    merge_archives(part_paths, output_dir, /*remap_genome_ids=*/false, cfg.build_cidx);

    // Per-part GSTX/GCOV/QUAL are partial aggregates (a genus can span parts), so
    // merge drops them. Re-derive once on the merged corpus so the -p archive is
    // complete and matches a single-process build.
    if (cfg.build_gstx) {
        spdlog::info("Re-deriving genus stats / quality on merged archive…");
        cmd_reindex(output_dir, /*force=*/false, /*build_txdb=*/false,
                    /*cidx_tsv=*/"", /*cidx_threads=*/8,
                    /*build_skch=*/false, /*skch_threads=*/8,
                    /*skch_kmer=*/-1, /*skch_size=*/-1, /*skch_syncmer=*/-1,
                    /*skch_kmers=*/{}, /*skip_gidx=*/true,
                    /*repack_threads=*/cfg.io_threads > 0 ? cfg.io_threads : 16,
                    /*build_gstx=*/true, /*build_qual=*/true,
                    /*build_gcov=*/cfg.build_gcov);
    }

    // Cleanup temp parts only after a fully successful merge.
    std::error_code ec;
    std::filesystem::remove_all(tmp_dir, ec);

    return 0;
}

// ── genopack merge ─────────────────────────────────────────────────────────────
static int cmd_merge(const std::vector<std::string>& inputs,
                     const std::string& list_file,
                     const std::string& output) {
    std::vector<std::filesystem::path> paths;
    for (const auto& s : inputs) paths.emplace_back(s);
    if (!list_file.empty()) {
        std::ifstream f(list_file);
        if (!f) throw std::runtime_error("Cannot open list file: " + list_file);
        std::string line;
        while (std::getline(f, line)) {
            if (!line.empty()) paths.emplace_back(line);
        }
    }
    if (paths.size() < 2)
        throw std::runtime_error("merge: need at least 2 input archives");
    merge_archives(paths, output);
    return 0;
}

// ── genopack extract ───────────────────────────────────────────────────────────
static int cmd_extract(const std::string& archive_dir,
                       const std::vector<std::string>& accessions,
                       const std::string& accessions_file,
                       float min_completeness, float max_contamination,
                       const std::string& out_fasta,
                       const std::string& out_dir) {
    ArchiveSetReader ar = open_archive_auto(archive_dir);

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

    // --output-dir: write one {accession}.fa per genome, open archive once
    if (!out_dir.empty()) {
        namespace fs = std::filesystem;
        fs::create_directories(out_dir);
        if (q.accessions.empty()) {
            spdlog::error("--output-dir requires --accession or --accessions-file");
            return 1;
        }
        size_t written = 0;
        ar.visit_shard_batches(q.accessions, [&](ArchiveSetReader::ShardBatch& batch) {
            for (auto& [idx, eg] : batch) {
                if (eg.fasta.empty()) {
                    spdlog::warn("Accession not found or deleted: {}", q.accessions[idx]);
                    continue;
                }
                auto p = fs::path(out_dir) / (q.accessions[idx] + ".fa");
                std::ofstream of(p);
                if (!of) { spdlog::error("Cannot write: {}", p.string()); return; }
                of << eg.fasta;
                ++written;
            }
        });
        spdlog::info("Extracted {} genomes to {}", written, out_dir);
        return 0;
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

// ── genopack slice ─────────────────────────────────────────────────────────────
static int cmd_slice(const std::string& archive_dir,
                     const std::string& accession,
                     uint64_t start,
                     uint64_t length,
                     bool fasta_header) {
    ArchiveSetReader ar = open_archive_auto(archive_dir);

    auto seq = ar.fetch_sequence_slice_by_accession(accession, start, length);
    if (!seq) {
        spdlog::error("Accession not found: {}", accession);
        return 1;
    }

    if (fasta_header)
        std::cout << ">" << accession << ":" << start << "+" << length << "\n";
    std::cout << *seq << "\n";
    return 0;
}

// ── genopack stat ──────────────────────────────────────────────────────────────
static int cmd_stat(const std::string& archive_dir, bool json) {
    ArchiveSetReader ar = open_archive_auto(archive_dir);

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

// ── genopack inspect ───────────────────────────────────────────────────────────
// Per-archive or per-directory sketch-layout report + preload-cost estimate.
// Accepts either a single .gpk archive or a directory containing one or more .gpk parts.
static int cmd_inspect(const std::string& path, bool json) {
    namespace fs = std::filesystem;
    std::vector<fs::path> archives;
    fs::path p = path;
    if (fs::is_directory(p)) {
        bool self_is_archive = fs::exists(p / "toc.bin") || p.extension() == ".gpk";
        if (self_is_archive) {
            archives.push_back(p);
        } else {
            for (auto& e : fs::directory_iterator(p)) {
                if (e.path().extension() != ".gpk") continue;
                if (e.is_directory() || e.is_regular_file())
                    archives.push_back(e.path());
            }
            std::sort(archives.begin(), archives.end());
        }
    } else {
        archives.push_back(p);
    }
    if (archives.empty()) {
        spdlog::error("No archive found at {}", path);
        return 1;
    }

    struct Row {
        std::string name;
        size_t n_genomes = 0;
        uint32_t sketch_size = 0;
        std::vector<uint32_t> kmer_sizes;
        uint32_t mask_words = 0;
        size_t bytes_per_sketch_per_k = 0;
        size_t bytes_per_genome = 0;  // summed over all k
        size_t preload_bytes = 0;
    };
    std::vector<Row> rows;
    rows.reserve(archives.size());

    size_t total_genomes = 0;
    size_t total_preload_bytes = 0;

    for (const auto& ap : archives) {
        Row r;
        r.name = ap.filename().string();
        ArchiveReader ar;
        ar.open(ap.string());
        auto stats = ar.archive_stats();
        r.n_genomes = stats.n_genomes_live;
        r.sketch_size = ar.sketch_sketch_size();
        r.kmer_sizes = ar.available_sketch_kmer_sizes();
        if (r.sketch_size > 0) {
            r.mask_words = (r.sketch_size + 63u) / 64u;
            // per k: sig1 + sig2 (uint16) + mask (uint64) + n_real_bins (uint32) + genome_length (uint64, shared once per genome but charge once)
            r.bytes_per_sketch_per_k = 2ull * 2ull * r.sketch_size   // sig1+sig2
                                     + 8ull * r.mask_words            // mask
                                     + 4ull;                          // n_real_bins (u32) per k
            const size_t n_k = r.kmer_sizes.empty() ? 1 : r.kmer_sizes.size();
            r.bytes_per_genome = r.bytes_per_sketch_per_k * n_k + 8ull; // + genome_length once
            r.preload_bytes = r.bytes_per_genome * r.n_genomes;
        }
        total_genomes += r.n_genomes;
        total_preload_bytes += r.preload_bytes;
        rows.push_back(std::move(r));
    }

    auto fmt_bytes = [](size_t b) {
        char buf[64];
        if (b >= (1ull << 40)) std::snprintf(buf, sizeof(buf), "%.2f TiB", b / double(1ull << 40));
        else if (b >= (1ull << 30)) std::snprintf(buf, sizeof(buf), "%.2f GiB", b / double(1ull << 30));
        else if (b >= (1ull << 20)) std::snprintf(buf, sizeof(buf), "%.2f MiB", b / double(1ull << 20));
        else std::snprintf(buf, sizeof(buf), "%zu B", b);
        return std::string(buf);
    };

    if (json) {
        std::cout << "{\n  \"archives\": [\n";
        for (size_t i = 0; i < rows.size(); ++i) {
            const auto& r = rows[i];
            std::cout << "    {\"name\":\"" << r.name << "\",\"n_genomes\":" << r.n_genomes
                      << ",\"sketch_size\":" << r.sketch_size
                      << ",\"mask_words\":" << r.mask_words
                      << ",\"kmer_sizes\":[";
            for (size_t j = 0; j < r.kmer_sizes.size(); ++j) {
                if (j) std::cout << ",";
                std::cout << r.kmer_sizes[j];
            }
            std::cout << "],\"bytes_per_genome\":" << r.bytes_per_genome
                      << ",\"preload_bytes\":" << r.preload_bytes << "}"
                      << (i + 1 < rows.size() ? "," : "") << "\n";
        }
        std::cout << "  ],\n  \"total_genomes\":" << total_genomes
                  << ",\n  \"total_preload_bytes\":" << total_preload_bytes << "\n}\n";
    } else {
        std::cout << "Inspect: " << path << "\n";
        for (const auto& r : rows) {
            std::cout << "  " << r.name << "\n"
                      << "    genomes (live):     " << r.n_genomes << "\n"
                      << "    sketch_size (bins): " << r.sketch_size << "\n"
                      << "    mask_words:         " << r.mask_words << "\n"
                      << "    kmer_sizes:         [";
            for (size_t j = 0; j < r.kmer_sizes.size(); ++j) {
                if (j) std::cout << ", ";
                std::cout << r.kmer_sizes[j];
            }
            std::cout << "]\n"
                      << "    bytes/sketch/k:     " << r.bytes_per_sketch_per_k << "\n"
                      << "    bytes/genome:       " << r.bytes_per_genome << "\n"
                      << "    preload size:       " << fmt_bytes(r.preload_bytes) << "\n";
        }
        if (rows.size() > 1) {
            std::cout << "  ─────────────────────────────\n"
                      << "  TOTAL genomes:        " << total_genomes << "\n"
                      << "  TOTAL preload size:   " << fmt_bytes(total_preload_bytes) << "\n";
        }
    }
    return 0;
}

// ── genopack dedup ─────────────────────────────────────────────────────────────
static int cmd_dedup(const std::string& archive_path, bool dry_run) {
    ArchiveReader reader;
    reader.open(archive_path);

    ExtractQuery q;
    auto all = reader.filter_meta(q);

    all.erase(std::remove_if(all.begin(), all.end(),
        [](const GenomeMeta& m){ return m.is_deleted(); }), all.end());

    std::sort(all.begin(), all.end(), [](const GenomeMeta& a, const GenomeMeta& b) {
        return a.oph_fingerprint < b.oph_fingerprint;
    });

    std::vector<GenomeId> to_tombstone;
    size_t n_groups = 0;

    size_t i = 0;
    while (i < all.size()) {
        size_t j = i + 1;
        while (j < all.size()
               && all[j].oph_fingerprint == all[i].oph_fingerprint
               && all[j].genome_length   == all[i].genome_length) {
            ++j;
        }
        if (j > i + 1) {
            ++n_groups;
            size_t keep = i;
            for (size_t k = i + 1; k < j; ++k) {
                if (all[k].completeness_x10 > all[keep].completeness_x10 ||
                    (all[k].completeness_x10 == all[keep].completeness_x10 &&
                     all[k].genome_id < all[keep].genome_id))
                    keep = k;
            }
            for (size_t k = i; k < j; ++k)
                if (k != keep)
                    to_tombstone.push_back(all[k].genome_id);
        }
        i = j;
    }

    spdlog::info("dedup: scanned {} genomes, found {} duplicate groups, {} to remove",
                 all.size(), n_groups, to_tombstone.size());

    if (dry_run || to_tombstone.empty()) {
        if (dry_run) spdlog::info("dedup: dry-run, no changes made");
        return 0;
    }

    ArchiveAppender appender(archive_path);
    for (GenomeId id : to_tombstone)
        appender.remove(id);
    appender.commit();

    spdlog::info("dedup: tombstoned {} duplicates", to_tombstone.size());
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

// Rebuild GCOV/FCOV/FMHR from the LIVE (non-tombstoned) genomes and write the three
// sections back into the archive, replacing any existing ones. Shared by the `gcov`
// subcommand and the per-iteration ref-cleaning in `decontaminate`. build_gcov_fcov_fmhr
// is tombstone-aware, so the rebuilt consensus excludes any genomes removed so far.
// Returns 0 on success, 1 if no genera were produced.
static int rebuild_gcov_fmhr(const std::filesystem::path& gp, int threads) {
    using namespace genopack;
    ArchiveReader ar;
    ar.open(gp);
    auto [gcov_w, fcov_w, fmhr_w] = build_gcov_fcov_fmhr(ar, threads > 0 ? threads : 8);
    ar.close();
    if (gcov_w.n_genera() == 0) {
        spdlog::warn("gcov: no genera found (archive has no TAXN section?)");
        return 1;
    }

    MmapFileReader mmap;
    mmap.open(gp);
    auto toc_r = TocReader::read(mmap);
    const auto* tail    = mmap.ptr_at<TailLocator>(mmap.size() - sizeof(TailLocator));
    const uint64_t prev = tail->toc_offset;
    const uint64_t gen  = toc_r.header.generation + 1;
    mmap.close();

    AppendWriter aw;
    aw.open_append(gp);
    uint64_t next_sid = toc_r.next_section_id();
    SectionDesc gcov_sd = gcov_w.finalize(aw, next_sid++);
    SectionDesc fcov_sd = fcov_w.n_genera() > 0
        ? fcov_w.finalize(aw, next_sid, SEC_FCOV) : SectionDesc{};
    ++next_sid;
    SectionDesc fmhr_sd = fmhr_w.n_genera() > 0
        ? fmhr_w.finalize(aw, next_sid) : SectionDesc{};

    TocWriter new_toc;
    for (const auto& sd : toc_r.sections)
        if (sd.type != SEC_GCOV && sd.type != SEC_FCOV && sd.type != SEC_FMHR)
            new_toc.add_section(sd);
    new_toc.add_section(gcov_sd);
    if (fcov_w.n_genera() > 0) new_toc.add_section(fcov_sd);
    if (fmhr_w.n_genera() > 0) new_toc.add_section(fmhr_sd);
    aw.flush();
    stamp_section_checksums(gp.c_str(), new_toc.sections(), /*only_if_zero=*/true);
    new_toc.finalize(aw, gen,
                     toc_r.header.live_genome_count,
                     toc_r.header.total_genome_count,
                     prev,
                     toc_r.header.catalog_root_section_id,
                     toc_r.header.accession_root_section_id,
                     toc_r.header.tombstone_root_section_id);
    spdlog::info("gcov: {} genera, {} families, {} fmhr genera rebuilt for {}",
                 gcov_w.n_genera(), fcov_w.n_genera(), fmhr_w.n_genera(), gp.string());
    return 0;
}

// Per-genus FMH null-background contamination scorer — the genus-band signal.
// Computes fmh_minority for EVERY live genome (per-contig FracMinHash containment
// vs its own genus FMHR ref + same-family candidate genus refs), then scores the
// DELTA from the genus's own baseline as a robust z (median + MAD). Absolute
// fmh_minority thresholds fail because closely-related genera share a conserved-
// core containment background (~11.6% E.coli/Salmonella), so only the per-genus
// delta is meaningful. Sketch-based (contigs >=10kb) -> works on fragmented MAGs
// where the per-contig TNF/CCO axis is dark. Builds on the validated
// score_bin_fmh_containment + SEC_FMHR refs; the baseline is robust to a
// contaminated minority because median/MAD ignore the tail.
static void compute_fmh_minority_z(const std::filesystem::path& gpk,
                                   std::unordered_map<std::string, float>& fmh_raw,
                                   std::unordered_map<std::string, float>& fmh_z,
                                   std::unordered_map<std::string, float>& genus_median) {
    ArchiveReader ar;
    ar.open(gpk);
    std::unordered_map<std::string, GenomeId> acc_gid;
    ar.scan_genome_accessions([&](std::string_view a, GenomeId g) {
        if (!ar.is_deleted(g)) acc_gid.emplace(std::string(a), g);
    });
    auto rank = [](std::string_view tax, const char* pre) -> std::string {
        auto p = tax.find(pre);
        if (p == std::string_view::npos) return {};
        auto e = tax.find(';', p);
        std::string s(tax.substr(p, e == std::string_view::npos ? tax.size() - p : e - p));
        return s == std::string(pre) ? std::string{} : s;
    };
    std::unordered_map<std::string, std::string> acc_genus, acc_family;
    std::unordered_map<std::string, std::set<std::string>> family_genera;
    ar.scan_taxonomy([&](std::string_view a, std::string_view tax) {
        std::string acc(a);
        if (!acc_gid.count(acc)) return;
        std::string g = rank(tax, "g__"), f = rank(tax, "f__");
        if (g.empty()) return;
        acc_genus[acc] = g;
        if (!f.empty()) { acc_family[acc] = f; family_genera[f].insert(g); }
    });
    std::unordered_map<std::string, FmhrView> host;
    for (auto& [acc, g] : acc_genus) {
        if (host.count(g)) continue;
        auto v = ar.fmhr_for_genus(g);
        if (v.valid()) host[g] = v;
    }
    std::unordered_map<std::string, std::vector<FmhrView>> fam_cands;
    for (auto& [f, genera] : family_genera)
        for (auto& g : genera) { auto it = host.find(g); if (it != host.end()) fam_cands[f].push_back(it->second); }

    std::vector<std::string> accs;
    for (auto& [acc, g] : acc_genus) if (host.count(g)) accs.push_back(acc);
    std::sort(accs.begin(), accs.end());
    auto egs = ar.batch_fetch_by_accessions(accs);

    std::unordered_map<std::string, std::vector<float>> by_genus;
    for (size_t i = 0; i < accs.size(); ++i) {
        if (!egs[i]) continue;
        const std::string& acc = accs[i];
        const FmhrView& h = host[acc_genus[acc]];
        std::vector<FmhrView> refs{ h };
        auto fit = acc_family.find(acc);
        if (fit != acc_family.end())
            for (auto& v : fam_cands[fit->second])
                if (v.genus_hash != h.genus_hash) refs.push_back(v);
        if (refs.size() < 2) continue;                       // no foreign candidate genus in pack
        auto contigs = genopack::check::parse_fasta(egs[i]->fasta);
        std::vector<std::string_view> seqs;
        seqs.reserve(contigs.size());
        for (auto& c : contigs) seqs.emplace_back(c.seq);
        float v = score_bin_fmh_containment(seqs, h.genus_hash, refs, 21, 125);
        if (!std::isnan(v)) { fmh_raw[acc] = v; by_genus[acc_genus[acc]].push_back(v); }
    }
    // Per-genus robust baseline: z = (fmh - median) / (1.4826 * MAD).
    std::unordered_map<std::string, std::pair<float, float>> base;   // genus -> (median, mad)
    for (auto& [g, vals] : by_genus) {
        if (vals.size() < 3) continue;
        std::vector<float> v = vals;
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        float med = v[v.size() / 2];
        std::vector<float> dev;
        dev.reserve(v.size());
        for (float x : vals) dev.push_back(std::fabs(x - med));
        std::nth_element(dev.begin(), dev.begin() + dev.size() / 2, dev.end());
        float mad = std::max(dev[dev.size() / 2], 1e-4f);
        base[g] = { med, mad };
    }
    for (auto& [acc, v] : fmh_raw) {
        auto bit = base.find(acc_genus[acc]);
        if (bit == base.end()) continue;
        genus_median[acc] = bit->second.first;
        fmh_z[acc] = (v - bit->second.first) / (1.4826f * bit->second.second);
    }
}

// Reference-internal CCO contamination scorer — the distant-band signal. For each
// live genome, computes the per-contig TNF outlier fraction (T²∪SPE Mahalanobis vs
// the genus GCOV covariance) via score_bin. Unlike FMH-minority it needs NO foreign
// reference, so it sees family/order/phylum contamination that the same-family FMH
// candidate set is structurally blind to. Needs contigs ≥20kb (TNF), so it complements
// FMH-minority (≥10kb sketch) rather than replacing it. Robust to a contaminated
// minority: the genus covariance centroid is dominated by the clean majority, and the
// iterative rebuild de-pollutes it further each round. cco_fraction is in [0,1].
static void compute_cco_fraction(const std::filesystem::path& gpk,
                                 std::unordered_map<std::string, float>& cco,
                                 std::unordered_map<std::string, float>& cco_z) {
    using namespace genopack;
    ArchiveReader ar;
    ar.open(gpk);
    const GcovReader* gcov_rd = ar.gcov_reader();
    const GcovReader* fcov_rd = ar.fcov_reader();
    if (!gcov_rd) return;

    std::unordered_map<std::string, GenomeId> acc_gid;
    ar.scan_genome_accessions([&](std::string_view a, GenomeId g) {
        if (!ar.is_deleted(g)) acc_gid.emplace(std::string(a), g);
    });
    auto rank = [](std::string_view tax, const char* pre) -> std::string {
        auto p = tax.find(pre);
        if (p == std::string_view::npos) return {};
        auto e = tax.find(';', p);
        std::string s(tax.substr(p, e == std::string_view::npos ? tax.size() - p : e - p));
        return s == std::string(pre) ? std::string{} : s;
    };
    std::unordered_map<std::string, std::string> acc_genus, acc_family;
    ar.scan_taxonomy([&](std::string_view a, std::string_view tax) {
        std::string acc(a);
        if (!acc_gid.count(acc)) return;
        std::string g = rank(tax, "g__");
        if (g.empty()) return;
        acc_genus[acc] = g;
        std::string f = rank(tax, "f__");
        if (!f.empty()) acc_family[acc] = f;
    });

    std::vector<std::string> accs;
    for (auto& [acc, g] : acc_genus) accs.push_back(acc);
    std::sort(accs.begin(), accs.end());
    auto egs = ar.batch_fetch_by_accessions(accs);

    BinScoreConfig cfg;
    std::unordered_map<std::string, std::vector<float>> by_genus;
    for (size_t i = 0; i < accs.size(); ++i) {
        if (!egs[i]) continue;
        const std::string& acc = accs[i];
        const GcovEntry* ge = ar.gcov_for_genus(acc_genus[acc]);
        if (!ge || !(ge->flags & GCOV_FLAG_VALID)) continue;
        const GcovEntry* fe = nullptr;
        auto fit = acc_family.find(acc);
        if (fit != acc_family.end() && fcov_rd) {
            fe = ar.fcov_for_family(fit->second);
            if (fe && !(fe->flags & GCOV_FLAG_VALID)) fe = nullptr;
        }
        auto contigs = genopack::check::parse_fasta(egs[i]->fasta);
        std::vector<std::string_view> seqs;
        seqs.reserve(contigs.size());
        for (auto& c : contigs) seqs.emplace_back(c.seq);
        BinScores sc = score_bin(seqs, ge, gcov_rd, fe, fcov_rd, cfg);
        if (sc.status == BinScores::Status::OK && !std::isnan(sc.cco_fraction)) {
            cco[acc] = sc.cco_fraction;
            by_genus[acc_genus[acc]].push_back(sc.cco_fraction);
        }
    }
    // Per-genus robust baseline — the same self-calibration the FMH axis uses, and the
    // answer to "why a fixed threshold?": a clean genome in a TNF-heterogeneous genus
    // (E. coli: pangenome, plasmids, prophages) reads a high ABSOLUTE cco, but so do its
    // genus-mates, so its delta-from-median z stays low. A foreign contaminant sits far
    // above its genus-mates. z = (cco - median) / (1.4826 * MAD). Note CCO is already
    // measured against the GENUS covariance, which contains the lineage's native rRNA/
    // plasmid spread — so native mobile elements are not flagged (they are not outliers
    // vs a genus that all carries them); only foreign content is.
    std::unordered_map<std::string, std::pair<float, float>> base;
    for (auto& [g, vals] : by_genus) {
        if (vals.size() < 3) continue;
        std::vector<float> v = vals;
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        float med = v[v.size() / 2];
        std::vector<float> dev; dev.reserve(v.size());
        for (float x : vals) dev.push_back(std::fabs(x - med));
        std::nth_element(dev.begin(), dev.begin() + dev.size() / 2, dev.end());
        float mad = std::max(dev[dev.size() / 2], 1e-4f);
        base[g] = { med, mad };
    }
    for (auto& [acc, v] : cco) {
        auto bit = base.find(acc_genus[acc]);
        if (bit == base.end()) continue;
        cco_z[acc] = (v - bit->second.first) / (1.4826f * bit->second.second);
    }
}

// ── genopack decontaminate ──────────────────────────────────────────────────────
// Iteratively remove contaminated genomes so neither the DB nor its consensus
// models carry contamination. Each round rebuilds the genus/family models from the
// LIVE (non-tombstoned) genomes, so the consensus gets cleaner every pass and the
// circularity — contaminated members inflating the model that scores them — is
// broken by construction. CCO (per-contig outlier vs the genus covariance) is the
// gate: it is reference-internal, so it flags foreign contigs whether or not the
// contaminant's genus is in the DB.
static int cmd_decontaminate(const std::string& archive_path, float max_fmh_z, float min_delta,
                             float cco_abs, int max_iters, int threads, bool dry_run) {
    std::filesystem::path gpk = archive_path;
    if (!std::filesystem::exists(gpk) && gpk.extension() != ".gpk")
        gpk = std::filesystem::path(gpk.string() + ".gpk");
    if (!std::filesystem::exists(gpk)) {
        spdlog::error("decontaminate: archive not found: {}", archive_path);
        return 1;
    }

    size_t total_removed = 0;
    bool converged = false;
    for (int it = 1; it <= max_iters; ++it) {
        // Score against the CURRENT FMHR refs. Round 1 uses the as-built refs (which
        // may be polluted if contaminants were in the build); each subsequent round
        // scores against the refs rebuilt at the end of the previous iteration from
        // the cleaned live set (see rebuild_gcov_fmhr below). The per-genus median+MAD
        // baseline is robust to a contaminated minority (the tail does not move the
        // median), which is what lets round 1 detect the strongest contamination even
        // against a polluted ref. compute_fmh_minority_z scores live genomes only.

        // Layered detection — a genome is flagged if EITHER axis fires:
        //   • FMH-minority (genus band): per-contig FracMinHash containment vs same-
        //     family candidate genus refs. Sees genus-distance contamination where TNF
        //     is information-theoretically blind, but needs the contaminant's sibling
        //     genus in the pack, so it is blind to family/order/phylum distance.
        //   • CCO (distant band): per-contig TNF outlier vs the genus GCOV covariance.
        //     Reference-internal — needs no foreign ref, so it catches family/order/
        //     phylum contamination FMH cannot. Needs contigs ≥20kb (TNF).
        std::unordered_map<std::string, float> fmh_raw, fmh_z, gmed;
        compute_fmh_minority_z(gpk, fmh_raw, fmh_z, gmed);
        std::unordered_map<std::string, float> cco, cco_z;   // raw fraction + per-genus robust z
        compute_cco_fraction(gpk, cco, cco_z);

        struct Flag { std::string acc; const char* axis; float fz; float cco_pct; };
        std::vector<Flag> flagged;
        std::set<std::string> seen;
        for (auto& [acc, z] : fmh_z)
            if (z > max_fmh_z && (fmh_raw[acc] - gmed[acc]) > min_delta) {
                flagged.push_back({ acc, "fmh", z, cco.count(acc) ? cco[acc] * 100.f : NAN });
                seen.insert(acc);
            }
        // CCO gate is an ABSOLUTE threshold, not a per-genus z. Empirically the per-genus
        // baseline poisons when contamination concentrates in a genus (median itself goes
        // contaminated), and CCO's clean baseline (~0.7%% TNF-outlier bp) is a near-universal
        // property of genome compositional homogeneity, so the absolute cut is the robust,
        // self-calibration-free choice here. (cco_z is computed for diagnostics only.)
        for (auto& [acc, c] : cco) {
            if (c * 100.f > cco_abs && !seen.count(acc)) {
                flagged.push_back({ acc, "cco", fmh_z.count(acc) ? fmh_z[acc] : NAN, c * 100.f });
                seen.insert(acc);
            }
        }

        spdlog::info("decontaminate iter {}: FMH-scored {}, CCO-scored {}; {} flagged "
                     "(FMH z>{:.1f} & delta>{:.3f}  OR  CCO>{:.2f}%)",
                     it, fmh_z.size(), cco.size(), flagged.size(), max_fmh_z, min_delta, cco_abs);
        if (flagged.empty()) { converged = true; break; }
        for (auto& f : flagged)
            spdlog::info("  iter {} {} [{}]: {}  fmh_z={:.1f} cco={:.2f}%",
                         it, dry_run ? "would-remove" : "remove", f.axis, f.acc, f.fz, f.cco_pct);
        if (dry_run) {
            spdlog::info("decontaminate: dry-run — {} genome(s) would be removed; no changes made",
                         flagged.size());
            return 0;
        }
        std::vector<std::string> rm_accs;
        rm_accs.reserve(flagged.size());
        for (auto& f : flagged) rm_accs.push_back(f.acc);
        cmd_rm(gpk.string(), rm_accs);
        total_removed += rm_accs.size();

        // Iterative reference cleaning: rebuild the genus/family consensus from the
        // now-cleaner live set so the NEXT round scores against a de-polluted FMHR.
        // This is the correct order — rebuild AFTER tombstoning, so the contaminants
        // just removed are excluded (build_gcov_fcov_fmhr is tombstone-aware). It
        // breaks the circularity (contaminated members inflating the model that
        // scores them): round 1 catches the strongest contamination even against a
        // polluted ref because the per-genus median is robust to a minority tail;
        // each rebuild de-pollutes the ref so the next round reaches closer-relative
        // contamination the previous ref had masked.
        rebuild_gcov_fmhr(gpk, threads);
    }
    spdlog::info("decontaminate: {} — removed {} contaminated genome(s)",
                 converged ? "converged" : "stopped at iteration cap", total_removed);
    return 0;
}

// ── genopack taxonomy ──────────────────────────────────────────────────────────
static int cmd_taxonomy(const std::string& archive_dir, const std::string& accession,
                        bool json) {
    ArchiveSetReader ar = open_archive_auto(archive_dir);

    if (!accession.empty()) {
        auto loc = ar.taxonomy_tree_for_accession(accession);
        if (!loc) {
            spdlog::error("Accession not found or no taxonomy data: {}", accession);
            return 1;
        }
        const auto& tree = *loc->tree;
        uint64_t taxid = tree.taxid_for_accession(accession);
        if (taxid == 0) {
            spdlog::error("Accession not found: {}", accession);
            return 1;
        }
        // Print lineage from accession up to root
        std::vector<std::pair<std::string, std::string>> lineage; // (rank, name)
        uint64_t cur = taxid;
        while (cur != 0 && cur != 1) {
            auto rk  = tree.rank(cur);
            auto nm  = tree.name(cur);
            const char* rank_str = "no_rank";
            switch (rk) {
                case TaxRank::DOMAIN:  rank_str = "domain";  break;
                case TaxRank::PHYLUM:  rank_str = "phylum";  break;
                case TaxRank::CLASS:   rank_str = "class";   break;
                case TaxRank::ORDER:   rank_str = "order";   break;
                case TaxRank::FAMILY:  rank_str = "family";  break;
                case TaxRank::GENUS:   rank_str = "genus";   break;
                case TaxRank::SPECIES: rank_str = "species"; break;
                default: break;
            }
            lineage.emplace_back(rank_str, std::string(nm));
            cur = tree.parent(cur);
        }
        std::reverse(lineage.begin(), lineage.end());
        if (json) {
            std::cout << "{\"accession\":\"" << accession << "\",\"lineage\":[";
            for (size_t i = 0; i < lineage.size(); ++i) {
                if (i) std::cout << ",";
                std::cout << "{\"rank\":\"" << lineage[i].first
                          << "\",\"name\":\"" << lineage[i].second << "\"}";
            }
            std::cout << "]}\n";
        } else {
            std::cout << "Accession: " << accession << "\n";
            for (const auto& [rank, name] : lineage)
                std::cout << "  " << rank << ": " << name << "\n";
        }
        return 0;
    }

    // Archive-wide summary
    auto sum = ar.taxonomy_summary();
    if (sum.n_parts_with_taxonomy == 0) {
        spdlog::error("No taxonomy data in archive");
        return 1;
    }
    if (json) {
        std::cout << "{\n"
                  << "  \"n_nodes\": "       << sum.n_nodes_union << ",\n"
                  << "  \"n_accessions\": "  << sum.n_accessions  << ",\n"
                  << "  \"n_parts\": "       << ar.part_count()   << ",\n"
                  << "  \"n_parts_with_taxonomy\": " << sum.n_parts_with_taxonomy << "\n"
                  << "}\n";
    } else {
        std::cout << "Taxonomy summary for: " << archive_dir << "\n"
                  << "  Parts:       " << ar.part_count() << " (with taxonomy: "
                                       << sum.n_parts_with_taxonomy << ")\n"
                  << "  Nodes:       " << sum.n_nodes_union << " (union across parts)\n"
                  << "  Accessions:  " << sum.n_accessions  << "\n";
    }
    return 0;
}

// ── genopack reindex ──────────────────────────────────────────────────────────
static int cmd_reindex(const std::string& archive_path, bool force, bool build_txdb,
                       const std::string& cidx_tsv, int cidx_threads,
                       bool build_skch, int skch_threads,
                       int skch_kmer, int skch_size, int skch_syncmer,
                       std::vector<int> skch_kmers = {},
                       bool skip_gidx = false,
                       int repack_threads = 16,
                       bool build_gstx = false,
                       bool build_qual = true,
                       bool build_gcov = false) {
    // Resolve .gpk path
    std::filesystem::path gpk = archive_path;
    if (!std::filesystem::exists(gpk) && gpk.extension() != ".gpk")
        gpk = std::filesystem::path(gpk.string() + ".gpk");
    if (!std::filesystem::exists(gpk) || !std::filesystem::is_regular_file(gpk)) {
        spdlog::error("Archive not found: {}", archive_path);
        return 1;
    }

    // Read existing TOC
    MmapFileReader mmap;
    mmap.open(gpk);
    if (mmap.size() < sizeof(FileHeader)) {
        spdlog::error("File too small to be a .gpk archive");
        return 1;
    }
    auto* fh = mmap.ptr_at<FileHeader>(0);
    if (fh->magic != GPK_MAGIC) {
        spdlog::error("Not a v2 .gpk file (bad magic)");
        return 1;
    }

    Toc toc = TocReader::read(mmap);
    const auto* tail = mmap.ptr_at<TailLocator>(mmap.size() - sizeof(TailLocator));
    spdlog::info("Archive: {} (gen {}, {} sections)",
                 gpk.string(), toc.header.generation, toc.sections.size());

    bool has_gidx = !toc.find_by_type(SEC_GIDX).empty();
    bool has_txdb = !toc.find_by_type(SEC_TXDB).empty();
    bool has_cidx = !toc.find_by_type(SEC_CIDX).empty();
    bool has_skch = !toc.find_by_type(SEC_SKCH).empty();
    bool need_gidx = (!has_gidx || force) && !skip_gidx;
    bool need_txdb = build_txdb && (!has_txdb || force);
    bool need_cidx = !cidx_tsv.empty() && (!has_cidx || force);
    bool need_skch = build_skch;

    bool has_gstx_sec = !toc.find_by_type(SEC_GSTX).empty();
    bool need_gstx_check = build_gstx && (!has_gstx_sec || force);
    bool has_qual_sec = !toc.find_by_type(SEC_QUAL).empty();
    bool need_qual_check = build_qual && build_gstx && (!has_qual_sec || force);
    bool has_gcov_sec = !toc.find_by_type(SEC_GCOV).empty();
    bool need_gcov = build_gcov && (!has_gcov_sec || force);
    bool need_metabundle = (tail->meta_bundle_offset() == 0);

    if (!need_gidx && !need_txdb && !need_cidx && !need_skch && !need_gstx_check && !need_qual_check && !need_gcov && !need_metabundle) {
        spdlog::info("GIDX: already present (use --force to rebuild)");
        if (build_txdb) spdlog::info("TXDB: already present (use --force to rebuild)");
        if (!cidx_tsv.empty()) spdlog::info("CIDX: already present (use --force to rebuild)");
        if (build_gstx) spdlog::info("GSTX: already present (use --force to rebuild)");
        if (build_qual) spdlog::info("QUAL: already present (use --force to rebuild)");
        spdlog::info("Nothing to do");
        return 0;
    }
    if (need_metabundle)
        spdlog::info("MetaBundle: missing — will inject after reindex");

    // ── Build GIDX ──────────────────────────────────────────────────────────

    // Step 1: scan all CATL sections to collect genome_id -> (shard_id, catl_row_index)
    struct CatlEntry {
        uint32_t shard_id;
        uint64_t catl_row_index;
    };
    std::unordered_map<GenomeId, CatlEntry> genome_map;
    std::vector<const SectionDesc*> shrd_sections;
    GidxWriter gidx_writer;

    if (need_gidx) {
        MergedCatalogReader catalog;
        auto catl_sections = toc.find_by_type(SEC_CATL);
        for (auto* sd : catl_sections)
            catalog.add_fragment(mmap.data(), sd->file_offset, sd->compressed_size);

        uint64_t catl_row = 0;
        catalog.scan([&](const GenomeMeta& m) {
            genome_map[m.genome_id] = {m.shard_id, catl_row};
            ++catl_row;
            return true;
        });

        spdlog::info("GIDX: scanned {} catalog rows", genome_map.size());

        // Step 2: for each SHRD, scan its genome directory to find dir_index per genome_id
        shrd_sections = toc.find_by_type(SEC_SHRD);
        size_t n_indexed = 0;

        for (auto* sd : shrd_sections) {
            // Open shard section to read its directory
            ShardReader shard;
            shard.open(mmap.data(), sd->file_offset, sd->compressed_size);

            uint32_t dir_idx = 0;
            for (auto* de = shard.dir_begin(); de != shard.dir_end(); ++de, ++dir_idx) {
                auto it = genome_map.find(de->genome_id);
                if (it == genome_map.end()) continue;

                gidx_writer.add(de->genome_id,
                                static_cast<uint32_t>(sd->section_id),
                                dir_idx,
                                it->second.catl_row_index);
                ++n_indexed;
            }
        }

        spdlog::info("GIDX: indexed {} genomes across {} shards", n_indexed, shrd_sections.size());
    }

    // ── Append GIDX section + rewrite TOC ───────────────────────────────────
    // Close mmap before appending (we'll write to the same file)
    uint64_t prev_generation     = toc.header.generation;
    uint64_t live_count          = toc.header.live_genome_count;
    uint64_t total_count         = toc.header.total_genome_count;
    uint64_t prev_toc_offset     = tail->toc_offset;
    uint64_t catalog_root_id     = toc.header.catalog_root_section_id;
    uint64_t accession_root_id   = toc.header.accession_root_section_id;
    uint64_t tombstone_root_id   = toc.header.tombstone_root_section_id;
    uint64_t next_section_id     = toc.next_section_id();

    // ── Collect TAXN entries for TXDB (if requested) ────────────────────────
    TxdbWriter txdb_writer;
    if (need_txdb) {
        auto taxn_sections = toc.find_by_type(SEC_TAXN);
        if (taxn_sections.empty()) {
            spdlog::warn("TXDB: no TAXN sections found — cannot build taxonomy tree");
            need_txdb = false;
        } else {
            size_t n_taxn = 0;
            for (auto* sd : taxn_sections) {
                TaxonomyIndexReader tir;
                tir.open(mmap.data(), sd->file_offset, sd->compressed_size);
                tir.scan([&](std::string_view acc, std::string_view tax) {
                    txdb_writer.add(acc, tax);
                    ++n_taxn;
                });
            }
            spdlog::info("TXDB: loaded {} accession→lineage entries from {} TAXN sections",
                         n_taxn, taxn_sections.size());
        }
    }

    // ── Build CIDX from FASTA TSV (if requested) ────────────────────────────
    CidxWriter cidx_writer;
    if (need_cidx) {
        // Build accession → genome_id map from all ACCX sections
        std::unordered_map<std::string, GenomeId> acc_to_gid;
        for (auto* sd : toc.find_by_type(SEC_ACCX)) {
            AccessionIndexReader ar;
            ar.open(mmap.data(), sd->file_offset, sd->compressed_size);
            ar.scan([&](std::string_view acc, GenomeId gid) {
                acc_to_gid.emplace(std::string(acc), gid);
            });
        }
        spdlog::info("CIDX: {} genome accessions loaded from ACCX", acc_to_gid.size());

        // Parse TSV: skip header line, columns = accession, taxonomy, file_path
        struct TsvRow { std::string accession; std::string file_path; };
        std::vector<TsvRow> rows;
        {
            std::ifstream f(cidx_tsv);
            if (!f) {
                spdlog::error("CIDX: cannot open TSV: {}", cidx_tsv);
                return 1;
            }
            std::string line;
            bool first = true;
            while (std::getline(f, line)) {
                if (first) { first = false; continue; } // skip header
                auto t1 = line.find('\t');
                if (t1 == std::string::npos) continue;
                auto t2 = line.find('\t', t1 + 1);
                if (t2 == std::string::npos) continue;
                rows.push_back({line.substr(0, t1), line.substr(t2 + 1)});
            }
        }
        spdlog::info("CIDX: {} genomes in TSV", rows.size());

        // Parallel FASTA scan
        int n_threads = std::max(1, cidx_threads);
        std::mutex mu;
        std::atomic<size_t> done{0}, n_contigs{0}, n_missing{0};
        size_t n_rows = rows.size();

        auto worker = [&](int tid) {
            CidxWriter local;
            for (size_t i = tid; i < n_rows; i += n_threads) {
                const auto& row = rows[i];
                auto it = acc_to_gid.find(row.accession);
                if (it == acc_to_gid.end()) { n_missing++; continue; }
                GenomeId gid = it->second;
                std::string fasta;
                try { fasta = decompress_gz(row.file_path); }
                catch (...) { n_missing++; continue; }
                size_t cnt = 0;
                parse_fasta_contig_accessions(fasta, [&](std::string_view acc) {
                    local.add(acc, gid);
                    cnt++;
                });
                n_contigs += cnt;
                size_t d = ++done;
                if (d % 10000 == 0)
                    spdlog::info("CIDX: {}/{} genomes scanned, {} contigs so far", d, n_rows, n_contigs.load());
            }
            std::lock_guard<std::mutex> lk(mu);
            // merge local into cidx_writer via add_hash
            // (access entries_ directly not possible; rebuild via re-add is fine for hashes)
            // Use add() on string is unavailable for batches; finalize local and merge entries
            // Actually, merge by re-feeding: cidx_writer has no merge API.
            // Workaround: accumulate all into a shared vector under mutex.
            // Better: use add_hash with pre-computed hashes stored in local.
            // CidxWriter stores CidxEntry{acc_hash, genome_id}; expose via finalize into temp buffer?
            // Simplest: just call cidx_writer.add() under the mutex (not hot path).
            // Reset: use the fact that local writer entries are accessible...
            // Actually the cleanest solution: finalize local into a temp AppendWriter and re-read.
            // But that's complex. Instead, keep a shared vector and merge at end.
            (void)local; // handled below via per-thread vectors
        };
        (void)worker; // replaced by simpler approach below

        // Simpler: per-thread vectors of entries, merge after
        struct Entry { uint64_t hash; GenomeId gid; };
        std::vector<std::vector<Entry>> thread_entries(n_threads);

        auto worker2 = [&](int tid) {
            auto& local = thread_entries[tid];
            for (size_t i = tid; i < n_rows; i += n_threads) {
                const auto& row = rows[i];
                auto it = acc_to_gid.find(row.accession);
                if (it == acc_to_gid.end()) { n_missing++; continue; }
                GenomeId gid = it->second;
                std::string fasta;
                try { fasta = decompress_gz(row.file_path); }
                catch (...) { n_missing++; continue; }
                size_t cnt = 0;
                parse_fasta_contig_accessions(fasta, [&](std::string_view acc) {
                    local.push_back({cidx_hash(acc), gid});
                    cnt++;
                });
                n_contigs += cnt;
                size_t d = ++done;
                if (d % 10000 == 0)
                    spdlog::info("CIDX: {}/{} genomes scanned, {} contigs so far", d, n_rows, n_contigs.load());
            }
        };

        std::vector<std::thread> thrs;
        thrs.reserve(n_threads);
        for (int t = 0; t < n_threads; t++) thrs.emplace_back(worker2, t);
        for (auto& t : thrs) t.join();

        for (auto& local : thread_entries)
            for (const auto& e : local)
                cidx_writer.add_hash(e.hash, static_cast<uint32_t>(e.gid));

        spdlog::info("CIDX: {} contigs indexed, {} genomes missing/skipped", n_contigs.load(), n_missing.load());
    }

    std::unique_ptr<SkchWriter>       skch_writer;
    std::unique_ptr<SkchWriterMultiK> skch_writer_mk;

    // ── Build SKCH for missing genomes (if --skch) ────────────────────────────
    size_t skch_n_missing = 0;

    // Normalise k list: --sketch-kmers takes precedence over --sketch-kmer.
    const bool multi_k = skch_kmers.size() > 1;

    if (need_skch) {
        // 1. Resolve target sketch parameters.
        //    Explicit CLI args take precedence; fall back to existing section params,
        //    then to built-in defaults (k=16, size=10000) for a fresh archive.
        uint32_t sk_sketch_size = (skch_size    > 0) ? static_cast<uint32_t>(skch_size)    : 10000;
        uint32_t sk_kmer_size   = (skch_kmer    > 0) ? static_cast<uint32_t>(skch_kmer)    : 16;
        uint32_t sk_syncmer_s;
        if (skch_syncmer > 0) {
            sk_syncmer_s = static_cast<uint32_t>(skch_syncmer);
        } else if (skch_syncmer == -1) {
            // auto: s = k/3, targeting ~6-8% syncmer density
            uint32_t ref_k = multi_k
                ? static_cast<uint32_t>(*std::min_element(skch_kmers.begin(), skch_kmers.end()))
                : sk_kmer_size;
            sk_syncmer_s = std::max(2u, ref_k / 3u);
        } else {
            sk_syncmer_s = 0;
        }
        uint64_t sk_seed1       = 42;
        uint64_t sk_seed2       = 43;

        auto skch_sections = toc.find_by_type(SEC_SKCH);
        std::vector<SkchReader> skch_readers;
        skch_readers.reserve(skch_sections.size());
        for (auto* sd : skch_sections) {
            skch_readers.emplace_back();
            skch_readers.back().open(mmap.data(), sd->file_offset, sd->compressed_size);
            // If no explicit override, inherit params from first existing section.
            if (skch_kmer <= 0 && !multi_k) sk_kmer_size  = skch_readers.back().kmer_size();
            if (skch_size <= 0)             sk_sketch_size = skch_readers.back().sketch_size();
        }

        if (multi_k) {
            std::string kmers_str;
            for (int k : skch_kmers) { if (!kmers_str.empty()) kmers_str += ','; kmers_str += std::to_string(k); }
            spdlog::info("SKCH: multi-k mode, k=[{}] size={} syncmer_s={}",
                         kmers_str, sk_sketch_size, sk_syncmer_s);
        } else {
            spdlog::info("SKCH: target params k={} size={} syncmer_s={}",
                         sk_kmer_size, sk_sketch_size, sk_syncmer_s);
        }

        // 2. Determine which genomes need sketching.
        //    Multi-k: covered if any existing v2 section has_kmer_size for ALL requested ks.
        //    Single-k: covered if exists in a matching section (same k, same size).
        auto params_match_single = [&](const SkchReader& r) {
            return r.kmer_size()  == sk_kmer_size
                && r.sketch_size() == sk_sketch_size;
        };
        std::vector<const SkchReader*> matching_readers;
        for (const auto& r : skch_readers) {
            if (force) continue;
            if (multi_k) {
                bool has_all = true;
                for (int k : skch_kmers) if (!r.has_kmer_size(static_cast<uint32_t>(k))) { has_all = false; break; }
                if (has_all && r.sketch_size() >= sk_sketch_size) matching_readers.push_back(&r);
            } else {
                if (params_match_single(r)) matching_readers.push_back(&r);
            }
        }

        struct MissingGenome { GenomeId gid; uint32_t shard_section_id; };
        std::vector<MissingGenome> missing;
        size_t total_genomes = 0;

        for (auto* sd : toc.find_by_type(SEC_SHRD)) {
            ShardReader shard;
            shard.open(mmap.data(), sd->file_offset, sd->compressed_size);
            for (auto* de = shard.dir_begin(); de != shard.dir_end(); ++de) {
                ++total_genomes;
                bool covered = false;
                for (const auto* r : matching_readers) {
                    if (r->contains(de->genome_id)) { covered = true; break; }
                }
                if (!covered)
                    missing.push_back({de->genome_id, static_cast<uint32_t>(sd->section_id)});
            }
        }
        skch_n_missing = missing.size();
        spdlog::info("SKCH: {}/{} genomes need sketching ({} matching sections exist)",
                     skch_n_missing, total_genomes, matching_readers.size());

        if (!missing.empty()) {
            // Initialise the appropriate writer.
            if (multi_k) {
                std::vector<uint32_t> ks_u32;
                for (int k : skch_kmers) ks_u32.push_back(static_cast<uint32_t>(k));
                skch_writer_mk = std::make_unique<SkchWriterMultiK>(
                    ks_u32, sk_sketch_size, sk_syncmer_s, sk_seed1, sk_seed2);
            } else {
                skch_writer = std::make_unique<SkchWriter>(
                    sk_sketch_size, sk_kmer_size, sk_syncmer_s, sk_seed1, sk_seed2);
            }

            // 3. Group missing genomes by shard section for sequential I/O.
            std::unordered_map<uint32_t, std::vector<GenomeId>> by_section;
            for (const auto& m : missing)
                by_section[m.shard_section_id].push_back(m.gid);

            size_t done = 0, failed = 0;
            for (auto* sd : toc.find_by_type(SEC_SHRD)) {
                auto it = by_section.find(static_cast<uint32_t>(sd->section_id));
                if (it == by_section.end()) continue;

                ShardReader shard;
                shard.open(mmap.data(), sd->file_offset, sd->compressed_size);

                const auto& gids = it->second;
                const int n_gids = static_cast<int>(gids.size());
                std::atomic<size_t> local_failed{0};

                if (multi_k) {
                    // Multi-k: sketch at each k per genome in parallel.
                    const size_t nk = skch_kmers.size();
                    struct MultiResult {
                        GenomeId gid;
                        uint64_t genome_length;
                        std::vector<std::vector<uint16_t>> sigs1;  // [ki]
                        std::vector<std::vector<uint16_t>> sigs2;  // [ki]
                        std::vector<uint32_t>              n_real;  // [ki]
                        std::vector<std::vector<uint64_t>> masks;   // [ki]
                    };
                    std::vector<MultiResult> results(static_cast<size_t>(n_gids));

                    #pragma omp parallel for schedule(dynamic, 1) num_threads(skch_threads)
                    for (int j = 0; j < n_gids; ++j) {
                        GenomeId gid = gids[static_cast<size_t>(j)];
                        try {
                            std::string fasta = shard.fetch_genome(gid);
                            if (fasta.empty()) { local_failed++; continue; }
                            MultiResult& mr = results[static_cast<size_t>(j)];
                            mr.gid = gid;
                            mr.sigs1.resize(nk);
                            mr.sigs2.resize(nk);
                            mr.n_real.resize(nk);
                            mr.masks.resize(nk);
                            for (size_t ki = 0; ki < nk; ++ki) {
                                auto sk = sketch_oph_dual_from_buffer(
                                    fasta.data(), fasta.size(),
                                    skch_kmers[ki],
                                    static_cast<int>(sk_sketch_size),
                                    static_cast<int>(sk_syncmer_s),
                                    sk_seed1, sk_seed2);
                                mr.genome_length = sk.genome_length;
                                mr.n_real[ki]    = sk.n_real_bins;
                                mr.masks[ki]     = sk.real_bins_bitmask;
                                const size_t ns = sk.signature1.size();
                                mr.sigs1[ki].resize(ns);
                                mr.sigs2[ki].resize(ns);
                                for (size_t si = 0; si < ns; ++si) {
                                    mr.sigs1[ki][si] = static_cast<uint16_t>(sk.signature1[si] >> 16);
                                    mr.sigs2[ki][si] = static_cast<uint16_t>(sk.signature2[si] >> 16);
                                }
                            }
                        } catch (...) {
                            local_failed++;
                        }
                    }
                    for (auto& mr : results) {
                        if (mr.sigs1.empty() || mr.sigs1[0].empty()) continue;
                        skch_writer_mk->add(mr.gid, mr.genome_length,
                                            mr.sigs1, mr.sigs2,
                                            mr.n_real, mr.masks);
                        ++done;
                    }
                } else {
                    // Single-k path (dual-seed).
                    std::vector<std::pair<GenomeId, OPHDualSketchResult>> results(static_cast<size_t>(n_gids));

                    #pragma omp parallel for schedule(dynamic, 1) num_threads(skch_threads)
                    for (int j = 0; j < n_gids; ++j) {
                        GenomeId gid = gids[static_cast<size_t>(j)];
                        try {
                            std::string fasta = shard.fetch_genome(gid);
                            if (fasta.empty()) { local_failed++; continue; }
                            auto sk = sketch_oph_dual_from_buffer(
                                fasta.data(), fasta.size(),
                                static_cast<int>(sk_kmer_size),
                                static_cast<int>(sk_sketch_size),
                                static_cast<int>(sk_syncmer_s),
                                sk_seed1, sk_seed2);
                            results[static_cast<size_t>(j)] = {gid, std::move(sk)};
                        } catch (...) {
                            local_failed++;
                        }
                    }
                    for (auto& [gid, sk] : results) {
                        if (sk.signature1.empty()) continue;
                        const size_t n = sk.signature1.size();
                        std::vector<uint16_t> sig1_16(n), sig2_16(n);
                        for (size_t si = 0; si < n; ++si) {
                            sig1_16[si] = static_cast<uint16_t>(sk.signature1[si] >> 16);
                            sig2_16[si] = static_cast<uint16_t>(sk.signature2[si] >> 16);
                        }
                        skch_writer->add(gid, sig1_16, sig2_16,
                                         sk.n_real_bins, sk.genome_length,
                                         sk.real_bins_bitmask);
                        ++done;
                    }
                }

                failed += local_failed.load();
                if (done % 10000 < static_cast<size_t>(n_gids))
                    spdlog::info("SKCH: {}/{} genomes sketched ({} failed)",
                                 done, skch_n_missing, failed);
            }
            spdlog::info("SKCH: {} genomes sketched, {} failed", done, failed);
            if (done == 0) { skch_writer.reset(); skch_writer_mk.reset(); }
        }
    }

    // ── Build GSTX (genus sketch stats) + QUAL (per-genome quality) ─────────────
    GstxWriter gstx_writer_new;
    QualWriter qual_writer_new;
    bool need_gstx = false;
    bool need_qual = false;
    if (need_gstx_check || need_qual_check) {
        if (toc.find_by_type(SEC_SKCH).empty()) {
            spdlog::warn("GSTX: no SKCH section — run reindex --skch first");
        } else if (toc.find_by_type(SEC_TAXN).empty()) {
            spdlog::warn("GSTX: no TAXN section found");
        } else {
            // Open a fresh ArchiveReader (its own mmap; the existing mmap above stays open)
            ArchiveReader ar;
            ar.open(gpk);

            const std::vector<uint32_t> avail_k = ar.available_sketch_kmer_sizes();
            const uint32_t sk_sz = ar.sketch_sketch_size();
            const int nk = static_cast<int>(avail_k.size());

            if (nk == 0 || sk_sz == 0 || sk_sz != GSTX_BINS) {
                spdlog::warn("GSTX: sketch_size={} != GSTX_BINS={} — skipping", sk_sz, GSTX_BINS);
            } else {
                // acc→gid
                std::unordered_map<std::string, GenomeId> acc_to_gid;
                ar.scan_genome_accessions([&](std::string_view acc, GenomeId gid) {
                    if (!ar.is_deleted(gid)) acc_to_gid.emplace(acc, gid);
                });

                // genus→members
                struct GMember { GenomeId gid; std::string acc; };
                std::vector<std::string>                      genus_keys;
                std::unordered_map<std::string, size_t>       genus_idx_map;
                std::vector<std::vector<GMember>>             genus_members;
                genus_keys.reserve(65536);
                genus_members.reserve(65536);

                ar.scan_taxonomy([&](std::string_view acc, std::string_view tax) {
                    // Key GSTX by species; fall back to genus if species is absent/empty.
                    constexpr std::string_view s_needle = ";s__";
                    constexpr std::string_view g_needle = ";g__";
                    std::string_view gsv;
                    auto try_extract = [&](std::string_view needle, std::string_view prefix) {
                        auto pos = tax.find(needle);
                        if (pos != std::string_view::npos) {
                            auto s = pos + 1;
                            auto e = tax.find(';', s);
                            gsv = tax.substr(s, e == std::string_view::npos ? e : e - s);
                        } else if (tax.starts_with(prefix)) {
                            auto e = tax.find(';', prefix.size());
                            gsv = tax.substr(0, e);
                        }
                    };
                    try_extract(s_needle, "s__");
                    if (gsv.empty() || gsv == "s__") { gsv = {}; try_extract(g_needle, "g__"); }
                    if (gsv.empty() || gsv == "g__") return;
                    auto it = acc_to_gid.find(std::string(acc));
                    if (it == acc_to_gid.end() || it->second == 0) return;
                    std::string gkey(gsv);
                    auto git = genus_idx_map.find(gkey);
                    if (git == genus_idx_map.end()) {
                        size_t idx = genus_keys.size();
                        genus_idx_map[gkey] = idx;
                        genus_keys.push_back(gkey);
                        genus_members.emplace_back();
                        git = genus_idx_map.find(gkey);
                    }
                    genus_members[git->second].push_back({it->second, std::string(acc)});
                });

                const size_t n_genera = genus_keys.size();
                spdlog::info("GSTX: {} genera, k=[{}], sketch_sz={}",
                             n_genera, [&]{ std::string s; for (size_t i=0;i<avail_k.size();++i){if(i)s+=',';s+=std::to_string(avail_k[i]);}return s;}(), sk_sz);

                // Sort each genus by gid for sequential SKCH access
                for (auto& members : genus_members)
                    std::sort(members.begin(), members.end(),
                              [](const GMember& a, const GMember& b) { return a.gid < b.gid; });

                // Dense re-index: only genera with ≥2 members get cand/vcnt/conts slots.
                // Singletons are skipped in GSTX/QUAL output anyway; allocating for them
                // was the source of the 100+ GB peak RSS.
                std::vector<size_t> multi_gi; // genus_members indices with size ≥ 2
                multi_gi.reserve(n_genera / 4);
                for (size_t gi = 0; gi < n_genera; ++gi)
                    if (genus_members[gi].size() >= 2) multi_gi.push_back(gi);
                const size_t n_multi = multi_gi.size();

                // gid → dense multi-genus index (singletons absent → ignored in callbacks)
                std::unordered_map<GenomeId, size_t> gid_to_gi;
                gid_to_gi.reserve(acc_to_gid.size());
                for (size_t mgi = 0; mgi < n_multi; ++mgi)
                    for (const auto& m : genus_members[multi_gi[mgi]])
                        gid_to_gi[m.gid] = mgi;

                // Only feed multi-genus gids to sketch_for_ids (saves I/O + memory)
                std::vector<GenomeId> all_sorted_gids;
                all_sorted_gids.reserve(gid_to_gi.size());
                for (const auto& [gid, mgi] : gid_to_gi) all_sorted_gids.push_back(gid);
                std::sort(all_sorted_gids.begin(), all_sorted_gids.end());

                // cand/vcnt only for multi-member genera.
                // vcnt uses int8_t: BM voting only needs sign, not magnitude (4× smaller).
                std::vector<std::vector<std::vector<uint16_t>>> cand(
                    n_multi, std::vector<std::vector<uint16_t>>(nk, std::vector<uint16_t>(sk_sz, 0)));
                std::vector<std::vector<std::vector<int8_t>>> vcnt(
                    n_multi, std::vector<std::vector<int8_t>>(nk, std::vector<int8_t>(sk_sz, 0)));

                // Pass 1: Boyer-Moore consensus — fused multi-k, one sequential scan
                ar.sketch_for_ids_multi_k(all_sorted_gids, sk_sz,
                    [&](size_t bidx, uint32_t ki, const SketchResult& sk) {
                        auto it = gid_to_gi.find(all_sorted_gids[bidx]);
                        if (it == gid_to_gi.end()) return;
                        const size_t mgi = it->second;
                        auto& c = cand[mgi][ki];
                        auto& v = vcnt[mgi][ki];
                        for (uint32_t b = 0; b < sk_sz; ++b) {
                            if (v[b] == 0)             { c[b] = sk.sig[b]; v[b] = 1; }
                            else if (sk.sig[b]==c[b])  { if (v[b] < 127) ++v[b]; }
                            else                       { --v[b]; }
                        }
                    }, 1);
                spdlog::info("GSTX: pass-1 done ({} k values fused)", nk);
                { decltype(vcnt) tmp; tmp.swap(vcnt); } // free vote counts — no longer needed

                ar.release_sketches(); // drop pass-1 SKCH pages before pass-2

                // Pass 2: containment vs consensus → p90 per multi-genus per k
                const bool have_kmrx = ar.kmer_profile(all_sorted_gids[0]) != nullptr;
                std::vector<std::vector<std::vector<float>>> conts(
                    n_multi, std::vector<std::vector<float>>(nk));
                for (size_t mgi = 0; mgi < n_multi; ++mgi)
                    for (int ki = 0; ki < nk; ++ki)
                        conts[mgi][ki].reserve(genus_members[multi_gi[mgi]].size());

                // Pass 2: containment vs consensus — fused multi-k, one sequential scan
                ar.sketch_for_ids_multi_k(all_sorted_gids, sk_sz,
                    [&](size_t bidx, uint32_t ki, const SketchResult& sk) {
                        auto it = gid_to_gi.find(all_sorted_gids[bidx]);
                        if (it == gid_to_gi.end()) return;
                        const size_t mgi = it->second;
                        uint32_t matches = 0;
                        const auto& c = cand[mgi][ki];
                        for (uint32_t b = 0; b < sk_sz; ++b)
                            if (sk.sig[b] == c[b]) ++matches;
                        conts[mgi][ki].push_back(
                            static_cast<float>(matches) / static_cast<float>(sk_sz));
                    }, 1);
                spdlog::info("GSTX: pass-2 done ({} k values fused)", nk);

                ar.release_sketches(); // drop pass-2 SKCH pages before TNF scan

                // Feed GstxWriter: p90 + TNF centroid + consensus
                for (size_t mgi = 0; mgi < n_multi; ++mgi) {
                    const size_t gi = multi_gi[mgi];
                    const auto& members = genus_members[gi];

                    // Exclude incomplete members from the p90 base — members with k=0
                    // containment < GSTX_P90_COMPLETE_FRAC * cluster_max are skipped so
                    // incomplete genomes don't deflate p90 and bias CCR at high completeness.
                    const auto& c0 = conts[mgi][0];
                    float max_c0 = 0.0f;
                    for (float v : c0) max_c0 = std::max(max_c0, v);
                    const float keep_thr = GSTX_P90_COMPLETE_FRAC * max_c0;

                    std::vector<float> p90(nk, 0.0f);
                    for (int ki = 0; ki < nk; ++ki) {
                        std::vector<float> cv;
                        cv.reserve(conts[mgi][ki].size());
                        for (size_t i = 0; i < conts[mgi][ki].size(); ++i)
                            if (i < c0.size() && c0[i] >= keep_thr)
                                cv.push_back(conts[mgi][ki][i]);
                        if (cv.empty()) cv = conts[mgi][ki]; // fallback: never empty p90
                        if (cv.empty()) continue;
                        std::sort(cv.begin(), cv.end());
                        size_t idx = static_cast<size_t>(
                            std::ceil(0.9f * static_cast<float>(cv.size())) - 1);
                        p90[ki] = cv[std::min(idx, cv.size() - 1)];
                    }

                    std::vector<float> tnf_sum(136, 0.0f);
                    uint32_t tnf_cnt = 0;
                    if (have_kmrx) {
                        for (const auto& m : members) {
                            const float* p = ar.kmer_profile(m.gid);
                            if (!p) continue;
                            for (int d = 0; d < 136; ++d) tnf_sum[d] += p[d];
                            ++tnf_cnt;
                        }
                    }
                    std::vector<float> tnf_mu(136, 0.0f);
                    if (tnf_cnt >= 2)
                        for (int d = 0; d < 136; ++d) tnf_mu[d] = tnf_sum[d] / tnf_cnt;

                    if (need_gstx_check)
                        gstx_writer_new.add_genus(
                            genus_keys[gi],
                            static_cast<uint32_t>(members.size()),
                            static_cast<uint32_t>(nk),
                            cand[mgi],
                            p90.data(),
                            tnf_cnt >= 2 ? tnf_mu.data() : nullptr,
                            avail_k.data());

                    if (need_qual_check) {
                        const float k0 = nk >= 1 ? static_cast<float>(avail_k[0]) : 0.0f;
                        const float k1 = nk >= 2 ? static_cast<float>(avail_k[1]) : 0.0f;
                        const float k2 = nk >= 3 ? static_cast<float>(avail_k[2]) : 0.0f;
                        for (size_t mi = 0; mi < members.size(); ++mi) {
                            auto qr           = QualRecord::make_empty(members[mi].gid);
                            qr.support_tier   = 0;
                            qr.interval_width = 0.05f;

                            auto meta = ar.genome_meta_by_accession(members[mi].acc);
                            uint32_t nc = meta ? meta->n_contigs : 0;
                            qr.completeness_fragmentation =
                                nc <= 1 ? 1.0f
                                : 1.0f / (1.0f + 0.333f * std::log2f(static_cast<float>(nc)));

                            float cont[GSTX_MAX_K] = {};
                            for (int ki = 0; ki < nk; ++ki)
                                if (mi < conts[mgi][ki].size())
                                    cont[ki] = conts[mgi][ki][mi];

                            if (p90[0] > 0.01f)
                                qr.completeness_cluster_relative =
                                    std::clamp(cont[0] / p90[0], 0.0f, 1.5f);

                            if (nk >= 3 && cont[0] >= 0.01f && cont[1] >= 0.01f && cont[2] >= 0.01f) {
                                const float lc0 = std::log(cont[0]), lc1 = std::log(cont[1]), lc2 = std::log(cont[2]);
                                const float beta    = (k0*lc0 + k1*lc1 + k2*lc2) / (k0*k0 + k1*k1 + k2*k2);
                                const float c2_pred = std::exp(beta * k2);
                                qr.contamination_leakage =
                                    std::max(0.0f, c2_pred - cont[2]) / std::max(c2_pred, 0.01f);
                                const float r0=lc0-beta*k0, r1=lc1-beta*k1, r2=lc2-beta*k2;
                                qr.leakage_residual = std::sqrt((r0*r0 + r1*r1 + r2*r2) / 2.0f);
                                // sketch_breadth field repurposed → fmh_minority_u8/contig_split_u8 (check-time only)
                            } else if (nk >= 2 && cont[0] >= 0.01f && cont[1] >= 0.01f) {
                                const float beta    = (k0*std::log(cont[0]) + k1*std::log(cont[1])) / (k0*k0 + k1*k1);
                                const float c1_pred = std::exp(beta * k1);
                                qr.contamination_leakage =
                                    std::max(0.0f, c1_pred - cont[1]) / std::max(c1_pred, 0.01f);
                            }

                            if (tnf_cnt >= 2) {
                                const float* tp = ar.kmer_profile(members[mi].gid);
                                if (tp) {
                                    float d2 = 0.0f, norm2 = 0.0f;
                                    for (int d = 0; d < 136; ++d) {
                                        float diff = tp[d] - tnf_mu[d];
                                        d2    += diff * diff;
                                        norm2 += tnf_mu[d] * tnf_mu[d];
                                    }
                                    const float dist = std::sqrt(d2);
                                    const float ref  = std::sqrt(norm2);
                                    if (ref > 1e-6f)
                                        qr.contamination_tnf_excess =
                                            std::max(0.0f, dist / ref - 1.0f);
                                }
                            }
                            qual_writer_new.add(qr);
                        }
                    }
                }

                { decltype(cand) tmp; tmp.swap(cand); }   // free consensus sketches
                { decltype(conts) tmp; tmp.swap(conts); } // free containment vectors

                // Singletons: no cluster → NaN for all cluster-derived fields,
                // but completeness_fragmentation is computable from n_contigs alone.
                if (need_qual_check) {
                    for (size_t gi = 0; gi < n_genera; ++gi) {
                        if (genus_members[gi].size() != 1) continue;
                        const auto& m = genus_members[gi][0];
                        auto qr       = QualRecord::make_empty(m.gid);
                        qr.support_tier   = static_cast<uint8_t>(2); // Singleton
                        qr.interval_width = 1.0f;
                        auto meta = ar.genome_meta_by_accession(m.acc);
                        if (meta) {
                            uint32_t nc = meta->n_contigs;
                            qr.completeness_fragmentation =
                                nc <= 1 ? 1.0f
                                : 1.0f / (1.0f + 0.333f * std::log2f(static_cast<float>(nc)));
                        }
                        qual_writer_new.add(qr);
                    }
                }

                if (need_gstx_check) {
                    need_gstx = true;
                    spdlog::info("GSTX: {} genera computed", gstx_writer_new.n_genera());
                }
                if (need_qual_check && qual_writer_new.size() > 0) {
                    need_qual = true;
                    spdlog::info("QUAL: {} genomes computed", qual_writer_new.size());
                }
                ar.close(); // drop all mmap pages (KMRX + residual) before write phase
            }
        }
    }

    // ── Rebuild GCOV/FCOV/FMHR (per-genus covariance + FMH) ─────────────────────
    // Full shard rescan — needed after merge, which drops these derived sections.
    // Computed before open_append so the reader sees the current, valid footer.
    std::optional<GcovFcovFmhrResult> gcov_res;
    if (need_gcov) {
        if (toc.find_by_type(SEC_TAXN).empty()) {
            spdlog::warn("GCOV: no TAXN section — skipping covariance rebuild");
            need_gcov = false;
        } else {
            ArchiveReader ar;
            ar.open(gpk);
            gcov_res = build_gcov_fcov_fmhr(ar, repack_threads > 0 ? repack_threads : 8);
            ar.close();
            if (gcov_res->gcov.n_genera() == 0) {
                spdlog::warn("GCOV: no genera produced — skipping");
                need_gcov = false;
                gcov_res.reset();
            } else {
                spdlog::info("GCOV: {} genera, {} families, {} fmhr genera",
                             gcov_res->gcov.n_genera(), gcov_res->fcov.n_genera(),
                             gcov_res->fmhr.n_genera());
            }
        }
    }

    // Copy existing sections into new TocWriter (excluding rebuilt types if --force)
    TocWriter new_toc;
    for (const auto& s : toc.sections) {
        if (force && s.type == SEC_GIDX) continue;
        if (force && need_txdb && s.type == SEC_TXDB) continue;
        if (force && need_cidx && s.type == SEC_CIDX) continue;
        if (force && need_gstx && s.type == SEC_GSTX) continue;
        if (force && need_qual && s.type == SEC_QUAL) continue;
        if (force && need_gcov && (s.type == SEC_GCOV || s.type == SEC_FCOV || s.type == SEC_FMHR)) continue;
        new_toc.add_section(s);
    }

    mmap.close();

    AppendWriter writer;
    writer.open_append(gpk);

    if (need_gidx) {
        SectionDesc gidx_sd = gidx_writer.finalize(writer, next_section_id++);
        new_toc.add_section(gidx_sd);
        spdlog::info("GIDX: written at offset {}, {} bytes",
                     gidx_sd.file_offset, gidx_sd.compressed_size);
    }

    if (need_txdb) {
        spdlog::info("TXDB: building taxonomy tree…");
        SectionDesc txdb_sd = txdb_writer.finalize(writer, next_section_id++);
        new_toc.add_section(txdb_sd);
        spdlog::info("TXDB: written at offset {}, {} bytes",
                     txdb_sd.file_offset, txdb_sd.compressed_size);
    }

    if (need_cidx && cidx_writer.size() > 0) {
        uint64_t batch_id = toc.next_section_id(); // use as monotonic batch id
        SectionDesc cidx_sd = cidx_writer.finalize(writer, next_section_id++, batch_id);
        new_toc.add_section(cidx_sd);
        spdlog::info("CIDX: written at offset {}, {} bytes ({} entries)",
                     cidx_sd.file_offset, cidx_sd.compressed_size, cidx_sd.item_count);
    }

    if (skch_writer) {
        SectionDesc skch_sd = skch_writer->finalize(writer, next_section_id++);
        new_toc.add_section(skch_sd);
        spdlog::info("SKCH: written at offset {}, {} bytes ({} sketches)",
                     skch_sd.file_offset, skch_sd.compressed_size, skch_sd.item_count);
    }
    if (skch_writer_mk) {
        SectionDesc skch_sd = skch_writer_mk->finalize(writer, next_section_id++);
        new_toc.add_section(skch_sd);
        spdlog::info("SKCH(multi-k): written at offset {}, {} bytes ({} sketches)",
                     skch_sd.file_offset, skch_sd.compressed_size, skch_sd.item_count);
    }

    if (need_gstx && gstx_writer_new.n_genera() > 0) {
        SectionDesc gstx_sd = gstx_writer_new.finalize(writer, next_section_id++);
        new_toc.add_section(gstx_sd);
        spdlog::info("GSTX: written at offset {}, {} bytes ({} genera)",
                     gstx_sd.file_offset, gstx_sd.compressed_size, gstx_sd.item_count);
    }

    if (need_qual && qual_writer_new.size() > 0) {
        SectionDesc qual_sd = qual_writer_new.finalize(writer, next_section_id++);
        new_toc.add_section(qual_sd);
        spdlog::info("QUAL: written at offset {}, {} bytes ({} genomes)",
                     qual_sd.file_offset, qual_sd.compressed_size, qual_sd.item_count);
    }

    if (need_gcov && gcov_res) {
        SectionDesc gcov_sd = gcov_res->gcov.finalize(writer, next_section_id++);
        new_toc.add_section(gcov_sd);
        spdlog::info("GCOV: written at offset {}, {} bytes ({} genera)",
                     gcov_sd.file_offset, gcov_sd.compressed_size, gcov_sd.item_count);
        if (gcov_res->fcov.n_genera() > 0)
            new_toc.add_section(gcov_res->fcov.finalize(writer, next_section_id++, SEC_FCOV));
        if (gcov_res->fmhr.n_genera() > 0)
            new_toc.add_section(gcov_res->fmhr.finalize(writer, next_section_id++));
    }

    // Stamp content checksums on the newly appended sections so verify can
    // validate them (mirrors build/merge). Read from the gpk itself — sections
    // were appended at absolute offsets; only_if_zero leaves the large
    // pre-existing SHRD/SKCH bodies (already checksummed) untouched.
    writer.flush();
    stamp_section_checksums(gpk.c_str(), new_toc.sections(), /*only_if_zero=*/true);

    // Write new TOC + TailLocator
    new_toc.finalize(writer,
                     prev_generation + 1,
                     live_count, total_count,
                     prev_toc_offset,
                     catalog_root_id, accession_root_id, tombstone_root_id);

    writer.flush();

    // If this reindex restored the genus-stats set that merge dropped (GSTX +
    // GCOV both present now), clear the FileHeader STATS_STALE marker so readers
    // stop warning about absent genus stats.
    {
        bool has_gstx_now = false, has_gcov_now = false;
        for (const auto& s : new_toc.sections()) {
            if (s.type == SEC_GSTX) has_gstx_now = true;
            if (s.type == SEC_GCOV) has_gcov_now = true;
        }
        if (has_gstx_now && has_gcov_now) {
            int ffd = ::open(gpk.c_str(), O_RDWR);
            if (ffd >= 0) {
                uint64_t flags = 0;
                if (::pread(ffd, &flags, sizeof(flags), offsetof(FileHeader, flags))
                        == static_cast<ssize_t>(sizeof(flags))
                    && (flags & FH_FLAG_STATS_STALE)) {
                    flags &= ~FH_FLAG_STATS_STALE;
                    ::pwrite(ffd, &flags, sizeof(flags), offsetof(FileHeader, flags));
                    ::fsync(ffd);
                }
                ::close(ffd);
            }
        }
    }

    spdlog::info("Reindex complete: gen {} -> {}", prev_generation, prev_generation + 1);
    return 0;
}

// ── genopack taxonomy patch ────────────────────────────────────────────────────

// Normalize a raw GTDB taxonomy string (7-rank, no l__/k__) to the canonical
// 10-rank format used by genopack (d,l,k,p,c,o,f,g,s,S). Already-normalized
// strings (containing l__) are returned unchanged.
static std::string normalize_gtdb_10rank(const std::string& tax, std::string_view acc) {
    // Already 10-rank?
    if (tax.find(";l__") != std::string::npos || tax.starts_with("l__")) return tax;

    // Derive stem for species/strain fallback
    if (acc.starts_with("RS_") || acc.starts_with("GB_")) acc = acc.substr(3);
    std::string stem;
    if (acc.starts_with("GCF_") || acc.starts_with("GCA_")) {
        auto dot = acc.rfind('.');
        stem = (dot != std::string_view::npos)
             ? std::string(acc.substr(0, dot)) : std::string(acc);
    } else {
        stem = std::string(acc);
    }

    static constexpr std::array<std::string_view, 10> kRanks =
        {"d__","l__","k__","p__","c__","o__","f__","g__","s__","S__"};

    std::unordered_map<std::string, std::string> rm;
    std::string_view sv(tax);
    while (!sv.empty()) {
        auto sep = sv.find(';');
        auto tok = sv.substr(0, sep);
        if (tok.size() >= 3 && tok[1] == '_' && tok[2] == '_')
            rm.emplace(std::string(tok.substr(0, 3)), std::string(tok));
        sv = (sep == std::string_view::npos) ? "" : sv.substr(sep + 1);
    }
    auto prop = [&](std::string_view c, std::string_view p) {
        std::string cs(c), ps(p);
        if (!rm.count(cs) || rm[cs] == cs)
            if (rm.count(ps)) rm[cs] = cs + rm[ps].substr(3);
    };
    prop("l__","d__"); prop("k__","l__"); prop("p__","k__"); prop("c__","p__");
    prop("o__","c__"); prop("f__","o__"); prop("g__","f__");
    auto& s = rm["s__"];
    if (s.empty() || s == "s__") s = "s__" + stem;
    rm["S__"] = "S__" + rm["s__"].substr(3);

    std::string res; res.reserve(tax.size() + 32);
    for (auto r : kRanks) {
        if (!res.empty()) res += ';';
        auto it = rm.find(std::string(r));
        res += (it != rm.end()) ? it->second : std::string(r);
    }
    return res;
}

// Patch the TAXN (and optionally TXDB) sections of a single .gpk file.
// patch_map: accession → new normalized taxonomy string.
// Returns the number of accessions actually patched in this file.
static size_t patch_gpk_taxn(const std::filesystem::path& gpk_path,
                              const std::unordered_map<std::string, std::string>& patch_map) {
    MmapFileReader mmap;
    mmap.open(gpk_path);
    if (mmap.size() < sizeof(FileHeader)) {
        spdlog::warn("Skipping {}: too small", gpk_path.string()); return 0;
    }
    auto* fh = mmap.ptr_at<FileHeader>(0);
    if (fh->magic != GPK_MAGIC) {
        spdlog::warn("Skipping {}: not a v2 .gpk", gpk_path.string()); return 0;
    }

    Toc toc = TocReader::read(mmap);
    auto taxn_sections = toc.find_by_type(SEC_TAXN);
    if (taxn_sections.empty()) {
        spdlog::debug("No TAXN sections in {}", gpk_path.string()); return 0;
    }

    // Load all existing TAXN entries into memory, apply patches.
    std::vector<std::pair<std::string, std::string>> entries; // acc → taxonomy
    size_t n_patched = 0;
    for (auto* sd : taxn_sections) {
        TaxonomyIndexReader tir;
        tir.open(mmap.data(), sd->file_offset, sd->compressed_size);
        tir.scan([&](std::string_view acc, std::string_view tax) {
            auto it = patch_map.find(std::string(acc));
            if (it != patch_map.end()) {
                entries.emplace_back(std::string(acc), it->second);
                ++n_patched;
            } else {
                entries.emplace_back(std::string(acc), std::string(tax));
            }
        });
    }
    if (n_patched == 0) return 0; // nothing to patch in this part

    // Capture TOC metadata before closing mmap.
    const auto* tail   = mmap.ptr_at<TailLocator>(mmap.size() - sizeof(TailLocator));
    uint64_t prev_gen  = toc.header.generation;
    uint64_t live      = toc.header.live_genome_count;
    uint64_t total     = toc.header.total_genome_count;
    uint64_t prev_toc  = tail->toc_offset;
    uint64_t cat_root  = toc.header.catalog_root_section_id;
    uint64_t acc_root  = toc.header.accession_root_section_id;
    uint64_t tomb_root = toc.header.tombstone_root_section_id;
    uint64_t next_sid  = toc.next_section_id();
    mmap.close();

    // Build new TOC: keep all non-TAXN/non-TXDB sections, add new TAXN + TXDB.
    TocWriter new_toc;
    for (const auto& s : toc.sections)
        if (s.type != SEC_TAXN && s.type != SEC_TXDB)
            new_toc.add_section(s);

    AppendWriter writer;
    writer.open_append(gpk_path);

    // Write new TAXN section.
    TaxonomyIndexWriter taxn_writer;
    for (const auto& [acc, tax] : entries)
        taxn_writer.add(acc, tax);
    auto taxn_sd = taxn_writer.finalize(writer, next_sid++);
    new_toc.add_section(taxn_sd);

    // Rebuild TXDB from the patched TAXN entries.
    TxdbWriter txdb_writer;
    for (const auto& [acc, tax] : entries)
        txdb_writer.add(acc, tax);
    auto txdb_sd = txdb_writer.finalize(writer, next_sid++);
    new_toc.add_section(txdb_sd);

    writer.flush();
    stamp_section_checksums(gpk_path.c_str(), new_toc.sections(), /*only_if_zero=*/true);
    new_toc.finalize(writer, prev_gen + 1, live, total, prev_toc,
                     cat_root, acc_root, tomb_root);
    writer.flush();
    return n_patched;
}

static int cmd_taxn_patch(
    const std::filesystem::path& archive_path,
    const std::filesystem::path& patch_tsv,
    const std::filesystem::path& input_tsv,
    const std::filesystem::path& tsv_out,
    bool gtdbtk_format,
    bool normalize)
{
    namespace fs = std::filesystem;

    // ── Step 1: load patch map ────────────────────────────────────────────
    // Auto-detect format from header line:
    //   GTDB-Tk: first field is "user_genome"
    //   taxonomy_changed.tsv: 3 fields, new taxonomy is col 2
    //   simple 2-col: accession TAB new_taxonomy
    spdlog::info("taxonomy patch: loading patch file {}", patch_tsv.string());
    std::ifstream fp(patch_tsv);
    if (!fp) { spdlog::error("Cannot open patch file: {}", patch_tsv.string()); return 1; }

    std::unordered_map<std::string, std::string> patch_map;  // acc → new normalized taxonomy
    std::string hdr, line;
    std::getline(fp, hdr);

    int tax_col = 1;  // default: simple 2-col
    bool skip_no_hit = false;
    if (gtdbtk_format || hdr.starts_with("user_genome")) {
        tax_col     = 1;    // classification column
        skip_no_hit = true; // GTDB-Tk reports "N/A" for unclassified
    } else {
        // Count tabs in header to detect 3-col taxonomy_changed.tsv
        size_t n_tabs = std::count(hdr.begin(), hdr.end(), '\t');
        if (n_tabs >= 2) tax_col = 2;
    }
    spdlog::info("taxonomy patch: format={}, taxonomy_col={}",
                 (tax_col == 1 && skip_no_hit) ? "gtdbtk" :
                 (tax_col == 2) ? "taxonomy_changed" : "simple",
                 tax_col);

    size_t n_loaded = 0;
    while (std::getline(fp, line)) {
        if (line.empty()) continue;
        std::vector<std::string_view> cols;
        std::string_view sv(line);
        while (!sv.empty()) {
            auto t = sv.find('\t');
            cols.push_back(sv.substr(0, t));
            if (t == std::string_view::npos) break;
            sv = sv.substr(t + 1);
        }
        if ((int)cols.size() <= tax_col) continue;
        std::string acc(cols[0]);
        std::string new_tax(cols[tax_col]);
        if (skip_no_hit && (new_tax == "N/A" || new_tax.empty())) continue;
        // Strip RS_/GB_ prefix from GTDB-Tk output accessions
        if (acc.starts_with("RS_") || acc.starts_with("GB_")) acc = acc.substr(3);
        // Normalize to 10-rank if requested or auto-detected as needed
        if (normalize || !new_tax.empty())
            new_tax = normalize_gtdb_10rank(new_tax, acc);
        patch_map.emplace(std::move(acc), std::move(new_tax));
        ++n_loaded;
    }
    spdlog::info("taxonomy patch: {} accessions to patch", n_loaded);
    if (n_loaded == 0) { spdlog::warn("Empty patch map — nothing to do"); return 0; }

    // ── Step 2: patch archive TAXN/TXDB sections ──────────────────────────
    if (!archive_path.empty()) {
        std::vector<fs::path> gpk_files;

        if (fs::is_directory(archive_path)) {
            // Multipart archive: collect all .gpk entries sorted.
            for (const auto& e : fs::directory_iterator(archive_path))
                if (e.path().extension() == ".gpk" &&
                    (e.is_directory() || e.is_regular_file() || e.is_symlink()))
                    gpk_files.push_back(e.path());
            std::sort(gpk_files.begin(), gpk_files.end());
            spdlog::info("taxonomy patch: {} archive parts in {}", gpk_files.size(), archive_path.string());
        } else {
            fs::path p = archive_path;
            if (p.extension() != ".gpk") p = fs::path(p.string() + ".gpk");
            gpk_files.push_back(p);
        }

        size_t total_patched = 0;
        for (const auto& gpk : gpk_files) {
            size_t n = patch_gpk_taxn(gpk, patch_map);
            spdlog::info("  {} → {} accessions patched", gpk.filename().string(), n);
            total_patched += n;
        }
        spdlog::info("taxonomy patch: {} total accessions patched across {} archive(s)",
                     total_patched, gpk_files.size());
    }

    // ── Step 3: patch flat input TSV (for geodesic) ───────────────────────
    if (!input_tsv.empty()) {
        std::ifstream fin(input_tsv);
        if (!fin) { spdlog::error("Cannot open input TSV: {}", input_tsv.string()); return 1; }

        fs::path out_path = tsv_out.empty() ? input_tsv : tsv_out;
        // Write to a temp file then rename to avoid clobbering on error.
        fs::path tmp_path = out_path.string() + ".patching";
        std::ofstream fout(tmp_path);
        if (!fout) { spdlog::error("Cannot write: {}", tmp_path.string()); return 1; }

        std::string hdr2, ln2;
        std::getline(fin, hdr2);
        fout << hdr2 << '\n';

        // Detect taxonomy column: header field "taxonomy" or col 1.
        int tsv_tax_col = 1;
        {
            std::string_view sv2(hdr2);
            int ci = 0;
            while (!sv2.empty()) {
                auto t = sv2.find('\t');
                if (sv2.substr(0, t) == "taxonomy") { tsv_tax_col = ci; break; }
                if (t == std::string_view::npos) break;
                sv2 = sv2.substr(t + 1); ++ci;
            }
        }

        size_t n_tsv_patched = 0, n_tsv_total = 0;
        while (std::getline(fin, ln2)) {
            if (ln2.empty()) { fout << '\n'; continue; }
            ++n_tsv_total;
            std::vector<std::string_view> cols;
            std::string_view sv2(ln2);
            while (!sv2.empty()) {
                auto t = sv2.find('\t');
                cols.push_back(sv2.substr(0, t));
                if (t == std::string_view::npos) break;
                sv2 = sv2.substr(t + 1);
            }
            // Strip prefix from accession for lookup
            std::string acc(cols[0]);
            if (acc.starts_with("RS_") || acc.starts_with("GB_")) acc = acc.substr(3);
            auto it = patch_map.find(acc);
            if (it != patch_map.end() && tsv_tax_col < (int)cols.size()) {
                // Reconstruct the line with the patched taxonomy column.
                for (int ci = 0; ci < (int)cols.size(); ++ci) {
                    if (ci > 0) fout << '\t';
                    if (ci == tsv_tax_col) fout << it->second;
                    else fout << cols[ci];
                }
                fout << '\n';
                ++n_tsv_patched;
            } else {
                fout << ln2 << '\n';
            }
            if (n_tsv_total % 500'000 == 0)
                spdlog::info("  TSV: {} rows processed...", n_tsv_total);
        }
        fout.close();
        fs::rename(tmp_path, out_path);
        spdlog::info("taxonomy patch: TSV — {} of {} rows patched → {}",
                     n_tsv_patched, n_tsv_total, out_path.string());
    }

    return 0;
}

// ── genopack taxdump ───────────────────────────────────────────────────────────
//
// Binary columnar file formats (standalone, no genopack dependency to read):
//
// acc2taxid.bin:
//   [Acc2TaxidHeader 32B] [Acc2TaxidEntry×n sorted by acc_hash]
//   Binary-search by FNV-1a-64(accession) → taxid
//
// taxnodes.bin:
//   [TaxnodesHeader 32B] [TaxnodeEntry×n sorted by taxid] [name_pool bytes]
//   Binary-search by taxid → (parent, rank, name)

struct Acc2TaxidHeader {
    uint32_t magic;       // 'G','P','A','T' = 0x54415047
    uint32_t version;     // 1
    uint64_t n_entries;
    uint8_t  reserved[16];
};
static_assert(sizeof(Acc2TaxidHeader) == 32);

struct Acc2TaxidEntry {
    uint64_t acc_hash;    // FNV-1a-64(accession)
    uint64_t taxid;       // uint64_t: NCBI taxids fit in low 32 bits; GTDB has bit 63 set
};
static_assert(sizeof(Acc2TaxidEntry) == 16);

struct TaxnodesHeader {
    uint32_t magic;           // 'G','P','T','N' = 0x4E545047
    uint32_t version;         // 1
    uint32_t n_nodes;
    uint32_t name_pool_size;
    uint8_t  reserved[16];
};
static_assert(sizeof(TaxnodesHeader) == 32);

struct TaxnodeEntry {
    uint64_t taxid;         // uint64_t: GTDB concept_ids have bit 63 set
    uint64_t parent_taxid;
    uint32_t name_offset;   // byte offset into name pool
    uint8_t  rank;          // TaxRank value
    uint8_t  flags;         // bit 0: synthetic node
    uint16_t name_len;
};
static_assert(sizeof(TaxnodeEntry) == 24);

static constexpr uint32_t TAXDUMP_ACC2TAXID_MAGIC = 0x54415047u;  // "GPAT"
static constexpr uint32_t TAXDUMP_TAXNODES_MAGIC  = 0x4E545047u;  // "GPTN"

static const char* rank_name_ncbi(TaxRank r) {
    switch (r) {
        case TaxRank::DOMAIN:  return "superkingdom";
        case TaxRank::PHYLUM:  return "phylum";
        case TaxRank::CLASS:   return "class";
        case TaxRank::ORDER:   return "order";
        case TaxRank::FAMILY:  return "family";
        case TaxRank::GENUS:   return "genus";
        case TaxRank::SPECIES: return "species";
        default:               return "no rank";
    }
}

static const char* rank_name_gtdb(TaxRank r) {
    switch (r) {
        case TaxRank::DOMAIN:  return "domain";
        case TaxRank::PHYLUM:  return "phylum";
        case TaxRank::CLASS:   return "class";
        case TaxRank::ORDER:   return "order";
        case TaxRank::FAMILY:  return "family";
        case TaxRank::GENUS:   return "genus";
        case TaxRank::SPECIES: return "species";
        default:               return "no_rank";
    }
}

static int cmd_taxdump(const std::string& archive_path,
                       const std::string& format,
                       const std::string& output_dir) {
    ArchiveReader ar;
    ar.open(archive_path);

    auto tree_opt = ar.taxonomy_tree();
    if (!tree_opt) {
        spdlog::error("No taxonomy data in archive");
        return 1;
    }
    auto& tree = *tree_opt;

    std::filesystem::create_directories(output_dir);
    std::filesystem::path out = output_dir;

    spdlog::info("Taxonomy: {} nodes, {} accessions", tree.n_nodes(), tree.n_accessions());

    if (format == "taxdump") {
        // ── NCBI taxdump ─────────────────────────────────────────────────────
        {
            std::ofstream f(out / "names.dmp");
            if (!f) throw std::runtime_error("Cannot open names.dmp");
            tree.scan_nodes([&](uint64_t taxid, uint64_t, TaxRank,
                                std::string_view name, bool) {
                f << taxid << "\t|\t" << name << "\t|\t\t|\tscientific name\t|\n";
            });
        }
        {
            std::ofstream f(out / "nodes.dmp");
            if (!f) throw std::runtime_error("Cannot open nodes.dmp");
            tree.scan_nodes([&](uint64_t taxid, uint64_t parent, TaxRank rank,
                                std::string_view, bool) {
                uint64_t p = (parent == 0) ? taxid : parent;
                f << taxid << "\t|\t" << p << "\t|\t" << rank_name_ncbi(rank)
                  << "\t|\t\t|\t0\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n";
            });
        }
        {
            std::ofstream f(out / "acc2taxid.dmp");
            if (!f) throw std::runtime_error("Cannot open acc2taxid.dmp");
            f << "accession\taccession.version\ttaxid\tgi\n";
            tree.scan_accessions([&](std::string_view acc, uint64_t taxid) {
                f << acc << "\t" << acc << "\t" << taxid << "\t0\n";
            });
        }
        // Stub files expected by tools like Kraken
        { std::ofstream(out / "merged.dmp");   }
        { std::ofstream(out / "delnodes.dmp"); }
        spdlog::info("taxdump written to {}", output_dir);

    } else if (format == "columnar") {
        // ── Binary columnar ──────────────────────────────────────────────────

        // acc2taxid.bin: sorted (FNV-1a-64(acc), taxid) entries
        {
            std::vector<Acc2TaxidEntry> entries;
            entries.reserve(tree.n_accessions());
            tree.scan_accessions([&](std::string_view acc, uint64_t taxid) {
                entries.push_back({cidx_hash(acc), taxid});
            });
            std::sort(entries.begin(), entries.end(),
                      [](const Acc2TaxidEntry& a, const Acc2TaxidEntry& b) {
                          return a.acc_hash < b.acc_hash;
                      });
            std::ofstream f(out / "acc2taxid.bin", std::ios::binary);
            if (!f) throw std::runtime_error("Cannot open acc2taxid.bin");
            Acc2TaxidHeader hdr{};
            hdr.magic     = TAXDUMP_ACC2TAXID_MAGIC;
            hdr.version   = 1;
            hdr.n_entries = entries.size();
            f.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));
            f.write(reinterpret_cast<const char*>(entries.data()),
                    static_cast<std::streamsize>(entries.size() * sizeof(Acc2TaxidEntry)));
            spdlog::info("acc2taxid.bin: {} entries, {} bytes",
                         entries.size(),
                         sizeof(hdr) + entries.size() * sizeof(Acc2TaxidEntry));
        }

        // taxnodes.bin: sorted TaxnodeEntry array + name pool
        {
            std::vector<TaxnodeEntry> nodes;
            std::string name_pool;
            nodes.reserve(tree.n_nodes());
            tree.scan_nodes([&](uint64_t taxid, uint64_t parent, TaxRank rank,
                                std::string_view name, bool synthetic) {
                TaxnodeEntry e{};
                e.taxid        = taxid;
                e.parent_taxid = (parent == 0) ? taxid : parent;
                e.name_offset  = static_cast<uint32_t>(name_pool.size());
                e.rank         = static_cast<uint8_t>(rank);
                e.flags        = synthetic ? 1u : 0u;
                e.name_len     = static_cast<uint16_t>(name.size());
                name_pool.append(name);
                name_pool.push_back('\0');
                nodes.push_back(e);
            });
            std::sort(nodes.begin(), nodes.end(),
                      [](const TaxnodeEntry& a, const TaxnodeEntry& b) {
                          return a.taxid < b.taxid;
                      });
            std::ofstream f(out / "taxnodes.bin", std::ios::binary);
            if (!f) throw std::runtime_error("Cannot open taxnodes.bin");
            TaxnodesHeader hdr{};
            hdr.magic          = TAXDUMP_TAXNODES_MAGIC;
            hdr.version        = 1;
            hdr.n_nodes        = static_cast<uint32_t>(nodes.size());
            hdr.name_pool_size = static_cast<uint32_t>(name_pool.size());
            f.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));
            f.write(reinterpret_cast<const char*>(nodes.data()),
                    static_cast<std::streamsize>(nodes.size() * sizeof(TaxnodeEntry)));
            f.write(name_pool.data(), static_cast<std::streamsize>(name_pool.size()));
            spdlog::info("taxnodes.bin: {} nodes, {} bytes name pool",
                         nodes.size(), name_pool.size());
        }

        // TSV sidecars for human inspection / pandas / polars
        {
            std::ofstream f(out / "acc2taxid.tsv");
            if (!f) throw std::runtime_error("Cannot open acc2taxid.tsv");
            f << "accession\ttaxid\n";
            tree.scan_accessions([&](std::string_view acc, uint64_t taxid) {
                f << acc << "\t" << taxid << "\n";
            });
        }
        {
            std::ofstream f(out / "taxonomy.tsv");
            if (!f) throw std::runtime_error("Cannot open taxonomy.tsv");
            f << "taxid\tparent_taxid\trank\tname\tis_synthetic\n";
            tree.scan_nodes([&](uint64_t taxid, uint64_t parent, TaxRank rank,
                                std::string_view name, bool synthetic) {
                uint64_t p = (parent == 0) ? taxid : parent;
                f << taxid << "\t" << p << "\t" << rank_name_gtdb(rank)
                  << "\t" << name << "\t" << (synthetic ? 1 : 0) << "\n";
            });
        }
        spdlog::info("columnar written to {}", output_dir);

    } else {
        spdlog::error("Unknown format '{}'. Use 'taxdump' or 'columnar'.", format);
        return 1;
    }

    return 0;
}

// ── genopack repack ────────────────────────────────────────────────────────────
static int cmd_repack(const std::string& input, const std::string& output,
                      int zstd_level, const std::string& taxonomy_rank,
                      int threads, int max_mem_gb, bool verbose) {
    RepackConfig cfg;
    cfg.taxonomy_rank              = taxonomy_rank.empty() ? 'g' : taxonomy_rank[0];
    cfg.shard_cfg.zstd_level       = zstd_level;
    cfg.shard_cfg.auto_codec       = false;
    cfg.shard_cfg.use_delta        = false;
    cfg.shard_cfg.compress_threads = static_cast<size_t>(std::max(1, threads));
    cfg.threads                    = static_cast<size_t>(std::max(1, threads));
    cfg.max_bucket_bytes           = static_cast<uint64_t>(std::max(1, max_mem_gb)) << 30;
    cfg.verbose                    = verbose;
    repack_archive(input, output, cfg);
    return 0;
}

// ── genopack verify ──────────────────────────────────────────────────────────
static int cmd_verify(const std::string& archive_path, bool verbose, bool strict) {
    MmapFileReader mmap;
    mmap.open(archive_path);
    if (mmap.size() < sizeof(TailLocator))
        throw std::runtime_error("verify: file too small");

    // Read TailLocator from end of file.
    uint64_t tail_off = mmap.size() - sizeof(TailLocator);
    const auto* tail = mmap.ptr_at<TailLocator>(tail_off);
    if (tail->magic != GPKT_MAGIC)
        throw std::runtime_error("verify: TailLocator magic mismatch");

    int failures = 0;

    static const uint8_t zeros16[16] = {};
    bool toc_checksums_present = std::memcmp(tail->toc_checksum, zeros16, 16) != 0;

    // --- Verify toc_checksum (covers TOC block bytes as written) ---
    if (toc_checksums_present) {
        const uint8_t* toc_bytes = mmap.data() + tail->toc_offset;
        XXH128_hash_t h = XXH3_128bits(toc_bytes, tail->toc_size);
        XXH128_canonical_t canon;
        XXH128_canonicalFromHash(&canon, h);
        bool ok = std::memcmp(tail->toc_checksum, canon.digest, 16) == 0;
        if (!ok) { spdlog::error("verify: TailLocator.toc_checksum MISMATCH"); ++failures; }
        else if (verbose) spdlog::info("verify: toc_checksum OK");
    } else {
        spdlog::warn("verify: TailLocator.toc_checksum is zero (pre-feature archive)");
    }

    // Read TOC to get section list.
    Toc toc = TocReader{}.read(mmap);

    // --- Verify TocHeader.checksum (covers header + SectionDescs with checksum zeroed) ---
    if (toc_checksums_present) {
        const uint8_t* toc_bytes = mmap.data() + tail->toc_offset;
        bool ok = verify_checksum(toc_bytes, tail->toc_size, offsetof(TocHeader, checksum));
        if (!ok) { spdlog::error("verify: TocHeader.checksum MISMATCH"); ++failures; }
        else if (verbose) spdlog::info("verify: TocHeader.checksum OK");
    }

    // --- Verify every section's content checksum (SectionDesc.checksum) ---
    // The descriptor checksum is the XXH128 of the whole section body. Sections
    // with an all-zero checksum predate the feature (or exceeded the build-time
    // hashing cap) and are skipped.
    size_t n_sections = toc.sections.size();
    size_t n_zero     = 0;
    size_t n_checked  = 0;
    static const uint8_t zeros[16] = {};
    for (const auto& sd : toc.sections) {
        if (std::memcmp(sd.checksum, zeros, 16) == 0) {
            ++n_zero;
            if (strict) {
                spdlog::error("verify: section type={:#x} id={} has no content checksum (--strict)",
                              sd.type, sd.section_id);
                ++failures;
            }
            continue;
        }
        if (sd.compressed_size > mmap.size() ||
            sd.file_offset > mmap.size() - sd.compressed_size) {
            spdlog::error("verify: section type={:#x} id={} out of bounds",
                          sd.type, sd.section_id);
            ++failures; continue;
        }
        const uint8_t* sec = mmap.data() + sd.file_offset;
        if (!checksum_matches(sec, sd.compressed_size, sd.checksum)) {
            spdlog::error("verify: section type={:#x} id={} checksum MISMATCH",
                          sd.type, sd.section_id);
            ++failures;
        } else {
            ++n_checked;
            if (verbose) spdlog::info("verify: section type={:#x} id={} deriv={:#018x} schema={:#018x} OK",
                                      sd.type, sd.section_id, sd.derivation_hash, sd.semantic_schema_hash);
        }
    }

    if (n_zero > 0) {
        if (strict)
            spdlog::error("verify: {} of {} sections have no content checksum "
                          "(--strict: counted as failures)", n_zero, n_sections);
        else
            spdlog::warn("verify: {} of {} sections have no content checksum "
                         "(index/metadata or >512MB section bodies — bytes unverified)",
                         n_zero, n_sections);
    }

    if (failures == 0)
        spdlog::info("verify: {} section(s) content-checked, {} skipped (no checksum)",
                     n_checked, n_zero);
    else
        spdlog::error("verify: {} checksum failure(s)", failures);

    return failures == 0 ? 0 : 1;
}

// ── main ──────────────────────────────────────────────────────────────────────
int main(int argc, char** argv) {
    CLI::App app{"genopack — genome archive"};
    app.require_subcommand(1);

    // genopack build
    auto* build = app.add_subcommand("build", "Build a new archive from a genome TSV (default codec mode auto-selects plain vs delta)");
    std::string build_input, build_output;
    int build_threads = 16, build_level = 6, build_parallel = 1;
    bool build_no_dict = false, build_ref_dict = false, build_delta = false;
    bool build_mem_delta = true, build_verbose = false, build_no_cidx = true;
    bool build_2bit = true, build_kmer_nn = true, build_taxon_group = true;
    bool build_sketch = true;
    int build_sketch_kmer = 16, build_sketch_size = 10000, build_sketch_syncmer = -1;
    std::string build_taxon_rank = "g";
    std::string build_sketch_kmers_str = "16,21,31";
    build->add_option("-i,--input",  build_input,  "Input TSV (accession, file_path, ...). Not required with --from-gpk.");
    build->add_option("-o,--output", build_output, "Output archive directory (.gpk)")->required();
    build->add_option("-t,--threads", build_threads, "I/O threads (decompression + compression)");
    build->add_option("-z,--zstd-level", build_level, "zstd compression level (1-22)");
    build->add_option("-p,--parallel", build_parallel, "Number of parallel build workers (auto-merge)");
    build->add_flag("--no-dict",    build_no_dict,    "Disable shared dictionary training");
    build->add_flag("--ref-dict",   build_ref_dict,   "Use first genome in each shard as reference content dictionary");
    build->add_flag("--delta",      build_delta,      "Compress non-reference blobs against first genome using zstd prefix");
    build->add_flag("--mem-delta,!--no-mem-delta", build_mem_delta,  "MEM-delta k-mer exact-match encoding (default: on; --no-mem-delta to disable)");
    build->add_flag("--no-cidx,!--cidx",          build_no_cidx,    "Skip CIDX contig index (default: on; --cidx to build it)");
    bool build_no_gstx = false;
    build->add_flag("--no-gstx,!--gstx",         build_no_gstx,    "Skip GSTX genus-stats index (default: on; --no-gstx to disable)");
    uint32_t build_micro_genus_threshold = 0;
    build->add_option("--micro-genus-threshold", build_micro_genus_threshold,
        "Min genomes for a genus to get its own shard + GSTX/GCOV/FMHR consensus model; smaller "
        "genera are bin-packed and unmodeled. 0 = auto (scales with corpus: 4 if <=50k genomes, "
        "8 if <=500k, 16 if <=2M, else 32). Lower = more genera modeled, more shards (default: auto)")
        ->default_val(0);
    build->add_flag("--2bit,!--no-2bit",          build_2bit,       "2-bit sequence packing before zstd (default: on; --no-2bit to disable)");
    build->add_flag("--kmer-sort,!--no-kmer-sort",   build_kmer_nn,    "Sort genomes within each shard by kmer4_profile NN chain (default: on; --no-kmer-sort to disable)");
    build->add_flag("--taxon-group,!--no-taxon-group",build_taxon_group,"Group genomes into per-taxon shards (default: on; --no-taxon-group to disable; requires taxonomy column)");
    build->add_option("--taxon-rank",build_taxon_rank,"Taxonomy rank for grouping: g=genus (default), f=family");
    build->add_flag("--sketch,!--no-sketch", build_sketch, "Compute OPH sketches (default: on; use --no-sketch to disable)");
    build->add_option("--sketch-kmer", build_sketch_kmer, "OPH sketch k-mer size (default: 16)");
    build->add_option("--sketch-kmers", build_sketch_kmers_str, "Comma-separated k-mer sizes for multi-k SKCH (default: 16,21,31)");
    build->add_option("--sketch-size", build_sketch_size, "Number of OPH bins (default: 10000)");
    build->add_option("--sketch-syncmer", build_sketch_syncmer, "Open syncmer prefilter s (0=disabled, -1=auto: s=k/3; default: auto)");
    build->add_flag("-v,--verbose", build_verbose, "Verbose progress");
    std::string build_markers;
    build->add_option("--markers", build_markers,
        "Path to markers.mrk; build-time marker completeness scoring. If omitted, "
        "auto-resolved from $GENOPACK_DATA/markers.mrk or the install share dir.");
    std::string build_contam_panel;
    build->add_option("--contam-panel", build_contam_panel,
        "Path to contamination .csp panel; build-time intra-genome duplication scoring. "
        "If omitted, auto-resolved from $GENOPACK_DATA/contamination.csp or the install "
        "share dir; absent → axis skipped.");
    std::string build_from_gpk;
    build->add_option("--from-gpk", build_from_gpk,
        "Rebuild from an existing .gpk: stream decoded sequence from its shards "
        "(sequential reads) instead of opening per-genome FASTA files on NFS. -i not required.");
    std::string build_from_stage;
    build->add_option("--from-stage", build_from_stage,
        "Rebuild from a local .gstage cache produced by `genopack stage`: stream "
        "decoded sequence sequentially (no per-genome NFS opens). -i not required.");
    std::string build_coordinator; // "manifest_dir:output.gpk" or empty
    build->add_option("--coordinator", build_coordinator,
        "NFS manifest coordinator: manifest_dir:/path/to/output.gpk. "
        "Build to a temp file then transfer sections via NFS manifest protocol. "
        "The legacy 'nfs:' prefix is accepted and stripped for backward compatibility.");
    build->callback([&]() {
        {
            int build_mode_count = (!build_input.empty()) + (!build_from_gpk.empty()) + (!build_from_stage.empty());
            if (build_mode_count != 1)
                throw std::runtime_error("build: provide exactly one of -i <TSV>, --from-gpk <source.gpk>, or --from-stage <cache.gstage>");
        }
        // Parse --sketch-kmers "16,21,31" into vector<int>
        std::vector<int> build_sketch_kmers;
        if (!build_sketch_kmers_str.empty()) {
            std::istringstream ss(build_sketch_kmers_str);
            std::string tok;
            while (std::getline(ss, tok, ','))
                if (!tok.empty()) build_sketch_kmers.push_back(std::stoi(tok));
        }
        // Resolve --sketch-syncmer auto (-1): s = min(k)/3
        if (build_sketch_syncmer == -1) {
            int ref_k = build_sketch_kmers.empty()
                ? build_sketch_kmer
                : *std::min_element(build_sketch_kmers.begin(), build_sketch_kmers.end());
            build_sketch_syncmer = std::max(2, ref_k / 3);
            spdlog::info("SKCH: syncmer auto → s={} (k={})", build_sketch_syncmer, ref_k);
        }
        if (!build_coordinator.empty()) {
            // Coordinator mode: build to temp file, then transfer via NFS manifest.
            std::string tmp_output = build_output + ".coord_tmp";
            int rc = cmd_build(build_input, tmp_output,
                                build_threads, build_level, build_no_dict, build_ref_dict,
                                build_delta, build_mem_delta, build_verbose, build_parallel,
                                build_no_cidx,
                                build_2bit, build_kmer_nn,
                                build_taxon_group, build_taxon_rank,
                                build_sketch, build_sketch_kmer, build_sketch_size,
                                build_sketch_syncmer, build_sketch_kmers,
                                build_no_gstx, build_markers, build_from_gpk,
                                build_micro_genus_threshold, build_from_stage,
                                build_contam_panel);
            if (rc != 0) std::exit(rc);
            std::string hostname = "worker";
            {
                char buf[256] = {};
                if (::gethostname(buf, sizeof(buf)) == 0) hostname = buf;
            }
            hostname += "_" + std::to_string(::getpid());
            std::string tmp_gpk = tmp_output + ".gpk";
            // Strip legacy "nfs:" prefix if present
            std::string coord = build_coordinator;
            if (coord.substr(0, 4) == "nfs:")
                coord = coord.substr(4);
            auto colon = coord.find(':');
            if (colon == std::string::npos)
                throw std::runtime_error("--coordinator: expected manifest_dir:/output.gpk");
            std::filesystem::path manifest_dir = coord.substr(0, colon);
            std::filesystem::path out_path     = coord.substr(colon + 1);
            genopack::transfer_nfs(tmp_gpk, manifest_dir, out_path, hostname);
            std::exit(0);
        }
        std::exit(cmd_build(build_input, build_output,
                             build_threads, build_level, build_no_dict, build_ref_dict,
                             build_delta, build_mem_delta, build_verbose, build_parallel,
                             build_no_cidx,
                             build_2bit, build_kmer_nn,
                             build_taxon_group, build_taxon_rank,
                             build_sketch, build_sketch_kmer, build_sketch_size,
                             build_sketch_syncmer, build_sketch_kmers,
                             build_no_gstx, build_markers, build_from_gpk,
                             build_micro_genus_threshold, build_from_stage,
                             build_contam_panel));
    });

    // genopack stage — durable local sequence cache for fast iterative rebuilds
    auto* stage_cmd = app.add_subcommand("stage",
        "Cache NFS genome FASTAs as a local zstd sequence store (.gstage) for fast "
        "rebuilds via `build --from-stage`");
    std::string stage_input, stage_output;
    int stage_threads = 48, stage_block_mb = 64;
    stage_cmd->add_option("-i,--input", stage_input,
        "Input TSV (accession, file_path, taxonomy, ...)")->required();
    stage_cmd->add_option("-o,--output", stage_output, "Output .gstage path")->required();
    stage_cmd->add_option("-t,--threads", stage_threads, "NFS reader threads (default: 48)");
    stage_cmd->add_option("--block-mb", stage_block_mb, "Uncompressed block size in MB (default: 64)");
    stage_cmd->callback([&]() {
        std::exit(genopack::cmd_stage(stage_input, stage_output, stage_threads, stage_block_mb));
    });

    // genopack extract
    auto* extract = app.add_subcommand("extract", "Extract genomes by quality or accession");
    std::string ext_archive, ext_out, ext_acc_file, ext_out_dir;
    std::vector<std::string> ext_accessions;
    float ext_min_comp = 0, ext_max_contam = 100;
    extract->add_option("archive", ext_archive, "Archive directory")->required();
    extract->add_option("--accession", ext_accessions, "Accession to extract (repeatable)");
    extract->add_option("--accessions-file", ext_acc_file, "File with one accession per line");
    extract->add_option("--min-completeness",  ext_min_comp,   "Minimum completeness %");
    extract->add_option("--max-contamination", ext_max_contam, "Maximum contamination %");
    extract->add_option("-o,--out", ext_out, "Output FASTA (default: stdout)");
    extract->add_option("--output-dir", ext_out_dir, "Write one {accession}.fa per genome to this directory");
    extract->callback([&]() {
        std::exit(cmd_extract(ext_archive, ext_accessions, ext_acc_file,
                              ext_min_comp, ext_max_contam, ext_out, ext_out_dir));
    });

    // genopack slice
    auto* slice = app.add_subcommand("slice", "Extract a subsequence by accession and sequence coordinates");
    std::string slice_archive, slice_accession;
    uint64_t slice_start = 0, slice_length = 0;
    bool slice_fasta = false;
    slice->add_option("archive", slice_archive, "Archive directory")->required();
    slice->add_option("accession", slice_accession, "Accession to slice")->required();
    slice->add_option("--start", slice_start, "0-based sequence start")->required();
    slice->add_option("--length", slice_length, "Number of bases to extract")->required();
    slice->add_flag("--fasta", slice_fasta, "Emit a FASTA-style header");
    slice->callback([&]() {
        std::exit(cmd_slice(slice_archive, slice_accession, slice_start, slice_length, slice_fasta));
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

    // genopack inspect
    auto* inspect = app.add_subcommand("inspect",
        "Report sketch layout and preload memory cost for a .gpk archive or directory of parts");
    std::string inspect_path;
    bool inspect_json = false;
    inspect->add_option("path", inspect_path,
        "Archive (.gpk) or directory containing one or more .gpk parts")->required();
    inspect->add_flag("--json", inspect_json, "Output JSON");
    inspect->callback([&]() {
        std::exit(cmd_inspect(inspect_path, inspect_json));
    });

    // genopack merge
    auto* merge_cmd = app.add_subcommand("merge", "Merge multiple .gpk archives into one");
    std::vector<std::string> merge_inputs;
    std::string merge_list;
    std::string merge_output;
    merge_cmd->add_option("inputs", merge_inputs, "Input .gpk archives")->expected(0, -1);
    merge_cmd->add_option("-l,--list", merge_list, "File with one .gpk path per line");
    merge_cmd->add_option("-o,--output", merge_output, "Output .gpk archive")->required();
    merge_cmd->callback([&]() {
        std::exit(cmd_merge(merge_inputs, merge_list, merge_output));
    });

    // genopack add
    auto* add_cmd = app.add_subcommand("add",
        "Append genomes to an existing archive. Updates CATL/ACCX/GIDX/TAXN/TXDB/CIDX/KMRX immediately.");
    std::string add_archive, add_input;
    add_cmd->add_option("archive", add_archive, "Archive directory")->required();
    add_cmd->add_option("-i,--input", add_input, "Input TSV")->required();
    add_cmd->callback([&]() {
        ArchiveAppender app2(add_archive);
        app2.add_from_tsv(add_input);
        app2.commit();
        std::exit(0);
    });

    // genopack dedup
    auto* dedup_cmd = app.add_subcommand("dedup", "Remove duplicate genomes (same sequence, different accession)");
    std::string dedup_archive;
    bool dedup_dry_run = false;
    dedup_cmd->add_option("archive", dedup_archive, "Path to .gpk archive")->required();
    dedup_cmd->add_flag("--dry-run", dedup_dry_run, "Report duplicates without tombstoning");
    dedup_cmd->callback([&]() {
        std::exit(cmd_dedup(dedup_archive, dedup_dry_run));
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

    // genopack taxonomy
    // genopack taxonomy  (parent — sub-subcommands below)
    auto* tax_cmd = app.add_subcommand("taxonomy",
        "Taxonomy utilities: lookup from archive, normalize TSV, partition for distributed build");
    tax_cmd->require_subcommand(1);

    // genopack taxonomy show <archive>
    auto* tax_show = tax_cmd->add_subcommand("show",
        "Show taxonomy lineage for an accession or archive summary");
    std::string tax_archive, tax_accession;
    bool tax_json = false;
    tax_show->add_option("archive", tax_archive, "Archive directory")->required();
    tax_show->add_option("--accession", tax_accession, "Accession to look up");
    tax_show->add_flag("--json", tax_json, "Output JSON");
    tax_show->callback([&]() {
        std::exit(cmd_taxonomy(tax_archive, tax_accession, tax_json));
    });

    // genopack taxonomy normalize -i <tsv> -o <tsv> [--ncbi-taxdump <dir>]
    auto* tax_norm = tax_cmd->add_subcommand("normalize",
        "Normalize TSV taxonomy to canonical 10-rank format (d,l,k,p,c,o,f,g,s,S). "
        "GTDB prokaryotes (d__-prefixed) are normalized by rank propagation. "
        "Eukaryotes and viruses are resolved from their taxonomy string via NCBI taxdump "
        "when --ncbi-taxdump is provided.");
    std::filesystem::path tax_norm_input, tax_norm_output, tax_norm_taxdump;
    tax_norm->add_option("-i,--input", tax_norm_input,
        "Input TSV (accession TAB taxonomy TAB file_path)")->required()->check(CLI::ExistingFile);
    tax_norm->add_option("-o,--output", tax_norm_output, "Output normalized TSV")->required();
    tax_norm->add_option("--ncbi-taxdump", tax_norm_taxdump,
        "Directory containing NCBI nodes.dmp + names.dmp. "
        "Downloaded automatically if absent or older than 30 days. "
        "Enables resolution of eukaryote/virus taxonomy by matching names against the NCBI tree.");
    tax_norm->callback([&]() {
        std::ifstream fin(tax_norm_input);
        if (!fin) { spdlog::error("Cannot open: {}", tax_norm_input.string()); std::exit(1); }
        std::ofstream fout(tax_norm_output);
        if (!fout) { spdlog::error("Cannot write: {}", tax_norm_output.string()); std::exit(1); }

        // Optionally load NCBI taxdb for eukaryote/virus rows
        std::optional<genopack::NcbiTaxdb> ncbi;
        if (!tax_norm_taxdump.empty()) {
            genopack::NcbiTaxdb::ensure_fresh(tax_norm_taxdump);
            ncbi = genopack::NcbiTaxdb::load(tax_norm_taxdump);
        }

        // 10-rank canonical order
        static constexpr std::array<std::string_view, 10> kRanks = {
            "d__","l__","k__","p__","c__","o__","f__","g__","s__","S__"
        };
        auto stem_fn = [](std::string_view acc) -> std::string {
            if (acc.starts_with("RS_") || acc.starts_with("GB_")) acc = acc.substr(3);
            if (acc.starts_with("GCF_") || acc.starts_with("GCA_")) {
                auto dot = acc.rfind('.');
                if (dot != std::string_view::npos) return std::string(acc.substr(0, dot));
            }
            return std::string(acc);
        };
        // GTDB path: propagate missing ranks, derive S__ from accession stem
        auto normalize_gtdb = [&](const std::string& tax, const std::string& acc) -> std::string {
            const std::string stem = stem_fn(acc);
            std::unordered_map<std::string,std::string> rm;
            std::string_view sv(tax);
            while (!sv.empty()) {
                auto sep = sv.find(';');
                auto tok = sv.substr(0, sep);
                if (tok.size() >= 3 && tok[1]=='_' && tok[2]=='_')
                    rm.emplace(std::string(tok.substr(0,3)), std::string(tok));
                sv = (sep == std::string_view::npos) ? "" : sv.substr(sep+1);
            }
            auto prop = [&](std::string_view c, std::string_view p){
                std::string cs(c), ps(p);
                if (!rm.count(cs)||rm[cs]==cs) if(rm.count(ps)) rm[cs]=cs+rm[ps].substr(3);
            };
            prop("l__","d__"); prop("k__","l__"); prop("p__","k__"); prop("c__","p__");
            prop("o__","c__"); prop("f__","o__"); prop("g__","f__");
            auto& s=rm["s__"]; if(s.empty()||s=="s__") s="s__"+stem;
            rm["S__"]="S__"+rm["s__"].substr(3);
            std::string res; res.reserve(tax.size()+30);
            for (auto r : kRanks) { if(!res.empty()) res+=';'; res+=rm[std::string(r)]; }
            return res;
        };

        std::string header, line;
        std::getline(fin, header);
        fout << header << '\n';
        size_t n_gtdb=0, n_ncbi=0, n_unresolved=0, total=0;
        while (std::getline(fin, line)) {
            if (line.empty()) continue; ++total;
            auto t1=line.find('\t');
            if (t1==std::string::npos) { fout<<line<<'\n'; continue; }
            auto t2=line.find('\t',t1+1);
            std::string acc=line.substr(0,t1);
            std::string tax=(t2==std::string::npos)?line.substr(t1+1):line.substr(t1+1,t2-t1-1);
            std::string rest=(t2==std::string::npos)?"":line.substr(t2);

            std::string norm;
            if (tax.starts_with("d__")) {
                norm = normalize_gtdb(tax, acc);
                ++n_gtdb;
            } else if (ncbi.has_value()) {
                norm = ncbi->taxonomy_for_string(tax, acc);
                if (norm.empty()) {
                    const std::string stem = stem_fn(acc);
                    norm = "d__Unclassified;l__" + stem + ";k__" + stem + ";"
                           "p__" + stem + ";c__" + stem + ";o__" + stem + ";"
                           "f__" + stem + ";g__" + stem + ";s__" + stem + ";S__" + stem;
                    ++n_unresolved;
                } else {
                    ++n_ncbi;
                }
            } else {
                const std::string stem = stem_fn(acc);
                norm = "d__Unclassified;l__" + stem + ";k__" + stem + ";"
                       "p__" + stem + ";c__" + stem + ";o__" + stem + ";"
                       "f__" + stem + ";g__" + stem + ";s__" + stem + ";S__" + stem;
                ++n_unresolved;
            }
            fout << acc << '\t' << norm << rest << '\n';
            if (total%500'000==0) spdlog::info("  {} rows processed...", total);
        }
        fout.close();
        spdlog::info("taxonomy normalize: {} rows → {}", total, tax_norm_output.string());
        spdlog::info("  GTDB prokaryotes:    {}", n_gtdb);
        spdlog::info("  NCBI euk/virus:      {}", n_ncbi);
        spdlog::info("  unresolved (stub):   {}", n_unresolved);
        std::exit(0);
    });

    // genopack taxonomy partition -i <tsv> -n N -o <dir> [-r g|f]
    auto* tax_part = tax_cmd->add_subcommand("partition",
        "Partition TSV into N genus-balanced parts using LPT bin-packing. "
        "All genomes of the same genus land in the same part, sorted by taxonomy within each part.");
    std::filesystem::path tax_part_input, tax_part_output;
    int tax_part_n = 1;
    std::string tax_part_rank = "g";
    tax_part->add_option("-i,--input", tax_part_input,
        "Input TSV (accession TAB taxonomy TAB file_path)")->required()->check(CLI::ExistingFile);
    tax_part->add_option("-n,--parts", tax_part_n, "Number of output parts")->required();
    tax_part->add_option("-o,--output-dir", tax_part_output,
        "Output directory for part_0.tsv ... part_N-1.tsv")->required();
    tax_part->add_option("-r,--rank", tax_part_rank,
        "Rank to partition by: g=genus (default), f=family")->default_val("g");
    tax_part->callback([&]() {
        if (tax_part_n <= 0) { spdlog::error("--parts must be >= 1"); std::exit(1); }
        std::filesystem::create_directories(tax_part_output);
        const std::string rp = tax_part_rank + "__";

        // Read header + group rows by rank key
        std::ifstream fin2(tax_part_input);
        if (!fin2) { spdlog::error("Cannot open: {}", tax_part_input.string()); std::exit(1); }
        std::string hdr2, ln2;
        std::getline(fin2, hdr2);
        std::map<std::string,std::vector<std::string>> rank_rows;
        size_t total2 = 0;
        auto extract_key = [&](std::string_view tax) -> std::string {
            std::string needle = ";" + rp;
            auto pos = tax.find(needle);
            if (pos != std::string_view::npos) {
                auto s = pos+1, e = tax.find(';', s);
                return std::string(tax.substr(s, e==std::string_view::npos?std::string_view::npos:e-s));
            }
            if (tax.starts_with(rp)) {
                auto e = tax.find(';');
                return std::string(tax.substr(0, e==std::string_view::npos?std::string_view::npos:e));
            }
            return "__unknown__";
        };
        while (std::getline(fin2, ln2)) {
            if (ln2.empty()) continue; ++total2;
            std::string_view sv2(ln2);
            auto t1=sv2.find('\t');
            if (t1==std::string_view::npos) { rank_rows["__unknown__"].push_back(ln2); continue; }
            auto t2=sv2.find('\t',t1+1);
            auto tax2=(t2==std::string_view::npos)?sv2.substr(t1+1):sv2.substr(t1+1,t2-t1-1);
            rank_rows[extract_key(tax2)].push_back(std::move(ln2));
        }
        fin2.close();
        spdlog::info("taxonomy partition: {} genomes, {} {} groups",
                     total2, rank_rows.size(), rp);

        // LPT bin-packing
        std::vector<std::pair<std::string,std::vector<std::string>*>> grps;
        grps.reserve(rank_rows.size());
        for (auto& [k,v] : rank_rows) grps.push_back({k,&v});
        std::sort(grps.begin(),grps.end(),[](auto& a,auto& b){ return a.second->size()>b.second->size(); });

        std::vector<std::vector<int>> bins(tax_part_n);
        std::vector<size_t> bcnt(tax_part_n, 0);
        for (int gi=0; gi<(int)grps.size(); ++gi) {
            int t=static_cast<int>(std::min_element(bcnt.begin(),bcnt.end())-bcnt.begin());
            bins[t].push_back(gi); bcnt[t]+=grps[gi].second->size();
        }

        for (int i=0; i<tax_part_n; ++i) {
            auto op = tax_part_output / ("part_"+std::to_string(i)+".tsv");
            std::ofstream fo(op);
            if (!fo) { spdlog::error("Cannot write: {}", op.string()); std::exit(1); }
            fo << hdr2 << '\n';
            std::sort(bins[i].begin(),bins[i].end(),[&](int a,int b){ return grps[a].first<grps[b].first; });
            for (int gi : bins[i]) for (auto& row : *grps[gi].second) fo << row << '\n';
            spdlog::info("  part_{}: {} genomes, {} {} groups → {}",
                         i, bcnt[i], bins[i].size(), rp, op.string());
        }
        auto [mn,mx]=std::minmax_element(bcnt.begin(),bcnt.end());
        spdlog::info("taxonomy partition: load balance min={} max={}", *mn, *mx);
        std::exit(0);
    });

    // genopack taxonomy assign-taxids -i <normalized.tsv> -o <registry.tsv> [--acc-map <acc2taxid.tsv>]
    auto* tax_assign = tax_cmd->add_subcommand("assign-taxids",
        "Assign stable GTDB concept_ids to every unique canonical path in a normalized taxonomy TSV. "
        "Concept_ids are FNV-64(canonical_path) | (1<<63) — deterministic, no central counter needed. "
        "Produces a registry TSV (path→concept_id) consumed by workers and embedded in the final GPK.");
    std::filesystem::path tax_assign_input, tax_assign_registry, tax_assign_accmap;
    tax_assign->add_option("-i,--input", tax_assign_input,
        "Normalized TSV (accession TAB taxonomy ...)")->required()->check(CLI::ExistingFile);
    tax_assign->add_option("-o,--registry", tax_assign_registry,
        "Output registry TSV: canonical_path TAB concept_id (sorted by path)")->required();
    tax_assign->add_option("--acc-map", tax_assign_accmap,
        "Optional: output per-accession TSV: accession TAB concept_id TAB taxonomy");
    tax_assign->callback([&]() {
        using genopack::TaxonomyTree;
        using genopack::TAXID_GTDB_BIT;

        std::ifstream fin(tax_assign_input);
        if (!fin) { spdlog::error("Cannot open: {}", tax_assign_input.string()); std::exit(1); }

        // path_to_id and id_to_path for collision resolution
        std::unordered_map<std::string, uint64_t> path_to_id;
        std::unordered_map<uint64_t, std::string> id_to_path;

        auto assign = [&](const std::string& path) -> uint64_t {
            auto it = path_to_id.find(path);
            if (it != path_to_id.end()) return it->second;
            uint64_t tid = TaxonomyTree::concept_id_for_path(path);
            while (true) {
                auto ex = id_to_path.find(tid);
                if (ex == id_to_path.end()) break;
                if (ex->second == path) break;
                ++tid;
                tid |= TAXID_GTDB_BIT;
                if ((tid & ~TAXID_GTDB_BIT) == 0) tid |= 1;
            }
            path_to_id[path] = tid;
            id_to_path[tid]  = path;
            return tid;
        };

        // Per-accession results
        struct AccEntry { std::string accession; uint64_t concept_id; std::string taxonomy; };
        std::vector<AccEntry> acc_entries;

        std::string header, line;
        std::getline(fin, header);
        size_t total = 0, n_collisions = 0;

        while (std::getline(fin, line)) {
            if (line.empty()) continue;
            ++total;
            auto t1 = line.find('\t');
            if (t1 == std::string::npos) continue;
            auto t2 = line.find('\t', t1 + 1);
            std::string acc = line.substr(0, t1);
            std::string tax = (t2 == std::string::npos)
                            ? line.substr(t1 + 1)
                            : line.substr(t1 + 1, t2 - t1 - 1);

            // Build cumulative paths at each rank and assign concept_ids
            // The deepest rank path is the genome's own concept_id
            std::string cur_path;
            uint64_t deepest_id = 0;
            std::string_view sv(tax);
            while (!sv.empty()) {
                auto sep = sv.find(';');
                auto tok = (sep == std::string_view::npos) ? sv : sv.substr(0, sep);
                sv = (sep == std::string_view::npos) ? "" : sv.substr(sep + 1);
                if (tok.size() < 3 || tok[1] != '_' || tok[2] != '_') continue;
                // Skip rank-prefix-only tokens (e.g. "s__" with empty name)
                if (tok.size() == 3) continue;

                cur_path = cur_path.empty()
                         ? std::string(tok)
                         : cur_path + ";" + std::string(tok);

                uint64_t id = assign(cur_path);
                deepest_id = id;
            }

            if (!tax_assign_accmap.empty() && deepest_id != 0)
                acc_entries.push_back({acc, deepest_id, tax});

            if (total % 500'000 == 0)
                spdlog::info("  {} rows processed, {} unique paths...", total, path_to_id.size());
        }

        // Check for collisions (paths that needed probing)
        for (const auto& [path, id] : path_to_id) {
            if (id != TaxonomyTree::concept_id_for_path(path)) ++n_collisions;
        }

        // Write registry TSV sorted by path
        {
            std::vector<std::pair<std::string_view, uint64_t>> rows;
            rows.reserve(path_to_id.size());
            for (const auto& [p, id] : path_to_id) rows.emplace_back(p, id);
            std::sort(rows.begin(), rows.end(),
                      [](const auto& a, const auto& b) { return a.first < b.first; });

            std::ofstream fo(tax_assign_registry);
            if (!fo) { spdlog::error("Cannot write: {}", tax_assign_registry.string()); std::exit(1); }
            fo << "canonical_path\tconcept_id\n";
            for (const auto& [p, id] : rows)
                fo << p << '\t' << id << '\n';
        }

        // Write optional per-accession map
        if (!tax_assign_accmap.empty()) {
            std::sort(acc_entries.begin(), acc_entries.end(),
                      [](const AccEntry& a, const AccEntry& b) {
                          return a.accession < b.accession;
                      });
            std::ofstream fa(tax_assign_accmap);
            if (!fa) { spdlog::error("Cannot write: {}", tax_assign_accmap.string()); std::exit(1); }
            fa << "accession\tconcept_id\ttaxonomy\n";
            for (const auto& e : acc_entries)
                fa << e.accession << '\t' << e.concept_id << '\t' << e.taxonomy << '\n';
        }

        spdlog::info("taxonomy assign-taxids: {} rows, {} unique paths, {} hash collisions resolved",
                     total, path_to_id.size(), n_collisions);
        spdlog::info("  registry:  {}", tax_assign_registry.string());
        if (!tax_assign_accmap.empty())
            spdlog::info("  acc-map:   {}", tax_assign_accmap.string());
        std::exit(0);
    });

    // genopack taxonomy diff
    // Diffs a current genome set against a new GTDB release taxonomy file.
    // Accepts either the original input TSV (accession/taxonomy/file) or the
    // genopack meta TSV (accession/genome_id/taxonomy/representative) as --current.
    // GTDB r232 uses 7-rank format (no l__/k__); comparison is done at s__ level.
    //
    // Output categories:
    //   unchanged        — same accession, same s__ in new release
    //   taxonomy_changed — same accession, different s__ → update string, no GTDB-Tk
    //   needs_gtdbtk     — custom MAG (non-GCF/GCA) whose assigned species no longer
    //                      exists in the new release → must re-run GTDB-Tk
    //   custom_stable    — custom MAG whose species still exists → no action needed
    //   dropped_from_gtdb — GCF/GCA accession absent from new release (GTDB removed it)
    //   new_in_release   — in new GTDB release but not in current set → needs FASTA + add
    auto* tax_diff = tax_cmd->add_subcommand("diff",
        "Diff current genome set against a new GTDB release. "
        "Classifies each genome as unchanged, taxonomy_changed, needs_gtdbtk, "
        "dropped_from_gtdb, or new_in_release. "
        "Accepts either the full input TSV (accession/taxonomy/file) or the "
        "genopack meta TSV (accession/genome_id/taxonomy/representative).");
    std::filesystem::path tax_diff_current, tax_diff_outdir;
    std::vector<std::string> tax_diff_gtdb_files;
    bool tax_diff_write_unchanged = false;
    tax_diff->add_option("--current", tax_diff_current,
        "Current genome set TSV. Accepts the original input TSV "
        "(accession TAB taxonomy TAB file) or the genopack meta TSV "
        "(accession TAB genome_id TAB taxonomy TAB representative)."
        )->required()->check(CLI::ExistingFile);
    tax_diff->add_option("--gtdb", tax_diff_gtdb_files,
        "New GTDB release taxonomy TSV(s). Pass bac120_taxonomy.tsv and "
        "ar53_taxonomy.tsv separately, or a pre-merged file. Repeatable."
        )->required()->expected(1, 10);
    tax_diff->add_option("-o,--output-dir", tax_diff_outdir,
        "Output directory for per-category TSV files and summary.txt")->required();
    tax_diff->add_flag("--write-unchanged", tax_diff_write_unchanged,
        "Also write unchanged.tsv (can be very large for full GTDB sets)");
    tax_diff->callback([&]() {
        namespace fs = std::filesystem;

        // Strip RS_/GB_ prefix used by GTDB on NCBI accessions.
        auto strip_prefix = [](std::string_view acc) -> std::string {
            if (acc.starts_with("RS_") || acc.starts_with("GB_")) acc = acc.substr(3);
            return std::string(acc);
        };

        // Extract a rank token (e.g. "s__Foo bar") from a semicolon-delimited
        // taxonomy string. The prefix must appear at the start or after a semicolon.
        auto extract_rank = [](std::string_view tax, std::string_view prefix) -> std::string {
            auto pos = tax.find(prefix);
            if (pos == std::string_view::npos) return {};
            if (pos > 0 && tax[pos - 1] != ';') return {};
            auto end = tax.find(';', pos);
            return std::string(end == std::string_view::npos
                               ? tax.substr(pos) : tax.substr(pos, end - pos));
        };

        // True for standard NCBI accessions (GCF_/GCA_) — these are dropped by GTDB,
        // not custom MAGs that need GTDB-Tk reclassification.
        auto is_ncbi_acc = [](std::string_view acc) -> bool {
            return acc.starts_with("GCF_") || acc.starts_with("GCA_");
        };

        fs::create_directories(tax_diff_outdir);

        // ── Step 1: load new GTDB taxonomy ───────────────────────────────────
        // Map: stripped_accession → raw_taxonomy (7-rank, no l__/k__).
        // Also collect all unique s__ tokens for the "species still exists?" check.
        spdlog::info("taxonomy diff: loading GTDB taxonomy file(s)...");
        std::unordered_map<std::string, std::string> r_new;
        std::unordered_set<std::string> new_species;
        r_new.reserve(600'000);

        for (const auto& gtdb_file : tax_diff_gtdb_files) {
            std::ifstream fg(gtdb_file);
            if (!fg) { spdlog::error("Cannot open GTDB file: {}", gtdb_file); std::exit(1); }
            std::string ln;
            size_t n_file = 0;
            while (std::getline(fg, ln)) {
                if (ln.empty()) continue;
                auto tab = ln.find('\t');
                if (tab == std::string::npos) continue;
                std::string_view tax_sv(ln.data() + tab + 1, ln.size() - tab - 1);
                // Skip header lines and non-GTDB entries.
                if (!tax_sv.starts_with("d__")) continue;
                std::string acc = strip_prefix(std::string_view(ln.data(), tab));
                std::string tax(tax_sv);
                auto sp = extract_rank(tax_sv, "s__");
                if (sp.size() > 3) new_species.insert(sp);
                r_new.emplace(std::move(acc), std::move(tax));
                ++n_file;
            }
            spdlog::info("  {} genomes from {}", n_file, gtdb_file);
        }
        spdlog::info("taxonomy diff: {} genomes in new release, {} unique species",
                     r_new.size(), new_species.size());

        // ── Step 2: detect column layout of current TSV ───────────────────
        // meta TSV:  accession | genome_id | taxonomy | representative
        // input TSV: accession | taxonomy  | file
        // We detect by checking whether the header contains "genome_id".
        std::ifstream fc(tax_diff_current);
        if (!fc) { spdlog::error("Cannot open: {}", tax_diff_current.string()); std::exit(1); }
        std::string hdr;
        std::getline(fc, hdr);
        const bool is_meta = (hdr.find("genome_id") != std::string::npos);
        const int col_acc  = 0;
        const int col_tax  = is_meta ? 2 : 1;
        const int col_file = is_meta ? -1 : 2;
        spdlog::info("taxonomy diff: format={} (tax_col={}, file_col={})",
                     is_meta ? "meta-tsv" : "input-tsv", col_tax, col_file);

        // ── Step 3: open output files ─────────────────────────────────────
        auto open_out = [&](const char* name) {
            std::ofstream f(tax_diff_outdir / name);
            if (!f) { spdlog::error("Cannot write: {}/{}", tax_diff_outdir.string(), name); std::exit(1); }
            return f;
        };
        std::ofstream f_changed  = open_out("taxonomy_changed.tsv");
        std::ofstream f_gtdbtk   = open_out("needs_gtdbtk.tsv");
        std::ofstream f_dropped  = open_out("dropped_from_gtdb.tsv");
        std::optional<std::ofstream> f_unch;
        if (tax_diff_write_unchanged) f_unch = open_out("unchanged.tsv");

        f_changed  << "accession\told_taxonomy\tnew_taxonomy\n";
        f_gtdbtk   << "accession\tcurrent_taxonomy" << (col_file >= 0 ? "\tfile\n" : "\n");
        f_dropped  << "accession\tcurrent_taxonomy\n";
        if (f_unch) *f_unch << "accession\ttaxonomy\n";

        // ── Step 4: classify each genome in current set ───────────────────
        size_t n_total=0, n_unchanged=0, n_changed=0,
               n_gtdbtk=0, n_custom_stable=0, n_dropped=0;
        // Only track accessions that matched in r_new (for new_in_release detection).
        std::unordered_set<std::string> found_in_new;
        found_in_new.reserve(r_new.size());

        std::string line;
        while (std::getline(fc, line)) {
            if (line.empty()) continue;
            ++n_total;

            // Split tab-delimited columns into string_views into `line`.
            std::array<std::string_view, 8> cols{};
            int n_cols = 0;
            std::string_view sv(line);
            while (!sv.empty() && n_cols < 8) {
                auto t = sv.find('\t');
                cols[n_cols++] = sv.substr(0, t);
                if (t == std::string_view::npos) break;
                sv = sv.substr(t + 1);
            }

            std::string_view acc_sv  = (col_acc  < n_cols) ? cols[col_acc]  : std::string_view{};
            std::string_view tax_sv  = (col_tax  < n_cols) ? cols[col_tax]  : std::string_view{};
            std::string_view file_sv = (col_file >= 0 && col_file < n_cols) ? cols[col_file] : std::string_view{};

            const std::string stripped = strip_prefix(acc_sv);
            auto it = r_new.find(stripped);

            if (it != r_new.end()) {
                // Accession exists in new release — compare at s__ level.
                found_in_new.insert(stripped);
                const std::string& new_raw = it->second;
                std::string cur_sp = extract_rank(tax_sv, "s__");
                std::string new_sp = extract_rank(new_raw, "s__");
                if (cur_sp == new_sp) {
                    ++n_unchanged;
                    if (f_unch) *f_unch << acc_sv << '\t' << tax_sv << '\n';
                } else {
                    ++n_changed;
                    f_changed << acc_sv << '\t' << tax_sv << '\t' << new_raw << '\n';
                }
            } else {
                // Not in new release.
                if (is_ncbi_acc(stripped)) {
                    // Standard NCBI accession that GTDB dropped.
                    ++n_dropped;
                    f_dropped << acc_sv << '\t' << tax_sv << '\n';
                } else {
                    // Custom MAG — check whether its assigned species still exists.
                    std::string cur_sp = extract_rank(tax_sv, "s__");
                    if (!cur_sp.empty() && cur_sp.size() > 3 && new_species.count(cur_sp)) {
                        ++n_custom_stable;  // species exists → no action needed
                    } else {
                        ++n_gtdbtk;
                        f_gtdbtk << acc_sv << '\t' << tax_sv;
                        if (col_file >= 0) f_gtdbtk << '\t' << file_sv;
                        f_gtdbtk << '\n';
                    }
                }
            }
            if (n_total % 500'000 == 0)
                spdlog::info("  {} genomes classified...", n_total);
        }

        // ── Step 5: genomes new in the release, absent from current set ───
        std::ofstream f_new = open_out("new_in_release.tsv");
        f_new << "accession\ttaxonomy\n";
        size_t n_new = 0;
        for (const auto& [acc, tax] : r_new) {
            if (!found_in_new.count(acc)) {
                f_new << acc << '\t' << tax << '\n';
                ++n_new;
            }
        }

        // ── Step 6: summary ───────────────────────────────────────────────
        std::ofstream f_sum = open_out("summary.txt");
        f_sum << "GTDB taxonomy diff summary\n"
              << "==========================\n"
              << "Current set:               " << n_total          << "\n"
              << "  unchanged:               " << n_unchanged       << "\n"
              << "  taxonomy_changed:        " << n_changed         << "  → taxonomy_changed.tsv\n"
              << "  custom_stable (no-op):   " << n_custom_stable   << "\n"
              << "  needs_gtdbtk:            " << n_gtdbtk          << "  → needs_gtdbtk.tsv\n"
              << "  dropped_from_gtdb:       " << n_dropped         << "  → dropped_from_gtdb.tsv\n"
              << "New in release:            " << n_new             << "  → new_in_release.tsv\n";

        // Close all streams before std::exit — exit() skips local destructors
        // so buffered data would be lost without explicit close().
        f_changed.close(); f_gtdbtk.close(); f_dropped.close();
        f_new.close();
        if (f_unch) f_unch->close();
        f_sum.close();

        spdlog::info("taxonomy diff: {} genomes classified", n_total);
        spdlog::info("  unchanged:              {}", n_unchanged);
        spdlog::info("  taxonomy_changed:       {} → taxonomy_changed.tsv", n_changed);
        spdlog::info("  custom_stable (no-op):  {}", n_custom_stable);
        spdlog::info("  needs_gtdbtk:           {} → needs_gtdbtk.tsv", n_gtdbtk);
        spdlog::info("  dropped_from_gtdb:      {} → dropped_from_gtdb.tsv", n_dropped);
        spdlog::info("  new_in_release:         {} → new_in_release.tsv", n_new);
        spdlog::info("  summary:                {}/summary.txt", tax_diff_outdir.string());
        std::exit(0);
    });

    // genopack taxonomy patch
    auto* tax_patch = tax_cmd->add_subcommand("patch",
        "Patch taxonomy strings in a .gpk archive and/or flat input TSV. "
        "Accepts three patch formats: simple 2-col (accession TAB new_taxonomy), "
        "taxonomy_changed.tsv from 'taxonomy diff' (3-col, uses col 2), "
        "or GTDB-Tk classify summary (--gtdbtk). "
        "Normalizes 7-rank GTDB strings to 10-rank automatically. "
        "Archive: rewrites TAXN + TXDB sections in-place (multipart-aware). "
        "TSV: atomic rename via temp file.");
    std::filesystem::path tax_patch_archive, tax_patch_file, tax_patch_tsv, tax_patch_tsv_out;
    bool tax_patch_gtdbtk  = false;
    bool tax_patch_normalize = true;
    tax_patch->add_option("--patch", tax_patch_file,
        "Patch TSV. Formats: (1) accession TAB new_taxonomy, "
        "(2) taxonomy_changed.tsv (accession/old_taxonomy/new_taxonomy — uses col 2), "
        "(3) GTDB-Tk summary with --gtdbtk flag."
        )->required()->check(CLI::ExistingFile);
    tax_patch->add_option("--archive", tax_patch_archive,
        "Archive to patch: single .gpk file or multipart directory. "
        "Rewrites TAXN and TXDB sections; all other sections preserved.");
    tax_patch->add_option("--tsv", tax_patch_tsv,
        "Flat input TSV to patch (accession/taxonomy[/file]). "
        "Written atomically via temp file."
        )->check(CLI::ExistingFile);
    tax_patch->add_option("--tsv-out", tax_patch_tsv_out,
        "Output path for patched TSV (default: overwrite --tsv in-place).");
    tax_patch->add_flag("--gtdbtk", tax_patch_gtdbtk,
        "Parse GTDB-Tk classify summary format "
        "(user_genome/classification/... — skips N/A entries).");
    tax_patch->add_flag("!--no-normalize", tax_patch_normalize,
        "Normalize 7-rank GTDB strings to 10-rank (default: on). "
        "Pass --no-normalize to skip if patch file is already normalized.");
    tax_patch->callback([&]() {
        if (tax_patch_archive.empty() && tax_patch_tsv.empty()) {
            spdlog::error("taxonomy patch: specify at least --archive or --tsv");
            std::exit(1);
        }
        std::exit(cmd_taxn_patch(tax_patch_archive, tax_patch_file,
                                 tax_patch_tsv, tax_patch_tsv_out,
                                 tax_patch_gtdbtk, tax_patch_normalize));
    });

    // genopack taxdump
    auto* taxdump_cmd = app.add_subcommand("taxdump", "Export taxonomy as NCBI taxdump or high-performance columnar binary");
    std::string taxdump_archive, taxdump_format = "columnar", taxdump_output;
    taxdump_cmd->add_option("archive", taxdump_archive, "Archive path")->required();
    taxdump_cmd->add_option("-f,--format", taxdump_format,
        "Output format: 'taxdump' (NCBI names/nodes/acc2taxid.dmp) or 'columnar' (binary + TSV, default)");
    taxdump_cmd->add_option("-o,--output", taxdump_output, "Output directory")->required();
    taxdump_cmd->callback([&]() {
        std::exit(cmd_taxdump(taxdump_archive, taxdump_format, taxdump_output));
    });

    // genopack reindex
    auto* reindex_cmd = app.add_subcommand("reindex", "Append or rebuild GIDX/TXDB/CIDX/SKCH sections on an existing archive");
    std::string reindex_archive;
    bool reindex_force = false;
    bool reindex_txdb = false;
    bool reindex_skch = false;
    std::string reindex_cidx_tsv;
    int reindex_cidx_threads = 8;
    int reindex_skch_threads = 8;
    int reindex_skch_kmer = -1, reindex_skch_size = -1, reindex_skch_syncmer = -1;
    std::string reindex_skch_kmers_str;  // comma-separated, e.g. "16,21,31"
    reindex_cmd->add_option("archive", reindex_archive, "Path to .gpk archive")->required();
    reindex_cmd->add_flag("--force", reindex_force, "Rebuild indexes even if already present");
    bool reindex_no_gidx = false;
    reindex_cmd->add_flag("--no-gidx", reindex_no_gidx, "Skip GIDX build");
    bool reindex_no_gstx = false;
    reindex_cmd->add_flag("--no-gstx", reindex_no_gstx, "Skip GSTX genus-stats index build");
    bool reindex_no_qual = false;
    reindex_cmd->add_flag("--no-qual", reindex_no_qual, "Skip QUAL per-genome quality score build");
    bool reindex_gcov = false;
    reindex_cmd->add_flag("--gcov", reindex_gcov, "Rebuild per-genus covariance (GCOV/FCOV/FMHR) — full shard rescan");
    reindex_cmd->add_flag("--txdb", reindex_txdb, "Build taxonomy tree (TXDB) from TAXN lineage strings");
    reindex_cmd->add_option("--cidx", reindex_cidx_tsv, "Build contig accession index (CIDX) from build TSV (accession<TAB>taxonomy<TAB>file_path)");
    reindex_cmd->add_option("--cidx-threads", reindex_cidx_threads, "Threads for parallel FASTA decompression (default: 8)");
    reindex_cmd->add_flag("--skch", reindex_skch, "Compute OPH sketches for genomes missing from existing SKCH sections");
    reindex_cmd->add_option("--skch-threads", reindex_skch_threads, "Threads for parallel sketch computation (default: 8)");
    reindex_cmd->add_option("--sketch-kmer", reindex_skch_kmer, "OPH k-mer size for single-k SKCH section (default: inherit from existing or 16)");
    reindex_cmd->add_option("--sketch-kmers", reindex_skch_kmers_str,
        "Comma-separated k-mer sizes for multi-k v2 SKCH section (e.g. 16,21,31). "
        "Takes precedence over --sketch-kmer when more than one value is given.");
    reindex_cmd->add_option("--sketch-size", reindex_skch_size, "OPH sketch size for new SKCH section (default: inherit from existing or 10000)");
    reindex_cmd->add_option("--sketch-syncmer", reindex_skch_syncmer, "Syncmer s for new SKCH section (0=disabled, default: inherit or 0)");
    reindex_cmd->callback([&]() {
        // Parse --sketch-kmers comma list into a vector<int>.
        std::vector<int> reindex_skch_kmers;
        if (!reindex_skch_kmers_str.empty()) {
            std::istringstream ss(reindex_skch_kmers_str);
            std::string tok;
            while (std::getline(ss, tok, ',')) {
                int k = std::stoi(tok);
                if (k > 0) reindex_skch_kmers.push_back(k);
            }
        }
        std::exit(cmd_reindex(reindex_archive, reindex_force, reindex_txdb,
                              reindex_cidx_tsv, reindex_cidx_threads,
                              reindex_skch, reindex_skch_threads,
                              reindex_skch_kmer, reindex_skch_size, reindex_skch_syncmer,
                              std::move(reindex_skch_kmers),
                              reindex_no_gidx, 16, !reindex_no_gstx, !reindex_no_qual,
                              reindex_gcov));
    });

    // genopack repack
    auto* repack_cmd = app.add_subcommand("repack",
        "Re-shard an archive by taxonomy (genus by default) for fast per-taxon access on NFS. "
        "One sequential pass through the source archive; genome IDs are preserved.");
    std::string repack_input, repack_output, repack_rank = "g";
    int repack_level = 6, repack_threads = 1, repack_max_mem = 32;
    bool repack_verbose = false;
    repack_cmd->add_option("input",  repack_input,  "Source .gpk archive")->required();
    repack_cmd->add_option("output", repack_output, "Output .gpk archive")->required();
    repack_cmd->add_option("-z,--zstd-level",  repack_level,   "zstd compression level (1-22)")->default_val(6);
    repack_cmd->add_option("-t,--threads",     repack_threads,
        "Threads for parallel decompression within each shard")->default_val(1);
    repack_cmd->add_option("-m,--max-memory",  repack_max_mem,
        "Max total bucket memory in GB before forcing a flush of all buckets (default: 32)")->default_val(32);
    repack_cmd->add_option("--taxonomy-rank",  repack_rank,
        "Rank for shard grouping: g=genus (default), f=family")->default_val("g");
    repack_cmd->add_flag("-v,--verbose", repack_verbose, "Log progress every 2000 shards");
    repack_cmd->callback([&]() {
        std::exit(cmd_repack(repack_input, repack_output, repack_level, repack_rank,
                             repack_threads, repack_max_mem, repack_verbose));
    });

    // genopack coordinator
    auto* coord_cmd = app.add_subcommand("coordinator",
        "Start an NFS manifest coordinator. Workers connect via 'genopack build --coordinator'. "
        "Creates the output .gpk file, allocates write offsets for workers, writes unified TOC.");
    std::string coord_output;
    std::filesystem::path coord_ntdb;
    std::string coord_nfs_dir;
    int coord_workers = 0;
    coord_cmd->add_option("-o,--output", coord_output, "Output .gpk file path")->required();
    coord_cmd->add_option("--workers", coord_workers,
        "Expected number of workers")->required();
    coord_cmd->add_option("--nfs-dir", coord_nfs_dir,
        "NFS manifest directory for coordination. "
        "Workers use '--coordinator this_dir:/output.gpk'.")->required();
    coord_cmd->add_option("--ntdb", coord_ntdb,
        "Directory containing NCBI nodes.dmp + names.dmp. "
        "If provided, the coordinator embeds the full NCBI tree as a NTDB section "
        "in the final archive before writing the TOC.");
    coord_cmd->callback([&]() {
        genopack::CoordinatorServer srv;
        srv.run_nfs(coord_nfs_dir, coord_output, coord_workers, coord_ntdb,
                    [](size_t n) { spdlog::info("coordinator-nfs: {} sections collected", n); });
        std::exit(0);
    });

    // genopack verify
    auto* verify_cmd = app.add_subcommand("verify",
        "Verify XXH128 checksums of all shards and TOC. Exits 0 if all pass.");
    std::string verify_path;
    bool verify_verbose = false;
    verify_cmd->add_option("archive", verify_path, "Path to .gpk archive or part directory")->required();
    verify_cmd->add_flag("-v,--verbose", verify_verbose, "Print OK lines for every shard");
    bool verify_strict = false;
    verify_cmd->add_flag("--strict", verify_strict,
        "Treat sections without a content checksum (index/metadata or >512MB) as failures");
    verify_cmd->callback([&]() {
        try {
            std::exit(cmd_verify(verify_path, verify_verbose, verify_strict));
        } catch (const std::exception& e) {
            spdlog::error("verify: {}", e.what());
            std::exit(1);
        }
    });

    // genopack check
    auto* check_cmd = app.add_subcommand("check",
        "Compute per-genome quality scores (completeness, contamination) and write QUAL section.");
    std::string check_pack;
    std::string check_genomes;
    std::string check_output;
    int check_threads = 8;
    int check_min_genus_size = 3;
    float check_leakage_threshold = 0.05f;
    bool check_recompute = false;
    bool check_scan_all = false;
    check_cmd->add_option("pack", check_pack, "Path to .gpk archive or directory of parts")->required();
    check_cmd->add_option("-g,--genomes", check_genomes, "Optional accession list (one per line); default: all in pack");
    check_cmd->add_option("-o,--output", check_output, "TSV output path for quality table");
    check_cmd->add_option("-t,--threads", check_threads, "Threads (default: 8)");
    check_cmd->add_option("--min-genus-size", check_min_genus_size, "Min genus members for saturated tier (default: 3)");
    check_cmd->add_option("--leakage-threshold", check_leakage_threshold, "Containment leakage threshold (default: 0.05)");
    check_cmd->add_flag("--recompute", check_recompute, "Ignore existing QUAL section and force full rescan");
    check_cmd->add_flag("--scan-all", check_scan_all, "Force every genome through FASTA-level analysis (intrinsic completeness for all, not just flagged)");
    std::string check_markers;
    check_cmd->add_option("--markers", check_markers,
        "Path to markers .mrk DB; enables marker-based completeness/redundancy scoring");
    check_cmd->callback([&]() {
        std::exit(genopack::check::cmd_check(
            check_pack,
            check_genomes.empty() ? std::filesystem::path{} : std::filesystem::path{check_genomes},
            check_threads,
            check_min_genus_size,
            check_leakage_threshold,
            check_output.empty() ? std::filesystem::path{} : std::filesystem::path{check_output},
            check_recompute,
            check_markers.empty() ? std::filesystem::path{} : std::filesystem::path{check_markers},
            check_scan_all));
    });

    // genopack ingest — attach external quality (CheckM2 / anvi'o) as SEC_XQAL
    auto* ingest_cmd = app.add_subcommand("ingest",
        "Ingest external quality (CheckM2/anvi'o) into the archive as provenance-carrying XQAL columns.");
    std::string ingest_pack, ingest_checkm2, ingest_anvio;
    ingest_cmd->add_option("pack", ingest_pack, "Path to .gpk archive or directory of parts")->required();
    ingest_cmd->add_option("--checkm2", ingest_checkm2, "CheckM2 quality_report.tsv (Name/Completeness/Contamination)");
    ingest_cmd->add_option("--anvio", ingest_anvio, "anvi'o completeness TSV (bin name/% completion/% redundancy)");
    ingest_cmd->callback([&]() {
        std::exit(genopack::cmd_ingest(
            ingest_pack,
            ingest_checkm2.empty() ? std::filesystem::path{} : std::filesystem::path{ingest_checkm2},
            ingest_anvio.empty()   ? std::filesystem::path{} : std::filesystem::path{ingest_anvio}));
    });

    // genopack report — unified per-genome quality via a named fusion profile
    auto* report_cmd = app.add_subcommand("report",
        "Emit a unified per-genome quality table resolved through a named profile. Each axis "
        "is sourced from exactly one column (built-in rule or stored policy) and the report "
        "carries that column's provenance (tool/method) — fusion is explicit, never silent.");
    std::string report_pack, report_profile = "best", report_output;
    bool report_list = false;
    report_cmd->add_option("pack", report_pack, "Path to .gpk archive or directory of parts")->required();
    report_cmd->add_option("-p,--profile", report_profile,
        "Profile name: built-in {intrinsic,external,calibrated,best} or a stored profile (default: best)");
    report_cmd->add_option("-o,--output", report_output, "TSV output path (default: stdout)");
    report_cmd->add_flag("--list", report_list,
        "List built-in + stored profiles and the available provenance columns, then exit");
    report_cmd->callback([&]() {
        std::exit(genopack::cmd_report(
            report_pack, report_profile,
            report_output.empty() ? std::filesystem::path{} : std::filesystem::path{report_output},
            report_list));
    });

    // genopack profile — author / list stored SEC_PROF fusion policies
    auto* profile_cmd = app.add_subcommand("profile",
        "Manage named reporting/fusion profiles stored in the archive (SEC_PROF).");
    profile_cmd->require_subcommand(1);
    auto* profile_add = profile_cmd->add_subcommand("add",
        "Author a named profile pinning each axis to an exact column identity present in the archive.");
    std::string prof_pack, prof_name, prof_desc;
    std::vector<std::string> prof_selects;
    profile_add->add_option("pack", prof_pack, "Path to .gpk archive or directory of parts")->required();
    profile_add->add_option("--name", prof_name, "Profile name (must not be a built-in name)")->required();
    profile_add->add_option("--description", prof_desc, "Human-facing description (cosmetic; excluded from policy hash)");
    profile_add->add_option("-s,--select", prof_selects,
        "Axis selector axis=tool:method[@cal][/tool:method[@cal]] (repeatable); the optional /… is a fallback")
        ->required();
    profile_add->callback([&]() {
        std::exit(genopack::cmd_profile_add(prof_pack, prof_name, prof_desc, prof_selects));
    });
    auto* profile_list = profile_cmd->add_subcommand("list", "List stored profiles and available columns.");
    std::string prof_list_pack;
    profile_list->add_option("pack", prof_list_pack, "Path to .gpk archive or directory of parts")->required();
    profile_list->callback([&]() {
        std::exit(genopack::cmd_profile_list(prof_list_pack));
    });

    // genopack qcontig — dump the per-contig quality overlay (SEC_QCONTIG)
    auto* qcontig_cmd = app.add_subcommand("qcontig",
        "Dump the per-contig quality overlay: one row per (genome, contig) with offset/length/"
        "TNF/leakage, so you can see which contigs drive a genome's contamination.");
    std::string qcontig_pack, qcontig_acc, qcontig_out;
    float qcontig_min_outlier = 0.0f, qcontig_min_leak = 0.0f;
    float qcontig_min_foreign = 0.0f, qcontig_min_lr = 0.0f;
    qcontig_cmd->add_option("pack", qcontig_pack, "Path to .gpk archive or directory of parts")->required();
    qcontig_cmd->add_option("-a,--accession", qcontig_acc, "Restrict to one genome accession (default: all)");
    qcontig_cmd->add_option("-o,--output", qcontig_out, "TSV output path (default: stdout)");
    qcontig_cmd->add_option("--min-outlier", qcontig_min_outlier,
        "TNF channel (LONG contigs): emit contigs whose GCOV T² or SPE percentile >= this (e.g. 0.99)");
    qcontig_cmd->add_option("--min-foreign", qcontig_min_foreign,
        "Protein-aamer channel (SMALL contigs): emit contigs whose foreign_fraction >= this (e.g. 0.30)");
    qcontig_cmd->add_option("--min-lr", qcontig_min_lr,
        "Protein-aamer channel: also require foreign log-LR >= this (default 0; pair with --min-foreign, e.g. 3.0)");
    qcontig_cmd->add_option("--min-leakage", qcontig_min_leak,
        "Also flag contigs with containment-leakage score >= this");
    qcontig_cmd->callback([&]() {
        std::exit(genopack::cmd_qcontig(qcontig_pack, qcontig_acc,
            qcontig_out.empty() ? std::filesystem::path{} : std::filesystem::path{qcontig_out},
            qcontig_min_outlier, qcontig_min_leak, qcontig_min_foreign, qcontig_min_lr));
    });

    // genopack decontaminate
    auto* decon_cmd = app.add_subcommand("decontaminate",
        "Iteratively remove contaminated genomes (per-contig CCO above a threshold), rebuilding "
        "genus/family models from the survivors each round so the DB and its consensus stay clean.");
    std::string decon_archive;
    float decon_max_fmh_z = 5.0f;    // robust SDs above the genus FMH baseline
    float decon_min_delta = 0.02f;   // absolute fmh_minority above baseline (guards tight genera)
    float decon_cco_abs   = 1.5f;    // CCO % absolute threshold (T²∪SPE outlier bp); clean genomes ~0.7%
    int   decon_iters     = 5;
    int   decon_threads   = 8;
    bool  decon_dry       = false;
    decon_cmd->add_option("archive", decon_archive, "Path to .gpk archive")->required();
    decon_cmd->add_option("--max-fmh-z", decon_max_fmh_z,
        "Remove genomes whose per-genus FMH-minority z (delta-from-genus-baseline, robust median+MAD) "
        "exceeds this — the genus-band signal, fragmentation-robust (default: 5.0)");
    decon_cmd->add_option("--min-delta", decon_min_delta,
        "Also require fmh_minority to exceed the genus baseline by this absolute amount (guards "
        "tight-baseline genera from z blow-up; default: 0.02)");
    decon_cmd->add_option("--max-cco", decon_cco_abs,
        "Remove genomes whose CCO contamination %% (per-contig T²∪SPE TNF-outlier bp vs the genus "
        "GCOV covariance) exceeds this — the distant-band signal (family/order/phylum), reference-"
        "internal so native rRNA/plasmids are not flagged (the genus covariance carries them). "
        "Clean genomes cap ~0.7%%; absolute (a per-genus z poisons when a genus is majority-"
        "contaminated). Needs contigs ≥20kb (default: 1.5)");
    decon_cmd->add_option("--max-iters", decon_iters, "Max decontamination rounds (default: 5)");
    decon_cmd->add_option("-t,--threads", decon_threads, "Threads (default: 8)");
    decon_cmd->add_flag("--dry-run", decon_dry,
        "Report what would be removed (with fmh_z/cco per axis) without modifying the archive");
    decon_cmd->callback([&]() {
        std::exit(cmd_decontaminate(decon_archive, decon_max_fmh_z, decon_min_delta,
                                    decon_cco_abs, decon_iters, decon_threads, decon_dry));
    });

    // genopack gcov
    auto* gcov_cmd = app.add_subcommand("gcov",
        "Build (or rebuild) the GCOV per-genus covariance section in an existing .gpk archive.");
    std::string gcov_pack;
    int gcov_threads = 8;
    gcov_cmd->add_option("archive", gcov_pack, "Path to .gpk archive")->required();
    gcov_cmd->add_option("-t,--threads", gcov_threads, "Parallel shard readers (default: 8)");
    gcov_cmd->callback([&]() {
        std::exit(rebuild_gcov_fmhr(std::filesystem::path(gcov_pack), gcov_threads));
    });

    // genopack fcore — build per-family prevalence cores (genus-core fallback)
    auto* fcore_cmd = app.add_subcommand("fcore",
        "Build SEC_FCORE per-family prevalence cores: the intrinsic-completeness fallback for "
        "novel/sparse genera whose genus has no (or too thin a) CORE. Reuses the CORE params so "
        "family-core coverage is comparable to genus-core coverage; check records both as distinct columns.");
    std::string fcore_pack, fcore_members;
    int   fcore_threads = 8;
    float fcore_theta   = 0.0f;
    unsigned fcore_min_members = 2;
    fcore_cmd->add_option("archive", fcore_pack, "Path to .gpk archive or directory of parts")->required();
    fcore_cmd->add_option("-t,--threads", fcore_threads, "Parallel shard readers (default: 8)");
    fcore_cmd->add_option("--theta", fcore_theta,
        "Prevalence threshold override (default: inherit the CORE section's theta)");
    fcore_cmd->add_option("--min-members", fcore_min_members,
        "Skip families with fewer than this many genomes (default: 2)");
    fcore_cmd->add_option("--members", fcore_members,
        "Reference accession list (one per line); only these genomes build the cores "
        "(excludes co-resident query/test genomes). Default: all live genomes");
    fcore_cmd->callback([&]() {
        std::exit(genopack::cmd_fcore(std::filesystem::path(fcore_pack),
                                      fcore_threads, fcore_theta, fcore_min_members,
                                      fcore_members.empty() ? std::filesystem::path{}
                                                            : std::filesystem::path{fcore_members}));
    });

    // genopack pcore — unified dense prevalence-annotated per-genus aamer reference
    auto* pcore_cmd = app.add_subcommand("pcore",
        "Build SEC_PCORE: the unified dense per-genus aamer reference (every aamer + prevalence). "
        "Dense enough to detect SMALL foreign contigs (the conserved-only CORE is too sparse for 1-2kb).");
    std::string pcore_pack, pcore_members;
    int pcore_threads = 8;
    pcore_cmd->add_option("archive", pcore_pack, "Path to .gpk archive or directory of parts")->required();
    pcore_cmd->add_option("-t,--threads", pcore_threads, "Parallel shard readers (default: 8)");
    pcore_cmd->add_option("--members", pcore_members,
        "Reference accession list (one per line); only these genomes build the reference. Default: all live");
    pcore_cmd->callback([&]() {
        std::exit(genopack::cmd_pcore(std::filesystem::path(pcore_pack), pcore_threads,
                                      pcore_members.empty() ? std::filesystem::path{}
                                                            : std::filesystem::path{pcore_members}));
    });


    // genopack calibrate
    auto* cal_cmd = app.add_subcommand("calibrate",
        "Calibrate intrinsic completeness (aamer genus-core) against ground truth and write a "
        "CQAL section of calibrated per-genome columns. Ground truth: --checkm2 TSV, or the "
        "ingested XQAL if omitted. Also writes an isotonic+OLS JSON model and prints RMSE.");
    std::string cal_archive;
    std::string cal_checkm2;
    std::string cal_output  = "calibration.json";
    int         cal_threads = 8;
    cal_cmd->add_option("archive", cal_archive,
        "Path to .gpk archive or directory of parts")->required();
    cal_cmd->add_option("--checkm2", cal_checkm2,
        "CheckM2 TSV ground truth (optional; falls back to ingested XQAL CheckM2 columns)");
    cal_cmd->add_option("-o,--output", cal_output,
        "Output JSON path for calibration model (default: calibration.json)");
    cal_cmd->add_option("-t,--threads", cal_threads,
        "Threads (default: 8)");
    cal_cmd->callback([&]() {
        std::exit(genopack::calibrate::run_calibrate(
            cal_archive, cal_checkm2, cal_output, cal_threads));
    });

    // genopack cladesplit — protein-aamer clade-split contamination tier
    auto* cs_cmd = app.add_subcommand("cladesplit",
        "Protein-aamer clade-split contamination tier (GUNC-style, fast). Build a panel "
        "of genus-diagnostic aamers from clean genomes, then score genomes by lineage-split.");
    auto* cs_build = cs_cmd->add_subcommand("build", "Build a .csp panel from clean reference genomes");
    std::string cs_b_tsv, cs_b_out, cs_b_mode = "aa"; int cs_b_c = 30, cs_b_minaa = 8, cs_b_threads = 0;
    cs_build->add_option("-i,--input", cs_b_tsv, "TSV (accession, file_path, taxonomy)")->required();
    cs_build->add_option("-o,--output", cs_b_out, "Output .csp panel")->required();
    cs_build->add_option("-t,--threads", cs_b_threads, "Worker threads (default: all cores)");
    cs_build->add_option("--frac-c", cs_b_c, "FracMinHash density 1/c on aamers (default: 30)");
    cs_build->add_option("--min-aa", cs_b_minaa, "Min inter-stop AA segment length (default: 8)");
    cs_build->add_option("--primitive", cs_b_mode, "Protein k-mer primitive: aamer | metamer | strobemer (default: aamer)");
    cs_build->callback([&]() {
        genopack::CladeSplitConfig cfg; cfg.frac_c = cs_b_c; cfg.min_aa = cs_b_minaa;
        cfg.mode = (cs_b_mode == "metamer" || cs_b_mode == "aadna")  ? 1
                 : (cs_b_mode == "strobemer" || cs_b_mode == "strobe") ? 2 : 0;
        genopack::cladesplit_build(cs_b_tsv, cs_b_out, cfg, cs_b_threads);
        std::exit(0);
    });
    auto* cs_score = cs_cmd->add_subcommand("score", "Score genomes against a .csp panel (per-genome minority_fraction)");
    std::string cs_s_tsv, cs_s_gpk, cs_s_panel, cs_s_out, cs_s_mode = "aa"; int cs_s_c = 30, cs_s_minaa = 8, cs_s_threads = 0;
    cs_score->add_option("-i,--input", cs_s_tsv, "TSV (accession, file_path) — or use --gpk");
    cs_score->add_option("--gpk", cs_s_gpk, "Score live genomes inside a .gpk archive directly (flag contamination in a genopack file)");
    cs_score->add_option("--panel", cs_s_panel, "Panel .csp")->required();
    cs_score->add_option("-o,--output", cs_s_out, "Output per-genome TSV (incl. redundancy_fraction)")->required();
    cs_score->add_option("-t,--threads", cs_s_threads, "Worker threads (default: all cores)");
    cs_score->add_option("--frac-c", cs_s_c, "FracMinHash density 1/c (must match build; default: 30)");
    cs_score->add_option("--min-aa", cs_s_minaa, "Min AA segment (must match build; default: 8)");
    cs_score->add_option("--primitive", cs_s_mode, "Protein k-mer primitive (must match build): aamer | metamer | strobemer");
    cs_score->callback([&]() {
        genopack::CladeSplitConfig cfg; cfg.frac_c = cs_s_c; cfg.min_aa = cs_s_minaa;
        cfg.mode = (cs_s_mode == "metamer" || cs_s_mode == "aadna")  ? 1
                 : (cs_s_mode == "strobemer" || cs_s_mode == "strobe") ? 2 : 0;
        if (!cs_s_gpk.empty())      genopack::cladesplit_score_gpk(cs_s_gpk, cs_s_panel, cs_s_out, cfg, cs_s_threads);
        else if (!cs_s_tsv.empty()) genopack::cladesplit_score_tsv(cs_s_tsv, cs_s_panel, cs_s_out, cfg, cs_s_threads);
        else throw std::runtime_error("cladesplit score: provide -i <TSV> or --gpk <archive.gpk>");
        std::exit(0);
    });
    auto* cs_dump = cs_cmd->add_subcommand("aamers", "Dump per-genome sorted-unique aamer sets (binary) for core/completeness R&D");
    std::string cs_d_tsv, cs_d_out, cs_d_mode = "aa"; int cs_d_c = 30, cs_d_minaa = 8, cs_d_threads = 0;
    cs_dump->add_option("-i,--input", cs_d_tsv, "TSV (accession, file_path)")->required();
    cs_dump->add_option("-o,--output", cs_d_out, "Output binary aamer dump")->required();
    cs_dump->add_option("-t,--threads", cs_d_threads, "Worker threads (default: all cores)");
    cs_dump->add_option("--frac-c", cs_d_c, "FracMinHash density 1/c (default: 30)");
    cs_dump->add_option("--min-aa", cs_d_minaa, "Min AA segment (default: 8)");
    cs_dump->add_option("--primitive", cs_d_mode, "Protein k-mer primitive: aamer | metamer | strobemer (default: aamer)");
    cs_dump->callback([&]() {
        genopack::CladeSplitConfig cfg; cfg.frac_c = cs_d_c; cfg.min_aa = cs_d_minaa;
        cfg.mode = (cs_d_mode == "metamer" || cs_d_mode == "aadna")  ? 1
                 : (cs_d_mode == "strobemer" || cs_d_mode == "strobe") ? 2 : 0;
        genopack::cladesplit_dump_aamers(cs_d_tsv, cs_d_out, cfg, cs_d_threads);
        std::exit(0);
    });

    // genopack bench-grid
    auto* bg_cmd = app.add_subcommand("bench-grid",
        "Heterogeneous spike-fraction × ANI-distance benchmark from a manifest TSV");
    std::string bg_archive, bg_manifest, bg_output, bg_completeness_str;
    int bg_threads = 4, bg_reps = 5;
    uint32_t bg_seed = 42;
    bg_cmd->add_option("archive",          bg_archive,          "Archive path (.gpk)")->required();
    bg_cmd->add_option("-o,--output",      bg_output,           "TSV output path")->required();
    bg_cmd->add_option("--manifest",       bg_manifest,         "Manifest TSV (host_genus, contam_genus, ani_label, intra_offset)")->required();
    bg_cmd->add_option("-t,--threads",     bg_threads,          "Thread count");
    bg_cmd->add_option("-r,--reps",        bg_reps,             "Replicates per cell");
    bg_cmd->add_option("--seed",           bg_seed,             "Random seed");
    bg_cmd->add_option("--completeness",   bg_completeness_str, "Comma-separated host completeness fractions 0.0-1.0 (default: 1.0)");
    bg_cmd->callback([&]() {
        std::vector<float> comp_fracs = {1.0f};
        if (!bg_completeness_str.empty()) {
            comp_fracs.clear();
            std::istringstream ss(bg_completeness_str);
            std::string tok;
            while (std::getline(ss, tok, ','))
                comp_fracs.push_back(std::stof(tok));
        }
        std::exit(genopack::bench::cmd_bench_grid(
            bg_archive, bg_manifest, bg_output, bg_threads, bg_reps, bg_seed,
            comp_fracs));
    });

    // genopack sim
    auto* sim_cmd = app.add_subcommand("sim",
        "Generate synthetic fragmented/contaminated genomes for completeness/"
        "contamination benchmarking (CheckM2-compatible 20kb chunk fragmentation).");
    std::vector<std::string> sim_refs, sim_tax;
    // --contam, --contam-label, --contam-taxonomy are parallel repeatable arrays.
    // --contam-self adds ref genome itself as a contam source (tests redundancy signal).
    std::vector<std::string> sim_contams, sim_contam_labels, sim_contam_tax;
    bool        sim_contam_self = false;
    std::string sim_comp_str = "0.1,0.3,0.5,0.7,0.9,1.0";
    std::string sim_cont_str = "0.0", sim_outdir, sim_out_tsv, sim_manifest_tsv;
    int sim_reps = 3, sim_chunk = 20000, sim_min_chunk = 1000, sim_threads = 4;
    uint64_t sim_seed = 42;
    sim_cmd->add_option("--ref",            sim_refs,         "Reference genome FASTA (repeatable)")->required();
    sim_cmd->add_option("--taxonomy",       sim_tax,          "GTDB taxonomy; one per --ref in order");
    sim_cmd->add_option("--contam",         sim_contams,      "Contamination source FASTA (repeatable; each adds a grid dimension)");
    sim_cmd->add_option("--contam-label",   sim_contam_labels,"Label for each --contam source (e.g. genus,family,phylum); one per --contam");
    sim_cmd->add_option("--contam-taxonomy",sim_contam_tax,   "GTDB taxonomy for each --contam source; one per --contam");
    sim_cmd->add_flag ("--contam-self",     sim_contam_self,  "Add each ref as its own contam source (tests marker_redundancy; label=self)");
    sim_cmd->add_option("--completeness",   sim_comp_str,     "Comma-separated completeness fractions 0.0-1.0");
    sim_cmd->add_option("--contamination",  sim_cont_str,     "Comma-separated contamination fractions 0.0-0.5");
    sim_cmd->add_option("--reps",           sim_reps,         "Replicates per combination (default: 3)");
    sim_cmd->add_option("--seed",           sim_seed,         "Base random seed (default: 42)");
    sim_cmd->add_option("--chunk-size",     sim_chunk,        "Fixed fragment size in bp (default: 20000; ignored if --contig-n50 set)");
    sim_cmd->add_option("--min-chunk",      sim_min_chunk,    "Min chunk size to keep in bp (default: 1000)");
    int    sim_contig_n50   = 0;
    double sim_contig_sigma = 1.2;
    sim_cmd->add_option("--contig-n50",    sim_contig_n50,   "N50 for lognormal contig length distribution (0=fixed --chunk-size)");
    sim_cmd->add_option("--contig-sigma",  sim_contig_sigma, "Lognormal sigma for contig lengths (default: 1.2)");
    sim_cmd->add_option("--output-dir",     sim_outdir,       "Output directory for FASTA files")->required();
    sim_cmd->add_option("--output-tsv",     sim_out_tsv,      "Ground-truth TSV (default: <output-dir>/sim_manifest.tsv)");
    sim_cmd->add_option("--manifest-tsv",   sim_manifest_tsv, "genopack-add manifest TSV (default: <output-dir>/add_manifest.tsv)");
    sim_cmd->add_option("-t,--threads",     sim_threads,      "Parallel workers (default: 4)");
    sim_cmd->callback([&]() {
        auto parse_fracs = [](const std::string& s) {
            std::vector<double> out;
            std::istringstream ss(s);
            std::string tok;
            while (std::getline(ss, tok, ','))
                if (!tok.empty()) out.push_back(std::stod(tok));
            return out;
        };
        if (!sim_tax.empty() && sim_tax.size() != sim_refs.size()) {
            spdlog::error("sim: --taxonomy count ({}) != --ref count ({})",
                          sim_tax.size(), sim_refs.size());
            std::exit(2);
        }
        if (!sim_contam_labels.empty() && sim_contam_labels.size() != sim_contams.size()) {
            spdlog::error("sim: --contam-label count ({}) != --contam count ({})",
                          sim_contam_labels.size(), sim_contams.size());
            std::exit(2);
        }
        genopack::SimConfig cfg;
        cfg.refs.reserve(sim_refs.size());
        for (size_t i = 0; i < sim_refs.size(); ++i)
            cfg.refs.push_back({sim_refs[i], sim_tax.empty() ? std::string{} : sim_tax[i]});
        // Build contam list from explicit --contam entries
        for (size_t k = 0; k < sim_contams.size(); ++k) {
            genopack::SimContam sc;
            sc.fasta    = sim_contams[k];
            sc.label    = sim_contam_labels.empty() ? std::string{} : sim_contam_labels[k];
            sc.taxonomy = sim_contam_tax.size() > k ? sim_contam_tax[k] : std::string{};
            cfg.contams.push_back(std::move(sc));
        }
        // --contam-self: for each ref, also add that ref as a contam source (shared label "self")
        // This is handled inside run_sim per-job; signal it via a sentinel SimContam with empty fasta.
        if (sim_contam_self) {
            for (size_t i = 0; i < sim_refs.size(); ++i) {
                genopack::SimContam sc;
                sc.fasta    = sim_refs[i];
                sc.label    = "self";
                sc.taxonomy = sim_tax.empty() ? std::string{} : sim_tax[i];
                cfg.contams.push_back(std::move(sc));
            }
        }
        cfg.completeness  = parse_fracs(sim_comp_str);
        cfg.contamination = parse_fracs(sim_cont_str);
        cfg.reps          = sim_reps;
        cfg.seed          = sim_seed;
        cfg.chunk_size    = sim_chunk;
        cfg.min_chunk     = sim_min_chunk;
        cfg.contig_n50    = sim_contig_n50;
        cfg.contig_sigma  = sim_contig_sigma;
        cfg.output_dir    = sim_outdir;
        cfg.output_tsv    = sim_out_tsv;
        cfg.manifest_tsv  = sim_manifest_tsv;
        cfg.threads       = sim_threads;
        std::exit(genopack::run_sim(cfg));
    });

    // genopack markers
    auto* markers_cmd = app.add_subcommand("markers",
        "Build or manage the SCG marker aamer panel (.mrk sidecar)");
    markers_cmd->require_subcommand(1);

    auto* markers_build_cmd = markers_cmd->add_subcommand("build",
        "Build markers.mrk from GTDB-Tk MSA files (bac120 + ar53)");
    std::string mb_gtdbtk_db, mb_output;
    float mb_threshold = 0.30f;
    int mb_threads = 1;
    markers_build_cmd->add_option("--gtdbtk-db", mb_gtdbtk_db,
        "GTDB-Tk reference database root")->required();
    markers_build_cmd->add_option("-o,--output", mb_output,
        "Output .mrk file")->required();
    markers_build_cmd->add_option("--threshold", mb_threshold,
        "Default presence threshold (default: 0.30)");
    markers_build_cmd->add_option("-t,--threads", mb_threads,
        "Thread count (default: 1)");
    int mb_scale = 1;
    markers_build_cmd->add_option("--scale", mb_scale,
        "FracMinHash scale factor: keep 1/N of hashes (default: 1 = keep all; ignored with --dayhoff6)");
    bool mb_dayhoff6 = false;
    markers_build_cmd->add_flag("--dayhoff6", mb_dayhoff6,
        "Build Dayhoff-6 k=12 syncmer profile pool (recommended; compact, robust to divergence)");
    float mb_ic_threshold = 0.25f;
    markers_build_cmd->add_option("--ic-threshold", mb_ic_threshold,
        "Min per-column IC fraction to include k-mer positions (default: 0.25; --dayhoff6 only)");
    float mb_expected_min_frac = 0.50f;
    markers_build_cmd->add_option("--expected-min-frac", mb_expected_min_frac,
        "A marker is counted as expected for a genus only if it is detectable (≥1 IC-passing "
        "syncmer) in at least this fraction of the genus's reference genomes. Mirrors CheckM2's "
        "~97%% single-copy universality criterion. (default: 0.95)");
    markers_build_cmd->callback([&]() {
        genopack::MarkersBuildConfig cfg;
        cfg.gtdbtk_db          = mb_gtdbtk_db;
        cfg.output             = mb_output;
        cfg.default_threshold  = mb_threshold;
        cfg.expected_min_frac  = mb_expected_min_frac;
        cfg.threads            = mb_threads;
        cfg.frac_scale         = static_cast<uint16_t>(mb_scale);
        cfg.dayhoff6           = mb_dayhoff6;
        cfg.ic_threshold       = mb_ic_threshold;
        genopack::build_markers_panel(cfg);
    });

    // genopack markers score
    auto* markers_score_cmd = markers_cmd->add_subcommand("score",
        "Score a FASTA file for SCG marker completeness");
    std::string ms_fasta, ms_mrk, ms_genus;
    int ms_min_hits = 1;
    markers_score_cmd->add_option("--fasta",   ms_fasta,  "Input FASTA file")->required();
    markers_score_cmd->add_option("--markers", ms_mrk,    "Markers .mrk file")->required();
    markers_score_cmd->add_option("--genus",   ms_genus,  "Target genus key (e.g. g__Escherichia)")->required();
    markers_score_cmd->add_option("--min-hits", ms_min_hits, "Min hits per marker to call present (default: 1)");
    markers_score_cmd->callback([&]() {
        genopack::MarkerReader mr;
        mr.open(ms_mrk);
        if (!mr.is_open()) throw std::runtime_error("cannot open markers: " + ms_mrk);
        mr.build_merged_pool();

        const uint64_t gh = genopack::GcovWriter::hash_genus(ms_genus);
        auto calib = mr.lookup_lineage(gh);
        if (!calib.valid()) {
            std::cout << "NA\tNA\tNA\n";
            return;
        }

        // Per-contig voting: each contig votes independently for each marker.
        // Dispatches to Dayhoff-6 k=12 syncmers or full-AA k=8 FracMinHash.
        const bool is_d6      = mr.is_dayhoff6();
        const uint64_t frac_max = mr.frac_max_hash();
        const int min_seg     = is_d6 ? 50 : genopack::AAMER_K;
        const uint32_t min_hits = static_cast<uint32_t>(ms_min_hits);

        const bool is_arc = (calib.header->domain == genopack::MRKR_DOMAIN_ARC);
        auto mh    = is_arc ? mr.merged_hashes_arc() : mr.merged_hashes_bac();
        auto mid   = is_arc ? mr.merged_ids_arc()    : mr.merged_ids_bac();
        const genopack::BlockedBloom& bloom =
            is_arc ? mr.merged_bloom_arc() : mr.merged_bloom_bac();
        const uint8_t n_markers = calib.header->n_markers;

        uint32_t contig_votes[173] = {};
        {
            std::ifstream fin(ms_fasta);
            if (!fin) throw std::runtime_error("cannot open FASTA: " + ms_fasta);
            std::string line, seq;
            std::vector<uint64_t> orf_mers;
            auto score_contig = [&]() {
                if (seq.empty()) return;
                uint32_t best[173] = {};
                if (is_d6) {
                    genopack::translate_6frame(seq, min_seg,
                        [&](int, const uint8_t* seg, int len, int, int) {
                            orf_mers.clear();
                            genopack::extract_d6_orf_syncmers(seg, len, orf_mers);
                            if (orf_mers.empty()) return;
                            std::sort(orf_mers.begin(), orf_mers.end());
                            orf_mers.erase(std::unique(orf_mers.begin(), orf_mers.end()), orf_mers.end());
                            uint32_t local[173] = {};
                            const uint64_t* mhp = mh.data(), *mhe = mhp + mh.size();
                            const uint8_t* mip = mid.data();
                            for (uint64_t h : orf_mers) {
                                if (!bloom.might_contain(h)) continue;
                                auto it = std::lower_bound(mhp, mhe, h);
                                if (it != mhe && *it == h)
                                    local[mip[it - mhp]]++;
                            }
                            const uint32_t q_sz = static_cast<uint32_t>(orf_mers.size());
                            const uint32_t thr = std::max(min_hits, std::max(1u, q_sz / 20u));
                            for (uint8_t mi = 0; mi < n_markers; ++mi)
                                if (local[mi] >= thr && local[mi] > best[mi])
                                    best[mi] = local[mi];
                        });
                } else {
                    orf_mers.clear();
                    genopack::extract_aamers_dna_into(seq, genopack::AAMER_K, min_seg,
                                                        frac_max, orf_mers);
                    if (!orf_mers.empty()) {
                        std::sort(orf_mers.begin(), orf_mers.end());
                        orf_mers.erase(std::unique(orf_mers.begin(), orf_mers.end()), orf_mers.end());
                        uint32_t hits[173] = {};
                                                const uint64_t* mhp = mh.data(), *mhe = mhp + mh.size();
                        const uint8_t* mip = mid.data();
                        for (uint64_t h : orf_mers) {
                            if (!bloom.might_contain(h)) continue;
                            auto it = std::lower_bound(mhp, mhe, h);
                            if (it != mhe && *it == h) hits[mip[it - mhp]]++;
                        }
                        const uint32_t q_sz = static_cast<uint32_t>(orf_mers.size());
                        const uint32_t thr = std::max(min_hits, std::max(1u, q_sz / 20u));
                        for (uint8_t mi = 0; mi < n_markers; ++mi)
                            if (hits[mi] >= thr && hits[mi] > best[mi])
                                best[mi] = hits[mi];
                    }
                }
                for (uint8_t mi = 0; mi < n_markers; ++mi)
                    if (best[mi] > 0) contig_votes[mi]++;
                seq.clear();
            };
            while (std::getline(fin, line)) {
                if (!line.empty() && line[0] == '>') score_contig();
                else seq += line;
            }
            score_contig();
        }

        int n_present = 0, n_expected = 0, n_redundant = 0;
        for (uint8_t mi = 0; mi < n_markers; ++mi) {
            if (!calib.marker_expected(mi)) continue;
            ++n_expected;
            if (contig_votes[mi] >= 1) ++n_present;
            if (contig_votes[mi] >= 2) ++n_redundant;
        }
        const float completeness = n_expected > 0
            ? static_cast<float>(n_present) / static_cast<float>(n_expected)
            : std::numeric_limits<float>::quiet_NaN();
        float redundancy = n_expected > 0
            ? static_cast<float>(n_redundant) / static_cast<float>(n_expected)
            : std::numeric_limits<float>::quiet_NaN();

        auto fmtf = [](float v) -> std::string {
            return std::isnan(v) ? "NA" : std::to_string(v);
        };
        std::cout << fmtf(completeness) << "\t" << fmtf(redundancy)
                  << "\t" << n_present << "\t" << n_expected << "\n";
    });

    CLI11_PARSE(app, argc, argv);
    return 0;
}
