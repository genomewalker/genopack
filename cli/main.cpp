#include <genopack/accx.hpp>
#include <genopack/archive.hpp>
#include <genopack/repack.hpp>
#include <genopack/catalog.hpp>
#include <genopack/cidx.hpp>
#include <genopack/format.hpp>
#include <genopack/gidx.hpp>
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
#include <thread>
#include <unordered_map>
#include <unordered_set>

using namespace genopack;

// ── genopack build ─────────────────────────────────────────────────────────────
static int cmd_build(const std::string& input_tsv, const std::string& output_dir,
                     int threads, int zstd_level, bool no_dict, bool ref_dict,
                     bool delta, bool mem_delta, bool verbose, int n_parallel,
                     bool no_cidx, bool use_2bit, bool kmer_nn_sort,
                     bool taxonomy_group, const std::string& taxonomy_rank,
                     bool sketch, int sketch_kmer, int sketch_size, int sketch_syncmer,
                     std::vector<int> sketch_kmers = {}) {
    ArchiveBuilder::Config cfg;
    const bool explicit_codec = no_dict || ref_dict || delta || mem_delta;
    cfg.io_threads                        = static_cast<size_t>(std::max(1, threads));
    cfg.shard_cfg.compress_threads        = static_cast<size_t>(std::max(1, threads));
    cfg.verbose                           = verbose;
    cfg.build_cidx                        = !no_cidx;
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

    std::vector<std::future<void>> futs;
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
            [part_cfg, part_path, slice = std::move(slice), p]() mutable {
                // Skip building if a valid complete archive already exists.
                // Validate by checking the TailLocator magic at end of file.
                if (std::filesystem::exists(part_path)) {
                    bool valid = false;
                    try {
                        MmapFileReader mm;
                        mm.open(part_path);
                        if (mm.size() >= sizeof(TailLocator)) {
                            const auto* tail = mm.ptr_at<TailLocator>(
                                mm.size() - sizeof(TailLocator));
                            valid = (tail->magic == GPKT_MAGIC);
                        }
                    } catch (...) {}
                    if (valid) {
                        spdlog::info("Part {}: skipping — valid archive exists ({})",
                                     p, part_path.string());
                        return;
                    }
                    spdlog::warn("Part {}: existing file invalid, rebuilding", p);
                    std::filesystem::remove(part_path);
                }
                spdlog::info("Part {}: {} genomes (gid_start={}) → {}", p, slice.size(),
                             part_cfg.starting_genome_id, part_path.string());
                ArchiveBuilder builder(part_path, part_cfg);
                for (auto& r : slice) builder.add(r);
                builder.finalize();
                spdlog::info("Part {} done", p);
            }));
    }

    for (auto& f : futs) f.get();
    spdlog::info("All {} parts built, merging…", part_paths.size());

    merge_archives(part_paths, output_dir, /*remap_genome_ids=*/false, cfg.build_cidx);

    // Cleanup temp parts
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

    // --output-dir: write one {accession}.fa per genome, open archive once
    if (!out_dir.empty()) {
        namespace fs = std::filesystem;
        fs::create_directories(out_dir);
        if (q.accessions.empty()) {
            spdlog::error("--output-dir requires --accession or --accessions-file");
            return 1;
        }
        auto results = ar.batch_fetch_by_accessions(q.accessions);
        size_t written = 0;
        for (size_t i = 0; i < q.accessions.size(); ++i) {
            if (!results[i]) { spdlog::warn("Accession not found: {}", q.accessions[i]); continue; }
            auto path = fs::path(out_dir) / (q.accessions[i] + ".fa");
            std::ofstream of(path);
            if (!of) { spdlog::error("Cannot write: {}", path.string()); return 1; }
            of << results[i]->fasta;
            ++written;
        }
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
    ArchiveReader ar;
    ar.open(archive_dir);

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

// ── genopack taxonomy ──────────────────────────────────────────────────────────
static int cmd_taxonomy(const std::string& archive_dir, const std::string& accession,
                        bool json) {
    ArchiveReader ar;
    ar.open(archive_dir);

    auto tree_opt = ar.taxonomy_tree();
    if (!tree_opt) {
        spdlog::error("No taxonomy data in archive");
        return 1;
    }
    auto& tree = *tree_opt;

    if (!accession.empty()) {
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
    size_t nn = tree.n_nodes();
    size_t na = tree.n_accessions();
    if (json) {
        std::cout << "{\n"
                  << "  \"n_nodes\": "       << nn << ",\n"
                  << "  \"n_accessions\": "  << na << "\n"
                  << "}\n";
    } else {
        std::cout << "Taxonomy summary for: " << archive_dir << "\n"
                  << "  Nodes:       " << nn << "\n"
                  << "  Accessions:  " << na << "\n";
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
                       bool repack_skch = false) {
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
    if (fh->magic != GPK2_MAGIC) {
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
    bool need_skch = build_skch || repack_skch;

    if (!need_gidx && !need_txdb && !need_cidx && !need_skch) {
        spdlog::info("GIDX: already present (use --force to rebuild)");
        if (build_txdb) spdlog::info("TXDB: already present (use --force to rebuild)");
        if (!cidx_tsv.empty()) spdlog::info("CIDX: already present (use --force to rebuild)");
        spdlog::info("Nothing to do");
        return 0;
    }

    // ── Build GIDX ──────────────────────────────────────────────────────────

    // Step 1: scan all CATL sections to collect genome_id -> (shard_id, catl_row_index)
    struct CatlEntry {
        uint32_t shard_id;
        uint64_t catl_row_index;
    };
    std::unordered_map<GenomeId, CatlEntry> genome_map;

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
    auto shrd_sections = toc.find_by_type(SEC_SHRD);
    GidxWriter gidx_writer;
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

    // Writers declared here so the repack block can assign into skch_writer_mk
    // before the normal sketching block runs (which is skipped when repack_skch).
    std::unique_ptr<SkchWriter>       skch_writer;
    std::unique_ptr<SkchWriterMultiK> skch_writer_mk;

    // ── Repack existing SKCH into seekable format (if --repack-skch) ────────────
    // Reads all sketches from the existing section and re-emits in seekable
    // multi-frame format without re-computing OPH signatures from FASTA.
    if (repack_skch) {
        auto skch_secs = toc.find_by_type(SEC_SKCH);
        if (skch_secs.empty()) {
            spdlog::error("--repack-skch: no SKCH section found in archive");
            return 1;
        }
        SkchReader src;
        src.open(mmap.data(), skch_secs[0]->file_offset, skch_secs[0]->compressed_size);

        if (src.version() == 3) {
            spdlog::info("SKCH: already seekable — nothing to repack");
            repack_skch = false;
            need_skch   = false;
        } else {
            src.ensure_loaded();

            const uint32_t nk = src.n_kmer_sizes();
            std::vector<uint32_t> ks;
            for (uint32_t i = 0; i < nk; ++i) ks.push_back(src.kmer_size_at(i));
            if (ks.empty()) ks.push_back(src.kmer_size());
            const uint32_t sketch_sz = src.sketch_size();
            const uint32_t mw        = (sketch_sz + 63) / 64;

            skch_writer_mk = std::make_unique<SkchWriterMultiK>(
                ks, sketch_sz, src.syncmer_s(), src.seed1(), src.seed2());

            size_t done = 0;
            for (GenomeId gid : src.genome_ids()) {
                std::vector<std::vector<uint16_t>> sigs(nk);
                std::vector<uint32_t>              nreal(nk);
                std::vector<std::vector<uint64_t>> masks(nk);
                uint64_t genome_length = 0;
                for (uint32_t ki = 0; ki < nk; ++ki) {
                    auto sk = src.sketch_for(gid, ks[ki], sketch_sz);
                    if (!sk) continue;
                    genome_length   = sk->genome_length;
                    sigs[ki].assign(sk->sig,  sk->sig  + sketch_sz);
                    nreal[ki]       = sk->n_real_bins;
                    masks[ki].assign(sk->mask, sk->mask + mw);
                }
                skch_writer_mk->add(gid, genome_length, sigs, nreal, masks);
                ++done;
                if (done % 100000 == 0)
                    spdlog::info("SKCH repack: {}/{} converted", done, src.n_genomes());
            }
            src.release();
            spdlog::info("SKCH repack: {} genomes read from existing section", done);
            need_skch = false;  // skip the normal FASTA-based sketching block
        }
    }

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
        uint32_t sk_syncmer_s   = (skch_syncmer > 0) ? static_cast<uint32_t>(skch_syncmer) : 0;
        uint64_t sk_seed1       = 42;
        uint64_t sk_seed2       = 1337;

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
                        std::vector<std::vector<uint16_t>> sigs;   // [ki]
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
                            mr.sigs.resize(nk);
                            mr.n_real.resize(nk);
                            mr.masks.resize(nk);
                            for (size_t ki = 0; ki < nk; ++ki) {
                                OPHSketchConfig sc;
                                sc.kmer_size   = skch_kmers[ki];
                                sc.sketch_size = static_cast<int>(sk_sketch_size);
                                sc.syncmer_s   = static_cast<int>(sk_syncmer_s);
                                sc.seed        = sk_seed1;
                                auto sk = sketch_oph_from_buffer(fasta.data(), fasta.size(), sc);
                                mr.genome_length = sk.genome_length;
                                mr.n_real[ki]    = sk.n_real_bins;
                                mr.masks[ki]     = sk.real_bins_bitmask;
                                mr.sigs[ki].resize(sk.signature.size());
                                for (size_t si = 0; si < sk.signature.size(); ++si)
                                    mr.sigs[ki][si] = static_cast<uint16_t>(sk.signature[si] >> 16);
                            }
                        } catch (...) {
                            local_failed++;
                        }
                    }
                    for (auto& mr : results) {
                        if (mr.sigs.empty() || mr.sigs[0].empty()) continue;
                        skch_writer_mk->add(mr.gid, mr.genome_length,
                                            mr.sigs, mr.n_real, mr.masks);
                        ++done;
                    }
                } else {
                    // Single-k path (unchanged).
                    OPHSketchConfig sc;
                    sc.kmer_size   = static_cast<int>(sk_kmer_size);
                    sc.sketch_size = static_cast<int>(sk_sketch_size);
                    sc.syncmer_s   = static_cast<int>(sk_syncmer_s);
                    sc.seed        = sk_seed1;
                    std::vector<std::pair<GenomeId, OPHSketchResult>> results(static_cast<size_t>(n_gids));

                    #pragma omp parallel for schedule(dynamic, 1) num_threads(skch_threads)
                    for (int j = 0; j < n_gids; ++j) {
                        GenomeId gid = gids[static_cast<size_t>(j)];
                        try {
                            std::string fasta = shard.fetch_genome(gid);
                            if (fasta.empty()) { local_failed++; continue; }
                            auto sk = sketch_oph_from_buffer(fasta.data(), fasta.size(), sc);
                            results[static_cast<size_t>(j)] = {gid, std::move(sk)};
                        } catch (...) {
                            local_failed++;
                        }
                    }
                    for (auto& [gid, sk] : results) {
                        if (sk.signature.empty()) continue;
                        std::vector<uint16_t> sig16(sk.signature.size());
                        for (size_t si = 0; si < sk.signature.size(); ++si)
                            sig16[si] = static_cast<uint16_t>(sk.signature[si] >> 16);
                        skch_writer->add(gid, sig16, sk.n_real_bins, sk.genome_length,
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

    // Copy existing sections into new TocWriter (excluding rebuilt types if --force)
    TocWriter new_toc;
    for (const auto& s : toc.sections) {
        if (force && s.type == SEC_GIDX) continue;
        if (force && need_txdb && s.type == SEC_TXDB) continue;
        if (force && need_cidx && s.type == SEC_CIDX) continue;
        if (repack_skch && s.type == SEC_SKCH) continue;  // drop old section; new one appended below
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

    // Write new TOC + TailLocator
    new_toc.finalize(writer,
                     prev_generation + 1,
                     live_count, total_count,
                     prev_toc_offset,
                     catalog_root_id, accession_root_id, tombstone_root_id);

    writer.flush();

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
    if (fh->magic != GPK2_MAGIC) {
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

// ── main ──────────────────────────────────────────────────────────────────────
int main(int argc, char** argv) {
    CLI::App app{"genopack — genome archive"};
    app.require_subcommand(1);

    // genopack build
    auto* build = app.add_subcommand("build", "Build a new archive from a genome TSV (default codec mode auto-selects plain vs delta)");
    std::string build_input, build_output;
    int build_threads = 16, build_level = 6, build_parallel = 1;
    bool build_no_dict = false, build_ref_dict = false, build_delta = false;
    bool build_mem_delta = false, build_verbose = false, build_no_cidx = false;
    bool build_2bit = false, build_kmer_nn = true, build_taxon_group = true;
    bool build_sketch = true;
    int build_sketch_kmer = 16, build_sketch_size = 10000, build_sketch_syncmer = 0;
    std::string build_taxon_rank = "g";
    std::string build_sketch_kmers_str;
    build->add_option("-i,--input",  build_input,  "Input TSV (accession, file_path, ...)")->required();
    build->add_option("-o,--output", build_output, "Output archive directory (.gpk)")->required();
    build->add_option("-t,--threads", build_threads, "I/O threads (decompression + compression)");
    build->add_option("-z,--zstd-level", build_level, "zstd compression level (1-22)");
    build->add_option("-p,--parallel", build_parallel, "Number of parallel build workers (auto-merge)");
    build->add_flag("--no-dict",    build_no_dict,    "Disable shared dictionary training");
    build->add_flag("--ref-dict",   build_ref_dict,   "Use first genome in each shard as reference content dictionary");
    build->add_flag("--delta",      build_delta,      "Compress non-reference blobs against first genome using zstd prefix");
    build->add_flag("--mem-delta",  build_mem_delta,  "Force MEM-delta: k=31 k-mer seeded exact-match encoding for highly similar shard groups");
    build->add_flag("--no-cidx",    build_no_cidx,    "Skip CIDX contig index (recommended for >1M genomes)");
    build->add_flag("--2bit",       build_2bit,       "Pack nucleotide sequence to 2 bits/base before zstd (~1.5-2x additional compression)");
    build->add_flag("--kmer-sort,!--no-kmer-sort",   build_kmer_nn,    "Sort genomes within each shard by kmer4_profile NN chain (default: on; --no-kmer-sort to disable)");
    build->add_flag("--taxon-group,!--no-taxon-group",build_taxon_group,"Group genomes into per-taxon shards (default: on; --no-taxon-group to disable; requires taxonomy column)");
    build->add_option("--taxon-rank",build_taxon_rank,"Taxonomy rank for grouping: g=genus (default), f=family");
    build->add_flag("--sketch,!--no-sketch", build_sketch, "Compute OPH sketches (default: on; use --no-sketch to disable)");
    build->add_option("--sketch-kmer", build_sketch_kmer, "OPH sketch k-mer size (default: 16)");
    build->add_option("--sketch-kmers", build_sketch_kmers_str, "Comma-separated k-mer sizes for multi-k SKCH v2 (e.g. 16,21,31)");
    build->add_option("--sketch-size", build_sketch_size, "Number of OPH bins (default: 10000)");
    build->add_option("--sketch-syncmer", build_sketch_syncmer, "Open syncmer prefilter s (0=disabled)");
    build->add_flag("-v,--verbose", build_verbose, "Verbose progress");
    std::string build_coordinator; // "manifest_dir:output.gpk" or empty
    build->add_option("--coordinator", build_coordinator,
        "NFS manifest coordinator: manifest_dir:/path/to/output.gpk. "
        "Build to a temp file then transfer sections via NFS manifest protocol. "
        "The legacy 'nfs:' prefix is accepted and stripped for backward compatibility.");
    build->callback([&]() {
        // Parse --sketch-kmers "16,21,31" into vector<int>
        std::vector<int> build_sketch_kmers;
        if (!build_sketch_kmers_str.empty()) {
            std::istringstream ss(build_sketch_kmers_str);
            std::string tok;
            while (std::getline(ss, tok, ','))
                if (!tok.empty()) build_sketch_kmers.push_back(std::stoi(tok));
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
                                build_sketch_syncmer, build_sketch_kmers);
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
                             build_sketch_syncmer, build_sketch_kmers));
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
    reindex_cmd->add_flag("--no-gidx", reindex_no_gidx, "Skip GIDX build (useful when only --skch is needed and GIDX is absent/unwanted)");
    reindex_cmd->add_flag("--txdb", reindex_txdb, "Build taxonomy tree (TXDB) from TAXN lineage strings");
    reindex_cmd->add_option("--cidx", reindex_cidx_tsv, "Build contig accession index (CIDX) from build TSV (accession<TAB>taxonomy<TAB>file_path)");
    reindex_cmd->add_option("--cidx-threads", reindex_cidx_threads, "Threads for parallel FASTA decompression (default: 8)");
    reindex_cmd->add_flag("--skch", reindex_skch, "Compute OPH sketches for genomes missing from existing SKCH sections");
    bool reindex_repack_skch = false;
    reindex_cmd->add_flag("--repack-skch", reindex_repack_skch,
        "Convert existing SKCH section to seekable multi-frame format without re-computing sketches. "
        "Reads all signatures from the current section and re-emits them as independent compressed frames "
        "so that per-taxon sketch access decompresses only the relevant rows.");
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
                              reindex_no_gidx,
                              reindex_repack_skch));
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

    CLI11_PARSE(app, argc, argv);
    return 0;
}
