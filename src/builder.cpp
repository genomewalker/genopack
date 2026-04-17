#include <genopack/archive.hpp>
#include <genopack/accx.hpp>
#include <genopack/oph_sketch.hpp>
#include <genopack/skch.hpp>
#include <genopack/cidx.hpp>
#include <genopack/gidx.hpp>
#include <genopack/kmrx.hpp>
#include <genopack/taxn.hpp>
#include <genopack/txdb.hpp>
#include <genopack/catalog.hpp>
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/shard.hpp>
#include <genopack/toc.hpp>
#include <genopack/util.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cassert>
#include <condition_variable>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <future>
#include <mutex>
#include <chrono>
#include <numeric>
#include <optional>
#include <queue>
#include <stdexcept>
#include <string_view>
#include <thread>
#include <unistd.h>
#include <fcntl.h>
#include <unordered_map>
#include <vector>

namespace genopack {

// ── Shard grouping helpers ────────────────────────────────────────────────────

// Extract taxonomy bucket key at the requested rank ('g'=genus, 'f'=family).
// Falls back to the next coarser rank if the requested rank is absent.
// Returns "__unclassified__" if taxonomy is empty.
static std::string extract_taxonomy_bucket(std::string_view taxonomy,
                                           char rank)
{
    if (taxonomy.empty()) return "__unclassified__";
    // Taxonomy format: "d__;p__;c__;o__;f__;g__;s__"
    // Each rank prefix: "g__", "f__", etc.
    char prefix[4] = {rank, '_', '_', ';'};
    auto pos = taxonomy.find(std::string_view(prefix, 3));
    if (pos == std::string_view::npos) {
        // Fallback: try family then order
        const char* fallbacks = (rank == 'g') ? "foc" : "oc";
        for (char fb : std::string_view(fallbacks)) {
            char fp[4] = {fb, '_', '_', ';'};
            pos = taxonomy.find(std::string_view(fp, 3));
            if (pos != std::string_view::npos) break;
        }
        if (pos == std::string_view::npos) return "__unclassified__";
    }
    size_t start = pos;
    size_t end   = taxonomy.find(';', start + 3);
    if (end == std::string_view::npos) end = taxonomy.size();
    return std::string(taxonomy.substr(start, end - start));
}

// Greedy nearest-neighbor chain over 136-dim kmer4_profiles.
// Returns permutation indices that order items for maximum within-shard similarity.
// Uses the centroid as the starting point.
template<typename Item>
static std::vector<size_t> greedy_nn_chain(const std::vector<Item>& items)
{
    const size_t n = items.size();
    if (n <= 1) {
        std::vector<size_t> idx(n);
        std::iota(idx.begin(), idx.end(), 0);
        return idx;
    }

    // Compute centroid
    std::array<float, 136> centroid{};
    for (const auto& item : items)
        for (int d = 0; d < 136; ++d)
            centroid[d] += item.stats.kmer4_profile[d];
    for (int d = 0; d < 136; ++d)
        centroid[d] /= static_cast<float>(n);

    auto dot = [](const float* a, const float* b) {
        float s = 0.f;
        for (int d = 0; d < 136; ++d) s += a[d] * b[d];
        return s;
    };

    // Find genome nearest to centroid as starting point
    size_t start = 0;
    float  best  = -2.f;
    for (size_t i = 0; i < n; ++i) {
        float sim = dot(items[i].stats.kmer4_profile.data(), centroid.data());
        if (sim > best) { best = sim; start = i; }
    }

    std::vector<size_t> order;
    order.reserve(n);
    std::vector<bool> visited(n, false);
    size_t cur = start;

    while (order.size() < n) {
        visited[cur] = true;
        order.push_back(cur);
        if (order.size() == n) break;
        // Find nearest unvisited
        size_t next = n;
        float  best_sim = -2.f;
        const float* cur_prof = items[cur].stats.kmer4_profile.data();
        for (size_t j = 0; j < n; ++j) {
            if (visited[j]) continue;
            float sim = dot(cur_prof, items[j].stats.kmer4_profile.data());
            if (sim > best_sim) { best_sim = sim; next = j; }
        }
        cur = next;
    }
    return order;
}

// ── Checkpoint binary format ──────────────────────────────────────────────────
//
// {gpk}.ckpt         — text key=value (updated atomically after each shard):
//   genome_count=N   total genomes written so far
//   byte_offset=B    file offset after last written shard
//   next_genome_id=G next genome_id to assign
//   current_shard_id=S next shard_id to open
//   next_section_id=I next section_id to assign
//
// {gpk}.ckpt_meta.bin — binary stream of shard blocks:
//   [ShardCkptHdr] [GenomeCkptFixed + accession bytes + taxonomy bytes] × n_genomes
//
// On resume: restore all in-memory state, truncate .gpk to byte_offset,
// open in append mode, skip first genome_count TSV rows.

#pragma pack(push, 1)
struct ShardCkptHdr {
    uint32_t shard_id;
    uint32_t n_genomes;
    uint64_t section_id;
    uint64_t file_offset;
    uint64_t compressed_size;
};  // 28 bytes

struct GenomeCkptFixed {
    uint64_t genome_id;
    uint64_t input_row_index;
    uint32_t shard_id;
    uint32_t dir_index;
    uint64_t oph_fingerprint;
    uint64_t genome_length;
    uint32_t n_contigs;
    uint32_t gc_pct_x100;
    uint16_t completeness_x10;
    uint16_t contamination_x10;
    uint32_t date_added;
    uint16_t accession_len;
    uint16_t taxonomy_len;
    uint32_t genome_type;   // GenomeType enum
    uint32_t _ckpt_pad0;
    float    kmer4[136];
};  // 612 bytes
#pragma pack(pop)

// ── ArchiveBuilder::Impl ──────────────────────────────────────────────────────

struct ArchiveBuilder::Impl {
    std::filesystem::path archive_dir;   // base path (no extension)
    std::filesystem::path gpk_path_;     // output .gpk file
    Config                cfg;

    std::vector<BuildRecord> pending;
    std::vector<uint64_t>    pending_input_rows;
    GenomeId next_genome_id = 1;  // overridden by cfg.starting_genome_id below

    // Pass-1 result: metadata only, FASTA discarded after stats computation
    struct GenomeMeta1 {
        BuildRecord record;
        GenomeId    genome_id;
        FastaStats  stats;
    };

    explicit Impl(const std::filesystem::path& dir, Config c)
        : archive_dir(dir), cfg(c), next_genome_id(c.starting_genome_id)
    {
        // Determine output .gpk path
        if (dir.extension() == ".gpk")
            gpk_path_ = dir;
        else
            gpk_path_ = std::filesystem::path(dir.string() + ".gpk");

        // Ensure parent directory exists (for meta.tsv sidecar)
        std::filesystem::create_directories(gpk_path_.parent_path());
    }

    void add(const BuildRecord& rec) {
        pending.push_back(rec);
        pending_input_rows.push_back(static_cast<uint64_t>(pending_input_rows.size()));
    }

    void add_from_tsv(const std::filesystem::path& tsv_path) {
        auto records = parse_tsv_records(tsv_path);
        uint64_t base = pending_input_rows.size();
        for (auto& r : records)
            pending.push_back(std::move(r));
        for (size_t i = 0; i < records.size(); ++i)
            pending_input_rows.push_back(base + i);
        spdlog::info("Loaded {} records from {}", pending.size(), tsv_path.string());
    }

    void finalize() {
        if (pending.empty()) { spdlog::warn("No records to build"); return; }

        const size_t original_total_records = pending.size();
        size_t total_records = pending.size();
        spdlog::info("Building archive: {} genomes, {} threads", total_records, cfg.io_threads);

        // ── Checkpoint paths ──────────────────────────────────────────────────
        std::filesystem::path ckpt_path      = std::filesystem::path(gpk_path_.string() + ".ckpt");
        std::filesystem::path ckpt_meta_path = std::filesystem::path(gpk_path_.string() + ".ckpt_meta.bin");

        AppendWriter app_writer;
        TocWriter    toc;
        uint64_t next_section_id = 1;
        ShardId  current_shard_id = 0;

        std::vector<GenomeMeta> catalog_rows;
        catalog_rows.reserve(total_records);
        std::vector<std::pair<std::string, GenomeId>> accession_pairs;
        accession_pairs.reserve(total_records);
        std::unordered_map<std::string, std::string> taxonomy_map;
        taxonomy_map.reserve(total_records);
        std::vector<std::pair<GenomeId, std::array<float, 136>>> kmer_pairs;
        kmer_pairs.reserve(total_records);

        // OPH sketch writer: spills sigs+masks to tmpfile immediately, no RAM accumulation
        std::unique_ptr<SkchWriter>       skch_writer;
        std::unique_ptr<SkchWriterMultiK> skch_writer_mk;
        const bool multi_k_sketch = cfg.build_sketch && cfg.sketch_kmer_sizes.size() > 1;
        if (cfg.build_sketch) {
            if (multi_k_sketch) {
                std::vector<uint32_t> ks;
                for (int k : cfg.sketch_kmer_sizes) ks.push_back(static_cast<uint32_t>(k));
                skch_writer_mk = std::make_unique<SkchWriterMultiK>(
                    ks, static_cast<uint32_t>(cfg.sketch_size),
                    static_cast<uint32_t>(cfg.sketch_syncmer_s),
                    cfg.sketch_seed, cfg.sketch_seed + 1,
                    gpk_path_.parent_path().string());
            } else {
                skch_writer = std::make_unique<SkchWriter>(
                    static_cast<uint32_t>(cfg.sketch_size),
                    static_cast<uint32_t>(cfg.sketch_kmer_size),
                    static_cast<uint32_t>(cfg.sketch_syncmer_s),
                    cfg.sketch_seed, cfg.sketch_seed + 1,
                    gpk_path_.parent_path().string());
            }
        }

        std::vector<uint64_t> source_row_indices;
        source_row_indices.reserve(original_total_records);

        struct GidxInfo { GenomeId genome_id; ShardId shard_id; uint32_t dir_index; };
        std::vector<GidxInfo> gidx_infos;
        gidx_infos.reserve(total_records);
        std::unordered_map<ShardId, uint64_t> shard_id_to_section_id;
        CidxWriter cidx_writer;

        // ── Try resume ───────────────────────────────────────────────────────
        bool     resuming      = false;

        if (std::filesystem::exists(ckpt_path) && std::filesystem::exists(ckpt_meta_path)
            && std::filesystem::exists(gpk_path_)) {

            size_t   ck_genome_count  = 0;
            uint64_t ck_byte_offset   = 0;
            uint64_t ck_next_gid      = 0;
            uint32_t ck_shard_id      = 0;
            uint64_t ck_section_id    = 0;

            // Parse .ckpt
            std::ifstream ck(ckpt_path);
            std::string line;
            bool parsed_ok = false;
            int parsed_fields = 0;
            while (std::getline(ck, line)) {
                auto eq = line.find('=');
                if (eq == std::string::npos) continue;
                std::string key = line.substr(0, eq);
                std::string val = line.substr(eq + 1);
                if      (key == "genome_count")    { ck_genome_count  = std::stoull(val); ++parsed_fields; }
                else if (key == "byte_offset")     { ck_byte_offset   = std::stoull(val); ++parsed_fields; }
                else if (key == "next_genome_id")  { ck_next_gid      = std::stoull(val); ++parsed_fields; }
                else if (key == "current_shard_id"){ ck_shard_id      = static_cast<uint32_t>(std::stoull(val)); ++parsed_fields; }
                else if (key == "next_section_id") { ck_section_id    = std::stoull(val); ++parsed_fields; }
            }
            parsed_ok = (parsed_fields == 5);

            // Read .ckpt_meta.bin and restore all in-memory state
            if (parsed_ok && ck_genome_count > 0 && ck_genome_count < total_records) {
                FILE* mf = std::fopen(ckpt_meta_path.c_str(), "rb");
                bool meta_ok = (mf != nullptr);
                size_t genomes_restored = 0;
                std::vector<char> processed_rows(total_records, 0);

                while (meta_ok && genomes_restored < ck_genome_count) {
                    ShardCkptHdr shdr{};
                    if (std::fread(&shdr, sizeof(shdr), 1, mf) != 1) { meta_ok = false; break; }

                    // Reconstruct TOC entry for this shard
                    SectionDesc sd{};
                    sd.type             = SEC_SHRD;
                    sd.version          = 4;
                    sd.flags            = 0;
                    sd.section_id       = shdr.section_id;
                    sd.file_offset      = shdr.file_offset;
                    sd.compressed_size  = shdr.compressed_size;
                    sd.uncompressed_size = 0;
                    sd.item_count       = shdr.n_genomes;
                    sd.aux0             = shdr.shard_id;
                    sd.aux1             = 0;
                    std::memset(sd.checksum, 0, sizeof(sd.checksum));
                    toc.add_section(sd);
                    shard_id_to_section_id[shdr.shard_id] = shdr.section_id;

                    for (uint32_t gi = 0; gi < shdr.n_genomes && meta_ok; ++gi) {
                        GenomeCkptFixed gf{};
                        if (std::fread(&gf, sizeof(gf), 1, mf) != 1) { meta_ok = false; break; }
                        if (gf.input_row_index >= processed_rows.size()) { meta_ok = false; break; }
                        processed_rows[gf.input_row_index] = 1;

                        std::string accession(gf.accession_len, '\0');
                        std::string taxonomy(gf.taxonomy_len, '\0');
                        if (gf.accession_len &&
                            std::fread(accession.data(), 1, gf.accession_len, mf) != gf.accession_len) {
                            meta_ok = false; break;
                        }
                        if (gf.taxonomy_len &&
                            std::fread(taxonomy.data(), 1, gf.taxonomy_len, mf) != gf.taxonomy_len) {
                            meta_ok = false; break;
                        }

                        GenomeMeta meta{};
                        meta.genome_id         = gf.genome_id;
                        meta.genome_type       = gf.genome_type;
                        meta.shard_id          = gf.shard_id;
                        meta.genome_length     = gf.genome_length;
                        meta.n_contigs         = gf.n_contigs;
                        meta.gc_pct_x100       = static_cast<uint16_t>(gf.gc_pct_x100);
                        meta.completeness_x10  = gf.completeness_x10;
                        meta.contamination_x10 = gf.contamination_x10;
                        meta.oph_fingerprint   = gf.oph_fingerprint;
                        meta.date_added        = gf.date_added;
                        catalog_rows.push_back(meta);

                        gidx_infos.push_back({gf.genome_id,
                                              static_cast<ShardId>(gf.shard_id),
                                              gf.dir_index});

                        accession_pairs.emplace_back(accession, gf.genome_id);

                        std::array<float, 136> kmer4;
                        std::memcpy(kmer4.data(), gf.kmer4, sizeof(gf.kmer4));
                        kmer_pairs.emplace_back(gf.genome_id, kmer4);
                        source_row_indices.push_back(gf.input_row_index);

                        if (gf.taxonomy_len)
                            taxonomy_map.emplace(accession, taxonomy);

                        ++genomes_restored;
                    }
                }
                if (mf) std::fclose(mf);

                if (meta_ok && genomes_restored == ck_genome_count) {
                    // Truncate .gpk to the checkpoint byte offset
                    if (::truncate(gpk_path_.c_str(), static_cast<off_t>(ck_byte_offset)) != 0)
                        throw std::runtime_error("checkpoint resume: truncate failed: " +
                                                 std::string(std::strerror(errno)));
                    app_writer.open_append(gpk_path_);
                    app_writer.seek_to(ck_byte_offset);

                    next_genome_id    = static_cast<GenomeId>(ck_next_gid);
                    current_shard_id  = static_cast<ShardId>(ck_shard_id);
                    next_section_id   = ck_section_id;

                    std::vector<BuildRecord> remaining_pending;
                    std::vector<uint64_t> remaining_rows;
                    remaining_pending.reserve(pending.size() - genomes_restored);
                    remaining_rows.reserve(pending_input_rows.size() - genomes_restored);
                    for (size_t i = 0; i < pending.size(); ++i) {
                        if (processed_rows[i]) continue;
                        remaining_pending.push_back(std::move(pending[i]));
                        remaining_rows.push_back(pending_input_rows[i]);
                    }
                    pending = std::move(remaining_pending);
                    pending_input_rows = std::move(remaining_rows);
                    total_records = pending.size();
                    resuming = true;

                    spdlog::info("Resuming from checkpoint: {}/{} genomes, {} shards, offset={}",
                                 genomes_restored, original_total_records, current_shard_id, ck_byte_offset);
                } else {
                    spdlog::warn("Checkpoint metadata corrupted (read {}/{}), starting fresh",
                                 genomes_restored, ck_genome_count);
                    catalog_rows.clear();
                    gidx_infos.clear();
                    accession_pairs.clear();
                    taxonomy_map.clear();
                    kmer_pairs.clear();
                    source_row_indices.clear();
                    shard_id_to_section_id.clear();
                    toc = TocWriter{};
                    next_section_id  = 1;
                    current_shard_id = 0;
                    next_genome_id   = cfg.starting_genome_id;
                    total_records    = pending.size();
                }
            }
        }

        // ── meta.tsv sidecar ─────────────────────────────────────────────────
        std::filesystem::path meta_tsv_path =
            gpk_path_.parent_path() / (gpk_path_.stem().string() + ".meta.tsv");

        if (!resuming) {
            app_writer.create(gpk_path_);

            // Write FileHeader (128B)
            {
                FileHeader fhdr{};
                fhdr.magic         = GPK2_MAGIC;
                fhdr.version_major = FORMAT_MAJOR;
                fhdr.version_minor = FORMAT_MINOR;
                uint64_t t = static_cast<uint64_t>(std::time(nullptr));
                fhdr.file_uuid_lo  = t ^ 0xdeadbeefcafe0001ULL;
                fhdr.file_uuid_hi  = (t << 17) ^ 0x1234567890abcdefULL;
                fhdr.created_at_unix = t;
                fhdr.flags         = 0;
                std::memset(fhdr.reserved, 0, sizeof(fhdr.reserved));
                app_writer.append(&fhdr, sizeof(fhdr));
            }

            std::ofstream meta_out(meta_tsv_path);
            meta_out << "accession\tgenome_id";
            if (!pending.empty())
                for (const auto& [k, v] : pending[0].extra_fields)
                    meta_out << "\t" << k;
            meta_out << "\n";
            // meta_out is a fresh stream; it will be reopened in append mode below
        }

        // Open meta_out in append mode (works for both fresh and resume)
        std::ofstream meta_out(meta_tsv_path, std::ios::app);

        uint32_t date = days_since_epoch();

        // ── IO writer thread: drains frozen shards and writes them to disk ─────
        struct WriteTask {
            std::future<FrozenShard> fut;
            uint64_t section_id;
            ShardId  shard_id;
            uint32_t n_genomes;
            size_t   catalog_start;
        };
        const size_t            write_q_max = 2;
        std::queue<WriteTask>   write_q;
        std::mutex              write_q_mx;
        std::condition_variable write_q_cv;
        std::condition_variable write_q_space_cv;
        bool                    writer_done = false;

        // ── Checkpoint write lambda ───────────────────────────────────────────
        FILE* ckpt_meta_file = std::fopen(ckpt_meta_path.c_str(),
                                          resuming ? "ab" : "wb");
        if (!ckpt_meta_file)
            throw std::runtime_error("Cannot open checkpoint meta file: " + ckpt_meta_path.string());

        struct DrainResult {
            ShardId  shard_id;
            uint64_t section_id;
            uint64_t file_offset;
            uint64_t compressed_size;
            uint32_t n_genomes;
            size_t   catalog_start;
        };

        auto write_checkpoint = [&](const DrainResult& dr) {
            // Append shard block to meta bin
            ShardCkptHdr shdr{};
            shdr.shard_id        = dr.shard_id;
            shdr.n_genomes       = dr.n_genomes;
            shdr.section_id      = dr.section_id;
            shdr.file_offset     = dr.file_offset;
            shdr.compressed_size = dr.compressed_size;
            std::fwrite(&shdr, sizeof(shdr), 1, ckpt_meta_file);

            for (size_t i = dr.catalog_start; i < dr.catalog_start + dr.n_genomes; ++i) {
                const auto& cm  = catalog_rows[i];
                const auto& acc = accession_pairs[i].first;
                std::string_view tax;
                auto tit = taxonomy_map.find(acc);
                if (tit != taxonomy_map.end()) tax = tit->second;

                const auto& kp = kmer_pairs[i];

                GenomeCkptFixed gf{};
                gf.genome_id          = cm.genome_id;
                gf.input_row_index    = source_row_indices[i];
                gf.shard_id           = cm.shard_id;
                gf.dir_index          = gidx_infos[i].dir_index;
                gf.oph_fingerprint    = cm.oph_fingerprint;
                gf.genome_length      = cm.genome_length;
                gf.n_contigs          = cm.n_contigs;
                gf.gc_pct_x100        = cm.gc_pct_x100;
                gf.completeness_x10   = cm.completeness_x10;
                gf.contamination_x10  = cm.contamination_x10;
                gf.date_added         = cm.date_added;
                gf.genome_type        = cm.genome_type;
                gf.accession_len      = static_cast<uint16_t>(acc.size());
                gf.taxonomy_len       = static_cast<uint16_t>(tax.size());
                std::memcpy(gf.kmer4, kp.second.data(), sizeof(gf.kmer4));

                std::fwrite(&gf, sizeof(gf), 1, ckpt_meta_file);
                if (gf.accession_len) std::fwrite(acc.data(),  1, gf.accession_len, ckpt_meta_file);
                if (gf.taxonomy_len)  std::fwrite(tax.data(),  1, gf.taxonomy_len,  ckpt_meta_file);
            }
            std::fflush(ckpt_meta_file);

            // Update scalar checkpoint (overwrite)
            std::ofstream ck(ckpt_path);
            ck << "genome_count="     << (dr.catalog_start + dr.n_genomes) << "\n"
               << "byte_offset="      << app_writer.current_offset()       << "\n"
               << "next_genome_id="   << next_genome_id                    << "\n"
               << "current_shard_id=" << current_shard_id                  << "\n"
               << "next_section_id="  << next_section_id                   << "\n";
        };

        // ── IO writer thread ──────────────────────────────────────────────────
        std::thread io_writer([&]() {
            while (true) {
                WriteTask wt;
                {
                    std::unique_lock lk(write_q_mx);
                    write_q_cv.wait(lk, [&]{ return !write_q.empty() || writer_done; });
                    if (write_q.empty()) break;
                    wt = std::move(write_q.front());
                    write_q.pop();
                }
                write_q_space_cv.notify_one();

                FrozenShard frozen = wt.fut.get();
                uint64_t shard_start = app_writer.current_offset();
                app_writer.append(frozen.bytes.data(), frozen.bytes.size());
                SectionDesc sd{};
                sd.type             = SEC_SHRD;
                sd.version          = 4;
                sd.flags            = 0;
                sd.section_id       = wt.section_id;
                sd.file_offset      = shard_start;
                sd.compressed_size  = static_cast<uint64_t>(frozen.bytes.size());
                sd.uncompressed_size = 0;
                sd.item_count       = wt.n_genomes;
                sd.aux0             = wt.shard_id;
                sd.aux1             = 0;
                std::memset(sd.checksum, 0, sizeof(sd.checksum));
                toc.add_section(sd);
                DrainResult dr{wt.shard_id, wt.section_id, shard_start,
                               static_cast<uint64_t>(frozen.bytes.size()),
                               wt.n_genomes, wt.catalog_start};
                write_checkpoint(dr);
            }
        });

        size_t catalog_start_for_current_shard = catalog_rows.size();

        std::unique_ptr<ShardWriter> shard_writer;
        uint32_t genomes_in_current_shard = 0;

        auto open_shard = [&]() {
            catalog_start_for_current_shard = catalog_rows.size();
            shard_writer = std::make_unique<ShardWriter>(
                current_shard_id, current_shard_id, cfg.shard_cfg);
            ++current_shard_id;
            genomes_in_current_shard = 0;
        };

        auto launch_shard_freeze = [&]() {
            if (!shard_writer || shard_writer->n_genomes() == 0) return;
            ShardId flushed_id = current_shard_id - 1;
            uint64_t sid       = next_section_id++;
            uint32_t ng        = static_cast<uint32_t>(shard_writer->n_genomes());
            shard_id_to_section_id[flushed_id] = sid; // pre-assign for GIDX
            WriteTask wt;
            wt.fut = std::async(std::launch::async,
                [sw = std::move(shard_writer)]() mutable { return sw->freeze(); });
            wt.section_id    = sid;
            wt.shard_id      = flushed_id;
            wt.n_genomes     = ng;
            wt.catalog_start = catalog_start_for_current_shard;
            {
                std::unique_lock lk(write_q_mx);
                write_q_space_cv.wait(lk, [&]{ return write_q.size() < write_q_max; });
                write_q.push(std::move(wt));
            }
            write_q_cv.notify_one();
            shard_writer.reset();
        };

        // ── Worker pool: io_threads persistent threads ─────────────────────────
        const size_t n_workers = std::max(size_t(1), cfg.io_threads);
        const size_t sort_buf  = n_workers * 4;

        struct ChunkItem {
            BuildRecord record;
            GenomeId    genome_id;
            uint64_t    input_row_index;
            FastaStats  stats;
            std::string fasta;
            OPHDualSketchResult sketch;                    // populated when build_sketch (single-k)
            std::vector<OPHDualSketchResult> sketches_mk;  // populated when multi_k_sketch
        };

        // Task queue (all tasks submitted upfront; poison-pill sentinel at end)
        // fd=-1 on poison pill; producer opens file + fadvise(WILLNEED) before queuing
        // so the kernel starts NFS prefetch 4*n_workers genomes ahead of workers.
        struct Task { BuildRecord* record; GenomeId gid; uint64_t input_row_index; int fd = -1; };
        std::queue<Task>        task_q;
        std::mutex              task_mx;
        std::condition_variable task_cv;

        // Completion queue — bounded to prevent OOM on large archives.
        const size_t            done_q_max = n_workers * 2;
        struct Done { std::optional<ChunkItem> item; };
        std::queue<Done>        done_q;
        std::mutex              done_mx;
        std::condition_variable done_cv;
        std::condition_variable done_push_cv;

        const size_t total = total_records;

        // Streaming producer: feeds tasks lazily to cap queue at 4*n_workers entries
        const size_t task_q_max = n_workers * 4;

        std::thread producer([&]() {
            for (size_t i = 0; i < total_records; ++i) {
                // Open + fadvise before acquiring the lock so kernel starts NFS
                // prefetch while we wait for queue space. Workers read from this fd.
                int fd = ::open(pending[i].file_path.c_str(), O_RDONLY);
#ifdef POSIX_FADV_WILLNEED
                if (fd >= 0)
                    ::posix_fadvise(fd, 0, 0, POSIX_FADV_WILLNEED | POSIX_FADV_SEQUENTIAL);
#endif
                std::unique_lock lk(task_mx);
                task_cv.wait(lk, [&]{ return task_q.size() < task_q_max; });
                task_q.push({&pending[i], next_genome_id++, pending_input_rows[i], fd});
                lk.unlock();
                task_cv.notify_one();
            }
            for (size_t i = 0; i < n_workers; ++i) {
                std::unique_lock lk(task_mx);
                task_cv.wait(lk, [&]{ return task_q.size() < task_q_max; });
                task_q.push({nullptr, 0, 0, -1});
                lk.unlock();
                task_cv.notify_one();
            }
        });

        std::vector<std::thread> workers;
        workers.reserve(n_workers);
        for (size_t i = 0; i < n_workers; ++i) {
            workers.emplace_back([&]() {
                while (true) {
                    Task t;
                    {
                        std::unique_lock lk(task_mx);
                        task_cv.wait(lk, [&]{ return !task_q.empty(); });
                        t = task_q.front();
                        task_q.pop();
                    }
                    task_cv.notify_one();
                    if (!t.record) return; // poison pill

                    Done d;
                    try {
                        std::string fasta = (t.fd >= 0)
                            ? decompress_gz_fd(t.fd, t.record->file_path)
                            : decompress_gz(t.record->file_path);
                        FastaStats stats  = compute_fasta_stats(fasta);
                        OPHDualSketchResult sk;
                        std::vector<OPHDualSketchResult> sks_mk;
                        if (cfg.build_sketch) {
                            const uint64_t seed1 = cfg.sketch_seed;
                            const uint64_t seed2 = cfg.sketch_seed + 1;
                            if (multi_k_sketch) {
                                for (int k : cfg.sketch_kmer_sizes) {
                                    sks_mk.push_back(sketch_oph_dual_from_buffer(
                                        fasta.data(), fasta.size(),
                                        k, cfg.sketch_size, cfg.sketch_syncmer_s,
                                        seed1, seed2));
                                }
                            } else {
                                sk = sketch_oph_dual_from_buffer(
                                    fasta.data(), fasta.size(),
                                    cfg.sketch_kmer_size, cfg.sketch_size,
                                    cfg.sketch_syncmer_s,
                                    seed1, seed2);
                            }
                        }
                        d.item = ChunkItem{*t.record, t.gid, t.input_row_index,
                                           stats, std::move(fasta), std::move(sk), std::move(sks_mk)};
                    } catch (const std::exception& ex) {
                        spdlog::warn("Skipping {}: {}", t.record->accession, ex.what());
                    }
                    {
                        std::unique_lock lk(done_mx);
                        done_push_cv.wait(lk, [&]{ return done_q.size() < done_q_max; });
                        done_q.push(std::move(d));
                    }
                    done_cv.notify_one();
                }
            });
        }
        task_cv.notify_all();

        // Staging-buffer flush: sort a full shard-sized batch before freezing.
        size_t n_done = 0, n_failed = 0;

        // When taxonomy_group is enabled: per-taxon buckets (flushed when bucket bytes >= threshold).
        // When disabled: single flat buffer (same as before, but sorted by kmer NN chain or oph).
        struct TaxonBucket { std::vector<ChunkItem> items; uint64_t raw_bytes = 0; };
        std::unordered_map<std::string, TaxonBucket> taxon_buckets;
        uint64_t taxon_total_bytes = 0; // sum of raw_bytes across all buckets
        // Global cap: flush the largest bucket whenever total in-memory genome
        // data exceeds this limit. Prevents unbounded accumulation across many
        // small-genus buckets that never individually hit max_shard_size_bytes.
        const uint64_t taxon_global_cap = 16ULL << 30; // 16 GB
        std::vector<ChunkItem> staging_buffer; // used when !cfg.taxonomy_group
        staging_buffer.reserve(sort_buf * 4);
        uint64_t staging_raw_bytes = 0;

        auto flush_staging_buf = [&](std::vector<ChunkItem>& buf) {
            if (buf.empty()) return;
            // Sort by kmer NN chain (if enabled) or oph_fingerprint
            if (cfg.kmer_nn_sort && buf.size() > 1) {
                auto order = greedy_nn_chain(buf);
                std::vector<ChunkItem> sorted;
                sorted.reserve(order.size());
                for (size_t idx : order) sorted.push_back(std::move(buf[idx]));
                buf = std::move(sorted);
            } else {
                std::sort(buf.begin(), buf.end(), [](const ChunkItem& a, const ChunkItem& b) {
                    return a.stats.oph_fingerprint < b.stats.oph_fingerprint;
                });
            }
            open_shard();
            for (auto& item : buf) {
                GenomeMeta meta{};
                meta.genome_id         = item.genome_id;
                meta.genome_type       = static_cast<uint32_t>(item.record.genome_type);
                meta.shard_id          = current_shard_id - 1;
                meta.genome_length     = item.record.genome_length > 0 ? item.record.genome_length : item.stats.genome_length;
                meta.n_contigs         = item.record.n_contigs > 0     ? item.record.n_contigs     : item.stats.n_contigs;
                meta.gc_pct_x100       = item.stats.gc_pct_x100;
                meta.completeness_x10  = static_cast<uint16_t>(item.record.completeness  * 10.0f);
                meta.contamination_x10 = static_cast<uint16_t>(item.record.contamination * 10.0f);
                meta.oph_fingerprint   = item.stats.oph_fingerprint;
                meta.date_added        = date;
                shard_writer->add_genome(item.genome_id, item.stats.oph_fingerprint,
                                         item.fasta.data(), item.fasta.size());
                gidx_infos.push_back({item.genome_id,
                                      static_cast<ShardId>(current_shard_id - 1),
                                      genomes_in_current_shard++});
                source_row_indices.push_back(item.input_row_index);
                catalog_rows.push_back(meta);
                accession_pairs.emplace_back(item.record.accession, item.genome_id);
                kmer_pairs.emplace_back(item.genome_id, item.stats.kmer4_profile);
                if (cfg.build_sketch) {
                    if (skch_writer_mk && !item.sketches_mk.empty()) {
                        std::vector<std::vector<uint16_t>> sigs1_per_k;
                        std::vector<std::vector<uint16_t>> sigs2_per_k;
                        std::vector<uint32_t>              n_real_bins_per_k;
                        std::vector<std::vector<uint64_t>> masks_per_k;
                        for (const auto& sk : item.sketches_mk) {
                            const size_t n = sk.signature1.size();
                            std::vector<uint16_t> sig1_16(n);
                            std::vector<uint16_t> sig2_16(n);
                            for (size_t si = 0; si < n; ++si) {
                                sig1_16[si] = static_cast<uint16_t>(sk.signature1[si] >> 16);
                                sig2_16[si] = static_cast<uint16_t>(sk.signature2[si] >> 16);
                            }
                            sigs1_per_k.push_back(std::move(sig1_16));
                            sigs2_per_k.push_back(std::move(sig2_16));
                            n_real_bins_per_k.push_back(sk.n_real_bins);
                            masks_per_k.push_back(sk.real_bins_bitmask);
                        }
                        skch_writer_mk->add(item.genome_id, item.stats.genome_length,
                                            sigs1_per_k, sigs2_per_k,
                                            n_real_bins_per_k, masks_per_k);
                    } else if (skch_writer && !item.sketch.signature1.empty()) {
                        const auto& sk = item.sketch;
                        const size_t n = sk.signature1.size();
                        std::vector<uint16_t> sig1_16(n);
                        std::vector<uint16_t> sig2_16(n);
                        for (size_t si = 0; si < n; ++si) {
                            sig1_16[si] = static_cast<uint16_t>(sk.signature1[si] >> 16);
                            sig2_16[si] = static_cast<uint16_t>(sk.signature2[si] >> 16);
                        }
                        skch_writer->add(item.genome_id, sig1_16, sig2_16,
                                         sk.n_real_bins, sk.genome_length,
                                         sk.real_bins_bitmask);
                    }
                }
                for (const auto& [k, v] : item.record.extra_fields)
                    if (k == "taxonomy") { taxonomy_map.emplace(item.record.accession, v); break; }
                parse_fasta_contig_accessions(item.fasta, [&](std::string_view contig_acc) {
                    cidx_writer.add(contig_acc, static_cast<uint32_t>(item.genome_id));
                });
                meta_out << item.record.accession << "\t" << item.genome_id;
                for (const auto& [k, v] : item.record.extra_fields) meta_out << "\t" << v;
                meta_out << "\n";
            }
            launch_shard_freeze();
            buf.clear();
            staging_raw_bytes = 0;
        };

        // Drain completion queue until all genomes processed
        for (size_t remaining = total; remaining > 0; --remaining) {
            Done d;
            {
                std::unique_lock lk(done_mx);
                done_cv.wait(lk, [&]{ return !done_q.empty(); });
                d = std::move(done_q.front());
                done_q.pop();
            }
            done_push_cv.notify_one();
            if (d.item) {
                ++n_done;
                if (cfg.taxonomy_group) {
                    // Find taxonomy string from extra_fields
                    std::string_view tax;
                    for (const auto& [k, v] : d.item->record.extra_fields)
                        if (k == "taxonomy") { tax = v; break; }
                    std::string key = extract_taxonomy_bucket(tax, cfg.taxonomy_rank.empty()
                                                              ? 'g' : cfg.taxonomy_rank[0]);
                    auto& bucket = taxon_buckets[key];
                    const uint64_t genome_bytes = d.item->fasta.size();
                    bucket.raw_bytes += genome_bytes;
                    taxon_total_bytes += genome_bytes;
                    bucket.items.push_back(std::move(*d.item));
                    if (bucket.raw_bytes >= cfg.shard_cfg.max_shard_size_bytes) {
                        taxon_total_bytes -= bucket.raw_bytes;
                        flush_staging_buf(bucket.items);
                        bucket.raw_bytes = 0;
                    } else if (taxon_total_bytes >= taxon_global_cap) {
                        // Find and flush the largest bucket to stay under cap
                        auto it = std::max_element(taxon_buckets.begin(), taxon_buckets.end(),
                            [](const auto& a, const auto& b){ return a.second.raw_bytes < b.second.raw_bytes; });
                        taxon_total_bytes -= it->second.raw_bytes;
                        flush_staging_buf(it->second.items);
                        it->second.raw_bytes = 0;
                    }
                } else {
                    staging_raw_bytes += d.item->fasta.size();
                    staging_buffer.push_back(std::move(*d.item));
                    if (staging_raw_bytes >= cfg.shard_cfg.max_shard_size_bytes)
                        flush_staging_buf(staging_buffer);
                }
            } else {
                ++n_failed;
            }
            if (cfg.verbose || n_done % 50000 == 0)
                spdlog::info("Build: {}/{} genomes ({:.1f}%) | {} shards | {} failed",
                             n_done, original_total_records,
                             100.0 * n_done / std::max<size_t>(size_t(1), original_total_records),
                             current_shard_id, n_failed);
        }
        // Final flush: remaining items in all buffers
        if (cfg.taxonomy_group) {
            for (auto& [key, bucket] : taxon_buckets)
                if (!bucket.items.empty()) { taxon_total_bytes -= bucket.raw_bytes; flush_staging_buf(bucket.items); bucket.raw_bytes = 0; }
        } else {
            if (!staging_buffer.empty()) flush_staging_buf(staging_buffer);
        }

        producer.join();
        for (auto& w : workers) w.join();
        pending.clear();
        pending_input_rows.clear();
        {
            std::unique_lock lk(write_q_mx);
            writer_done = true;
        }
        write_q_cv.notify_one();
        io_writer.join();

        std::fclose(ckpt_meta_file);
        ckpt_meta_file = nullptr;

        spdlog::info("Build: {}/{} genomes ({:.1f}%) | {} shards | {} failed",
                     n_done, original_total_records,
                     100.0 * n_done / std::max<size_t>(size_t(1), original_total_records),
                     current_shard_id, n_failed);
        spdlog::info("Wrote {} shards, {} genomes ({} failed)", current_shard_id,
                     catalog_rows.size(), n_failed);

        // Flush all shard writes to NFS, then write all metadata to a local
        // temp file. NFS write-back + ENOSPC can silently corrupt metadata,
        // leaving a garbage TailLocator that makes the archive unreadable by
        // merge_archives. Writing locally guarantees integrity.
        app_writer.flush();

        const uint64_t meta_base = app_writer.current_offset();
        auto meta_tmp_path = gpk_path_.parent_path() / "gpk_bld_meta_XXXXXX";
        std::string meta_tmp_str = meta_tmp_path.string();
        std::vector<char> bld_meta_tmp(meta_tmp_str.begin(), meta_tmp_str.end());
        bld_meta_tmp.push_back('\0');
        {
            int fd = ::mkstemp(bld_meta_tmp.data());
            if (fd < 0) throw std::runtime_error("Cannot create temp metadata file");
            ::close(fd);
        }
        AppendWriter mw;
        mw.create(bld_meta_tmp.data());
        mw.seek_to(meta_base);  // so section.file_offset values match the NFS file

        // Write CATL section
        uint64_t catalog_root_id = 0;
        {
            CatalogSectionWriter csw;
            for (const auto& m : catalog_rows)
                csw.add(m);
            SectionDesc catl_sd = csw.finalize(mw, next_section_id++);
            catalog_root_id = catl_sd.section_id;
            toc.add_section(catl_sd);
        }

        // Write GIDX section (O(1) genome_id lookup index)
        {
            GidxWriter gw;
            for (uint64_t i = 0; i < gidx_infos.size(); ++i) {
                const auto& gi = gidx_infos[i];
                auto it = shard_id_to_section_id.find(gi.shard_id);
                uint32_t sec_id = (it != shard_id_to_section_id.end())
                    ? static_cast<uint32_t>(it->second) : 0;
                gw.add(gi.genome_id, sec_id, gi.dir_index, i);
            }
            SectionDesc gidx_sd = gw.finalize(mw, next_section_id++);
            toc.add_section(gidx_sd);
        }

        // Write ACCX section
        uint64_t accession_root_id = 0;
        {
            AccessionIndexWriter aiw;
            for (const auto& [accession, genome_id] : accession_pairs)
                aiw.add(accession, genome_id);
            SectionDesc accx_sd = aiw.finalize(mw, next_section_id++);
            accession_root_id = accx_sd.section_id;
            toc.add_section(accx_sd);
        }

        // Write CIDX section (contig accession → genome_id index)
        if (cfg.build_cidx && cidx_writer.size() > 0) {
            spdlog::info("Writing CIDX: {} contig accessions", cidx_writer.size());
            SectionDesc cidx_sd = cidx_writer.finalize(mw, next_section_id++, /*batch_id=*/0);
            toc.add_section(cidx_sd);
        }

        // Write TAXN section (if taxonomy data available)
        if (!taxonomy_map.empty()) {
            TaxonomyIndexWriter tiw;
            for (const auto& [accession, taxonomy] : taxonomy_map)
                tiw.add(accession, taxonomy);
            SectionDesc taxn_sd = tiw.finalize(mw, next_section_id++);
            toc.add_section(taxn_sd);

            std::unordered_map<std::string, std::array<float, 136>> acc_kmer_profiles;
            if (!kmer_pairs.empty()) {
                std::unordered_map<GenomeId, const std::string*> gid_to_acc;
                gid_to_acc.reserve(accession_pairs.size());
                for (const auto& [acc, gid] : accession_pairs)
                    gid_to_acc[gid] = &acc;
                acc_kmer_profiles.reserve(kmer_pairs.size());
                for (const auto& [gid, prof] : kmer_pairs) {
                    auto it = gid_to_acc.find(gid);
                    if (it != gid_to_acc.end())
                        acc_kmer_profiles[*it->second] = prof;
                }
            }

            TxdbWriter txw;
            for (const auto& [accession, taxonomy] : taxonomy_map)
                txw.add(accession, taxonomy);
            if (!acc_kmer_profiles.empty())
                txw.set_kmer_profiles(acc_kmer_profiles);
            SectionDesc txdb_sd = txw.finalize(mw, next_section_id++);
            toc.add_section(txdb_sd);
        }

        // Write KMRX section
        if (!kmer_pairs.empty()) {
            KmrxWriter kw;
            for (auto& [gid, prof] : kmer_pairs)
                kw.add(gid, prof);
            SectionDesc kmrx_sd = kw.finalize(mw, next_section_id++);
            toc.add_section(kmrx_sd);
        }

        // Write SKCH section (OPH sketches)
        if (cfg.build_sketch) {
            if (skch_writer_mk) {
                std::string klist;
                for (size_t i = 0; i < cfg.sketch_kmer_sizes.size(); ++i) {
                    if (i) klist += ',';
                    klist += std::to_string(cfg.sketch_kmer_sizes[i]);
                }
                spdlog::info("Writing SKCH (v2 multi-k): {} sketches (k=[{}], m={})",
                             catalog_rows.size(), klist, cfg.sketch_size);
                SectionDesc skch_sd = skch_writer_mk->finalize(mw, next_section_id++);
                toc.add_section(skch_sd);
                skch_writer_mk.reset();
            } else if (skch_writer) {
                spdlog::info("Writing SKCH: {} sketches (k={}, m={})",
                             catalog_rows.size(), cfg.sketch_kmer_size, cfg.sketch_size);
                SectionDesc skch_sd = skch_writer->finalize(mw, next_section_id++);
                toc.add_section(skch_sd);
                skch_writer.reset();
            }
        }

        // Write TOC + TailLocator to local file
        toc.finalize(mw,
                     /*generation=*/1,
                     /*live_count=*/catalog_rows.size(),
                     /*total_count=*/catalog_rows.size(),
                     /*prev_toc_offset=*/0,
                     catalog_root_id,
                     accession_root_id,
                     /*tombstone_root_id=*/0);
        mw.flush();  // local flush — always reliable

        // Copy metadata from local file to NFS using O_SYNC in 64MB chunks
        {
            const uint64_t meta_size = mw.current_offset() - meta_base;
            const size_t   CHUNK     = 64 * 1024 * 1024;
            int local_fd = ::open(bld_meta_tmp.data(), O_RDONLY);
            int nfs_fd   = ::open(gpk_path_.c_str(), O_WRONLY | O_SYNC);
            if (local_fd < 0 || nfs_fd < 0) {
                if (local_fd >= 0) ::close(local_fd);
                if (nfs_fd   >= 0) ::close(nfs_fd);
                throw std::runtime_error("Cannot open files for metadata copy");
            }
            std::vector<uint8_t> buf(CHUNK);
            uint64_t copied = 0;
            while (copied < meta_size) {
                size_t  chunk = static_cast<size_t>(std::min<uint64_t>(CHUNK, meta_size - copied));
                off_t   soff  = static_cast<off_t>(meta_base + copied);
                ssize_t nr    = ::pread(local_fd, buf.data(), chunk, soff);
                if (nr <= 0) {
                    ::close(local_fd); ::close(nfs_fd);
                    throw std::runtime_error("Local metadata read failed");
                }
                size_t  wr = 0;
                int     retries = 0;
                while (wr < static_cast<size_t>(nr)) {
                    ssize_t nw = ::pwrite(nfs_fd, buf.data() + wr,
                                          static_cast<size_t>(nr) - wr,
                                          soff + static_cast<off_t>(wr));
                    if (nw < 0) {
                        if ((errno == ENOSPC || errno == EIO) && retries < 60) {
                            spdlog::warn("NFS metadata write ENOSPC, retry {}…", ++retries);
                            std::this_thread::sleep_for(std::chrono::seconds(5));
                            continue;
                        }
                        ::close(local_fd); ::close(nfs_fd);
                        throw std::runtime_error("NFS metadata write failed: " +
                                                 std::string(std::strerror(errno)));
                    }
                    wr += static_cast<size_t>(nw);
                }
                copied += static_cast<uint64_t>(nr);
            }
            ::fsync(nfs_fd);
            ::close(local_fd);
            ::close(nfs_fd);
        }
        ::unlink(bld_meta_tmp.data());

        // Remove checkpoint files — build completed successfully
        std::error_code ec;
        std::filesystem::remove(ckpt_path, ec);
        std::filesystem::remove(ckpt_meta_path, ec);

        spdlog::info("genopack archive written: {}", gpk_path_.string());
    }
};

ArchiveBuilder::ArchiveBuilder(const std::filesystem::path& dir, Config cfg)
    : impl_(std::make_unique<Impl>(dir, cfg))
{}
ArchiveBuilder::~ArchiveBuilder() = default;

void ArchiveBuilder::add_from_tsv(const std::filesystem::path& tsv_path) {
    impl_->add_from_tsv(tsv_path);
}
void ArchiveBuilder::add(const BuildRecord& rec) { impl_->add(rec); }
void ArchiveBuilder::finalize() { impl_->finalize(); }

} // namespace genopack
