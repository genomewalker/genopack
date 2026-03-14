#include <genopack/archive.hpp>
#include <genopack/shard.hpp>
#include <genopack/catalog.hpp>
#include <genopack/util.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <future>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <thread>
#include <vector>

namespace genopack {

// ── ArchiveBuilder::Impl ──────────────────────────────────────────────────────

struct ArchiveBuilder::Impl {
    std::filesystem::path archive_dir;
    Config                cfg;

    std::vector<BuildRecord> pending;
    GenomeId next_genome_id = 1;

    // Pass-1 result: metadata only, FASTA discarded after stats computation
    struct GenomeMeta1 {
        BuildRecord record;
        GenomeId    genome_id;
        FastaStats  stats;
    };

    explicit Impl(const std::filesystem::path& dir, Config c)
        : archive_dir(dir), cfg(c)
    {
        std::filesystem::create_directories(dir / "shards");
    }

    void add(const BuildRecord& rec) {
        pending.push_back(rec);
    }

    void add_from_tsv(const std::filesystem::path& tsv_path) {
        auto records = parse_tsv_records(tsv_path);
        for (auto& r : records)
            pending.push_back(std::move(r));
        spdlog::info("Loaded {} records from {}", pending.size(), tsv_path.string());
    }

    void finalize() {
        if (pending.empty()) { spdlog::warn("No records to build"); return; }

        // ── Pass 1: parallel stats computation ────────────────────────────────
        // Read each FASTA, compute MinHash + stats, immediately discard FASTA.
        // Only metadata (~300 bytes/genome) is retained.

        std::vector<GenomeMeta1> work;
        work.reserve(pending.size());

        size_t batch = std::max(size_t(1), cfg.io_threads);
        for (size_t i = 0; i < pending.size(); i += batch) {
            size_t end = std::min(i + batch, pending.size());
            std::vector<std::future<std::optional<GenomeMeta1>>> futures;
            futures.reserve(end - i);
            for (size_t j = i; j < end; ++j) {
                auto& r = pending[j];
                GenomeId gid = next_genome_id++;
                futures.push_back(std::async(std::launch::async,
                    [&r, gid]() -> std::optional<GenomeMeta1> {
                        std::string fasta;
                        try {
                            fasta = decompress_gz(r.file_path);
                        } catch (const std::exception& ex) {
                            spdlog::warn("Skipping {}: {}", r.accession, ex.what());
                            return std::nullopt;
                        }
                        FastaStats stats = compute_fasta_stats(fasta);
                        // fasta discarded here
                        return GenomeMeta1{r, gid, stats};
                    }));
            }
            for (auto& fut : futures) {
                auto result = fut.get();
                if (result) work.push_back(std::move(*result));
            }
            if (cfg.verbose)
                spdlog::info("Stats pass: {}/{}", std::min(i + batch, pending.size()), pending.size());
        }
        pending.clear();

        // ── Sort by MinHash (oph_fingerprint) ─────────────────────────────────
        std::sort(work.begin(), work.end(), [](const GenomeMeta1& a, const GenomeMeta1& b) {
            return a.stats.oph_fingerprint < b.stats.oph_fingerprint;
        });
        spdlog::info("MinHash sorted {} genomes", work.size());

        // ── Pass 2: re-read FASTAs in sorted order, write shards ──────────────
        // Only one shard-worth of FASTAs is in RAM at a time.

        std::ofstream meta_out(archive_dir / "meta.tsv");
        if (!work.empty()) {
            meta_out << "accession\tgenome_id";
            for (const auto& [k, v] : work[0].record.extra_fields)
                meta_out << "\t" << k;
            meta_out << "\n";
        }

        ShardId current_shard_id = 0;
        std::unique_ptr<ShardWriter> shard_writer;
        size_t shard_compressed_size = 0;
        std::vector<GenomeMeta> catalog_rows;
        catalog_rows.reserve(work.size());
        uint32_t date = days_since_epoch();

        auto flush_shard = [&]() {
            if (shard_writer) {
                shard_writer->finalize();
                shard_writer.reset();
                shard_compressed_size = 0;
            }
        };

        auto open_shard = [&]() {
            char fname[64];
            snprintf(fname, sizeof(fname), "shard_%05u.gpks", current_shard_id);
            shard_writer = std::make_unique<ShardWriter>(
                archive_dir / "shards" / fname, current_shard_id, current_shard_id, cfg.shard_cfg);
            ++current_shard_id;
        };

        size_t n_done = 0;
        for (auto& gw : work) {
            bool shard_full = shard_writer &&
                (shard_compressed_size >= cfg.shard_cfg.max_shard_size_bytes);
            if (!shard_writer || shard_full) { flush_shard(); open_shard(); }

            // Re-read FASTA (only one genome at a time in RAM)
            std::string fasta;
            try {
                fasta = decompress_gz(gw.record.file_path);
            } catch (const std::exception& ex) {
                spdlog::warn("Skipping {} on write pass: {}", gw.record.accession, ex.what());
                continue;
            }

            GenomeMeta meta{};
            meta.genome_id         = gw.genome_id;
            meta._reserved0        = 0;
            meta.shard_id          = current_shard_id - 1;
            meta.genome_length     = gw.record.genome_length > 0 ? gw.record.genome_length : gw.stats.genome_length;
            meta.n_contigs         = gw.record.n_contigs > 0     ? gw.record.n_contigs     : gw.stats.n_contigs;
            meta.gc_pct_x100       = gw.stats.gc_pct_x100;
            meta.completeness_x10  = static_cast<uint16_t>(gw.record.completeness  * 10.0f);
            meta.contamination_x10 = static_cast<uint16_t>(gw.record.contamination * 10.0f);
            meta.oph_fingerprint   = gw.stats.oph_fingerprint;
            meta.date_added        = date;

            shard_writer->add_genome(gw.genome_id, gw.stats.oph_fingerprint,
                                     fasta.data(), fasta.size());
            shard_compressed_size = shard_writer->compressed_size();
            catalog_rows.push_back(meta);

            meta_out << gw.record.accession << "\t" << gw.genome_id;
            for (const auto& [k, v] : gw.record.extra_fields) meta_out << "\t" << v;
            meta_out << "\n";

            if (cfg.verbose && ++n_done % 1000 == 0)
                spdlog::info("Write pass: {}/{}", n_done, work.size());
        }
        flush_shard();

        spdlog::info("Wrote {} shards, {} genomes", current_shard_id, catalog_rows.size());

        // Write catalog
        {
            CatalogWriter cw(archive_dir / "catalog.gpkc");
            for (const auto& m : catalog_rows)
                cw.add(m);
            cw.finalize();
        }

        write_manifest(current_shard_id, catalog_rows.size(), date);

        spdlog::info("genopack archive built: {}", archive_dir.string());
    }

    void write_manifest(uint32_t n_shards, uint64_t n_genomes, uint32_t /*date*/) {
        ManifestHeader hdr{};
        hdr.magic          = GPKM_MAGIC;
        hdr.version        = FORMAT_VERSION;
        hdr.generation     = 1;
        hdr.created_at     = static_cast<uint64_t>(std::time(nullptr));
        hdr.n_shards       = n_shards;
        hdr._reserved0     = 0;
        hdr.n_genomes      = n_genomes;
        hdr.n_genomes_live = n_genomes;

        std::ofstream f(archive_dir / "MANIFEST.bin", std::ios::binary | std::ios::trunc);
        if (!f) throw std::runtime_error("Cannot write MANIFEST.bin");
        f.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));

        for (uint32_t i = 0; i < n_shards; ++i) {
            ShardDescriptor sd{};
            sd.shard_id   = i;
            sd.cluster_id = i;
            char fname[64];
            snprintf(fname, sizeof(fname), "shard_%05u.gpks", i);
            std::strncpy(sd.filename, fname, sizeof(sd.filename) - 1);
            auto fpath = archive_dir / "shards" / fname;
            if (std::filesystem::exists(fpath))
                sd.file_size = std::filesystem::file_size(fpath);
            f.write(reinterpret_cast<const char*>(&sd), sizeof(sd));
        }
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
