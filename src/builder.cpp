#include <genopack/archive.hpp>
#include <genopack/accx.hpp>
#include <genopack/catalog.hpp>
#include <genopack/format_v2.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/shard.hpp>
#include <genopack/toc.hpp>
#include <genopack/util.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <ctime>
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
    std::filesystem::path archive_dir;   // base path (no extension)
    std::filesystem::path gpk_path_;     // output .gpk file
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

        // ── Pass 2: write v2 single-file archive ──────────────────────────────

        // Open AppendWriter on the output .gpk file
        AppendWriter app_writer;
        app_writer.create(gpk_path_);

        // Write FileHeader (128B)
        {
            FileHeader fhdr{};
            fhdr.magic         = GPK2_MAGIC;
            fhdr.version_major = FORMAT_V2_MAJOR;
            fhdr.version_minor = FORMAT_V2_MINOR;
            // Simple UUID: XOR of time and a counter
            uint64_t t = static_cast<uint64_t>(std::time(nullptr));
            fhdr.file_uuid_lo  = t ^ 0xdeadbeefcafe0001ULL;
            fhdr.file_uuid_hi  = (t << 17) ^ 0x1234567890abcdefULL;
            fhdr.created_at_unix = t;
            fhdr.flags         = 0;
            std::memset(fhdr.reserved, 0, sizeof(fhdr.reserved));
            app_writer.append(&fhdr, sizeof(fhdr));
        }

        // meta.tsv sidecar (alongside the .gpk)
        std::filesystem::path meta_tsv_path =
            gpk_path_.parent_path() / (gpk_path_.stem().string() + ".meta.tsv");
        std::ofstream meta_out(meta_tsv_path);
        if (!work.empty()) {
            meta_out << "accession\tgenome_id";
            for (const auto& [k, v] : work[0].record.extra_fields)
                meta_out << "\t" << k;
            meta_out << "\n";
        }

        TocWriter toc;
        uint64_t next_section_id = 1;

        ShardId current_shard_id = 0;
        std::unique_ptr<ShardWriterV2> shard_writer;
        std::vector<GenomeMeta> catalog_rows;
        catalog_rows.reserve(work.size());
        std::vector<std::pair<std::string, GenomeId>> accession_pairs;
        accession_pairs.reserve(work.size());
        uint32_t date = days_since_epoch();

        auto flush_shard = [&]() {
            if (!shard_writer || shard_writer->n_genomes() == 0)
                return;

            uint64_t shard_offset_before = app_writer.current_offset();
            uint64_t shard_start = shard_writer->finalize(app_writer);
            uint64_t shard_end   = app_writer.current_offset();

            SectionDesc sd{};
            sd.type             = SEC_SHRD;
            sd.version          = 2;
            sd.flags            = 0;
            sd.section_id       = next_section_id++;
            sd.file_offset      = shard_start;
            sd.compressed_size  = shard_end - shard_start;
            sd.uncompressed_size = 0;  // not tracked at this level
            sd.item_count       = shard_writer->n_genomes();
            sd.aux0             = current_shard_id - 1;  // shard_id
            sd.aux1             = 0;
            std::memset(sd.checksum, 0, sizeof(sd.checksum));
            toc.add_section(sd);

            shard_writer.reset();

            (void)shard_offset_before;
        };

        auto open_shard = [&]() {
            shard_writer = std::make_unique<ShardWriterV2>(
                current_shard_id, current_shard_id, cfg.shard_cfg);
            ++current_shard_id;
        };

        size_t n_done = 0;
        for (auto& gw : work) {
            // Open shard if needed
            if (!shard_writer) open_shard();

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
            // blob_offset/blob_len_cmp left as 0 — shard section has correct values

            shard_writer->add_genome(gw.genome_id, gw.stats.oph_fingerprint,
                                     fasta.data(), fasta.size());
            catalog_rows.push_back(meta);
            accession_pairs.emplace_back(gw.record.accession, gw.genome_id);

            meta_out << gw.record.accession << "\t" << gw.genome_id;
            for (const auto& [k, v] : gw.record.extra_fields) meta_out << "\t" << v;
            meta_out << "\n";

            if (cfg.verbose && ++n_done % 1000 == 0)
                spdlog::info("Write pass: {}/{}", n_done, work.size());

            // Check shard fullness after adding genome
            if (shard_writer->n_bytes_raw() >= cfg.shard_cfg.max_shard_size_bytes) {
                flush_shard();
            }
        }
        // Flush any remaining genomes in the current shard
        flush_shard();

        spdlog::info("Wrote {} shards, {} genomes", current_shard_id, catalog_rows.size());

        // Write CATL section
        uint64_t catalog_root_id = 0;
        {
            CatalogSectionWriter csw;
            for (const auto& m : catalog_rows)
                csw.add(m);
            SectionDesc catl_sd = csw.finalize(app_writer, next_section_id++);
            catalog_root_id = catl_sd.section_id;
            toc.add_section(catl_sd);
        }

        // Write ACCX section
        uint64_t accession_root_id = 0;
        {
            AccessionIndexWriter aiw;
            for (const auto& [accession, genome_id] : accession_pairs)
                aiw.add(accession, genome_id);
            SectionDesc accx_sd = aiw.finalize(app_writer, next_section_id++);
            accession_root_id = accx_sd.section_id;
            toc.add_section(accx_sd);
        }

        // Write TOC + TailLocator
        toc.finalize(app_writer,
                     /*generation=*/1,
                     /*live_count=*/catalog_rows.size(),
                     /*total_count=*/catalog_rows.size(),
                     /*prev_toc_offset=*/0,
                     catalog_root_id,
                     accession_root_id,
                     /*tombstone_root_id=*/0);

        app_writer.flush();

        spdlog::info("genopack v2 archive written: {}", gpk_path_.string());
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
