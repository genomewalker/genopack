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

        const size_t total = pending.size();
        spdlog::info("Building archive: {} genomes, {} io_threads", total, cfg.io_threads);

        // ── Single-pass streaming build ────────────────────────────────────────
        // Process genomes in chunks of io_threads in parallel:
        //   read FASTA → compute OPH+stats → sort chunk by OPH → write to shard.
        // No global sort, no dedup (run 'genopack dedup' post-build).
        // One NFS read pass total.

        // Open AppendWriter on the output .gpk file
        AppendWriter app_writer;
        app_writer.create(gpk_path_);

        // Write FileHeader (128B)
        {
            FileHeader fhdr{};
            fhdr.magic         = GPK2_MAGIC;
            fhdr.version_major = FORMAT_V2_MAJOR;
            fhdr.version_minor = FORMAT_V2_MINOR;
            uint64_t t = static_cast<uint64_t>(std::time(nullptr));
            fhdr.file_uuid_lo  = t ^ 0xdeadbeefcafe0001ULL;
            fhdr.file_uuid_hi  = (t << 17) ^ 0x1234567890abcdefULL;
            fhdr.created_at_unix = t;
            fhdr.flags         = 0;
            std::memset(fhdr.reserved, 0, sizeof(fhdr.reserved));
            app_writer.append(&fhdr, sizeof(fhdr));
        }

        // meta.tsv sidecar
        std::filesystem::path meta_tsv_path =
            gpk_path_.parent_path() / (gpk_path_.stem().string() + ".meta.tsv");
        std::ofstream meta_out(meta_tsv_path);
        {
            // Write header from first record's extra_fields
            meta_out << "accession\tgenome_id";
            if (!pending.empty())
                for (const auto& [k, v] : pending[0].extra_fields)
                    meta_out << "\t" << k;
            meta_out << "\n";
        }

        TocWriter toc;
        uint64_t next_section_id = 1;

        ShardId current_shard_id = 0;
        std::unique_ptr<ShardWriter> shard_writer;
        std::vector<GenomeMeta> catalog_rows;
        catalog_rows.reserve(total);
        std::vector<std::pair<std::string, GenomeId>> accession_pairs;
        accession_pairs.reserve(total);
        uint32_t date = days_since_epoch();

        auto flush_shard = [&]() {
            if (!shard_writer || shard_writer->n_genomes() == 0)
                return;
            uint64_t shard_start = shard_writer->finalize(app_writer);
            uint64_t shard_end   = app_writer.current_offset();
            SectionDesc sd{};
            sd.type             = SEC_SHRD;
            sd.version          = 2;
            sd.flags            = 0;
            sd.section_id       = next_section_id++;
            sd.file_offset      = shard_start;
            sd.compressed_size  = shard_end - shard_start;
            sd.uncompressed_size = 0;
            sd.item_count       = shard_writer->n_genomes();
            sd.aux0             = current_shard_id - 1;
            sd.aux1             = 0;
            std::memset(sd.checksum, 0, sizeof(sd.checksum));
            toc.add_section(sd);
            shard_writer.reset();
        };

        auto open_shard = [&]() {
            shard_writer = std::make_unique<ShardWriter>(
                current_shard_id, current_shard_id, cfg.shard_cfg);
            ++current_shard_id;
        };

        // Chunk: read io_threads genomes in parallel, then sort chunk by OPH and write.
        // Chunk size = io_threads so we keep all NFS readers busy at once.
        const size_t chunk_sz = std::max(size_t(1), cfg.io_threads);

        // Result of one genome read
        struct ChunkItem {
            BuildRecord  record;
            GenomeId     genome_id;
            FastaStats   stats;
            std::string  fasta;
        };

        size_t n_done = 0, n_failed = 0;

        for (size_t i = 0; i < total; i += chunk_sz) {
            const size_t end = std::min(i + chunk_sz, total);

            // Parallel read + OPH
            std::vector<std::future<std::optional<ChunkItem>>> futures;
            futures.reserve(end - i);
            for (size_t j = i; j < end; ++j) {
                auto& r = pending[j];
                GenomeId gid = next_genome_id++;
                futures.push_back(std::async(std::launch::async,
                    [&r, gid]() -> std::optional<ChunkItem> {
                        std::string fasta;
                        try {
                            fasta = decompress_gz(r.file_path);
                        } catch (const std::exception& ex) {
                            spdlog::warn("Skipping {}: {}", r.accession, ex.what());
                            return std::nullopt;
                        }
                        FastaStats stats = compute_fasta_stats(fasta);
                        return ChunkItem{r, gid, stats, std::move(fasta)};
                    }));
            }

            std::vector<ChunkItem> chunk;
            chunk.reserve(end - i);
            for (auto& fut : futures) {
                auto res = fut.get();
                if (res) chunk.push_back(std::move(*res));
                else     ++n_failed;
            }

            // Local OPH sort for better intra-shard compression
            std::sort(chunk.begin(), chunk.end(), [](const ChunkItem& a, const ChunkItem& b) {
                return a.stats.oph_fingerprint < b.stats.oph_fingerprint;
            });

            // Write chunk to shard(s)
            for (auto& item : chunk) {
                if (!shard_writer) open_shard();

                GenomeMeta meta{};
                meta.genome_id         = item.genome_id;
                meta._reserved0        = 0;
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
                catalog_rows.push_back(meta);
                accession_pairs.emplace_back(item.record.accession, item.genome_id);

                meta_out << item.record.accession << "\t" << item.genome_id;
                for (const auto& [k, v] : item.record.extra_fields) meta_out << "\t" << v;
                meta_out << "\n";

                if (shard_writer->n_bytes_raw() >= cfg.shard_cfg.max_shard_size_bytes)
                    flush_shard();
            }

            n_done += chunk.size();
            if (cfg.verbose || n_done % 50000 == 0 || n_done == total)
                spdlog::info("Build: {}/{} genomes ({:.1f}%) | {} shards | {} failed",
                             n_done, total, 100.0 * n_done / total,
                             current_shard_id, n_failed);
        }
        pending.clear();
        flush_shard();

        spdlog::info("Wrote {} shards, {} genomes ({} failed)", current_shard_id, catalog_rows.size(), n_failed);

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
