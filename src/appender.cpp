#include <genopack/archive.hpp>
#include <genopack/catalog.hpp>
#include <genopack/shard.hpp>
#include <genopack/util.hpp>
#include <genopack/format_v2.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/toc.hpp>
#include <genopack/accx.hpp>
#include <genopack/tombstone.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cstring>
#include <ctime>
#include <fstream>
#include <stdexcept>
#include <unordered_map>

namespace genopack {

// ── Load meta.tsv ─────────────────────────────────────────────────────────────

static void load_meta_tsv(const std::filesystem::path& path,
                           std::unordered_map<std::string, GenomeId>& acc_map,
                           GenomeId& max_genome_id)
{
    std::ifstream f(path);
    if (!f) return;
    std::string line;
    std::getline(f, line); // skip header
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto tab1 = line.find('\t');
        if (tab1 == std::string::npos) continue;
        auto tab2 = line.find('\t', tab1 + 1);
        std::string acc    = line.substr(0, tab1);
        std::string gidstr = line.substr(tab1 + 1, tab2 - tab1 - 1);
        try {
            GenomeId gid = std::stoull(gidstr);
            acc_map[acc] = gid;
            if (gid > max_genome_id) max_genome_id = gid;
        } catch (...) {}
    }
}

// ── ArchiveAppender::Impl ─────────────────────────────────────────────────────

struct ArchiveAppender::Impl {
    // v1 fields
    std::filesystem::path    archive_dir;
    std::vector<std::string> tombstone_accessions;
    std::vector<GenomeId>    tombstone_ids;
    std::vector<BuildRecord> pending;

    // v2 fields
    std::filesystem::path gpk_path_;
    std::vector<BuildRecord>  pending_add_;
    std::vector<std::string>  pending_remove_accessions_;
    Toc                       existing_toc_;
    bool                      is_v2_ = false;
    GenomeId                  next_genome_id_ = 1;
    ShardWriterConfig         shard_cfg_;

    explicit Impl(const std::filesystem::path& path) {
        // Try v2 .gpk file first
        std::filesystem::path gpk = path;
        if (!std::filesystem::exists(gpk) && gpk.extension() != ".gpk")
            gpk = std::filesystem::path(path.string() + ".gpk");

        if (std::filesystem::exists(gpk) && std::filesystem::is_regular_file(gpk)) {
            MmapFileReader mmap;
            mmap.open(gpk);
            if (mmap.size() >= sizeof(FileHeader)) {
                auto* fh = mmap.ptr_at<FileHeader>(0);
                if (fh->magic == GPK2_MAGIC) {
                    gpk_path_ = gpk;
                    existing_toc_ = TocReader::read(mmap);
                    is_v2_ = true;
                    next_genome_id_ = existing_toc_.header.total_genome_count + 1;
                    return;
                }
            }
        }

        // Fall back to v1 directory
        if (std::filesystem::is_directory(path)) {
            if (!std::filesystem::exists(path / "MANIFEST.bin"))
                throw std::runtime_error("Not a genopack archive: " + path.string());
            archive_dir = path;
            is_v2_ = false;
            return;
        }

        throw std::runtime_error("Archive not found: " + path.string());
    }

    // ── v1 commit ─────────────────────────────────────────────────────────────

    void commit_v1() {
        // ── 1. Load existing manifest ────────────────────────────────────────
        ManifestHeader manifest{};
        {
            std::ifstream mf(archive_dir / "MANIFEST.bin", std::ios::binary);
            if (!mf) throw std::runtime_error("Cannot read MANIFEST.bin");
            mf.read(reinterpret_cast<char*>(&manifest), sizeof(manifest));
            if (manifest.magic != GPKM_MAGIC)
                throw std::runtime_error("Invalid MANIFEST.bin magic");
        }

        // ── 2. Load existing catalog rows ────────────────────────────────────
        std::vector<GenomeMeta> all_rows;
        {
            auto catalog_path = archive_dir / "catalog.gpkc";
            if (std::filesystem::exists(catalog_path)) {
                CatalogReader cr;
                cr.open(catalog_path);
                cr.scan([&](const GenomeMeta& m) { all_rows.push_back(m); return true; });
                cr.close();
            }
        }

        // ── 3. Load meta.tsv → accession map + max genome_id ────────────────
        std::unordered_map<std::string, GenomeId> acc_map;
        GenomeId max_id = 0;
        load_meta_tsv(archive_dir / "meta.tsv", acc_map, max_id);
        GenomeId next_id = max_id + 1;

        // ── 4. Apply tombstones ──────────────────────────────────────────────
        for (const auto& acc : tombstone_accessions) {
            auto it = acc_map.find(acc);
            if (it == acc_map.end()) {
                spdlog::warn("Tombstone: accession not found: {}", acc);
                continue;
            }
            GenomeId gid = it->second;
            for (auto& row : all_rows) {
                if (row.genome_id == gid) {
                    row.flags |= GenomeMeta::FLAG_DELETED;
                    break;
                }
            }
        }
        for (GenomeId gid : tombstone_ids) {
            for (auto& row : all_rows) {
                if (row.genome_id == gid) {
                    row.flags |= GenomeMeta::FLAG_DELETED;
                    break;
                }
            }
        }

        if (!tombstone_accessions.empty() || !tombstone_ids.empty()) {
            std::ofstream tlog(archive_dir / "tombstones.txt", std::ios::app);
            for (const auto& a : tombstone_accessions) tlog << a << "\n";
            for (GenomeId id : tombstone_ids) tlog << id << "\n";
            spdlog::info("Tombstoned {} accessions, {} ids",
                         tombstone_accessions.size(), tombstone_ids.size());
        }

        // ── 5. Process pending records ───────────────────────────────────────
        ShardId next_shard_id = manifest.n_shards;

        struct GenomeWork {
            BuildRecord record;
            FastaStats  stats;
            std::string fasta;
            GenomeMeta  meta;
        };

        std::vector<GenomeWork> work;
        work.reserve(pending.size());

        for (auto& r : pending) {
            GenomeWork gw;
            gw.record = std::move(r);
            try {
                gw.fasta = decompress_gz(gw.record.file_path);
            } catch (const std::exception& ex) {
                spdlog::warn("Skipping {}: {}", gw.record.accession, ex.what());
                continue;
            }
            gw.stats = compute_fasta_stats(gw.fasta);

            gw.meta = GenomeMeta{};
            gw.meta.genome_id         = next_id++;
            gw.meta._reserved0        = 0;
            gw.meta.shard_id          = INVALID_SHARD_ID;
            gw.meta.genome_length     = gw.record.genome_length > 0 ? gw.record.genome_length : gw.stats.genome_length;
            gw.meta.n_contigs         = gw.record.n_contigs > 0     ? gw.record.n_contigs     : gw.stats.n_contigs;
            gw.meta.gc_pct_x100       = gw.stats.gc_pct_x100;
            gw.meta.completeness_x10  = static_cast<uint16_t>(gw.record.completeness  * 10.0f);
            gw.meta.contamination_x10 = static_cast<uint16_t>(gw.record.contamination * 10.0f);
            gw.meta.oph_fingerprint   = gw.stats.oph_fingerprint;
            gw.meta.date_added        = days_since_epoch();
            gw.meta.flags             = 0;

            work.push_back(std::move(gw));
        }
        pending.clear();

        // ── 6. Sort by oph_fingerprint ───────────────────────────────────────
        std::sort(work.begin(), work.end(), [](const GenomeWork& a, const GenomeWork& b) {
            return a.meta.oph_fingerprint < b.meta.oph_fingerprint;
        });

        // ── 7. Write new shard(s) ────────────────────────────────────────────
        std::vector<std::pair<std::string, GenomeId>> new_accessions;
        std::vector<GenomeMeta>                       new_rows;

        if (!work.empty()) {
            std::filesystem::create_directories(archive_dir / "shards");

            ShardWriterConfig shard_cfg{};
            std::unique_ptr<ShardWriterV1> sw;
            size_t sw_compressed = 0;

            auto open_shard = [&]() {
                char fname[64];
                snprintf(fname, sizeof(fname), "shard_%05u.gpks", next_shard_id);
                sw = std::make_unique<ShardWriterV1>(
                    archive_dir / "shards" / fname,
                    next_shard_id, next_shard_id, shard_cfg);
                ++next_shard_id;
                sw_compressed = 0;
            };
            auto flush_shard = [&]() {
                if (sw) { sw->finalize(); sw.reset(); sw_compressed = 0; }
            };

            for (auto& gw : work) {
                if (!sw || sw_compressed >= shard_cfg.max_shard_size_bytes) {
                    flush_shard();
                    open_shard();
                }
                gw.meta.shard_id = next_shard_id - 1;
                sw->add_genome(gw.meta.genome_id, gw.meta.oph_fingerprint,
                               gw.fasta.data(), gw.fasta.size());
                sw_compressed = sw->compressed_size();

                new_rows.push_back(gw.meta);
                new_accessions.emplace_back(gw.record.accession, gw.meta.genome_id);
            }
            flush_shard();

            spdlog::info("Appender wrote {} new shard(s), {} genomes",
                         next_shard_id - manifest.n_shards, new_rows.size());
        }

        // ── 8. Merge + rewrite catalog ────────────────────────────────────────
        for (const auto& r : new_rows)
            all_rows.push_back(r);

        std::sort(all_rows.begin(), all_rows.end(), [](const GenomeMeta& a, const GenomeMeta& b) {
            return a.oph_fingerprint < b.oph_fingerprint;
        });

        {
            CatalogWriter cw(archive_dir / "catalog.gpkc");
            for (const auto& m : all_rows)
                cw.add(m);
            cw.finalize();
        }

        // ── 9. Append to meta.tsv ─────────────────────────────────────────────
        if (!new_accessions.empty()) {
            bool meta_exists = std::filesystem::exists(archive_dir / "meta.tsv");
            std::ofstream mf(archive_dir / "meta.tsv", std::ios::app);
            if (!meta_exists)
                mf << "accession\tgenome_id\n";
            for (const auto& [acc, gid] : new_accessions)
                mf << acc << "\t" << gid << "\n";
        }

        // ── 10. Write new MANIFEST.bin ────────────────────────────────────────
        uint64_t n_live = std::count_if(all_rows.begin(), all_rows.end(),
            [](const GenomeMeta& m) { return !m.is_deleted(); });

        ManifestHeader new_manifest{};
        new_manifest.magic          = GPKM_MAGIC;
        new_manifest.version        = FORMAT_VERSION;
        new_manifest.generation     = manifest.generation + 1;
        new_manifest.created_at     = static_cast<uint64_t>(std::time(nullptr));
        new_manifest.n_shards       = next_shard_id;
        new_manifest._reserved0     = 0;
        new_manifest.n_genomes      = static_cast<uint64_t>(all_rows.size());
        new_manifest.n_genomes_live = n_live;

        {
            std::ofstream mf_out(archive_dir / "MANIFEST.bin",
                                 std::ios::binary | std::ios::trunc);
            if (!mf_out) throw std::runtime_error("Cannot write MANIFEST.bin");
            mf_out.write(reinterpret_cast<const char*>(&new_manifest), sizeof(new_manifest));

            for (uint32_t i = 0; i < next_shard_id; ++i) {
                ShardDescriptor sd{};
                sd.shard_id   = i;
                sd.cluster_id = i;
                char fname[64];
                snprintf(fname, sizeof(fname), "shard_%05u.gpks", i);
                std::memcpy(sd.filename, fname, std::min(sizeof(fname), sizeof(sd.filename)));
                auto fpath = archive_dir / "shards" / fname;
                if (std::filesystem::exists(fpath))
                    sd.file_size = std::filesystem::file_size(fpath);
                mf_out.write(reinterpret_cast<const char*>(&sd), sizeof(sd));
            }
        }

        spdlog::info("Commit complete: gen={}, shards={}, genomes={}, live={}",
                     new_manifest.generation, next_shard_id,
                     all_rows.size(), n_live);
    }

    // ── v2 commit ─────────────────────────────────────────────────────────────

    void commit_v2() {
        AppendWriter writer;
        writer.open_append(gpk_path_);

        TocWriter new_toc;
        for (const auto& sd : existing_toc_.sections)
            new_toc.add_section(sd);

        uint64_t next_section_id   = existing_toc_.next_section_id();
        uint64_t new_live_count    = existing_toc_.header.live_genome_count;
        uint64_t new_total_count   = existing_toc_.header.total_genome_count;
        uint64_t catalog_root_id   = existing_toc_.header.catalog_root_section_id;
        uint64_t accession_root_id = existing_toc_.header.accession_root_section_id;
        uint64_t tombstone_root_id = existing_toc_.header.tombstone_root_section_id;

        // ── Process new genomes ───────────────────────────────────────────────
        if (!pending_add_.empty()) {
            struct GenomeMeta1 {
                BuildRecord record;
                GenomeId    genome_id;
                FastaStats  stats;
            };

            std::vector<GenomeMeta1> work;
            work.reserve(pending_add_.size());

            GenomeId gid = next_genome_id_;
            for (auto& r : pending_add_) {
                std::string fasta;
                try { fasta = decompress_gz(r.file_path); }
                catch (const std::exception& ex) {
                    spdlog::warn("Skipping {}: {}", r.accession, ex.what());
                    continue;
                }
                FastaStats stats = compute_fasta_stats(fasta);
                work.push_back({r, gid++, stats});
            }
            pending_add_.clear();

            std::sort(work.begin(), work.end(), [](const GenomeMeta1& a, const GenomeMeta1& b) {
                return a.stats.oph_fingerprint < b.stats.oph_fingerprint;
            });

            // Determine next shard_id from existing sections
            uint32_t new_shard_id = 0;
            for (const auto& sd : existing_toc_.sections) {
                if (sd.type == SEC_SHRD)
                    new_shard_id = std::max(new_shard_id, static_cast<uint32_t>(sd.aux0) + 1);
            }

            std::vector<GenomeMeta>                       new_catalog_rows;
            std::vector<std::pair<std::string, GenomeId>> new_accessions;

            std::unique_ptr<ShardWriter> shard_writer;
            uint32_t current_shard_id = new_shard_id;

            auto flush_shard = [&]() {
                if (!shard_writer) return;
                uint64_t shard_offset = shard_writer->finalize(writer);
                uint64_t shard_size   = writer.current_offset() - shard_offset;
                SectionDesc sd{};
                sd.type            = SEC_SHRD;
                sd.section_id      = next_section_id++;
                sd.file_offset     = shard_offset;
                sd.compressed_size = shard_size;
                sd.item_count      = static_cast<uint64_t>(shard_writer->n_genomes());
                sd.aux0            = current_shard_id;
                new_toc.add_section(sd);
                shard_writer.reset();
            };

            auto open_shard = [&]() {
                shard_writer = std::make_unique<ShardWriter>(
                    current_shard_id, current_shard_id, shard_cfg_);
            };

            uint32_t date = days_since_epoch();
            for (auto& gw : work) {
                if (!shard_writer) open_shard();
                if (shard_writer->n_bytes_raw() >= shard_cfg_.max_shard_size_bytes) {
                    flush_shard();
                    ++current_shard_id;
                    open_shard();
                }

                std::string fasta;
                try { fasta = decompress_gz(gw.record.file_path); }
                catch (const std::exception& ex) {
                    spdlog::warn("Skipping {} on write pass: {}", gw.record.accession, ex.what());
                    continue;
                }

                GenomeMeta meta{};
                meta.genome_id         = gw.genome_id;
                meta.shard_id          = current_shard_id;
                meta.genome_length     = gw.record.genome_length > 0 ? gw.record.genome_length : gw.stats.genome_length;
                meta.n_contigs         = gw.record.n_contigs > 0     ? gw.record.n_contigs     : gw.stats.n_contigs;
                meta.gc_pct_x100       = gw.stats.gc_pct_x100;
                meta.completeness_x10  = static_cast<uint16_t>(gw.record.completeness * 10.0f);
                meta.contamination_x10 = static_cast<uint16_t>(gw.record.contamination * 10.0f);
                meta.oph_fingerprint   = gw.stats.oph_fingerprint;
                meta.date_added        = date;

                shard_writer->add_genome(gw.genome_id, gw.stats.oph_fingerprint,
                                         fasta.data(), fasta.size());
                new_catalog_rows.push_back(meta);
                new_accessions.emplace_back(gw.record.accession, gw.genome_id);
            }
            flush_shard();

            // Write new CATL fragment
            if (!new_catalog_rows.empty()) {
                CatalogSectionWriter csw;
                for (const auto& m : new_catalog_rows) csw.add(m);
                SectionDesc catl_sd = csw.finalize(writer, next_section_id++);
                new_toc.add_section(catl_sd);
                catalog_root_id  = catl_sd.section_id;
                new_live_count  += new_catalog_rows.size();
                new_total_count += new_catalog_rows.size();
            }

            // Write new ACCX fragment
            if (!new_accessions.empty()) {
                AccessionIndexWriter aiw;
                for (const auto& [acc, id2] : new_accessions) aiw.add(acc, id2);
                SectionDesc accx_sd = aiw.finalize(writer, next_section_id++);
                new_toc.add_section(accx_sd);
                accession_root_id = accx_sd.section_id;
            }

            next_genome_id_ = gid;
        }

        // ── Process tombstones ────────────────────────────────────────────────
        if (!pending_remove_accessions_.empty()) {
            MmapFileReader mmap;
            mmap.open(gpk_path_);

            std::unordered_map<std::string, GenomeId> acc_map;
            for (auto* sd : existing_toc_.find_by_type(SEC_ACCX)) {
                AccessionIndexReader reader;
                reader.open(mmap.data(), sd->file_offset, sd->compressed_size);
                reader.scan([&](std::string_view acc, GenomeId id2) {
                    acc_map[std::string(acc)] = id2;
                });
            }

            TombstoneWriter tw;
            uint64_t n_removed = 0;
            for (const auto& acc : pending_remove_accessions_) {
                auto it = acc_map.find(acc);
                if (it == acc_map.end()) {
                    spdlog::warn("remove: accession not found: {}", acc);
                    continue;
                }
                tw.add(it->second);
                ++n_removed;
            }
            pending_remove_accessions_.clear();

            if (n_removed > 0) {
                SectionDesc tomb_sd = tw.finalize(writer, next_section_id++);
                new_toc.add_section(tomb_sd);
                tombstone_root_id = tomb_sd.section_id;
                new_live_count = (new_live_count >= n_removed)
                    ? new_live_count - n_removed : 0;
            }
        }

        // ── Write new TOCB + TailLocator ──────────────────────────────────────
        new_toc.finalize(writer,
                         existing_toc_.header.generation + 1,
                         new_live_count, new_total_count,
                         existing_toc_.header.prev_toc_offset,
                         catalog_root_id, accession_root_id, tombstone_root_id);

        writer.flush();
        spdlog::info("genopack v2 committed: gen {}, {} live genomes",
                     existing_toc_.header.generation + 1, new_live_count);
    }
};

// ── Public interface ──────────────────────────────────────────────────────────

ArchiveAppender::ArchiveAppender(const std::filesystem::path& dir)
    : impl_(std::make_unique<Impl>(dir))
{}
ArchiveAppender::~ArchiveAppender() = default;

void ArchiveAppender::add_from_tsv(const std::filesystem::path& tsv_path) {
    if (impl_->is_v2_) {
        size_t before = impl_->pending_add_.size();
        auto records = parse_tsv_records(tsv_path);
        for (auto& r : records)
            impl_->pending_add_.push_back(std::move(r));
        spdlog::info("Appender loaded {} records from {}",
                     impl_->pending_add_.size() - before, tsv_path.string());
    } else {
        size_t before = impl_->pending.size();
        auto records = parse_tsv_records(tsv_path);
        for (auto& r : records)
            impl_->pending.push_back(std::move(r));
        spdlog::info("Appender loaded {} records from {}",
                     impl_->pending.size() - before, tsv_path.string());
    }
}

void ArchiveAppender::add(const BuildRecord& rec) {
    if (impl_->is_v2_)
        impl_->pending_add_.push_back(rec);
    else
        impl_->pending.push_back(rec);
}

void ArchiveAppender::remove(GenomeId id) {
    impl_->tombstone_ids.push_back(id);
}

void ArchiveAppender::remove_by_accession(std::string_view accession) {
    if (impl_->is_v2_) {
        impl_->pending_remove_accessions_.emplace_back(accession);
    } else {
        impl_->tombstone_accessions.emplace_back(accession);
    }
    spdlog::info("Queued tombstone: {}", accession);
}

void ArchiveAppender::commit() {
    if (impl_->is_v2_) {
        impl_->commit_v2();
    } else {
        impl_->commit_v1();
    }
}

} // namespace genopack
