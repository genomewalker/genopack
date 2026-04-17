#include <genopack/archive.hpp>
#include <genopack/catalog.hpp>
#include <genopack/cidx.hpp>
#include <genopack/gidx.hpp>
#include <genopack/kmrx.hpp>
#include <genopack/shard.hpp>
#include <genopack/taxn.hpp>
#include <genopack/txdb.hpp>
#include <genopack/util.hpp>
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/toc.hpp>
#include <genopack/accx.hpp>
#include <genopack/tombstone.hpp>
#include <tuple>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace genopack {

// ── ArchiveAppender::Impl ─────────────────────────────────────────────────────

struct ArchiveAppender::Impl {
    std::filesystem::path     gpk_path_;
    std::vector<BuildRecord>  pending_add_;
    std::vector<std::string>  pending_remove_accessions_;
    std::vector<GenomeId>     tombstone_ids;
    Toc                       existing_toc_;
    uint64_t                  current_toc_offset_ = 0;
    GenomeId                  next_genome_id_ = 1;
    ShardWriterConfig         shard_cfg_;

    explicit Impl(const std::filesystem::path& path) {
        std::filesystem::path gpk = path;
        if (!std::filesystem::exists(gpk) && gpk.extension() != ".gpk")
            gpk = std::filesystem::path(path.string() + ".gpk");

        if (!std::filesystem::exists(gpk) || !std::filesystem::is_regular_file(gpk))
            throw std::runtime_error("Archive not found: " + path.string());

        MmapFileReader mmap;
        mmap.open(gpk);
        if (mmap.size() < sizeof(FileHeader))
            throw std::runtime_error("File too small to be a .gpk: " + gpk.string());

        auto* fh = mmap.ptr_at<FileHeader>(0);
        if (fh->magic != GPK2_MAGIC)
            throw std::runtime_error("Not a .gpk file: " + gpk.string());

        gpk_path_ = gpk;
        existing_toc_ = TocReader::read(mmap);
        const auto* tail = mmap.ptr_at<TailLocator>(mmap.size() - sizeof(TailLocator));
        current_toc_offset_ = tail->toc_offset;
        next_genome_id_ = existing_toc_.header.total_genome_count + 1;
    }

    void commit() {
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
            // Build accession set from existing ACCX sections for dedup check
            MmapFileReader mmap_dedup;
            mmap_dedup.open(gpk_path_);
            std::unordered_set<std::string> existing_accessions;
            for (auto* sd : existing_toc_.find_by_type(SEC_ACCX)) {
                AccessionIndexReader reader;
                reader.open(mmap_dedup.data(), sd->file_offset, sd->compressed_size);
                reader.scan([&](std::string_view acc, GenomeId) {
                    existing_accessions.emplace(acc);
                });
            }

            struct GenomeMeta1 {
                BuildRecord record;
                GenomeId    genome_id;
                FastaStats  stats;
            };

            std::vector<GenomeMeta1> work;
            work.reserve(pending_add_.size());

            size_t n_dup_accession = 0;
            GenomeId gid = next_genome_id_;
            for (auto& r : pending_add_) {
                if (existing_accessions.count(r.accession)) {
                    ++n_dup_accession;
                    continue;
                }
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
            if (n_dup_accession > 0)
                spdlog::info("Dedup: skipped {} duplicate accessions already in archive", n_dup_accession);

            std::sort(work.begin(), work.end(), [](const GenomeMeta1& a, const GenomeMeta1& b) {
                return a.stats.oph_fingerprint < b.stats.oph_fingerprint;
            });

            // ── Sequence dedup ──────────────────────────────────────────────────
            // Step 1: build set of (oph_fingerprint, genome_length) from existing catalog
            struct SeqKey { uint64_t oph; uint64_t len;
                bool operator==(const SeqKey& o) const { return oph==o.oph && len==o.len; } };
            struct SeqKeyHash { size_t operator()(const SeqKey& k) const {
                return std::hash<uint64_t>{}(k.oph) ^ (std::hash<uint64_t>{}(k.len) * 0x9e3779b97f4a7c15ULL); } };
            std::unordered_set<SeqKey, SeqKeyHash> existing_seqs;
            {
                MergedCatalogReader cat;
                for (auto* sd : existing_toc_.find_by_type(SEC_CATL))
                    cat.add_fragment(mmap_dedup.data(), sd->file_offset, sd->compressed_size);
                cat.scan([&](const GenomeMeta& m) {
                    if (!m.is_deleted())
                        existing_seqs.insert({m.oph_fingerprint, m.genome_length});
                    return true;
                });
            }

            // Step 2: remove new genomes already in catalog (cross-catalog dedup)
            {
                size_t n_before = work.size();
                work.erase(std::remove_if(work.begin(), work.end(), [&](const GenomeMeta1& g) {
                    return existing_seqs.count({g.stats.oph_fingerprint, g.stats.genome_length}) > 0;
                }), work.end());
                size_t n_removed = n_before - work.size();
                if (n_removed > 0)
                    spdlog::info("Dedup: removed {} cross-catalog sequence duplicates", n_removed);
            }

            // Step 3: intra-batch dedup (identical sequences sort adjacent)
            {
                size_t n_before = work.size();
                size_t out = 0;
                for (size_t i = 0; i < work.size(); ) {
                    size_t j = i + 1;
                    while (j < work.size()
                           && work[j].stats.oph_fingerprint == work[i].stats.oph_fingerprint
                           && work[j].stats.genome_length   == work[i].stats.genome_length)
                        ++j;
                    size_t keep = i;
                    for (size_t k = i + 1; k < j; ++k)
                        if (work[k].record.completeness > work[keep].record.completeness)
                            keep = k;
                    if (out != keep) work[out] = std::move(work[keep]);
                    ++out;
                    i = j;
                }
                work.resize(out);
                if (work.size() < n_before)
                    spdlog::info("Dedup: removed {} intra-batch duplicates", n_before - work.size());
            }

            // Determine next shard_id from existing sections
            uint32_t new_shard_id = 0;
            for (const auto& sd : existing_toc_.sections) {
                if (sd.type == SEC_SHRD)
                    new_shard_id = std::max(new_shard_id, static_cast<uint32_t>(sd.aux0) + 1);
            }

            std::vector<GenomeMeta>                         new_catalog_rows;
            std::vector<std::pair<std::string, GenomeId>>   new_accessions;
            std::vector<std::pair<std::string, std::string>> new_taxonomies; // (accession, taxonomy)
            std::vector<std::pair<GenomeId, std::array<float,136>>> new_kmer_profiles;
            // (genome_id, new_shard_section_id, dir_index_in_shard)
            std::vector<std::tuple<GenomeId,uint32_t,uint32_t>> new_gidx_entries;
            std::unordered_map<uint32_t, uint64_t> shard_id_to_new_section_id; // shard_id → section_id

            // For contig index
            std::vector<std::pair<uint64_t, uint32_t>> new_cidx_hashes; // (acc_hash, genome_id)

            std::unique_ptr<ShardWriter> shard_writer;
            uint32_t current_shard_id = new_shard_id;

            auto flush_shard = [&]() {
                if (!shard_writer) return;
                uint64_t shard_offset = shard_writer->finalize(writer);
                uint64_t shard_size   = writer.current_offset() - shard_offset;
                SectionDesc sd{};
                sd.type            = SEC_SHRD;
                sd.version         = 4;
                sd.section_id      = next_section_id++;
                sd.file_offset     = shard_offset;
                sd.compressed_size = shard_size;
                sd.item_count      = static_cast<uint64_t>(shard_writer->n_genomes());
                sd.aux0            = current_shard_id;
                new_toc.add_section(sd);
                shard_id_to_new_section_id[current_shard_id] = sd.section_id;
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

                uint32_t dir_idx = shard_writer->n_genomes(); // index before add
                shard_writer->add_genome(gw.genome_id, gw.stats.oph_fingerprint,
                                         fasta.data(), fasta.size());
                new_catalog_rows.push_back(meta);
                new_accessions.emplace_back(gw.record.accession, gw.genome_id);
                new_gidx_entries.emplace_back(gw.genome_id, current_shard_id, dir_idx);
                new_kmer_profiles.emplace_back(gw.genome_id, gw.stats.kmer4_profile);
                for (const auto& [k, v] : gw.record.extra_fields)
                    if (k == "taxonomy") { new_taxonomies.emplace_back(gw.record.accession, v); break; }
                parse_fasta_contig_accessions(fasta, [&](std::string_view contig_acc) {
                    new_cidx_hashes.emplace_back(cidx_hash(contig_acc),
                                                 static_cast<uint32_t>(gw.genome_id));
                });
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

            // Write GIDX fragment for new genomes
            if (!new_gidx_entries.empty()) {
                // Build catl_row_index for new genomes (appended after existing)
                uint64_t base_catl_row = existing_toc_.header.total_genome_count;
                GidxWriter gw_idx;
                for (uint64_t i = 0; i < new_gidx_entries.size(); ++i) {
                    auto [gid_e, shard_id_e, dir_idx_e] = new_gidx_entries[i];
                    auto it = shard_id_to_new_section_id.find(shard_id_e);
                    if (it != shard_id_to_new_section_id.end())
                        gw_idx.add(gid_e, static_cast<uint32_t>(it->second), dir_idx_e,
                                   base_catl_row + i);
                }
                SectionDesc gidx_sd = gw_idx.finalize(writer, next_section_id++);
                new_toc.add_section(gidx_sd);
            }

            // Write TAXN fragment for new genomes (if taxonomy available)
            if (!new_taxonomies.empty()) {
                TaxonomyIndexWriter tiw;
                for (const auto& [acc, tax] : new_taxonomies)
                    tiw.add(acc, tax);
                SectionDesc taxn_sd = tiw.finalize(writer, next_section_id++);
                new_toc.add_section(taxn_sd);

                TxdbWriter txw;
                for (const auto& [acc, tax] : new_taxonomies)
                    txw.add(acc, tax);
                SectionDesc txdb_sd = txw.finalize(writer, next_section_id++);
                new_toc.add_section(txdb_sd);
            }

            // Write CIDX fragment for new genomes
            if (!new_cidx_hashes.empty()) {
                CidxWriter cw;
                for (const auto& [h, gid_e] : new_cidx_hashes)
                    cw.add_hash(h, gid_e);
                SectionDesc cidx_sd = cw.finalize(writer, next_section_id++, /*batch_id=*/0);
                new_toc.add_section(cidx_sd);
            }

            // Write KMRX fragment for new genomes
            if (!new_kmer_profiles.empty()) {
                KmrxWriter kw;
                for (const auto& [gid_e, prof] : new_kmer_profiles)
                    kw.add(gid_e, prof);
                SectionDesc kmrx_sd = kw.finalize(writer, next_section_id++);
                new_toc.add_section(kmrx_sd);
            }
            // HNSW is a global ANN graph and cannot be incrementally updated.
            // After append:
            //   - similar ACC_NEW works immediately (ACC_NEW has a KMRX profile and
            //     can query the existing HNSW to find older genomes).
            //   - similar ACC_BASE does NOT return ACC_NEW until the index is rebuilt.
            // Run 'genopack reindex --hnsw' to rebuild the full similarity index.
            // This is intentional: cheap incremental writes, explicit expensive reindex.

            next_genome_id_ = gid;
        }

        // ── Process tombstones ────────────────────────────────────────────────
        if (!pending_remove_accessions_.empty() || !tombstone_ids.empty()) {
            TombstoneWriter tw;
            uint64_t n_removed = 0;

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
            }

            for (GenomeId gid : tombstone_ids) {
                tw.add(gid);
                ++n_removed;
            }
            tombstone_ids.clear();

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
                         current_toc_offset_,
                         catalog_root_id, accession_root_id, tombstone_root_id);

        writer.flush();
        spdlog::info("genopack committed: gen {}, {} live genomes",
                     existing_toc_.header.generation + 1, new_live_count);
    }
};

// ── Public interface ──────────────────────────────────────────────────────────

ArchiveAppender::ArchiveAppender(const std::filesystem::path& dir)
    : impl_(std::make_unique<Impl>(dir))
{}
ArchiveAppender::~ArchiveAppender() = default;

void ArchiveAppender::add_from_tsv(const std::filesystem::path& tsv_path) {
    size_t before = impl_->pending_add_.size();
    auto records = parse_tsv_records(tsv_path);
    for (auto& r : records)
        impl_->pending_add_.push_back(std::move(r));
    spdlog::info("Appender loaded {} records from {}",
                 impl_->pending_add_.size() - before, tsv_path.string());
}

void ArchiveAppender::add(const BuildRecord& rec) {
    impl_->pending_add_.push_back(rec);
}

void ArchiveAppender::remove(GenomeId id) {
    impl_->tombstone_ids.push_back(id);
}

void ArchiveAppender::remove_by_accession(std::string_view accession) {
    impl_->pending_remove_accessions_.emplace_back(accession);
    spdlog::info("Queued tombstone: {}", accession);
}

void ArchiveAppender::commit() {
    impl_->commit();
}

} // namespace genopack
