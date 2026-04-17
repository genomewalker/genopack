#include <genopack/archive.hpp>
#include <genopack/cidx.hpp>
#include <genopack/format.hpp>
#include <genopack/gidx.hpp>
#include <genopack/kmrx.hpp>
#include <genopack/skch.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/toc.hpp>
#include <genopack/accx.hpp>
#include <genopack/taxn.hpp>
#include <genopack/txdb.hpp>
#include <genopack/tombstone.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <cstring>
#include <future>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <utility>

namespace genopack {

// ── ShardBox ─────────────────────────────────────────────────────────────────
// ShardReader() = default is inline in shard.hpp, which forces the compiler
// to instantiate the destructor of unique_ptr<ShardReader::Impl> here — but
// ShardReader::Impl is incomplete in this TU. Workaround: zero-initialize
// aligned storage for a ShardReader, then call only the out-of-line methods
// open() and ~ShardReader(). A zero-initialized unique_ptr<Impl> is nullptr
// on all platforms, so open() (which does `if (!impl_) impl_ = make_unique<Impl>()`)
// initialises it correctly in shard.cpp where Impl is complete.

struct ShardBox {
    alignas(ShardReader) unsigned char storage_[sizeof(ShardReader)];
    bool live_ = false;

    ShardBox() { std::memset(storage_, 0, sizeof(storage_)); }

    ~ShardBox() {
        if (live_) reinterpret_cast<ShardReader*>(storage_)->~ShardReader();
    }

    ShardBox(const ShardBox&)            = delete;
    ShardBox& operator=(const ShardBox&) = delete;

    ShardReader& reader() {
        return *reinterpret_cast<ShardReader*>(storage_);
    }
    const ShardReader& reader() const {
        return *reinterpret_cast<const ShardReader*>(storage_);
    }

    void open(const uint8_t* base, uint64_t offset, uint64_t size) {
        reader().open(base, offset, size);
        live_ = true;
    }
};

// ── ArchiveReader::Impl ───────────────────────────────────────────────────────

struct ArchiveReader::Impl {
    MmapFileReader      mmap_;
    Toc                 toc_;
    MergedCatalogReader catalog_;

    // shard_id -> SectionDesc pointer (into toc_.sections)
    std::unordered_map<uint32_t, const SectionDesc*> shard_descs_;
    // lazily-opened shard readers keyed by logical shard_id
    mutable std::unordered_map<uint32_t, ShardBox> shards_;
    std::unordered_map<uint64_t, uint32_t> shard_section_to_id_;
    mutable std::mutex shard_open_mx_;

    // accession <-> genome_id maps built from all ACCX sections
    std::unordered_map<std::string, GenomeId> accession_map_;
    std::unordered_map<GenomeId, std::string> genome_accession_map_;
    std::unordered_map<GenomeId, const GenomeMeta*> genome_meta_map_;

    // accession -> taxonomy string (loaded lazily on first taxonomy query)
    mutable std::unordered_map<std::string, std::string> taxonomy_map_;
    void ensure_taxonomy_loaded() const;

    // TXDB section descriptors (section_id -> SectionDesc pointer into toc_.sections)
    std::unordered_map<uint64_t, const SectionDesc*> txdb_descs_;

    // KMRX readers — all sections (merged archives have one per original part)
    std::vector<KmrxReader> kmrx_readers_;
    bool                    has_kmrx_ = false;

    // GIDX reader for O(1) genome_id -> shard position lookup
    GidxReader gidx_;
    bool       has_gidx_ = false;

    // CIDX reader for contig accession -> genome_id lookup
    MergedCidxReader cidx_;

    // SKCH section descriptors (loaded lazily on first sketch_for call).
    struct SkchSectionDesc {
        uint64_t              file_offset;
        uint64_t              compressed_size;
        std::vector<uint32_t> kmer_sizes;  // peek_params() result; v1 has 1 element
        uint32_t              sketch_size; // from SectionDesc.aux0
    };
    std::vector<SkchSectionDesc>    skch_descs_;
    mutable std::vector<SkchReader> skch_readers_;
    mutable std::atomic<bool>       skch_loaded_ = false;

    // cached taxonomy tree (lazy-built)
    mutable bool                       tree_built_ = false;
    mutable std::optional<TaxonomyTree> cached_tree_;

    // merged tombstones from all TOMB sections
    std::vector<TombstoneReader> tombstones_;

    bool     open_        = false;
    uint64_t live_count_  = 0;
    uint64_t total_count_ = 0;
    uint64_t generation_  = 0;
    uint32_t n_shards_    = 0;

    // ── taxonomy tree (lazy) ──────────────────────────────────────────────────

    std::optional<TaxonomyTree> get_tree() const {
        if (tree_built_) return cached_tree_;
        tree_built_ = true;

        // Prefer highest section_id TXDB section (most recent generation)
        if (!txdb_descs_.empty()) {
            uint64_t best_id = 0;
            const SectionDesc* best_sd = nullptr;
            for (const auto& [sid, sd] : txdb_descs_) {
                if (sid > best_id) { best_id = sid; best_sd = sd; }
            }
            if (best_sd) {
                TxdbReader reader;
                reader.open(mmap_.data(), best_sd->file_offset, best_sd->compressed_size);
                cached_tree_ = reader.tree();
                return cached_tree_;
            }
        }

        cached_tree_ = std::nullopt;
        return std::nullopt;
    }

    void open_gpk(const std::filesystem::path& path) {
        mmap_.open(path);

        auto* fh = mmap_.ptr_at<FileHeader>(0);
        if (fh->magic != GPK2_MAGIC)
            throw std::runtime_error("Not a .gpk file");

        toc_ = TocReader::read(mmap_);

        // Load catalog fragments (all CATL sections)
        for (auto* sd : toc_.find_by_type(SEC_CATL)) {
            catalog_.add_fragment(mmap_.data(), sd->file_offset, sd->compressed_size);
            CatalogSectionReader cat;
            cat.open(mmap_.data(), sd->file_offset, sd->compressed_size);
            for (uint32_t i = 0; i < cat.n_rows(); ++i) {
                const GenomeMeta* row = cat.row_at(i);
                if (row) genome_meta_map_[row->genome_id] = row;
            }
        }

        // Index shard sections by shard_id (aux0)
        n_shards_ = 0;
        for (auto& sd : toc_.sections) {
            if (sd.type == SEC_SHRD) {
                shard_descs_[static_cast<uint32_t>(sd.aux0)] = &sd;
                shard_section_to_id_[sd.section_id] = static_cast<uint32_t>(sd.aux0);
                ++n_shards_;
            }
        }

        // Load GIDX section (highest section_id wins)
        {
            uint64_t best_id = 0;
            const SectionDesc* best_sd = nullptr;
            for (auto* sd : toc_.find_by_type(SEC_GIDX)) {
                if (sd->section_id > best_id) {
                    best_id = sd->section_id;
                    best_sd = sd;
                }
            }
            if (best_sd) {
                gidx_.open(mmap_.data(), best_sd->file_offset, best_sd->compressed_size);
                has_gidx_ = true;
            }
        }

        // Load CIDX sections (contig accession → genome_id)
        for (auto* sd : toc_.find_by_type(SEC_CIDX))
            cidx_.add_section(mmap_.data(), sd->file_offset, sd->compressed_size);

        // Load accession index from all ACCX sections
        for (auto* sd : toc_.find_by_type(SEC_ACCX)) {
            AccessionIndexReader reader;
            reader.open(mmap_.data(), sd->file_offset, sd->compressed_size);
            reader.scan([&](std::string_view acc, GenomeId gid) {
                accession_map_[std::string(acc)] = gid;
                genome_accession_map_[gid]       = std::string(acc);
            });
        }

        // Defer taxonomy loading — only loaded on first taxonomy_for_accession() call.
        // 4.7M × ~150-char taxonomy strings cause significant heap pressure when loaded
        // eagerly alongside DuckDB; lazy loading avoids fragmentation-induced corruption.
        // (Genome fetch, KMRX, and HNSW access do not need taxonomy_map_.)

        // Index TXDB sections by section_id
        for (auto* sd : toc_.find_by_type(SEC_TXDB)) {
            txdb_descs_[sd->section_id] = sd;
        }

        // Load tombstones from all TOMB sections
        for (auto* sd : toc_.find_by_type(SEC_TOMB)) {
            tombstones_.emplace_back();
            tombstones_.back().open(mmap_.data(), sd->file_offset, sd->compressed_size);
        }

        // Load all KMRX sections (merged archives have one per original part).
        // kmer_profile() searches all readers so merged archives work correctly.
        {
            for (auto* sd : toc_.find_by_type(SEC_KMRX)) {
                KmrxReader r;
                r.open(mmap_.data(), sd->file_offset, sd->compressed_size);
                kmrx_readers_.push_back(std::move(r));
            }
            has_kmrx_ = !kmrx_readers_.empty();
        }

        // Record SKCH section locations for lazy loading on first sketch_for call.
        // Use peek_params() to read kmer_sizes[] without full decompression.
        {
            for (auto* sd : toc_.find_by_type(SEC_SKCH)) {
                auto [ver, ks] = SkchReader::peek_params(
                    mmap_.data(), sd->file_offset, sd->compressed_size);
                skch_descs_.push_back({sd->file_offset, sd->compressed_size,
                                       std::move(ks),
                                       static_cast<uint32_t>(sd->aux0)});
            }
            if (!skch_descs_.empty())
                spdlog::info("SKCH: {} section(s) deferred", skch_descs_.size());
        }

        live_count_  = toc_.header.live_genome_count;
        total_count_ = toc_.header.total_genome_count;
        generation_  = toc_.header.generation;

        open_  = true;
        spdlog::info("genopack opened: {} genomes, {} shards, gen {}",
                     live_count_, n_shards_, generation_);
    }

    void open(const std::filesystem::path& path) {
        // Reset all state before loading the new archive so that reusing the same
        // ArchiveReader instance across multiple open() calls never leaks stale data.
        if (open_) {
            mmap_.close();
            toc_                  = Toc{};
            catalog_              = MergedCatalogReader{};
            shard_descs_.clear();
            shards_.clear();
            shard_section_to_id_.clear();
            accession_map_.clear();
            genome_accession_map_.clear();
            genome_meta_map_.clear();
            taxonomy_map_.clear();
            txdb_descs_.clear();
            gidx_                 = GidxReader{};
            has_gidx_             = false;
            kmrx_readers_.clear();
            has_kmrx_             = false;
            cidx_                 = MergedCidxReader{};
            skch_descs_.clear();
            skch_readers_.clear();
            skch_loaded_          = false;
            tombstones_.clear();
            tree_built_           = false;
            cached_tree_.reset();
            live_count_           = 0;
            total_count_          = 0;
            generation_           = 0;
            n_shards_             = 0;
            open_                 = false;
        }

        auto gpk = path;
        if (!std::filesystem::exists(gpk) && gpk.extension() != ".gpk")
            gpk = std::filesystem::path(path.string() + ".gpk");

        open_gpk(gpk);
    }

    // ── shard access ────────────────────────────────────────────────────

    const ShardReader& get_shard(uint32_t shard_id) const {
        auto desc_it = shard_descs_.find(shard_id);
        if (desc_it == shard_descs_.end())
            throw std::runtime_error("Shard " + std::to_string(shard_id) + " not found");
        std::lock_guard<std::mutex> lk(shard_open_mx_);
        auto it = shards_.find(shard_id);
        if (it == shards_.end()) {
            auto [inserted, _] = shards_.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(shard_id),
                std::forward_as_tuple());
            inserted->second.open(mmap_.data(), desc_it->second->file_offset, desc_it->second->compressed_size);
            it = inserted;
        }
        return it->second.reader();
    }

    const ShardReader& get_shard_by_section_id(uint64_t section_id) const {
        auto it = shard_section_to_id_.find(section_id);
        if (it == shard_section_to_id_.end())
            throw std::runtime_error("Shard section " + std::to_string(section_id) + " not found");
        return get_shard(it->second);
    }

    // ── tombstone check ───────────────────────────────────────────────────────

    bool is_deleted(GenomeId id) const {
        for (const auto& t : tombstones_)
            if (t.is_deleted(id)) return true;
        return false;
    }

    std::optional<ExtractedGenome> fetch_genome_impl(GenomeId id) const {
        auto meta_it = genome_meta_map_.find(id);
        const GenomeMeta* meta = (meta_it != genome_meta_map_.end()) ? meta_it->second : nullptr;
        if (!meta || is_deleted(meta->genome_id)) return std::nullopt;

        ExtractedGenome eg;
        eg.meta  = *meta;
        if (has_gidx_) {
            const GidxEntry* ge = gidx_.lookup(id);
            if (ge) {
                const auto& shard = get_shard_by_section_id(ge->shard_section_id);
                const GenomeDirEntry* de = shard.dir_entry(ge->dir_index);
                if (de && de->genome_id == id)
                    eg.fasta = shard.fetch_genome_at(ge->dir_index);
                else
                    eg.fasta = get_shard(meta->shard_id).fetch_genome(id);
            } else {
                eg.fasta = get_shard(meta->shard_id).fetch_genome(id);
            }
        } else {
            eg.fasta = get_shard(meta->shard_id).fetch_genome(id);
        }

        auto acc_it = genome_accession_map_.find(id);
        if (acc_it != genome_accession_map_.end())
            eg.accession = acc_it->second;

        return eg;
    }

    std::optional<std::string> fetch_sequence_slice_impl(GenomeId id,
                                                         uint64_t start,
                                                         uint64_t length) const
    {
        auto meta_it = genome_meta_map_.find(id);
        const GenomeMeta* meta = (meta_it != genome_meta_map_.end()) ? meta_it->second : nullptr;
        if (!meta || is_deleted(meta->genome_id)) return std::nullopt;
        if (start >= meta->genome_length || length == 0) return std::string{};
        uint64_t clamped = std::min<uint64_t>(length, meta->genome_length - start);

        if (has_gidx_) {
            const GidxEntry* ge = gidx_.lookup(id);
            if (ge) {
                const auto& shard = get_shard_by_section_id(ge->shard_section_id);
                const GenomeDirEntry* de = shard.dir_entry(ge->dir_index);
                if (de && de->genome_id == id)
                    return shard.fetch_sequence_slice_at(ge->dir_index, start, clamped);
            }
        }
        return get_shard(meta->shard_id).fetch_sequence_slice(id, start, clamped);
    }

    std::vector<ExtractedGenome> extract_impl(const ExtractQuery& q) const {
        auto ptrs = catalog_.filter(q);

        std::vector<ExtractedGenome> out;
        out.reserve(ptrs.size());

        if (has_gidx_) {
            struct IndexedMeta {
                const GenomeMeta* meta;
                uint32_t          dir_index;
            };
            std::unordered_map<uint32_t, std::vector<IndexedMeta>> by_shard;
            for (const auto* p : ptrs) {
                if (is_deleted(p->genome_id)) continue;
                const GidxEntry* ge = gidx_.lookup(p->genome_id);
                if (!ge) {
                    by_shard[p->shard_id].push_back({p, UINT32_MAX});
                    continue;
                }
                auto sid_it = shard_section_to_id_.find(ge->shard_section_id);
                uint32_t shard_id = (sid_it != shard_section_to_id_.end()) ? sid_it->second : p->shard_id;
                by_shard[shard_id].push_back({p, ge->dir_index});
            }

            for (auto& [shard_id, metas] : by_shard) {
                const auto& shard = get_shard(shard_id);
                for (const auto& im : metas) {
                    ExtractedGenome eg;
                    eg.meta = *im.meta;
                    if (im.dir_index != UINT32_MAX) {
                        const GenomeDirEntry* de = shard.dir_entry(im.dir_index);
                        eg.fasta = (de && de->genome_id == im.meta->genome_id)
                            ? shard.fetch_genome_at(im.dir_index)
                            : shard.fetch_genome(im.meta->genome_id);
                    } else {
                        eg.fasta = shard.fetch_genome(im.meta->genome_id);
                    }
                    auto acc_it = genome_accession_map_.find(im.meta->genome_id);
                    if (acc_it != genome_accession_map_.end())
                        eg.accession = acc_it->second;
                    out.push_back(std::move(eg));
                }
            }
        } else {
            std::unordered_map<uint32_t, std::vector<const GenomeMeta*>> by_shard;
            for (const auto* p : ptrs) {
                if (!is_deleted(p->genome_id))
                    by_shard[p->shard_id].push_back(p);
            }

            for (auto& [shard_id, metas] : by_shard) {
                const auto& shard = get_shard(shard_id);
                for (const auto* m : metas) {
                    ExtractedGenome eg;
                    eg.meta  = *m;
                    eg.fasta = shard.fetch_genome(m->genome_id);
                    auto acc_it = genome_accession_map_.find(m->genome_id);
                    if (acc_it != genome_accession_map_.end())
                        eg.accession = acc_it->second;
                    out.push_back(std::move(eg));
                }
            }
        }
        return out;
    }
};

// ── ArchiveReader public API ──────────────────────────────────────────────────

ArchiveReader::ArchiveReader()  : impl_(std::make_unique<Impl>()) {}
ArchiveReader::~ArchiveReader() = default;

void ArchiveReader::open(const std::filesystem::path& path) { impl_->open(path); }

void ArchiveReader::close() {
    impl_->mmap_.close();
    impl_->toc_                  = Toc{};
    impl_->catalog_              = MergedCatalogReader{};
    impl_->shard_descs_.clear();
    impl_->shards_.clear();
    impl_->shard_section_to_id_.clear();
    impl_->accession_map_.clear();
    impl_->genome_accession_map_.clear();
    impl_->genome_meta_map_.clear();
    impl_->taxonomy_map_.clear();
    impl_->txdb_descs_.clear();
    impl_->gidx_                 = GidxReader{};
    impl_->has_gidx_             = false;
    impl_->kmrx_readers_.clear();
    impl_->has_kmrx_             = false;
    impl_->cidx_                 = MergedCidxReader{};
    impl_->tombstones_.clear();
    impl_->tree_built_           = false;
    impl_->cached_tree_.reset();
    impl_->live_count_           = 0;
    impl_->total_count_          = 0;
    impl_->generation_           = 0;
    impl_->n_shards_             = 0;
    impl_->open_                 = false;
}

bool ArchiveReader::is_open() const { return impl_->open_; }

std::optional<GenomeMeta> ArchiveReader::genome_meta(GenomeId id) const {
    auto it = impl_->genome_meta_map_.find(id);
    const GenomeMeta* p = (it != impl_->genome_meta_map_.end()) ? it->second : nullptr;
    if (!p) return std::nullopt;
    return *p;
}

size_t ArchiveReader::count(const ExtractQuery& q) const {
    return impl_->catalog_.filter(q).size();
}

std::vector<GenomeMeta> ArchiveReader::filter_meta(const ExtractQuery& q) const {
    auto ptrs = impl_->catalog_.filter(q);
    std::vector<GenomeMeta> out;
    out.reserve(ptrs.size());
    for (const auto* p : ptrs) out.push_back(*p);
    return out;
}

std::optional<ExtractedGenome> ArchiveReader::fetch_genome(GenomeId id) const {
    return impl_->fetch_genome_impl(id);
}

std::optional<std::string> ArchiveReader::fetch_sequence_slice(GenomeId id,
                                                               uint64_t start,
                                                               uint64_t length) const {
    return impl_->fetch_sequence_slice_impl(id, start, length);
}

std::optional<GenomeMeta> ArchiveReader::genome_meta_by_accession(
    std::string_view accession) const
{
    auto it = impl_->accession_map_.find(std::string(accession));
    if (it == impl_->accession_map_.end()) return std::nullopt;
    auto meta_it = impl_->genome_meta_map_.find(it->second);
    if (meta_it == impl_->genome_meta_map_.end()) return std::nullopt;
    return *meta_it->second;
}

std::optional<ExtractedGenome> ArchiveReader::fetch_by_accession(
    std::string_view accession) const
{
    auto it = impl_->accession_map_.find(std::string(accession));
    if (it == impl_->accession_map_.end()) return std::nullopt;
    return fetch_genome(it->second);
}

std::vector<std::optional<ExtractedGenome>>
ArchiveReader::batch_fetch_by_accessions(
    const std::vector<std::string>& accessions) const
{
    const size_t n = accessions.size();
    std::vector<std::optional<ExtractedGenome>> results(n);

    // Group requests by shard_id to read each shard exactly once.
    struct Req { size_t out_idx; GenomeId gid; };
    std::unordered_map<uint32_t, std::vector<Req>> by_shard;
    by_shard.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        auto it = impl_->accession_map_.find(accessions[i]);
        if (it == impl_->accession_map_.end()) continue;
        GenomeId gid = it->second;
        auto meta_it = impl_->genome_meta_map_.find(gid);
        if (meta_it == impl_->genome_meta_map_.end()) continue;
        uint32_t shard_id = meta_it->second->shard_id;
        by_shard[shard_id].push_back({i, gid});
    }

    // Fetch all genomes from each shard, then release that shard's pages.
    for (const auto& [shard_id, reqs] : by_shard) {
        const ShardReader& shard = impl_->get_shard(shard_id);
        for (const auto& req : reqs) {
            ExtractedGenome eg;
            eg.meta = *impl_->genome_meta_map_.at(req.gid);
            eg.fasta = shard.fetch_genome(req.gid);
            auto acc_it = impl_->genome_accession_map_.find(req.gid);
            if (acc_it != impl_->genome_accession_map_.end())
                eg.accession = acc_it->second;
            results[req.out_idx] = std::move(eg);
        }
        // Release this shard's mmap pages now that all genomes are extracted.
        shard.release_pages();
    }

    return results;
}

void ArchiveReader::visit_by_shard(
    const std::vector<std::string>& accessions,
    const std::function<void(size_t idx, ExtractedGenome)>& cb) const
{
    // Same grouping as batch_fetch, but calls cb inline per-shard and immediately
    // releases that shard's pages. Peak memory = one shard's genomes at a time.
    struct Req { size_t out_idx; GenomeId gid; };
    std::unordered_map<uint32_t, std::vector<Req>> by_shard;
    by_shard.reserve(accessions.size());

    for (size_t i = 0; i < accessions.size(); ++i) {
        auto it = impl_->accession_map_.find(accessions[i]);
        if (it == impl_->accession_map_.end()) continue;
        GenomeId gid = it->second;
        auto meta_it = impl_->genome_meta_map_.find(gid);
        if (meta_it == impl_->genome_meta_map_.end()) continue;
        by_shard[meta_it->second->shard_id].push_back({i, gid});
    }

    for (const auto& [shard_id, reqs] : by_shard) {
        const ShardReader& shard = impl_->get_shard(shard_id);
        for (const auto& req : reqs) {
            ExtractedGenome eg;
            eg.meta = *impl_->genome_meta_map_.at(req.gid);
            eg.fasta = shard.fetch_genome(req.gid);
            auto acc_it = impl_->genome_accession_map_.find(req.gid);
            if (acc_it != impl_->genome_accession_map_.end())
                eg.accession = acc_it->second;
            cb(req.out_idx, std::move(eg));
        }
        shard.release_pages();
    }
}

void ArchiveReader::visit_shard_batches(
    const std::vector<std::string>& accessions,
    const std::function<void(ShardBatch&)>& cb) const
{
    struct Req { size_t out_idx; GenomeId gid; };
    std::unordered_map<uint32_t, std::vector<Req>> by_shard;
    by_shard.reserve(accessions.size() / 200 + 1);

    for (size_t i = 0; i < accessions.size(); ++i) {
        auto it = impl_->accession_map_.find(accessions[i]);
        if (it == impl_->accession_map_.end()) continue;
        GenomeId gid = it->second;
        auto meta_it = impl_->genome_meta_map_.find(gid);
        if (meta_it == impl_->genome_meta_map_.end()) continue;
        by_shard[meta_it->second->shard_id].push_back({i, gid});
    }

    // Build shard list sorted by file_offset for sequential NFS access.
    // Sorted order enables NFS client read-ahead and makes MADV_WILLNEED effective.
    std::vector<uint32_t> shard_order;
    shard_order.reserve(by_shard.size());
    for (const auto& [shard_id, _] : by_shard)
        shard_order.push_back(shard_id);
    std::sort(shard_order.begin(), shard_order.end(),
        [&](uint32_t a, uint32_t b) {
            auto da = impl_->shard_descs_.find(a);
            auto db = impl_->shard_descs_.find(b);
            if (da == impl_->shard_descs_.end()) return false;
            if (db == impl_->shard_descs_.end()) return true;
            return da->second->file_offset < db->second->file_offset;
        });
    // Double-buffered pread: explicit large I/O per shard instead of mmap page faults.
    // On NFS, MADV_WILLNEED is often silently ignored and mmap faults arrive one 4KB page
    // at a time (N RPCs per shard). A single pread() for the full compressed shard (~1 MB)
    // issues one large NFS RPC and overlaps I/O with genome embedding in a background thread.
    const int fd = impl_->mmap_.fd();

    auto do_pread = [fd](uint64_t offset, uint64_t size, std::vector<uint8_t>& out) {
        out.resize(size);
        uint8_t* p   = out.data();
        off_t    pos = static_cast<off_t>(offset);
        size_t   rem = size;
        while (rem > 0) {
            ssize_t n = ::pread(fd, p, rem, pos);
            if (n < 0) {
                if (errno == EINTR) continue;
                throw std::runtime_error("genopack pread failed: " + std::string(strerror(errno)));
            }
            if (n == 0)
                throw std::runtime_error("genopack pread: unexpected EOF at offset " +
                                         std::to_string(pos) + " (" + std::to_string(rem) + " bytes remaining)");
            p   += static_cast<size_t>(n);
            pos += static_cast<off_t>(n);
            rem -= static_cast<size_t>(n);
        }
    };

    std::array<std::vector<uint8_t>, 2> bufs;
    std::array<ShardBox, 2> boxes;
    int cur = 0;
    std::future<void> bg;

    // Synchronously load first shard before entering the loop.
    if (!shard_order.empty()) {
        auto d = impl_->shard_descs_.find(shard_order[0]);
        if (d != impl_->shard_descs_.end()) {
            do_pread(d->second->file_offset, d->second->compressed_size, bufs[0]);
            boxes[0].open(bufs[0].data(), 0, bufs[0].size());
        }
    }

    ShardBatch batch;
    for (size_t s = 0; s < shard_order.size(); ++s) {
        // Wait for background pread of this shard (skipped for s=0 which was loaded above).
        if (s > 0 && bg.valid()) {
            bg.get();
            boxes[cur].open(bufs[cur].data(), 0, bufs[cur].size());
        }

        // Start background pread of shard s+1 into the other buffer.
        const int nxt = 1 - cur;
        if (s + 1 < shard_order.size()) {
            auto d_nxt = impl_->shard_descs_.find(shard_order[s + 1]);
            if (d_nxt != impl_->shard_descs_.end()) {
                const uint64_t off = d_nxt->second->file_offset;
                const uint64_t sz  = d_nxt->second->compressed_size;
                bg = std::async(std::launch::async,
                    [&bufs, nxt, off, sz, &do_pread]() { do_pread(off, sz, bufs[nxt]); });
            }
        }

        // Parallel decompress + deliver: each blob is independently decompressible.
        // Previously decompression was serial (one genome at a time), wasting 23/24 threads.
        uint32_t shard_id = shard_order[s];
        const auto& reqs  = by_shard.at(shard_id);
        const ShardReader& shard = boxes[cur].reader();
        const int n_reqs = static_cast<int>(reqs.size());

        batch.clear();
        batch.resize(static_cast<size_t>(n_reqs));

        #pragma omp parallel for schedule(dynamic, 1) num_threads(std::min(n_reqs, 8))
        for (int j = 0; j < n_reqs; ++j) {
            const auto& req = reqs[static_cast<size_t>(j)];
            ExtractedGenome eg;
            eg.meta = *impl_->genome_meta_map_.at(req.gid);
            try {
                eg.fasta = shard.fetch_genome(req.gid);
            } catch (const std::exception&) {
                // genome_id in archive metadata but absent from shard directory —
                // return empty FASTA so the caller can handle it gracefully.
                eg.fasta.clear();
            }
            auto acc_it = impl_->genome_accession_map_.find(req.gid);
            if (acc_it != impl_->genome_accession_map_.end())
                eg.accession = acc_it->second;
            batch[static_cast<size_t>(j)] = {req.out_idx, std::move(eg)};
        }
        cb(batch);
        // No release_pages() — buffer is heap-owned, not mmap

        cur = nxt;
    }
    if (bg.valid()) bg.get();
}

std::optional<std::string> ArchiveReader::fetch_sequence_slice_by_accession(
    std::string_view accession, uint64_t start, uint64_t length) const
{
    auto it = impl_->accession_map_.find(std::string(accession));
    if (it == impl_->accession_map_.end()) return std::nullopt;
    return fetch_sequence_slice(it->second, start, length);
}

uint32_t ArchiveReader::find_contig_genome_id(std::string_view contig_acc) const {
    return impl_->cidx_.find(contig_acc);
}

void ArchiveReader::batch_find_contig_genome_ids(const std::string_view* accs,
                                                  uint32_t*               out_genome_ids,
                                                  size_t                  n,
                                                  size_t                  n_threads) const {
    impl_->cidx_.batch_find(accs, out_genome_ids, n, n_threads);
}

// Lazily populate taxonomy_map_ from TAXN sections on first use.
void ArchiveReader::Impl::ensure_taxonomy_loaded() const {
    if (!taxonomy_map_.empty()) return;
    std::unique_lock<std::mutex> lk(shard_open_mx_);
    if (!taxonomy_map_.empty()) return;
    for (auto* sd : toc_.find_by_type(SEC_TAXN)) {
        TaxonomyIndexReader tir;
        tir.open(mmap_.data(), sd->file_offset, sd->compressed_size);
        tir.scan([&](std::string_view acc, std::string_view tax) {
            taxonomy_map_[std::string(acc)] = std::string(tax);
        });
    }
}

std::optional<std::string> ArchiveReader::taxonomy_for_accession(
    std::string_view accession) const
{
    impl_->ensure_taxonomy_loaded();
    auto it = impl_->taxonomy_map_.find(std::string(accession));
    if (it == impl_->taxonomy_map_.end()) return std::nullopt;
    return it->second;
}

void ArchiveReader::scan_taxonomy(
    const std::function<void(std::string_view, std::string_view)>& cb) const
{
    impl_->ensure_taxonomy_loaded();
    for (const auto& [acc, tax] : impl_->taxonomy_map_)
        cb(acc, tax);
}

std::string ArchiveReader::accession_for_genome_id(GenomeId id) const {
    auto it = impl_->genome_accession_map_.find(id);
    if (it == impl_->genome_accession_map_.end()) return {};
    return it->second;
}

void ArchiveReader::scan_genome_accessions(
    const std::function<void(std::string_view, GenomeId)>& cb) const {
    for (const auto& [acc, gid] : impl_->accession_map_)
        cb(acc, gid);
}

std::optional<TaxonomyTree> ArchiveReader::taxonomy_tree() const {
    return impl_->get_tree();
}

std::vector<ExtractedGenome> ArchiveReader::extract(const ExtractQuery& q) const {
    return impl_->extract_impl(q);
}

const float* ArchiveReader::kmer_profile(GenomeId genome_id) const {
    if (!impl_->has_kmrx_) return nullptr;
    for (const auto& r : impl_->kmrx_readers_) {
        const float* p = r.profile_for(genome_id);
        if (p) return p;
    }
    return nullptr;
}

const float* ArchiveReader::kmer_profile_by_accession(std::string_view accession) const {
    if (!impl_->has_kmrx_) return nullptr;
    auto it = impl_->accession_map_.find(std::string(accession));
    if (it == impl_->accession_map_.end()) return nullptr;
    return kmer_profile(it->second);
}

bool ArchiveReader::has_sketches() const {
    return !impl_->skch_descs_.empty();
}

std::optional<SketchResult> ArchiveReader::sketch_for(GenomeId genome_id) const {
    if (impl_->skch_descs_.empty()) return std::nullopt;

    // Lazy-load SKCH sections on first call (thread-safe)
    if (!impl_->skch_loaded_) {
        std::lock_guard<std::mutex> lk(impl_->shard_open_mx_);
        if (!impl_->skch_loaded_) {
            impl_->skch_readers_.reserve(impl_->skch_descs_.size());
            for (const auto& desc : impl_->skch_descs_) {
                impl_->skch_readers_.emplace_back();
                impl_->skch_readers_.back().open(impl_->mmap_.data(),
                                                 desc.file_offset, desc.compressed_size);
            }
            impl_->skch_loaded_ = true;
        }
    }

    for (const auto& r : impl_->skch_readers_) {
        auto res = r.sketch_for(genome_id);
        if (res) return res;
    }
    return std::nullopt;
}

void ArchiveReader::release_sketches() const {
    // V4 readers don't hold decompressed buffers between calls
    // (sketch_for / sketch_for_ids decompress on demand), so there's
    // nothing to release. Kept for API compatibility.
}

size_t ArchiveReader::sketch_memory_bytes() const {
    // V4 readers hold only the lightweight id/frame index; per-frame
    // buffers are transient. Report the index size across readers.
    size_t total = 0;
    for (const auto& r : impl_->skch_readers_) {
        total += r.genome_ids().size() * sizeof(uint64_t);
    }
    return total;
}

bool ArchiveReader::has_sig2() const {
    return !impl_->skch_descs_.empty();  // V4 always has sig2
}

std::optional<SketchResult> ArchiveReader::sketch_for(GenomeId genome_id,
                                                       uint32_t k,
                                                       uint32_t sz) const {
    if (impl_->skch_descs_.empty()) return std::nullopt;

    // Trigger lazy load (reuses same lock and readers as the unparameterised overload).
    if (!impl_->skch_loaded_) {
        std::lock_guard<std::mutex> lk(impl_->shard_open_mx_);
        if (!impl_->skch_loaded_) {
            impl_->skch_readers_.reserve(impl_->skch_descs_.size());
            for (const auto& desc : impl_->skch_descs_) {
                impl_->skch_readers_.emplace_back();
                impl_->skch_readers_.back().open(impl_->mmap_.data(),
                                                 desc.file_offset, desc.compressed_size);
            }
            impl_->skch_loaded_ = true;
        }
    }

    for (const auto& r : impl_->skch_readers_) {
        if (!r.has_kmer_size(k)) continue;
        if (r.sketch_size() < sz) continue;
        auto res = r.sketch_for(genome_id, k, sz);
        if (res) return res;
    }
    return std::nullopt;
}

void ArchiveReader::sketch_for_ids(const std::vector<GenomeId>& sorted_ids,
                                    uint32_t k, uint32_t sz,
                                    const SketchCallback& cb) const
{
    if (sorted_ids.empty() || impl_->skch_descs_.empty()) return;

    if (!impl_->skch_loaded_) {
        std::lock_guard<std::mutex> lk(impl_->shard_open_mx_);
        if (!impl_->skch_loaded_) {
            impl_->skch_readers_.reserve(impl_->skch_descs_.size());
            for (const auto& desc : impl_->skch_descs_) {
                impl_->skch_readers_.emplace_back();
                impl_->skch_readers_.back().open(impl_->mmap_.data(),
                                                 desc.file_offset, desc.compressed_size);
            }
            impl_->skch_loaded_ = true;
        }
    }

    for (const auto& r : impl_->skch_readers_) {
        if (k > 0 && !r.has_kmer_size(k)) continue;
        if (sz > 0 && r.sketch_size() < sz) continue;
        r.sketch_for_ids(sorted_ids, k, sz, cb);
    }
}

uint32_t ArchiveReader::sketch_kmer_size() const {
    if (impl_->skch_descs_.empty()) return 0;
    const auto& ks = impl_->skch_descs_[0].kmer_sizes;
    return ks.empty() ? 0 : ks[0];
}

std::vector<uint32_t> ArchiveReader::available_sketch_kmer_sizes() const {
    std::vector<uint32_t> all;
    for (const auto& desc : impl_->skch_descs_)
        for (uint32_t k : desc.kmer_sizes)
            if (std::find(all.begin(), all.end(), k) == all.end())
                all.push_back(k);
    std::sort(all.begin(), all.end());
    return all;
}

uint32_t ArchiveReader::sketch_sketch_size() const {
    if (impl_->skch_descs_.empty()) return 0;
    return impl_->skch_descs_[0].sketch_size;
}

void ArchiveReader::scan_shards(
    const std::function<void(const uint8_t* data,
                             uint64_t offset,
                             uint64_t compressed_size,
                             uint32_t shard_id)>& cb) const
{
    // Collect shard sections sorted by file_offset for sequential access
    struct ShardInfo {
        uint64_t offset;
        uint64_t size;
        uint32_t shard_id;
    };
    std::vector<ShardInfo> shards;
    shards.reserve(impl_->shard_descs_.size());
    for (const auto& [sid, sd] : impl_->shard_descs_) {
        shards.push_back({sd->file_offset, sd->compressed_size, sid});
    }
    std::sort(shards.begin(), shards.end(),
              [](const ShardInfo& a, const ShardInfo& b) { return a.offset < b.offset; });

    for (const auto& s : shards) {
        cb(impl_->mmap_.data(), s.offset, s.size, s.shard_id);
    }
}

int ArchiveReader::fd() const {
    return impl_->mmap_.fd();
}

ArchiveReader::ArchiveStats ArchiveReader::archive_stats() const {
    ArchiveStats s{};
    s.generation      = impl_->generation_;
    s.n_shards        = impl_->n_shards_;
    s.n_genomes_total = impl_->total_count_;
    s.n_genomes_live  = impl_->live_count_;

    impl_->catalog_.scan([&](const GenomeMeta& m) {
        if (!m.is_deleted() && !impl_->is_deleted(m.genome_id))
            s.total_raw_bp += m.genome_length;
        return true;
    });
    for (const auto& sd : impl_->toc_.sections)
        if (sd.type == SEC_SHRD)
            s.total_compressed_bytes += sd.compressed_size;

    s.compression_ratio = (s.total_compressed_bytes > 0)
        ? static_cast<double>(s.total_raw_bp) / s.total_compressed_bytes
        : 0.0;
    return s;
}

} // namespace genopack
