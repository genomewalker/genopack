#include <genopack/archive.hpp>
#include <genopack/cidx.hpp>
#include <genopack/format.hpp>
#include <genopack/meta.hpp>
#include <genopack/gidx.hpp>
#include <genopack/gcov.hpp>
#include <genopack/gstx.hpp>
#include <genopack/qual.hpp>
#include <genopack/kmrx.hpp>
#include <genopack/skch.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/toc.hpp>
#include <genopack/accx.hpp>
#include <genopack/taxn.hpp>
#include <genopack/txdb.hpp>
#include <genopack/tombstone.hpp>
#include <spdlog/spdlog.h>
#include <unistd.h>
#include <sys/mman.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <cstring>
#include <future>
#include <mutex>
#include <stdexcept>
#ifdef __unix__
#include <sys/mman.h>
#endif
#include <string>
#include <fcntl.h>
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

    // ACCX readers kept alive so accession lookups are zero-copy into the mmap.
    // genome_accession_map_ stores const char* into the ACCX string area (valid
    // while mmap_ is open) — avoids duplicating 2.3M accession strings on heap.
    std::vector<AccessionIndexReader>         accx_readers_;
    std::unordered_map<GenomeId, const char*> genome_accession_map_;
    std::unordered_map<GenomeId, const GenomeMeta*> genome_meta_map_;

    // TAXN readers kept alive — open() just sets pointers into the mmap, O(1),
    // zero heap allocation. Replaces the old lazy taxonomy_map_ copy (~900 MB).
    std::vector<TaxonomyIndexReader> taxn_readers_;

    // TXDB section descriptors (section_id -> SectionDesc pointer into toc_.sections)
    std::unordered_map<uint64_t, const SectionDesc*> txdb_descs_;

    // KMRX readers — all sections (merged archives have one per original part)
    std::vector<KmrxReader> kmrx_readers_;
    bool                    has_kmrx_ = false;

    // GIDX reader for O(1) genome_id -> shard position lookup
    GidxReader gidx_;
    bool       has_gidx_ = false;

    // GSTX reader: genus sketch stats (consensus + p90 + TNF centroid)
    std::vector<uint8_t> gstx_buf_;  // pread heap buffer — avoids NFS mmap page faults
    GstxReader gstx_;
    bool       has_gstx_ = false;

    // GCOV reader: per-genus biological covariance eigenbasis + calibrated quantiles
    std::vector<uint8_t> gcov_buf_;
    GcovReader gcov_;
    bool       has_gcov_ = false;

    // FCOV reader: per-family biological covariance (same layout as GCOV)
    std::vector<uint8_t> fcov_buf_;
    GcovReader fcov_;
    bool       has_fcov_ = false;

    // FMHR reader: per-genus FracMinHash reference sketches (k=21, c=125)
    FmhrReader fmhr_;
    bool       has_fmhr_ = false;

    // QUAL reader: per-genome quality scores
    std::vector<uint8_t> qual_buf_;  // pread heap buffer — avoids NFS mmap page faults
    QualReader qual_;
    bool       has_qual_ = false;

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

    // O(1) accession → GenomeId lookup across all ACCX readers (mmap'd hash tables).
    std::optional<GenomeId> find_accession(std::string_view acc) const {
        for (const auto& r : accx_readers_) {
            auto gid = r.find(acc);
            if (gid) return gid;
        }
        return std::nullopt;
    }

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

    // pread a section into a fresh heap buffer — one NFS bulk transfer,
    // no mmap page-fault storm. Returns the buffer; throws on short read.
    static std::vector<uint8_t> pread_section(int fd, uint64_t offset, uint64_t size) {
        std::vector<uint8_t> buf(size);
        uint8_t* p = buf.data();
        uint64_t remaining = size;
        while (remaining > 0) {
            ssize_t n = ::pread(fd, p, remaining, static_cast<off_t>(offset));
            if (n <= 0) throw std::runtime_error("pread_section: short read at offset "
                                                  + std::to_string(offset));
            p         += n;
            offset    += static_cast<uint64_t>(n);
            remaining -= static_cast<uint64_t>(n);
        }
        return buf;
    }

    void open_gpk(const std::filesystem::path& path) {
        mmap_.open(path);

        // Read FileHeader from mmap offset 0 — always a warm page, no fault.
        auto* fh = mmap_.ptr_at<FileHeader>(0);
        if (fh->magic != GPK2_MAGIC)
            throw std::runtime_error("Not a .gpk file");
        // Reject archives whose major format version is newer than we support,
        // rather than silently misparsing a future layout as garbage (P15).
        if (fh->version_major > FORMAT_MAJOR)
            throw std::runtime_error(
                "genopack: archive format major v" + std::to_string(fh->version_major) +
                " > supported v" + std::to_string(FORMAT_MAJOR) +
                " — upgrade genopack to read this archive");

        // pread TailLocator from EOF — avoids mmap page-fault at multi-TB offset.
        TailLocator tail{};
        {
            const uint64_t tail_off = mmap_.size() - sizeof(TailLocator);
            const ssize_t n = ::pread(mmap_.fd(), &tail, sizeof(tail),
                                      static_cast<off_t>(tail_off));
            if (n != static_cast<ssize_t>(sizeof(tail)) || tail.magic != GPKT_MAGIC)
                throw std::runtime_error("Not a .gpk file or TailLocator corrupt");
        }

        // If MetaBundle is present: one more pread gives the full section directory.
        // This is the fast path for any archive written by genopack >= this version.
        const uint64_t mb_off  = tail.meta_bundle_offset();
        const uint64_t mb_size = tail.meta_bundle_size();
        MetaBundleReader meta;
        if (mb_off != 0 && mb_size != 0 && meta.open(mmap_.fd(), mb_off, mb_size)) {
            toc_ = meta.to_toc();
            spdlog::debug("genopack: MetaBundle fast-open ({} sections)", toc_.sections.size());
        } else {
            // Legacy archive (MetaBundle absent): fall back to pread-based TOC read.
            // Slower on NFS (page-faults at multi-TB offset) but correct.
            spdlog::info("genopack: legacy archive (no MetaBundle), reading TOC via pread — "
                         "run 'genopack reindex' to add MetaBundle for fast open");
            const uint64_t toc_off  = tail.toc_offset;
            const uint64_t toc_size = tail.toc_size;
            std::vector<uint8_t> toc_buf(toc_size);
            const ssize_t nr = ::pread(mmap_.fd(), toc_buf.data(), toc_size,
                                       static_cast<off_t>(toc_off));
            if (nr != static_cast<ssize_t>(toc_size))
                throw std::runtime_error("genopack: legacy archive TOC pread failed");
            toc_ = TocReader::read_at(toc_buf.data(), 0, toc_size);
        }

        // Validate every section descriptor against the file size so a truncated
        // or corrupt directory fails loudly here instead of faulting (SIGBUS/OOB)
        // when a section is lazily touched (P5).
        {
            const uint64_t fsz = mmap_.size();
            for (const auto& sd : toc_.sections) {
                if (sd.compressed_size > fsz || sd.file_offset > fsz - sd.compressed_size)
                    throw std::runtime_error(
                        "genopack: section descriptor out of bounds (type=" +
                        std::to_string(sd.type) + " offset=" + std::to_string(sd.file_offset) +
                        " size=" + std::to_string(sd.compressed_size) +
                        " file=" + std::to_string(fsz) + ") — archive truncated/corrupt");
            }
        }

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

        // Load GSTX section (highest section_id wins)
        {
            uint64_t best_id = 0;
            const SectionDesc* best_sd = nullptr;
            for (auto* sd : toc_.find_by_type(SEC_GSTX)) {
                if (sd->section_id > best_id) {
                    best_id = sd->section_id;
                    best_sd = sd;
                }
            }
            if (best_sd) {
                gstx_buf_ = pread_section(mmap_.fd(), best_sd->file_offset,
                                          best_sd->compressed_size);
                gstx_.open(gstx_buf_.data(), 0, best_sd->compressed_size);
                has_gstx_ = true;
                // Release mmap pages for this region — data is now in gstx_buf_.
                mmap_.advise(best_sd->file_offset, best_sd->compressed_size, MADV_DONTNEED);
            }
        }

        // Load GCOV section (highest section_id wins)
        {
            uint64_t best_id = 0;
            const SectionDesc* best_sd = nullptr;
            for (auto* sd : toc_.find_by_type(SEC_GCOV)) {
                if (sd->section_id > best_id) {
                    best_id = sd->section_id;
                    best_sd = sd;
                }
            }
            if (best_sd) {
                try {
                    gcov_buf_ = pread_section(mmap_.fd(), best_sd->file_offset,
                                              best_sd->compressed_size);
                    gcov_.open(gcov_buf_.data(), 0, best_sd->compressed_size);
                    has_gcov_ = true;
                    mmap_.advise(best_sd->file_offset, best_sd->compressed_size, MADV_DONTNEED);
                } catch (const std::exception& ex) {
                    spdlog::warn("GCOV section unreadable ({}); rebuild with 'genopack gcov'", ex.what());
                }
            }
        }

        // Load FCOV section (per-family covariance; highest section_id wins)
        {
            uint64_t best_id = 0;
            const SectionDesc* best_sd = nullptr;
            for (auto* sd : toc_.find_by_type(SEC_FCOV)) {
                if (sd->section_id > best_id) { best_id = sd->section_id; best_sd = sd; }
            }
            if (best_sd) {
                try {
                    fcov_buf_ = pread_section(mmap_.fd(), best_sd->file_offset,
                                              best_sd->compressed_size);
                    fcov_.open(fcov_buf_.data(), 0, best_sd->compressed_size);
                    has_fcov_ = true;
                    mmap_.advise(best_sd->file_offset, best_sd->compressed_size, MADV_DONTNEED);
                } catch (const std::exception& ex) {
                    spdlog::warn("FCOV section unreadable ({}); rebuild with 'genopack gcov'", ex.what());
                }
            }
        }

        // Load FMHR section (highest section_id wins; zero-copy — points into mmap)
        {
            uint64_t best_id = 0;
            const SectionDesc* best_sd = nullptr;
            for (auto* sd : toc_.find_by_type(SEC_FMHR)) {
                if (sd->section_id > best_id) { best_id = sd->section_id; best_sd = sd; }
            }
            if (best_sd) {
                try {
                    fmhr_.open(mmap_.data(), best_sd->file_offset, best_sd->compressed_size);
                    has_fmhr_ = true;
                } catch (const std::exception& ex) {
                    spdlog::warn("FMHR section unreadable ({}); rebuild with 'genopack gcov'", ex.what());
                }
            }
        }

        // Load QUAL section (highest section_id wins)
        {
            uint64_t best_id = 0;
            const SectionDesc* best_sd = nullptr;
            for (auto* sd : toc_.find_by_type(SEC_QUAL)) {
                if (sd->section_id > best_id) {
                    best_id = sd->section_id;
                    best_sd = sd;
                }
            }
            if (best_sd) {
                try {
                    qual_buf_ = pread_section(mmap_.fd(), best_sd->file_offset,
                                              best_sd->compressed_size);
                    qual_.open(qual_buf_.data(), 0, best_sd->compressed_size);
                    has_qual_ = true;
                    mmap_.advise(best_sd->file_offset, best_sd->compressed_size, MADV_DONTNEED);
                } catch (const std::exception& e) {
                    spdlog::warn("QUAL section unreadable ({}); rebuild with reindex", e.what());
                }
            }
        }

        // Load CIDX sections (contig accession → genome_id)
        for (auto* sd : toc_.find_by_type(SEC_CIDX))
            cidx_.add_section(mmap_.data(), sd->file_offset, sd->compressed_size);

        // Load accession index from all ACCX sections.
        // Keep readers alive so lookups use the mmap directly (no string copies).
        // genome_accession_map_ stores const char* into the mmap string area.
        for (auto* sd : toc_.find_by_type(SEC_ACCX)) {
            accx_readers_.emplace_back();
            accx_readers_.back().open(mmap_.data(), sd->file_offset, sd->compressed_size);
            accx_readers_.back().scan([&](std::string_view acc, GenomeId gid) {
                genome_accession_map_[gid] = acc.data(); // points into mmap string area
            });
        }

        // Open TAXN readers eagerly — each open() sets 5 pointers into the existing mmap,
        // O(1) work, zero heap allocation. Replaces the old lazy ~900 MB taxonomy_map_ copy.

        for (auto* sd : toc_.find_by_type(SEC_TAXN)) {
            taxn_readers_.emplace_back();
            taxn_readers_.back().open(mmap_.data(), sd->file_offset, sd->compressed_size);
        }

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
            accx_readers_.clear();
            genome_accession_map_.clear();
            genome_meta_map_.clear();
            taxn_readers_.clear();
            txdb_descs_.clear();
            gidx_                 = GidxReader{};
            has_gidx_             = false;
            gstx_                 = GstxReader{};
            has_gstx_             = false;
            gcov_                 = GcovReader{};
            has_gcov_             = false;
            gcov_buf_.clear();
            fcov_                 = GcovReader{};
            has_fcov_             = false;
            fcov_buf_.clear();
            fmhr_                 = FmhrReader{};
            has_fmhr_             = false;
            qual_                 = QualReader{};
            has_qual_             = false;
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
            eg.accession = std::string(acc_it->second);

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
                        eg.accession = std::string(acc_it->second);
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
                        eg.accession = std::string(acc_it->second);
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
    impl_->accx_readers_.clear();
    impl_->genome_accession_map_.clear();
    impl_->genome_meta_map_.clear();
    impl_->taxn_readers_.clear();
    impl_->txdb_descs_.clear();
    impl_->gidx_                 = GidxReader{};
    impl_->has_gidx_             = false;
    impl_->gstx_                 = GstxReader{};
    impl_->has_gstx_             = false;
    impl_->gcov_                 = GcovReader{};
    impl_->has_gcov_             = false;
    impl_->gcov_buf_.clear();
    impl_->fcov_                 = GcovReader{};
    impl_->has_fcov_             = false;
    impl_->fcov_buf_.clear();
    impl_->qual_                 = QualReader{};
    impl_->has_qual_             = false;
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
    auto gid = impl_->find_accession(accession);
    if (!gid) return std::nullopt;
    auto meta_it = impl_->genome_meta_map_.find(*gid);
    if (meta_it == impl_->genome_meta_map_.end()) return std::nullopt;
    return *meta_it->second;
}

std::optional<ExtractedGenome> ArchiveReader::fetch_by_accession(
    std::string_view accession) const
{
    auto gid = impl_->find_accession(accession);
    if (!gid) return std::nullopt;
    return fetch_genome(*gid);
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
        auto gid = impl_->find_accession(accessions[i]);
        if (!gid) continue;
        auto meta_it = impl_->genome_meta_map_.find(*gid);
        if (meta_it == impl_->genome_meta_map_.end()) continue;
        by_shard[meta_it->second->shard_id].push_back({i, *gid});
    }

    // Sort shards by file_offset for sequential NFS reads — avoids random seeks
    // across the pack file and enables OS read-ahead on contiguous access.
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

    // Fetch all genomes from each shard in offset order, then release that shard's pages.
    for (uint32_t shard_id : shard_order) {
        const auto& reqs = by_shard.at(shard_id);
        const ShardReader& shard = impl_->get_shard(shard_id);
        for (const auto& req : reqs) {
            ExtractedGenome eg;
            eg.meta = *impl_->genome_meta_map_.at(req.gid);
            eg.fasta = shard.fetch_genome(req.gid);
            auto acc_it = impl_->genome_accession_map_.find(req.gid);
            if (acc_it != impl_->genome_accession_map_.end())
                eg.accession = std::string(acc_it->second);
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
        auto gid = impl_->find_accession(accessions[i]);
        if (!gid) continue;
        auto meta_it = impl_->genome_meta_map_.find(*gid);
        if (meta_it == impl_->genome_meta_map_.end()) continue;
        by_shard[meta_it->second->shard_id].push_back({i, *gid});
    }

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

    for (uint32_t shard_id : shard_order) {
        const auto& reqs = by_shard.at(shard_id);
        const ShardReader& shard = impl_->get_shard(shard_id);
        for (const auto& req : reqs) {
            ExtractedGenome eg;
            eg.meta = *impl_->genome_meta_map_.at(req.gid);
            eg.fasta = shard.fetch_genome(req.gid);
            auto acc_it = impl_->genome_accession_map_.find(req.gid);
            if (acc_it != impl_->genome_accession_map_.end())
                eg.accession = std::string(acc_it->second);
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
        auto gid = impl_->find_accession(accessions[i]);
        if (!gid) continue;
        auto meta_it = impl_->genome_meta_map_.find(*gid);
        if (meta_it == impl_->genome_meta_map_.end()) continue;
        by_shard[meta_it->second->shard_id].push_back({i, *gid});
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
                eg.accession = std::string(acc_it->second);
            batch[static_cast<size_t>(j)] = {req.out_idx, std::move(eg)};
        }
        cb(batch);
        // No release_pages() — buffer is heap-owned, not mmap

        cur = nxt;
    }
    if (bg.valid()) bg.get();
}

void ArchiveReader::visit_shard_batches_parallel(
    const std::vector<std::string>& accessions,
    int n_readers,
    const std::function<void(ShardBatch&)>& cb) const
{
    if (n_readers <= 1) { visit_shard_batches(accessions, cb); return; }

    struct Req { size_t out_idx; GenomeId gid; };
    std::unordered_map<uint32_t, std::vector<Req>> by_shard;
    by_shard.reserve(accessions.size() / 200 + 1);

    for (size_t i = 0; i < accessions.size(); ++i) {
        auto gid = impl_->find_accession(accessions[i]);
        if (!gid) continue;
        auto meta_it = impl_->genome_meta_map_.find(*gid);
        if (meta_it == impl_->genome_meta_map_.end()) continue;
        by_shard[meta_it->second->shard_id].push_back({i, *gid});
    }

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

    const int src_fd = impl_->mmap_.fd();
    const int n = static_cast<int>(shard_order.size());
    n_readers = std::min(n_readers, std::max(1, n));

    std::exception_ptr first_error;
    std::mutex err_mtx;
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(n_readers));

    for (int t = 0; t < n_readers; ++t) {
        const int band_start = (t * n) / n_readers;
        const int band_end   = ((t + 1) * n) / n_readers;
        if (band_start >= band_end) continue;

        const int tfd = ::dup(src_fd);
        if (tfd < 0) throw std::runtime_error("visit_shard_batches_parallel: dup failed");
        ::posix_fadvise(tfd, 0, 0, POSIX_FADV_SEQUENTIAL);

        workers.emplace_back([&, band_start, band_end, tfd]() {
            try {
                auto do_pread = [tfd](uint64_t offset, uint64_t size, std::vector<uint8_t>& out) {
                    out.resize(size);
                    uint8_t* p   = out.data();
                    off_t    pos = static_cast<off_t>(offset);
                    size_t   rem = size;
                    while (rem > 0) {
                        ssize_t r = ::pread(tfd, p, rem, pos);
                        if (r < 0) {
                            if (errno == EINTR) continue;
                            throw std::runtime_error("genopack pread failed: " + std::string(strerror(errno)));
                        }
                        if (r == 0)
                            throw std::runtime_error("genopack pread: unexpected EOF");
                        p   += static_cast<size_t>(r);
                        pos += static_cast<off_t>(r);
                        rem -= static_cast<size_t>(r);
                    }
                };

                std::array<std::vector<uint8_t>, 2> bufs;
                std::array<ShardBox, 2> boxes;
                int cur = 0;
                std::future<void> bg;

                {
                    auto d = impl_->shard_descs_.find(shard_order[static_cast<size_t>(band_start)]);
                    if (d != impl_->shard_descs_.end()) {
                        do_pread(d->second->file_offset, d->second->compressed_size, bufs[0]);
                        boxes[0].open(bufs[0].data(), 0, bufs[0].size());
                    }
                }

                ShardBatch batch;
                for (int s = band_start; s < band_end; ++s) {
                    if (s > band_start && bg.valid()) {
                        bg.get();
                        boxes[cur].open(bufs[cur].data(), 0, bufs[cur].size());
                    }

                    const int nxt = 1 - cur;
                    if (s + 1 < band_end) {
                        auto d_nxt = impl_->shard_descs_.find(shard_order[static_cast<size_t>(s + 1)]);
                        if (d_nxt != impl_->shard_descs_.end()) {
                            const uint64_t off = d_nxt->second->file_offset;
                            const uint64_t sz  = d_nxt->second->compressed_size;
                            bg = std::async(std::launch::async,
                                [&bufs, nxt, off, sz, &do_pread]() { do_pread(off, sz, bufs[nxt]); });
                        }
                    }

                    uint32_t shard_id     = shard_order[static_cast<size_t>(s)];
                    const auto& reqs      = by_shard.at(shard_id);
                    const ShardReader& shard = boxes[cur].reader();
                    const int n_reqs      = static_cast<int>(reqs.size());

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
                            eg.fasta.clear();
                        }
                        auto acc_it = impl_->genome_accession_map_.find(req.gid);
                        if (acc_it != impl_->genome_accession_map_.end())
                            eg.accession = std::string(acc_it->second);
                        batch[static_cast<size_t>(j)] = {req.out_idx, std::move(eg)};
                    }
                    cb(batch);

                    cur = nxt;
                }
                if (bg.valid()) bg.get();
            } catch (...) {
                std::lock_guard<std::mutex> lk(err_mtx);
                if (!first_error) first_error = std::current_exception();
            }
            ::close(tfd);
        });
    }

    for (auto& w : workers) w.join();
    if (first_error) std::rethrow_exception(first_error);
}

std::optional<std::string> ArchiveReader::fetch_sequence_slice_by_accession(
    std::string_view accession, uint64_t start, uint64_t length) const
{
    auto gid = impl_->find_accession(accession);
    if (!gid) return std::nullopt;
    return fetch_sequence_slice(*gid, start, length);
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

std::optional<std::string> ArchiveReader::taxonomy_for_accession(
    std::string_view accession) const
{
    for (const auto& tir : impl_->taxn_readers_) {
        auto t = tir.find(accession);
        if (t) return std::string(*t);
    }
    return std::nullopt;
}

void ArchiveReader::scan_taxonomy(
    const std::function<void(std::string_view, std::string_view)>& cb) const
{
    for (const auto& tir : impl_->taxn_readers_)
        tir.scan(cb);
}

std::string ArchiveReader::accession_for_genome_id(GenomeId id) const {
    auto it = impl_->genome_accession_map_.find(id);
    if (it == impl_->genome_accession_map_.end()) return {};
    return std::string(it->second);
}

void ArchiveReader::scan_genome_accessions(
    const std::function<void(std::string_view, GenomeId)>& cb) const {
    for (const auto& r : impl_->accx_readers_)
        r.scan(cb);
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
    auto gid = impl_->find_accession(accession);
    if (!gid) return nullptr;
    return kmer_profile(*gid);
}

bool ArchiveReader::has_sketches() const {
    return !impl_->skch_descs_.empty();
}

bool ArchiveReader::has_gstx() const { return impl_->has_gstx_; }

const GstxEntry* ArchiveReader::gstx_for_genus(std::string_view genus) const {
    if (!impl_->has_gstx_) return nullptr;
    return impl_->gstx_.lookup(genus);
}

std::vector<uint32_t> ArchiveReader::gstx_kmer_sizes() const {
    if (!impl_->has_gstx_) return {};
    std::vector<uint32_t> ks;
    for (uint32_t ki = 0; ki < impl_->gstx_.n_k(); ++ki)
        ks.push_back(impl_->gstx_.kmer_size(ki));
    return ks;
}

const GstxReader* ArchiveReader::gstx_reader() const {
    return impl_->has_gstx_ ? &impl_->gstx_ : nullptr;
}

bool ArchiveReader::has_gcov() const { return impl_->has_gcov_; }

const GcovEntry* ArchiveReader::gcov_for_genus(std::string_view genus) const {
    if (!impl_->has_gcov_) return nullptr;
    return impl_->gcov_.lookup(GcovWriter::hash_genus(genus));
}

const GcovReader* ArchiveReader::gcov_reader() const {
    return impl_->has_gcov_ ? &impl_->gcov_ : nullptr;
}

bool ArchiveReader::has_fcov() const { return impl_->has_fcov_; }

const GcovEntry* ArchiveReader::fcov_for_family(std::string_view family) const {
    if (!impl_->has_fcov_) return nullptr;
    return impl_->fcov_.lookup(GcovWriter::hash_genus(family));
}

const GcovReader* ArchiveReader::fcov_reader() const {
    return impl_->has_fcov_ ? &impl_->fcov_ : nullptr;
}

bool ArchiveReader::has_fmhr() const { return impl_->has_fmhr_; }

FmhrView ArchiveReader::fmhr_for_genus(std::string_view genus) const {
    if (!impl_->has_fmhr_) return {};
    return impl_->fmhr_.lookup(GcovWriter::hash_genus(genus));
}

const FmhrReader* ArchiveReader::fmhr_reader() const {
    return impl_->has_fmhr_ ? &impl_->fmhr_ : nullptr;
}

bool ArchiveReader::has_qual() const { return impl_->has_qual_; }

void ArchiveReader::scan_qual(const std::function<void(const QualRecord&)>& cb) const {
    if (impl_->has_qual_) impl_->qual_.scan(cb);
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
    // Advise the kernel it can evict SKCH section pages from the page cache.
    // Without this, mmap pages touched during preload accumulate in RSS across
    // waves — hundreds of GB on nodes with ample RAM — because the kernel has
    // no signal to reclaim them. The preloaded KStore buffer (heap) already
    // holds all needed sketch data, so these file-backed pages are redundant.
    for (const auto& desc : impl_->skch_descs_) {
#ifdef MADV_DONTNEED
        impl_->mmap_.advise(desc.file_offset, desc.compressed_size, MADV_DONTNEED);
#endif
    }
    // Also drop the in-memory id/frame index (small heap; rebuilt on next load).
    impl_->skch_readers_.clear();
    impl_->skch_loaded_ = false;
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
                                    const SketchCallback& cb,
                                    int num_threads) const
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

    // Hint sequential access so NFS/kernel prefetches frames in order.
    // Only effective when num_threads==1 (frames read sequentially).
    if (num_threads == 1) {
        for (const auto& desc : impl_->skch_descs_)
            impl_->mmap_.advise(desc.file_offset, desc.compressed_size, MADV_SEQUENTIAL);
    }

    for (const auto& r : impl_->skch_readers_) {
        if (k > 0 && !r.has_kmer_size(k)) continue;
        if (sz > 0 && r.sketch_size() < sz) continue;
        r.sketch_for_ids(sorted_ids, k, sz, cb, num_threads);
    }
}

void ArchiveReader::sketch_for_ids_multi_k(const std::vector<GenomeId>& sorted_ids,
                                            uint32_t sz,
                                            const SketchCallbackMultiK& cb,
                                            int num_threads) const
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

    if (num_threads == 1) {
        for (const auto& desc : impl_->skch_descs_)
            impl_->mmap_.advise(desc.file_offset, desc.compressed_size, MADV_SEQUENTIAL);
    }

    for (const auto& r : impl_->skch_readers_) {
        if (sz > 0 && r.sketch_size() < sz) continue;
        r.sketch_for_ids_multi_k(sorted_ids, sz, cb, num_threads);
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
