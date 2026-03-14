#include <genopack/archive.hpp>
#include <genopack/format_v2.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/toc.hpp>
#include <genopack/accx.hpp>
#include <genopack/tombstone.hpp>
#include <spdlog/spdlog.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace genopack {

// ── ShardBox ─────────────────────────────────────────────────────────────────
// ShardReaderV2() = default is inline in shard.hpp, which forces the compiler
// to instantiate the destructor of unique_ptr<ShardReaderV2::Impl> here — but
// ShardReaderV2::Impl is incomplete in this TU. Workaround: zero-initialize
// aligned storage for a ShardReaderV2, then call only the out-of-line methods
// open() and ~ShardReaderV2(). A zero-initialized unique_ptr<Impl> is nullptr
// on all platforms, so open() (which does `if (!impl_) impl_ = make_unique<Impl>()`)
// initialises it correctly in shard.cpp where Impl is complete.

struct ShardBox {
    alignas(ShardReaderV2) unsigned char storage_[sizeof(ShardReaderV2)];
    bool live_ = false;

    ShardBox() { std::memset(storage_, 0, sizeof(storage_)); }

    ~ShardBox() {
        if (live_) reinterpret_cast<ShardReaderV2*>(storage_)->~ShardReaderV2();
    }

    ShardBox(const ShardBox&)            = delete;
    ShardBox& operator=(const ShardBox&) = delete;

    ShardReaderV2& reader() {
        return *reinterpret_cast<ShardReaderV2*>(storage_);
    }
    const ShardReaderV2& reader() const {
        return *reinterpret_cast<const ShardReaderV2*>(storage_);
    }

    void open(const uint8_t* base, uint64_t offset, uint64_t size) {
        reader().open(base, offset, size);
        live_ = true;
    }
};

// ── ArchiveReader::Impl ───────────────────────────────────────────────────────

struct ArchiveReader::Impl {
    // ── v2 state ──────────────────────────────────────────────────────────────
    MmapFileReader      mmap_;
    Toc                 toc_;
    MergedCatalogReader catalog_v2_;

    // shard_id -> SectionDesc pointer (into toc_.sections)
    std::unordered_map<uint32_t, const SectionDesc*> shard_descs_;
    // lazy-opened v2 shard readers (ShardBox avoids the incomplete-Impl issue)
    mutable std::unordered_map<uint32_t, ShardBox> shards_v2_;

    // accession <-> genome_id maps built from all ACCX sections
    std::unordered_map<std::string, GenomeId> accession_map_;
    std::unordered_map<GenomeId, std::string> genome_accession_map_;

    // merged tombstones from all TOMB sections
    std::vector<TombstoneReader> tombstones_;

    // ── v1 state ──────────────────────────────────────────────────────────────
    std::filesystem::path archive_dir_v1_;
    ManifestHeader        manifest_v1_{};
    CatalogReader         catalog_v1_;
    mutable std::unordered_map<ShardId, ShardReader> shards_v1_;

    // ── common state ──────────────────────────────────────────────────────────
    bool     is_v2_       = false;
    bool     open_        = false;
    uint64_t live_count_  = 0;
    uint64_t total_count_ = 0;
    uint64_t generation_  = 0;
    uint32_t n_shards_    = 0;

    // ── v1 helpers ────────────────────────────────────────────────────────────

    void load_meta_tsv(const std::filesystem::path& path) {
        std::ifstream f(path);
        if (!f) {
            spdlog::debug("genopack: meta.tsv not found at {}, accession lookup unavailable",
                          path.string());
            return;
        }

        std::string line;
        if (!std::getline(f, line)) return;  // skip header

        while (std::getline(f, line)) {
            if (line.empty()) continue;
            std::istringstream ss(line);
            std::string accession, genome_id_str;
            if (!std::getline(ss, accession, '\t')) continue;
            if (!std::getline(ss, genome_id_str, '\t')) continue;

            GenomeId gid = std::stoull(genome_id_str);
            accession_map_[accession]  = gid;
            genome_accession_map_[gid] = accession;
        }
        spdlog::debug("genopack: loaded {} accession mappings from meta.tsv",
                      accession_map_.size());
    }

    void open_v1(const std::filesystem::path& dir) {
        archive_dir_v1_ = dir;

        std::ifstream mf(dir / "MANIFEST.bin", std::ios::binary);
        if (!mf) throw std::runtime_error("Cannot open MANIFEST.bin in " + dir.string());
        mf.read(reinterpret_cast<char*>(&manifest_v1_), sizeof(manifest_v1_));
        if (manifest_v1_.magic != GPKM_MAGIC)
            throw std::runtime_error("Invalid MANIFEST.bin magic");

        catalog_v1_.open(dir / "catalog.gpkc");
        load_meta_tsv(dir / "meta.tsv");

        live_count_  = manifest_v1_.n_genomes_live;
        total_count_ = manifest_v1_.n_genomes;
        generation_  = manifest_v1_.generation;
        n_shards_    = manifest_v1_.n_shards;

        is_v2_ = false;
        open_  = true;
        spdlog::info("genopack v1 opened: {} genomes, {} shards",
                     live_count_, n_shards_);
    }

    void open_v2(const std::filesystem::path& path) {
        mmap_.open(path);

        auto* fh = mmap_.ptr_at<FileHeader>(0);
        if (fh->magic != GPK2_MAGIC)
            throw std::runtime_error("Not a v2 .gpk file");

        toc_ = TocReader::read(mmap_);

        // Load catalog fragments (all CATL sections)
        for (auto* sd : toc_.find_by_type(SEC_CATL)) {
            catalog_v2_.add_fragment(mmap_.data(), sd->file_offset, sd->compressed_size);
        }

        // Index shard sections by shard_id (aux0)
        n_shards_ = 0;
        for (auto& sd : toc_.sections) {
            if (sd.type == SEC_SHRD) {
                shard_descs_[static_cast<uint32_t>(sd.aux0)] = &sd;
                ++n_shards_;
            }
        }

        // Load accession index from all ACCX sections
        for (auto* sd : toc_.find_by_type(SEC_ACCX)) {
            AccessionIndexReader reader;
            reader.open(mmap_.data(), sd->file_offset, sd->compressed_size);
            reader.scan([&](std::string_view acc, GenomeId gid) {
                accession_map_[std::string(acc)] = gid;
                genome_accession_map_[gid]       = std::string(acc);
            });
        }

        // Load tombstones from all TOMB sections
        for (auto* sd : toc_.find_by_type(SEC_TOMB)) {
            tombstones_.emplace_back();
            tombstones_.back().open(mmap_.data(), sd->file_offset, sd->compressed_size);
        }

        live_count_  = toc_.header.live_genome_count;
        total_count_ = toc_.header.total_genome_count;
        generation_  = toc_.header.generation;

        is_v2_ = true;
        open_  = true;
        spdlog::info("genopack v2 opened: {} genomes, {} shards, gen {}",
                     live_count_, n_shards_, generation_);
    }

    void open(const std::filesystem::path& path) {
        if (std::filesystem::is_directory(path)) {
            open_v1(path);
            return;
        }

        std::ifstream probe(path, std::ios::binary);
        if (!probe) throw std::runtime_error("Cannot open: " + path.string());
        uint32_t magic = 0;
        probe.read(reinterpret_cast<char*>(&magic), 4);
        probe.close();

        if (magic == GPK2_MAGIC) {
            open_v2(path);
        } else if (magic == GPKM_MAGIC) {
            throw std::runtime_error(
                "v1 file format not supported for single-file open; use directory path");
        } else {
            auto gpk = path;
            if (gpk.extension() != ".gpk") gpk += ".gpk";
            if (std::filesystem::exists(gpk)) {
                open_v2(gpk);
            } else {
                throw std::runtime_error("Unknown file format at: " + path.string());
            }
        }
    }

    // ── v1 shard access ───────────────────────────────────────────────────────

    ShardReader& get_shard_v1(ShardId shard_id) const {
        auto it = shards_v1_.find(shard_id);
        if (it != shards_v1_.end()) return it->second;

        for (const auto& entry :
             std::filesystem::directory_iterator(archive_dir_v1_ / "shards")) {
            if (entry.path().extension() != ".gpks") continue;
            ShardReader tmp;
            tmp.open(entry.path());
            if (tmp.shard_id() == shard_id) {
                auto [ins_it, ok] = shards_v1_.emplace(shard_id, std::move(tmp));
                return ins_it->second;
            }
        }
        throw std::runtime_error("Shard " + std::to_string(shard_id) + " not found");
    }

    // ── v2 shard access ───────────────────────────────────────────────────────

    ShardReaderV2& get_shard_v2(uint32_t shard_id) const {
        auto it = shards_v2_.find(shard_id);
        if (it != shards_v2_.end()) return it->second.reader();

        auto desc_it = shard_descs_.find(shard_id);
        if (desc_it == shard_descs_.end())
            throw std::runtime_error("Shard " + std::to_string(shard_id) + " not found");

        const SectionDesc* sd = desc_it->second;
        auto [ins, ok] = shards_v2_.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(shard_id),
            std::forward_as_tuple());
        ins->second.open(mmap_.data(), sd->file_offset, sd->compressed_size);
        return ins->second.reader();
    }

    // ── tombstone check ───────────────────────────────────────────────────────

    bool is_deleted(GenomeId id) const {
        for (const auto& t : tombstones_)
            if (t.is_deleted(id)) return true;
        return false;
    }

    // ── v2 fetch helpers ──────────────────────────────────────────────────────

    std::optional<ExtractedGenome> fetch_genome_v2(GenomeId id) const {
        const GenomeMeta* meta = catalog_v2_.find_genome(id);
        if (!meta || is_deleted(meta->genome_id)) return std::nullopt;

        auto& shard = get_shard_v2(meta->shard_id);
        ExtractedGenome eg;
        eg.meta  = *meta;
        eg.fasta = shard.fetch_genome(id);

        auto acc_it = genome_accession_map_.find(id);
        if (acc_it != genome_accession_map_.end())
            eg.accession = acc_it->second;

        return eg;
    }

    std::vector<ExtractedGenome> extract_v2(const ExtractQuery& q) const {
        auto ptrs = catalog_v2_.filter(q);

        std::unordered_map<uint32_t, std::vector<const GenomeMeta*>> by_shard;
        for (const auto* p : ptrs) {
            if (!is_deleted(p->genome_id))
                by_shard[p->shard_id].push_back(p);
        }

        std::vector<ExtractedGenome> out;
        out.reserve(ptrs.size());

        for (auto& [shard_id, metas] : by_shard) {
            auto& shard = get_shard_v2(shard_id);
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
        return out;
    }
};

// ── ArchiveReader public API ──────────────────────────────────────────────────

ArchiveReader::ArchiveReader()  : impl_(std::make_unique<Impl>()) {}
ArchiveReader::~ArchiveReader() = default;

void ArchiveReader::open(const std::filesystem::path& path) { impl_->open(path); }

void ArchiveReader::close() {
    impl_->mmap_.close();
    impl_->shards_v2_.clear();
    impl_->accession_map_.clear();
    impl_->genome_accession_map_.clear();
    impl_->tombstones_.clear();
    impl_->catalog_v1_.close();
    impl_->shards_v1_.clear();
    impl_->open_  = false;
    impl_->is_v2_ = false;
}

bool ArchiveReader::is_open() const { return impl_->open_; }

std::optional<GenomeMeta> ArchiveReader::genome_meta(GenomeId id) const {
    if (impl_->is_v2_) {
        const GenomeMeta* p = impl_->catalog_v2_.find_genome(id);
        if (!p) return std::nullopt;
        return *p;
    }
    const GenomeMeta* p = impl_->catalog_v1_.find_genome(id);
    if (!p) return std::nullopt;
    return *p;
}

size_t ArchiveReader::count(const ExtractQuery& q) const {
    if (impl_->is_v2_) return impl_->catalog_v2_.filter(q).size();
    return impl_->catalog_v1_.filter(q).size();
}

std::vector<GenomeMeta> ArchiveReader::filter_meta(const ExtractQuery& q) const {
    std::vector<const GenomeMeta*> ptrs;
    if (impl_->is_v2_) {
        ptrs = impl_->catalog_v2_.filter(q);
    } else {
        ptrs = impl_->catalog_v1_.filter(q);
    }
    std::vector<GenomeMeta> out;
    out.reserve(ptrs.size());
    for (const auto* p : ptrs) out.push_back(*p);
    return out;
}

std::optional<ExtractedGenome> ArchiveReader::fetch_genome(GenomeId id) const {
    if (impl_->is_v2_) return impl_->fetch_genome_v2(id);

    const GenomeMeta* meta = impl_->catalog_v1_.find_genome(id);
    if (!meta || meta->is_deleted()) return std::nullopt;

    auto& shard = impl_->get_shard_v1(meta->shard_id);
    ExtractedGenome eg;
    eg.meta  = *meta;
    eg.fasta = shard.fetch_genome(id);

    auto acc_it = impl_->genome_accession_map_.find(id);
    if (acc_it != impl_->genome_accession_map_.end())
        eg.accession = acc_it->second;

    return eg;
}

std::optional<ExtractedGenome> ArchiveReader::fetch_by_accession(
    std::string_view accession) const
{
    auto it = impl_->accession_map_.find(std::string(accession));
    if (it == impl_->accession_map_.end()) return std::nullopt;
    return fetch_genome(it->second);
}

std::vector<ExtractedGenome> ArchiveReader::extract(const ExtractQuery& q) const {
    if (impl_->is_v2_) return impl_->extract_v2(q);

    auto ptrs = impl_->catalog_v1_.filter(q);
    std::vector<ExtractedGenome> out;
    out.reserve(ptrs.size());

    std::unordered_map<ShardId, std::vector<const GenomeMeta*>> by_shard;
    for (const auto* p : ptrs) by_shard[p->shard_id].push_back(p);

    for (auto& [shard_id, metas] : by_shard) {
        auto& shard = impl_->get_shard_v1(shard_id);
        for (const auto* m : metas) {
            ExtractedGenome eg;
            eg.meta  = *m;
            eg.fasta = shard.fetch_genome(m->genome_id);

            auto acc_it = impl_->genome_accession_map_.find(m->genome_id);
            if (acc_it != impl_->genome_accession_map_.end())
                eg.accession = acc_it->second;

            out.push_back(std::move(eg));
        }
    }
    return out;
}

ArchiveReader::ArchiveStats ArchiveReader::archive_stats() const {
    ArchiveStats s{};
    s.generation      = impl_->generation_;
    s.n_shards        = impl_->n_shards_;
    s.n_genomes_total = impl_->total_count_;
    s.n_genomes_live  = impl_->live_count_;

    if (impl_->is_v2_) {
        impl_->catalog_v2_.scan([&](const GenomeMeta& m) {
            if (!m.is_deleted() && !impl_->is_deleted(m.genome_id)) {
                s.total_raw_bp           += m.genome_length;
                s.total_compressed_bytes += m.blob_len_cmp;
            }
            return true;
        });
    } else {
        impl_->catalog_v1_.scan([&](const GenomeMeta& m) {
            if (!m.is_deleted()) {
                s.total_raw_bp           += m.genome_length;
                s.total_compressed_bytes += m.blob_len_cmp;
            }
            return true;
        });
    }

    s.compression_ratio = (s.total_compressed_bytes > 0)
        ? static_cast<double>(s.total_raw_bp) / s.total_compressed_bytes
        : 0.0;
    return s;
}

} // namespace genopack
