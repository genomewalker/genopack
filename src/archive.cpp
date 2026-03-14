#include <genopack/archive.hpp>
#include <spdlog/spdlog.h>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace genopack {

// ── ArchiveReader::Impl ───────────────────────────────────────────────────────

struct ArchiveReader::Impl {
    std::filesystem::path archive_dir;
    ManifestHeader        manifest{};
    CatalogReader         catalog;

    // Accession <-> GenomeId maps loaded from meta.tsv sidecar (optional)
    std::unordered_map<std::string, GenomeId> accession_map;
    std::unordered_map<GenomeId, std::string> genome_accession_map;

    // Lazy-loaded shard readers (cache by shard_id)
    mutable std::unordered_map<ShardId, ShardReader> shards;

    bool open_ = false;

    void load_meta_tsv(const std::filesystem::path& path) {
        std::ifstream f(path);
        if (!f) {
            spdlog::debug("genopack: meta.tsv not found at {}, accession lookup unavailable",
                          path.string());
            return;
        }

        std::string line;
        // Skip header
        if (!std::getline(f, line)) return;

        while (std::getline(f, line)) {
            if (line.empty()) continue;
            std::istringstream ss(line);
            std::string accession, genome_id_str;
            if (!std::getline(ss, accession, '\t')) continue;
            if (!std::getline(ss, genome_id_str, '\t')) continue;

            GenomeId gid = std::stoull(genome_id_str);
            accession_map[accession]   = gid;
            genome_accession_map[gid]  = accession;
        }
        spdlog::debug("genopack: loaded {} accession mappings from meta.tsv",
                      accession_map.size());
    }

    void open(const std::filesystem::path& dir) {
        archive_dir = dir;

        // Read manifest
        std::ifstream mf(dir / "MANIFEST.bin", std::ios::binary);
        if (!mf) throw std::runtime_error("Cannot open MANIFEST.bin in " + dir.string());
        mf.read(reinterpret_cast<char*>(&manifest), sizeof(manifest));
        if (manifest.magic != GPKM_MAGIC)
            throw std::runtime_error("Invalid MANIFEST.bin magic");

        // Open catalog
        catalog.open(dir / "catalog.gpkc");

        // Load meta.tsv sidecar if present
        load_meta_tsv(dir / "meta.tsv");

        open_ = true;
        spdlog::info("genopack opened: {} genomes, {} shards",
                     manifest.n_genomes_live, manifest.n_shards);
    }

    ShardReader& get_shard(ShardId shard_id) const {
        auto it = shards.find(shard_id);
        if (it != shards.end()) return it->second;

        // Scan shards/ directory for any .gpks file with matching shard_id
        for (const auto& entry : std::filesystem::directory_iterator(archive_dir / "shards")) {
            if (entry.path().extension() != ".gpks") continue;
            ShardReader tmp;
            tmp.open(entry.path());
            if (tmp.shard_id() == shard_id) {
                auto [ins_it, ok] = shards.emplace(shard_id, std::move(tmp));
                return ins_it->second;
            }
        }
        throw std::runtime_error("Shard " + std::to_string(shard_id) + " not found");
    }
};

ArchiveReader::ArchiveReader()  : impl_(std::make_unique<Impl>()) {}
ArchiveReader::~ArchiveReader() = default;

void ArchiveReader::open(const std::filesystem::path& dir) { impl_->open(dir); }
void ArchiveReader::close() {
    impl_->catalog.close();
    impl_->shards.clear();
    impl_->accession_map.clear();
    impl_->genome_accession_map.clear();
    impl_->open_ = false;
}
bool ArchiveReader::is_open() const { return impl_->open_; }

std::optional<GenomeMeta> ArchiveReader::genome_meta(GenomeId id) const {
    const GenomeMeta* p = impl_->catalog.find_genome(id);
    if (!p) return std::nullopt;
    return *p;
}

size_t ArchiveReader::count(const ExtractQuery& q) const {
    return impl_->catalog.filter(q).size();
}

std::vector<GenomeMeta> ArchiveReader::filter_meta(const ExtractQuery& q) const {
    auto ptrs = impl_->catalog.filter(q);
    std::vector<GenomeMeta> out;
    out.reserve(ptrs.size());
    for (const auto* p : ptrs) out.push_back(*p);
    return out;
}

std::optional<ExtractedGenome> ArchiveReader::fetch_genome(GenomeId id) const {
    const GenomeMeta* meta = impl_->catalog.find_genome(id);
    if (!meta || meta->is_deleted()) return std::nullopt;

    auto& shard = impl_->get_shard(meta->shard_id);
    ExtractedGenome eg;
    eg.meta  = *meta;
    eg.fasta = shard.fetch_genome(id);

    auto acc_it = impl_->genome_accession_map.find(id);
    if (acc_it != impl_->genome_accession_map.end())
        eg.accession = acc_it->second;

    return eg;
}

std::optional<ExtractedGenome> ArchiveReader::fetch_by_accession(std::string_view accession) const {
    auto it = impl_->accession_map.find(std::string(accession));
    if (it == impl_->accession_map.end()) return std::nullopt;
    return fetch_genome(it->second);
}

std::vector<ExtractedGenome> ArchiveReader::extract(const ExtractQuery& q) const {
    auto ptrs = impl_->catalog.filter(q);
    std::vector<ExtractedGenome> out;
    out.reserve(ptrs.size());

    std::unordered_map<ShardId, std::vector<const GenomeMeta*>> by_shard;
    for (const auto* p : ptrs) by_shard[p->shard_id].push_back(p);

    for (auto& [shard_id, metas] : by_shard) {
        auto& shard = impl_->get_shard(shard_id);
        for (const auto* m : metas) {
            ExtractedGenome eg;
            eg.meta  = *m;
            eg.fasta = shard.fetch_genome(m->genome_id);

            auto acc_it = impl_->genome_accession_map.find(m->genome_id);
            if (acc_it != impl_->genome_accession_map.end())
                eg.accession = acc_it->second;

            out.push_back(std::move(eg));
        }
    }
    return out;
}

ArchiveReader::ArchiveStats ArchiveReader::archive_stats() const {
    ArchiveStats s{};
    s.generation      = impl_->manifest.generation;
    s.n_shards        = impl_->manifest.n_shards;
    s.n_genomes_total = impl_->manifest.n_genomes;
    s.n_genomes_live  = impl_->manifest.n_genomes_live;

    impl_->catalog.scan([&](const GenomeMeta& m) {
        if (!m.is_deleted()) {
            s.total_raw_bp           += m.genome_length;
            s.total_compressed_bytes += m.blob_len_cmp;
        }
        return true;
    });
    s.compression_ratio = (s.total_compressed_bytes > 0)
        ? static_cast<double>(s.total_raw_bp) / s.total_compressed_bytes
        : 0.0;
    return s;
}

} // namespace genopack
