#include <genopack/archive.hpp>
#include "taxon_registry.hpp"
#include <genopack/shard.hpp>
#include <genopack/catalog.hpp>
#include <zlib.h>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <vector>

namespace genopack {

// Days since 2024-01-01
static uint32_t days_since_epoch() {
    using namespace std::chrono;
    auto now = system_clock::now();
    auto tp  = floor<days>(now);
    auto ref = sys_days{year{2024}/January/1};
    return static_cast<uint32_t>((tp - ref).count());
}

// Decompress gzipped file into string
static std::string decompress_gz(const std::filesystem::path& path) {
    gzFile gz = gzopen(path.string().c_str(), "rb");
    if (!gz) {
        // Try as plain FASTA
        std::ifstream f(path, std::ios::binary | std::ios::ate);
        if (!f) throw std::runtime_error("Cannot open: " + path.string());
        size_t sz = static_cast<size_t>(f.tellg());
        f.seekg(0);
        std::string out(sz, '\0');
        f.read(out.data(), sz);
        return out;
    }
    std::string out;
    char buf[65536];
    int  n;
    while ((n = gzread(gz, buf, sizeof(buf))) > 0)
        out.append(buf, n);
    gzclose(gz);
    return out;
}

// Quick genome stats from FASTA string
struct FastaStats {
    uint64_t genome_length = 0;
    uint32_t n_contigs     = 0;
    uint16_t gc_pct_x100   = 0;
    uint64_t oph_fingerprint = 0;  // FNV-1a 64-bit of first 1000 chars (cheap dedup hint)
};

static FastaStats compute_fasta_stats(const std::string& fasta) {
    FastaStats s;
    uint64_t gc = 0;
    uint64_t at = 0;
    uint64_t h  = 14695981039346656037ULL; // FNV offset basis
    size_t fp_chars = 0;

    for (char c : fasta) {
        if (c == '>') {
            ++s.n_contigs;
        } else if (c != '\n' && c != '\r') {
            ++s.genome_length;
            char u = static_cast<char>(c & ~0x20);  // uppercase
            if (u == 'G' || u == 'C') ++gc;
            else if (u == 'A' || u == 'T') ++at;
            if (fp_chars < 1000) {
                h = (h ^ static_cast<uint8_t>(c)) * 1099511628211ULL;
                ++fp_chars;
            }
        }
    }
    s.oph_fingerprint = h;
    uint64_t total = gc + at;
    s.gc_pct_x100 = (total > 0)
        ? static_cast<uint16_t>(gc * 10000 / total)
        : 0;
    return s;
}

// ── ArchiveBuilder::Impl ──────────────────────────────────────────────────────

struct ArchiveBuilder::Impl {
    std::filesystem::path archive_dir;
    Config                cfg;

    TaxonRegistry         registry;
    std::vector<BuildRecord> pending;
    GenomeId next_genome_id = 1;

    explicit Impl(const std::filesystem::path& dir, Config c)
        : archive_dir(dir), cfg(c)
    {
        std::filesystem::create_directories(dir / "shards");
    }

    void add(const BuildRecord& rec) {
        pending.push_back(rec);
    }

    void add_from_tsv(const std::filesystem::path& tsv_path) {
        std::ifstream f(tsv_path);
        if (!f) throw std::runtime_error("Cannot open TSV: " + tsv_path.string());
        std::string header_line;
        std::getline(f, header_line);

        std::string line;
        while (std::getline(f, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            BuildRecord r;
            std::string completeness, contamination;
            std::string path_str;
            if (!(ss >> r.accession >> r.taxonomy >> path_str)) continue;
            r.file_path = path_str;
            if (ss >> completeness) r.completeness = std::stof(completeness);
            if (ss >> contamination) r.contamination = std::stof(contamination);
            pending.push_back(std::move(r));
        }
        spdlog::info("Loaded {} records from {}", pending.size(), tsv_path.string());
    }

    void finalize() {
        if (pending.empty()) {
            spdlog::warn("No records to build");
            return;
        }

        // Register all taxa
        for (const auto& r : pending)
            registry.insert_lineage(r.taxonomy);
        registry.finalize();
        spdlog::info("Taxonomy: {} unique taxa", registry.size());

        // Assign genome IDs and resolve taxon IDs
        struct GenomeWork {
            BuildRecord  record;
            GenomeId     genome_id;
            TaxonId      taxon_id;
            TaxonId      genus_taxon_id;
        };
        std::vector<GenomeWork> work;
        work.reserve(pending.size());
        for (auto& r : pending) {
            TaxonId species_id = registry.find_by_lineage(r.taxonomy);
            // Find genus (parent of species)
            TaxonId genus_id = species_id;
            if (species_id != INVALID_TAXON_ID) {
                const auto& info = registry.info(species_id);
                genus_id = info.parent_id;
            }
            work.push_back({std::move(r), next_genome_id++, species_id, genus_id});
        }
        pending.clear();

        // Sort by taxon_id (genus first, then species, then oph_fingerprint after stats)
        std::sort(work.begin(), work.end(), [](const auto& a, const auto& b) {
            if (a.genus_taxon_id != b.genus_taxon_id) return a.genus_taxon_id < b.genus_taxon_id;
            return a.taxon_id < b.taxon_id;
        });

        // Group into genus-level shard batches
        // Each shard targets cfg.shard_cfg.max_shard_size_bytes
        std::vector<GenomeMeta> catalog_rows;
        catalog_rows.reserve(work.size());

        ShardId current_shard_id = 0;
        std::unique_ptr<ShardWriter> shard_writer;
        size_t shard_compressed_size = 0;
        TaxonId current_genus = INVALID_TAXON_ID;
        uint32_t date = days_since_epoch();

        auto flush_shard = [&]() {
            if (shard_writer) {
                shard_writer->finalize();
                shard_writer.reset();
                shard_compressed_size = 0;
            }
        };

        auto open_shard = [&](TaxonId genus_id) {
            char fname[128];
            snprintf(fname, sizeof(fname), "genus_%08u_%05u.gpks",
                     genus_id, current_shard_id);
            auto shard_path = archive_dir / "shards" / fname;
            shard_writer = std::make_unique<ShardWriter>(
                shard_path, current_shard_id, genus_id, cfg.shard_cfg);
            ++current_shard_id;
            current_genus = genus_id;
        };

        size_t n_done = 0;
        for (auto& gw : work) {
            // Open new shard when genus changes or shard is too large
            bool new_genus = (gw.genus_taxon_id != current_genus);
            bool shard_full = (shard_compressed_size >= cfg.shard_cfg.max_shard_size_bytes);

            if (!shard_writer || new_genus || shard_full) {
                flush_shard();
                open_shard(gw.genus_taxon_id);
            }

            // Decompress FASTA
            std::string fasta;
            try {
                fasta = decompress_gz(gw.record.file_path);
            } catch (const std::exception& ex) {
                spdlog::warn("Skipping {}: {}", gw.record.accession, ex.what());
                continue;
            }

            FastaStats stats = compute_fasta_stats(fasta);

            // Record metadata
            GenomeMeta meta{};
            meta.genome_id         = gw.genome_id;
            meta.taxon_id          = gw.taxon_id;
            meta.shard_id          = shard_writer ? (current_shard_id - 1) : INVALID_SHARD_ID;
            meta.genome_length     = (gw.record.genome_length > 0)
                                     ? gw.record.genome_length : stats.genome_length;
            meta.n_contigs         = (gw.record.n_contigs > 0)
                                     ? gw.record.n_contigs : stats.n_contigs;
            meta.gc_pct_x100       = stats.gc_pct_x100;
            meta.completeness_x10  = static_cast<uint16_t>(gw.record.completeness * 10.0f);
            meta.contamination_x10 = static_cast<uint16_t>(gw.record.contamination * 10.0f);
            meta.oph_fingerprint   = stats.oph_fingerprint;
            meta.date_added        = date;

            // Add to shard writer (blob_offset set after all genomes written)
            shard_writer->add_genome(gw.genome_id, gw.taxon_id, stats.oph_fingerprint,
                                     fasta.data(), fasta.size());
            shard_compressed_size = shard_writer->compressed_size();

            catalog_rows.push_back(meta);

            if (cfg.verbose && ++n_done % 1000 == 0)
                spdlog::info("Built {}/{} genomes", n_done, work.size());
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

        // Write manifest
        write_manifest(current_shard_id, catalog_rows.size(), date);

        spdlog::info("genopack archive built: {}", archive_dir.string());
    }

    void write_manifest(uint32_t n_shards, uint64_t n_genomes, uint32_t date) {
        ManifestHeader hdr{};
        hdr.magic          = GPKM_MAGIC;
        hdr.version        = FORMAT_VERSION;
        hdr.generation     = 1;
        hdr.created_at     = static_cast<uint64_t>(std::time(nullptr));
        hdr.n_shards       = n_shards;
        hdr.n_taxons       = static_cast<uint32_t>(registry.size());
        hdr.n_genomes      = n_genomes;
        hdr.n_genomes_live = n_genomes;

        std::ofstream f(archive_dir / "MANIFEST.bin", std::ios::binary | std::ios::trunc);
        if (!f) throw std::runtime_error("Cannot write MANIFEST.bin");
        f.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));

        // Write shard table (filenames)
        // v1: sequential genus_XXXXXXXX_XXXXX.gpks files
        for (uint32_t i = 0; i < n_shards; ++i) {
            ShardDescriptor sd{};
            sd.shard_id = i;
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
