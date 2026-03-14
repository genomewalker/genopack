#include <genopack/archive.hpp>
#include <genopack/shard.hpp>
#include <genopack/catalog.hpp>
#include <zlib.h>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <numeric>
#include <set>
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

// MinHash minimum over k=21 k-mers. Locality-sensitive: similar genomes → nearby values.
static uint64_t genome_minhash(const std::string& fasta, int k = 21) {
    uint64_t min_hash = UINT64_MAX;
    uint64_t kmer = 0;
    int kmer_len = 0;
    for (char c : fasta) {
        if (c == '>' || c == '\n' || c == '\r') { kmer_len = 0; kmer = 0; continue; }
        char u = static_cast<char>(c & ~0x20);
        uint8_t base;
        if      (u == 'A') base = 0;
        else if (u == 'C') base = 1;
        else if (u == 'G') base = 2;
        else if (u == 'T') base = 3;
        else { kmer_len = 0; kmer = 0; continue; }
        kmer = (kmer << 2) | base;
        if (++kmer_len >= k) {
            uint64_t v = kmer & ((1ULL << (2 * k)) - 1);
            // MurmurHash3 finalizer mix
            v ^= v >> 33;
            v *= 0xff51afd7ed558ccdULL;
            v ^= v >> 33;
            v *= 0xc4ceb9fe1a85ec53ULL;
            v ^= v >> 33;
            if (v < min_hash) min_hash = v;
        }
    }
    return min_hash;
}

// Quick genome stats from FASTA string
struct FastaStats {
    uint64_t genome_length   = 0;
    uint32_t n_contigs       = 0;
    uint16_t gc_pct_x100     = 0;
    uint64_t oph_fingerprint = 0;  // MinHash minimum (k=21)
};

static FastaStats compute_fasta_stats(const std::string& fasta) {
    FastaStats s;
    uint64_t gc = 0;
    uint64_t at = 0;

    for (char c : fasta) {
        if (c == '>') {
            ++s.n_contigs;
        } else if (c != '\n' && c != '\r') {
            ++s.genome_length;
            char u = static_cast<char>(c & ~0x20);  // uppercase
            if (u == 'G' || u == 'C') ++gc;
            else if (u == 'A' || u == 'T') ++at;
        }
    }
    s.oph_fingerprint = genome_minhash(fasta);
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

        // Parse tab-separated header
        std::vector<std::string> cols;
        {
            std::istringstream ss(header_line);
            std::string tok;
            while (std::getline(ss, tok, '\t')) cols.push_back(tok);
        }

        // Find required column indices
        auto find_col = [&](std::initializer_list<const char*> names) -> int {
            for (const char* name : names)
                for (int i = 0; i < (int)cols.size(); ++i)
                    if (cols[i] == name) return i;
            return -1;
        };

        int idx_acc  = find_col({"accession", "acc"});
        int idx_path = find_col({"file_path", "path", "fasta_path", "fasta"});
        int idx_comp = find_col({"completeness"});
        int idx_cont = find_col({"contamination"});
        int idx_glen = find_col({"genome_length"});
        int idx_ctg  = find_col({"n_contigs"});

        if (idx_acc  < 0) throw std::runtime_error("TSV missing 'accession' column");
        if (idx_path < 0) throw std::runtime_error("TSV missing 'file_path' column");

        // Identify extra column indices (all except the known ones)
        std::vector<int> extra_indices;
        std::vector<std::string> extra_names;
        std::set<int> known = {};
        if (idx_acc  >= 0) known.insert(idx_acc);
        if (idx_path >= 0) known.insert(idx_path);
        if (idx_comp >= 0) known.insert(idx_comp);
        if (idx_cont >= 0) known.insert(idx_cont);
        if (idx_glen >= 0) known.insert(idx_glen);
        if (idx_ctg  >= 0) known.insert(idx_ctg);
        for (int i = 0; i < (int)cols.size(); ++i) {
            if (!known.count(i)) { extra_indices.push_back(i); extra_names.push_back(cols[i]); }
        }

        std::string line;
        while (std::getline(f, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::vector<std::string> fields;
            std::istringstream ss(line);
            std::string tok;
            while (std::getline(ss, tok, '\t')) fields.push_back(tok);
            if ((int)fields.size() <= std::max(idx_acc, idx_path)) continue;

            BuildRecord r;
            r.accession = fields[idx_acc];
            r.file_path = fields[idx_path];
            if (idx_comp >= 0 && idx_comp < (int)fields.size()) r.completeness  = std::stof(fields[idx_comp]);
            if (idx_cont >= 0 && idx_cont < (int)fields.size()) r.contamination = std::stof(fields[idx_cont]);
            if (idx_glen >= 0 && idx_glen < (int)fields.size()) r.genome_length  = std::stoull(fields[idx_glen]);
            if (idx_ctg  >= 0 && idx_ctg  < (int)fields.size()) r.n_contigs      = std::stoul(fields[idx_ctg]);
            for (int j = 0; j < (int)extra_indices.size(); ++j) {
                int ci = extra_indices[j];
                r.extra_fields.emplace_back(extra_names[j],
                    ci < (int)fields.size() ? fields[ci] : "");
            }
            pending.push_back(std::move(r));
        }
        spdlog::info("Loaded {} records from {}", pending.size(), tsv_path.string());
    }

    void finalize() {
        if (pending.empty()) { spdlog::warn("No records to build"); return; }

        struct GenomeWork {
            BuildRecord  record;
            GenomeId     genome_id;
            std::string  fasta;
            FastaStats   stats;
        };

        // Phase 1: Decompress + compute MinHash stats
        std::vector<GenomeWork> work;
        work.reserve(pending.size());
        for (auto& r : pending) {
            GenomeWork gw;
            gw.record    = std::move(r);
            gw.genome_id = next_genome_id++;
            try {
                gw.fasta = decompress_gz(gw.record.file_path);
            } catch (const std::exception& ex) {
                spdlog::warn("Skipping {}: {}", gw.record.accession, ex.what());
                continue;
            }
            gw.stats = compute_fasta_stats(gw.fasta);
            work.push_back(std::move(gw));
        }
        pending.clear();

        // Phase 2: Sort by MinHash (oph_fingerprint) — locality-sensitive, approximates phylogenetic order
        std::sort(work.begin(), work.end(), [](const auto& a, const auto& b) {
            return a.stats.oph_fingerprint < b.stats.oph_fingerprint;
        });

        spdlog::info("MinHash sorted {} genomes", work.size());

        // Phase 3: Write shards by accumulated compressed size
        uint32_t date = days_since_epoch();
        std::vector<GenomeMeta> catalog_rows;
        catalog_rows.reserve(work.size());

        ShardId current_shard_id = 0;
        std::unique_ptr<ShardWriter> shard_writer;
        size_t shard_compressed_size = 0;

        // Write metadata sidecar
        std::ofstream meta_out(archive_dir / "meta.tsv");
        if (!work.empty()) {
            meta_out << "accession\tgenome_id";
            for (const auto& [k, v] : work[0].record.extra_fields)
                meta_out << "\t" << k;
            meta_out << "\n";
        }

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
            auto shard_path = archive_dir / "shards" / fname;
            shard_writer = std::make_unique<ShardWriter>(
                shard_path, current_shard_id, current_shard_id, cfg.shard_cfg);
            ++current_shard_id;
        };

        size_t n_done = 0;
        for (auto& gw : work) {
            bool shard_full = shard_writer &&
                (shard_compressed_size >= cfg.shard_cfg.max_shard_size_bytes);

            if (!shard_writer || shard_full) {
                flush_shard();
                open_shard();
            }

            GenomeMeta meta{};
            meta.genome_id         = gw.genome_id;
            meta._reserved0        = 0;
            meta.shard_id          = current_shard_id - 1;
            meta.genome_length     = gw.record.genome_length > 0 ? gw.record.genome_length : gw.stats.genome_length;
            meta.n_contigs         = gw.record.n_contigs > 0 ? gw.record.n_contigs : gw.stats.n_contigs;
            meta.gc_pct_x100       = gw.stats.gc_pct_x100;
            meta.completeness_x10  = static_cast<uint16_t>(gw.record.completeness  * 10.0f);
            meta.contamination_x10 = static_cast<uint16_t>(gw.record.contamination * 10.0f);
            meta.oph_fingerprint   = gw.stats.oph_fingerprint;
            meta.date_added        = date;

            shard_writer->add_genome(gw.genome_id, gw.stats.oph_fingerprint,
                                     gw.fasta.data(), gw.fasta.size());
            shard_compressed_size = shard_writer->compressed_size();

            catalog_rows.push_back(meta);

            // Write metadata sidecar row
            meta_out << gw.record.accession << "\t" << gw.genome_id;
            for (const auto& [k, v] : gw.record.extra_fields)
                meta_out << "\t" << v;
            meta_out << "\n";

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
            sd.shard_id = i;
            char fname[64];
            snprintf(fname, sizeof(fname), "shard_%05u.gpks", i);
            std::memcpy(sd.filename, fname, std::min(sizeof(fname), sizeof(sd.filename)));
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
