#include <genopack/archive.hpp>
#include <genopack/catalog.hpp>
#include <genopack/shard.hpp>
#include <spdlog/spdlog.h>
#include <zlib.h>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <ctime>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace genopack {

// ── Helpers ───────────────────────────────────────────────────────────────────
// These replicate builder.cpp's helpers; consolidate in a shared header later.

static uint32_t days_since_epoch_app() {
    using namespace std::chrono;
    auto now = system_clock::now();
    auto tp  = floor<days>(now);
    auto ref = sys_days{year{2024}/January/1};
    return static_cast<uint32_t>((tp - ref).count());
}

static std::string decompress_gz_app(const std::filesystem::path& path) {
    gzFile gz = gzopen(path.string().c_str(), "rb");
    if (!gz) {
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
    int n;
    while ((n = gzread(gz, buf, sizeof(buf))) > 0)
        out.append(buf, n);
    gzclose(gz);
    return out;
}

static uint64_t genome_minhash_app(const std::string& fasta, int k = 21) {
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
            v ^= v >> 33; v *= 0xff51afd7ed558ccdULL;
            v ^= v >> 33; v *= 0xc4ceb9fe1a85ec53ULL;
            v ^= v >> 33;
            if (v < min_hash) min_hash = v;
        }
    }
    return min_hash;
}

struct AppFastaStats {
    uint64_t genome_length   = 0;
    uint32_t n_contigs       = 0;
    uint16_t gc_pct_x100     = 0;
    uint64_t oph_fingerprint = 0;
};

static AppFastaStats compute_app_stats(const std::string& fasta) {
    AppFastaStats s;
    uint64_t gc = 0, at = 0;
    for (char c : fasta) {
        if (c == '>') { ++s.n_contigs; continue; }
        if (c == '\n' || c == '\r') continue;
        ++s.genome_length;
        char u = static_cast<char>(c & ~0x20);
        if      (u == 'G' || u == 'C') ++gc;
        else if (u == 'A' || u == 'T') ++at;
    }
    uint64_t total = gc + at;
    s.gc_pct_x100    = total > 0 ? static_cast<uint16_t>(gc * 10000 / total) : 0;
    s.oph_fingerprint = genome_minhash_app(fasta);
    return s;
}

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
    std::filesystem::path    archive_dir;
    std::vector<std::string> tombstone_accessions;
    std::vector<GenomeId>    tombstone_ids;
    std::vector<BuildRecord> pending;

    explicit Impl(const std::filesystem::path& dir) : archive_dir(dir) {
        if (!std::filesystem::exists(dir / "MANIFEST.bin"))
            throw std::runtime_error("Not a genopack archive: " + dir.string());
    }

    void commit() {
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
            BuildRecord  record;
            AppFastaStats stats;
            std::string  fasta;
            GenomeMeta   meta;
        };

        std::vector<GenomeWork> work;
        work.reserve(pending.size());

        for (auto& r : pending) {
            GenomeWork gw;
            gw.record = std::move(r);
            try {
                gw.fasta = decompress_gz_app(gw.record.file_path);
            } catch (const std::exception& ex) {
                spdlog::warn("Skipping {}: {}", gw.record.accession, ex.what());
                continue;
            }
            gw.stats = compute_app_stats(gw.fasta);

            gw.meta = GenomeMeta{};
            gw.meta.genome_id         = next_id++;
            gw.meta._reserved0        = 0;
            gw.meta.shard_id          = INVALID_SHARD_ID; // filled in during write pass
            gw.meta.genome_length     = gw.record.genome_length > 0 ? gw.record.genome_length : gw.stats.genome_length;
            gw.meta.n_contigs         = gw.record.n_contigs > 0     ? gw.record.n_contigs     : gw.stats.n_contigs;
            gw.meta.gc_pct_x100       = gw.stats.gc_pct_x100;
            gw.meta.completeness_x10  = static_cast<uint16_t>(gw.record.completeness  * 10.0f);
            gw.meta.contamination_x10 = static_cast<uint16_t>(gw.record.contamination * 10.0f);
            gw.meta.oph_fingerprint   = gw.stats.oph_fingerprint;
            gw.meta.date_added        = days_since_epoch_app();
            gw.meta.flags             = 0;

            work.push_back(std::move(gw));
        }
        pending.clear();

        // ── 6. Sort by oph_fingerprint (maximises zstd LDM reuse) ───────────
        std::sort(work.begin(), work.end(), [](const GenomeWork& a, const GenomeWork& b) {
            return a.meta.oph_fingerprint < b.meta.oph_fingerprint;
        });

        // ── 7. Write new shard(s) ────────────────────────────────────────────
        std::vector<std::pair<std::string, GenomeId>> new_accessions; // for meta.tsv
        std::vector<GenomeMeta>                       new_rows;

        if (!work.empty()) {
            std::filesystem::create_directories(archive_dir / "shards");

            ShardWriterConfig shard_cfg{};
            std::unique_ptr<ShardWriter> sw;
            size_t sw_compressed = 0;

            auto open_shard = [&]() {
                char fname[64];
                snprintf(fname, sizeof(fname), "shard_%05u.gpks", next_shard_id);
                sw = std::make_unique<ShardWriter>(
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

        // CatalogWriter::finalize() re-sorts internally, but we sort here too
        // so the merge is correct before handing off.
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
};

// ── Public interface ──────────────────────────────────────────────────────────

ArchiveAppender::ArchiveAppender(const std::filesystem::path& dir)
    : impl_(std::make_unique<Impl>(dir))
{}
ArchiveAppender::~ArchiveAppender() = default;

void ArchiveAppender::add_from_tsv(const std::filesystem::path& tsv_path) {
    std::ifstream f(tsv_path);
    if (!f) throw std::runtime_error("Cannot open TSV: " + tsv_path.string());

    std::string header_line;
    std::getline(f, header_line);

    std::vector<std::string> cols;
    {
        std::istringstream ss(header_line);
        std::string tok;
        while (std::getline(ss, tok, '\t')) cols.push_back(tok);
    }

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

    // Track extra column indices (all beyond the six known ones)
    std::set<int> known;
    if (idx_acc  >= 0) known.insert(idx_acc);
    if (idx_path >= 0) known.insert(idx_path);
    if (idx_comp >= 0) known.insert(idx_comp);
    if (idx_cont >= 0) known.insert(idx_cont);
    if (idx_glen >= 0) known.insert(idx_glen);
    if (idx_ctg  >= 0) known.insert(idx_ctg);

    std::vector<int>         extra_indices;
    std::vector<std::string> extra_names;
    for (int i = 0; i < (int)cols.size(); ++i) {
        if (!known.count(i)) { extra_indices.push_back(i); extra_names.push_back(cols[i]); }
    }

    std::string line;
    size_t before = impl_->pending.size();
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> fields;
        {
            std::istringstream ss(line);
            std::string tok;
            while (std::getline(ss, tok, '\t')) fields.push_back(tok);
        }
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
        impl_->pending.push_back(std::move(r));
    }

    spdlog::info("Appender loaded {} records from {}",
                 impl_->pending.size() - before, tsv_path.string());
}

void ArchiveAppender::add(const BuildRecord& rec) {
    impl_->pending.push_back(rec);
}

void ArchiveAppender::remove(GenomeId id) {
    impl_->tombstone_ids.push_back(id);
}

void ArchiveAppender::remove_by_accession(std::string_view accession) {
    impl_->tombstone_accessions.emplace_back(accession);
    spdlog::info("Queued tombstone: {}", accession);
}

void ArchiveAppender::commit() {
    impl_->commit();
}

} // namespace genopack
