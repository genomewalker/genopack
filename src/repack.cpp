#include <genopack/repack.hpp>
#include <genopack/accx.hpp>
#include <genopack/catalog.hpp>
#include <genopack/cidx.hpp>
#include <genopack/format.hpp>
#include <genopack/gidx.hpp>
#include <genopack/kmrx.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/shard.hpp>
#include <genopack/taxn.hpp>
#include <genopack/toc.hpp>
#include <genopack/txdb.hpp>
#include <spdlog/spdlog.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <condition_variable>
#include <cstring>
#include <ctime>
#include <fcntl.h>
#include <filesystem>
#include <future>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <string>
#include <thread>
#include <unistd.h>
#include <unordered_map>
#include <vector>

namespace genopack {

static std::string taxonomy_bucket_key(std::string_view taxonomy, char rank)
{
    if (taxonomy.empty()) return "__unclassified__";
    char prefix[4] = {rank, '_', '_', '\0'};
    auto pos = taxonomy.find(prefix);
    if (pos == std::string_view::npos) {
        const char* fallbacks = (rank == 'g') ? "foc" : "oc";
        for (char fb : std::string_view(fallbacks)) {
            char fp[4] = {fb, '_', '_', '\0'};
            pos = taxonomy.find(fp);
            if (pos != std::string_view::npos) break;
        }
        if (pos == std::string_view::npos) return "__unclassified__";
    }
    size_t end = taxonomy.find(';', pos + 3);
    if (end == std::string_view::npos) end = taxonomy.size();
    return std::string(taxonomy.substr(pos, end - pos));
}

void repack_archive(const std::filesystem::path& input_gpk,
                    const std::filesystem::path& output_gpk,
                    const RepackConfig& cfg)
{
    // ── Open source archive ─────────────────────────────────────────────────
    MmapFileReader src_mmap;
    src_mmap.open(input_gpk);
    if (src_mmap.size() < sizeof(FileHeader))
        throw std::runtime_error("repack: file too small: " + input_gpk.string());
    if (src_mmap.ptr_at<FileHeader>(0)->magic != GPK2_MAGIC)
        throw std::runtime_error("repack: not a v2 .gpk archive: " + input_gpk.string());

    Toc src_toc = TocReader::read(src_mmap);
    auto src_shards = src_toc.find_by_type(SEC_SHRD);

    // Sort by file offset for sequential NFS access
    std::sort(src_shards.begin(), src_shards.end(), [](const SectionDesc* a, const SectionDesc* b) {
        return a->file_offset < b->file_offset;
    });

    spdlog::info("Repack: source {} ({} shards)", input_gpk.string(), src_shards.size());

    // ── Load catalog, ACCX, TAXN ────────────────────────────────────────────
    std::unordered_map<GenomeId, GenomeMeta> gid_to_meta;
    {
        MergedCatalogReader cat;
        for (auto* sd : src_toc.find_by_type(SEC_CATL))
            cat.add_fragment(src_mmap.data(), sd->file_offset, sd->compressed_size);
        cat.scan([&](const GenomeMeta& m) {
            if (!m.is_deleted()) gid_to_meta.emplace(m.genome_id, m);
            return true;
        });
    }
    spdlog::info("Repack: {} live genomes", gid_to_meta.size());

    std::unordered_map<GenomeId, std::string> gid_to_acc;
    gid_to_acc.reserve(gid_to_meta.size());
    for (auto* sd : src_toc.find_by_type(SEC_ACCX)) {
        AccessionIndexReader air;
        air.open(src_mmap.data(), sd->file_offset, sd->compressed_size);
        air.scan([&](std::string_view acc, GenomeId gid) {
            if (gid_to_meta.count(gid)) gid_to_acc.emplace(gid, std::string(acc));
        });
    }

    std::unordered_map<std::string, std::string> acc_to_tax;
    for (auto* sd : src_toc.find_by_type(SEC_TAXN)) {
        TaxonomyIndexReader tir;
        tir.open(src_mmap.data(), sd->file_offset, sd->compressed_size);
        tir.scan([&](std::string_view acc, std::string_view tax) {
            acc_to_tax.emplace(std::string(acc), std::string(tax));
        });
    }
    spdlog::info("Repack: {} taxonomy entries loaded", acc_to_tax.size());

    // ── Phase 1: directory scan (NO FASTA decompression) ────────────────────
    // Reads only GenomeDirEntry headers (~300 MB) instead of full blobs (~3.1 TB).
    struct GenomeRecord {
        GenomeId  genome_id;
        uint64_t  oph;
        uint32_t  src_shard_idx;
        uint32_t  dir_idx;
        GenomeMeta meta;
        uint32_t  tax_key_idx;
        uint32_t  _pad;
    };

    std::vector<std::string>  tax_keys;
    std::unordered_map<std::string, uint32_t> tax_key_map;
    std::vector<GenomeRecord> records;
    records.reserve(gid_to_meta.size());

    for (size_t s = 0; s < src_shards.size(); ++s) {
        const SectionDesc* sd = src_shards[s];
        ShardReader shard;
        shard.open(src_mmap.data(), sd->file_offset, sd->compressed_size);

        uint32_t dir_idx = 0;
        for (auto* de = shard.dir_begin(); de != shard.dir_end(); ++de, ++dir_idx) {
            GenomeId gid = static_cast<GenomeId>(de->genome_id);
            auto meta_it = gid_to_meta.find(gid);
            if (meta_it == gid_to_meta.end()) continue;

            std::string key = "__unclassified__";
            auto acc_it = gid_to_acc.find(gid);
            if (acc_it != gid_to_acc.end()) {
                auto tax_it = acc_to_tax.find(acc_it->second);
                if (tax_it != acc_to_tax.end())
                    key = taxonomy_bucket_key(tax_it->second, cfg.taxonomy_rank);
            }

            auto [it, inserted] = tax_key_map.emplace(key, static_cast<uint32_t>(tax_keys.size()));
            if (inserted) tax_keys.push_back(key);

            records.push_back({gid, de->oph_fingerprint,
                               static_cast<uint32_t>(s), dir_idx,
                               meta_it->second, it->second, 0u});
        }

        if ((s + 1) % 5000 == 0)
            spdlog::info("Phase 1: {}/{} shards scanned, {} genomes indexed",
                         s + 1, src_shards.size(), records.size());
    }
    spdlog::info("Phase 1 complete: {} genomes, {} taxa", records.size(), tax_keys.size());

    // ── Phase 2: sort by (taxonomy, oph) ────────────────────────────────────
    std::sort(records.begin(), records.end(), [](const GenomeRecord& a, const GenomeRecord& b) {
        if (a.tax_key_idx != b.tax_key_idx) return a.tax_key_idx < b.tax_key_idx;
        return a.oph < b.oph;
    });

    // Build reverse index: src_shard_idx → sorted list of record positions.
    // Because we iterate i in order after the sort, src_shard_records[s] is
    // automatically in global (tax_key, oph) order — genomes are routed to their
    // taxon writers in the correct sorted sequence.
    std::vector<std::vector<uint32_t>> src_shard_records(src_shards.size());
    for (uint32_t i = 0; i < static_cast<uint32_t>(records.size()); ++i)
        src_shard_records[records[i].src_shard_idx].push_back(i);

    spdlog::info("Phase 2 complete: sorted by ({}, oph), ready for Phase 3",
                 cfg.taxonomy_rank == 'g' ? "genus" : "family");

    // ── Open output file ────────────────────────────────────────────────────
    std::filesystem::path out_path = output_gpk;
    if (out_path.extension() != ".gpk")
        out_path = std::filesystem::path(out_path.string() + ".gpk");
    std::filesystem::create_directories(out_path.parent_path());

    AppendWriter app_writer;
    app_writer.create(out_path);

    {
        FileHeader ofh{};
        ofh.magic           = GPK2_MAGIC;
        ofh.version_major   = FORMAT_MAJOR;
        ofh.version_minor   = FORMAT_MINOR;
        uint64_t t = static_cast<uint64_t>(std::time(nullptr));
        ofh.file_uuid_lo    = t ^ 0xdeadbeefcafe0002ULL;
        ofh.file_uuid_hi    = (t << 17) ^ 0xfedcba9876543210ULL;
        ofh.created_at_unix = t;
        ofh.flags           = 0;
        std::memset(ofh.reserved, 0, sizeof(ofh.reserved));
        app_writer.append(&ofh, sizeof(ofh));
    }

    // ── Background IO writer thread ─────────────────────────────────────────
    TocWriter new_toc;
    uint64_t next_section_id = 1;
    ShardId  current_shard_id = 0;

    struct WriteTask {
        std::future<FrozenShard> fut;
        uint64_t section_id;
        ShardId  shard_id;
    };
    std::queue<WriteTask>   write_q;
    std::mutex              write_q_mx;
    std::condition_variable write_q_cv;
    std::condition_variable write_q_space_cv;
    bool                    writer_done = false;

    std::unordered_map<ShardId, uint64_t> shard_id_to_section_id;

    std::thread io_writer([&]() {
        while (true) {
            WriteTask wt;
            {
                std::unique_lock lk(write_q_mx);
                write_q_cv.wait(lk, [&]{ return !write_q.empty() || writer_done; });
                if (write_q.empty()) break;
                wt = std::move(write_q.front());
                write_q.pop();
            }
            write_q_space_cv.notify_one();

            FrozenShard frozen = wt.fut.get();
            uint64_t shard_start = app_writer.current_offset();
            app_writer.append(frozen.bytes.data(), frozen.bytes.size());

            SectionDesc sd{};
            sd.type              = SEC_SHRD;
            sd.version           = 4;
            sd.flags             = 0;
            sd.section_id        = wt.section_id;
            sd.file_offset       = shard_start;
            sd.compressed_size   = static_cast<uint64_t>(frozen.bytes.size());
            sd.uncompressed_size = 0;
            sd.item_count        = frozen.n_genomes;
            sd.aux0              = wt.shard_id;
            sd.aux1              = 0;
            std::memset(sd.checksum, 0, sizeof(sd.checksum));
            new_toc.add_section(sd);
        }
    });

    // ── Per-genome output data (CATL / GIDX) ────────────────────────────────
    struct GidxEntry {
        GenomeId genome_id;
        ShardId  shard_id;
        uint32_t dir_index;
        uint64_t catl_row_index;
    };
    std::vector<GenomeMeta> new_catalog;
    std::vector<GidxEntry>  new_gidx;
    new_catalog.reserve(records.size());
    new_gidx.reserve(records.size());

    // Freeze a shard writer and enqueue it for the IO thread.
    auto launch_freeze = [&](std::unique_ptr<ShardWriter>& sw, ShardId sid) {
        if (!sw || sw->n_genomes() == 0) { sw.reset(); return; }
        uint64_t sec_id = next_section_id++;
        shard_id_to_section_id[sid] = sec_id;
        WriteTask wt;
        wt.section_id = sec_id;
        wt.shard_id   = sid;
        wt.fut = std::async(std::launch::async,
            [s = std::move(sw)]() mutable { return s->freeze(); });
        {
            std::unique_lock lk(write_q_mx);
            write_q_space_cv.wait(lk, [&]{ return write_q.size() < 2; });
            write_q.push(std::move(wt));
        }
        write_q_cv.notify_one();
    };

    // ── Phase 3: decompress + route + write ─────────────────────────────────
    // Each taxon gets a ShardWriter. Flush naturally when full (512 MB).
    // When total_buffered exceeds cap: evict LARGEST writer only — not all writers.
    // This minimises shard fragmentation vs the naive "flush everything" strategy.
    struct TaxonWriter {
        std::unique_ptr<ShardWriter> writer;
        ShardId  shard_id       = 0;
        uint64_t bytes_buffered = 0;
    };
    std::unordered_map<uint32_t, TaxonWriter> taxon_writers;
    uint64_t total_buffered_bytes = 0;

    static constexpr int PREFETCH_AHEAD = 4;
    auto prefetch_shard = [&]([[maybe_unused]] size_t s) {
#ifdef MADV_WILLNEED
        if (s >= src_shards.size()) return;
        const uint8_t* ptr = src_mmap.data() + src_shards[s]->file_offset;
        size_t len = static_cast<size_t>(src_shards[s]->compressed_size);
        static const size_t PAGE = 4096;
        uintptr_t addr    = reinterpret_cast<uintptr_t>(ptr);
        uintptr_t aligned = addr & ~(PAGE - 1);
        ::madvise(reinterpret_cast<void*>(aligned), (addr + len) - aligned, MADV_WILLNEED);
#endif
    };

    for (int k = 0; k < std::min<int>(PREFETCH_AHEAD, static_cast<int>(src_shards.size())); ++k)
        prefetch_shard(static_cast<size_t>(k));

    const int n_threads = static_cast<int>(cfg.threads > 0 ? cfg.threads : 1);
    size_t n_shards_done   = 0;
    size_t n_genomes_done  = 0;

    for (size_t s = 0; s < src_shards.size(); ++s) {
        const auto& rec_idxs = src_shard_records[s];
        if (rec_idxs.empty()) { ++n_shards_done; continue; }

        prefetch_shard(s + PREFETCH_AHEAD);

        ShardReader shard;
        shard.open(src_mmap.data(), src_shards[s]->file_offset, src_shards[s]->compressed_size);

        // OMP parallel decompress — fetch_genome_at is read-only on mmap
        std::vector<std::string> fastas(rec_idxs.size());
        #pragma omp parallel for schedule(dynamic, 4) num_threads(n_threads)
        for (int j = 0; j < static_cast<int>(rec_idxs.size()); ++j)
            fastas[j] = shard.fetch_genome_at(records[rec_idxs[j]].dir_idx);

        // Route in sorted (tax_key, oph) order — rec_idxs is already sorted by
        // global position in records[], which is the (tax_key, oph) sort order.
        for (size_t j = 0; j < rec_idxs.size(); ++j) {
            const GenomeRecord& rec = records[rec_idxs[j]];
            auto& tw = taxon_writers[rec.tax_key_idx];

            if (!tw.writer) {
                tw.shard_id = current_shard_id++;
                tw.writer   = std::make_unique<ShardWriter>(tw.shard_id, tw.shard_id, cfg.shard_cfg);
            }

            tw.writer->add_genome(rec.genome_id, rec.oph, fastas[j].data(), fastas[j].size());
            uint64_t fsize = fastas[j].size();
            tw.bytes_buffered    += fsize;
            total_buffered_bytes += fsize;

            uint64_t catl_idx = new_catalog.size();
            GenomeMeta m = rec.meta;
            m.shard_id = tw.shard_id;
            new_catalog.push_back(m);
            new_gidx.push_back({rec.genome_id, tw.shard_id,
                                 static_cast<uint32_t>(tw.writer->n_genomes() - 1),
                                 catl_idx});
            ++n_genomes_done;

            // Natural flush: shard is full
            if (tw.bytes_buffered >= cfg.shard_cfg.max_shard_size_bytes) {
                total_buffered_bytes -= tw.bytes_buffered;
                launch_freeze(tw.writer, tw.shard_id);
                tw.bytes_buffered = 0;
            }
        }

        shard.release_pages();
        ++n_shards_done;

        // Smart eviction: flush only the LARGEST writer when cap is hit.
        // Compared to "flush all", this creates far fewer partial shards.
        while (total_buffered_bytes >= cfg.max_bucket_bytes) {
            TaxonWriter* largest = nullptr;
            for (auto& [kid, tw] : taxon_writers) {
                if (tw.writer && (!largest || tw.bytes_buffered > largest->bytes_buffered))
                    largest = &tw;
            }
            if (!largest) break;

            spdlog::info("Repack: cap ({:.1f} GB), evicting shard {} ({:.0f} MB)",
                         total_buffered_bytes / double(1ULL << 30),
                         largest->shard_id,
                         largest->bytes_buffered / 1e6);
            total_buffered_bytes -= largest->bytes_buffered;
            launch_freeze(largest->writer, largest->shard_id);
            largest->bytes_buffered = 0;
        }

        if (cfg.verbose || n_shards_done % 2000 == 0)
            spdlog::info("Phase 3: {}/{} shards, {} genomes, {:.1f} GB buffered",
                         n_shards_done, src_shards.size(), n_genomes_done,
                         total_buffered_bytes / double(1ULL << 30));
    }

    // Flush remaining partial writers
    for (auto& [kid, tw] : taxon_writers) {
        if (tw.writer) {
            launch_freeze(tw.writer, tw.shard_id);
            tw.bytes_buffered = 0;
        }
    }

    // Signal IO writer and wait
    {
        std::unique_lock lk(write_q_mx);
        writer_done = true;
    }
    write_q_cv.notify_one();
    io_writer.join();

    spdlog::info("Phase 3 complete: {} genomes → {} shards", n_genomes_done, current_shard_id);
    app_writer.flush();

    // ── Write metadata (CATL, GIDX, ACCX, TAXN, KMRX, TXDB/CIDX/HNSW, TOC) ──
    const uint64_t meta_base = app_writer.current_offset();
    auto meta_tmp_fs = std::filesystem::path(out_path).parent_path() / "gpk_repack_XXXXXX";
    std::string meta_tmp_s = meta_tmp_fs.string();
    std::vector<char> meta_tmp(meta_tmp_s.begin(), meta_tmp_s.end());
    meta_tmp.push_back('\0');
    {
        int fd = ::mkstemp(meta_tmp.data());
        if (fd < 0) throw std::runtime_error("Cannot create temp metadata file");
        ::close(fd);
    }
    AppendWriter mw;
    mw.create(meta_tmp.data());
    mw.seek_to(meta_base);

    uint64_t catalog_root_id   = 0;
    uint64_t accession_root_id = 0;

    // CATL
    {
        CatalogSectionWriter csw;
        for (const auto& m : new_catalog) csw.add(m);
        SectionDesc sd = csw.finalize(mw, next_section_id++);
        catalog_root_id = sd.section_id;
        new_toc.add_section(sd);
    }

    // GIDX
    {
        GidxWriter gw;
        for (const auto& e : new_gidx) {
            auto it = shard_id_to_section_id.find(e.shard_id);
            uint32_t sec_id = (it != shard_id_to_section_id.end())
                ? static_cast<uint32_t>(it->second) : 0;
            gw.add(e.genome_id, sec_id, e.dir_index, e.catl_row_index);
        }
        new_toc.add_section(gw.finalize(mw, next_section_id++));
    }

    // ACCX
    {
        AccessionIndexWriter aiw;
        for (const auto& [gid, acc] : gid_to_acc) aiw.add(acc, gid);
        SectionDesc sd = aiw.finalize(mw, next_section_id++);
        accession_root_id = sd.section_id;
        new_toc.add_section(sd);
    }

    // TAXN
    if (!acc_to_tax.empty()) {
        TaxonomyIndexWriter tiw;
        for (auto* sd : src_toc.find_by_type(SEC_TAXN)) {
            TaxonomyIndexReader tir;
            tir.open(src_mmap.data(), sd->file_offset, sd->compressed_size);
            tir.scan([&](std::string_view acc, std::string_view tax) {
                tiw.add(std::string(acc), std::string(tax));
            });
        }
        new_toc.add_section(tiw.finalize(mw, next_section_id++));
    }

    // KMRX
    {
        auto kmrx_secs = src_toc.find_by_type(SEC_KMRX);
        if (!kmrx_secs.empty()) {
            std::vector<KmrxReader> readers;
            readers.reserve(kmrx_secs.size());
            for (auto* sd : kmrx_secs) {
                KmrxReader kr;
                kr.open(src_mmap.data(), sd->file_offset, sd->compressed_size);
                readers.push_back(std::move(kr));
            }
            KmrxWriter kw;
            for (const auto& [gid, _] : gid_to_meta) {
                for (const auto& kr : readers) {
                    const float* p = kr.profile_for(gid);
                    if (p) {
                        std::array<float, 136> arr;
                        std::copy(p, p + 136, arr.begin());
                        kw.add(gid, arr);
                        break;
                    }
                }
            }
            new_toc.add_section(kw.finalize(mw, next_section_id++));
        }
    }

    // TXDB / CIDX / HNSW — copy raw bytes (genome_id indexed, unaffected by reshard)
    for (uint32_t type : {SEC_TXDB, SEC_CIDX, SEC_HNSW}) {
        for (auto* sd : src_toc.find_by_type(type)) {
            uint64_t new_offset = mw.current_offset();
            mw.append(src_mmap.data() + sd->file_offset, sd->compressed_size);
            SectionDesc new_sd  = *sd;
            new_sd.section_id   = next_section_id++;
            new_sd.file_offset  = new_offset;
            new_toc.add_section(new_sd);
        }
    }

    // TOC + TailLocator
    new_toc.finalize(mw,
                     /*generation=*/1,
                     /*live_count=*/static_cast<uint64_t>(n_genomes_done),
                     /*total_count=*/static_cast<uint64_t>(n_genomes_done),
                     /*prev_toc_offset=*/0,
                     catalog_root_id,
                     accession_root_id,
                     /*tombstone_root_id=*/0);
    mw.flush();

    // Copy metadata from local temp to NFS in 64 MB chunks
    {
        const uint64_t meta_size = mw.current_offset() - meta_base;
        const size_t   CHUNK     = 64 * 1024 * 1024;
        int local_fd = ::open(meta_tmp.data(), O_RDONLY);
        int nfs_fd   = ::open(out_path.c_str(), O_WRONLY | O_SYNC);
        if (local_fd < 0 || nfs_fd < 0) {
            if (local_fd >= 0) ::close(local_fd);
            if (nfs_fd   >= 0) ::close(nfs_fd);
            throw std::runtime_error("Cannot open files for metadata copy");
        }
        std::vector<uint8_t> buf(CHUNK);
        uint64_t copied = 0;
        while (copied < meta_size) {
            size_t  chunk = static_cast<size_t>(std::min<uint64_t>(CHUNK, meta_size - copied));
            off_t   soff  = static_cast<off_t>(meta_base + copied);
            ssize_t nr    = ::pread(local_fd, buf.data(), chunk, soff);
            if (nr <= 0) { ::close(local_fd); ::close(nfs_fd);
                throw std::runtime_error("Local metadata read failed"); }
            size_t wr = 0;
            while (wr < static_cast<size_t>(nr)) {
                ssize_t nw = ::pwrite(nfs_fd, buf.data() + wr,
                                      static_cast<size_t>(nr) - wr,
                                      soff + static_cast<off_t>(wr));
                if (nw <= 0) { ::close(local_fd); ::close(nfs_fd);
                    throw std::runtime_error("NFS metadata write failed"); }
                wr += static_cast<size_t>(nw);
            }
            copied += static_cast<uint64_t>(nr);
        }
        ::fsync(nfs_fd);
        ::close(local_fd);
        ::close(nfs_fd);
    }
    ::unlink(meta_tmp.data());

    spdlog::info("Repack complete: {} → {}", input_gpk.string(), out_path.string());
}

} // namespace genopack
