#include <genopack/subset.hpp>
#include <genopack/accx.hpp>
#include <genopack/catalog.hpp>
#include <genopack/cidx.hpp>
#include <genopack/format.hpp>
#include <genopack/gidx.hpp>
#include <genopack/kmrx.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/qual.hpp>
#include <genopack/qual_columns.hpp>
#include <genopack/section_checksum.hpp>
#include <genopack/shard.hpp>
#include <genopack/taxn.hpp>
#include <genopack/toc.hpp>
#include <genopack/txdb.hpp>
#include <spdlog/spdlog.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <cctype>
#include <condition_variable>
#include <cstring>
#include <ctime>
#include <fcntl.h>
#include <sys/sendfile.h>
#include <filesystem>
#include <fstream>
#include <future>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <string>
#include <thread>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace genopack {

// ── Accession file loader ────────────────────────────────────────────────────

std::unordered_set<std::string> load_accession_set(const std::filesystem::path& path)
{
    std::ifstream f(path);
    if (!f) throw std::runtime_error("subset: cannot open accessions file: " + path.string());

    std::unordered_set<std::string> result;
    std::string line;

    // Peek at first line to detect TSV vs plain-text
    if (!std::getline(f, line)) return result;

    if (line.find('\t') != std::string::npos) {
        // TSV: find column index of "accession" or "name" (case-insensitive)
        auto lower = [](std::string s) {
            for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            return s;
        };

        std::vector<std::string> headers;
        {
            std::istringstream ss(line);
            std::string col;
            while (std::getline(ss, col, '\t')) headers.push_back(col);
        }

        int col_idx = -1;
        for (int i = 0; i < static_cast<int>(headers.size()); ++i) {
            std::string h = lower(headers[i]);
            if (h == "accession" || h == "name") { col_idx = i; break; }
        }
        if (col_idx < 0)
            throw std::runtime_error(
                "subset: TSV header has no 'accession' or 'name' column in: " + path.string());

        while (std::getline(f, line)) {
            if (line.empty()) continue;
            std::istringstream ss(line);
            std::string field;
            for (int i = 0; i <= col_idx; ++i) {
                if (!std::getline(ss, field, '\t')) { field.clear(); break; }
            }
            if (!field.empty()) result.insert(field);
        }
    } else {
        // Plain-text: first line is already an accession
        if (!line.empty()) result.insert(line);
        while (std::getline(f, line))
            if (!line.empty()) result.insert(line);
    }

    return result;
}

// ── subset_archive ────────────────────────────────────────────────────────────

void subset_archive(const std::filesystem::path& input_gpk,
                    const std::filesystem::path& output_gpk,
                    const SubsetConfig& cfg)
{
    // ── Open source archive ──────────────────────────────────────────────────
    MmapFileReader src_mmap;
    src_mmap.open(input_gpk);
    if (src_mmap.size() < sizeof(FileHeader))
        throw std::runtime_error("subset: file too small: " + input_gpk.string());
    if (src_mmap.ptr_at<FileHeader>(0)->magic != GPK_MAGIC)
        throw std::runtime_error("subset: not a v2 .gpk archive: " + input_gpk.string());

    Toc src_toc = TocReader::read(src_mmap);
    auto src_shards = src_toc.find_by_type(SEC_SHRD);

    // Sort shards by file offset for sequential NFS access
    std::sort(src_shards.begin(), src_shards.end(), [](const SectionDesc* a, const SectionDesc* b) {
        return a->file_offset < b->file_offset;
    });

    spdlog::info("Subset: source {} ({} shards, {} target accessions)",
                 input_gpk.string(), src_shards.size(), cfg.accessions.size());

    // ── Load catalog, ACCX, TAXN ─────────────────────────────────────────────
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
    spdlog::info("Subset: {} live genomes in source", gid_to_meta.size());

    std::unordered_map<GenomeId, std::string> gid_to_acc;
    gid_to_acc.reserve(gid_to_meta.size());
    for (auto* sd : src_toc.find_by_type(SEC_ACCX)) {
        AccessionIndexReader air;
        air.open(src_mmap.data(), sd->file_offset, sd->compressed_size);
        air.scan([&](std::string_view acc, GenomeId gid) {
            if (gid_to_meta.count(gid)) gid_to_acc.emplace(gid, std::string(acc));
        });
    }

    // Build kept genome set: gids whose accession is in cfg.accessions
    std::unordered_set<GenomeId> kept_gids;
    kept_gids.reserve(cfg.accessions.size());
    for (const auto& [gid, acc] : gid_to_acc) {
        if (cfg.accessions.count(acc)) kept_gids.insert(gid);
    }
    spdlog::info("Subset: {} genomes match target set", kept_gids.size());

    // Filtered gid_to_meta and gid_to_acc restricted to kept set
    {
        std::unordered_map<GenomeId, GenomeMeta> kept_meta;
        kept_meta.reserve(kept_gids.size());
        for (const auto& gid : kept_gids) {
            auto it = gid_to_meta.find(gid);
            if (it != gid_to_meta.end()) kept_meta.emplace(gid, it->second);
        }
        gid_to_meta = std::move(kept_meta);
    }

    std::unordered_map<std::string, std::string> acc_to_tax;
    for (auto* sd : src_toc.find_by_type(SEC_TAXN)) {
        TaxonomyIndexReader tir;
        tir.open(src_mmap.data(), sd->file_offset, sd->compressed_size);
        tir.scan([&](std::string_view acc, std::string_view tax) {
            if (cfg.accessions.count(std::string(acc)))
                acc_to_tax.emplace(std::string(acc), std::string(tax));
        });
    }
    spdlog::info("Subset: {} taxonomy entries for kept genomes", acc_to_tax.size());

    // ── Phase 1: directory scan — skip genomes not in kept_gids ─────────────
    struct GenomeRecord {
        GenomeId  genome_id;
        uint64_t  oph;
        uint32_t  src_shard_idx;
        uint32_t  dir_idx;
        GenomeMeta meta;
        uint32_t  _pad0;
        uint32_t  _pad1;
    };

    std::vector<GenomeRecord> records;
    records.reserve(kept_gids.size());

    for (size_t s = 0; s < src_shards.size(); ++s) {
        const SectionDesc* sd = src_shards[s];
        ShardReader shard;
        shard.open(src_mmap.data(), sd->file_offset, sd->compressed_size);

        uint32_t dir_idx = 0;
        for (auto* de = shard.dir_begin(); de != shard.dir_end(); ++de, ++dir_idx) {
            GenomeId gid = static_cast<GenomeId>(de->genome_id);
            if (!kept_gids.count(gid)) continue;

            auto meta_it = gid_to_meta.find(gid);
            if (meta_it == gid_to_meta.end()) continue;

            records.push_back({gid, de->oph_fingerprint,
                               static_cast<uint32_t>(s), dir_idx,
                               meta_it->second, 0u, 0u});
        }

        if ((s + 1) % 5000 == 0)
            spdlog::info("Phase 1: {}/{} shards scanned, {} genomes indexed",
                         s + 1, src_shards.size(), records.size());
    }
    spdlog::info("Phase 1 complete: {} genomes", records.size());

    // ── Phase 2: sort by oph (no taxonomy regrouping for subset) ────────────
    std::sort(records.begin(), records.end(), [](const GenomeRecord& a, const GenomeRecord& b) {
        return a.oph < b.oph;
    });

    // Build reverse index: src_shard_idx → sorted record positions
    std::vector<std::vector<uint32_t>> src_shard_records(src_shards.size());
    for (uint32_t i = 0; i < static_cast<uint32_t>(records.size()); ++i)
        src_shard_records[records[i].src_shard_idx].push_back(i);

    spdlog::info("Phase 2 complete: sorted by oph, ready for Phase 3");

    // ── Open output file ─────────────────────────────────────────────────────
    std::filesystem::path out_path = output_gpk;
    if (out_path.extension() != ".gpk")
        out_path = std::filesystem::path(out_path.string() + ".gpk");
    std::filesystem::create_directories(out_path.parent_path());

    AppendWriter app_writer;
    app_writer.create(out_path);

    {
        FileHeader ofh{};
        ofh.magic           = GPK_MAGIC;
        ofh.version_major   = FORMAT_MAJOR;
        ofh.endian_abi_tag  = ENDIAN_ABI_TAG;
        ofh.version_minor   = FORMAT_MINOR;
        uint64_t t = static_cast<uint64_t>(std::time(nullptr));
        ofh.file_uuid_lo    = t ^ 0xdeadbeefcafe0003ULL;
        ofh.file_uuid_hi    = (t << 17) ^ 0xfedcba9876543210ULL;
        ofh.created_at_unix = t;
        ofh.flags           = 0;
        std::memset(ofh.reserved, 0, sizeof(ofh.reserved));
        app_writer.append(&ofh, sizeof(ofh));
    }

    // ── Background IO writer thread ──────────────────────────────────────────
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
            uint64_t shard_bytes;
            if (!frozen.tmp_path.empty()) {
                FrozenShardReader rdr(frozen); // RAII: closes + unlinks on scope exit
                app_writer.append_from_fd(rdr.fd, frozen.section_size);
                shard_bytes = frozen.section_size;
            } else {
                app_writer.append(frozen.bytes.data(), frozen.bytes.size());
                shard_bytes = static_cast<uint64_t>(frozen.bytes.size());
            }

            SectionDesc sd{};
            sd.type              = SEC_SHRD;
            sd.version           = 4;
            sd.flags             = 0;
            sd.section_id        = wt.section_id;
            sd.file_offset       = shard_start;
            sd.compressed_size   = shard_bytes;
            sd.uncompressed_size = 0;
            sd.item_count        = frozen.n_genomes;
            sd.aux0              = wt.shard_id;
            sd.aux1              = 0;
            if (!frozen.tmp_path.empty()) {
                if (!checksum_of_fd(app_writer.fd(), shard_start, shard_bytes, sd.checksum))
                    throw std::runtime_error("subset drain: checksum_of_fd failed");
            } else {
                checksum_of(frozen.bytes.data(), frozen.bytes.size(), sd.checksum);
            }
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

    // Single shard writer (no taxonomy bucketing for subset — just one sequential shard stream)
    std::unique_ptr<ShardWriter> cur_writer;
    ShardId cur_shard_id = 0;
    uint64_t cur_bytes   = 0;
    uint64_t total_buffered_bytes = 0;

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
    size_t n_shards_done  = 0;
    size_t n_genomes_done = 0;

    // ── Phase 3: decompress + write ──────────────────────────────────────────
    for (size_t s = 0; s < src_shards.size(); ++s) {
        const auto& rec_idxs = src_shard_records[s];
        if (rec_idxs.empty()) { ++n_shards_done; continue; }

        prefetch_shard(s + PREFETCH_AHEAD);

        ShardReader shard;
        shard.open(src_mmap.data(), src_shards[s]->file_offset, src_shards[s]->compressed_size);

        // Raw-blob fast path: when no partial output shard is open, every genome
        // in this source shard is kept, and there are no deleted entries, copy the
        // compressed bytes directly without decompress→recompress.
        {
            const ShardHeader* src_hdr =
                reinterpret_cast<const ShardHeader*>(
                    src_mmap.data() + src_shards[s]->file_offset);
            if (!cur_writer &&
                src_hdr->n_deleted == 0 &&
                rec_idxs.size() == static_cast<size_t>(shard.n_genomes()))
            {
                ShardId raw_sid = current_shard_id++;
                uint64_t sec_id = next_section_id++;
                shard_id_to_section_id[raw_sid] = sec_id;

                for (size_t j = 0; j < rec_idxs.size(); ++j) {
                    const GenomeRecord& rec = records[rec_idxs[j]];
                    uint64_t catl_idx = new_catalog.size();
                    GenomeMeta m = rec.meta;
                    m.shard_id = raw_sid;
                    new_catalog.push_back(m);
                    new_gidx.push_back({rec.genome_id, raw_sid, rec.dir_idx, catl_idx});
                    ++n_genomes_done;
                }

                const uint8_t* blob_ptr = src_mmap.data() + src_shards[s]->file_offset;
                size_t blob_len = static_cast<size_t>(src_shards[s]->compressed_size);
                FrozenShard fs;
                fs.bytes.assign(blob_ptr, blob_ptr + blob_len);
                fs.shard_id  = raw_sid;
                fs.n_genomes = src_hdr->n_genomes;
                fs.raw_bytes = src_hdr->shard_raw_bp;
                WriteTask wt;
                wt.section_id = sec_id;
                wt.shard_id   = raw_sid;
                wt.fut = std::async(std::launch::deferred,
                    [fs = std::move(fs)]() mutable { return std::move(fs); });
                {
                    std::unique_lock lk(write_q_mx);
                    write_q_space_cv.wait(lk, [&]{ return write_q.size() < 2; });
                    write_q.push(std::move(wt));
                }
                write_q_cv.notify_one();
                shard.release_pages();
                ++n_shards_done;
                if (cfg.verbose || n_shards_done % 2000 == 0)
                    spdlog::info("Phase 3: {}/{} shards (raw), {} genomes",
                                 n_shards_done, src_shards.size(), n_genomes_done);
                continue;
            }
        }

        std::vector<std::string> fastas(rec_idxs.size());
        #pragma omp parallel for schedule(dynamic, 4) num_threads(n_threads)
        for (int j = 0; j < static_cast<int>(rec_idxs.size()); ++j)
            fastas[j] = shard.fetch_genome_at(records[rec_idxs[j]].dir_idx);

        for (size_t j = 0; j < rec_idxs.size(); ++j) {
            const GenomeRecord& rec = records[rec_idxs[j]];

            if (!cur_writer) {
                cur_shard_id = current_shard_id++;
                cur_writer   = std::make_unique<ShardWriter>(cur_shard_id, cur_shard_id, cfg.shard_cfg);
            }

            cur_writer->add_genome(rec.genome_id, rec.oph, fastas[j].data(), fastas[j].size());
            uint64_t fsize        = fastas[j].size();
            cur_bytes            += fsize;
            total_buffered_bytes += fsize;

            uint64_t catl_idx = new_catalog.size();
            GenomeMeta m = rec.meta;
            m.shard_id = cur_shard_id;
            new_catalog.push_back(m);
            new_gidx.push_back({rec.genome_id, cur_shard_id,
                                 static_cast<uint32_t>(cur_writer->n_genomes() - 1),
                                 catl_idx});
            ++n_genomes_done;

            // Natural flush: shard is full
            if (cur_bytes >= cfg.shard_cfg.max_shard_size_bytes) {
                total_buffered_bytes -= cur_bytes;
                launch_freeze(cur_writer, cur_shard_id);
                cur_bytes = 0;
            }
        }

        shard.release_pages();
        ++n_shards_done;

        // Evict if over cap
        if (total_buffered_bytes >= cfg.max_bucket_bytes && cur_writer) {
            spdlog::info("Subset: cap ({:.1f} GB), flushing current shard {}",
                         total_buffered_bytes / double(1ULL << 30), cur_shard_id);
            total_buffered_bytes -= cur_bytes;
            launch_freeze(cur_writer, cur_shard_id);
            cur_bytes = 0;
        }

        if (cfg.verbose || n_shards_done % 2000 == 0)
            spdlog::info("Phase 3: {}/{} shards, {} genomes, {:.1f} GB buffered",
                         n_shards_done, src_shards.size(), n_genomes_done,
                         total_buffered_bytes / double(1ULL << 30));
    }

    // Flush remaining writer
    if (cur_writer) {
        launch_freeze(cur_writer, cur_shard_id);
        cur_bytes = 0;
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

    // ── Write metadata sections ──────────────────────────────────────────────
    const uint64_t meta_base = app_writer.current_offset();
    auto meta_tmp_fs = std::filesystem::path(out_path).parent_path() / "gpk_subset_XXXXXX";
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

    // ACCX — only kept genomes
    {
        AccessionIndexWriter aiw;
        for (const auto& [gid, acc] : gid_to_acc) {
            if (kept_gids.count(gid)) aiw.add(acc, gid);
        }
        SectionDesc sd = aiw.finalize(mw, next_section_id++);
        accession_root_id = sd.section_id;
        new_toc.add_section(sd);
    }

    // TAXN — only kept accessions
    if (!acc_to_tax.empty()) {
        TaxonomyIndexWriter tiw;
        for (const auto& [acc, tax] : acc_to_tax) tiw.add(acc, tax);
        new_toc.add_section(tiw.finalize(mw, next_section_id++));
    }

    // KMRX — only kept genomes
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
            for (const GenomeId gid : kept_gids) {
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

    // Quality — only kept genomes; re-emitted as QCOL (reads QCOL or legacy QUAL).
    {
        std::vector<QualRecord> kept;
        for (auto* sd : src_toc.find_by_type(SEC_QCOL)) {
            ColStoreReader cr;
            cr.open(src_mmap.data(), sd->file_offset, sd->compressed_size);
            qcol_scan(cr, [&](const QualRecord& rec) {
                if (kept_gids.count(static_cast<GenomeId>(rec.genome_id))) kept.push_back(rec);
            });
        }
        for (auto* sd : src_toc.find_by_type(SEC_QUAL)) {
            QualReader qr;
            qr.open(src_mmap.data(), sd->file_offset, sd->compressed_size);
            qr.scan([&](const QualRecord& rec) {
                if (kept_gids.count(static_cast<GenomeId>(rec.genome_id))) kept.push_back(rec);
            });
        }
        if (!kept.empty())
            new_toc.add_section(qcol_write(mw, next_section_id++, std::move(kept)));
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

    // Content-checksum the metadata sections (SHRD hashed inline) before finalize.
    mw.flush();
    stamp_section_checksums(meta_tmp.data(), new_toc.sections());

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

    // Copy metadata from local temp to output in 64 MB chunks
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

    spdlog::info("Subset complete: {} → {} ({} genomes kept of {})",
                 input_gpk.string(), out_path.string(),
                 n_genomes_done, gid_to_acc.size());
}

} // namespace genopack
