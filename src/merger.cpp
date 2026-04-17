#include <genopack/merger.hpp>
#include <genopack/string_collision_heap.hpp>
#include <genopack/accx.hpp>
#include <genopack/cidx.hpp>
#include <genopack/catalog.hpp>
#include <genopack/format.hpp>
#include <genopack/gidx.hpp>
#include <genopack/hnsw_section.hpp>
#include <genopack/kmrx.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/shard.hpp>
#include <genopack/skch.hpp>
#include <genopack/taxn.hpp>
#include <genopack/toc.hpp>
#include <genopack/tombstone.hpp>
#include <genopack/txdb.hpp>
#include <genopack/types.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <ctime>
#include <fcntl.h>
#include <filesystem>
#include <future>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <unistd.h>
#include <zstd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace genopack {

void merge_archives(const std::vector<std::filesystem::path>& inputs,
                    const std::filesystem::path& output,
                    bool remap_genome_ids,
                    bool merge_cidx)
{
    if (inputs.empty())
        throw std::runtime_error("merge_archives: no inputs");

    std::filesystem::path out_path = output;
    if (out_path.extension() != ".gpk")
        out_path = std::filesystem::path(out_path.string() + ".gpk");
    if (!out_path.parent_path().empty())
        std::filesystem::create_directories(out_path.parent_path());

    spdlog::info("Merging {} archives -> {}", inputs.size(), out_path.string());

    struct InputArchive {
        std::filesystem::path path;
        MmapFileReader        mmap;
        Toc                   toc;
    };

    std::vector<InputArchive> archives;
    archives.reserve(inputs.size());
    for (const auto& p : inputs) {
        InputArchive a;
        a.path = p;
        a.mmap.open(p);
        if (a.mmap.size() < sizeof(FileHeader))
            throw std::runtime_error("merge: too small: " + p.string());
        const auto* fhdr = a.mmap.ptr_at<FileHeader>(0);
        if (fhdr->magic != GPK2_MAGIC)
            throw std::runtime_error("merge: not a v2 .gpk file: " + p.string());
        // Validate TailLocator before reading the TOC.
        // A corrupt TailLocator (from NFS write-back failure) would cause
        // TocReader to read garbage SectionDescs with invalid compressed_size values,
        // leading to wildly wrong write offsets in the merged output.
        if (a.mmap.size() < sizeof(TailLocator))
            throw std::runtime_error("merge: file too small for TailLocator: " + p.string());
        const auto* tail = a.mmap.ptr_at<TailLocator>(a.mmap.size() - sizeof(TailLocator));
        if (tail->magic != GPKT_MAGIC)
            throw std::runtime_error("merge: corrupt TailLocator in " + p.string() +
                                     " — rebuild this archive before merging");
        a.toc = TocReader::read(a.mmap);
        spdlog::info("  Input: {} ({} sections)", p.string(), a.toc.sections.size());
        archives.push_back(std::move(a));
    }

    std::vector<GenomeId> gid_offsets(archives.size(), 0);
    std::vector<uint32_t> sid_offsets(archives.size(), 0);
    {
        GenomeId next_gid = 0;
        uint32_t next_sid = 0;
        for (size_t i = 0; i < archives.size(); ++i) {
            gid_offsets[i] = remap_genome_ids ? next_gid : 0;
            sid_offsets[i] = next_sid;

            GenomeId local_max_gid = 0;
            MergedCatalogReader merged;
            for (const auto* sd : archives[i].toc.find_by_type(SEC_CATL))
                merged.add_fragment(archives[i].mmap.data(), sd->file_offset, sd->compressed_size);
            merged.scan([&](const GenomeMeta& m) {
                local_max_gid = std::max(local_max_gid, m.genome_id);
                return true;
            });
            if (remap_genome_ids)
                next_gid = gid_offsets[i] + local_max_gid + 1;

            uint32_t local_max_sid = 0;
            for (const auto& sd : archives[i].toc.sections) {
                if (sd.type == SEC_SHRD) {
                    uint32_t shard_id = static_cast<uint32_t>(sd.aux0);
                    local_max_sid = std::max(local_max_sid, shard_id + 1);
                }
            }
            next_sid = sid_offsets[i] + local_max_sid;
        }
    }

    struct ShrdJob {
        size_t            archive_idx;
        const SectionDesc* sd;
        uint64_t          out_offset;
        SectionDesc       new_sd;
    };

    struct BulkCopyJob {
        size_t            archive_idx;
        const SectionDesc* sd;
        uint64_t          out_offset;
        SectionDesc       new_sd;
        GenomeId          gid_offset;
    };

    std::vector<ShrdJob>     shrd_jobs;
    std::vector<BulkCopyJob> hnsw_jobs;
    std::vector<BulkCopyJob> kmrx_jobs;
    uint64_t shrd_end_offset = sizeof(FileHeader);
    uint64_t next_section_id = 1;

    for (size_t i = 0; i < archives.size(); ++i) {
        uint32_t sid_off = sid_offsets[i];
        for (const auto* sd : archives[i].toc.find_by_type(SEC_SHRD)) {
            ShrdJob job;
            job.archive_idx         = i;
            job.sd                  = sd;
            job.out_offset          = shrd_end_offset;
            job.new_sd              = *sd;
            job.new_sd.section_id   = next_section_id++;
            job.new_sd.file_offset  = shrd_end_offset;
            job.new_sd.aux0         = static_cast<uint32_t>(sd->aux0) + sid_off;
            shrd_jobs.push_back(job);
            shrd_end_offset += sd->compressed_size;
        }
    }

    // Plan HNSW+KMRX bulk copies (placed immediately after SHRD sections)
    // The label_map / genome_ids arrays are patched in-place with the gid_offset.
    uint64_t hnsw_pos = shrd_end_offset;
    for (size_t i = 0; i < archives.size(); ++i) {
        const SectionDesc* best = nullptr;
        for (const auto* sd : archives[i].toc.find_by_type(SEC_HNSW))
            if (!best || sd->section_id > best->section_id) best = sd;
        if (!best) continue;
        BulkCopyJob job;
        job.archive_idx   = i;
        job.sd            = best;
        job.out_offset    = hnsw_pos;
        job.new_sd        = *best;
        job.new_sd.section_id  = next_section_id++;
        job.new_sd.file_offset = hnsw_pos;
        job.gid_offset    = gid_offsets[i];
        hnsw_jobs.push_back(job);
        hnsw_pos += best->compressed_size;
    }
    uint64_t kmrx_pos = hnsw_pos;
    for (size_t i = 0; i < archives.size(); ++i) {
        const SectionDesc* best = nullptr;
        for (const auto* sd : archives[i].toc.find_by_type(SEC_KMRX))
            if (!best || sd->section_id > best->section_id) best = sd;
        if (!best) continue;
        BulkCopyJob job;
        job.archive_idx   = i;
        job.sd            = best;
        job.out_offset    = kmrx_pos;
        job.new_sd        = *best;
        job.new_sd.section_id  = next_section_id++;
        job.new_sd.file_offset = kmrx_pos;
        job.gid_offset    = gid_offsets[i];
        kmrx_jobs.push_back(job);
        kmrx_pos += best->compressed_size;
    }
    uint64_t meta_start_offset = kmrx_pos;

    // Build per-archive map: old section_id → new section_id (for GIDX remapping).
    // Used by futures to read existing GIDX sections instead of re-scanning shards.
    std::vector<std::unordered_map<uint32_t, uint32_t>> old_to_new_secid(archives.size());
    for (const auto& job : shrd_jobs)
        old_to_new_secid[job.archive_idx][static_cast<uint32_t>(job.sd->section_id)] =
            static_cast<uint32_t>(job.new_sd.section_id);

    AppendWriter writer;
    writer.create(out_path);
    {
        FileHeader fhdr{};
        fhdr.magic           = GPK2_MAGIC;
        fhdr.version_major   = FORMAT_MAJOR;
        fhdr.version_minor   = FORMAT_MINOR;
        uint64_t t           = static_cast<uint64_t>(std::time(nullptr));
        fhdr.file_uuid_lo    = t ^ 0xdeadbeefcafe0001ULL;
        fhdr.file_uuid_hi    = (t << 17) ^ 0x1234567890abcdefULL;
        fhdr.created_at_unix = t;
        fhdr.flags           = 0;
        std::memset(fhdr.reserved, 0, sizeof(fhdr.reserved));
        writer.append(&fhdr, sizeof(fhdr));
    }

    struct ArchiveResult {
        std::vector<GenomeMeta>                              metas;
        std::vector<std::pair<std::string, GenomeId>>        accessions;  // sorted by accession
        std::vector<std::pair<std::string, std::string>>     taxonomies;  // sorted by accession
        std::vector<std::tuple<GenomeId, uint32_t, uint32_t>> gidx_entries;
        std::vector<GenomeId>                                tombstones;
    };

    std::vector<std::future<ArchiveResult>> futures;
    futures.reserve(archives.size());

    for (size_t i = 0; i < archives.size(); ++i) {
        futures.push_back(std::async(std::launch::async, [&, i]() -> ArchiveResult {
            const auto& a       = archives[i];
            GenomeId gid_off    = gid_offsets[i];
            uint32_t sid_off    = sid_offsets[i];
            int      in_fd      = a.mmap.fd();

            // Copy SHRD sections using pread (not memcpy-from-mmap).
            // pread returns EIO on NFS error instead of delivering SIGBUS.
            // POSIX_FADV_DONTNEED after each shard evicts NFS page cache,
            // keeping RSS bounded regardless of archive size.
            for (const auto& job : shrd_jobs) {
                if (job.archive_idx != i) continue;

                std::vector<uint8_t> shard_bytes(job.sd->compressed_size);
                size_t total = 0;
                while (total < job.sd->compressed_size) {
                    ssize_t r = ::pread(in_fd,
                                        shard_bytes.data() + total,
                                        job.sd->compressed_size - total,
                                        static_cast<off_t>(job.sd->file_offset + total));
                    if (r <= 0)
                        throw std::runtime_error("pread failed reading shard from archive "
                                                 + std::to_string(i));
                    total += static_cast<size_t>(r);
                }
#ifdef POSIX_FADV_DONTNEED
                ::posix_fadvise(in_fd,
                                static_cast<off_t>(job.sd->file_offset),
                                static_cast<off_t>(job.sd->compressed_size),
                                POSIX_FADV_DONTNEED);
#endif

                auto* hdr = reinterpret_cast<ShardHeader*>(shard_bytes.data());
                hdr->shard_id   += sid_off;
                hdr->cluster_id += sid_off;

                auto* dir = reinterpret_cast<GenomeDirEntry*>(shard_bytes.data() + hdr->genome_dir_offset);
                for (uint32_t di = 0; di < hdr->n_genomes; ++di) {
                    if (remap_genome_ids)
                        dir[di].genome_id += gid_off;
                }

                writer.write_at(job.out_offset, shard_bytes.data(), job.sd->compressed_size);
            }

            ArchiveResult res;

            MergedCatalogReader merged_cat;
            for (const auto* sd : a.toc.find_by_type(SEC_CATL))
                merged_cat.add_fragment(a.mmap.data(), sd->file_offset, sd->compressed_size);
            merged_cat.scan([&](const GenomeMeta& m) {
                GenomeMeta rm = m;
                rm.genome_id += gid_off;
                rm.shard_id  += sid_off;
                res.metas.push_back(rm);
                return true;
            });

            // Read GIDX entries from existing GIDX section (no shard re-scan).
            // gidx_entries stores (new_genome_id, new_section_id, dir_index).
            const auto& secid_map = old_to_new_secid[i];
            for (const auto* sd : a.toc.find_by_type(SEC_GIDX)) {
                const uint8_t* sec = a.mmap.data() + sd->file_offset;
                const auto* hdr = reinterpret_cast<const GidxHeader*>(sec);
                const auto* entries = reinterpret_cast<const GidxEntry*>(sec + hdr->entries_offset);
                for (uint32_t j = 0; j < hdr->n_entries; ++j) {
                    const auto& e = entries[j];
                    auto it = secid_map.find(e.shard_section_id);
                    if (it == secid_map.end()) continue;
                    res.gidx_entries.emplace_back(e.genome_id + gid_off, it->second, e.dir_index);
                }
            }

            for (const auto* sd : a.toc.find_by_type(SEC_ACCX)) {
                AccessionIndexReader air;
                air.open(a.mmap.data(), sd->file_offset, sd->compressed_size);
                // scan() returns entries in storage order (sorted by accession).
                air.scan([&](std::string_view acc, GenomeId gid) {
                    res.accessions.emplace_back(std::string(acc), gid + gid_off);
                });
            }

            for (const auto* sd : a.toc.find_by_type(SEC_TAXN)) {
                TaxonomyIndexReader tir;
                tir.open(a.mmap.data(), sd->file_offset, sd->compressed_size);
                // scan() returns entries in storage order; we re-sort by accession
                // after collecting so the per-archive list is sorted by accession key
                // (required for the merge-heap invariant).
                size_t tax_start = res.taxonomies.size();
                tir.scan([&](std::string_view acc, std::string_view tax) {
                    res.taxonomies.emplace_back(std::string(acc), std::string(tax));
                });
                std::sort(res.taxonomies.begin() + static_cast<ptrdiff_t>(tax_start),
                          res.taxonomies.end(),
                          [](const std::pair<std::string,std::string>& a,
                             const std::pair<std::string,std::string>& b) {
                              return a.first < b.first;
                          });
            }

            for (const auto* sd : a.toc.find_by_type(SEC_TOMB)) {
                TombstoneReader tr;
                tr.open(a.mmap.data(), sd->file_offset, sd->compressed_size);
                tr.scan([&](GenomeId id) {
                    res.tombstones.push_back(id + gid_off);
                });
            }

            spdlog::info("  Finished archive {}: {} genomes, {} accessions, {} tombstones",
                         i, res.metas.size(), res.accessions.size(), res.tombstones.size());
            return res;
        }));
    }

    std::vector<GenomeMeta> all_metas;
    std::vector<std::pair<std::string, GenomeId>>      all_accessions;
    std::vector<std::pair<std::string, std::string>>   all_taxonomies;
    std::vector<std::tuple<GenomeId, uint32_t, uint32_t>> all_gidx_entries;
    std::unordered_set<GenomeId> all_tombstones;

    std::vector<std::vector<std::pair<std::string, GenomeId>>>    acc_lists;
    std::vector<std::vector<std::pair<std::string, std::string>>> tax_lists;
    acc_lists.reserve(futures.size());
    tax_lists.reserve(futures.size());

    for (auto& f : futures) {
        ArchiveResult res = f.get();
        for (auto& m : res.metas) all_metas.push_back(std::move(m));
        for (auto& ge : res.gidx_entries) all_gidx_entries.push_back(std::move(ge));
        for (GenomeId gid : res.tombstones) all_tombstones.insert(gid);
        acc_lists.push_back(std::move(res.accessions));
        tax_lists.push_back(std::move(res.taxonomies));
    }

    // Merge T sorted lists using Myers string heap: O(M log T + S).
    string_heap_merge(acc_lists, [&](const std::string& k, GenomeId v) {
        all_accessions.emplace_back(k, v);
    });
    string_heap_merge(tax_lists, [&](const std::string& k, const std::string& v) {
        all_taxonomies.emplace_back(k, v);
    });

    std::sort(all_metas.begin(), all_metas.end(), [](const GenomeMeta& a, const GenomeMeta& b) {
        if (a.oph_fingerprint != b.oph_fingerprint)
            return a.oph_fingerprint < b.oph_fingerprint;
        return a.genome_id < b.genome_id;
    });

    spdlog::info("Merge: {} shards, {} genomes collected", shrd_jobs.size(), all_metas.size());

    // Copy HNSW sections — patch label_map genome IDs in-place, no rebuild.
    for (const auto& job : hnsw_jobs) {
        const auto& a = archives[job.archive_idx];
        std::vector<uint8_t> buf(job.sd->compressed_size);
        std::memcpy(buf.data(), a.mmap.data() + job.sd->file_offset, job.sd->compressed_size);
        if (job.gid_offset > 0) {
            auto* hdr = reinterpret_cast<HnswSectionHeader*>(buf.data());
            auto* lm  = reinterpret_cast<uint64_t*>(buf.data() + hdr->label_map_offset);
            for (uint32_t j = 0; j < hdr->n_elements; ++j)
                lm[j] += job.gid_offset;
        }
        writer.write_at(job.out_offset, buf.data(), buf.size());
        spdlog::info("  Copied HNSW from archive {}: {} elements", job.archive_idx, job.sd->item_count);
    }

    // Copy KMRX sections — patch sorted genome_ids array in-place, no rebuild.
    for (const auto& job : kmrx_jobs) {
        const auto& a = archives[job.archive_idx];
        std::vector<uint8_t> buf(job.sd->compressed_size);
        std::memcpy(buf.data(), a.mmap.data() + job.sd->file_offset, job.sd->compressed_size);
        if (job.gid_offset > 0) {
            auto* hdr = reinterpret_cast<KmrxHeader*>(buf.data());
            auto* ids = reinterpret_cast<uint64_t*>(buf.data() + hdr->index_offset);
            for (uint32_t j = 0; j < hdr->n_genomes; ++j)
                ids[j] += job.gid_offset;
        }
        writer.write_at(job.out_offset, buf.data(), buf.size());
        spdlog::info("  Copied KMRX from archive {}: {} genomes", job.archive_idx, job.sd->item_count);
    }

    TocWriter toc_out;
    for (const auto& job : shrd_jobs)
        toc_out.add_section(job.new_sd);
    for (const auto& job : hnsw_jobs)
        toc_out.add_section(job.new_sd);
    for (const auto& job : kmrx_jobs)
        toc_out.add_section(job.new_sd);

    // Flush all async shard/HNSW/KMRX writes to NFS before writing metadata.
    writer.flush();

    // Write all metadata to a LOCAL temp file first.
    // NFS write-back + ENOSPC can silently lose metadata regardless of O_DSYNC.
    // Writing locally guarantees integrity; we then copy to NFS in small O_SYNC chunks.
    auto meta_tmp_fs = out_path.parent_path() / "gpk_meta_XXXXXX";
    std::string meta_tmp_s = meta_tmp_fs.string();
    std::vector<char> meta_tmp_path(meta_tmp_s.begin(), meta_tmp_s.end());
    meta_tmp_path.push_back('\0');
    {
        int fd = ::mkstemp(meta_tmp_path.data());
        if (fd < 0)
            throw std::runtime_error("Cannot create temp metadata file: " +
                                     std::string(std::strerror(errno)));
        ::close(fd);
    }
    AppendWriter mw;   // local metadata writer
    mw.create(meta_tmp_path.data());
    mw.seek_to(meta_start_offset);  // so all section.file_offset values are NFS-relative

    uint64_t catalog_root_id = 0;
    {
        CatalogSectionWriter csw;
        for (const auto& m : all_metas)
            csw.add(m);
        SectionDesc catl_sd = csw.finalize(mw, next_section_id++);
        catalog_root_id = catl_sd.section_id;
        toc_out.add_section(catl_sd);
    }

    {
        std::unordered_map<GenomeId, uint64_t> gid_to_catl_row;
        gid_to_catl_row.reserve(all_metas.size());
        for (uint64_t i = 0; i < all_metas.size(); ++i)
            gid_to_catl_row[all_metas[i].genome_id] = i;

        GidxWriter gw;
        for (const auto& [gid, new_section_id, dir_idx] : all_gidx_entries) {
            auto catl_it = gid_to_catl_row.find(gid);
            if (catl_it == gid_to_catl_row.end()) continue;
            gw.add(gid, new_section_id, dir_idx, catl_it->second);
        }
        SectionDesc gidx_sd = gw.finalize(mw, next_section_id++);
        toc_out.add_section(gidx_sd);
    }

    uint64_t accession_root_id = 0;
    {
        AccessionIndexWriter aiw;
        for (const auto& [acc, gid] : all_accessions)
            aiw.add(acc, gid);
        SectionDesc accx_sd = aiw.finalize_presorted(mw, next_section_id++);
        accession_root_id = accx_sd.section_id;
        toc_out.add_section(accx_sd);
    }

    {
        // Merge CIDX sections per-archive, applying each archive's gid_offset directly.
        // Using a flat genome_id → offset map is wrong when archives reuse IDs (1,2,3…);
        // iterating per-archive guarantees the right offset for every entry.
        {
            CidxWriter cw;
            uint64_t total_entries = 0;
            for (size_t ai = 0; ai < archives.size(); ++ai) {
                GenomeId gid_off = gid_offsets[ai];
                for (const auto* sd : archives[ai].toc.find_by_type(SEC_CIDX)) {
                    CidxReader cr;
                    cr.open(archives[ai].mmap.data(), sd->file_offset, sd->compressed_size);
                    cr.scan([&](uint64_t h, uint32_t old_gid) {
                        cw.add_hash(h, static_cast<uint32_t>(old_gid + gid_off));
                    });
                    total_entries += cr.n_entries();
                }
            }
            if (total_entries > 0 && merge_cidx) {
                spdlog::info("Merging CIDX: {} total contig accessions", total_entries);
                SectionDesc cidx_sd = cw.finalize(mw, next_section_id++, /*batch_id=*/0);
                toc_out.add_section(cidx_sd);
            }
        }
    }

    if (!all_taxonomies.empty()) {
        TaxonomyIndexWriter tiw;
        for (const auto& [acc, tax] : all_taxonomies)
            tiw.add(acc, tax);
        SectionDesc taxn_sd = tiw.finalize_presorted(mw, next_section_id++);
        toc_out.add_section(taxn_sd);

        TxdbWriter txw;
        for (const auto& [acc, tax] : all_taxonomies)
            txw.add(acc, tax);
        SectionDesc txdb_sd = txw.finalize(mw, next_section_id++);
        toc_out.add_section(txdb_sd);
    }

    uint64_t tombstone_root_id = 0;
    if (!all_tombstones.empty()) {
        TombstoneWriter tw;
        for (GenomeId gid : all_tombstones)
            tw.add(gid);
        SectionDesc tomb_sd = tw.finalize(mw, next_section_id++);
        tombstone_root_id = tomb_sd.section_id;
        toc_out.add_section(tomb_sd);
    }

    // SKCH: merge OPH sketch sections from all archives.
    // V4-only: every source archive must have SKC4 sections. Each sketch
    // is read via SkchReader::sketch_for_ids and re-emitted carrying both
    // sig1 and sig2. If ANY archive has a SKCH section that isn't V4 we
    // abort — there is no sig2 to transit.
    {
        uint32_t skch_sketch_size = 0, skch_kmer_size = 0, skch_syncmer_s = 0;
        uint64_t skch_seed1 = 42, skch_seed2 = 43;
        bool     skch_params_set = false;
        size_t   total_sketches  = 0;

        for (size_t ai = 0; ai < archives.size(); ++ai) {
            for (const auto* sd : archives[ai].toc.find_by_type(SEC_SKCH)) {
                // peek_params throws if the section is not SKC4.
                auto [ver, ks] = SkchReader::peek_params(
                    archives[ai].mmap.data(), sd->file_offset, sd->compressed_size);
                if (ver != 4) {
                    throw std::runtime_error(
                        "SKC4 required; source archive lacks sig2 — rebuild first");
                }
                SkchReader reader;
                reader.open(archives[ai].mmap.data(), sd->file_offset, sd->compressed_size);
                if (!skch_params_set) {
                    skch_sketch_size = reader.sketch_size();
                    skch_kmer_size   = reader.kmer_size();
                    skch_syncmer_s   = reader.syncmer_s();
                    skch_seed1       = reader.seed1();
                    skch_seed2       = reader.seed2();
                    skch_params_set  = true;
                } else {
                    if (reader.sketch_size() != skch_sketch_size ||
                        reader.kmer_size()   != skch_kmer_size   ||
                        reader.seed1()       != skch_seed1       ||
                        reader.seed2()       != skch_seed2) {
                        throw std::runtime_error(
                            "SKCH merge: mismatched parameters across source archives");
                    }
                }
                total_sketches += reader.n_genomes();
            }
        }

        if (total_sketches > 0) {
            spdlog::info("Merging SKCH (V4): {} sketches from {} archives",
                         total_sketches, archives.size());
            SkchWriter skch_out(skch_sketch_size, skch_kmer_size,
                                skch_syncmer_s, skch_seed1, skch_seed2);

            const uint32_t mask_words = (skch_sketch_size + 63) / 64;

            for (size_t ai = 0; ai < archives.size(); ++ai) {
                GenomeId gid_off = gid_offsets[ai];
                for (const auto* sd : archives[ai].toc.find_by_type(SEC_SKCH)) {
                    SkchReader reader;
                    reader.open(archives[ai].mmap.data(),
                                sd->file_offset, sd->compressed_size);

                    const auto& src_ids = reader.genome_ids(); // sorted ascending

                    // Pre-allocate per-row remapped output since cb may fire
                    // concurrently from multiple threads — serialise through a
                    // mutex-protected SkchWriter::add.
                    std::vector<GenomeId>              out_ids(src_ids.size());
                    std::vector<std::vector<uint16_t>> out_sig1(src_ids.size());
                    std::vector<std::vector<uint16_t>> out_sig2(src_ids.size());
                    std::vector<std::vector<uint64_t>> out_mask(src_ids.size());
                    std::vector<uint32_t>              out_nrb(src_ids.size());
                    std::vector<uint64_t>              out_glen(src_ids.size());

                    reader.sketch_for_ids(src_ids, skch_kmer_size, skch_sketch_size,
                        [&](size_t row, const SketchResult& sr) {
                            out_ids[row]  = static_cast<GenomeId>(src_ids[row]) + gid_off;
                            out_sig1[row].assign(sr.sig,  sr.sig  + sr.sketch_size);
                            out_sig2[row].assign(sr.sig2, sr.sig2 + sr.sketch_size);
                            out_mask[row].assign(sr.mask, sr.mask + sr.mask_words);
                            out_nrb[row]  = sr.n_real_bins;
                            out_glen[row] = sr.genome_length;
                        });

                    for (size_t j = 0; j < src_ids.size(); ++j) {
                        if (out_sig1[j].empty()) continue;
                        skch_out.add(out_ids[j],
                                     out_sig1[j], out_sig2[j],
                                     out_nrb[j], out_glen[j],
                                     out_mask[j]);
                    }
                    (void)mask_words;
                }
            }

            SectionDesc skch_sd = skch_out.finalize(mw, next_section_id++);
            toc_out.add_section(skch_sd);
        }
    }

    uint64_t live_count = 0;
    for (const auto& m : all_metas) {
        if (!m.is_deleted() && !all_tombstones.contains(m.genome_id))
            ++live_count;
    }

    toc_out.finalize(mw,
                     /*generation=*/1,
                     live_count,
                     static_cast<uint64_t>(all_metas.size()),
                     /*prev_toc_offset=*/0,
                     catalog_root_id,
                     accession_root_id,
                     tombstone_root_id);
    mw.flush();

    // Validate TailLocator in the local file before copying to NFS.
    {
        MmapFileReader local_mm;
        local_mm.open(meta_tmp_path.data());
        if (local_mm.size() < meta_start_offset + sizeof(TailLocator))
            throw std::runtime_error("Local metadata file too small — internal error");
        const uint64_t tail_off = mw.current_offset() - sizeof(TailLocator);
        const auto* tail = local_mm.ptr_at<TailLocator>(tail_off);
        if (tail->magic != GPKT_MAGIC)
            throw std::runtime_error("Local metadata TailLocator invalid — internal error");
        spdlog::info("Local metadata validated ({} bytes), copying to NFS...",
                     mw.current_offset() - meta_start_offset);
    }

    // Copy metadata from local file to NFS in 64 MB O_SYNC chunks.
    // O_SYNC at open() time is respected by NFS; small chunks bound retry granularity.
    {
        const uint64_t meta_size   = mw.current_offset() - meta_start_offset;
        const size_t   CHUNK       = 64 * 1024 * 1024;
        int local_fd = ::open(meta_tmp_path.data(), O_RDONLY);
        if (local_fd < 0)
            throw std::runtime_error("Cannot re-open local metadata file");
        int nfs_fd = ::open(out_path.c_str(), O_WRONLY | O_SYNC);
        if (nfs_fd < 0) {
            ::close(local_fd);
            throw std::runtime_error("Cannot open NFS output for O_SYNC write");
        }

        std::vector<uint8_t> buf(CHUNK);
        uint64_t copied = 0;
        while (copied < meta_size) {
            size_t  chunk     = static_cast<size_t>(std::min<uint64_t>(CHUNK, meta_size - copied));
            off_t   src_off   = static_cast<off_t>(meta_start_offset + copied);
            ssize_t nr        = ::pread(local_fd, buf.data(), chunk, src_off);
            if (nr <= 0) {
                ::close(local_fd); ::close(nfs_fd);
                throw std::runtime_error("Local metadata read failed at offset " +
                                         std::to_string(src_off));
            }
            size_t  written = 0;
            int     retries = 0;
            while (written < static_cast<size_t>(nr)) {
                ssize_t nw = ::pwrite(nfs_fd, buf.data() + written,
                                      static_cast<size_t>(nr) - written,
                                      src_off + static_cast<off_t>(written));
                if (nw < 0) {
                    if ((errno == ENOSPC || errno == EIO) && retries < 60) {
                        spdlog::warn("NFS metadata copy ENOSPC at offset {}, retry {}...",
                                     src_off + written, ++retries);
                        std::this_thread::sleep_for(std::chrono::seconds(5));
                        continue;
                    }
                    ::close(local_fd); ::close(nfs_fd);
                    throw std::runtime_error("NFS metadata write failed: " +
                                             std::string(std::strerror(errno)));
                }
                written += static_cast<size_t>(nw);
            }
            copied += static_cast<uint64_t>(nr);
        }
        if (::fsync(nfs_fd) != 0)
            spdlog::warn("Final fsync of metadata failed ({}), continuing anyway",
                         std::strerror(errno));
        ::close(local_fd);
        ::close(nfs_fd);
    }

    ::unlink(meta_tmp_path.data());
    spdlog::info("Merge complete: {} inputs -> {}", inputs.size(), out_path.string());
}

} // namespace genopack
