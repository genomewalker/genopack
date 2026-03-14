#include <genopack/merger.hpp>
#include <genopack/accx.hpp>
#include <genopack/catalog.hpp>
#include <genopack/format_v2.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/toc.hpp>
#include <genopack/types.hpp>
#include <spdlog/spdlog.h>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <stdexcept>
#include <vector>

namespace genopack {

void merge_archives(const std::vector<std::filesystem::path>& inputs,
                    const std::filesystem::path& output)
{
    if (inputs.empty())
        throw std::runtime_error("merge_archives: no inputs");

    // Determine output path
    std::filesystem::path out_path = output;
    if (out_path.extension() != ".gpk")
        out_path = std::filesystem::path(out_path.string() + ".gpk");
    std::filesystem::create_directories(out_path.parent_path());

    spdlog::info("Merging {} archives → {}", inputs.size(), out_path.string());

    // ── Open and validate all inputs ─────────────────────────────────────────
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
        // Validate magic
        if (a.mmap.size() < sizeof(FileHeader))
            throw std::runtime_error("merge: too small: " + p.string());
        const auto* fhdr = a.mmap.ptr_at<FileHeader>(0);
        if (fhdr->magic != GPK2_MAGIC)
            throw std::runtime_error("merge: not a v2 .gpk file: " + p.string());
        a.toc = TocReader::read(a.mmap);
        spdlog::info("  Input: {} ({} sections)", p.string(), a.toc.sections.size());
        archives.push_back(std::move(a));
    }

    // ── Compute per-archive genome_id and shard_id offsets ───────────────────
    // Each input's genome_ids and shard_ids are renumbered starting after the
    // previous input's maximum values.
    std::vector<GenomeId> gid_offsets(archives.size(), 0);
    std::vector<uint32_t> sid_offsets(archives.size(), 0);
    {
        GenomeId next_gid = 0;  // next available genome_id after all renumbered ids so far
        uint32_t next_sid = 0;
        for (size_t i = 0; i < archives.size(); ++i) {
            gid_offsets[i] = next_gid;
            sid_offsets[i] = next_sid;
            // Find max genome_id in this archive (original, before renumbering)
            GenomeId local_max_gid = 0;
            for (const auto* sd : archives[i].toc.find_by_type(SEC_CATL)) {
                CatalogSectionReader cr;
                cr.open(archives[i].mmap.data(), sd->file_offset, sd->compressed_size);
                cr.scan([&](const GenomeMeta& m) {
                    if (m.genome_id > local_max_gid) local_max_gid = m.genome_id;
                    return true;
                });
            }
            // Find max shard_id in this archive's TOC
            uint32_t local_max_sid = 0;
            for (const auto& sd : archives[i].toc.sections) {
                if (sd.type == SEC_SHRD) {
                    uint32_t shard_id = static_cast<uint32_t>(sd.aux0);
                    if (shard_id + 1 > local_max_sid) local_max_sid = shard_id + 1;
                }
            }
            // Next archive's ids start after this archive's renumbered maximum
            next_gid = gid_offsets[i] + local_max_gid + 1;
            next_sid = sid_offsets[i] + local_max_sid;
        }
    }

    // ── Open output and write FileHeader ─────────────────────────────────────
    AppendWriter writer;
    writer.create(out_path);
    {
        FileHeader fhdr{};
        fhdr.magic           = GPK2_MAGIC;
        fhdr.version_major   = FORMAT_V2_MAJOR;
        fhdr.version_minor   = FORMAT_V2_MINOR;
        uint64_t t           = static_cast<uint64_t>(std::time(nullptr));
        fhdr.file_uuid_lo    = t ^ 0xdeadbeefcafe0001ULL;
        fhdr.file_uuid_hi    = (t << 17) ^ 0x1234567890abcdefULL;
        fhdr.created_at_unix = t;
        fhdr.flags           = 0;
        std::memset(fhdr.reserved, 0, sizeof(fhdr.reserved));
        writer.append(&fhdr, sizeof(fhdr));
    }

    TocWriter toc_out;
    uint64_t next_section_id = 1;

    // Accumulators for merged CATL and ACCX
    std::vector<GenomeMeta> all_metas;
    std::vector<std::pair<std::string, GenomeId>> all_accessions;

    // ── Copy shard sections from each input ───────────────────────────────────
    size_t total_shards = 0;
    size_t total_genomes = 0;

    for (size_t i = 0; i < archives.size(); ++i) {
        const auto& a        = archives[i];
        GenomeId    gid_off  = gid_offsets[i];
        uint32_t    sid_off  = sid_offsets[i];

        // Copy each SHRD section verbatim
        for (const auto* sd : a.toc.find_by_type(SEC_SHRD)) {
            uint64_t new_offset = writer.current_offset();
            writer.append(a.mmap.data() + sd->file_offset, sd->compressed_size);

            SectionDesc new_sd = *sd;
            new_sd.section_id  = next_section_id++;
            new_sd.file_offset = new_offset;
            new_sd.aux0        = static_cast<uint32_t>(sd->aux0) + sid_off;  // renumber shard_id
            toc_out.add_section(new_sd);
            ++total_shards;
        }

        // Scan CATL — collect renumbered GenomeMeta rows
        for (const auto* sd : a.toc.find_by_type(SEC_CATL)) {
            CatalogSectionReader cr;
            cr.open(a.mmap.data(), sd->file_offset, sd->compressed_size);
            cr.scan([&](const GenomeMeta& m) {
                GenomeMeta rm = m;
                rm.genome_id = m.genome_id + gid_off;
                rm.shard_id  = m.shard_id  + sid_off;
                all_metas.push_back(rm);
                ++total_genomes;
                return true;
            });
        }

        // Scan ACCX — collect renumbered accession→genome_id pairs
        for (const auto* sd : a.toc.find_by_type(SEC_ACCX)) {
            AccessionIndexReader air;
            air.open(a.mmap.data(), sd->file_offset, sd->compressed_size);
            air.scan([&](std::string_view acc, GenomeId gid) {
                all_accessions.emplace_back(std::string(acc), gid + gid_off);
            });
        }
    }

    spdlog::info("Merge: {} shards, {} genomes collected", total_shards, total_genomes);

    // ── Write merged CATL ─────────────────────────────────────────────────────
    uint64_t catalog_root_id = 0;
    {
        CatalogSectionWriter csw;
        for (const auto& m : all_metas)
            csw.add(m);
        SectionDesc catl_sd = csw.finalize(writer, next_section_id++);
        catalog_root_id = catl_sd.section_id;
        toc_out.add_section(catl_sd);
    }

    // ── Write merged ACCX ─────────────────────────────────────────────────────
    uint64_t accession_root_id = 0;
    {
        AccessionIndexWriter aiw;
        for (const auto& [acc, gid] : all_accessions)
            aiw.add(acc, gid);
        SectionDesc accx_sd = aiw.finalize(writer, next_section_id++);
        accession_root_id = accx_sd.section_id;
        toc_out.add_section(accx_sd);
    }

    // ── Write TOC + TailLocator ───────────────────────────────────────────────
    toc_out.finalize(writer,
                     /*generation=*/1,
                     /*live_count=*/total_genomes,
                     /*total_count=*/total_genomes,
                     /*prev_toc_offset=*/0,
                     catalog_root_id,
                     accession_root_id,
                     /*tombstone_root_id=*/0);
    writer.flush();

    spdlog::info("Merge complete: {} → {}", inputs.size(), out_path.string());
}

} // namespace genopack
