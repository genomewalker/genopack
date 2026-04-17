#include "genopack/coordinator.hpp"
#include "genopack/mmap_file.hpp"
#include "genopack/ntdb.hpp"
#include "genopack/toc.hpp"

#include <spdlog/spdlog.h>

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstring>
#include <ctime>
#include <fstream>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace genopack {

namespace {

// Encode/decode a SectionDesc as a fixed-size hex blob (160 hex chars = 80 bytes).
static std::string encode_section(const SectionDesc& sd) {
    const uint8_t* p = reinterpret_cast<const uint8_t*>(&sd);
    static constexpr char hex[] = "0123456789abcdef";
    std::string out(sizeof(SectionDesc) * 2, '\0');
    for (size_t i = 0; i < sizeof(SectionDesc); ++i) {
        out[2*i]   = hex[p[i] >> 4];
        out[2*i+1] = hex[p[i] & 0xf];
    }
    return out;
}

static SectionDesc decode_section(const std::string& hex) {
    if (hex.size() != sizeof(SectionDesc) * 2)
        throw std::runtime_error("coordinator: bad SECTION payload length");
    SectionDesc sd{};
    uint8_t* p = reinterpret_cast<uint8_t*>(&sd);
    auto unhex = [](char c) -> uint8_t {
        if (c >= '0' && c <= '9') return c - '0';
        if (c >= 'a' && c <= 'f') return c - 'a' + 10;
        if (c >= 'A' && c <= 'F') return c - 'A' + 10;
        throw std::runtime_error("coordinator: bad hex char");
    };
    for (size_t i = 0; i < sizeof(SectionDesc); ++i)
        p[i] = (unhex(hex[2*i]) << 4) | unhex(hex[2*i+1]);
    return sd;
}

} // anonymous namespace


// ── NFS manifest transfer (worker side) ───────────────────────────────────────

void transfer_nfs(const std::filesystem::path& tmp_gpk,
                  const std::filesystem::path& manifest_dir,
                  const std::filesystem::path& output_path,
                  const std::string& worker_name) {
    MmapFileReader mmap;
    mmap.open(tmp_gpk);
    Toc toc = TocReader::read(mmap);

    uint64_t total_bytes = 0;
    for (auto& sd : toc.sections) total_bytes += sd.compressed_size;

    // Write .pending atomically (coordinator reads total_bytes to allocate space)
    auto pending_tmp = manifest_dir / (worker_name + ".pending.tmp");
    auto pending     = manifest_dir / (worker_name + ".pending");
    {
        std::ofstream f(pending_tmp, std::ios::trunc);
        f << total_bytes << "\n";
    }
    std::filesystem::rename(pending_tmp, pending);
    spdlog::info("nfs-transfer [{}]: pending manifest written ({} bytes, {} sections)",
                 worker_name, total_bytes, toc.sections.size());

    // Poll for .alloc (coordinator writes this once it has allocated our range)
    auto alloc_file = manifest_dir / (worker_name + ".alloc");
    while (!std::filesystem::exists(alloc_file))
        std::this_thread::sleep_for(std::chrono::seconds(5));

    uint64_t start_offset = 0;
    {
        std::ifstream f(alloc_file);
        if (!(f >> start_offset))
            throw std::runtime_error("nfs-transfer: cannot read alloc file: " + alloc_file.string());
    }
    spdlog::info("nfs-transfer [{}]: allocated at offset {}", worker_name, start_offset);

    // pwrite every section into the shared output file at the allocated range
    int out_fd = ::open(output_path.c_str(), O_WRONLY);
    if (out_fd < 0)
        throw std::runtime_error("nfs-transfer: cannot open output: " + output_path.string());

    uint64_t cur = start_offset;
    std::vector<SectionDesc> out_sds;
    for (auto& sd : toc.sections) {
        const uint8_t* src = mmap.data() + sd.file_offset;
        uint64_t rem = sd.compressed_size;
        uint64_t dst = cur;
        while (rem > 0) {
            ssize_t n = ::pwrite(out_fd, src, static_cast<size_t>(rem), static_cast<off_t>(dst));
            if (n <= 0) { ::close(out_fd); throw std::runtime_error("nfs-transfer: pwrite failed"); }
            src += n; dst += n; rem -= n;
        }
        SectionDesc out_sd  = sd;
        out_sd.file_offset  = cur;
        out_sds.push_back(out_sd);
        cur += sd.compressed_size;
    }
    ::fsync(out_fd);
    ::close(out_fd);

    // Write .done atomically (hex-encoded SectionDescs, one per line)
    auto done_tmp  = manifest_dir / (worker_name + ".done.tmp");
    auto done_file = manifest_dir / (worker_name + ".done");
    {
        std::ofstream f(done_tmp, std::ios::trunc);
        for (auto& sd : out_sds) f << encode_section(sd) << "\n";
    }
    std::filesystem::rename(done_tmp, done_file);
    spdlog::info("nfs-transfer [{}]: complete, {} sections written", worker_name, out_sds.size());

    std::filesystem::remove(tmp_gpk);
}

// ── NFS manifest coordinator (server side) ────────────────────────────────────

void CoordinatorServer::run_nfs(
    const std::filesystem::path& manifest_dir,
    const std::filesystem::path& output_path,
    int expected_workers,
    const std::filesystem::path& taxdump_dir,
    const std::function<void(size_t)>& on_progress)
{
    if (expected_workers <= 0)
        throw std::runtime_error("coordinator-nfs: --workers must be > 0");

    // Create output file and write FileHeader
    AppendWriter writer;
    writer.create(output_path);

    FileHeader fhdr{};
    fhdr.magic           = GPK2_MAGIC;
    fhdr.version_major   = FORMAT_MAJOR;
    fhdr.version_minor   = FORMAT_MINOR;
    fhdr.created_at_unix = static_cast<uint64_t>(std::time(nullptr));
    if (FILE* ufd = fopen("/dev/urandom", "rb")) {
        fread(&fhdr.file_uuid_lo, 8, 1, ufd);
        fread(&fhdr.file_uuid_hi, 8, 1, ufd);
        fclose(ufd);
    }
    writer.append(&fhdr, sizeof(fhdr));
    writer.flush();
    uint64_t next_offset = writer.current_offset();

    spdlog::info("coordinator-nfs: output {} created", output_path.string());
    spdlog::info("coordinator-nfs: manifest dir {}, expecting {} workers",
                 manifest_dir.string(), expected_workers);

    // Phase 1: watch for .pending files, allocate contiguous ranges
    std::vector<std::string> worker_names;
    std::map<std::string, bool> allocated;

    while ((int)worker_names.size() < expected_workers) {
        for (auto& entry : std::filesystem::directory_iterator(manifest_dir)) {
            std::string fname = entry.path().filename().string();
            if (!fname.ends_with(".pending")) continue;
            std::string wname = fname.substr(0, fname.size() - 8);
            if (allocated.count(wname)) continue;

            uint64_t total_bytes = 0;
            {
                std::ifstream f(entry.path());
                if (!(f >> total_bytes)) continue; // incomplete write, retry
            }

            uint64_t alloc_start = next_offset;
            next_offset += total_bytes;
            allocated[wname] = true;
            worker_names.push_back(wname);

            auto alloc_tmp  = manifest_dir / (wname + ".alloc.tmp");
            auto alloc_file = manifest_dir / (wname + ".alloc");
            { std::ofstream f(alloc_tmp); f << alloc_start << "\n"; }
            std::filesystem::rename(alloc_tmp, alloc_file);

            spdlog::info("coordinator-nfs: worker {} → {} bytes at offset {}",
                         wname, total_bytes, alloc_start);
        }
        if ((int)worker_names.size() < expected_workers)
            std::this_thread::sleep_for(std::chrono::seconds(5));
    }

    // Pre-allocate output file to cover all worker data (avoids sparse file / NFS holes)
    {
        int fd = ::open(output_path.c_str(), O_WRONLY);
        if (fd >= 0) {
            ::posix_fallocate(fd, 0, static_cast<off_t>(next_offset));
            ::close(fd);
        }
    }
    spdlog::info("coordinator-nfs: all {} workers allocated, waiting for .done files...",
                 expected_workers);

    // Phase 2: watch for .done files, collect SectionDescs
    std::vector<SectionDesc> all_sections;
    std::map<std::string, bool> collected;

    while ((int)collected.size() < expected_workers) {
        for (auto& wname : worker_names) {
            if (collected.count(wname)) continue;
            auto done_file = manifest_dir / (wname + ".done");
            if (!std::filesystem::exists(done_file)) continue;

            std::ifstream f(done_file);
            std::string line;
            while (std::getline(f, line)) {
                if (!line.empty()) all_sections.push_back(decode_section(line));
            }
            collected[wname] = true;
            if (on_progress) on_progress(all_sections.size());
            spdlog::info("coordinator-nfs: worker {} done ({} total sections)",
                         wname, all_sections.size());
        }
        if ((int)collected.size() < expected_workers)
            std::this_thread::sleep_for(std::chrono::seconds(5));
    }

    spdlog::info("coordinator-nfs: all {} workers done, {} sections collected",
                 expected_workers, all_sections.size());

    // Finalize: optional NTDB + unified TOC
    if (!taxdump_dir.empty() &&
        std::filesystem::exists(taxdump_dir / "nodes.dmp") &&
        std::filesystem::exists(taxdump_dir / "names.dmp"))
    {
        writer.seek_to(next_offset);
        writer.align(8);
        uint64_t ntdb_section_id = all_sections.size() + 1;
        NtdbWriter ntdb;
        ntdb.load(taxdump_dir);
        SectionDesc ntdb_desc = ntdb.finalize(writer, ntdb_section_id);
        all_sections.push_back(ntdb_desc);
        next_offset = writer.current_offset();
        spdlog::info("coordinator-nfs: NTDB written ({} bytes compressed)", ntdb_desc.compressed_size);
    } else if (!taxdump_dir.empty()) {
        spdlog::warn("coordinator-nfs: --ntdb dir missing nodes.dmp/names.dmp, skipping NTDB");
    }

    writer.seek_to(next_offset);
    writer.align(8);

    TocWriter toc;
    uint64_t catalog_root_id = 0, accession_root_id = 0, tombstone_root_id = 0;
    uint64_t live_genomes = 0;
    for (auto& sd : all_sections) {
        toc.add_section(sd);
        if (sd.type == SEC_CATL) { catalog_root_id  = sd.section_id; live_genomes += sd.item_count; }
        if (sd.type == SEC_ACCX) accession_root_id = sd.section_id;
        if (sd.type == SEC_TOMB) tombstone_root_id  = sd.section_id;
    }
    toc.finalize(writer, /*generation=*/1, live_genomes, live_genomes,
                 /*prev_toc_offset=*/0,
                 catalog_root_id, accession_root_id, tombstone_root_id);
    writer.flush();

    spdlog::info("coordinator-nfs: {} written successfully", output_path.string());
}

} // namespace genopack
