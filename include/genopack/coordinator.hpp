#pragma once

#include <genopack/format.hpp>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

namespace genopack {

// ── CoordinatorServer ─────────────────────────────────────────────────────────
// Started by `genopack coordinator -o <file> --nfs-dir <dir> --workers N`.
// Writes FileHeader immediately. Allocates offsets via NFS manifest protocol.
// After all expected workers complete, writes unified TOC + TailLocator.

class CoordinatorServer {
public:
    // NFS manifest mode: workers use 'genopack build --coordinator manifest_dir:output.gpk'.
    // manifest_dir: shared NFS directory for .pending/.alloc/.done coordination files.
    // expected_workers must be > 0.
    void run_nfs(const std::filesystem::path& manifest_dir,
                 const std::filesystem::path& output_path,
                 int      expected_workers,
                 const std::filesystem::path& taxdump_dir = {},
                 const std::function<void(size_t)>& on_progress = {});
};

// NFS worker side: transfer a temp .gpk into the shared output via manifest protocol.
// Writes manifest_dir/{worker_name}.pending, polls for .alloc, pwrite()s sections,
// writes .done, removes tmp_gpk.
void transfer_nfs(const std::filesystem::path& tmp_gpk,
                  const std::filesystem::path& manifest_dir,
                  const std::filesystem::path& output_path,
                  const std::string& worker_name);

} // namespace genopack
