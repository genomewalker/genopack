#pragma once
#include <cstdint>
#include <filesystem>
#include <string>

namespace genopack::bench {

int cmd_bench_grid(
    const std::filesystem::path& archive_path,
    const std::filesystem::path& manifest_path,
    const std::filesystem::path& output,
    int      threads,
    int      reps,
    uint32_t seed);

} // namespace genopack::bench
