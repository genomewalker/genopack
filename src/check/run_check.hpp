#pragma once
#include <filesystem>
#include <string>

namespace genopack::check {

int cmd_check(const std::filesystem::path& pack_path,
              const std::filesystem::path& genomes_file,
              int threads,
              int min_genus_size,
              float leakage_threshold,
              const std::filesystem::path& output,
              bool recompute = false);

} // namespace genopack::check
