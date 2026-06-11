#pragma once
#include <filesystem>

namespace genopack {
// Build and append a SEC_GAMI v2 section to each .gpk in pack_path.
// Reads PCORE, builds the global multiplicity index, serialises it.
int cmd_gami_build(const std::filesystem::path& pack_path, int threads);
} // namespace genopack
