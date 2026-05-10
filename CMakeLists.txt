cmake_minimum_required(VERSION 3.20)
project(genopack VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
add_compile_options(-Wall -Wextra -O3 -mavx2 -mfma -Wno-dangling-reference)

# OpenMP
find_package(OpenMP)

# ── Dependencies ──────────────────────────────────────────────────────────────

# zstd
find_package(PkgConfig REQUIRED)
pkg_check_modules(ZSTD REQUIRED libzstd)

# zlib (for gzip decompression)
find_package(ZLIB REQUIRED)

# xxhash — required by derep_view (.gpd reader)
find_path(Xxhash_INCLUDE_DIR xxhash.h
  HINTS /usr/local/include /usr/include
        /maps/projects/fernandezguerra/apps/opt/conda/envs/bioinfo/include
        ENV CONDA_PREFIX
  PATH_SUFFIXES include)
if(NOT Xxhash_INCLUDE_DIR)
  message(FATAL_ERROR "xxhash.h not found. Install xxhash or set Xxhash_INCLUDE_DIR.")
endif()
message(STATUS "xxhash include: ${Xxhash_INCLUDE_DIR}")

# libdeflate (fast gzip decompression, ~3x faster than zlib)
pkg_check_modules(LIBDEFLATE libdeflate)
if(NOT LIBDEFLATE_FOUND)
  find_library(LIBDEFLATE_LIBRARIES NAMES deflate libdeflate
    HINTS ${CMAKE_PREFIX_PATH}/lib /usr/lib /usr/local/lib)
  find_path(LIBDEFLATE_INCLUDE_DIRS libdeflate.h
    HINTS ${CMAKE_PREFIX_PATH}/include /usr/include /usr/local/include)
  if(LIBDEFLATE_LIBRARIES AND LIBDEFLATE_INCLUDE_DIRS)
    set(LIBDEFLATE_FOUND TRUE)
  endif()
endif()

# spdlog — always fetch v1.14.1 to match geodesic's fmt version (fmt v10)
include(FetchContent)
FetchContent_Declare(spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG        v1.14.1)
FetchContent_MakeAvailable(spdlog)

# CLI11
find_path(CLI11_INCLUDE_DIR CLI/CLI.hpp
  HINTS /usr/include /usr/local/include)
if(NOT CLI11_INCLUDE_DIR)
  include(FetchContent)
  FetchContent_Declare(CLI11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG        v2.4.2)
  FetchContent_MakeAvailable(CLI11)
endif()

# ── Library ───────────────────────────────────────────────────────────────────

add_library(genopack_lib STATIC
  src/catalog.cpp
  src/shard.cpp
  src/builder.cpp
  src/archive.cpp
  src/archive_set_reader.cpp
  src/derep_view.cpp
  src/appender.cpp
  src/mmap_file.cpp
  src/util.cpp
  src/accx.cpp
  src/cidx.cpp
  src/taxn.cpp
  src/txdb.cpp
  src/tombstone.cpp
  src/toc.cpp
  src/merger.cpp
  src/repack.cpp
  src/kmrx.cpp
  src/skch.cpp
  src/gidx.cpp
  src/oph_sketch.cpp
  src/scan_engine.cpp
  src/mem_delta.cpp
  src/ncbi_taxdb.cpp
  src/coordinator.cpp
  src/ntdb.cpp
  src/gtax.cpp
)

target_include_directories(genopack_lib PUBLIC
  include
  ${ZSTD_INCLUDE_DIRS}
  ${Xxhash_INCLUDE_DIR}
)
if(SPDLOG_INCLUDE_DIR)
  target_include_directories(genopack_lib PUBLIC ${SPDLOG_INCLUDE_DIR})
endif()

target_link_libraries(genopack_lib PUBLIC
  ${ZSTD_LIBRARIES}
  ZLIB::ZLIB
  pthread
)
if(OpenMP_CXX_FOUND)
  target_link_libraries(genopack_lib PUBLIC OpenMP::OpenMP_CXX)
  message(STATUS "OpenMP found: parallel shard decompression enabled")
endif()
if(LIBDEFLATE_FOUND)
  target_compile_definitions(genopack_lib PUBLIC GENOPACK_HAVE_LIBDEFLATE=1)
  target_include_directories(genopack_lib PUBLIC ${LIBDEFLATE_INCLUDE_DIRS})
  target_link_libraries(genopack_lib PUBLIC ${LIBDEFLATE_LIBRARIES})
  message(STATUS "libdeflate found: fast gzip decompression enabled")
else()
  message(STATUS "libdeflate not found: using zlib for gzip decompression")
endif()
if(TARGET spdlog::spdlog)
  target_link_libraries(genopack_lib PUBLIC spdlog::spdlog)
endif()

# ── CLI binary ────────────────────────────────────────────────────────────────

add_executable(genopack cli/main.cpp)
target_link_libraries(genopack PRIVATE genopack_lib)
if(CLI11_INCLUDE_DIR)
  target_include_directories(genopack PRIVATE ${CLI11_INCLUDE_DIR})
endif()
if(TARGET CLI11::CLI11)
  target_link_libraries(genopack PRIVATE CLI11::CLI11)
endif()

# ── Tests ─────────────────────────────────────────────────────────────────────

enable_testing()

add_executable(genopack_smoke_test tests/smoke.cpp)
target_link_libraries(genopack_smoke_test PRIVATE genopack_lib)
target_compile_definitions(genopack_smoke_test PRIVATE GENOPACK_BIN="$<TARGET_FILE:genopack>")
add_test(NAME genopack-smoke COMMAND genopack_smoke_test)

add_executable(genopack_mem_delta_test tests/mem_delta.cpp)
target_link_libraries(genopack_mem_delta_test PRIVATE genopack_lib)
target_compile_definitions(genopack_mem_delta_test PRIVATE GENOPACK_BIN="$<TARGET_FILE:genopack>")
add_test(NAME genopack-mem-delta COMMAND genopack_mem_delta_test)

add_executable(genopack_skch_dual_seed_test tests/skch_dual_seed.cpp)
target_link_libraries(genopack_skch_dual_seed_test PRIVATE genopack_lib)
add_test(NAME genopack-skch-dual-seed COMMAND genopack_skch_dual_seed_test)
set_tests_properties(genopack-skch-dual-seed PROPERTIES LABELS "unit")

add_executable(genopack_derep_view_test tests/derep_view_test.cpp)
target_link_libraries(genopack_derep_view_test PRIVATE genopack_lib)
add_test(NAME genopack-derep-view COMMAND genopack_derep_view_test)
set_tests_properties(genopack-derep-view PROPERTIES LABELS "unit")

add_executable(genopack_benchmark tests/benchmark.cpp)
target_link_libraries(genopack_benchmark PRIVATE genopack_lib)
target_compile_definitions(genopack_benchmark PRIVATE GENOPACK_BIN="$<TARGET_FILE:genopack>")

# ── Install ───────────────────────────────────────────────────────────────────
install(TARGETS genopack DESTINATION bin)
install(TARGETS genopack_lib DESTINATION lib)
install(DIRECTORY include/genopack DESTINATION include)
