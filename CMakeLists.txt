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
  src/score_bin.cpp
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
  src/gstx.cpp
  src/gcov.cpp
  src/fmhr.cpp
  src/markers.cpp
  src/markers_build.cpp
  src/cladesplit.cpp
  src/colstore.cpp
  src/core_section.cpp
  src/pcore.cpp
  src/gami.cpp
  src/xqual_columns.cpp
  src/prof.cpp
  src/qcontig.cpp
  src/build_genus_stats.cpp
  src/meta.cpp
  src/oph_sketch.cpp
  src/skani.cpp
  src/scan_engine.cpp
  src/mem_delta.cpp
  src/ncbi_taxdb.cpp
  src/coordinator.cpp
  src/ntdb.cpp
  src/gtax.cpp
  src/sim.cpp
  src/calibrate/calibrate.cpp
)

target_include_directories(genopack_lib PUBLIC
  include
  ${ZSTD_INCLUDE_DIRS}
  ${Xxhash_INCLUDE_DIR}
)
target_include_directories(genopack_lib PRIVATE src)
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

# Eigen3 for check module (TNF/GC-skew contig scoring)
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
  FetchContent_Declare(Eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0
    GIT_SHALLOW    TRUE)
  set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
  set(BUILD_TESTING   OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(Eigen3)
endif()

# Eigen3 is also needed by genopack_lib (transition-operator eigendecomposition in util.cpp)
if(TARGET Eigen3::Eigen)
  target_link_libraries(genopack_lib PUBLIC Eigen3::Eigen)
elseif(DEFINED eigen3_SOURCE_DIR)
  target_include_directories(genopack_lib PUBLIC "${eigen3_SOURCE_DIR}")
endif()

# ── CLI binary ────────────────────────────────────────────────────────────────

add_executable(genopack
  cli/main.cpp
  src/check/pass_a.cpp
  src/check/pass_b.cpp
  src/check/run_check.cpp
  src/run_ingest.cpp
  src/run_report.cpp
  src/run_fcore.cpp
  src/run_pcore.cpp
  src/run_gami.cpp
  src/run_tier.cpp
  src/bench/bench_grid.cpp
  src/calibrate/calibrate.cpp
)
target_include_directories(genopack PRIVATE src)
target_link_libraries(genopack PRIVATE genopack_lib)
if(TARGET Eigen3::Eigen)
  target_link_libraries(genopack PRIVATE Eigen3::Eigen)
elseif(DEFINED eigen3_SOURCE_DIR)
  target_include_directories(genopack PRIVATE "${eigen3_SOURCE_DIR}")
endif()
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

add_executable(genopack_check_roundtrip_test tests/check_roundtrip.cpp)
target_link_libraries(genopack_check_roundtrip_test PRIVATE genopack_lib)
target_compile_definitions(genopack_check_roundtrip_test PRIVATE GENOPACK_BIN="$<TARGET_FILE:genopack>")
add_test(NAME genopack-check-roundtrip COMMAND genopack_check_roundtrip_test)
set_tests_properties(genopack-check-roundtrip PROPERTIES LABELS "unit")

add_executable(genopack_colstore_test tests/colstore_test.cpp)
target_link_libraries(genopack_colstore_test PRIVATE genopack_lib)
add_test(NAME genopack-colstore COMMAND genopack_colstore_test)
set_tests_properties(genopack-colstore PROPERTIES LABELS "unit")

add_executable(genopack_core_section_test tests/core_section_test.cpp)
target_link_libraries(genopack_core_section_test PRIVATE genopack_lib)
add_test(NAME genopack-core-section COMMAND genopack_core_section_test)
set_tests_properties(genopack-core-section PROPERTIES LABELS "unit")

add_executable(genopack_xqual_test tests/xqual_test.cpp)
target_link_libraries(genopack_xqual_test PRIVATE genopack_lib)
add_test(NAME genopack-xqual COMMAND genopack_xqual_test)
set_tests_properties(genopack-xqual PROPERTIES LABELS "unit")

add_executable(genopack_prof_test tests/prof_test.cpp)
target_link_libraries(genopack_prof_test PRIVATE genopack_lib)
add_test(NAME genopack-prof COMMAND genopack_prof_test)
set_tests_properties(genopack-prof PROPERTIES LABELS "unit")

add_executable(genopack_qcontig_test tests/qcontig_test.cpp)
target_link_libraries(genopack_qcontig_test PRIVATE genopack_lib)
add_test(NAME genopack-qcontig COMMAND genopack_qcontig_test)
set_tests_properties(genopack-qcontig PROPERTIES LABELS "unit")

add_executable(genopack_pcore_test tests/pcore_test.cpp)
target_link_libraries(genopack_pcore_test PRIVATE genopack_lib)
add_test(NAME genopack-pcore COMMAND genopack_pcore_test)
set_tests_properties(genopack-pcore PROPERTIES LABELS "unit")
set_tests_properties(genopack-xqual PROPERTIES LABELS "unit")

add_executable(genopack_skch_merge_roundtrip_test tests/skch_merge_roundtrip.cpp)
target_link_libraries(genopack_skch_merge_roundtrip_test PRIVATE genopack_lib)
target_compile_definitions(genopack_skch_merge_roundtrip_test PRIVATE GENOPACK_BIN="$<TARGET_FILE:genopack>")
add_test(NAME genopack-skch-merge-roundtrip COMMAND genopack_skch_merge_roundtrip_test)
set_tests_properties(genopack-skch-merge-roundtrip PROPERTIES LABELS "unit")

add_executable(genopack_repack_roundtrip_test tests/repack_roundtrip.cpp)
target_link_libraries(genopack_repack_roundtrip_test PRIVATE genopack_lib)
target_compile_definitions(genopack_repack_roundtrip_test PRIVATE GENOPACK_BIN="$<TARGET_FILE:genopack>")
add_test(NAME genopack-repack-roundtrip COMMAND genopack_repack_roundtrip_test)
set_tests_properties(genopack-repack-roundtrip PROPERTIES LABELS "unit")

add_executable(genopack_compression_guard_test tests/compression_guard.cpp)
target_link_libraries(genopack_compression_guard_test PRIVATE genopack_lib)
target_compile_definitions(genopack_compression_guard_test PRIVATE GENOPACK_BIN="$<TARGET_FILE:genopack>")
add_test(NAME genopack-compression-guard COMMAND genopack_compression_guard_test)
set_tests_properties(genopack-compression-guard PROPERTIES LABELS "unit")

add_executable(genopack_benchmark tests/benchmark.cpp)
target_link_libraries(genopack_benchmark PRIVATE genopack_lib)
target_compile_definitions(genopack_benchmark PRIVATE GENOPACK_BIN="$<TARGET_FILE:genopack>")

# ── Shared library (.so) ──────────────────────────────────────────────────────
# All symbols pulled from genopack_lib via --whole-archive; stub provides the
# required translation unit for CMake.
set_target_properties(genopack_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(genopack_shared SHARED src/genopack_so_stub.cpp)
set_target_properties(genopack_shared PROPERTIES OUTPUT_NAME "genopack")
target_include_directories(genopack_shared PUBLIC include)
target_link_options(genopack_shared PRIVATE
    -Wl,--whole-archive $<TARGET_FILE:genopack_lib> -Wl,--no-whole-archive)
add_dependencies(genopack_shared genopack_lib)
target_link_libraries(genopack_shared PUBLIC
    ${ZSTD_LIBRARIES} ZLIB::ZLIB pthread)
if(OpenMP_CXX_FOUND)
  target_link_libraries(genopack_shared PUBLIC OpenMP::OpenMP_CXX)
endif()
if(TARGET spdlog::spdlog)
  target_link_libraries(genopack_shared PUBLIC spdlog::spdlog)
endif()
if(TARGET Eigen3::Eigen)
  target_link_libraries(genopack_shared PUBLIC Eigen3::Eigen)
elseif(DEFINED eigen3_SOURCE_DIR)
  target_include_directories(genopack_shared PUBLIC "${eigen3_SOURCE_DIR}")
endif()
if(LIBDEFLATE_FOUND)
  target_link_libraries(genopack_shared PUBLIC ${LIBDEFLATE_LIBRARIES})
endif()

# ── Install ───────────────────────────────────────────────────────────────────
install(TARGETS genopack DESTINATION bin)
install(TARGETS genopack_lib genopack_shared DESTINATION lib)
install(DIRECTORY include/genopack DESTINATION include)
