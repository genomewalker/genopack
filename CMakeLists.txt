cmake_minimum_required(VERSION 3.20)
project(genopack VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
add_compile_options(-Wall -Wextra -O3 -march=native -Wno-dangling-reference)

# ── Dependencies ──────────────────────────────────────────────────────────────

# zstd
find_package(PkgConfig REQUIRED)
pkg_check_modules(ZSTD REQUIRED libzstd)

# zlib (for gzip decompression)
find_package(ZLIB REQUIRED)

# spdlog
find_package(spdlog QUIET)
if(NOT spdlog_FOUND)
  find_path(SPDLOG_INCLUDE_DIR spdlog/spdlog.h
    HINTS /usr/include /usr/local/include)
  if(NOT SPDLOG_INCLUDE_DIR)
    # Try building from source via FetchContent
    include(FetchContent)
    FetchContent_Declare(spdlog
      GIT_REPOSITORY https://github.com/gabime/spdlog.git
      GIT_TAG        v1.13.0)
    FetchContent_MakeAvailable(spdlog)
  endif()
endif()

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
  src/appender.cpp
)

target_include_directories(genopack_lib PUBLIC
  include
  ${ZSTD_INCLUDE_DIRS}
)
if(SPDLOG_INCLUDE_DIR)
  target_include_directories(genopack_lib PUBLIC ${SPDLOG_INCLUDE_DIR})
endif()

target_link_libraries(genopack_lib PUBLIC
  ${ZSTD_LIBRARIES}
  ZLIB::ZLIB
)
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

# ── Install ───────────────────────────────────────────────────────────────────
install(TARGETS genopack DESTINATION bin)
install(TARGETS genopack_lib DESTINATION lib)
install(DIRECTORY include/genopack DESTINATION include)
