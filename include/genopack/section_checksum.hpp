#pragma once
#include "format.hpp"
#include "checksum.hpp"
#include <cstdint>
#include <cstring>
#include <span>
#include <vector>
#include <fcntl.h>
#include <unistd.h>

namespace genopack {

// Stamp SectionDesc.checksum (content_xxh128) for sections whose bytes live in
// the file at `path` at their absolute file_offset. This is the canonical P1
// integrity pass shared by build/merge/repack/reindex so their metadata bodies
// are content-hashed (otherwise verify reports them as "no checksum" and skips
// them, masking corruption).
//
// Hashing is streamed (checksum_of_fd) so multi-GB sections are covered without
// a full in-RAM copy — there is no size cap. Skips SEC_SHRD (hashed inline at
// write time) and empty sections. When `only_if_zero` is set, sections already
// carrying a non-zero checksum are left untouched, so large already-hashed
// sections (e.g. SHRD/SKCH copied through reindex) are not needlessly re-read.
inline void stamp_section_checksums(const char* path,
                                    std::span<SectionDesc> secs,
                                    bool only_if_zero = false) {
    int fd = ::open(path, O_RDONLY);
    if (fd < 0) return;
    static const uint8_t zeros[16] = {};
    for (auto& sd : secs) {
        if (sd.type == SEC_SHRD) continue;                  // hashed inline at write time
        if (sd.compressed_size == 0) continue;
        if (only_if_zero && std::memcmp(sd.checksum, zeros, 16) != 0) continue;
        checksum_of_fd(fd, sd.file_offset, sd.compressed_size, sd.checksum);
    }
    ::close(fd);
}

} // namespace genopack
