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
// integrity pass: build and merge inline their own equivalent loop; repack and
// reindex call this so their metadata bodies are content-hashed too (otherwise
// verify reports them as "no checksum" and skips them, masking corruption).
//
// Skips SEC_SHRD (hashed inline at write time), empty sections, and sections
// larger than `cap` (left zero to avoid multi-GB re-reads — the same bound
// build/merge use). When `only_if_zero` is set, sections already carrying a
// non-zero checksum are left untouched, so large already-hashed sections (e.g.
// SHRD/SKCH copied through reindex) are not needlessly re-read.
inline void stamp_section_checksums(const char* path,
                                    std::span<SectionDesc> secs,
                                    uint64_t cap = (512ull << 20),
                                    bool only_if_zero = false) {
    int fd = ::open(path, O_RDONLY);
    if (fd < 0) return;
    static const uint8_t zeros[16] = {};
    std::vector<uint8_t> sbuf;
    for (auto& sd : secs) {
        if (sd.type == SEC_SHRD) continue;                  // hashed inline at write time
        if (sd.compressed_size == 0) continue;
        if (sd.compressed_size > cap) continue;             // over-cap: left zero like build/merge
        if (only_if_zero && std::memcmp(sd.checksum, zeros, 16) != 0) continue;
        sbuf.resize(static_cast<size_t>(sd.compressed_size));
        ssize_t nr = ::pread(fd, sbuf.data(), sbuf.size(),
                             static_cast<off_t>(sd.file_offset));
        if (nr == static_cast<ssize_t>(sbuf.size()))
            checksum_of(sbuf.data(), sbuf.size(), sd.checksum);
    }
    ::close(fd);
}

} // namespace genopack
