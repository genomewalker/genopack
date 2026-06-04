#pragma once
#include <cstdint>
#include <cstring>
#include <span>
#include <algorithm>
#include <vector>
#include <unistd.h>
#include <xxhash.h>

namespace genopack {

// Compute XXH128 over a byte buffer with the 16-byte checksum field at
// `checksum_offset` zeroed during hashing. Store the result back at that offset.
// This is the canonical write-side operation: assemble bytes, call this, done.
inline void compute_checksum(uint8_t* data, size_t size, size_t checksum_offset) {
    uint8_t saved[16];
    std::memcpy(saved, data + checksum_offset, 16);
    std::memset(data + checksum_offset, 0, 16);
    XXH128_hash_t h = XXH3_128bits(data, size);
    // Restore original bytes first in case caller wants to reuse buffer,
    // then overwrite with the computed hash in little-endian order.
    XXH128_canonical_t canon;
    XXH128_canonicalFromHash(&canon, h);
    std::memcpy(data + checksum_offset, canon.digest, 16);
}

// Verify: zero the checksum field, recompute, compare. Returns true if matches.
inline bool verify_checksum(const uint8_t* data, size_t size, size_t checksum_offset) {
    // Copy buffer to zero the checksum field without modifying the original mmap.
    // For large shards this would be expensive — callers should use the streaming
    // variant below when operating on mmap'd data.
    std::vector<uint8_t> buf(data, data + size);
    std::memset(buf.data() + checksum_offset, 0, 16);
    XXH128_hash_t h = XXH3_128bits(buf.data(), size);
    XXH128_canonical_t canon;
    XXH128_canonicalFromHash(&canon, h);
    return std::memcmp(data + checksum_offset, canon.digest, 16) == 0;
}

// Compute XXH128 over a standalone section body (no embedded checksum field) and
// store the 16-byte canonical digest into out16. This is the SectionDesc.checksum
// operation: the digest lives in the descriptor, separate from the body.
inline void checksum_of(const uint8_t* data, size_t size, uint8_t out16[16]) {
    XXH128_hash_t h = XXH3_128bits(data, size);
    XXH128_canonical_t canon;
    XXH128_canonicalFromHash(&canon, h);
    std::memcpy(out16, canon.digest, 16);
}

// Streaming XXH3-128 over a file region [offset, offset+size), read in chunks so
// multi-GB sections are hashed without a full in-RAM copy. Writes the 16-byte
// canonical digest to out16; returns false on a short/failed read (out16 left
// untouched). Equivalent to checksum_of() over the same bytes.
inline bool checksum_of_fd(int fd, uint64_t offset, uint64_t size, uint8_t out16[16]) {
    XXH3_state_t* st = XXH3_createState();
    if (!st) return false;
    XXH3_128bits_reset(st);
    constexpr size_t CHUNK = 8u << 20;   // 8 MB
    std::vector<uint8_t> buf(static_cast<size_t>(std::min<uint64_t>(CHUNK, size ? size : 1)));
    uint64_t done = 0;
    bool ok = true;
    while (done < size) {
        size_t want = static_cast<size_t>(std::min<uint64_t>(buf.size(), size - done));
        ssize_t nr = ::pread(fd, buf.data(), want, static_cast<off_t>(offset + done));
        if (nr <= 0) { ok = false; break; }
        XXH3_128bits_update(st, buf.data(), static_cast<size_t>(nr));
        done += static_cast<uint64_t>(nr);
    }
    if (ok) {
        XXH128_canonical_t canon;
        XXH128_canonicalFromHash(&canon, XXH3_128bits_digest(st));
        std::memcpy(out16, canon.digest, 16);
    }
    XXH3_freeState(st);
    return ok;
}

// True if XXH128(data,size) equals the 16-byte canonical digest in expect16.
inline bool checksum_matches(const uint8_t* data, size_t size, const uint8_t expect16[16]) {
    XXH128_hash_t h = XXH3_128bits(data, size);
    XXH128_canonical_t canon;
    XXH128_canonicalFromHash(&canon, h);
    return std::memcmp(expect16, canon.digest, 16) == 0;
}

} // namespace genopack
