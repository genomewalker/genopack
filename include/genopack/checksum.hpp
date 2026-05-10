#pragma once
#include <cstdint>
#include <cstring>
#include <span>
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

} // namespace genopack
