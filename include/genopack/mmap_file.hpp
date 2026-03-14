#pragma once
#include <cstdint>
#include <filesystem>
#include <span>
#include <stdexcept>

namespace genopack {

// ── MmapFileReader ────────────────────────────────────────────────────────────
// RAII read-only memory-mapped file. Not copyable, movable.

class MmapFileReader {
public:
    MmapFileReader() = default;
    ~MmapFileReader() { close(); }

    MmapFileReader(const MmapFileReader&) = delete;
    MmapFileReader& operator=(const MmapFileReader&) = delete;

    MmapFileReader(MmapFileReader&& o) noexcept
        : data_(o.data_), size_(o.size_), fd_(o.fd_)
    {
        o.data_ = nullptr; o.size_ = 0; o.fd_ = -1;
    }
    MmapFileReader& operator=(MmapFileReader&& o) noexcept {
        if (this != &o) {
            close();
            data_ = o.data_; size_ = o.size_; fd_ = o.fd_;
            o.data_ = nullptr; o.size_ = 0; o.fd_ = -1;
        }
        return *this;
    }

    void open(const std::filesystem::path& path);
    void close();

    const uint8_t* data() const { return data_; }
    uint64_t       size() const { return size_; }

    std::span<const uint8_t> view(uint64_t offset, uint64_t len) const {
        if (offset + len > size_)
            throw std::out_of_range("MmapFileReader::view out of bounds");
        return {data_ + offset, static_cast<size_t>(len)};
    }

    template<typename T>
    const T* ptr_at(uint64_t offset) const {
        if (offset + sizeof(T) > size_)
            throw std::out_of_range("MmapFileReader::ptr_at out of bounds");
        return reinterpret_cast<const T*>(data_ + offset);
    }

private:
    uint8_t* data_ = nullptr;
    uint64_t size_ = 0;
    int      fd_   = -1;
};

// ── AppendWriter ──────────────────────────────────────────────────────────────
// Sequential append writer with random-write patch support.

class AppendWriter {
public:
    AppendWriter() = default;
    ~AppendWriter() { close(); }

    AppendWriter(const AppendWriter&) = delete;
    AppendWriter& operator=(const AppendWriter&) = delete;

    void create(const std::filesystem::path& path);
    void open_append(const std::filesystem::path& path);
    void close();

    // Append bytes; returns start offset of this write.
    uint64_t append(const void* data, uint64_t len);

    // Pad to alignment boundary with zero bytes; returns new offset.
    uint64_t align(uint64_t alignment = 8);

    // Patch already-written region (pwrite).
    void write_at(uint64_t offset, const void* data, uint64_t len);

    uint64_t current_offset() const { return offset_; }

    void flush();

private:
    int      fd_     = -1;
    uint64_t offset_ = 0;
};

} // namespace genopack
