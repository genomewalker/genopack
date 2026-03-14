#include <genopack/mmap_file.hpp>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cerrno>
#include <cstring>
#include <stdexcept>

namespace genopack {

// ── MmapFileReader ────────────────────────────────────────────────────────────

void MmapFileReader::open(const std::filesystem::path& path) {
    close();

    fd_ = ::open(path.c_str(), O_RDONLY);
    if (fd_ < 0)
        throw std::runtime_error("MmapFileReader: cannot open " + path.string()
                                 + ": " + std::strerror(errno));

    struct stat st{};
    if (::fstat(fd_, &st) < 0) {
        ::close(fd_); fd_ = -1;
        throw std::runtime_error("MmapFileReader: fstat failed: " + std::string(std::strerror(errno)));
    }
    size_ = static_cast<uint64_t>(st.st_size);

    if (size_ == 0) {
        // Empty file: leave data_ as nullptr, keep fd_ open so close() is clean.
        return;
    }

    void* p = ::mmap(nullptr, size_, PROT_READ, MAP_SHARED, fd_, 0);
    if (p == MAP_FAILED) {
        ::close(fd_); fd_ = -1; size_ = 0;
        throw std::runtime_error("MmapFileReader: mmap failed: " + std::string(std::strerror(errno)));
    }
    data_ = static_cast<uint8_t*>(p);
}

void MmapFileReader::close() {
    if (data_) {
        ::munmap(data_, size_);
        data_ = nullptr;
    }
    size_ = 0;
    if (fd_ >= 0) {
        ::close(fd_);
        fd_ = -1;
    }
}

// ── AppendWriter ──────────────────────────────────────────────────────────────

void AppendWriter::create(const std::filesystem::path& path) {
    close();
    fd_ = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd_ < 0)
        throw std::runtime_error("AppendWriter: cannot create " + path.string()
                                 + ": " + std::strerror(errno));
    offset_ = 0;
}

void AppendWriter::open_append(const std::filesystem::path& path) {
    close();
    fd_ = ::open(path.c_str(), O_WRONLY, 0644);
    if (fd_ < 0)
        throw std::runtime_error("AppendWriter: cannot open " + path.string()
                                 + ": " + std::strerror(errno));
    off_t end = ::lseek(fd_, 0, SEEK_END);
    if (end < 0) {
        ::close(fd_); fd_ = -1;
        throw std::runtime_error("AppendWriter: lseek failed: " + std::string(std::strerror(errno)));
    }
    offset_ = static_cast<uint64_t>(end);
}

void AppendWriter::close() {
    if (fd_ >= 0) {
        ::fsync(fd_);
        ::close(fd_);
        fd_ = -1;
    }
    offset_ = 0;
}

uint64_t AppendWriter::append(const void* data, uint64_t len) {
    uint64_t start = offset_;
    const uint8_t* p = static_cast<const uint8_t*>(data);
    uint64_t remaining = len;
    while (remaining > 0) {
        ssize_t n = ::write(fd_, p, remaining);
        if (n < 0) {
            if (errno == EINTR) continue;
            throw std::runtime_error("AppendWriter: write failed: " + std::string(std::strerror(errno)));
        }
        p         += n;
        remaining -= static_cast<uint64_t>(n);
        offset_   += static_cast<uint64_t>(n);
    }
    return start;
}

uint64_t AppendWriter::align(uint64_t alignment) {
    uint64_t rem = offset_ % alignment;
    if (rem == 0) return offset_;
    uint64_t pad = alignment - rem;
    static const uint8_t zeros[64] = {};
    while (pad > 0) {
        uint64_t chunk = (pad < sizeof(zeros)) ? pad : sizeof(zeros);
        append(zeros, chunk);
        pad -= chunk;
    }
    return offset_;
}

void AppendWriter::write_at(uint64_t offset, const void* data, uint64_t len) {
    const uint8_t* p = static_cast<const uint8_t*>(data);
    uint64_t remaining = len;
    off_t pos = static_cast<off_t>(offset);
    while (remaining > 0) {
        ssize_t n = ::pwrite(fd_, p, remaining, pos);
        if (n < 0) {
            if (errno == EINTR) continue;
            throw std::runtime_error("AppendWriter: pwrite failed: " + std::string(std::strerror(errno)));
        }
        p         += n;
        remaining -= static_cast<uint64_t>(n);
        pos       += n;
    }
}

void AppendWriter::flush() {
    if (fd_ >= 0) ::fsync(fd_);
}

} // namespace genopack
