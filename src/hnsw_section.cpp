#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sys/mman.h>
#include <unistd.h>

#include <genopack/hnsw_section.hpp>
#include <hnswlib/hnswlib.h>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <thread>
#include <vector>

namespace genopack {

// ── HnswSectionWriter ────────────────────────────────────────────────────────

void HnswSectionWriter::build(
    const std::vector<std::pair<GenomeId, std::array<float, 136>>>& profiles,
    uint32_t M, uint32_t ef_construction)
{
    M_  = M;
    ef_ = ef_construction;
    n_elements_ = static_cast<uint32_t>(profiles.size());

    if (n_elements_ == 0) return;

    // InnerProductSpace: for L2-normalised vectors, IP = cosine similarity.
    // hnswlib maximises IP, so nearest = highest cosine similarity.
    hnswlib::InnerProductSpace space(136);
    hnswlib::HierarchicalNSW<float> index(&space, n_elements_, M_, ef_construction);

    label_map_.resize(n_elements_);
    for (uint32_t i = 0; i < n_elements_; ++i)
        label_map_[i] = profiles[i].first;

    // hnswlib addPoint is thread-safe for concurrent adds (not with queries/save).
    // 24 threads balances parallelism vs lock contention on the shared HNSW graph.
    const size_t n_threads = 24;
    std::vector<std::exception_ptr> worker_exceptions(n_threads);
    std::vector<std::thread> workers;
    workers.reserve(n_threads);
    const size_t chunk = (n_elements_ + n_threads - 1) / n_threads;
    for (size_t t = 0; t < n_threads; ++t) {
        const size_t start = t * chunk;
        const size_t end   = std::min(start + chunk, static_cast<size_t>(n_elements_));
        if (start >= end) break;
        workers.emplace_back([&, t, start, end]() {
            try {
                for (size_t i = start; i < end; ++i)
                    index.addPoint(profiles[i].second.data(), static_cast<hnswlib::labeltype>(i));
            } catch (...) {
                worker_exceptions[t] = std::current_exception();
            }
        });
    }
    for (auto& w : workers) w.join();
    for (auto& ex : worker_exceptions)
        if (ex) std::rethrow_exception(ex);

    // Serialize index to in-memory blob via a temp stringstream
    // hnswlib only supports saveIndex(path), so we write to a pipe-style buffer.
    // Use a temporary file approach with a stringstream-backed custom streambuf.
    // Actually, hnswlib v0.8.0 supports saveIndex to std::ostream.
    // Let's use the byte-level serialization.
    {
        // hnswlib::HierarchicalNSW::saveIndex writes to a file path.
        // We'll serialize to a temporary buffer by capturing the file output.
        // The cleanest approach: write to a pipe or use getDataSize + manual serialization.
        // hnswlib exposes: index.saveIndex(const std::string& location)
        // We need to compute the size and serialize. Use a tmpfile approach with memfd.

        // Calculate serialized size using hnswlib internals:
        // offset_level0_ = 0 at label 0 position
        // We'll use a pragmatic approach: serialize to a memory buffer via a
        // temporary file descriptor with memfd_create.
        // Simpler: just use /dev/shm or serialize to a std::ostringstream.

        // hnswlib v0.8.0 has saveIndex(std::ostream&) overload? No, only path.
        // Use memfd_create to get an anonymous in-memory file descriptor.
        int memfd = memfd_create("hnsw_serialize", 0);
        if (memfd < 0)
            throw std::runtime_error("HnswSectionWriter: memfd_create failed");

        // Convert fd to a path hnswlib can open
        std::string fd_path = "/proc/self/fd/" + std::to_string(memfd);
        index.saveIndex(fd_path);

        // Read back the serialized bytes
        off_t blob_size = lseek(memfd, 0, SEEK_END);
        lseek(memfd, 0, SEEK_SET);
        index_blob_.resize(static_cast<size_t>(blob_size));
        size_t total_read = 0;
        while (total_read < index_blob_.size()) {
            ssize_t r = read(memfd, index_blob_.data() + total_read,
                             index_blob_.size() - total_read);
            if (r <= 0) break;
            total_read += static_cast<size_t>(r);
        }
        ::close(memfd);

        if (total_read != index_blob_.size())
            throw std::runtime_error("HnswSectionWriter: incomplete read from memfd");
    }

    spdlog::info("HNSW index built: {} elements, {} bytes, M={}, ef_construction={}",
                 n_elements_, index_blob_.size(), M_, ef_);
}

SectionDesc HnswSectionWriter::finalize(AppendWriter& w, uint64_t section_id) {
    uint64_t section_start = w.current_offset();

    // Layout:
    //   [HnswSectionHeader (64)]
    //   [hnswlib serialized blob]
    //   [uint64_t label_map[n_elements]]

    uint64_t idx_offset = sizeof(HnswSectionHeader);
    uint64_t idx_size   = index_blob_.size();
    uint64_t lm_offset  = idx_offset + idx_size;

    HnswSectionHeader hdr{};
    hdr.magic           = SEC_HNSW;
    hdr.version         = 1;
    hdr.n_elements      = n_elements_;
    hdr.dim             = 136;
    hdr.M               = M_;
    hdr.ef_construction = ef_;
    hdr.index_offset    = idx_offset;
    hdr.index_size      = idx_size;
    hdr.label_map_offset = lm_offset;
    std::memset(hdr.reserved, 0, sizeof(hdr.reserved));

    w.append(&hdr, sizeof(hdr));
    w.append(index_blob_.data(), index_blob_.size());
    w.append(label_map_.data(), label_map_.size() * sizeof(uint64_t));

    uint64_t section_end = w.current_offset();

    SectionDesc desc{};
    desc.type              = SEC_HNSW;
    desc.version           = 1;
    desc.flags             = 0;
    desc.section_id        = section_id;
    desc.file_offset       = section_start;
    desc.compressed_size   = section_end - section_start;
    desc.uncompressed_size = section_end - section_start;
    desc.item_count        = n_elements_;
    desc.aux0              = M_;
    desc.aux1              = ef_;
    std::memset(desc.checksum, 0, sizeof(desc.checksum));
    return desc;
}

// ── HnswSectionReader ────────────────────────────────────────────────────────

struct HnswSectionReader::Impl {
    std::unique_ptr<hnswlib::InnerProductSpace> space;
    std::unique_ptr<hnswlib::HierarchicalNSW<float>> index;
    const uint64_t* label_map = nullptr;
    uint32_t n_elements = 0;
};

void HnswSectionReader::open(const uint8_t* data, uint64_t offset, uint64_t size) {
    if (size < sizeof(HnswSectionHeader))
        throw std::runtime_error("HnswSectionReader: section too small");

    const uint8_t* payload = data + offset;
    const auto* hdr = reinterpret_cast<const HnswSectionHeader*>(payload);

    if (hdr->magic != SEC_HNSW)
        throw std::runtime_error("HnswSectionReader: bad magic");
    if (hdr->dim != 136)
        throw std::runtime_error("HnswSectionReader: unexpected dim");

    impl_ = std::make_shared<Impl>();
    impl_->n_elements = hdr->n_elements;
    impl_->label_map  = reinterpret_cast<const uint64_t*>(payload + hdr->label_map_offset);

    // Load hnswlib index from the serialized blob in mmap'd memory.
    // hnswlib's loadIndex expects a file path, so we use memfd_create to provide
    // the blob as a readable fd.
    impl_->space = std::make_unique<hnswlib::InnerProductSpace>(136);

    int memfd = memfd_create("hnsw_load", 0);
    if (memfd < 0)
        throw std::runtime_error("HnswSectionReader: memfd_create failed");

    const uint8_t* blob_ptr = payload + hdr->index_offset;
    size_t blob_size = hdr->index_size;
    size_t written = 0;
    while (written < blob_size) {
        ssize_t w = write(memfd, blob_ptr + written, blob_size - written);
        if (w <= 0) { ::close(memfd); throw std::runtime_error("HnswSectionReader: write to memfd failed"); }
        written += static_cast<size_t>(w);
    }
    lseek(memfd, 0, SEEK_SET);

    std::string fd_path = "/proc/self/fd/" + std::to_string(memfd);
    impl_->index = std::make_unique<hnswlib::HierarchicalNSW<float>>(
        impl_->space.get(), fd_path, false, hdr->n_elements);
    ::close(memfd);

    open_ = true;
    spdlog::debug("HNSW index loaded: {} elements", hdr->n_elements);
}

std::vector<std::pair<GenomeId, float>>
HnswSectionReader::find_similar(const float* query_profile_136, size_t k) const {
    if (!open_ || !impl_ || impl_->n_elements == 0) return {};

    size_t actual_k = std::min(k, static_cast<size_t>(impl_->n_elements));

    // Set ef to max(64, 2*k) for recall quality
    impl_->index->setEf(std::max(size_t(64), 2 * actual_k));

    auto result = impl_->index->searchKnn(query_profile_136, actual_k);

    // result is a max-heap of (distance, label). For InnerProductSpace,
    // hnswlib returns 1 - IP as distance (since it minimises distance).
    // So cosine_similarity = 1 - distance.
    std::vector<std::pair<GenomeId, float>> out;
    out.reserve(result.size());
    while (!result.empty()) {
        auto [dist, label] = result.top();
        result.pop();
        GenomeId gid = impl_->label_map[label];
        float similarity = 1.0f - dist;
        out.emplace_back(gid, similarity);
    }

    // Sort descending by similarity
    std::sort(out.begin(), out.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    return out;
}

} // namespace genopack
