#pragma once
#include <genopack/markers.hpp>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <vector>

namespace genopack::check {

// Thread-safe LRU cache of GenusMarkerIndex, keyed by genus_hash.
// Builds are done outside the lock; reads via shared_ptr are lock-free after construction.
class GenusIndexCache {
    struct Entry {
        uint64_t                          genus_hash = 0;
        std::shared_ptr<GenusMarkerIndex> idx;
        size_t                            lru_tick = 0;
    };
    std::vector<Entry>                  entries_;
    std::unordered_map<uint64_t,size_t> index_;
    std::mutex                          mu_;
    size_t                              tick_    = 0;
    size_t                              max_size_;
public:
    explicit GenusIndexCache(size_t max_size = 32) : max_size_(max_size) {}

    // Get or build the index for a genus. Builder is called outside the lock.
    std::shared_ptr<const GenusMarkerIndex>
    get_or_build(uint64_t genus_hash,
                 std::function<GenusMarkerIndex()> builder)
    {
        {
            std::lock_guard<std::mutex> lk(mu_);
            auto it = index_.find(genus_hash);
            if (it != index_.end()) {
                entries_[it->second].lru_tick = ++tick_;
                return entries_[it->second].idx;
            }
        }

        // Build outside the lock (expensive).
        auto built = std::make_shared<GenusMarkerIndex>(builder());

        std::lock_guard<std::mutex> lk(mu_);
        // Double-check: another thread may have built it concurrently.
        auto it = index_.find(genus_hash);
        if (it != index_.end()) {
            entries_[it->second].lru_tick = ++tick_;
            return entries_[it->second].idx;
        }
        // Evict LRU entry if at capacity.
        if (entries_.size() >= max_size_) {
            size_t victim = 0;
            for (size_t i = 1; i < entries_.size(); ++i)
                if (entries_[i].lru_tick < entries_[victim].lru_tick) victim = i;
            index_.erase(entries_[victim].genus_hash);
            entries_[victim] = Entry{genus_hash, built, ++tick_};
            index_[genus_hash] = victim;
        } else {
            index_[genus_hash] = entries_.size();
            entries_.push_back(Entry{genus_hash, built, ++tick_});
        }
        return built;
    }
};

} // namespace genopack::check
