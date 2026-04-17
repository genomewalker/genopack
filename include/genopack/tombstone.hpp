#pragma once
#include <genopack/types.hpp>
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <cstdint>
#include <functional>
#include <vector>

namespace genopack {

// ── TombHeader: 16 bytes ──────────────────────────────────────────────────────
// 4+4+8 = 16
struct TombHeader {
    uint32_t magic;       // SEC_TOMB = 0x424D4F54
    uint32_t n_entries;
    uint64_t reserved;
};
static_assert(sizeof(TombHeader) == 16);

// ── TombstoneWriter ───────────────────────────────────────────────────────────

class TombstoneWriter {
public:
    void add(GenomeId id);

    // Sort ids, write TombHeader + sorted GenomeId array.
    SectionDesc finalize(AppendWriter& writer, uint64_t section_id);

private:
    std::vector<GenomeId> ids_;
};

// ── TombstoneReader ───────────────────────────────────────────────────────────

class TombstoneReader {
public:
    void open(const uint8_t* base, uint64_t offset, uint64_t size);

    bool is_deleted(GenomeId id) const;
    void scan(const std::function<void(GenomeId)>& cb) const;

    size_t size() const { return count_; }

private:
    const GenomeId* ids_   = nullptr;
    size_t          count_ = 0;
};

} // namespace genopack
