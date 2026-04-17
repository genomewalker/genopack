#include <genopack/gtax.hpp>
#include <algorithm>
#include <stdexcept>

namespace genopack {

// ── GtaxWriter ────────────────────────────────────────────────────────────────

void GtaxWriter::add(uint64_t old_taxid, uint64_t new_taxid,
                     GtaxAliasType type, uint32_t source_ver) {
    GtaxEntry e{};
    e.old_taxid  = old_taxid;
    e.new_taxid  = new_taxid;
    e.source_ver = source_ver;
    e.alias_type = static_cast<uint8_t>(type);
    entries_.push_back(e);
}

SectionDesc GtaxWriter::finalize(AppendWriter& writer, uint64_t section_id) {
    // Sort by old_taxid for binary search; stable sort preserves insertion order
    // for equal old_taxids (multiple successors in a split).
    std::stable_sort(entries_.begin(), entries_.end(),
                     [](const GtaxEntry& a, const GtaxEntry& b) {
                         return a.old_taxid < b.old_taxid;
                     });

    GtaxHeader hdr{};
    hdr.magic    = SEC_GTAX;
    hdr.version  = 1;
    hdr.n_entries = entries_.size();

    uint64_t section_start = writer.current_offset();

    writer.append(&hdr, sizeof(hdr));
    if (!entries_.empty())
        writer.append(entries_.data(), sizeof(GtaxEntry) * entries_.size());

    uint64_t payload_size = writer.current_offset() - section_start;

    SectionDesc desc{};
    desc.type              = SEC_GTAX;
    desc.version           = 1;
    desc.section_id        = section_id;
    desc.file_offset       = section_start;
    desc.compressed_size   = payload_size;
    desc.uncompressed_size = payload_size;
    desc.item_count        = entries_.size();
    return desc;
}

// ── GtaxReader ────────────────────────────────────────────────────────────────

void GtaxReader::open(const uint8_t* base, uint64_t offset, uint64_t size) {
    if (size < sizeof(GtaxHeader))
        throw std::runtime_error("GtaxReader: section too small");

    const uint8_t* payload = base + offset;
    hdr_ = reinterpret_cast<const GtaxHeader*>(payload);

    if (hdr_->magic != SEC_GTAX)
        throw std::runtime_error("GtaxReader: bad magic");

    entries_ = reinterpret_cast<const GtaxEntry*>(payload + sizeof(GtaxHeader));
}

size_t GtaxReader::lower_bound(uint64_t old_taxid) const {
    if (!hdr_) return 0;
    size_t n = static_cast<size_t>(hdr_->n_entries);
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (entries_[mid].old_taxid < old_taxid) lo = mid + 1;
        else                                      hi = mid;
    }
    return lo;
}

GtaxResult GtaxReader::resolve(uint64_t old_taxid) const {
    GtaxResult result;

    if (!hdr_ || hdr_->n_entries == 0) {
        result.kind  = GtaxResult::Kind::NOT_FOUND;
        return result;
    }

    size_t n   = static_cast<size_t>(hdr_->n_entries);
    size_t idx = lower_bound(old_taxid);

    if (idx >= n || entries_[idx].old_taxid != old_taxid) {
        result.kind = GtaxResult::Kind::NOT_FOUND;
        return result;
    }

    // Count how many consecutive entries share this old_taxid
    size_t count = 0;
    for (size_t i = idx; i < n && entries_[i].old_taxid == old_taxid; ++i)
        ++count;

    const GtaxEntry& first = entries_[idx];
    auto type = static_cast<GtaxAliasType>(first.alias_type);

    if (type == GtaxAliasType::TOMBSTONED) {
        result.kind  = GtaxResult::Kind::TOMBSTONED;
        result.taxid = 0;
        return result;
    }

    if (type == GtaxAliasType::RETIRED_SPLIT || count > 1) {
        // Multiple successors → AMBIGUOUS (never auto-follow)
        result.kind = GtaxResult::Kind::AMBIGUOUS;
        for (size_t i = idx; i < idx + count; ++i)
            result.successors.push_back(entries_[i].new_taxid);
        return result;
    }

    // Single successor: MERGED or RENAMED
    result.taxid = first.new_taxid;
    result.kind  = (type == GtaxAliasType::RENAMED)
                 ? GtaxResult::Kind::RENAMED
                 : GtaxResult::Kind::MERGED;
    return result;
}

} // namespace genopack
