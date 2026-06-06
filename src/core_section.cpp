#include <genopack/core_section.hpp>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <xxhash.h>

namespace genopack {

void CoreWriter::add(uint64_t genus_hash, std::vector<uint64_t> core_aamers,
                     std::vector<uint64_t> member_ids) {
    if (core_aamers.empty()) return;
    std::sort(member_ids.begin(), member_ids.end());
    entries_.push_back({genus_hash, std::move(core_aamers), std::move(member_ids)});
}

void CoreWriter::add_from_members(uint64_t genus_hash,
                                  const std::vector<std::vector<uint64_t>>& member_qmers,
                                  std::vector<uint64_t> member_ids) {
    const auto n_mem = static_cast<uint32_t>(member_qmers.size());
    if (n_mem == 0) return;
    const auto thr = static_cast<uint32_t>(std::ceil(static_cast<double>(theta_) * n_mem));

    size_t est = 0;
    for (const auto& q : member_qmers) est = std::max(est, q.size());
    std::unordered_map<uint64_t, uint32_t> cnt;
    cnt.reserve(est * 2 + 16);
    for (const auto& q : member_qmers)
        for (uint64_t h : q) ++cnt[h];

    std::vector<uint64_t> core;
    core.reserve(est);
    for (const auto& [h, c] : cnt)
        if (c >= thr) core.push_back(h);
    std::sort(core.begin(), core.end());
    add(genus_hash, std::move(core), std::move(member_ids));
}

void CoreWriter::add_from_counts(uint64_t genus_hash,
                                 const std::unordered_map<uint64_t, uint32_t>& cnt,
                                 uint32_t n_members, std::vector<uint64_t> member_ids) {
    if (n_members == 0 || cnt.empty()) return;
    const auto thr = static_cast<uint32_t>(std::ceil(static_cast<double>(theta_) * n_members));
    std::vector<uint64_t> core;
    core.reserve(cnt.size());
    for (const auto& [h, c] : cnt)
        if (c >= thr) core.push_back(h);
    std::sort(core.begin(), core.end());
    add(genus_hash, std::move(core), std::move(member_ids));
}

uint64_t CoreWriter::core_model_hash() const {
    std::vector<const Entry*> ord;
    ord.reserve(entries_.size());
    for (const auto& e : entries_) ord.push_back(&e);
    std::sort(ord.begin(), ord.end(),
              [](const Entry* a, const Entry* b) { return a->genus_hash < b->genus_hash; });

    XXH3_state_t* st = XXH3_createState();
    XXH3_64bits_reset(st);
    uint32_t theta_bits;
    std::memcpy(&theta_bits, &theta_, 4);
    XXH3_64bits_update(st, &k_, sizeof(k_));
    XXH3_64bits_update(st, &min_seg_aa_, sizeof(min_seg_aa_));
    XXH3_64bits_update(st, &theta_bits, sizeof(theta_bits));
    for (const Entry* e : ord) {
        XXH3_64bits_update(st, &e->genus_hash, sizeof(e->genus_hash));
        uint32_t n = static_cast<uint32_t>(e->aamers.size());
        XXH3_64bits_update(st, &n, sizeof(n));
        if (n) XXH3_64bits_update(st, e->aamers.data(), n * sizeof(uint64_t));
    }
    uint64_t h = XXH3_64bits_digest(st);
    XXH3_freeState(st);
    return h;
}

SectionDesc CoreWriter::finalize(AppendWriter& w, uint64_t section_id, uint64_t frac_max_hash) {
    const auto n = static_cast<uint32_t>(entries_.size());
    const uint32_t min_buckets = (n == 0) ? 1
        : static_cast<uint32_t>(static_cast<double>(n) / 0.7 + 1.0);
    const uint32_t n_buckets = next_pow2(min_buckets);
    const uint32_t mask      = n_buckets - 1;

    std::vector<uint32_t> buckets(n_buckets, CORE_EMPTY);
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t slot = static_cast<uint32_t>(entries_[i].genus_hash) & mask;
        while (buckets[slot] != CORE_EMPTY) slot = (slot + 1) & mask;
        buckets[slot] = i;
    }

    const uint64_t entries_off = sizeof(CoreHeader);
    const uint64_t buckets_off = entries_off + static_cast<uint64_t>(n) * sizeof(CoreEntry);
    uint64_t pool_off = buckets_off + static_cast<uint64_t>(n_buckets) * sizeof(uint32_t);
    pool_off = (pool_off + 7) & ~uint64_t{7};   // 8-align the u64 pool

    // Assign each genus's pool offsets (aamers then member ids), in entry order.
    std::vector<CoreEntry> ents(n);
    uint64_t cursor = pool_off;
    for (uint32_t i = 0; i < n; ++i) {
        const Entry& e = entries_[i];
        ents[i].genus_hash    = e.genus_hash;
        ents[i].aamers_offset = cursor;
        ents[i].n_aamers      = static_cast<uint32_t>(e.aamers.size());
        cursor += static_cast<uint64_t>(e.aamers.size()) * sizeof(uint64_t);
        ents[i].ids_offset    = cursor;
        ents[i].n_members     = static_cast<uint32_t>(e.members.size());
        cursor += static_cast<uint64_t>(e.members.size()) * sizeof(uint64_t);
    }

    w.align(8);
    const uint64_t section_start = w.current_offset();

    CoreHeader hdr{};
    hdr.magic           = section_type_;
    hdr.n_genera        = n;
    hdr.n_buckets       = n_buckets;
    hdr.k               = k_;
    hdr.min_seg_aa      = min_seg_aa_;
    hdr.theta           = theta_;
    hdr.frac_max_hash   = frac_max_hash;
    hdr.core_model_hash = core_model_hash();
    hdr.entries_offset  = entries_off;
    hdr.buckets_offset  = buckets_off;

    w.append(&hdr, sizeof(hdr));
    if (n) w.append(ents.data(), static_cast<uint64_t>(n) * sizeof(CoreEntry));
    w.append(buckets.data(), static_cast<uint64_t>(n_buckets) * sizeof(uint32_t));
    w.align(8);   // reach pool_off
    for (uint32_t i = 0; i < n; ++i) {
        const Entry& e = entries_[i];
        if (!e.aamers.empty())
            w.append(e.aamers.data(), e.aamers.size() * sizeof(uint64_t));
        if (!e.members.empty())
            w.append(e.members.data(), e.members.size() * sizeof(uint64_t));
    }

    const uint64_t section_end = w.current_offset();
    SectionDesc sd{};
    sd.type              = section_type_;
    sd.version           = 1;
    sd.flags             = 0;
    sd.section_id        = section_id;
    sd.file_offset       = section_start;
    sd.compressed_size   = section_end - section_start;
    sd.uncompressed_size = section_end - section_start;
    sd.item_count        = n;
    sd.aux0              = hdr.core_model_hash;   // quick directory-level model id
    sd.header_size       = sizeof(CoreHeader);
    return sd;
}

} // namespace genopack
