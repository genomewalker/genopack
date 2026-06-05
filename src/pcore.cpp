#include <genopack/pcore.hpp>
#include <algorithm>
#include <cstring>
#include <unordered_map>
#include <xxhash.h>

namespace genopack {

void PcoreWriter::add_from_members(uint64_t key_hash,
                                   const std::vector<std::vector<uint64_t>>& member_qmers) {
    const auto n_mem = static_cast<uint32_t>(member_qmers.size());
    if (n_mem == 0) return;
    std::unordered_map<uint64_t, uint32_t> cnt;
    size_t est = 0;
    for (const auto& q : member_qmers) est += q.size();
    cnt.reserve(est);
    for (const auto& q : member_qmers)
        for (uint64_t h : q) ++cnt[h];

    Entry e;
    e.key_hash = key_hash;
    e.n_members = n_mem;
    e.aamers.reserve(cnt.size());
    for (const auto& [h, _] : cnt) e.aamers.push_back(h);
    std::sort(e.aamers.begin(), e.aamers.end());
    e.prev.resize(e.aamers.size());
    for (size_t i = 0; i < e.aamers.size(); ++i) {
        const uint32_t c = cnt[e.aamers[i]];
        e.prev[i] = static_cast<uint8_t>(c > 255u ? 255u : c);  // member count (cap 255)
    }
    if (!e.aamers.empty()) entries_.push_back(std::move(e));
}

uint64_t PcoreWriter::model_hash() const {
    std::vector<const Entry*> ord;
    ord.reserve(entries_.size());
    for (const auto& e : entries_) ord.push_back(&e);
    std::sort(ord.begin(), ord.end(),
              [](const Entry* a, const Entry* b) { return a->key_hash < b->key_hash; });
    XXH3_state_t* st = XXH3_createState();
    XXH3_64bits_reset(st);
    XXH3_64bits_update(st, &k_, sizeof(k_));
    XXH3_64bits_update(st, &min_seg_aa_, sizeof(min_seg_aa_));
    XXH3_64bits_update(st, &theta_, sizeof(theta_));
    for (const Entry* e : ord) {
        XXH3_64bits_update(st, &e->key_hash, sizeof(e->key_hash));
        uint32_t n = static_cast<uint32_t>(e->aamers.size());
        XXH3_64bits_update(st, &n, sizeof(n));
        if (n) XXH3_64bits_update(st, e->aamers.data(), n * sizeof(uint64_t));
        if (n) XXH3_64bits_update(st, e->prev.data(), n);
    }
    uint64_t h = XXH3_64bits_digest(st);
    XXH3_freeState(st);
    return h;
}

SectionDesc PcoreWriter::finalize(AppendWriter& w, uint64_t section_id, uint64_t frac_max_hash) {
    const auto n = static_cast<uint32_t>(entries_.size());
    const uint32_t min_buckets = (n == 0) ? 1
        : static_cast<uint32_t>(static_cast<double>(n) / 0.7 + 1.0);
    const uint32_t n_buckets = next_pow2(min_buckets);
    const uint32_t mask      = n_buckets - 1;

    std::vector<uint32_t> buckets(n_buckets, PCORE_EMPTY);
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t slot = static_cast<uint32_t>(entries_[i].key_hash) & mask;
        while (buckets[slot] != PCORE_EMPTY) slot = (slot + 1) & mask;
        buckets[slot] = i;
    }

    const uint64_t entries_off = sizeof(PcoreHeader);
    const uint64_t buckets_off = entries_off + static_cast<uint64_t>(n) * sizeof(PcoreEntry);
    uint64_t pool_off = buckets_off + static_cast<uint64_t>(n_buckets) * sizeof(uint32_t);
    pool_off = (pool_off + 7) & ~uint64_t{7};

    std::vector<PcoreEntry> ents(n);
    uint64_t cursor = pool_off;
    for (uint32_t i = 0; i < n; ++i) {
        const Entry& e = entries_[i];
        const uint32_t na = static_cast<uint32_t>(e.aamers.size());
        cursor = (cursor + 7) & ~uint64_t{7};        // 8-align the u64 aamer array
        ents[i].key_hash      = e.key_hash;
        ents[i].aamers_offset = cursor;
        ents[i].n_aamers      = na;
        ents[i].n_members     = e.n_members;
        cursor += static_cast<uint64_t>(na) * sizeof(uint64_t);
        ents[i].prev_offset   = cursor;
        cursor += na;                                 // u8 prevalence
    }

    w.align(8);
    const uint64_t section_start = w.current_offset();

    PcoreHeader hdr{};
    hdr.magic          = SEC_PCORE;
    hdr.n_entries      = n;
    hdr.n_buckets      = n_buckets;
    hdr.k              = k_;
    hdr.min_seg_aa     = min_seg_aa_;
    hdr.theta          = theta_;
    hdr.frac_max_hash  = frac_max_hash;
    hdr.model_hash     = model_hash();
    hdr.entries_offset = entries_off;
    hdr.buckets_offset = buckets_off;

    w.append(&hdr, sizeof(hdr));
    if (n) w.append(ents.data(), static_cast<uint64_t>(n) * sizeof(PcoreEntry));
    w.append(buckets.data(), static_cast<uint64_t>(n_buckets) * sizeof(uint32_t));
    for (uint32_t i = 0; i < n; ++i) {
        const Entry& e = entries_[i];
        w.align(8);
        if (!e.aamers.empty()) w.append(e.aamers.data(), e.aamers.size() * sizeof(uint64_t));
        if (!e.prev.empty())   w.append(e.prev.data(), e.prev.size());
    }

    const uint64_t section_end = w.current_offset();
    SectionDesc sd{};
    sd.type              = SEC_PCORE;
    sd.version           = 1;
    sd.section_id        = section_id;
    sd.file_offset       = section_start;
    sd.compressed_size   = section_end - section_start;
    sd.uncompressed_size = section_end - section_start;
    sd.item_count        = n;
    sd.aux0              = hdr.model_hash;
    sd.header_size       = sizeof(PcoreHeader);
    return sd;
}

} // namespace genopack
