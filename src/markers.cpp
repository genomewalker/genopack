#include <genopack/markers.hpp>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <numeric>
#include <stdexcept>

namespace genopack {

// ── MarkerWriter::finalize ────────────────────────────────────────────────────

void MarkerWriter::finalize(const std::filesystem::path& path,
                            uint8_t n_bac, uint8_t n_arc, uint8_t k,
                            uint16_t frac_scale, uint8_t alphabet) const {
    // Sort lineages by hash for binary search at read time.
    std::vector<size_t> order(lineages_.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](size_t a, size_t b) { return lineages_[a].hash < lineages_[b].hash; });

    // ── Pass 1: compute layout offsets ───────────────────────────────────────

    const uint32_t n_lin    = static_cast<uint32_t>(lineages_.size());
    const uint32_t n_pool   = static_cast<uint32_t>(n_bac) + static_cast<uint32_t>(n_arc);
    const uint32_t hdr_size = sizeof(MarkerHeader);

    const uint32_t lookup_off   = hdr_size;
    const uint32_t pool_idx_off = lookup_off   + n_lin  * sizeof(MarkerLookupEntry);
    const uint32_t calib_off    = pool_idx_off + n_pool * sizeof(MarkerPoolEntry);

    // Build packed calib data and per-lineage offsets.
    std::vector<uint8_t> calib_data;
    std::vector<uint32_t> calib_offsets(n_lin);

    for (size_t si = 0; si < n_lin; ++si) {
        const auto& e = lineages_[order[si]];
        const uint8_t n_markers = static_cast<uint8_t>(e.slots.size());
        const uint32_t bitset_bytes = (n_markers + 7u) / 8u;

        calib_offsets[si] = static_cast<uint32_t>(calib_data.size());

        CalibEntryHeader ch{};
        ch.n_markers = n_markers;
        ch.ref_count = e.ref_count;
        ch.domain    = e.domain;
        const uint8_t* chp = reinterpret_cast<const uint8_t*>(&ch);
        calib_data.insert(calib_data.end(), chp, chp + sizeof(ch));

        for (const auto& slot : e.slots) {
            const uint8_t* sp = reinterpret_cast<const uint8_t*>(&slot);
            calib_data.insert(calib_data.end(), sp, sp + sizeof(slot));
        }

        // Pack expected_markers bitset.
        std::vector<uint8_t> bits(bitset_bytes, 0);
        for (int i = 0; i < (int)e.expected.size() && i < n_markers; ++i)
            if (e.expected[i]) bits[i / 8] |= (1u << (i % 8));
        calib_data.insert(calib_data.end(), bits.begin(), bits.end());
    }

    // Align calib section end to 8 bytes.
    while (calib_data.size() % 8 != 0) calib_data.push_back(0);

    const uint32_t pool_off = calib_off + static_cast<uint32_t>(calib_data.size());

    // Compute pool section offsets and total hash bytes.
    std::vector<uint64_t> pool_hash_offsets(n_pool, 0);
    uint64_t pool_bytes = 0;
    for (uint32_t mi = 0; mi < n_pool; ++mi) {
        pool_hash_offsets[mi] = pool_bytes;
        const size_t nh = (mi < pool_.size()) ? pool_[mi].size() : 0;
        pool_bytes += nh * sizeof(uint64_t);
    }

    // Pre-merged pool section: k-way merge done at build time so loading is O(1).
    // Layout: [N_bac uint64 hashes][N_bac uint8 ids] then [N_arc uint64][N_arc uint8]
    // Aligned to 8 bytes. Stored after pool hash arrays.
    const uint64_t pool_end = static_cast<uint64_t>(pool_off) + pool_bytes;
    // Align merged section to 8 bytes.
    const uint64_t merged_bac_off_u64 = (pool_end + 7) & ~uint64_t{7};

    // Build merged bac pool (k-way merge of the per-marker arrays).
    auto build_merged_section = [&](uint8_t base, uint8_t count)
        -> std::pair<std::vector<uint64_t>, std::vector<uint8_t>> {
        std::vector<std::pair<uint64_t,uint8_t>> pairs;
        size_t total = 0;
        for (uint8_t i = 0; i < count; ++i)
            if (base + i < pool_.size()) total += pool_[base + i].size();
        pairs.reserve(total);
        for (uint8_t i = 0; i < count; ++i) {
            if (base + i >= pool_.size()) continue;
            for (uint64_t h : pool_[base + i]) pairs.emplace_back(h, i);
        }
        std::stable_sort(pairs.begin(), pairs.end(),
                         [](const auto& a, const auto& b) { return a.first < b.first; });
        std::vector<uint64_t> hashes(pairs.size());
        std::vector<uint8_t>  ids(pairs.size());
        for (size_t i = 0; i < pairs.size(); ++i) {
            hashes[i] = pairs[i].first;
            ids[i]    = pairs[i].second;
        }
        return {std::move(hashes), std::move(ids)};
    };

    auto [merged_bac_h, merged_bac_id] = build_merged_section(0,     n_bac);
    auto [merged_arc_h, merged_arc_id] = build_merged_section(n_bac, n_arc);

    // Merged arc section follows immediately after bac.
    const uint64_t bac_section_bytes = merged_bac_h.size() * 8 + merged_bac_h.size();
    const uint64_t bac_section_padded = (bac_section_bytes + 7) & ~uint64_t{7};
    const uint64_t merged_arc_off_u64 = merged_arc_h.empty() ? 0
                                       : merged_bac_off_u64 + bac_section_padded;

    const uint32_t merged_bac_off = static_cast<uint32_t>(merged_bac_off_u64);
    const uint32_t merged_arc_off = static_cast<uint32_t>(merged_arc_off_u64);

    // ── Pass 2: write ─────────────────────────────────────────────────────────

    std::ofstream f(path, std::ios::binary | std::ios::trunc);
    if (!f) throw std::runtime_error("markers: cannot open output: " + path.string());

    auto write = [&](const void* p, size_t n) {
        f.write(reinterpret_cast<const char*>(p), static_cast<std::streamsize>(n));
    };

    MarkerHeader hdr{};
    hdr.magic         = MRKR_MAGIC;
    hdr.version       = MRKR_VERSION;
    hdr.k             = k;
    hdr.alphabet      = alphabet;
    hdr.n_lineages    = n_lin;
    hdr.n_bac_markers = n_bac;
    hdr.n_arc_markers = n_arc;
    hdr.lookup_off    = lookup_off;
    hdr.pool_idx_off  = pool_idx_off;
    hdr.calib_off     = calib_off;
    hdr.pool_off      = pool_off;
    hdr.frac_scale      = frac_scale;
    hdr.merged_bac_off  = merged_bac_off;
    hdr.merged_arc_off  = merged_arc_off;
    write(&hdr, sizeof(hdr));

    // Sorted lookup table.
    for (size_t si = 0; si < n_lin; ++si) {
        MarkerLookupEntry le{};
        le.lineage_hash   = lineages_[order[si]].hash;
        le.calib_data_off = calib_offsets[si];
        le.domain         = lineages_[order[si]].domain;
        write(&le, sizeof(le));
    }

    // Pool index.
    for (uint32_t mi = 0; mi < n_pool; ++mi) {
        MarkerPoolEntry pe{};
        pe.hashes_off = pool_hash_offsets[mi];
        pe.n_hashes   = static_cast<uint32_t>((mi < pool_.size()) ? pool_[mi].size() : 0);
        write(&pe, sizeof(pe));
    }

    // Calibration data (packed).
    write(calib_data.data(), calib_data.size());

    // Hash arrays (sorted uint64, concatenated).
    for (uint32_t mi = 0; mi < n_pool; ++mi) {
        if (mi < pool_.size() && !pool_[mi].empty())
            write(pool_[mi].data(), pool_[mi].size() * sizeof(uint64_t));
    }

    // Pre-merged pool (bac then arc): enables O(1) loading via mmap pointer.
    // Align to 8 bytes.
    {
        const uint64_t cur = static_cast<uint64_t>(f.tellp());
        const uint64_t pad = (merged_bac_off_u64 > cur) ? merged_bac_off_u64 - cur : 0;
        if (pad) { const std::vector<uint8_t> zeros(pad, 0); write(zeros.data(), pad); }
    }
    write(merged_bac_h.data(),  merged_bac_h.size()  * sizeof(uint64_t));
    write(merged_bac_id.data(), merged_bac_id.size());
    if (!merged_arc_h.empty()) {
        // Pad to 8-byte alignment between bac and arc sections.
        const uint64_t cur = static_cast<uint64_t>(f.tellp());
        const uint64_t pad = (merged_arc_off_u64 > cur) ? merged_arc_off_u64 - cur : 0;
        if (pad) { const std::vector<uint8_t> zeros(pad, 0); write(zeros.data(), pad); }
        write(merged_arc_h.data(),  merged_arc_h.size()  * sizeof(uint64_t));
        write(merged_arc_id.data(), merged_arc_id.size());
    }

    if (!f) throw std::runtime_error("markers: write failed: " + path.string());
}

// Append the pre-merged pool to an existing panel + patch the header offsets so
// MarkerReader takes its zero-copy mmap fast-path. No source genomes needed.
bool markers_add_premerged(const std::filesystem::path& path) {
    std::vector<uint64_t> bh, ah;
    std::vector<uint8_t>  bi, ai;
    uint64_t fsz = 0;
    {
        MarkerReader mr;
        mr.open(path);
        if (mr.has_premerged()) return false;          // already has it
        mr.build_merged_pool();                        // k-way merge into owned arrays
        const auto bhs = mr.merged_hashes_bac(); const auto bis = mr.merged_ids_bac();
        const auto ahs = mr.merged_hashes_arc(); const auto ais = mr.merged_ids_arc();
        bh.assign(bhs.begin(), bhs.end()); bi.assign(bis.begin(), bis.end());
        ah.assign(ahs.begin(), ahs.end()); ai.assign(ais.begin(), ais.end());
        fsz = mr.file_size();
    }

    const uint64_t mbac = (fsz + 7) & ~uint64_t{7};
    const uint64_t bac_bytes = bh.size() * sizeof(uint64_t) + bi.size();
    const uint64_t marc = ah.empty() ? 0 : ((mbac + bac_bytes + 7) & ~uint64_t{7});
    if (mbac > UINT32_MAX || marc > UINT32_MAX)
        throw std::runtime_error("markers remerge: panel too large for 32-bit offsets");

    std::fstream f(path, std::ios::in | std::ios::out | std::ios::binary);
    if (!f) throw std::runtime_error("markers remerge: cannot open " + path.string());
    auto wr = [&](const void* p, size_t n) { f.write(reinterpret_cast<const char*>(p),
                                                     static_cast<std::streamsize>(n)); };
    auto pad_to = [&](uint64_t target) {
        const uint64_t cur = static_cast<uint64_t>(f.tellp());
        if (target > cur) { std::vector<char> z(target - cur, 0); wr(z.data(), z.size()); }
    };
    f.seekp(0, std::ios::end);
    pad_to(mbac);
    wr(bh.data(), bh.size() * sizeof(uint64_t)); wr(bi.data(), bi.size());
    if (!ah.empty()) {
        pad_to(marc);
        wr(ah.data(), ah.size() * sizeof(uint64_t)); wr(ai.data(), ai.size());
    }
    const uint32_t mb = static_cast<uint32_t>(mbac), ma = static_cast<uint32_t>(marc);
    f.seekp(offsetof(MarkerHeader, merged_bac_off), std::ios::beg); wr(&mb, 4);
    f.seekp(offsetof(MarkerHeader, merged_arc_off), std::ios::beg); wr(&ma, 4);
    f.flush();
    if (!f) throw std::runtime_error("markers remerge: write failed: " + path.string());
    return true;
}

} // namespace genopack
