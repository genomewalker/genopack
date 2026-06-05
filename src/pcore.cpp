#include <genopack/pcore.hpp>
#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <unordered_map>
#include <vector>
#include <xxhash.h>
#include <unistd.h>

namespace genopack {

PcoreWriter::PcoreWriter(uint32_t k, uint32_t min_seg_aa, float theta)
    : k_(k), min_seg_aa_(min_seg_aa), theta_(theta) {}

PcoreWriter::~PcoreWriter() {
    if (spill_) { std::fclose(spill_); std::remove(spill_path_.c_str()); }
}

void PcoreWriter::open_spill_() {
    if (spill_) return;
    const char* env = std::getenv("GENOPACK_SPILL_DIR");
    std::filesystem::path dir = (env && *env) ? std::filesystem::path(env)
                                              : std::filesystem::temp_directory_path();
    static std::atomic<uint64_t> ctr{0};
    spill_path_ = (dir / ("pcore_spill_" + std::to_string(::getpid()) + "_" +
                          std::to_string(ctr.fetch_add(1)) + ".tmp")).string();
    spill_ = std::fopen(spill_path_.c_str(), "w+b");
    if (!spill_) throw std::runtime_error("PcoreWriter: cannot open spill file " + spill_path_);
}

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

    std::vector<uint64_t> aamers;
    aamers.reserve(cnt.size());
    for (const auto& [h, _] : cnt) aamers.push_back(h);
    if (aamers.empty()) return;
    std::sort(aamers.begin(), aamers.end());
    std::vector<uint8_t> prev(aamers.size());
    for (size_t i = 0; i < aamers.size(); ++i) {
        const uint32_t c = cnt[aamers[i]];
        prev[i] = static_cast<uint8_t>(c > 255u ? 255u : c);   // member count (cap 255)
    }
    const uint32_t na = static_cast<uint32_t>(aamers.size());

    // Spill this genus's pool (8-aligned aamers, then prev). Only metadata stays in RAM.
    open_spill_();
    static const char zeros[8] = {0};
    while (pool_cursor_ & 7u) { std::fwrite(zeros, 1, 1, spill_); ++pool_cursor_; }
    const uint64_t aoff = pool_cursor_;
    std::fwrite(aamers.data(), sizeof(uint64_t), na, spill_); pool_cursor_ += static_cast<uint64_t>(na) * 8;
    const uint64_t poff = pool_cursor_;
    std::fwrite(prev.data(), 1, na, spill_);                  pool_cursor_ += na;
    meta_.push_back({key_hash, aoff, poff, na, n_mem});

    // Order-independent incremental model hash: XOR of per-entry digests.
    XXH3_state_t* st = XXH3_createState();
    XXH3_64bits_reset(st);
    XXH3_64bits_update(st, &key_hash, sizeof(key_hash));
    XXH3_64bits_update(st, &na, sizeof(na));
    XXH3_64bits_update(st, aamers.data(), static_cast<size_t>(na) * sizeof(uint64_t));
    XXH3_64bits_update(st, prev.data(), na);
    hash_fold_ ^= XXH3_64bits_digest(st);
    XXH3_freeState(st);
}

uint64_t PcoreWriter::model_hash() const {
    uint32_t theta_bits;
    std::memcpy(&theta_bits, &theta_, 4);
    uint64_t buf[3] = { (static_cast<uint64_t>(k_) << 32) | min_seg_aa_, theta_bits, 0 };
    return XXH3_64bits(buf, sizeof(buf)) ^ hash_fold_;
}

SectionDesc PcoreWriter::finalize(AppendWriter& w, uint64_t section_id, uint64_t frac_max_hash) {
    if (spill_) std::fflush(spill_);
    const auto n = static_cast<uint32_t>(meta_.size());
    const uint32_t min_buckets = (n == 0) ? 1
        : static_cast<uint32_t>(static_cast<double>(n) / 0.7 + 1.0);
    const uint32_t n_buckets = next_pow2(min_buckets);
    const uint32_t mask      = n_buckets - 1;

    std::vector<uint32_t> buckets(n_buckets, PCORE_EMPTY);
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t slot = static_cast<uint32_t>(meta_[i].key_hash) & mask;
        while (buckets[slot] != PCORE_EMPTY) slot = (slot + 1) & mask;
        buckets[slot] = i;
    }

    const uint64_t entries_off = sizeof(PcoreHeader);
    const uint64_t buckets_off = entries_off + static_cast<uint64_t>(n) * sizeof(PcoreEntry);
    uint64_t pool_off = buckets_off + static_cast<uint64_t>(n_buckets) * sizeof(uint32_t);
    pool_off = (pool_off + 7) & ~uint64_t{7};

    // Section-relative offsets = pool_off + spilled-pool-relative offset (both 8-aligned).
    std::vector<PcoreEntry> ents(n);
    for (uint32_t i = 0; i < n; ++i) {
        ents[i].key_hash      = meta_[i].key_hash;
        ents[i].aamers_offset = pool_off + meta_[i].aamers_off;
        ents[i].prev_offset   = pool_off + meta_[i].prev_off;
        ents[i].n_aamers      = meta_[i].n_aamers;
        ents[i].n_members     = meta_[i].n_members;
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
    w.align(8);   // reach pool_off (section_start is 8-aligned)

    // Stream the spilled pool verbatim into the section (bounded memory).
    if (spill_ && pool_cursor_ > 0) {
        std::rewind(spill_);
        std::vector<char> buf(1u << 20);
        size_t r;
        while ((r = std::fread(buf.data(), 1, buf.size(), spill_)) > 0)
            w.append(buf.data(), r);
        std::fclose(spill_);
        std::remove(spill_path_.c_str());
        spill_ = nullptr;
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
