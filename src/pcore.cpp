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

// ── .ptier side-channel binary format ─────────────────────────────────────────
// magic  = 0x5054494552000001ULL  ("PTIER\0\0\1")
// header: { uint64 magic; uint64 n_pairs; uint64 fmh_seed; uint64 frac_max_hash;
//           uint32 k; uint32 _pad; }   — 48 bytes total
// body:   { uint64 aamer_hash; uint64 genus_hash; }[] — n_pairs records, unsorted

struct PtierFileHeader {
    uint64_t magic         = 0x5054494552000001ULL;
    uint64_t n_pairs       = 0;
    uint64_t fmh_seed      = 0;
    uint64_t frac_max_hash = 0;
    uint32_t k             = 0;
    uint32_t _pad0         = 0;
    uint64_t _pad1         = 0;   // reserved — pad to 48 bytes
};
static_assert(sizeof(PtierFileHeader) == 48);

PcoreWriter::PcoreWriter(uint32_t k, uint32_t min_seg_aa, float theta, std::string spill_dir)
    : k_(k), min_seg_aa_(min_seg_aa), theta_(theta), spill_dir_(std::move(spill_dir)) {}

PcoreWriter::~PcoreWriter() {
    if (spill_) { std::fclose(spill_); std::remove(spill_path_.c_str()); }
    if (tier_sc_) { std::fclose(tier_sc_); tier_sc_ = nullptr; }
}

void PcoreWriter::set_tier_sidechannel(const std::string& path) {
    if (path.empty()) return;
    tier_sc_path_ = path;
    tier_sc_ = std::fopen(path.c_str(), "w+b");
    if (!tier_sc_) throw std::runtime_error("PcoreWriter: cannot open tier sidechannel " + path);
    // Write placeholder header — n_pairs and frac_max_hash filled in finalize().
    PtierFileHeader hdr{};
    hdr.fmh_seed = fmh_seed_;
    hdr.k        = k_;
    std::fwrite(&hdr, sizeof(hdr), 1, tier_sc_);
    tier_sc_n_pairs_ = 0;
}

void PcoreWriter::open_spill_() {
    if (spill_) return;
    std::filesystem::path dir;
    if (!spill_dir_.empty()) {
        dir = std::filesystem::path(spill_dir_);
    } else if (const char* env = std::getenv("GENOPACK_SKETCH_SPILL_DIR"); env && *env) {
        dir = std::filesystem::path(env);
    } else if (const char* env = std::getenv("GENOPACK_SPILL_DIR"); env && *env) {
        dir = std::filesystem::path(env);
    } else {
        dir = std::filesystem::temp_directory_path();
    }
    std::error_code ec; std::filesystem::create_directories(dir, ec);
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
        for (uint64_t h : q) ++cnt[h];                 // member sets are sorted-unique → counts members
    add_from_counts(key_hash, cnt, n_mem);
}

void PcoreWriter::add_from_counts(uint64_t key_hash,
                                  const std::unordered_map<uint64_t, uint32_t>& cnt,
                                  uint32_t n_mem) {
    if (n_mem == 0 || cnt.empty()) return;
    std::vector<uint64_t> aamers; aamers.reserve(cnt.size());
    for (const auto& [h, c] : cnt) { (void)c; aamers.push_back(h); }
    std::sort(aamers.begin(), aamers.end());
    std::vector<uint32_t> counts(aamers.size());
    for (size_t i = 0; i < aamers.size(); ++i) counts[i] = cnt.at(aamers[i]);
    add_sorted(key_hash, aamers, counts, n_mem);
}

void PcoreWriter::add_sorted(uint64_t key_hash,
                             const std::vector<uint64_t>& aamers,
                             const std::vector<uint32_t>& counts,
                             uint32_t n_mem) {
    if (n_mem == 0 || aamers.empty()) return;

    const uint32_t C = pcore_core_threshold(n_mem, theta_);   // ⌈θ·n⌉

    // aamers is globally sorted ascending → each run is a sorted subsequence, so no
    // per-run re-sort. Stratify by TRUE count: core (c≥C) | singleton (c==1) | multi.
    std::vector<uint64_t> single, multi_h, core_h;
    std::vector<uint8_t>  multi_q, core_q;
    for (size_t i = 0; i < aamers.size(); ++i) {
        const uint64_t h = aamers[i]; const uint32_t c = counts[i];
        if (c >= C)      { core_h.push_back(h);  core_q.push_back(pcore_quant_prev(c, n_mem)); }
        else if (c <= 1) { single.push_back(h); }
        else             { multi_h.push_back(h); multi_q.push_back(pcore_quant_prev(c, n_mem)); }
    }

    const std::vector<uint8_t> es = pfor::encode_sorted(single);
    const std::vector<uint8_t> em = pfor::encode_sorted(multi_h);
    const std::vector<uint8_t> ec = pfor::encode_sorted(core_h);

    open_spill_();
    static const char zeros[8] = {0};
    while (pool_cursor_ & 7u) { std::fwrite(zeros, 1, 1, spill_); ++pool_cursor_; }
    const uint64_t run_off = pool_cursor_;
    auto spill = [&](const void* p, size_t n) { if (n) std::fwrite(p, 1, n, spill_); pool_cursor_ += n; };
    spill(es.data(), es.size());
    spill(em.data(), em.size());
    spill(ec.data(), ec.size());
    spill(multi_q.data(), multi_q.size());
    spill(core_q.data(),  core_q.size());

    // Emit (aamer_hash, genus_hash) pairs to side-channel if open.
    if (tier_sc_) {
        struct PtierPair { uint64_t aamer; uint64_t genus; };
        for (uint64_t h : aamers) {
            PtierPair p{h, key_hash};
            std::fwrite(&p, sizeof(p), 1, tier_sc_);
        }
        tier_sc_n_pairs_ += static_cast<uint64_t>(aamers.size());
    }

    EntryMeta m{};
    m.key_hash      = key_hash;
    m.run_off       = run_off;
    m.n_singleton   = static_cast<uint32_t>(single.size());
    m.n_multi       = static_cast<uint32_t>(multi_h.size());
    m.n_core        = static_cast<uint32_t>(core_h.size());
    m.n_members     = n_mem;
    m.enc_singleton = static_cast<uint32_t>(es.size());
    m.enc_multi     = static_cast<uint32_t>(em.size());
    m.enc_core      = static_cast<uint32_t>(ec.size());
    meta_.push_back(m);

    // Order-independent incremental model hash (XOR of per-entry digests).
    XXH3_state_t* st = XXH3_createState();
    XXH3_64bits_reset(st);
    XXH3_64bits_update(st, &key_hash, sizeof(key_hash));
    const uint32_t cnts[4] = { m.n_singleton, m.n_multi, m.n_core, m.n_members };
    XXH3_64bits_update(st, cnts, sizeof(cnts));
    if (!es.empty()) XXH3_64bits_update(st, es.data(), es.size());
    if (!em.empty()) XXH3_64bits_update(st, em.data(), em.size());
    if (!ec.empty()) XXH3_64bits_update(st, ec.data(), ec.size());
    if (!multi_q.empty()) XXH3_64bits_update(st, multi_q.data(), multi_q.size());
    if (!core_q.empty())  XXH3_64bits_update(st, core_q.data(),  core_q.size());
    hash_fold_ ^= XXH3_64bits_digest(st);
    XXH3_freeState(st);
}

uint64_t PcoreWriter::model_hash() const {
    uint32_t theta_bits;
    std::memcpy(&theta_bits, &theta_, 4);
    uint64_t buf[4] = { (static_cast<uint64_t>(k_) << 32) | min_seg_aa_, theta_bits,
                        fmh_seed_, taxonomy_hash_ };
    return XXH3_64bits(buf, sizeof(buf)) ^ hash_fold_;
}

SectionDesc PcoreWriter::finalize(AppendWriter& w, uint64_t section_id, uint64_t frac_max_hash) {
    // Finalise the .ptier side-channel: rewrite the header with actual n_pairs + frac_max.
    if (tier_sc_) {
        std::fflush(tier_sc_);
        std::rewind(tier_sc_);
        PtierFileHeader hdr{};
        hdr.fmh_seed      = fmh_seed_;
        hdr.n_pairs       = tier_sc_n_pairs_;
        hdr.frac_max_hash = frac_max_hash;
        hdr.k             = k_;
        std::fwrite(&hdr, sizeof(hdr), 1, tier_sc_);
        std::fclose(tier_sc_);
        tier_sc_ = nullptr;
    }
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

    const uint64_t entries_off = sizeof(PcoreHeader) + sizeof(PcoreHeaderExtV1);   // 128
    const uint64_t buckets_off = entries_off + static_cast<uint64_t>(n) * sizeof(PcoreEntryV1);
    uint64_t pool_off = buckets_off + static_cast<uint64_t>(n_buckets) * sizeof(uint32_t);
    pool_off = (pool_off + 7) & ~uint64_t{7};

    std::vector<PcoreEntryV1> ents(n);
    for (uint32_t i = 0; i < n; ++i) {
        const auto& m = meta_[i];
        ents[i].key_hash     = m.key_hash;
        ents[i].run_offset   = pool_off + m.run_off;
        ents[i].n_singleton  = m.n_singleton;
        ents[i].n_multi      = m.n_multi;
        ents[i].n_core       = m.n_core;
        ents[i].n_members    = m.n_members;
        ents[i].enc_singleton= m.enc_singleton;
        ents[i].enc_multi    = m.enc_multi;
        ents[i].enc_core     = m.enc_core;
        ents[i]._pad         = 0;
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
    hdr.codec          = PCORE_CODEC_V1;
    hdr.prev_quant     = PCORE_PREVQ_LOGN;
    hdr.frame_table_id = frame_table_id_;
    hdr.header_bytes   = 128;

    PcoreHeaderExtV1 ext{};
    ext.fmh_seed      = fmh_seed_;
    ext.taxonomy_hash = taxonomy_hash_;
    ext.k_frame_pin   = (static_cast<uint64_t>(k_) << 8) | frame_table_id_;

    w.append(&hdr, sizeof(hdr));
    w.append(&ext, sizeof(ext));
    if (n) w.append(ents.data(), static_cast<uint64_t>(n) * sizeof(PcoreEntryV1));
    w.append(buckets.data(), static_cast<uint64_t>(n_buckets) * sizeof(uint32_t));
    w.align(8);

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
    sd.version           = 2;                          // v1 stratified-PFOR codec
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
