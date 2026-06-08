#include "run_tier.hpp"
#include <genopack/run_tier.hpp>
#include <genopack/pcore.hpp>
#include <genopack/archive.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/toc.hpp>
#include <genopack/format.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace genopack {

// ── .ptier on-disk layout (written by PcoreWriter::set_tier_sidechannel) ─────
// header: PtierFileHeader (48 bytes)
// body:   { uint64 aamer_hash; uint64 genus_hash; }[]  n_pairs records, unsorted

struct PtierFileHeader {
    uint64_t magic         = 0x5054494552000001ULL;
    uint64_t n_pairs       = 0;
    uint64_t fmh_seed      = 0;
    uint64_t frac_max_hash = 0;
    uint32_t k             = 0;
    uint32_t _pad0         = 0;
    uint64_t _pad1         = 0;
};
static_assert(sizeof(PtierFileHeader) == 48);

struct PtierPair {
    uint64_t aamer;
    uint64_t genus;
};

// ── .gtier on-disk layout ─────────────────────────────────────────────────────
// header: GtierFileHeader (48 bytes)
// body:   { uint64 aamer_hash; uint8 tier; }[] — n_aamers records, sorted by aamer_hash

static constexpr uint64_t GTIER_MAGIC = 0x4754494552000001ULL;  // "GTIER\0\0\1"

struct GtierFileHeader {
    uint64_t magic          = GTIER_MAGIC;
    uint64_t n_aamers       = 0;
    uint64_t total_genera_G = 0;
    uint64_t fmh_seed       = 0;
    uint64_t frac_max_hash  = 0;
    uint32_t k              = 0;
    uint32_t _pad           = 0;
};
static_assert(sizeof(GtierFileHeader) == 48);

struct GtierEntry {
    uint64_t aamer_hash;
    uint8_t  tier;
} __attribute__((packed));
static_assert(sizeof(GtierEntry) == 9);

// ── cmd_tier_merge ────────────────────────────────────────────────────────────

int cmd_tier_merge(const std::vector<std::string>& ptier_files,
                   const std::string& out_path) {
    if (ptier_files.empty()) {
        spdlog::error("tier merge: no input .ptier files");
        return 1;
    }

    // 1. Read all headers and verify consistency.
    uint64_t ref_fmh_seed = 0, ref_frac_max = 0;
    uint32_t ref_k = 0;
    bool     ref_set = false;
    std::vector<PtierFileHeader> hdrs(ptier_files.size());
    for (size_t fi = 0; fi < ptier_files.size(); ++fi) {
        FILE* f = std::fopen(ptier_files[fi].c_str(), "rb");
        if (!f) { spdlog::error("tier merge: cannot open {}", ptier_files[fi]); return 1; }
        if (std::fread(&hdrs[fi], sizeof(PtierFileHeader), 1, f) != 1) {
            std::fclose(f); spdlog::error("tier merge: truncated header in {}", ptier_files[fi]); return 1;
        }
        std::fclose(f);
        if (hdrs[fi].magic != 0x5054494552000001ULL) {
            spdlog::error("tier merge: bad magic in {}", ptier_files[fi]); return 1;
        }
        if (!ref_set) {
            ref_fmh_seed = hdrs[fi].fmh_seed;
            ref_frac_max = hdrs[fi].frac_max_hash;
            ref_k        = hdrs[fi].k;
            ref_set      = true;
        } else {
            if (hdrs[fi].fmh_seed != ref_fmh_seed || hdrs[fi].frac_max_hash != ref_frac_max || hdrs[fi].k != ref_k) {
                spdlog::error("tier merge: mismatched fmh_seed/frac_max/k between {} and first file", ptier_files[fi]);
                return 1;
            }
        }
    }

    // 2. Single pass: for each aamer, accumulate distinct genus hashes.
    //    aamer_hash → set<genus_hash>  (use unordered_set<uint64_t> per aamer)
    //    For large corpora this fits in RAM: a 9M-genus build with frac=100 has
    //    ~600M aamers → ~4.8 GB at 8 bytes/aamer key. Acceptable.
    std::unordered_map<uint64_t, std::unordered_set<uint64_t>> aamer_to_genera;
    aamer_to_genera.reserve(64'000'000ULL);

    uint64_t total_pairs = 0;
    for (size_t fi = 0; fi < ptier_files.size(); ++fi) {
        FILE* f = std::fopen(ptier_files[fi].c_str(), "rb");
        if (!f) { spdlog::error("tier merge: cannot open {} for reading", ptier_files[fi]); return 1; }
        std::fseek(f, static_cast<long>(sizeof(PtierFileHeader)), SEEK_SET);
        const uint64_t n = hdrs[fi].n_pairs;
        constexpr size_t kBuf = 1u << 16;
        std::vector<PtierPair> buf(kBuf);
        uint64_t remaining = n;
        while (remaining > 0) {
            const size_t batch = std::min(remaining, static_cast<uint64_t>(kBuf));
            if (std::fread(buf.data(), sizeof(PtierPair), batch, f) != batch) {
                std::fclose(f);
                spdlog::error("tier merge: short read in {}", ptier_files[fi]);
                return 1;
            }
            for (size_t i = 0; i < batch; ++i)
                aamer_to_genera[buf[i].aamer].insert(buf[i].genus);
            remaining -= batch;
            total_pairs += batch;
        }
        std::fclose(f);
        spdlog::info("tier merge: read {} pairs from {}", n, ptier_files[fi]);
    }

    // 3. Count total distinct genera across all aamers' genus sets → G.
    std::unordered_set<uint64_t> all_genera;
    for (const auto& [a, gs] : aamer_to_genera)
        for (uint64_t g : gs) all_genera.insert(g);
    const uint64_t G = static_cast<uint64_t>(all_genera.size());
    all_genera.clear();  // free RAM

    spdlog::info("tier merge: {} unique aamers, {} total genera (G={})", aamer_to_genera.size(), G, G);

    // 4. Build sorted (aamer_hash, tier) table.
    std::vector<GtierEntry> entries;
    entries.reserve(aamer_to_genera.size());
    for (const auto& [aamer, genera] : aamer_to_genera) {
        GtierEntry e;
        e.aamer_hash = aamer;
        e.tier       = pcore_tier_from_count(static_cast<uint32_t>(
                           std::min(genera.size(), static_cast<size_t>(UINT32_MAX))));
        entries.push_back(e);
    }
    aamer_to_genera.clear();
    std::sort(entries.begin(), entries.end(),
              [](const GtierEntry& a, const GtierEntry& b) { return a.aamer_hash < b.aamer_hash; });

    // 5. Write global .gtier file.
    FILE* out = std::fopen(out_path.c_str(), "wb");
    if (!out) { spdlog::error("tier merge: cannot write {}", out_path); return 1; }

    GtierFileHeader ghdr{};
    ghdr.magic          = GTIER_MAGIC;
    ghdr.n_aamers       = static_cast<uint64_t>(entries.size());
    ghdr.total_genera_G = G;
    ghdr.fmh_seed       = ref_fmh_seed;
    ghdr.frac_max_hash  = ref_frac_max;
    ghdr.k              = ref_k;
    std::fwrite(&ghdr, sizeof(ghdr), 1, out);

    // Write in chunks to avoid one giant syscall.
    constexpr size_t kChunk = 1u << 16;
    for (size_t off = 0; off < entries.size(); off += kChunk) {
        const size_t n = std::min(kChunk, entries.size() - off);
        std::fwrite(entries.data() + off, sizeof(GtierEntry), n, out);
    }
    std::fclose(out);

    spdlog::info("tier merge: wrote {} aamers, G={} → {}", entries.size(), G, out_path);
    return 0;
}

// ── GtierView — binary search over a mmap'd .gtier body ──────────────────────

struct GtierView {
    const GtierEntry*  data    = nullptr;
    uint64_t           n       = 0;
    uint64_t           G       = 0;

    bool valid() const noexcept { return data != nullptr && n > 0; }

    // Returns tier byte for aamer_hash, or 0xFF if not found.
    uint8_t lookup(uint64_t aamer_hash) const noexcept {
        if (!valid()) return 0xFF;
        // Binary search over packed 9-byte entries — use index arithmetic.
        uint64_t lo = 0, hi = n;
        while (lo < hi) {
            const uint64_t mid = lo + (hi - lo) / 2;
            const uint64_t h   = data[mid].aamer_hash;
            if (h == aamer_hash) return data[mid].tier;
            if (h  < aamer_hash) lo = mid + 1;
            else                 hi = mid;
        }
        return 0xFF;  // absent
    }
};

// ── cmd_tier_stamp ────────────────────────────────────────────────────────────

int cmd_tier_stamp(const std::string& gpk_path,
                   const std::string& gtier_path,
                   const std::string& out_path_arg) {
    const std::string out_path = out_path_arg.empty()
        ? (gpk_path.substr(0, gpk_path.size() - (gpk_path.size() > 4 ? 4 : 0)) + ".tier.gpk")
        : out_path_arg;

    // 1. Mmap .gtier and build the lookup view.
    MmapFileReader gtier_mmap;
    gtier_mmap.open(gtier_path);
    if (gtier_mmap.size() < sizeof(GtierFileHeader))
        throw std::runtime_error("tier stamp: .gtier file too small");
    const auto* ghdr = reinterpret_cast<const GtierFileHeader*>(gtier_mmap.data());
    if (ghdr->magic != GTIER_MAGIC)
        throw std::runtime_error("tier stamp: bad .gtier magic");
    GtierView tv;
    tv.data = reinterpret_cast<const GtierEntry*>(gtier_mmap.data() + sizeof(GtierFileHeader));
    tv.n    = ghdr->n_aamers;
    tv.G    = ghdr->total_genera_G;
    spdlog::info("tier stamp: .gtier loaded: {} aamers, G={}", tv.n, tv.G);

    // 2. Open the source .gpk and locate its PCORE section.
    MmapFileReader src_mmap;
    src_mmap.open(gpk_path);
    Toc src_toc = TocReader::read(src_mmap);

    const SectionDesc* pcore_sd = nullptr;
    {
        uint64_t best_id = 0;
        for (auto* sd : src_toc.find_by_type(SEC_PCORE))
            if (sd->section_id > best_id) { best_id = sd->section_id; pcore_sd = sd; }
    }
    if (!pcore_sd) { spdlog::error("tier stamp: no PCORE section in {}", gpk_path); return 1; }

    // Open the existing PCORE reader to iterate entries.
    PcoreReader pr;
    pr.open(src_mmap.data(), pcore_sd->file_offset, pcore_sd->compressed_size);
    if (pr.codec() != PCORE_CODEC_V1 && pr.codec() != PCORE_CODEC_V2) {
        spdlog::error("tier stamp: PCORE codec {} not V1/V2 — cannot stamp", pr.codec());
        return 1;
    }
    const uint32_t n_entries = pr.n_entries();
    spdlog::info("tier stamp: PCORE has {} entries (codec {})", n_entries, pr.codec());

    // 3. Build a new .gpk copying all sections verbatim except PCORE (rewritten as V2).
    const std::string partial = out_path + ".partial";
    std::filesystem::create_directories(
        std::filesystem::path(out_path).parent_path().empty()
        ? std::filesystem::path(".")
        : std::filesystem::path(out_path).parent_path());

    AppendWriter w;
    w.create(partial);

    // Copy FileHeader verbatim (new UUID/generation handled by TocWriter).
    w.append(src_mmap.data(), sizeof(FileHeader));

    TocWriter tw;
    (void)src_toc.next_section_id();  // unused here — TOC writer assigns its own ids

    // Copy all non-PCORE sections verbatim.
    for (const auto& sd : src_toc.sections) {
        if (sd.type == SEC_PCORE) continue;
        if (sd.compressed_size == 0) { tw.add_section(sd); continue; }
        w.align(8);
        SectionDesc nd = sd;
        nd.file_offset = w.current_offset();
        w.append(src_mmap.data() + sd.file_offset, sd.compressed_size);
        tw.add_section(nd);
    }

    // 4. Build new V2 PCORE section: copy V1 streams + append tier bytes per entry.
    //    Layout per entry in spill: enc_s | enc_m | enc_c | prevq_m | prevq_c | tier_s | tier_m | tier_c
    //    We build a new PcoreWriter-like byte stream by hand to avoid re-encoding.
    w.align(8);
    const uint64_t pcore_start = w.current_offset();

    // Read the original PCORE header+ext verbatim (128 bytes).
    // We'll patch codec and total_genera_G then write.
    // Copy raw header bytes first so layout matches.
    std::vector<uint8_t> hdr_bytes(128);
    std::memcpy(hdr_bytes.data(), src_mmap.data() + pcore_sd->file_offset, 128);

    auto* new_hdr = reinterpret_cast<PcoreHeader*>(hdr_bytes.data());
    auto* new_ext = reinterpret_cast<PcoreHeaderExtV1*>(hdr_bytes.data() + sizeof(PcoreHeader));
    new_hdr->codec       = PCORE_CODEC_V2;
    new_ext->total_genera_G  = tv.G;

    // We need to recompute all entry run_offsets because tier bytes expand each entry.
    // Strategy: collect (new_run_off, entry_bytes) pairs, write header+entries+buckets, then pool.
    struct NewEntry {
        uint64_t key_hash;
        uint64_t new_run_off;
        uint32_t n_singleton, n_multi, n_core, n_members;
        uint32_t enc_singleton, enc_multi, enc_core;
        std::vector<uint8_t> tier_s, tier_m, tier_c;
    };

    std::vector<NewEntry> new_entries;
    new_entries.reserve(n_entries);

    // Spill the new run pool to a temp buffer.
    std::vector<uint8_t> pool;

    std::vector<uint64_t> aamers_tmp;
    std::vector<float>    prevf_tmp;

    uint64_t n_stamped = 0, n_absent = 0;

    for (uint32_t i = 0; i < n_entries; ++i) {
        const PcoreView v = pr.lookup(pr.key_hash_at(i));
        if (!v.valid() || !v.is_v1) continue;

        NewEntry ne{};
        ne.key_hash     = v.key_hash;
        ne.n_singleton  = v.n_singleton;
        ne.n_multi      = v.n_multi;
        ne.n_core       = v.n_core;
        ne.n_members    = v.n_members;
        ne.enc_singleton = v.enc_s_bytes;
        ne.enc_multi     = v.enc_m_bytes;
        ne.enc_core      = v.enc_c_bytes;

        // Decode aamer lists to look up tiers.
        pfor::decode_into(v.enc_s, v.enc_s_bytes, v.n_singleton, aamers_tmp);
        ne.tier_s.resize(v.n_singleton);
        for (uint32_t j = 0; j < v.n_singleton; ++j) {
            const uint8_t t = tv.lookup(aamers_tmp[j]);
            ne.tier_s[j] = (t == 0xFF) ? 0 : t;
            if (t == 0xFF) ++n_absent; else ++n_stamped;
        }
        pfor::decode_into(v.enc_m, v.enc_m_bytes, v.n_multi, aamers_tmp);
        ne.tier_m.resize(v.n_multi);
        for (uint32_t j = 0; j < v.n_multi; ++j) {
            const uint8_t t = tv.lookup(aamers_tmp[j]);
            ne.tier_m[j] = (t == 0xFF) ? 0 : t;
            if (t == 0xFF) ++n_absent; else ++n_stamped;
        }
        pfor::decode_into(v.enc_c, v.enc_c_bytes, v.n_core, aamers_tmp);
        ne.tier_c.resize(v.n_core);
        for (uint32_t j = 0; j < v.n_core; ++j) {
            const uint8_t t = tv.lookup(aamers_tmp[j]);
            ne.tier_c[j] = (t == 0xFF) ? 0 : t;
            if (t == 0xFF) ++n_absent; else ++n_stamped;
        }

        // Align to 8 bytes, record pool offset.
        while (pool.size() & 7u) pool.push_back(0);
        ne.new_run_off = static_cast<uint64_t>(pool.size());

        // Append: enc_s | enc_m | enc_c | prevq_m | prevq_c | tier_s | tier_m | tier_c
        pool.insert(pool.end(), v.enc_s,   v.enc_s   + v.enc_s_bytes);
        pool.insert(pool.end(), v.enc_m,   v.enc_m   + v.enc_m_bytes);
        pool.insert(pool.end(), v.enc_c,   v.enc_c   + v.enc_c_bytes);
        pool.insert(pool.end(), v.prevq_m, v.prevq_m + v.n_multi);
        pool.insert(pool.end(), v.prevq_c, v.prevq_c + v.n_core);
        pool.insert(pool.end(), ne.tier_s.begin(), ne.tier_s.end());
        pool.insert(pool.end(), ne.tier_m.begin(), ne.tier_m.end());
        pool.insert(pool.end(), ne.tier_c.begin(), ne.tier_c.end());

        new_entries.push_back(std::move(ne));
    }
    spdlog::info("tier stamp: {} aamers stamped, {} absent in gtier", n_stamped, n_absent);

    // Rebuild hash buckets.
    const uint32_t n = static_cast<uint32_t>(new_entries.size());
    const uint32_t min_bkts = (n == 0) ? 1 : static_cast<uint32_t>(static_cast<double>(n) / 0.7 + 1.0);
    uint32_t n_bkts = 1; while (n_bkts < min_bkts) n_bkts <<= 1;
    const uint32_t mask = n_bkts - 1;
    std::vector<uint32_t> buckets(n_bkts, PCORE_EMPTY);
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t slot = static_cast<uint32_t>(new_entries[i].key_hash) & mask;
        while (buckets[slot] != PCORE_EMPTY) slot = (slot + 1) & mask;
        buckets[slot] = i;
    }

    // Compute section-relative offsets.
    const uint64_t ent_off  = 128;  // after hdr+ext
    const uint64_t bkt_off  = ent_off + static_cast<uint64_t>(n) * sizeof(PcoreEntryV1);
    uint64_t pool_off = bkt_off + static_cast<uint64_t>(n_bkts) * sizeof(uint32_t);
    pool_off = (pool_off + 7) & ~uint64_t{7};

    // Patch header fields.
    new_hdr->n_entries      = n;
    new_hdr->n_buckets      = n_bkts;
    new_hdr->entries_offset = ent_off;
    new_hdr->buckets_offset = bkt_off;

    // Build entry array with updated run_offset.
    std::vector<PcoreEntryV1> ents(n);
    for (uint32_t i = 0; i < n; ++i) {
        const auto& ne = new_entries[i];
        ents[i].key_hash      = ne.key_hash;
        ents[i].run_offset    = pool_off + ne.new_run_off;
        ents[i].n_singleton   = ne.n_singleton;
        ents[i].n_multi       = ne.n_multi;
        ents[i].n_core        = ne.n_core;
        ents[i].n_members     = ne.n_members;
        ents[i].enc_singleton = ne.enc_singleton;
        ents[i].enc_multi     = ne.enc_multi;
        ents[i].enc_core      = ne.enc_core;
        ents[i]._pad          = 0;
    }

    // Write PCORE section.
    w.append(hdr_bytes.data(), 128);
    if (n) w.append(ents.data(), static_cast<uint64_t>(n) * sizeof(PcoreEntryV1));
    w.append(buckets.data(), static_cast<uint64_t>(n_bkts) * sizeof(uint32_t));
    w.align(8);
    if (!pool.empty()) w.append(pool.data(), pool.size());

    const uint64_t pcore_end = w.current_offset();

    SectionDesc pcore_nd = *pcore_sd;
    pcore_nd.file_offset       = pcore_start;
    pcore_nd.compressed_size   = pcore_end - pcore_start;
    pcore_nd.uncompressed_size = pcore_end - pcore_start;
    pcore_nd.version           = 3;  // V2 codec
    tw.add_section(pcore_nd);

    // 5. Write TOC.
    const auto& src_toc_hdr = src_toc.header;
    tw.finalize(w,
                src_toc_hdr.generation + 1,
                src_toc_hdr.live_genome_count,
                src_toc_hdr.total_genome_count,
                /*prev_toc_offset=*/0,
                src_toc_hdr.catalog_root_section_id,
                src_toc_hdr.accession_root_section_id,
                src_toc_hdr.tombstone_root_section_id);
    w.flush();
    w.close();

    std::filesystem::rename(partial, out_path);
    spdlog::info("tier stamp: wrote {}", out_path);
    return 0;
}

} // namespace genopack
