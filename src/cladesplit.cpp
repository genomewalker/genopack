#include "genopack/cladesplit.hpp"
#include "genopack/aamer.hpp"
#include "genopack/util.hpp"
#include "genopack/types.hpp"
#include "genopack/archive.hpp"
#include "check/pass_b.hpp"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <fstream>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>

namespace genopack {

namespace {

struct CspHeader {
    uint64_t magic;
    uint64_t n_aamers;
    uint64_t n_clades;
    uint64_t names_off;
    uint16_t frac_c;     // build config — scoring must reproduce these exactly
    uint16_t min_aa;
    uint16_t mode;
    uint16_t _reserved;
    uint64_t scc_off;    // v3 only: byte offset of the per-genus SCC section (0 in v2; first 40 bytes match v2)
};
static_assert(sizeof(CspHeader) == 48, "csp v3 header must be 48 bytes (first 40 match v2)");
constexpr size_t kCspHdrV2 = 40;   // v2 on-disk header size (no scc_off)
constexpr size_t kCspHdrV3 = 48;

std::string load_fasta(const std::filesystem::path& p) {
    if (p.extension() == ".gz") return decompress_gz(p);
    std::ifstream f(p, std::ios::binary);
    if (!f) throw std::runtime_error("cladesplit: cannot read " + p.string());
    return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

std::string extract_rank(std::string_view tax, const char* pre) {
    auto p = tax.find(pre);
    if (p == std::string_view::npos) return {};
    auto e = tax.find(';', p);
    std::string s(tax.substr(p, e == std::string_view::npos ? tax.size() - p : e - p));
    return s == std::string(pre) ? std::string{} : s;
}

std::string taxonomy_of(const BuildRecord& r) {
    for (const auto& [k, v] : r.extra_fields)
        if (k == "taxonomy" || k == "gtdb_taxonomy") return v;
    return {};
}

inline uint64_t splitmix(uint64_t v) {
    v ^= v >> 30; v *= 0xbf58476d1ce4e5b9ULL;
    v ^= v >> 27; v *= 0x94d049bb133111ebULL;
    v ^= v >> 31; return v;
}
inline uint64_t aa4_hash(const uint8_t* a) {  // hash one 4-AA strobe
    uint64_t v = (uint64_t)a[0] | ((uint64_t)a[1] << 5) | ((uint64_t)a[2] << 10) | ((uint64_t)a[3] << 15);
    v ^= v >> 17; v *= 0x9E3779B97F4A7C15ULL; v ^= v >> 29; return v;
}

// Emit aamers from one inter-stop segment: aa[i] = amino-acid 0-19, wob[i] = 3rd
// codon-position base 0-3. Mode selects the primitive.
// Emit (FracMinHash-subsampled, 1/c) aamers from one segment. Subsampling is INLINE
// so the caller never materializes the full ~30× larger raw stream.
void emit_seg(const uint8_t* aa, const uint8_t* wob, int len, int mode, uint64_t c, std::vector<uint64_t>& out) {
    if (mode == PK_STROBE) {                                   // order-2 protein randstrobe
        constexpr int s = 4, wmin = 5, wmax = 20;
        for (int i = 0; i + s <= len; ++i) {
            const uint64_t h1 = aa4_hash(aa + i);
            const int lo = i + s + wmin, hi = std::min(i + s + wmax, len - s);
            if (lo > hi) continue;
            uint64_t best = UINT64_MAX; int bj = lo;
            for (int j = lo; j <= hi; ++j) {                  // pick 2nd strobe minimizing the link
                const uint64_t link = h1 ^ aa4_hash(aa + j);
                if (link < best) { best = link; bj = j; }
            }
            const uint64_t h = splitmix(h1 ^ (aa4_hash(aa + bj) * 0x9E3779B97F4A7C15ULL));
            if (h % c == 0) out.push_back(h);
        }
        return;
    }
    for (int i = 0; i + AAMER_K <= len; ++i) {              // PK_AAMER / PK_METAMER: sliding k=8
        if (aamer_is_low_complexity(aa + i, AAMER_K)) continue;
        uint64_t h;
        if (mode == PK_AAMER) {
            h = aamer_hash(aa + i);
        } else {                                              // PK_METAMER: 40 AA bits + 16 wobble bits
            uint64_t v = 0;
            for (int j = 0; j < AAMER_K; ++j) v |= (uint64_t)aa[i + j]        << (5 * j);
            for (int j = 0; j < AAMER_K; ++j) v |= (uint64_t)(wob[i + j] & 3) << (40 + 2 * j);
            h = splitmix(v);
        }
        if (h % c == 0) out.push_back(h);
    }
}

// 6-frame walk that yields AA + 3rd-codon-position base per inter-stop segment.
void walk_codons(std::string_view seq, int min_aa, int mode, uint64_t c, std::vector<uint64_t>& out) {
    const int n = (int)seq.size();
    const char* s = seq.data();
    std::vector<uint8_t> aa, wob;
    aa.reserve(n / 3 + 4); wob.reserve(n / 3 + 4);
    auto flush = [&]() {
        if ((int)aa.size() >= min_aa) emit_seg(aa.data(), wob.data(), (int)aa.size(), mode, c, out);
        aa.clear(); wob.clear();
    };
    for (int f = 0; f < 3; ++f) {                             // forward frames
        aa.clear(); wob.clear();
        for (int i = f; i + 2 < n; i += 3) {
            const uint8_t b0 = DNA_ENC[(uint8_t)s[i]], b1 = DNA_ENC[(uint8_t)s[i+1]], b2 = DNA_ENC[(uint8_t)s[i+2]];
            if (b0 == 0xFF || b1 == 0xFF || b2 == 0xFF) { flush(); continue; }
            const uint8_t a = CODON11[b0 * 16 + b1 * 4 + b2];
            if (a == AA_STOP) { flush(); continue; }
            aa.push_back(a); wob.push_back(b2);
        }
        flush();
    }
    for (int f = 0; f < 3; ++f) {                             // reverse-complement frames
        aa.clear(); wob.clear();
        for (int pos = n - 1 - f; pos - 2 >= 0; pos -= 3) {
            const uint8_t e0 = DNA_ENC[(uint8_t)s[pos]], e1 = DNA_ENC[(uint8_t)s[pos-1]], e2 = DNA_ENC[(uint8_t)s[pos-2]];
            if (e0 == 0xFF || e1 == 0xFF || e2 == 0xFF) { flush(); continue; }
            const uint8_t b0 = BASE_COMP[e0], b1 = BASE_COMP[e1], b2 = BASE_COMP[e2];
            const uint8_t a = CODON11[b0 * 16 + b1 * 4 + b2];
            if (a == AA_STOP) { flush(); continue; }
            aa.push_back(a); wob.push_back(b2);
        }
        flush();
    }
}

} // namespace

std::vector<uint64_t> cladesplit_aamers(std::string_view fasta, const CladeSplitConfig& cfg) {
    const uint64_t cc = cfg.frac_c > 0 ? (uint64_t)cfg.frac_c : 1;
    std::vector<uint64_t> out;
    auto contigs = genopack::check::parse_fasta(fasta);
    for (const auto& ctg : contigs)
        walk_codons(ctg.seq, cfg.min_aa, cfg.mode, cc, out);   // subsample is inline now
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

namespace {
// Per-genome aamer MULTIPLICITY (needed for the SCC single-copy criterion, which
// cladesplit_aamers discards by uniquing). Same primitive/config as the panel.
std::unordered_map<uint64_t, uint32_t> aamer_counts(std::string_view fasta, const CladeSplitConfig& cfg) {
    const uint64_t cc = cfg.frac_c > 0 ? (uint64_t)cfg.frac_c : 1;
    std::vector<uint64_t> raw;
    auto contigs = genopack::check::parse_fasta(fasta);
    for (const auto& ctg : contigs)
        walk_codons(ctg.seq, cfg.min_aa, cfg.mode, cc, raw);
    std::unordered_map<uint64_t, uint32_t> cnt;
    cnt.reserve(raw.size());
    for (uint64_t h : raw) ++cnt[h];
    return cnt;
}
// core_dup of one genome's count map over a genus SCC hash set.
float core_dup_of(const std::unordered_map<uint64_t, uint32_t>& cnt,
                  const std::unordered_set<uint64_t>& scc) {
    uint32_t present = 0, dup = 0;
    for (const auto& [h, c] : cnt) if (scc.count(h)) { ++present; if (c >= 2) ++dup; }
    return present ? (float)dup / (float)present : 0.0f;
}
} // namespace

void cladesplit_build(const std::filesystem::path& tsv,
                      const std::filesystem::path& out,
                      const CladeSplitConfig& cfg, int threads) {
    auto recs = parse_tsv_records(tsv);
    // Pre-assign genus clade ids (sequential, cheap).
    std::unordered_map<std::string, int32_t> clade_id;
    std::vector<std::string>                 clade_names;
    std::vector<int32_t>                      rec_cid(recs.size(), -1);
    for (size_t i = 0; i < recs.size(); ++i) {
        std::string genus = extract_rank(taxonomy_of(recs[i]), "g__");
        if (genus.empty()) continue;
        auto it = clade_id.find(genus);
        if (it == clade_id.end()) {
            rec_cid[i] = (int32_t)clade_names.size();
            clade_id.emplace(genus, rec_cid[i]);
            clade_names.push_back(genus);
        } else rec_cid[i] = it->second;
    }
    const unsigned nt = threads > 0 ? (unsigned)threads : std::max(1u, std::thread::hardware_concurrency());
    spdlog::info("cladesplit build: {} genomes, {} genera, {} threads", recs.size(), clade_names.size(), nt);

    // Sharded global aamer diagnosticity map, 256-way by low 8 hash bits. Each shard
    // is a compact open-addressing table {uint64 key, uint32 val} at <=70% load —
    // ~17 B/entry vs ~50 B for std::unordered_map, the dominant Phase-A footprint at
    // GTDB scale. val packing: 0 = empty; bit31 = multi-genus (non-diagnostic);
    // bits0..30 = cid+1 (single-genus diagnostic while bit31 clear).
    constexpr int NSH = 256;
    struct DiagShard {
        std::vector<uint64_t> keys; std::vector<uint32_t> vals;
        uint32_t size = 0, mask = 0, occ = 0;
        std::mutex mu;
        void reset(uint32_t cap) { size = cap; mask = cap - 1; keys.assign(cap, 0); vals.assign(cap, 0); occ = 0; }
        inline void reput_(uint64_t h, uint32_t v) {          // reinsert on grow, no update logic
            uint32_t s = (uint32_t)(h >> 8) & mask;
            while (vals[s]) s = (s + 1) & mask;
            keys[s] = h; vals[s] = v; ++occ;
        }
        void grow() {
            std::vector<uint64_t> ok = std::move(keys); std::vector<uint32_t> ov = std::move(vals);
            const uint32_t osz = size; reset(size << 1);
            for (uint32_t i = 0; i < osz; ++i) if (ov[i]) reput_(ok[i], ov[i]);
        }
        inline void upsert(uint64_t h, int32_t cid) {
            if (size == 0) reset(256);
            if ((occ + 1) * 10 >= size * 7) grow();
            uint32_t s = (uint32_t)(h >> 8) & mask;
            while (vals[s]) {
                if (keys[s] == h) {                            // seen before → mark multi if new genus
                    const uint32_t v = vals[s];
                    if (!(v & 0x80000000u) && (int32_t)((v & 0x7fffffffu) - 1) != cid) vals[s] = v | 0x80000000u;
                    return;
                }
                s = (s + 1) & mask;
            }
            keys[s] = h; vals[s] = (uint32_t)(cid + 1); ++occ;  // first sighting
        }
    };
    std::vector<DiagShard> shards(NSH);
    std::atomic<size_t> next{0}, done{0};
    auto worker = [&]() {
        std::vector<std::vector<uint64_t>> buckets(NSH);
        size_t i;
        while ((i = next.fetch_add(1)) < recs.size()) {
            if (rec_cid[i] < 0) continue;
            std::string fasta;
            try { fasta = load_fasta(recs[i].file_path); } catch (const std::exception&) { continue; }
            auto ms = cladesplit_aamers(fasta, cfg);
            const int32_t cid = rec_cid[i];
            for (auto& b : buckets) b.clear();
            for (uint64_t h : ms) buckets[h & (NSH - 1)].push_back(h);
            for (int sidx = 0; sidx < NSH; ++sidx) {
                if (buckets[sidx].empty()) continue;
                std::lock_guard<std::mutex> lk(shards[sidx].mu);
                for (uint64_t h : buckets[sidx]) shards[sidx].upsert(h, cid);
            }
            if (size_t d = ++done; d % 1000 == 0) spdlog::info("cladesplit build: {}/{}", d, recs.size());
        }
    };
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < nt; ++t) pool.emplace_back(worker);
    for (auto& th : pool) th.join();

    std::vector<std::pair<uint64_t, int32_t>> diag;
    size_t total_mm = 0;
    for (auto& sh : shards) {
        total_mm += sh.occ;
        for (uint32_t s = 0; s < sh.size; ++s) {
            const uint32_t v = sh.vals[s];
            if (v && !(v & 0x80000000u)) diag.push_back({ sh.keys[s], (int32_t)((v & 0x7fffffffu) - 1) });
        }
    }
    std::sort(diag.begin(), diag.end());
    spdlog::info("cladesplit build: {} clades, {} aamers, {} genus-diagnostic ({:.1f}%)",
                 clade_names.size(), total_mm, diag.size(),
                 total_mm == 0 ? 0.0 : 100.0 * diag.size() / total_mm);

    // ── v3: per-genus single-copy-core (SCC) index + clean core_dup ceiling ──
    // SCC is NOT the diagnostic set (a prevalence-core aamer may be shared across genera).
    // It is derived per genus from its own genomes' aamer MULTIPLICITY: prevalent AND
    // single-copy. Only genera with >= kSccMinRefs get an SCC set; others -> core_dup NA.
    // Streamed one genus at a time, parallel across genera: peak RAM is
    // O(refs-per-genus x threads), not O(all genomes). This replaces the previous
    // gcounts(recs.size()) materialisation that held every participating genome's
    // count map resident and did not scale to full GTDB.
    std::vector<std::vector<size_t>> genus_recs(clade_names.size());
    for (size_t i = 0; i < recs.size(); ++i) if (rec_cid[i] >= 0) genus_recs[rec_cid[i]].push_back(i);

    struct SccOut { int32_t clade; float ceiling; std::vector<uint64_t> hashes; };
    std::vector<SccOut> scc_out;
    std::mutex scc_mu;
    std::atomic<size_t> gnext{0}, gdone{0};
    auto scc_worker = [&]() {
        size_t cid;
        while ((cid = gnext.fetch_add(1)) < genus_recs.size()) {
            const auto& idxs = genus_recs[cid];
            if (idxs.size() < kSccMinRefs) continue;
            // Load only this genus's genomes — resident for this genus alone, then freed.
            std::vector<std::unordered_map<uint64_t, uint32_t>> cg;
            cg.reserve(idxs.size());
            for (size_t i : idxs) {
                std::string fasta;
                try { fasta = load_fasta(recs[i].file_path); } catch (const std::exception&) { continue; }
                cg.push_back(aamer_counts(fasta, cfg));
            }
            std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>> pm;  // hash -> (present, multi)
            uint32_t N = 0;
            for (auto& c : cg) {
                if (c.empty()) continue;
                ++N;
                for (auto& [h, cnt] : c) { auto& e = pm[h]; ++e.first; if (cnt >= 2) ++e.second; }
            }
            if (N < kSccMinRefs) continue;
            const uint32_t need_prev = (uint32_t)std::ceil(kSccPrevalence * N);
            std::unordered_set<uint64_t> scc;
            for (auto& [h, e] : pm)
                if (e.first >= need_prev && e.second <= (uint32_t)std::floor(kSccMultiFrac * e.first)) scc.insert(h);
            if (scc.empty()) continue;
            std::vector<float> cds;
            for (auto& c : cg) if (!c.empty()) cds.push_back(core_dup_of(c, scc));
            std::sort(cds.begin(), cds.end());
            const float med = cds[cds.size() / 2];
            std::vector<float> ad(cds.size());
            for (size_t k = 0; k < cds.size(); ++k) ad[k] = std::fabs(cds[k] - med);
            std::sort(ad.begin(), ad.end());
            const float mad = ad[ad.size() / 2];
            const float ceiling = med + 20.0f * mad;
            std::vector<uint64_t> sv(scc.begin(), scc.end());
            std::sort(sv.begin(), sv.end());
            std::lock_guard<std::mutex> lk(scc_mu);
            spdlog::info("cladesplit scc: genus '{}' refs={} SCC={} ceiling={:.5f} (median={:.5f} MAD={:.5f})",
                         clade_names[cid], N, sv.size(), ceiling, med, mad);
            scc_out.push_back({ (int32_t)cid, ceiling, std::move(sv) });
            if (size_t d = ++gdone; d % 200 == 0) spdlog::info("cladesplit scc: {} genera done", d);
        }
    };
    {
        std::vector<std::thread> gp;
        for (unsigned t = 0; t < nt; ++t) gp.emplace_back(scc_worker);
        for (auto& th : gp) th.join();
    }
    // Deterministic on-disk order (parallel completion order is nondeterministic).
    std::sort(scc_out.begin(), scc_out.end(),
              [](const SccOut& a, const SccOut& b) { return a.clade < b.clade; });

    AppendWriter w;
    w.create(out);
    CspHeader h{};
    w.append(&h, sizeof(h));   // placeholder, patched below
    for (const auto& [hash, cid] : diag) w.append(&hash, sizeof(hash));
    for (const auto& [hash, cid] : diag) w.append(&cid, sizeof(cid));
    uint64_t names_off = w.current_offset();
    for (const auto& nm : clade_names) {
        uint16_t len = (uint16_t)nm.size();
        w.append(&len, sizeof(len));
        w.append(nm.data(), nm.size());
    }
    uint64_t scc_off = w.current_offset();
    uint32_t n_genus = (uint32_t)scc_out.size();
    w.append(&n_genus, sizeof(n_genus));
    for (const auto& e : scc_out) {
        w.append(&e.clade, sizeof(e.clade));
        w.append(&e.ceiling, sizeof(e.ceiling));
        uint32_t ns = (uint32_t)e.hashes.size();
        w.append(&ns, sizeof(ns));
        if (ns) w.append(e.hashes.data(), (size_t)ns * sizeof(uint64_t));
    }
    h.magic = kCspMagicV3;
    h.n_aamers = diag.size();
    h.n_clades = clade_names.size();
    h.names_off = names_off;
    h.frac_c = (uint16_t)cfg.frac_c;
    h.min_aa    = (uint16_t)cfg.min_aa;
    h.mode      = (uint16_t)cfg.mode;
    h.scc_off   = scc_off;
    w.write_at(0, &h, sizeof(h));
    w.flush();
    w.close();
    spdlog::info("cladesplit build: wrote {}", out.string());
}

void CladeSplitPanel::open(const std::filesystem::path& path) {
    m_.open(path);
    if (m_.size() < sizeof(CspHeader)) throw std::runtime_error("cladesplit: file too small");
    CspHeader h{};
    // v3 header (48B) shares its first 40 bytes with v2, so a 48-byte read is safe for
    // both (the file has diagnostic hashes after the header); scc_off is only used for v3.
    std::memcpy(&h, m_.data(), sizeof(h));
    if (h.magic != kCspMagic && h.magic != kCspMagicV3)
        throw std::runtime_error("cladesplit: bad magic (not a .csp v2/v3 panel — rebuild)");
    const bool v3 = (h.magic == kCspMagicV3);
    const size_t hdr = v3 ? kCspHdrV3 : kCspHdrV2;   // v2 hashes start at 40, v3 at 48
    n_aamers_   = h.n_aamers;
    panel_c_      = h.frac_c ? h.frac_c : 30;
    panel_min_aa_ = h.min_aa    ? h.min_aa    : 8;
    panel_mode_   = h.mode;
    hashes_    = reinterpret_cast<const uint64_t*>(m_.data() + hdr);
    clade_ids_ = reinterpret_cast<const uint32_t*>(m_.data() + hdr + n_aamers_ * sizeof(uint64_t));
    const uint8_t* p   = m_.data() + h.names_off;
    const uint8_t* end = m_.data() + m_.size();
    clade_names_.reserve(h.n_clades);
    for (uint64_t i = 0; i < h.n_clades && p + 2 <= end; ++i) {
        uint16_t len; std::memcpy(&len, p, 2); p += 2;
        if (p + len > end) break;
        clade_names_.emplace_back(reinterpret_cast<const char*>(p), len);
        p += len;
    }
    // v3: per-genus SCC index (clade_id, ceiling, sorted SCC hashes).
    if (v3 && h.scc_off > 0 && h.scc_off + 4 <= m_.size()) {
        const uint8_t* q    = m_.data() + h.scc_off;
        const uint8_t* qend = m_.data() + m_.size();
        uint32_t ng; std::memcpy(&ng, q, 4); q += 4;
        for (uint32_t g = 0; g < ng && q + 12 <= qend; ++g) {
            int32_t  clade;   std::memcpy(&clade, q, 4);   q += 4;
            float    ceiling; std::memcpy(&ceiling, q, 4); q += 4;
            uint32_t ns;      std::memcpy(&ns, q, 4);      q += 4;
            if (q + (size_t)ns * 8 > qend) break;
            SccEntry se;
            se.ceiling = ceiling;
            se.hashes.reserve(ns);
            for (uint32_t k = 0; k < ns; ++k) { uint64_t hh; std::memcpy(&hh, q, 8); q += 8; se.hashes.insert(hh); }
            scc_[clade] = std::move(se);
        }
    }
    // Blocked Bloom prefilter over the diagnostic aamers (10 bits/elem, k=4).
    if (n_aamers_ > 0) {
        bloom_blocks_ = std::max<uint32_t>(1, (uint32_t)((double)n_aamers_ * 10.0 / 512.0));
        bloom_.assign((size_t)bloom_blocks_ * 8, 0);
        for (uint64_t i = 0; i < n_aamers_; ++i) {
            const uint32_t b = (uint32_t)(hashes_[i] % bloom_blocks_) * 8;
            uint64_t mix = hashes_[i];
            for (int k = 0; k < 4; ++k) {
                mix = mix * 0x9e3779b97f4a7c15ULL ^ (mix >> 32);
                bloom_[b + ((mix >> 9) & 7)] |= (uint64_t)1 << (mix & 63);
            }
        }
    }
}

std::string_view CladeSplitPanel::clade_name(int32_t id) const {
    if (id < 0 || (size_t)id >= clade_names_.size()) return {};
    return clade_names_[id];
}

namespace {
std::string fmt_line(const std::string& acc, const CladeSplitScore& s, const CladeSplitPanel& panel) {
    char num[96];
    std::string l = acc; l += '\t';
    std::snprintf(num, sizeof(num), "%.6g\t%.6g\t", s.minority_fraction, s.redundancy_fraction); l += num;
    l += (s.majority_clade >= 0 ? std::string(panel.clade_name(s.majority_clade)) : "-"); l += '\t';
    l += (s.minority_clade >= 0 ? std::string(panel.clade_name(s.minority_clade)) : "-"); l += '\t';
    l += std::to_string(s.n_votes); l += '\t';
    // core_dup / core_dup_excess: "NA" when the majority genus has no SCC set (v2 panel or underpopulated genus).
    if (std::isnan(s.core_dup)) l += "NA\tNA";
    else { std::snprintf(num, sizeof(num), "%.6g\t%.6g", s.core_dup, s.core_dup_excess); l += num; }
    l += '\n';
    return l;
}
unsigned resolve_threads(int t) { return t > 0 ? (unsigned)t : std::max(1u, std::thread::hardware_concurrency()); }
} // namespace

void cladesplit_score_tsv(const std::filesystem::path& tsv,
                          const std::filesystem::path& panel_path,
                          const std::filesystem::path& out,
                          const CladeSplitConfig& cfg, int threads) {
    auto recs = parse_tsv_records(tsv);
    CladeSplitPanel panel;
    panel.open(panel_path);
    const unsigned nt = resolve_threads(threads);
    spdlog::info("cladesplit score: {} genomes vs panel ({} aamers, {} clades), {} threads",
                 recs.size(), panel.n_aamers(), panel.n_clades(), nt);
    std::vector<std::string> lines(recs.size());
    std::atomic<size_t> next{0}, done{0};
    auto worker = [&]() {
        size_t i;
        while ((i = next.fetch_add(1)) < recs.size()) {
            std::string fasta;
            try { fasta = load_fasta(recs[i].file_path); } catch (const std::exception&) { continue; }
            lines[i] = fmt_line(recs[i].accession, panel.score(fasta, cfg), panel);
            if (size_t d = ++done; d % 500 == 0) spdlog::info("cladesplit score: {}/{}", d, recs.size());
        }
    };
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < nt; ++t) pool.emplace_back(worker);
    for (auto& th : pool) th.join();
    std::ofstream o(out);
    o << "accession\tminority_fraction\tredundancy_fraction\tmajority_clade\tminority_clade\tn_votes\tcore_dup\tcore_dup_excess\n";
    for (const auto& l : lines) if (!l.empty()) o << l;
    spdlog::info("cladesplit score: wrote {}", out.string());
}

void cladesplit_dump_aamers(const std::filesystem::path& tsv,
                            const std::filesystem::path& out,
                            const CladeSplitConfig& cfg, int threads) {
    auto recs = parse_tsv_records(tsv);
    const unsigned nt = resolve_threads(threads);
    spdlog::info("cladesplit aamers: extracting from {} genomes, {} threads", recs.size(), nt);
    std::vector<std::vector<uint64_t>> sets(recs.size());
    std::atomic<size_t> next{0}, done{0};
    auto worker = [&]() {
        size_t i;
        while ((i = next.fetch_add(1)) < recs.size()) {
            std::string fasta;
            try { fasta = load_fasta(recs[i].file_path); } catch (const std::exception&) { continue; }
            sets[i] = cladesplit_aamers(fasta, cfg);
            if (size_t d = ++done; d % 500 == 0) spdlog::info("cladesplit aamers: {}/{}", d, recs.size());
        }
    };
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < nt; ++t) pool.emplace_back(worker);
    for (auto& th : pool) th.join();
    std::ofstream o(out, std::ios::binary);
    for (size_t i = 0; i < recs.size(); ++i) {
        if (sets[i].empty()) continue;
        const std::string& acc = recs[i].accession;
        uint16_t al = (uint16_t)acc.size();
        uint32_t nh = (uint32_t)sets[i].size();
        o.write((const char*)&al, 2); o.write(acc.data(), al);
        o.write((const char*)&nh, 4);
        o.write((const char*)sets[i].data(), (std::streamsize)nh * 8);
    }
    spdlog::info("cladesplit aamers: wrote {} ({} genomes)", out.string(), recs.size());
}

void cladesplit_score_gpk(const std::filesystem::path& gpk,
                          const std::filesystem::path& panel_path,
                          const std::filesystem::path& out,
                          const CladeSplitConfig& cfg, int threads) {
    CladeSplitPanel panel;
    panel.open(panel_path);
    ArchiveReader ar;
    ar.open(gpk);
    std::vector<std::string> accs;
    ar.scan_genome_accessions([&](std::string_view a, GenomeId g) {
        if (!ar.is_deleted(g)) accs.push_back(std::string(a));
    });
    std::sort(accs.begin(), accs.end());
    const unsigned nt = resolve_threads(threads);
    spdlog::info("cladesplit score-gpk: {} live genomes vs panel ({} aamers, {} clades), {} threads",
                 accs.size(), panel.n_aamers(), panel.n_clades(), nt);
    std::vector<std::string> lines(accs.size());
    std::atomic<size_t> done{0};
    ar.visit_shard_batches_parallel(accs, static_cast<int>(nt),
        [&](ArchiveReader::ShardBatch& batch) {
            for (auto& [idx, eg] : batch) {
                lines[idx] = fmt_line(accs[idx], panel.score(eg.fasta, cfg), panel);
                if (size_t d = ++done; d % 500 == 0)
                    spdlog::info("cladesplit score-gpk: {}/{}", d, accs.size());
            }
        });
    std::ofstream o(out);
    o << "accession\tminority_fraction\tredundancy_fraction\tmajority_clade\tminority_clade\tn_votes\tcore_dup\tcore_dup_excess\n";
    for (const auto& l : lines) if (!l.empty()) o << l;
    spdlog::info("cladesplit score-gpk: wrote {}", out.string());
}

CladeSplitScore CladeSplitPanel::score(std::string_view fasta, const CladeSplitConfig&) const {
    CladeSplitScore s;
    // Reproduce the panel's own build config exactly (the cfg arg is ignored): subsample
    // density, segment length and aamer mode must match how the panel was constructed.
    const uint64_t cc = panel_c_ > 0 ? (uint64_t)panel_c_ : 1;
    // Inline-subsampled aamer multiset (no giant raw stream).
    std::vector<uint64_t> ms;
    {
        auto contigs = genopack::check::parse_fasta(fasta);
        for (const auto& ctg : contigs) walk_codons(ctg.seq, panel_min_aa_, panel_mode_, cc, ms);
    }
    // Blocked-Bloom prefilter rejects ~99% of non-panel aamers in one cache line; only
    // the candidates land in a small local count map. No per-genome sort of the full set.
    std::unordered_map<uint64_t, uint32_t> hit_cnt;
    hit_cnt.reserve(1024);
    for (uint64_t h : ms)
        if (bloom_might(h)) ++hit_cnt[h];

    std::vector<uint32_t> votes(clade_names_.size(), 0);
    uint32_t total = 0, diag_dup = 0;
    for (const auto& [h, cnt] : hit_cnt) {
        const uint64_t* it = std::lower_bound(hashes_, hashes_ + n_aamers_, h);
        if (it != hashes_ + n_aamers_ && *it == h) {        // confirm (Bloom can false-positive)
            votes[clade_ids_[it - hashes_]]++;
            ++total;
            if (cnt >= 2) ++diag_dup;
        }
    }
    s.n_votes = total;
    if (total > 0) s.redundancy_fraction = (float)diag_dup / (float)total;
    if (total == 0) return s;
    uint32_t maj_v = 0, min_v = 0; int32_t maj_c = -1, min_c = -1;
    for (size_t c = 0; c < votes.size(); ++c) {
        if (votes[c] > maj_v) { min_v = maj_v; min_c = maj_c; maj_v = votes[c]; maj_c = (int32_t)c; }
        else if (votes[c] > min_v) { min_v = votes[c]; min_c = (int32_t)c; }
    }
    s.majority_clade = maj_c;
    s.minority_clade = min_c;
    s.minority_fraction = 1.0f - (float)maj_v / (float)total;

    // v3: SCC-restricted core duplication for the placed (majority) genus. Structurally
    // excludes accessory/multi-copy families, so it flags same-genus/strain redundancy
    // that redundancy_fraction (diagnostic-only) and CheckM/CheckM2 miss. NAN when the
    // majority genus has no SCC set (underpopulated/novel) -> defers to observability cap.
    if (auto sit = scc_.find(maj_c); sit != scc_.end() && !sit->second.hashes.empty()) {
        const auto& sset = sit->second.hashes;
        std::unordered_map<uint64_t, uint32_t> sc_cnt;
        for (uint64_t h : ms) if (sset.count(h)) ++sc_cnt[h];
        uint32_t present = (uint32_t)sc_cnt.size(), dup = 0;
        for (const auto& [h, c] : sc_cnt) if (c >= 2) ++dup;
        s.core_dup = present ? (float)dup / (float)present : 0.0f;
        const float ceil_ = sit->second.ceiling;
        float ex = ceil_ > 0.0f ? (s.core_dup - ceil_) / ceil_ : 0.0f;
        s.core_dup_excess = ex < 0.0f ? 0.0f : (ex > 1.0f ? 1.0f : ex);
    }
    return s;
}

} // namespace genopack
