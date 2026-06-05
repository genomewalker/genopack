#pragma once
// ── SEC_QCONTIG — per-contig quality overlay ──────────────────────────────────
// The per-genome quality columns (QCOL) collapse a genome to one contamination
// number, but contamination is a per-CONTIG phenomenon — a genome scores dirty
// because SOME contigs are foreign. QCONTIG persists the per-contig signal that
// `check` already computes (pass_b's ContigFlag: offset, length, TNF/GCOV-outlier
// score, containment-leakage score), keyed by (genome_id, contig_index), so a
// consumer can see WHICH contigs drive a score and decontaminate by dropping them.
//
// Records are sorted by (genome_id, contig_index); a reader binary-searches the
// first record of a genome and scans forward. Layout: [QcontigHeader][Record[]].
#include "format.hpp"
#include "mmap_file.hpp"
#include <algorithm>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace genopack {

static constexpr uint16_t QCONTIG_FORMAT_VERSION = 1;

struct QcontigHeader {           // 32 bytes
    uint32_t magic;              //  0  SEC_QCONTIG
    uint16_t format_version;     //  4
    uint16_t flags;              //  6
    uint64_t n_rows;             //  8  number of contig records
    uint64_t records_off;        // 16  == sizeof(QcontigHeader)
    uint64_t _reserved;          // 24
};
static_assert(sizeof(QcontigHeader) == 32);

struct QcontigRecord {           // 64 bytes
    uint64_t genome_id;          //  0
    uint32_t contig_index;       //  8  ordinal within the genome (0-based)
    uint32_t contig_offset;      // 12  byte offset within the genome's sequence
    uint32_t contig_length;      // 16  bp
    float    tnf_score;          // 20  coarse centroid-distance ratio (NaN = not scored)
    float    leakage_score;      // 24  per-contig containment-leakage score (NaN = not scored)
    float    gcov_t2_pct;        // 28  GCOV Hotelling-T² percentile (within-genus Mahalanobis)
    float    gcov_spe_pct;       // 32  GCOV SPE percentile (residual; catches cross-genus/foreign LONG contigs)
    // Protein-aamer foreign-containment channel — the SMALL-contig signal (TNF-independent).
    float    prot_foreign_frac;  // 36  foreign_specific / (host_specific + foreign_specific)
    float    prot_loglr;         // 40  log((foreign_specific+1)/(host_specific+1))  [v1 presence-only]
    uint16_t prot_classifiable;  // 44  host_specific + foreign_specific aamers (abstain if < ~32)
    uint16_t prot_host_specific; // 46
    uint16_t prot_foreign_specific;//48
    uint16_t prot_family_hits;   // 50
    uint32_t prot_best_genus;    // 52  22-bit tag of the dominant foreign genus
    uint8_t  prot_flags;         // 56  bit0 ABSTAIN_LOW_N, bit1 MOBILE_NATIVE, bit2 NOVEL_GENUS_FALLBACK
    uint8_t  _pad0;              // 57
    uint16_t _pad1;              // 58
    uint32_t _pad2;              // 60
};
static_assert(sizeof(QcontigRecord) == 64);

// ── Writer ────────────────────────────────────────────────────────────────────
class QcontigWriter {
public:
    // Append a genome's per-contig flags (ordinal index assigned in input order).
    template <class ContigFlagT>
    void add_genome(uint64_t genome_id, const std::vector<ContigFlagT>& flags) {
        for (uint32_t i = 0; i < flags.size(); ++i) {
            const auto& f = flags[i];
            recs_.push_back(QcontigRecord{genome_id, i, f.contig_offset, f.contig_length,
                                          f.tnf_score, f.leakage_score, f.gcov_t2_pct, f.gcov_spe_pct,
                                          f.prot_foreign_frac, f.prot_loglr, f.prot_classifiable,
                                          f.prot_host_specific, f.prot_foreign_specific, f.prot_family_hits,
                                          f.prot_best_genus, f.prot_flags, 0, 0, 0});
        }
    }
    size_t n_rows() const { return recs_.size(); }

    SectionDesc finalize(AppendWriter& w, uint64_t section_id);

private:
    std::vector<QcontigRecord> recs_;
};

// ── Reader ────────────────────────────────────────────────────────────────────
class QcontigReader {
public:
    void open(const uint8_t* data, uint64_t offset, uint64_t size) {
        if (size < sizeof(QcontigHeader)) throw std::runtime_error("QCONTIG section too small");
        data_   = data + offset;
        header_ = reinterpret_cast<const QcontigHeader*>(data_);
        if (header_->magic != SEC_QCONTIG) throw std::runtime_error("QCONTIG: bad magic");
        const uint64_t end = header_->records_off + header_->n_rows * sizeof(QcontigRecord);
        if (end > size) throw std::runtime_error("QCONTIG section truncated");
        recs_ = reinterpret_cast<const QcontigRecord*>(data_ + header_->records_off);
    }

    bool     is_open() const noexcept { return data_ != nullptr; }
    uint64_t n_rows() const noexcept { return header_ ? header_->n_rows : 0; }

    // Contiguous range of records for one genome (empty if none). Records are
    // sorted by (genome_id, contig_index), so a genome's contigs are contiguous.
    std::span<const QcontigRecord> contigs_for(uint64_t genome_id) const noexcept {
        if (!header_ || header_->n_rows == 0) return {};
        const QcontigRecord* begin = recs_;
        const QcontigRecord* end   = recs_ + header_->n_rows;
        auto lo = std::lower_bound(begin, end, genome_id,
            [](const QcontigRecord& r, uint64_t g) { return r.genome_id < g; });
        auto hi = std::upper_bound(lo, end, genome_id,
            [](uint64_t g, const QcontigRecord& r) { return g < r.genome_id; });
        return {lo, static_cast<size_t>(hi - lo)};
    }

    const QcontigRecord* records() const noexcept { return recs_; }

private:
    const uint8_t*       data_   = nullptr;
    const QcontigHeader* header_ = nullptr;
    const QcontigRecord* recs_   = nullptr;
};

} // namespace genopack
