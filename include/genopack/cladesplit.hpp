#pragma once
#include "mmap_file.hpp"
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace genopack {

// ── Clade-split contamination tier (.csp) ─────────────────────────────────────
// A fast protein-aamer contamination signal (GUNC-style, k-mer form). Each genome
// is 6-frame translated; inter-stop AA≥8 segments yield k=8 amino-acid aamers
// (reusing aamer.hpp). A panel built from CLEAN reference genomes maps each
// "genus-diagnostic" aamer (one carried by exactly one genus) → that genus. At
// check time a query's aamers vote their diagnostic genus; a clean genome resolves
// to one genus, a contaminated one splits — minority_fraction = votes to non-majority
// genus / total diagnostic votes. No per-genus model, no contig length floor that
// kills TNF/FMH — works on fragmented MAGs and across all taxonomic distances.

constexpr uint64_t kCspMagic = 0x4753435350000002ULL; // "GCSP" v2: header carries build config

// Protein k-mer primitive (build and score must match):
//   AAMER   = 8 amino acids only (the production signal; not a Metabuli metamer)
//   METAMER = 8 AA + 8 third-codon-position bases (joint AA+DNA, true Metabuli metamer)
//   STROBE  = order-2 protein randstrobe (gapped, longer span, mutation-robust)
enum ProteinKmer : int { PK_AAMER = 0, PK_METAMER = 1, PK_STROBE = 2 };

struct CladeSplitConfig {
    int frac_c = 30;          // FracMinHash density 1/c on the primitive
    int min_aa = 8;           // inter-stop AA segment minimum length
    int mode   = PK_AAMER;    // protein k-mer primitive (build and score must match)
};

struct CladeSplitScore {
    float    minority_fraction   = 0.0f; // votes to non-majority clade / total (inter-lineage)
    float    redundancy_fraction = 0.0f; // diagnostic aamers seen ≥2× / distinct diagnostic (intra-lineage dup)
    uint32_t n_votes             = 0;
    int32_t  majority_clade      = -1;
    int32_t  minority_clade      = -1;
};

// Build a panel from a TSV (accession, file_path, taxonomy[, ...]); writes `out`.
// Parallel over genomes (sharded global aamer→genus map). threads 0 = all cores.
void cladesplit_build(const std::filesystem::path& tsv,
                      const std::filesystem::path& out,
                      const CladeSplitConfig& cfg = {},
                      int threads = 0);

class CladeSplitPanel {
public:
    void open(const std::filesystem::path& path);
    // Scoring always uses the panel's own build config (stored in the .csp header), so
    // the result is identical to how the panel was validated — the cfg arg is ignored
    // for aamer extraction and kept only for source compatibility.
    CladeSplitScore score(std::string_view fasta, const CladeSplitConfig& cfg = {}) const;
    std::string_view clade_name(int32_t id) const;
    uint64_t n_aamers() const { return n_aamers_; }
    uint64_t n_clades()   const { return clade_names_.size(); }
    CladeSplitConfig config() const { return { panel_c_, panel_min_aa_, panel_mode_ }; }

    // Cache-line blocked Bloom prefilter over the diagnostic aamers: rejects ~99%
    // of non-panel aamers in ONE cache line, so the binary search almost never runs.
    bool bloom_might(uint64_t h) const noexcept {
        if (bloom_.empty()) return true;
        const uint32_t b = (uint32_t)(h % bloom_blocks_) * 8;
        uint64_t mix = h;
        for (int i = 0; i < 4; ++i) {
            mix = mix * 0x9e3779b97f4a7c15ULL ^ (mix >> 32);
            if (!((bloom_[b + ((mix >> 9) & 7)] >> (mix & 63)) & 1)) return false;
        }
        return true;
    }

private:
    MmapFileReader           m_;
    const uint64_t*          hashes_    = nullptr;  // sorted, n_aamers_
    const uint32_t*          clade_ids_ = nullptr;  // parallel to hashes_
    uint64_t                 n_aamers_ = 0;
    std::vector<std::string> clade_names_;
    std::vector<uint64_t>    bloom_;                 // blocked Bloom, 8 u64 per 64-byte block
    uint32_t                 bloom_blocks_ = 0;
    uint16_t                 panel_c_      = 30;     // build config — score uses these, not the cfg arg
    uint16_t                 panel_min_aa_ = 8;
    uint16_t                 panel_mode_   = PK_AAMER;
};

// Extract per-genome subsampled, sorted-unique k=8 protein aamers (6-frame).
std::vector<uint64_t> cladesplit_aamers(std::string_view fasta, const CladeSplitConfig& cfg);

// Score every genome in a TSV (accession, file_path, ...) against a panel; write a
// per-genome TSV (accession, minority_fraction, redundancy_fraction, ...).
void cladesplit_score_tsv(const std::filesystem::path& tsv,
                          const std::filesystem::path& panel,
                          const std::filesystem::path& out,
                          const CladeSplitConfig& cfg = {},
                          int threads = 0);   // 0 = hardware_concurrency

// Score every live genome IN a .gpk archive directly (sequence streamed from shards),
// writing the same per-genome TSV. This is the "flag contamination in a genopack file" path.
void cladesplit_score_gpk(const std::filesystem::path& gpk,
                          const std::filesystem::path& panel,
                          const std::filesystem::path& out,
                          const CladeSplitConfig& cfg = {},
                          int threads = 0);

// Dump per-genome sorted-unique subsampled aamer sets (for genus-core / completeness
// R&D). Binary stream of records: u16 acc_len, acc bytes, u32 n_hashes, u64 hashes[].
void cladesplit_dump_aamers(const std::filesystem::path& tsv,
                            const std::filesystem::path& out,
                            const CladeSplitConfig& cfg = {},
                            int threads = 0);

} // namespace genopack
