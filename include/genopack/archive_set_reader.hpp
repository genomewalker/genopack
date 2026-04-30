#pragma once
#include "archive.hpp"
#include "types.hpp"
#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace genopack {

// ── ArchiveSetReader ─────────────────────────────────────────────────────────
// Read facade over one or more genopack archives. Transparently handles:
//   single archive:  <dir> containing toc.bin (or path with .gpk extension)
//   multipart set:   <dir> containing part_*.gpk subdirs (or .gpk files)
//
// Discovery for multipart layouts is done from part_*.gpk entries only;
// part_N.meta.tsv sidecars are ignored (they may be stale relative to the
// archives' canonical metadata).
//
// Methods route by accession when applicable. For multi-part queries that
// produce GenomeMeta, results are wrapped in LocatedGenomeMeta to preserve
// (part_index, part_path) since local GenomeId values can collide across parts.

struct LocatedGenomeMeta {
    size_t                 part_index;
    std::filesystem::path  part_path;
    GenomeMeta             meta;  // local genome_id within its part
};

class ArchiveSetReader {
public:
    ArchiveSetReader();
    ~ArchiveSetReader();
    ArchiveSetReader(const ArchiveSetReader&) = delete;
    ArchiveSetReader& operator=(const ArchiveSetReader&) = delete;
    ArchiveSetReader(ArchiveSetReader&&) noexcept;
    ArchiveSetReader& operator=(ArchiveSetReader&&) noexcept;

    void open(const std::filesystem::path& path);
    void close();
    bool is_open() const;

    bool   is_multipart() const;
    size_t part_count() const;
    const std::vector<std::filesystem::path>& part_paths() const;

    // Aggregate stats: sums n_shards, n_genomes_total, n_genomes_live, raw_bp,
    // compressed_bytes; takes max(generation); recomputes compression_ratio.
    ArchiveReader::ArchiveStats archive_stats() const;

    // Sum across parts.
    size_t count(const ExtractQuery& q) const;

    // Per-part filter, wrapped with part info.
    std::vector<LocatedGenomeMeta> filter_meta(const ExtractQuery& q) const;

    // If q.accessions present, routes per-accession to the owning part.
    // Otherwise delegates to each part and concatenates (q.limit enforced
    // globally if > 0).
    std::vector<ExtractedGenome> extract(const ExtractQuery& q) const;

    // Probes parts in order. First hit wins.
    std::optional<ExtractedGenome> fetch_by_accession(std::string_view accession) const;

    // Groups accessions by owning part, calls each part's batch fetch once,
    // reassembles in original order. Misses → nullopt at that index.
    std::vector<std::optional<ExtractedGenome>>
    batch_fetch_by_accessions(const std::vector<std::string>& accessions) const;

    std::optional<std::string>
    fetch_sequence_slice_by_accession(std::string_view accession,
                                      uint64_t start,
                                      uint64_t length) const;

    std::optional<std::string>
    taxonomy_for_accession(std::string_view accession) const;

    // Iterates each part's TAXN section in part order.
    void scan_taxonomy(const std::function<void(std::string_view accession,
                                                std::string_view taxonomy)>& cb) const;

    // Locates the part owning `accession` and returns its TaxonomyTree
    // alongside the part index. Caller should query `taxid_for_accession`
    // on the returned tree. Returns nullopt if accession unknown or no
    // taxonomy data.
    struct LocatedTaxonomy {
        size_t                  part_index;
        const TaxonomyTree*     tree;     // owned by part; valid while reader is open
    };
    std::optional<LocatedTaxonomy>
    taxonomy_tree_for_accession(std::string_view accession) const;

    // Aggregate node/accession counts across all parts.
    struct TaxonomySummary {
        size_t n_nodes_union  = 0;  // unique taxids across parts
        size_t n_accessions   = 0;  // sum across parts
        size_t n_parts_with_taxonomy = 0;
    };
    TaxonomySummary taxonomy_summary() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// Detects layout, opens accordingly. Throws std::runtime_error on failure.
ArchiveSetReader open_archive_auto(const std::filesystem::path& path);

} // namespace genopack
