#pragma once
#include "types.hpp"
#include "catalog.hpp"
#include "fmhr.hpp"
#include "core_section.hpp"
#include "colstore.hpp"
#include "bprm.hpp"
#include "gcov.hpp"
#include "gstx.hpp"
#include "pcore.hpp"
#include "prof.hpp"
#include "qcontig.hpp"
#include "qual.hpp"
#include "shard.hpp"
#include "skch.hpp"
#include "txdb.hpp"
#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace genopack {

// ── ArchiveReader ─────────────────────────────────────────────────────────────

class ArchiveReader {
public:
    ArchiveReader();
    ~ArchiveReader();
    ArchiveReader(const ArchiveReader&) = delete;
    ArchiveReader& operator=(const ArchiveReader&) = delete;

    void open(const std::filesystem::path& archive_dir);
    void close();
    bool is_open() const;

    // Metadata without sequence decompression
    std::optional<GenomeMeta> genome_meta(GenomeId id) const;

    // Count genomes matching query (no decompression)
    size_t count(const ExtractQuery& q) const;

    // Metadata-only filter result
    std::vector<GenomeMeta> filter_meta(const ExtractQuery& q) const;

    // Full extraction with FASTA decompression
    std::vector<ExtractedGenome> extract(const ExtractQuery& q) const;

    // Fetch one genome by id
    std::optional<ExtractedGenome> fetch_genome(GenomeId id) const;
    std::optional<std::string> fetch_sequence_slice(GenomeId id,
                                                    uint64_t start,
                                                    uint64_t length) const;

    // Fetch metadata only (no FASTA decompression) — O(1) lookup via ACCX+CATL
    std::optional<GenomeMeta> genome_meta_by_accession(std::string_view accession) const;

    // Fetch one genome by accession (uses ACCX section)
    std::optional<ExtractedGenome> fetch_by_accession(std::string_view accession) const;

    // Batch fetch: groups accessions by shard, reads each shard once, then
    // calls MADV_DONTNEED on it. Dramatically faster than N individual fetches
    // when genomes span many shards (avoids repeated NFS page faults).
    // results[i] = nullopt if accession not found.
    std::vector<std::optional<ExtractedGenome>>
    batch_fetch_by_accessions(const std::vector<std::string>& accessions) const;

    // Streaming shard visitor: calls cb(original_index, genome) for each found genome,
    // grouped by shard (each shard is read once then its pages released).
    // Peak memory: O(max_genomes_per_shard) instead of O(total_genomes).
    // Use this instead of batch_fetch_by_accessions when n_genomes is large.
    void visit_by_shard(
        const std::vector<std::string>& accessions,
        const std::function<void(size_t idx, ExtractedGenome)>& cb) const;

    // Shard-batch visitor: delivers all genomes from one shard at a time as a vector.
    // Allows the caller to parallelize embedding within each shard batch using OMP/threads.
    // Peak memory: O(max_genomes_per_shard). Shard pages released after each callback.
    using ShardBatch = std::vector<std::pair<size_t, ExtractedGenome>>; // (idx, genome)
    void visit_shard_batches(
        const std::vector<std::string>& accessions,
        const std::function<void(ShardBatch&)>& cb) const;
    // Like visit_shard_batches but uses n_readers independent fds for parallel NFS I/O.
    // cb may be called concurrently from multiple threads; caller must be thread-safe.
    void visit_shard_batches_parallel(
        const std::vector<std::string>& accessions,
        int n_readers,
        const std::function<void(ShardBatch&)>& cb) const;
    std::optional<std::string> fetch_sequence_slice_by_accession(std::string_view accession,
                                                                 uint64_t start,
                                                                 uint64_t length) const;

    // Contig accession → genome_id lookup (uses CIDX section; UINT32_MAX if not found)
    uint32_t find_contig_genome_id(std::string_view contig_acc) const;

    // High-performance batch lookup: maps N contig accessions → genome_ids.
    // out_genome_ids[i] = genome_id for accs[i], or UINT32_MAX if not found.
    // n_threads=0 → hardware_concurrency()
    void batch_find_contig_genome_ids(const std::string_view* accs,
                                      uint32_t*               out_genome_ids,
                                      size_t                  n,
                                      size_t                  n_threads = 0) const;

    // Taxonomy lookup by accession (uses TAXN section; nullopt if not present)
    std::optional<std::string> taxonomy_for_accession(std::string_view accession) const;

    // Iterate all (accession, taxonomy) pairs from TAXN section
    void scan_taxonomy(const std::function<void(std::string_view accession,
                                                std::string_view taxonomy)>& cb) const;

    // Reverse lookup: genome_id → accession (from ACCX section; empty if not found)
    std::string accession_for_genome_id(GenomeId id) const;

    // Iterate all (accession, genome_id) pairs from ACCX sections
    void scan_genome_accessions(const std::function<void(std::string_view, GenomeId)>& cb) const;

    // True if `id` has been tombstoned (soft-deleted). The raw ACCX/TAXN/QUAL
    // scans do NOT filter on this, so live-only consumers must check it.
    bool is_deleted(GenomeId id) const;

    // Parsed taxonomy tree from TXDB section (or auto-built from TAXN)
    std::optional<TaxonomyTree> taxonomy_tree() const;

    // Returns pointer to float[136] L2-normalised k=4 profile, nullptr if not stored.
    const float* kmer_profile(GenomeId genome_id) const;
    const float* kmer_profile_by_accession(std::string_view accession) const;

    // OPH sketch access (SKCH section)
    bool has_sketches() const;
    // V4 archives always carry both sig1 (seed1) and sig2 (seed2); returns
    // false only if no SKCH section is present.
    bool has_sig2() const;
    // Returns the first available sketch for genome_id, regardless of parameters.
    std::optional<SketchResult> sketch_for(GenomeId genome_id) const;
    // Param-aware: finds the section where kmer_size==k and sketch_size>=sz, then
    // returns a (possibly sliced) result. Returns nullopt if no matching section
    // exists or the genome is not found.
    std::optional<SketchResult> sketch_for(GenomeId genome_id,
                                           uint32_t k, uint32_t sz) const;
    // k-mer size of the first SKCH section (0 if no sketches present).
    uint32_t sketch_kmer_size() const;
    // sketch_size of the first SKCH section (0 if no sketches present).
    uint32_t sketch_sketch_size() const;
    // All k-mer sizes available across all SKCH sections (sorted ascending).
    std::vector<uint32_t> available_sketch_kmer_sizes() const;

    // Batch lookup by sorted genome_ids. For V3 archives, decompresses only the frames
    // that contain the requested ids. sorted_ids must be sorted ascending.
    // cb(idx_in_sorted_ids, SketchResult) — pointers are valid only during the callback.
    using SketchCallback = SkchReader::SketchCallback;
    // num_threads=0 → omp_get_max_threads(). Pass 1 for sequential NFS access.
    void sketch_for_ids(const std::vector<GenomeId>& sorted_ids,
                        uint32_t k, uint32_t sz,
                        const SketchCallback& cb,
                        int num_threads = 0) const;

    // Fused multi-k: one sequential pass yields all k values per genome.
    // cb(row_idx, ki, result) — ki is index into available_sketch_kmer_sizes().
    using SketchCallbackMultiK = SkchReader::SketchCallbackMultiK;
    void sketch_for_ids_multi_k(const std::vector<GenomeId>& sorted_ids,
                                 uint32_t sz,
                                 const SketchCallbackMultiK& cb,
                                 int num_threads = 0) const;

    // Release decompressed SKCH buffers to free memory. Sections auto-reload on next sketch_for().
    void release_sketches() const;

    // Memory used by decompressed SKCH buffers.
    size_t sketch_memory_bytes() const;

    // Iterate all live shard sections in file order.
    // Callback receives: shard section data pointer, file offset, compressed size, shard_id.
    // Thread-safe: mmap is read-only; each invocation of cb gets independent data.
    void scan_shards(const std::function<void(const uint8_t* data,
                                              uint64_t offset,
                                              uint64_t compressed_size,
                                              uint32_t shard_id)>& cb) const;

    // GSTX: genus sketch stats index (consensus sig + p90 + TNF centroid per genus).
    // Built during archive construction; enables O(1) check queries without loading members.
    bool has_gstx() const;
    // Returns nullptr if genus not found or no GSTX section present.
    const GstxEntry* gstx_for_genus(std::string_view genus) const;
    // k-mer sizes stored in GSTX (same as sketch k-values); empty if no GSTX.
    std::vector<uint32_t> gstx_kmer_sizes() const;
    // Raw reader pointer (nullptr if no GSTX section); valid for archive lifetime.
    const GstxReader* gstx_reader() const;

    // GCOV: per-genus biological covariance eigenbasis + calibrated quantiles.
    bool has_gcov() const;
    // Returns nullptr if genus not found or no GCOV section present.
    const GcovEntry* gcov_for_genus(std::string_view genus) const;
    // Raw reader pointer (nullptr if no GCOV section); valid for archive lifetime.
    const GcovReader* gcov_reader() const;

    // FCOV: per-family biological covariance (same layout as GCOV, keyed by family hash).
    bool has_fcov() const;
    // Returns nullptr if family not found or no FCOV section present.
    const GcovEntry* fcov_for_family(std::string_view family) const;
    // Raw reader pointer (nullptr if no FCOV section); valid for archive lifetime.
    const GcovReader* fcov_reader() const;

    // FMHR: per-genus FracMinHash reference sketches (k=21, c=125).
    bool has_fmhr() const;
    // Returns invalid FmhrView (valid()==false) if genus not found or no FMHR section.
    FmhrView fmhr_for_genus(std::string_view genus) const;
    // Raw reader pointer (nullptr if no FMHR section); valid for archive lifetime.
    const FmhrReader* fmhr_reader() const;

    // SEC_CORE — per-genus prevalence cores (intrinsic completeness reference).
    bool has_core() const;
    CoreView core_for_genus(std::string_view genus) const;
    const CoreReader* core_reader() const;

    // SEC_FCORE — per-family prevalence cores. The genus-core fallback: when a query
    // genome's genus has no (or too sparse a) CORE, family-core coverage is a coarser
    // but more universally-available intrinsic completeness reference. Built post-hoc
    // by `genopack fcore`. core_for_family returns invalid (valid()==false) if absent.
    bool has_fcore() const;
    CoreView core_for_family(std::string_view family) const;
    const CoreReader* fcore_reader() const;

    // SEC_PCORE — unified prevalence-annotated per-genus aamer reference (dense: every
    // aamer + u8 prevalence). Supersedes CORE/FCORE; dense enough for small-contig
    // foreign detection. Built by `genopack pcore`. Invalid view if genus absent.
    bool has_pcore() const;
    PcoreView pcore_for_genus(std::string_view genus) const;
    PcoreView pcore_for_family(std::string_view family) const;
    const PcoreReader* pcore_reader() const;

    // BPRM: self-describing build parameters (one per archive).
    bool has_bprm() const;
    // Pointer to the build-params header (nullptr if no BPRM section).
    const BprmHeader* build_params() const;
    // O(1) cross-section consistency key; 0 if no BPRM section.
    uint64_t build_params_hash() const;

    // QUAL: per-genome quality scores (completeness, contamination, consistency).
    // Written during build; enables O(1) quality lookup without re-computing.
    bool has_qual() const;
    void scan_qual(const std::function<void(const QualRecord&)>& cb) const;
    // QCOL columnar reader (nullptr if only legacy SEC_QUAL or no quality). Lets a
    // consumer read columns that have no QualRecord field (e.g. AAMER_GENUS_CORE).
    const ColStoreReader* qcol_reader() const;

    // XQAL: external quality (CheckM2 / anvi'o), columnar, keyed by genome_id.
    // Ingested via `genopack ingest`; each tool's measurements are distinct columns
    // (tool != "genopack"). Reader pointer is nullptr if no XQAL section.
    bool has_xqual() const;
    const ColStoreReader* xqual_reader() const;

    // CQAL: calibrated / fused quality, columnar, keyed by genome_id. Written by
    // `genopack calibrate`; columns carry a calibration_hash. nullptr if no section.
    bool has_cqal() const;
    const ColStoreReader* cqal_reader() const;

    // PROF: named reporting/fusion profiles (user-authored), each pinning per-axis
    // column identities. Written by `genopack profile`; consumed by `genopack
    // report --profile <name>`. Reader pointer is nullptr if no PROF section.
    bool has_prof() const;
    const ProfReader* prof_reader() const;

    // QCONTIG: per-contig quality overlay (offset/length/TNF/leakage per contig),
    // keyed by genome_id. Written by `check`; lets a consumer see which contigs
    // drive a genome's contamination. Reader pointer is nullptr if no section.
    bool has_qcontig() const;
    const QcontigReader* qcontig_reader() const;

    // File descriptor for the underlying mmap (for posix_fadvise).
    int fd() const;

    // Archive-wide stats
    struct ArchiveStats {
        uint64_t generation;
        size_t   n_shards;
        size_t   n_genomes_total;
        size_t   n_genomes_live;
        uint64_t total_raw_bp;
        uint64_t total_compressed_bytes;
        double   compression_ratio;
    };
    ArchiveStats archive_stats() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// ── ArchiveBuilder ────────────────────────────────────────────────────────────
// Builds a new archive from scratch. Not thread-safe (single writer).

struct ArchiveBuilderConfig {
    ShardWriterConfig shard_cfg;
    size_t io_threads = 16;     // parallel FASTA decompression
    bool   verbose    = false;
    size_t batch_size = 1000;   // genomes buffered before shard flush
    GenomeId starting_genome_id = 1;  // first genome_id assigned; set per-part in parallel builds
    bool     build_cidx = true;       // build CIDX contig index; disable with --no-cidx for large archives
    bool     kmer_nn_sort        = true;   // sort genomes within each shard by kmer4_profile NN chain
                                           // instead of oph_fingerprint; improves zstd compression
    bool     taxonomy_group      = true;   // bucket genomes by taxonomy before shard formation
    std::string taxonomy_rank    = "g";    // "g"=genus, "f"=family; genus with family fallback
    bool     build_sketch        = false;  // compute OPH sketches and write SKCH section
    bool     build_gstx          = true;   // build GSTX genus-stats index (needs taxonomy_group + build_sketch)
    bool     build_gcov          = true;   // build GCOV per-genus covariance eigenbasis (needs build_gstx)
    bool     build_core          = true;   // build SEC_CORE per-genus prevalence cores (needs markers_path)
    bool     build_pcore         = false;  // also build SEC_PCORE dense per-genus reference (small-contig
                                           // contamination). OFF by default: the dense union is ~10-50x CORE
                                           // and held in memory until finalize — opt in per build (--pcore).
    float    core_theta          = 0.90f;  // prevalence threshold for genus-core membership
    uint32_t micro_genus_threshold = 0;    // min genus members for a dedicated shard + GSTX/GCOV/FMHR
                                           // model; smaller genera are bin-packed and unmodeled.
                                           // 0 = auto (scales with corpus: 4 if <=50k genomes, 8 if
                                           // <=500k, 16 if <=2M, else 32). Lower = more genera modeled,
                                           // more shards
    int      sketch_kmer_size    = 16;     // k-mer size for OPH sketching (single-k path)
    std::vector<int> sketch_kmer_sizes;   // if size > 1: write multi-k SKCH v2 (overrides sketch_kmer_size)
    int      sketch_size         = 10000;  // number of OPH bins
    int      sketch_syncmer_s    = 0;      // 0 = disabled; >0 = open syncmer prefilter
    uint64_t sketch_seed         = 42;
    int      fmh_k               = 21;     // FracMinHash k for FMHR reference sketches
    int      fmh_c               = 125;    // FracMinHash density (1/c) for FMHR sketches
    std::string markers_path;           // path to .mrk file; empty = skip build-time marker scoring
    int         marker_min_hits  = 5;  // min pool hits to call a marker present (matches check default)
    std::string contam_panel_path;      // path to .csp contamination panel; empty = skip duplication scoring
    std::filesystem::path from_gpk_source;  // non-empty: rebuild by streaming decoded sequence from this .gpk
    std::filesystem::path from_stage_source; // non-empty: rebuild by streaming decoded sequence from this .gstage cache
};

// `build --from-gpk` fast-path: if the source's build params equal what cfg
// would produce, repack it verbatim (no decompress/recompute) into output and
// return true. Returns false if params differ — caller does a streaming rebuild.
bool from_gpk_try_verbatim_reuse(const std::filesystem::path& source,
                                 const std::filesystem::path& output,
                                 const ArchiveBuilderConfig& cfg);

class ArchiveBuilder {
public:
    using Config = ArchiveBuilderConfig;

    explicit ArchiveBuilder(const std::filesystem::path& archive_dir,
                            Config cfg = Config{});
    ~ArchiveBuilder();

    // Add all genomes from a TSV (accession, file_path [, completeness, contamination])
    void add_from_tsv(const std::filesystem::path& tsv_path);

    // Rebuild from an existing .gpk: stream decoded sequence from its shards
    // (sequential reads) instead of opening per-genome FASTA files on NFS.
    void add_from_gpk(const std::filesystem::path& source);

    // Rebuild from a .gstage cache: stream decoded sequence sequentially from a
    // local store produced by `genopack stage` (no per-genome NFS file opens).
    void add_from_stage(const std::filesystem::path& source);

    // Add individual record
    void add(const BuildRecord& rec);

    // Write all shards, catalog, manifest, meta.tsv sidecar
    void finalize();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// ── ArchiveAppender ───────────────────────────────────────────────────────────
// Append genomes to an existing archive (new shard generation).

class ArchiveAppender {
public:
    explicit ArchiveAppender(const std::filesystem::path& archive_dir);
    ~ArchiveAppender();

    void add_from_tsv(const std::filesystem::path& tsv_path);
    void add(const BuildRecord& rec);

    // Tombstone genomes (mark deleted; does not remove from shards)
    void remove(GenomeId id);
    void remove_by_accession(std::string_view accession);

    // Write new shard(s), update catalog, increment generation
    void commit();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace genopack
