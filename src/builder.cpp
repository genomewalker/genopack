#include <genopack/archive.hpp>
#include <genopack/aamer_counter.hpp>
#include <genopack/accx.hpp>
#include <genopack/oph_sketch.hpp>
#include <genopack/skch.hpp>
#include <genopack/cidx.hpp>
#include <genopack/gidx.hpp>
#include <genopack/gstx.hpp>
#include <genopack/build_genus_stats.hpp>
#include <genopack/skani.hpp>
#include <genopack/markers.hpp>
#include <genopack/aamer.hpp>
#include <genopack/cladesplit.hpp>
#include <genopack/core_section.hpp>
#include <genopack/pcore_spiller.hpp>
#include <genopack/qual.hpp>
#include <genopack/qual_columns.hpp>
#include <genopack/kmrx.hpp>
#include <genopack/taxn.hpp>
#include <genopack/txdb.hpp>
#include <genopack/catalog.hpp>
#include <genopack/format.hpp>
#include <genopack/bprm.hpp>
#include <genopack/derivation.hpp>
#include <genopack/section_checksum.hpp>
#include <genopack/mmap_file.hpp>

#include <genopack/shard.hpp>
#include <genopack/toc.hpp>
#include <genopack/util.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cassert>
#include <condition_variable>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <future>
#include <mutex>
#include <chrono>
#include <numeric>
#include <optional>
#include <queue>
#include <random>
#include <stdexcept>
#include <span>
#include <string_view>
#include <thread>
#include <unistd.h>
#include <fcntl.h>
#include <sys/sendfile.h>
#include <genopack/checksum.hpp>
#include <unordered_map>
#include <vector>

namespace genopack {

namespace {
// LSD byte-radix sort + dedup for uint64_t. O(8n), branch-free inner loop.
// Output is ascending with duplicates removed — bit-identical to std::sort+unique.
// `scratch` is a caller-supplied reusable buffer (use thread_local to amortise alloc).
static void radix_sort_dedup_u64(std::vector<uint64_t>& v,
                                 std::vector<uint64_t>& scratch) {
    const size_t n = v.size();
    if (n < 2) return;
    scratch.resize(n);
    uint64_t* src = v.data();
    uint64_t* dst = scratch.data();
    for (int shift = 0; shift < 64; shift += 8) {
        size_t cnt[257] = {};
        for (size_t i = 0; i < n; ++i) ++cnt[((src[i] >> shift) & 0xFF) + 1];
        bool single = false;
        for (int b = 0; b < 256; ++b) if (cnt[b + 1] == n) { single = true; break; }
        if (single) continue;
        for (int b = 0; b < 256; ++b) cnt[b + 1] += cnt[b];
        for (size_t i = 0; i < n; ++i) dst[cnt[(src[i] >> shift) & 0xFF]++] = src[i];
        std::swap(src, dst);
    }
    if (src != v.data()) std::memcpy(v.data(), src, n * sizeof(uint64_t));
    size_t w = 1;
    for (size_t r = 1; r < n; ++r) if (v[r] != v[r - 1]) v[w++] = v[r];
    v.resize(w);
}
} // namespace

// ── Build-params (BPRM) ───────────────────────────────────────────────────────
// Build the canonical BprmHeader from the config. Pure function of cfg, so it is
// computed once at header-write time (for FileHeader.build_params_hash) and again
// at section-write time (for the BPRM body) — guaranteeing the two agree.
// Resolve the micro-genus threshold: genera with fewer members than this are
// bin-packed into combined shards and get no GSTX/GCOV/FMHR consensus model.
// configured==0 means auto — scale with corpus size so small corpora model
// nearly every genus while huge corpora bound shard fragmentation.
static uint32_t resolve_micro_threshold(uint32_t configured, uint64_t n_genomes) {
    if (configured != 0) return configured;
    if (n_genomes <=   50000) return 4;
    if (n_genomes <=  500000) return 8;
    if (n_genomes <= 2000000) return 16;
    return 32;
}

static BprmHeader make_bprm_header_from_cfg(const ArchiveBuilderConfig& cfg) {
    BprmHeader bp{};
    bp.sketch_kmer_size = static_cast<uint32_t>(cfg.sketch_kmer_size);
    bp.sketch_size      = static_cast<uint32_t>(cfg.sketch_size);
    bp.sketch_syncmer_s = static_cast<uint32_t>(cfg.sketch_syncmer_s);
    const uint32_t nk = static_cast<uint32_t>(
        std::min<size_t>(cfg.sketch_kmer_sizes.size(), 8));
    bp.n_kmer_sizes = nk;
    for (uint32_t i = 0; i < nk; ++i)
        bp.kmer_sizes[i] = static_cast<uint32_t>(cfg.sketch_kmer_sizes[i]);
    bp.sketch_seed        = cfg.sketch_seed;
    bp.fmh_k              = static_cast<uint32_t>(cfg.fmh_k);
    bp.fmh_c              = static_cast<uint32_t>(cfg.fmh_c);
    bp.tnf_order          = TNF_ORDER;
    bp.gstx_model_version = GSTX_MODEL_VERSION;
    bp.gcov_model_version = GCOV_MODEL_VERSION;
    bp.taxonomy_rank      = cfg.taxonomy_rank.empty()
        ? static_cast<uint8_t>('g')
        : static_cast<uint8_t>(cfg.taxonomy_rank[0]);
    uint32_t bf = 0;
    if (cfg.build_cidx)            bf |= BPRM_F_CIDX;
    if (cfg.build_sketch)          bf |= BPRM_F_SKETCH;
    if (cfg.build_gstx)            bf |= BPRM_F_GSTX;
    if (cfg.build_gcov)            bf |= BPRM_F_GCOV;
    if (cfg.taxonomy_group)        bf |= BPRM_F_TAXGROUP;
    if (!cfg.markers_path.empty()) bf |= BPRM_F_MARKERS;
    if (cfg.thin_archive)          bf |= BPRM_F_THIN;
    bp.build_flags = bf;
    bp.micro_genus_threshold = cfg.micro_genus_threshold;
    bp.core_theta            = cfg.build_pcore ? cfg.core_theta : 0.0f;
    bp.magic       = SEC_BPRM;
    bp.version     = 1;
    bp.header_size = sizeof(BprmHeader);
    bp.params_hash = bprm_compute_params_hash(bp);
    return bp;
}

// ── --from-gpk verbatim reuse fast-path (§5.3, the reuse oracle) ──────────────
// If the source archive's build_params_hash equals the params hash cfg would
// produce, params + section set are identical, so the whole archive can be
// repacked VERBATIM — no decompress, no recompute. Corpus/content identity is
// NOT proven by the hash (it is corpus-blind, see bprm.hpp); it holds here only
// because this path repacks a SINGLE source archive and copies its own bodies,
// so every derived section travels with the content it was derived from. Do NOT
// reuse this oracle to justify copying derived sections across a corpus-changing
// op (e.g. merge), where source and destination corpora differ. Returns false
// when params differ; the caller then falls back to a full streaming rebuild.
bool from_gpk_try_verbatim_reuse(const std::filesystem::path& source,
                                 const std::filesystem::path& output,
                                 const ArchiveBuilderConfig& cfg) {
    MmapFileReader src;
    src.open(source);
    if (src.size() < sizeof(FileHeader)) return false;
    const auto* sfh = src.ptr_at<FileHeader>(0);
    if (sfh->magic != GPK_MAGIC) return false;
    // The source must record its build params for an equivalence proof.
    if (sfh->build_params_hash == 0) return false;
    const uint64_t new_params = make_bprm_header_from_cfg(cfg).params_hash;
    if (sfh->build_params_hash != new_params) return false;   // params/section-set changed

    Toc toc = TocReader::read(src);

    std::filesystem::path partial = output;
    partial += ".partial";
    if (!output.parent_path().empty())
        std::filesystem::create_directories(output.parent_path());

    std::random_device rd;
    const uint64_t new_uuid_lo = (static_cast<uint64_t>(rd()) << 32) | rd();
    const uint64_t new_uuid_hi = (static_cast<uint64_t>(rd()) << 32) | rd();
    const uint64_t new_gen     = sfh->generation + 1;

    AppendWriter w;
    w.create(partial);
    {
        FileHeader fh{};
        fh.magic             = GPK_MAGIC;
        fh.version_major     = FORMAT_MAJOR;
        fh.version_minor     = FORMAT_MINOR;
        fh.file_uuid_lo      = new_uuid_lo;
        fh.file_uuid_hi      = new_uuid_hi;
        fh.created_at_unix   = static_cast<uint64_t>(std::time(nullptr));
        fh.endian_abi_tag    = ENDIAN_ABI_TAG;
        fh.build_params_hash = sfh->build_params_hash;
        fh.generation        = new_gen;
        fh.shard_set_uuid_lo = sfh->shard_set_uuid_lo;
        fh.shard_set_uuid_hi = sfh->shard_set_uuid_hi;
        fh.shard_id          = sfh->shard_id;
        w.append(&fh, sizeof(fh));
    }

    // Copy every section body verbatim; bodies use section-relative offsets, so a
    // new 8-aligned start is valid (the directory is the only absolute reference,
    // and it is rewritten below). content_xxh128 / derivation_hash /
    // semantic_schema_hash carry over unchanged because the bytes are unchanged.
    TocWriter tw;
    for (const auto& sd : toc.sections) {
        if (sd.compressed_size == 0) { tw.add_section(sd); continue; }
        if (sd.file_offset > src.size() ||
            sd.compressed_size > src.size() - sd.file_offset)
            throw std::runtime_error("from-gpk reuse: section out of bounds in source");
        w.align(8);
        SectionDesc nd = sd;
        nd.file_offset = w.current_offset();
        w.append(src.data() + sd.file_offset, sd.compressed_size);
        tw.add_section(nd);
    }

    tw.finalize(w, new_gen,
                toc.header.live_genome_count, toc.header.total_genome_count,
                /*prev_toc_offset=*/0,
                toc.header.catalog_root_section_id,
                toc.header.accession_root_section_id,
                toc.header.tombstone_root_section_id,
                new_uuid_lo, new_uuid_hi);
    w.flush();
    w.close();
    std::filesystem::rename(partial, output);
    if (!output.parent_path().empty()) {
        int dfd = ::open(output.parent_path().c_str(), O_RDONLY | O_DIRECTORY);
        if (dfd >= 0) { ::fsync(dfd); ::close(dfd); }
    }
    spdlog::info("from-gpk: params unchanged — verbatim repack of {} section(s), no recompute",
                 toc.sections.size());
    return true;
}

// ── Shard grouping helpers ────────────────────────────────────────────────────

// Extract taxonomy bucket key at the requested rank ('g'=genus, 'f'=family).
// Falls back to the next coarser rank if the requested rank is absent.
// Returns "__unclassified__" if taxonomy is empty.
static std::string extract_taxonomy_bucket(std::string_view taxonomy,
                                           char rank)
{
    if (taxonomy.empty()) return "__unclassified__";
    // Taxonomy format: "d__;p__;c__;o__;f__;g__;s__"
    // Each rank prefix: "g__", "f__", etc.
    char prefix[4] = {rank, '_', '_', ';'};
    auto pos = taxonomy.find(std::string_view(prefix, 3));
    if (pos == std::string_view::npos) {
        // Fallback: try family then order
        const char* fallbacks = (rank == 'g') ? "foc" : "oc";
        for (char fb : std::string_view(fallbacks)) {
            char fp[4] = {fb, '_', '_', ';'};
            pos = taxonomy.find(std::string_view(fp, 3));
            if (pos != std::string_view::npos) break;
        }
        if (pos == std::string_view::npos) return "__unclassified__";
    }
    size_t start = pos;
    size_t end   = taxonomy.find(';', start + 3);
    if (end == std::string_view::npos) end = taxonomy.size();
    if (end - start <= 3) return "__unclassified__";  // prefix only, no name
    return std::string(taxonomy.substr(start, end - start));
}

// Approximate nearest-neighbor ordering over 136-dim kmer4_profiles.
// O(n log n): projects each feature vector onto a fixed random direction then
// sorts by that scalar.  Deterministic (fixed seed — no std::random_device).
template<typename Item>
static std::vector<size_t> greedy_nn_chain(const std::vector<Item>& items)
{
    const size_t n = items.size();
    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    if (n <= 1) return idx;

    // Build a fixed random projection vector (136-dim).
    std::mt19937_64 rng(0x5eed5eed5eed5eedULL);
    std::uniform_real_distribution<float> ud(-1.f, 1.f);
    std::array<float, 136> proj{};
    for (int d = 0; d < 136; ++d) proj[d] = ud(rng);

    // Project each genome and sort indices by scalar projection value.
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
        float pa = 0.f, pb = 0.f;
        const float* va = items[a].stats.kmer4_profile.data();
        const float* vb = items[b].stats.kmer4_profile.data();
        for (int d = 0; d < 136; ++d) { pa += va[d] * proj[d]; pb += vb[d] * proj[d]; }
        return pa < pb;
    });
    return idx;
}

// ── Checkpoint binary format ──────────────────────────────────────────────────
//
// {gpk}.ckpt         — text key=value (updated atomically after each shard):
//   genome_count=N   total genomes written so far
//   byte_offset=B    file offset after last written shard
//   next_genome_id=G next genome_id to assign
//   current_shard_id=S next shard_id to open
//   next_section_id=I next section_id to assign
//
// {gpk}.ckpt_meta.bin — binary stream of shard blocks:
//   [ShardCkptHdr] [GenomeCkptFixed + accession bytes + taxonomy bytes] × n_genomes
//
// On resume: restore all in-memory state, truncate .gpk to byte_offset,
// open in append mode, skip first genome_count TSV rows.

#pragma pack(push, 1)
struct ShardCkptHdr {
    uint32_t shard_id;
    uint32_t n_genomes;
    uint64_t section_id;
    uint64_t file_offset;
    uint64_t compressed_size;
};  // 28 bytes

struct GenomeCkptFixed {
    uint64_t genome_id;
    uint64_t input_row_index;
    uint32_t shard_id;
    uint32_t dir_index;
    uint64_t oph_fingerprint;
    uint64_t genome_length;
    uint32_t n_contigs;
    uint32_t gc_pct_x100;
    uint16_t completeness_x10;
    uint16_t contamination_x10;
    uint32_t date_added;
    uint16_t accession_len;
    uint16_t taxonomy_len;
    uint32_t genome_type;   // GenomeType enum
    uint32_t _ckpt_pad0;
    float    kmer4[136];
};  // 612 bytes
#pragma pack(pop)

// ── ArchiveBuilder::Impl ──────────────────────────────────────────────────────

struct ArchiveBuilder::Impl {
    std::filesystem::path archive_dir;   // base path (no extension)
    std::filesystem::path gpk_path_;     // output .gpk file
    std::unique_ptr<ArchiveReader> src_reader_;  // from-gpk source (streams decoded sequence)
    Config                cfg;

    std::vector<BuildRecord> pending;
    std::vector<uint64_t>    pending_input_rows;
    GenomeId next_genome_id = 1;  // overridden by cfg.starting_genome_id below

    // Pass-1 result: metadata only, FASTA discarded after stats computation
    struct GenomeMeta1 {
        BuildRecord record;
        GenomeId    genome_id;
        FastaStats  stats;
    };

    explicit Impl(const std::filesystem::path& dir, Config c)
        : archive_dir(dir), cfg(c), next_genome_id(c.starting_genome_id)
    {
        // Determine output .gpk path
        if (dir.extension() == ".gpk")
            gpk_path_ = dir;
        else
            gpk_path_ = std::filesystem::path(dir.string() + ".gpk");

        // Ensure parent directory exists (for meta.tsv sidecar). A relative
        // output like "out.gpk" has an empty parent_path() — create_directories("")
        // throws "Invalid argument", so guard it.
        if (!gpk_path_.parent_path().empty())
            std::filesystem::create_directories(gpk_path_.parent_path());
    }

    void add(const BuildRecord& rec) {
        pending.push_back(rec);
        pending_input_rows.push_back(static_cast<uint64_t>(pending_input_rows.size()));
    }

    void add_from_gpk(const std::filesystem::path& source) {
        cfg.from_gpk_source = source;
        src_reader_ = std::make_unique<ArchiveReader>();
        src_reader_->open(source);
        // Collect accession -> taxonomy from the source TAXN section.
        std::unordered_map<std::string, std::string> tax;
        src_reader_->scan_taxonomy([&](std::string_view acc, std::string_view t) {
            tax.emplace(std::string(acc), std::string(t));
        });
        // One BuildRecord per source genome; the sequence is streamed in finalize().
        src_reader_->scan_genome_accessions([&](std::string_view acc, GenomeId) {
            BuildRecord r;
            r.accession = std::string(acc);
            auto it = tax.find(r.accession);
            if (it != tax.end())
                r.extra_fields.emplace_back("taxonomy", it->second);
            add(r);
        });
        spdlog::info("from-gpk: {} genomes from {}", pending.size(), source.string());
    }

    void add_from_tsv(const std::filesystem::path& tsv_path) {
        auto records = parse_tsv_records(tsv_path);
        uint64_t base = pending_input_rows.size();
        for (auto& r : records)
            pending.push_back(std::move(r));
        for (size_t i = 0; i < records.size(); ++i)
            pending_input_rows.push_back(base + i);
        spdlog::info("Loaded {} records from {}", pending.size(), tsv_path.string());
    }

    void finalize() {
        if (pending.empty()) { spdlog::warn("No records to build"); return; }

        const size_t original_total_records = pending.size();
        size_t total_records = pending.size();
        spdlog::info("Building archive: {} genomes, {} threads", total_records, cfg.io_threads);

        // Pre-scan per-genus member totals so CORE/PCORE aggregate ACROSS the 512MB
        // shard-chunks a mega-genus is split into, emitting one full-density model
        // per genus when all members have arrived (genus-contiguous input ⇒ bounded
        // active combiners). Only modeled genera (≥ micro threshold) participate.
        std::unordered_map<std::string, uint32_t> genus_total;
        uint32_t cc_micro_threshold = 0;
        if (cfg.taxonomy_group && cfg.build_pcore) {
            const char grank = cfg.taxonomy_rank.empty() ? 'g' : cfg.taxonomy_rank[0];
            for (const auto& rec : pending)
                for (const auto& [k, v] : rec.extra_fields)
                    if (k == "taxonomy") { ++genus_total[extract_taxonomy_bucket(v, grank)]; break; }
            cc_micro_threshold = resolve_micro_threshold(cfg.micro_genus_threshold, original_total_records);
        }

        // ── Checkpoint paths ──────────────────────────────────────────────────
        std::filesystem::path ckpt_path      = std::filesystem::path(gpk_path_.string() + ".ckpt");
        std::filesystem::path ckpt_meta_path = std::filesystem::path(gpk_path_.string() + ".ckpt_meta.bin");
        // Build into a .partial sidecar and atomically rename to the final path
        // only after a durable footer. A crash leaves .partial + checkpoints for
        // resume; the final .gpk only ever appears fully sealed.
        std::filesystem::path partial_path   = std::filesystem::path(gpk_path_.string() + ".partial");

        AppendWriter app_writer;
        TocWriter    toc;
        uint64_t next_section_id = 1;
        ShardId  current_shard_id = 0;

        std::vector<GenomeMeta> catalog_rows;
        catalog_rows.reserve(total_records);
        std::vector<std::pair<std::string, GenomeId>> accession_pairs;
        accession_pairs.reserve(total_records);
        std::unordered_map<std::string, std::string> taxonomy_map;
        taxonomy_map.reserve(total_records);
        std::vector<std::pair<GenomeId, std::array<float, 136>>> kmer_pairs;
        kmer_pairs.reserve(total_records);

        // OPH sketch writer: spills sigs+masks to tmpfile immediately, no RAM accumulation
        std::unique_ptr<SkchWriter>       skch_writer;
        std::unique_ptr<SkchWriterMultiK> skch_writer_mk;
        const bool multi_k_sketch = cfg.build_sketch && cfg.sketch_kmer_sizes.size() > 1;
        if (cfg.build_sketch) {
            if (multi_k_sketch) {
                std::vector<uint32_t> ks;
                for (int k : cfg.sketch_kmer_sizes) ks.push_back(static_cast<uint32_t>(k));
                skch_writer_mk = std::make_unique<SkchWriterMultiK>(
                    ks, static_cast<uint32_t>(cfg.sketch_size),
                    static_cast<uint32_t>(cfg.sketch_syncmer_s),
                    cfg.sketch_seed, cfg.sketch_seed + 1,
                    cfg.tmpdir.empty() ? gpk_path_.parent_path().string() : cfg.tmpdir);
            } else {
                skch_writer = std::make_unique<SkchWriter>(
                    static_cast<uint32_t>(cfg.sketch_size),
                    static_cast<uint32_t>(cfg.sketch_kmer_size),
                    static_cast<uint32_t>(cfg.sketch_syncmer_s),
                    cfg.sketch_seed, cfg.sketch_seed + 1,
                    cfg.tmpdir.empty() ? gpk_path_.parent_path().string() : cfg.tmpdir);
            }
        }

        std::vector<uint64_t> source_row_indices;
        source_row_indices.reserve(original_total_records);

        struct GidxInfo { GenomeId genome_id; ShardId shard_id; uint32_t dir_index; };
        std::vector<GidxInfo> gidx_infos;
        gidx_infos.reserve(total_records);
        std::unordered_map<ShardId, uint64_t> shard_id_to_section_id;
        CidxWriter  cidx_writer;
        GstxWriter  gstx_writer;
        QualWriter  qual_writer;

        // ── Try resume ───────────────────────────────────────────────────────
        bool     resuming      = false;

        if (std::filesystem::exists(ckpt_path) && std::filesystem::exists(ckpt_meta_path)
            && std::filesystem::exists(partial_path)) {

            size_t   ck_genome_count  = 0;
            uint64_t ck_byte_offset   = 0;
            uint64_t ck_next_gid      = 0;
            uint32_t ck_shard_id      = 0;
            uint64_t ck_section_id    = 0;

            // Parse .ckpt
            std::ifstream ck(ckpt_path);
            std::string line;
            bool parsed_ok = false;
            int parsed_fields = 0;
            while (std::getline(ck, line)) {
                auto eq = line.find('=');
                if (eq == std::string::npos) continue;
                std::string key = line.substr(0, eq);
                std::string val = line.substr(eq + 1);
                if      (key == "genome_count")    { ck_genome_count  = std::stoull(val); ++parsed_fields; }
                else if (key == "byte_offset")     { ck_byte_offset   = std::stoull(val); ++parsed_fields; }
                else if (key == "next_genome_id")  { ck_next_gid      = std::stoull(val); ++parsed_fields; }
                else if (key == "current_shard_id"){ ck_shard_id      = static_cast<uint32_t>(std::stoull(val)); ++parsed_fields; }
                else if (key == "next_section_id") { ck_section_id    = std::stoull(val); ++parsed_fields; }
            }
            parsed_ok = (parsed_fields == 5);

            // Read .ckpt_meta.bin and restore all in-memory state
            // ck_genome_count == total_records is valid: all SHRDs written, finalization interrupted
            if (parsed_ok && ck_genome_count > 0 && ck_genome_count <= total_records) {
                FILE* mf = std::fopen(ckpt_meta_path.c_str(), "rb");
                bool meta_ok = (mf != nullptr);
                size_t genomes_restored = 0;
                std::vector<char> processed_rows(total_records, 0);
                // Restored genomes grouped by their shard section, so their sketches can be
                // rebuilt from the archive's stored sequences (the interrupted run's SKCH spill
                // is gone). Without this they ship resolvable-but-sketchless.
                std::vector<std::pair<SectionDesc, std::vector<GenomeId>>> restored_shards;

                while (meta_ok && genomes_restored < ck_genome_count) {
                    ShardCkptHdr shdr{};
                    if (std::fread(&shdr, sizeof(shdr), 1, mf) != 1) { meta_ok = false; break; }

                    // Reconstruct TOC entry for this shard
                    SectionDesc sd{};
                    sd.type             = SEC_SHRD;
                    sd.version          = 4;
                    sd.flags            = 0;
                    sd.section_id       = shdr.section_id;
                    sd.file_offset      = shdr.file_offset;
                    sd.compressed_size  = shdr.compressed_size;
                    sd.uncompressed_size = 0;
                    sd.item_count       = shdr.n_genomes;
                    sd.aux0             = shdr.shard_id;
                    sd.aux1             = 0;
                    std::memset(sd.checksum, 0, sizeof(sd.checksum));
                    toc.add_section(sd);
                    if (cfg.build_sketch) restored_shards.push_back({sd, {}});
                    shard_id_to_section_id[shdr.shard_id] = shdr.section_id;

                    for (uint32_t gi = 0; gi < shdr.n_genomes && meta_ok; ++gi) {
                        GenomeCkptFixed gf{};
                        if (std::fread(&gf, sizeof(gf), 1, mf) != 1) { meta_ok = false; break; }
                        if (gf.input_row_index >= processed_rows.size()) { meta_ok = false; break; }
                        processed_rows[gf.input_row_index] = 1;

                        std::string accession(gf.accession_len, '\0');
                        std::string taxonomy(gf.taxonomy_len, '\0');
                        if (gf.accession_len &&
                            std::fread(accession.data(), 1, gf.accession_len, mf) != gf.accession_len) {
                            meta_ok = false; break;
                        }
                        if (gf.taxonomy_len &&
                            std::fread(taxonomy.data(), 1, gf.taxonomy_len, mf) != gf.taxonomy_len) {
                            meta_ok = false; break;
                        }

                        GenomeMeta meta{};
                        meta.genome_id         = gf.genome_id;
                        meta.genome_type       = gf.genome_type;
                        meta.shard_id          = gf.shard_id;
                        meta.genome_length     = gf.genome_length;
                        meta.n_contigs         = gf.n_contigs;
                        meta.gc_pct_x100       = static_cast<uint16_t>(gf.gc_pct_x100);
                        meta.completeness_x10  = gf.completeness_x10;
                        meta.contamination_x10 = gf.contamination_x10;
                        meta.oph_fingerprint   = gf.oph_fingerprint;
                        meta.date_added        = gf.date_added;
                        catalog_rows.push_back(meta);
                        if (cfg.build_sketch) restored_shards.back().second.push_back(gf.genome_id);

                        gidx_infos.push_back({gf.genome_id,
                                              static_cast<ShardId>(gf.shard_id),
                                              gf.dir_index});

                        accession_pairs.emplace_back(accession, gf.genome_id);

                        std::array<float, 136> kmer4;
                        std::memcpy(kmer4.data(), gf.kmer4, sizeof(gf.kmer4));
                        kmer_pairs.emplace_back(gf.genome_id, kmer4);
                        source_row_indices.push_back(gf.input_row_index);

                        if (gf.taxonomy_len)
                            taxonomy_map.emplace(accession, taxonomy);

                        ++genomes_restored;
                    }
                }
                if (mf) std::fclose(mf);

                if (meta_ok && genomes_restored == ck_genome_count) {
                    // Rebuild sketches for restored genomes from the archive's stored sequences,
                    // matching the main build's params (k, size, syncmer s, seeds). Their catalog
                    // rows exist but their SKCH rows were lost with the interrupted run's spill.
                    if (cfg.build_sketch) {
                        MmapFileReader pmm;
                        pmm.open(partial_path);
                        const uint64_t seed1 = cfg.sketch_seed, seed2 = cfg.sketch_seed + 1;
                        for (const auto& [sd, gids] : restored_shards) {
                            if (gids.empty()) continue;
                            ShardReader shard;
                            shard.open(pmm.data(), sd.file_offset, sd.compressed_size);
                            const int ng = static_cast<int>(gids.size());
                            std::vector<std::vector<OPHDualSketchResult>> mk(ng);
                            std::vector<OPHDualSketchResult>              s1(ng);
                            #pragma omp parallel for schedule(dynamic, 1)
                            for (int j = 0; j < ng; ++j) {
                                std::string fa = shard.fetch_genome(gids[j]);
                                if (fa.empty()) continue;
                                if (multi_k_sketch)
                                    mk[j] = sketch_oph_dual_multik(fa.data(), fa.size(),
                                        std::span<const int>(cfg.sketch_kmer_sizes),
                                        cfg.sketch_size, cfg.sketch_syncmer_s, seed1, seed2);
                                else
                                    s1[j] = sketch_oph_dual_from_buffer(fa.data(), fa.size(),
                                        cfg.sketch_kmer_size, cfg.sketch_size,
                                        cfg.sketch_syncmer_s, seed1, seed2);
                            }
                            for (int j = 0; j < ng; ++j) {
                                if (multi_k_sketch) {
                                    auto& sks = mk[j];
                                    if (sks.empty())
                                        throw std::runtime_error("resume re-sketch failed for genome "
                                                                 + std::to_string(gids[j]));
                                    std::vector<std::vector<uint16_t>> g1, g2;
                                    std::vector<std::vector<uint64_t>> gm;
                                    std::vector<uint32_t>              gn;
                                    for (auto& s : sks) {
                                        g1.push_back(std::move(s.signature1));
                                        g2.push_back(std::move(s.signature2));
                                        gm.push_back(std::move(s.real_bins_bitmask));
                                        gn.push_back(s.n_real_bins);
                                    }
                                    skch_writer_mk->add(gids[j], sks[0].genome_length, g1, g2, gn, gm);
                                } else {
                                    auto& s = s1[j];
                                    if (s.signature1.empty())
                                        throw std::runtime_error("resume re-sketch failed for genome "
                                                                 + std::to_string(gids[j]));
                                    skch_writer->add(gids[j], s.signature1, s.signature2,
                                                     s.n_real_bins, s.genome_length, s.real_bins_bitmask);
                                }
                            }
                        }
                    }
                    // Truncate .gpk to the checkpoint byte offset
                    if (::truncate(partial_path.c_str(), static_cast<off_t>(ck_byte_offset)) != 0)
                        throw std::runtime_error("checkpoint resume: truncate failed: " +
                                                 std::string(std::strerror(errno)));
                    app_writer.open_append(partial_path);
                    app_writer.seek_to(ck_byte_offset);

                    next_genome_id    = static_cast<GenomeId>(ck_next_gid);
                    current_shard_id  = static_cast<ShardId>(ck_shard_id);
                    next_section_id   = ck_section_id;

                    std::vector<BuildRecord> remaining_pending;
                    std::vector<uint64_t> remaining_rows;
                    remaining_pending.reserve(pending.size() - genomes_restored);
                    remaining_rows.reserve(pending_input_rows.size() - genomes_restored);
                    for (size_t i = 0; i < pending.size(); ++i) {
                        if (processed_rows[i]) continue;
                        remaining_pending.push_back(std::move(pending[i]));
                        remaining_rows.push_back(pending_input_rows[i]);
                    }
                    pending = std::move(remaining_pending);
                    pending_input_rows = std::move(remaining_rows);
                    total_records = pending.size();
                    resuming = true;

                    spdlog::info("Resuming from checkpoint: {}/{} genomes, {} shards, offset={}",
                                 genomes_restored, original_total_records, current_shard_id, ck_byte_offset);
                } else {
                    spdlog::warn("Checkpoint metadata corrupted (read {}/{}), starting fresh",
                                 genomes_restored, ck_genome_count);
                    catalog_rows.clear();
                    gidx_infos.clear();
                    accession_pairs.clear();
                    taxonomy_map.clear();
                    kmer_pairs.clear();
                    source_row_indices.clear();
                    shard_id_to_section_id.clear();
                    toc = TocWriter{};
                    next_section_id  = 1;
                    current_shard_id = 0;
                    next_genome_id   = cfg.starting_genome_id;
                    total_records    = pending.size();
                }
            }
        }

        // ── meta.tsv sidecar ─────────────────────────────────────────────────
        std::filesystem::path meta_tsv_path =
            gpk_path_.parent_path() / (gpk_path_.stem().string() + ".meta.tsv");

        if (!resuming) {
            app_writer.create(partial_path);

            // Write FileHeader (256B, GPK3)
            {
                FileHeader fhdr{};
                fhdr.magic         = GPK_MAGIC;
                fhdr.version_major = FORMAT_MAJOR;
                fhdr.version_minor = FORMAT_MINOR;
                // Real random 128-bit file UUID (P23).
                std::random_device rd;
                fhdr.file_uuid_lo  = (static_cast<uint64_t>(rd()) << 32) | rd();
                fhdr.file_uuid_hi  = (static_cast<uint64_t>(rd()) << 32) | rd();
                fhdr.created_at_unix = static_cast<uint64_t>(std::time(nullptr));
                fhdr.flags         = 0;
                fhdr.endian_abi_tag   = ENDIAN_ABI_TAG;                     // P26
                fhdr.build_params_hash = make_bprm_header_from_cfg(cfg).params_hash;
                fhdr.generation    = 1;
                app_writer.append(&fhdr, sizeof(fhdr));
            }

            std::ofstream meta_out(meta_tsv_path);
            meta_out << "accession\tgenome_id";
            if (!pending.empty())
                for (const auto& [k, v] : pending[0].extra_fields)
                    meta_out << "\t" << k;
            meta_out << "\n";
            // meta_out is a fresh stream; it will be reopened in append mode below
        }

        // Open meta_out in append mode (works for both fresh and resume)
        std::ofstream meta_out(meta_tsv_path, std::ios::app);

        uint32_t date = days_since_epoch();

        // ── IO writer thread: drains frozen shards and writes them to disk ─────
        struct WriteTask {
            std::future<FrozenShard> fut;
            uint64_t section_id;
            ShardId  shard_id;
            uint32_t n_genomes;
            size_t   catalog_start;
        };
        const size_t            write_q_max = 2;
        std::queue<WriteTask>   write_q;
        std::mutex              write_q_mx;
        std::condition_variable write_q_cv;
        std::condition_variable write_q_space_cv;
        bool                    writer_done = false;

        // ── Checkpoint write lambda ───────────────────────────────────────────
        FILE* ckpt_meta_file = std::fopen(ckpt_meta_path.c_str(),
                                          resuming ? "ab" : "wb");
        if (!ckpt_meta_file)
            throw std::runtime_error("Cannot open checkpoint meta file: " + ckpt_meta_path.string());

        struct DrainResult {
            ShardId  shard_id;
            uint64_t section_id;
            uint64_t file_offset;
            uint64_t compressed_size;
            uint32_t n_genomes;
            size_t   catalog_start;
        };

        auto write_checkpoint = [&](const DrainResult& dr) {
            // Append shard block to meta bin
            ShardCkptHdr shdr{};
            shdr.shard_id        = dr.shard_id;
            shdr.n_genomes       = dr.n_genomes;
            shdr.section_id      = dr.section_id;
            shdr.file_offset     = dr.file_offset;
            shdr.compressed_size = dr.compressed_size;
            std::fwrite(&shdr, sizeof(shdr), 1, ckpt_meta_file);

            for (size_t i = dr.catalog_start; i < dr.catalog_start + dr.n_genomes; ++i) {
                const auto& cm  = catalog_rows[i];
                const auto& acc = accession_pairs[i].first;
                std::string_view tax;
                auto tit = taxonomy_map.find(acc);
                if (tit != taxonomy_map.end()) tax = tit->second;

                const auto& kp = kmer_pairs[i];

                GenomeCkptFixed gf{};
                gf.genome_id          = cm.genome_id;
                gf.input_row_index    = source_row_indices[i];
                gf.shard_id           = cm.shard_id;
                gf.dir_index          = gidx_infos[i].dir_index;
                gf.oph_fingerprint    = cm.oph_fingerprint;
                gf.genome_length      = cm.genome_length;
                gf.n_contigs          = cm.n_contigs;
                gf.gc_pct_x100        = cm.gc_pct_x100;
                gf.completeness_x10   = cm.completeness_x10;
                gf.contamination_x10  = cm.contamination_x10;
                gf.date_added         = cm.date_added;
                gf.genome_type        = cm.genome_type;
                gf.accession_len      = static_cast<uint16_t>(acc.size());
                gf.taxonomy_len       = static_cast<uint16_t>(tax.size());
                std::memcpy(gf.kmer4, kp.second.data(), sizeof(gf.kmer4));

                std::fwrite(&gf, sizeof(gf), 1, ckpt_meta_file);
                if (gf.accession_len) std::fwrite(acc.data(),  1, gf.accession_len, ckpt_meta_file);
                if (gf.taxonomy_len)  std::fwrite(tax.data(),  1, gf.taxonomy_len,  ckpt_meta_file);
            }
            std::fflush(ckpt_meta_file);

            // Update scalar checkpoint atomically (tmp + rename). Non-fatal:
            // a transient NFS error here just means we lose the resume point
            // for this shard but the archive data is already durable.
            try {
                std::filesystem::path ckpt_tmp =
                    std::filesystem::path(ckpt_path.string() + ".tmp");
                {
                    std::ofstream ck(ckpt_tmp, std::ios::trunc);
                    ck << "genome_count="     << (dr.catalog_start + dr.n_genomes) << "\n"
                       << "byte_offset="      << app_writer.current_offset()       << "\n"
                       << "next_genome_id="   << next_genome_id                    << "\n"
                       << "current_shard_id=" << current_shard_id                  << "\n"
                       << "next_section_id="  << next_section_id                   << "\n";
                }
                std::filesystem::rename(ckpt_tmp, ckpt_path);
            } catch (const std::exception& e) {
                spdlog::warn("checkpoint write skipped (transient error): {}", e.what());
            }
        };

        // ── IO writer thread ──────────────────────────────────────────────────
        std::thread io_writer([&]() {
            while (true) {
                WriteTask wt;
                {
                    std::unique_lock lk(write_q_mx);
                    write_q_cv.wait(lk, [&]{ return !write_q.empty() || writer_done; });
                    if (write_q.empty()) break;
                    wt = std::move(write_q.front());
                    write_q.pop();
                }
                write_q_space_cv.notify_one();

                FrozenShard frozen = wt.fut.get();
                uint64_t shard_start = app_writer.current_offset();
                uint64_t shard_bytes;
                if (!frozen.tmp_path.empty()) {
                    FrozenShardReader rdr(frozen); // RAII: closes + unlinks on scope exit
                    app_writer.append_from_fd(rdr.fd, frozen.section_size);
                    shard_bytes = frozen.section_size;
                } else {
                    app_writer.append(frozen.bytes.data(), frozen.bytes.size());
                    shard_bytes = static_cast<uint64_t>(frozen.bytes.size());
                }
                // Make the shard bytes server-durable BEFORE recording the resume
                // offset in the checkpoint, so a crash cannot leave a checkpoint
                // pointing past bytes the NFS server never committed (P13).
                app_writer.flush();
                SectionDesc sd{};
                sd.type             = SEC_SHRD;
                sd.version          = 4;
                sd.flags            = 0;
                sd.section_id       = wt.section_id;
                sd.file_offset      = shard_start;
                sd.compressed_size  = shard_bytes;
                sd.uncompressed_size = 0;
                sd.item_count       = wt.n_genomes;
                sd.aux0             = wt.shard_id;
                sd.aux1             = 0;
                // XXH128 content hash of the whole shard section body (P1) —
                // checked by `genopack verify`.
                if (!frozen.tmp_path.empty()) {
                    if (!checksum_of_fd(app_writer.fd(), shard_start, shard_bytes, sd.checksum))
                        throw std::runtime_error("builder drain: checksum_of_fd failed");
                } else {
                    checksum_of(frozen.bytes.data(), frozen.bytes.size(), sd.checksum);
                }
                toc.add_section(sd);
                DrainResult dr{wt.shard_id, wt.section_id, shard_start, shard_bytes,
                               wt.n_genomes, wt.catalog_start};
                write_checkpoint(dr);
            }
        });

        size_t catalog_start_for_current_shard = catalog_rows.size();

        std::unique_ptr<ShardWriter> shard_writer;
        uint32_t genomes_in_current_shard = 0;

        auto open_shard = [&]() {
            catalog_start_for_current_shard = catalog_rows.size();
            shard_writer = std::make_unique<ShardWriter>(
                current_shard_id, current_shard_id, cfg.shard_cfg);
            ++current_shard_id;
            genomes_in_current_shard = 0;
        };

        auto launch_shard_freeze = [&]() {
            if (!shard_writer || shard_writer->n_genomes() == 0) return;
            ShardId flushed_id = current_shard_id - 1;
            uint64_t sid       = next_section_id++;
            uint32_t ng        = static_cast<uint32_t>(shard_writer->n_genomes());
            shard_id_to_section_id[flushed_id] = sid; // pre-assign for GIDX
            WriteTask wt;
            wt.fut = std::async(std::launch::async,
                [sw = std::move(shard_writer)]() mutable { return sw->freeze(); });
            wt.section_id    = sid;
            wt.shard_id      = flushed_id;
            wt.n_genomes     = ng;
            wt.catalog_start = catalog_start_for_current_shard;
            {
                std::unique_lock lk(write_q_mx);
                write_q_space_cv.wait(lk, [&]{ return write_q.size() < write_q_max; });
                write_q.push(std::move(wt));
            }
            write_q_cv.notify_one();
            shard_writer.reset();
        };

        // ── Worker pool: io_threads persistent threads ─────────────────────────
        const size_t n_workers = std::max(size_t(1), cfg.io_threads);
        const size_t sort_buf  = n_workers * 4;

        struct ChunkItem {
            BuildRecord record;
            GenomeId    genome_id;
            uint64_t    input_row_index;
            FastaStats  stats;
            std::string fasta;
            OPHDualSketchResult sketch;                    // populated when build_sketch (single-k)
            std::vector<OPHDualSketchResult> sketches_mk;  // populated when multi_k_sketch
            std::vector<ContigProfile> long_profiles;      // long-contig TNF profiles for GCOV (worker-computed)
            std::vector<uint64_t> aamers;                  // sorted-unique k=8 protein aamers (worker-computed)
            std::vector<uint64_t> fmh;                     // FracMinHash hashes, ≥1kb contigs (worker-computed)
            float contam_dup_excess = NAN;                 // SCC core_dup excess (worker-computed; fasta is freed before finalize)
            float contam_dup_mass   = NAN;                 // Phase-1: non-saturating SCC dup mass (worker-computed; not serialized in Phase-1)
        };

        // Task queue (all tasks submitted upfront; poison-pill sentinel at end)
        // fd=-1 on poison pill; producer opens file + fadvise(WILLNEED) before queuing
        // so the kernel starts NFS prefetch 4*n_workers genomes ahead of workers.
        struct Task { BuildRecord* record; GenomeId gid; uint64_t input_row_index; int fd = -1;
                      std::shared_ptr<std::string> seq; };  // from-gpk: in-memory decoded FASTA
        std::queue<Task>        task_q;
        std::mutex              task_mx;
        std::condition_variable task_cv;

        // Completion queue — sized to absorb n_workers genomes during a full shard
        // flush (~500 ms) so workers never stall on done_push_cv.
        const size_t            done_q_max = n_workers * 4;
        struct Done { std::optional<ChunkItem> item; };
        std::queue<Done>        done_q;
        std::mutex              done_mx;
        std::condition_variable done_cv;
        std::condition_variable done_push_cv;

        const size_t total = total_records;

        // Sort pending by genus so same-genus genomes arrive consecutively at
        // the drain loop. Without this, random input order causes the 16 GB
        // global cap to fire while genus buckets are still sparse, fragmenting
        // each genus across many tiny shards and making sketch I/O expensive.
        if (cfg.taxonomy_group) {
            const char rank = cfg.taxonomy_rank.empty() ? 'g' : cfg.taxonomy_rank[0];
            auto genus_key = [&](const BuildRecord& r) {
                for (const auto& [k, v] : r.extra_fields)
                    if (k == "taxonomy") return extract_taxonomy_bucket(v, rank);
                return std::string("__unclassified__");
            };
            std::vector<size_t> idx(pending.size());
            std::iota(idx.begin(), idx.end(), 0);
            std::stable_sort(idx.begin(), idx.end(),
                [&](size_t a, size_t b){ return genus_key(pending[a]) < genus_key(pending[b]); });
            std::vector<BuildRecord> sp;
            std::vector<uint64_t>    sr;
            sp.reserve(pending.size());
            sr.reserve(pending_input_rows.size());
            for (size_t i : idx) {
                sp.push_back(std::move(pending[i]));
                sr.push_back(pending_input_rows[i]);
            }
            pending = std::move(sp);
            pending_input_rows = std::move(sr);
        }

        // Streaming producer: fadvise(WILLNEED) is issued before the lock, so the queue
        // depth determines how far ahead NFS prefetch runs. n_workers*16 = 256 ahead
        // gives ~4s of prefetch runway at 250ms/genome × 16 workers, covering NFS cold-
        // miss latency of 100-400ms/file for random-order inputs.
        const size_t task_q_max = n_workers * 16;

        std::thread producer([&]() {
            if (src_reader_) {
                // from-gpk: stream decoded sequence shard-by-shard from the source
                // archive (sequential reads, no per-genome NFS file opens). Each
                // genome's FASTA is handed to a worker via the in-memory Task.seq;
                // the drain loop re-buckets by genus regardless of arrival order.
                std::vector<std::string> accs;
                accs.reserve(total_records);
                for (size_t i = 0; i < total_records; ++i) accs.push_back(pending[i].accession);
                src_reader_->visit_shard_batches(accs,
                    [&](ArchiveReader::ShardBatch& batch) {
                        for (auto& pr : batch) {
                            size_t idx = pr.first;
                            auto seq = std::make_shared<std::string>(std::move(pr.second.fasta));
                            std::unique_lock lk(task_mx);
                            task_cv.wait(lk, [&]{ return task_q.size() < task_q_max; });
                            Task t;
                            t.record          = &pending[idx];
                            t.gid             = next_genome_id++;
                            t.input_row_index = pending_input_rows[idx];
                            t.fd              = -1;
                            t.seq             = std::move(seq);
                            task_q.push(std::move(t));
                            lk.unlock();
                            task_cv.notify_one();
                        }
                    });
            } else {
                for (size_t i = 0; i < total_records; ++i) {
                    // Open + fadvise before acquiring the lock so kernel starts NFS
                    // prefetch while we wait for queue space. Workers read from this fd.
                    int fd = ::open(pending[i].file_path.c_str(), O_RDONLY);
#ifdef POSIX_FADV_WILLNEED
                    if (fd >= 0)
                        ::posix_fadvise(fd, 0, 0, POSIX_FADV_WILLNEED | POSIX_FADV_SEQUENTIAL);
#endif
                    std::unique_lock lk(task_mx);
                    task_cv.wait(lk, [&]{ return task_q.size() < task_q_max; });
                    task_q.push({&pending[i], next_genome_id++, pending_input_rows[i], fd});
                    lk.unlock();
                    task_cv.notify_one();
                }
            }
            for (size_t i = 0; i < n_workers; ++i) {
                std::unique_lock lk(task_mx);
                task_cv.wait(lk, [&]{ return task_q.size() < task_q_max; });
                task_q.push({nullptr, 0, 0, -1});
                lk.unlock();
                task_cv.notify_one();
            }
        });

        // Contamination (SCC core_dup) at build time — own 6-frame pass per genome via the
        // .csp panel. Declared BEFORE the worker pool and scored INSIDE the worker (while the
        // FASTA is alive): finalize_genus runs after flush_staging_buf frees item.fasta, so the
        // result is carried on ChunkItem.contam_dup_excess. panel.score is const/thread-safe.
        std::unique_ptr<CladeSplitPanel> build_contam_panel;
        if (!cfg.contam_panel_path.empty()) {
            build_contam_panel = std::make_unique<CladeSplitPanel>();
            try {
                build_contam_panel->open(cfg.contam_panel_path);
                spdlog::info("build: contamination panel {} aamers, {} clades",
                             build_contam_panel->n_aamers(), build_contam_panel->n_clades());
            } catch (const std::exception& ex) {
                spdlog::warn("build: cannot open contamination panel ({}); skipping", ex.what());
                build_contam_panel.reset();
            }
        }

        std::vector<std::thread> workers;
        workers.reserve(n_workers);
        for (size_t i = 0; i < n_workers; ++i) {
            workers.emplace_back([&]() {
                while (true) {
                    Task t;
                    {
                        std::unique_lock lk(task_mx);
                        task_cv.wait(lk, [&]{ return !task_q.empty(); });
                        t = task_q.front();
                        task_q.pop();
                    }
                    task_cv.notify_one();
                    if (!t.record) return; // poison pill

                    Done d;
                    try {
                        std::string fasta = t.seq
                            ? std::move(*t.seq)
                            : ((t.fd >= 0)
                                ? decompress_gz_fd(t.fd, t.record->file_path)
                                : decompress_gz(t.record->file_path));
                        FastaStats stats  = compute_fasta_stats(fasta);
                        OPHDualSketchResult sk;
                        std::vector<OPHDualSketchResult> sks_mk;
                        if (cfg.build_sketch) {
                            const uint64_t seed1 = cfg.sketch_seed;
                            const uint64_t seed2 = cfg.sketch_seed + 1;
                            if (multi_k_sketch) {
                                sks_mk = sketch_oph_dual_multik(
                                    fasta.data(), fasta.size(),
                                    std::span<const int>(cfg.sketch_kmer_sizes),
                                    cfg.sketch_size, cfg.sketch_syncmer_s,
                                    seed1, seed2);
                            } else {
                                sk = sketch_oph_dual_from_buffer(
                                    fasta.data(), fasta.size(),
                                    cfg.sketch_kmer_size, cfg.sketch_size,
                                    cfg.sketch_syncmer_s,
                                    seed1, seed2);
                            }
                        }
                        // Long-contig TNF profiles are consensus-INDEPENDENT — compute them HERE
                        // in the parallel worker (this was the single-threaded drain bottleneck).
                        std::vector<ContigProfile> lprof =
                            compute_long_contig_profiles(fasta, GCOV_MIN_LONG_BP);
                        // One per-contig dewrap walk → aamers (sort+unique) + FMH hashes. Both
                        // consensus-INDEPENDENT, so computed HERE in the parallel worker (these
                        // were the last per-genome FASTA walks in the serial drain). The dewrap
                        // must match the prior in-drain walks exactly for model-hash parity.
                        std::vector<uint64_t> aamers, fmh;
                        {
                            const bool want_aamers = cfg.build_pcore;
                            uint64_t fmax = UINT64_MAX;
                            if (want_aamers) {
                                const char* af = std::getenv("GENOPACK_AAMER_FRAC");
                                const uint64_t div = (af && *af) ? std::strtoull(af, nullptr, 10)
                                    : static_cast<uint64_t>(cfg.pcore_frac > 0 ? cfg.pcore_frac : 1u);
                                if (div > 1) fmax = UINT64_MAX / div;
                            }
                            const int fmh_c_vals[1] = { cfg.fmh_c };
                            const char* wp = fasta.data(), *we = wp + fasta.size();
                            std::string wseq;
                            while (wp < we) {
                                while (wp < we && *wp != '>') ++wp;
                                if (wp >= we) break;
                                while (wp < we && *wp != '\n') ++wp;
                                if (wp < we) ++wp;
                                const char* ws = wp; size_t wl = 0;
                                while (wp < we && *wp != '>') {
                                    while (wp < we && *wp != '\n' && *wp != '\r') { ++wl; ++wp; }
                                    while (wp < we && (*wp == '\n' || *wp == '\r')) ++wp;
                                }
                                if (wl == 0) continue;
                                wseq.clear(); wseq.reserve(wl);
                                for (const char* sp = ws; sp < wp; ++sp)
                                    if (*sp != '\n' && *sp != '\r') wseq.push_back(*sp);
                                if (want_aamers) extract_aamers_dna_into(wseq, AAMER_K, 8, fmax, aamers);
                                if (wl >= 1000) {   // kFmhMinBp: FMH on ≥1kb contigs
                                    auto vecs = fmh_multi_c(wseq, cfg.fmh_k, fmh_c_vals, 1);
                                    if (!vecs.empty()) fmh.insert(fmh.end(), vecs[0].begin(), vecs[0].end());
                                }
                            }
                            if (want_aamers) {
                                thread_local std::vector<uint64_t> radix_scratch;
                                radix_sort_dedup_u64(aamers, radix_scratch);
                            }
                        }
                        d.item = ChunkItem{*t.record, t.gid, t.input_row_index,
                                           stats, std::move(fasta), std::move(sk), std::move(sks_mk),
                                           std::move(lprof),
                                           std::move(aamers), std::move(fmh)};
                        // SCC core_dup excess — score here (FASTA alive); finalize_genus reads
                        // the stored value (item.fasta is freed by flush_staging_buf before it runs).
                        if (build_contam_panel) {
                            const auto cds = build_contam_panel->score(d.item->fasta);
                            d.item->contam_dup_excess = cds.core_dup_excess;
                            d.item->contam_dup_mass   = cds.core_dup_mass;
                        }
                    } catch (const std::exception& ex) {
                        spdlog::warn("Skipping {}: {}", t.record->accession, ex.what());
                    }
                    {
                        std::unique_lock lk(done_mx);
                        done_push_cv.wait(lk, [&]{ return done_q.size() < done_q_max; });
                        done_q.push(std::move(d));
                    }
                    done_cv.notify_one();
                }
            });
        }
        task_cv.notify_all();

        // Staging-buffer flush: sort a full shard-sized batch before freezing.
        size_t n_done = 0, n_failed = 0;

        // When taxonomy_group is enabled: per-taxon buckets (flushed when bucket bytes >= threshold).
        // When disabled: single flat buffer (same as before, but sorted by kmer NN chain or oph).
        struct TaxonBucket {
            std::vector<ChunkItem> items;
            uint64_t raw_bytes = 0;        // FASTA bytes (drives 512MB shard sizing)
            uint64_t aamer_bytes = 0;      // worker aamer arrays — ~10x the FASTA at keep-all,
                                           // so they, not the FASTA, dominate buffered RAM
            // Cross-chunk CORE/PCORE aggregation: a mega-genus splits across many
            // 512MB shard-chunks; fold each chunk's aamer counts here and emit one
            // full-density CORE/PCORE once all the genus's members have arrived.
            std::unique_ptr<GenusAamerCounter> genus_counter;   // lazily allocated on first fold
            std::vector<uint64_t> member_ids;
            uint32_t arrived = 0;
            bool model_emitted = false;
        };
        std::unordered_map<std::string, TaxonBucket> taxon_buckets;
        uint64_t taxon_total_bytes = 0; // sum of bucket FASTA bytes
        uint64_t taxon_total_aamer = 0; // sum of bucket aamer-array bytes (the RAM driver)
        // Global cap: flush the largest bucket whenever total in-memory genome
        // data exceeds this limit. Prevents unbounded accumulation across many
        // small-genus buckets that never individually hit max_shard_size_bytes.
        const uint64_t taxon_global_cap = 2ULL << 30;   // 2 GB FASTA buffered
        const uint64_t taxon_aamer_cap  = 512ULL << 20; // 512 MB aamer arrays buffered (FMH-subsampled)
        // Cross-genus bin-pack buffer for micro-genera flushed immediately at completion.
        std::vector<ChunkItem> micro_combine;
        uint64_t micro_combine_bytes = 0;
        std::vector<ChunkItem> staging_buffer; // used when !cfg.taxonomy_group
        staging_buffer.reserve(sort_buf * 4);
        uint64_t staging_raw_bytes = 0;

        // In-memory GCOV/FCOV/FMHR writers populated during flush_staging_buf (one-pass).
        GcovWriter build_gcov_w, build_fcov_w;
        // FmhrWriter spills per-genus FMH hash sets to disk to avoid accumulating
        // O(n_genomes × hashes_per_genome) in RAM throughout the build.
        FmhrWriter build_fmhr_w(cfg.tmpdir.empty()
            ? std::filesystem::temp_directory_path().string()
            : cfg.tmpdir);

        // ── Genus-finalization futures (Steps 4-5) ─────────────────────────────
        // Each genus's GSTX/GCOV/QUAL block is dispatched to the thread pool so
        // all n_workers threads stay busy while the main thread drains the queue.
        std::vector<std::future<void>> genus_futs;
        genus_futs.reserve(256);
        auto reap_futs = [&]() {
            genus_futs.erase(
                std::remove_if(genus_futs.begin(), genus_futs.end(),
                    [](std::future<void>& f) {
                        return f.wait_for(std::chrono::seconds(0))
                               == std::future_status::ready;
                    }),
                genus_futs.end());
        };
        // Mutexes protecting each writer from concurrent finalize_genus calls.
        // core_writer_mx guards the (shared) PCORE writer during finalize.
        std::mutex gstx_writer_mx, gcov_writer_mx, fmhr_writer_mx,
                   core_writer_mx, qual_writer_mx;
        // SEC_PCORE: unified DENSE per-genus aamer reference (every aamer + member
        // count). The prev>=ceil(theta*n) slice reproduces CORE; the full set is the
        // dense reference the small-contig foreign-contamination channel needs. Built
        // inline from the same chunk_aamers (zero extra passes); genus tier only here
        // (the family tier is added post-hoc by `genopack pcore`, like FCORE).
        PcoreWriter build_pcore_w(AAMER_K, /*min_seg_aa=*/8, cfg.core_theta, cfg.tmpdir);
        PcoreSpillWriter pcore_spiller(cfg.tmpdir.empty()
            ? std::filesystem::temp_directory_path()
            : std::filesystem::path(cfg.tmpdir),
            /*buf_tuples=*/4096);
        if (cfg.build_tier && cfg.build_pcore) {
            const std::string tier_path = gpk_path_.string() + ".ptier";
            build_pcore_w.set_tier_sidechannel(tier_path);
            spdlog::info("build: tier side-channel → {}", tier_path);
        }
        static constexpr float  kGcovOutlierPct = 0.99f;
        const int               kFmhK           = cfg.fmh_k, kFmhC = cfg.fmh_c;
        static constexpr uint32_t kFmhMinBp     = 1000;

        // Marker scoring at build time — same FASTA pass as GCOV/FMH.
        std::unique_ptr<MarkerReader> build_mrk_rd;
        // build_pcore_frac_max: FMH threshold for PCORE/CORE aamer extraction.
        // NOT overridden by the marker panel — PCORE and markers may use different rates.
        // Priority: GENOPACK_AAMER_FRAC env > cfg.pcore_frac (default 100) > 1.
        uint64_t build_pcore_frac_max = UINT64_MAX;
        {
            const char* af = std::getenv("GENOPACK_AAMER_FRAC");
            const uint64_t div = (af && *af) ? std::strtoull(af, nullptr, 10)
                                             : static_cast<uint64_t>(cfg.pcore_frac > 0 ? cfg.pcore_frac : 1u);
            if (div > 1) { build_pcore_frac_max = UINT64_MAX / div;
                spdlog::info("build: aamer subsample 1/{} (CORE/PCORE)", div); }
        }
        // build_mrk_frac_max: FMH threshold for marker scoring only.
        // Overridden by the marker panel's own frac_max when --markers is used.
        uint64_t build_mrk_frac_max = build_pcore_frac_max;
        if (!cfg.markers_path.empty()) {
            build_mrk_rd = std::make_unique<MarkerReader>();
            try {
                build_mrk_rd->open(cfg.markers_path);
                build_mrk_rd->build_merged_pool();   // copies pools to owned vectors
                build_mrk_frac_max = build_mrk_rd->frac_max_hash();
                spdlog::info("build: marker panel {} bac + {} arc, {} lineages",
                             build_mrk_rd->n_bac(), build_mrk_rd->n_arc(),
                             build_mrk_rd->n_lineages());
            } catch (const std::exception& ex) {
                spdlog::warn("build: cannot open markers ({}); skipping", ex.what());
                build_mrk_rd.reset();
            }
        }


        // ── finalize_genus ──────────────────────────────────────────────────────
        // Computes GSTX consensus, GCOV covariance, QUAL records, CORE/PCORE for
        // one homogeneous genus chunk. Designed to run off the main thread: all
        // reads of `buf` are local, all writer calls are mutex-guarded.
        // `buf` is taken by value so callers can std::move it in (Step 5).
        auto finalize_genus = [&](std::vector<ChunkItem> buf,
                                  std::unique_ptr<GenusAamerCounter>* accum_holder,
                                  std::vector<uint64_t>* accum_member_ids) {
            if (buf.size() < 2) return;

            // genus_homogeneous check (same logic as flush_staging_buf)
            {
                std::string g0;
                for (const auto& [k, v] : buf[0].record.extra_fields)
                    if (k == "taxonomy") { g0 = extract_taxonomy_bucket(v, 'g'); break; }
                for (size_t gi = 1; gi < buf.size(); ++gi)
                    for (const auto& [k, v] : buf[gi].record.extra_fields)
                        if (k == "taxonomy") {
                            if (extract_taxonomy_bucket(v, 'g') != g0) return;
                            break;
                        }
            }

            if (!cfg.build_gstx || !cfg.taxonomy_group || !multi_k_sketch) return;

            std::string genus_key;
            std::string species_key;
            for (const auto& [k, v] : buf[0].record.extra_fields)
                if (k == "taxonomy") {
                    genus_key   = extract_taxonomy_bucket(v, 'g');
                    species_key = extract_taxonomy_bucket(v, 's');
                    break;
                }
            const std::string& gstx_key = (!species_key.empty() && species_key != "__unclassified__")
                                          ? species_key : genus_key;

            const int nk   = std::min((int)buf[0].sketches_mk.size(), (int)GSTX_MAX_K);
            const int bins = cfg.sketch_size;

            if (nk <= 0 || genus_key.empty() || genus_key == "__unclassified__"
                || bins != static_cast<int>(GSTX_BINS)) return;

            // Pass 1: Boyer-Moore majority vote → consensus signature per k
            std::vector<std::vector<uint16_t>> cand(nk, std::vector<uint16_t>(bins, 0));
            std::vector<std::vector<int32_t>>  vote(nk, std::vector<int32_t>(bins, 0));
            for (const auto& item : buf) {
                for (int ki = 0; ki < nk && ki < (int)item.sketches_mk.size(); ++ki) {
                    const auto& sk = item.sketches_mk[ki];
                    for (int b = 0; b < bins; ++b) {
                        uint16_t v = sk.signature1[b];
                        if      (vote[ki][b] == 0) { cand[ki][b] = v; vote[ki][b] = 1; }
                        else if (v == cand[ki][b]) { ++vote[ki][b]; }
                        else                       { --vote[ki][b]; }
                    }
                }
            }

            // Pass 2: containment vs consensus → exact p90 (all data in RAM)
            float p90[GSTX_MAX_K] = {};
            // Robust dispersion of the same member-containment vector p90 comes from:
            // median and 1.4826·MAD per k. NaN sentinel = not computed (empty genus).
            float med_c[GSTX_MAX_K] = {NAN, NAN, NAN};
            float mad_c[GSTX_MAX_K] = {NAN, NAN, NAN};
            std::vector<float> cont0_v;
            cont0_v.reserve(buf.size());
            if (nk > 0) {
                for (const auto& item : buf) {
                    if (item.sketches_mk.empty()) { cont0_v.push_back(0.0f); continue; }
                    const auto& sk = item.sketches_mk[0];
                    int match = 0;
                    for (int b = 0; b < bins; ++b)
                        if (sk.signature1[b] == cand[0][b]) ++match;
                    cont0_v.push_back(static_cast<float>(match) / bins);
                }
            }
            const float max_c0 = cont0_v.empty() ? 0.0f
                : *std::max_element(cont0_v.begin(), cont0_v.end());
            const float keep_thr = GSTX_P90_COMPLETE_FRAC * max_c0;

            for (int ki = 0; ki < nk; ++ki) {
                std::vector<float> c;
                c.reserve(buf.size());
                for (size_t ii = 0; ii < buf.size(); ++ii) {
                    if (ki >= (int)buf[ii].sketches_mk.size()) continue;
                    if (!cont0_v.empty() && cont0_v[ii] < keep_thr) continue;
                    const auto& sk = buf[ii].sketches_mk[ki];
                    int match = 0;
                    for (int b = 0; b < bins; ++b)
                        if (sk.signature1[b] == cand[ki][b]) ++match;
                    c.push_back(static_cast<float>(match) / bins);
                }
                if (c.empty()) {
                    for (const auto& item : buf) {
                        if (ki >= (int)item.sketches_mk.size()) continue;
                        const auto& sk = item.sketches_mk[ki];
                        int match = 0;
                        for (int b = 0; b < bins; ++b)
                            if (sk.signature1[b] == cand[ki][b]) ++match;
                        c.push_back(static_cast<float>(match) / bins);
                    }
                }
                if (!c.empty()) {
                    // median + 1.4826·MAD (compute before nth_element reorders for p90).
                    std::vector<float> t = c;
                    const size_t mid = t.size() / 2;
                    std::nth_element(t.begin(), t.begin() + mid, t.end());
                    const float med = t[mid];
                    med_c[ki] = med;
                    for (float& v : t) v = std::fabs(v - med);
                    std::nth_element(t.begin(), t.begin() + mid, t.end());
                    mad_c[ki] = 1.4826f * t[mid];

                    size_t idx = static_cast<size_t>(0.9f * c.size());
                    if (idx >= c.size()) idx = c.size() - 1;
                    std::nth_element(c.begin(), c.begin() + idx, c.end());
                    p90[ki] = c[idx];
                }
            }

            // TNF centroid
            float tnf_mu[136] = {};
            for (const auto& item : buf)
                for (int d = 0; d < 136; ++d) tnf_mu[d] += item.stats.kmer4_profile[d];
            const bool has_tnf = buf.size() >= 2;
            if (has_tnf)
                for (int d = 0; d < 136; ++d) tnf_mu[d] /= static_cast<float>(buf.size());

            uint32_t ksizes[GSTX_MAX_K] = {};
            for (int ki = 0; ki < nk; ++ki)
                ksizes[ki] = static_cast<uint32_t>(cfg.sketch_kmer_sizes[ki]);

            // nrb_p90: p90 of n_real_bins at k=0 — used for sketch-fill completeness
            float nrb_p90_gstx = 0.0f;
            if (nk > 0) {
                std::vector<float> nrb_v;
                nrb_v.reserve(buf.size());
                for (const auto& item : buf)
                    if (!item.sketches_mk.empty())
                        nrb_v.push_back(static_cast<float>(item.sketches_mk[0].n_real_bins));
                if (!nrb_v.empty()) {
                    size_t idx = static_cast<size_t>(0.9f * static_cast<float>(nrb_v.size()));
                    if (idx >= nrb_v.size()) idx = nrb_v.size() - 1;
                    std::nth_element(nrb_v.begin(), nrb_v.begin() + static_cast<ptrdiff_t>(idx), nrb_v.end());
                    nrb_p90_gstx = nrb_v[idx];
                }
            }

            {
                std::lock_guard<std::mutex> lk(gstx_writer_mx);
                gstx_writer.add_genus(gstx_key,
                                      static_cast<uint32_t>(buf.size()),
                                      static_cast<uint32_t>(nk),
                                      cand, p90,
                                      has_tnf ? tnf_mu : nullptr,
                                      ksizes,
                                      nrb_p90_gstx,
                                      med_c, mad_c);
            }

            // Per-genome quality records + long-contig profiles for GCOV scoring.
            const float k0 = nk >= 1 ? static_cast<float>(cfg.sketch_kmer_sizes[0]) : 0.0f;
            const float k1 = nk >= 2 ? static_cast<float>(cfg.sketch_kmer_sizes[1]) : 0.0f;
            const float k2 = nk >= 3 ? static_cast<float>(cfg.sketch_kmer_sizes[2]) : 0.0f;

            std::vector<QualRecord>                 pending_qrs;
            std::vector<std::vector<ContigProfile>> all_profiles;
            pending_qrs.reserve(buf.size());
            all_profiles.reserve(buf.size());

            for (const auto& item : buf) {
                auto qr               = QualRecord::make_empty(item.genome_id);
                qr.support_tier       = 0; // GenusSaturated

                float cont[GSTX_MAX_K] = {};
                for (int ki = 0; ki < nk && ki < (int)item.sketches_mk.size(); ++ki) {
                    const auto& sk = item.sketches_mk[ki];
                    int match = 0;
                    for (int b = 0; b < bins; ++b)
                        if (sk.signature1[b] == cand[ki][b])
                            ++match;
                    cont[ki] = static_cast<float>(match) / bins;
                }

                if (p90[0] > 0.01f)
                    qr.completeness_cluster_relative =
                        std::clamp(cont[0] / p90[0], 0.0f, 1.5f);

                if (nk >= 3 && cont[0] >= 0.01f && cont[1] >= 0.01f && cont[2] >= 0.01f) {
                    const float lc0 = std::log(cont[0]), lc1 = std::log(cont[1]), lc2 = std::log(cont[2]);
                    const float beta = (k0*lc0 + k1*lc1 + k2*lc2) / (k0*k0 + k1*k1 + k2*k2);
                    const float c2_pred = std::exp(beta * k2);
                    qr.contamination_leakage =
                        std::max(0.0f, c2_pred - cont[2]) / std::max(c2_pred, 0.01f);
                } else if (nk >= 2 && cont[0] >= 0.01f && cont[1] >= 0.01f) {
                    const float beta = (k0*std::log(cont[0]) + k1*std::log(cont[1])) / (k0*k0 + k1*k1);
                    const float c1_pred = std::exp(beta * k1);
                    qr.contamination_leakage =
                        std::max(0.0f, c1_pred - cont[1]) / std::max(c1_pred, 0.01f);
                }

                if (has_tnf) {
                    float d2 = 0.0f, norm2 = 0.0f;
                    for (int di = 0; di < 136; ++di) {
                        float diff = item.stats.kmer4_profile[di] - tnf_mu[di];
                        d2    += diff * diff;
                        norm2 += tnf_mu[di] * tnf_mu[di];
                    }
                    const float dist = std::sqrt(d2);
                    const float ref  = std::sqrt(norm2);
                    if (ref > 1e-6f)
                        qr.contamination_tnf_excess = std::max(0.0f, dist / ref - 1.0f);

                    qr.completeness_post_decontam =
                        compute_completeness_post_decontam(item.fasta, tnf_mu);
                }

                qr.qual_flags |= QualRecord::QUAL_FLAG_BUILD_SIGNALS;

                all_profiles.push_back(std::move(item.long_profiles));
                pending_qrs.push_back(std::move(qr));
            }

            // ── One-pass GCOV/FCOV/FMHR ──────────────────────────────────────
            if (cfg.build_gcov) {
                std::string family_key;
                for (const auto& [k, v] : buf[0].record.extra_fields)
                    if (k == "taxonomy") {
                        auto fp = v.find("f__");
                        if (fp != std::string::npos) {
                            auto fe = v.find(';', fp + 3);
                            family_key = v.substr(fp, fe == std::string::npos
                                                       ? v.size() - fp : fe - fp);
                            if (family_key == "f__") family_key.clear();
                        }
                        break;
                    }

                GenusAccum genus_acc, family_acc;
                std::vector<uint64_t> fmh_all;
                const bool do_markers = build_mrk_rd && build_mrk_rd->has_merged_pool();
                const bool need_aamers = do_markers || cfg.build_pcore;
                const bool aamers_capped = (cfg.genus_member_cap > 0) && accum_member_ids
                    && accum_member_ids->size() >= cfg.genus_member_cap;
                (void)need_aamers;
                (void)accum_holder;

                for (size_t ii = 0; ii < buf.size(); ++ii) {
                    genus_acc.add_genome(all_profiles[ii]);
                    if (!family_key.empty()) family_acc.add_genome(all_profiles[ii]);
                    fmh_all.insert(fmh_all.end(), buf[ii].fmh.begin(), buf[ii].fmh.end());
                }

                if (!aamers_capped && cfg.build_pcore) {
                    const uint64_t gkh_spill = GcovWriter::hash_genus(genus_key);
                    for (size_t ii = 0; ii < buf.size(); ++ii) {
                        if (buf[ii].aamers.empty()) continue;
                        pcore_spiller.add_member(gkh_spill);
                        pcore_spiller.add_genome(gkh_spill, buf[ii].aamers);
                    }
                }

                const uint32_t n_mem = static_cast<uint32_t>(buf.size());
                GcovEntry ge{}, fe{};
                bool fcov_ok = false;
                {
                    std::lock_guard<std::mutex> lk(gcov_writer_mx);
                    ge = finalize_and_add_genus(genus_key, n_mem, genus_acc, build_gcov_w);
                    if (!family_key.empty()) {
                        fe = finalize_and_add_genus(family_key, n_mem, family_acc, build_fcov_w);
                        fcov_ok = (fe.flags & GCOV_FLAG_VALID) != 0;
                    }
                }

                if (!fmh_all.empty()) {
                    std::sort(fmh_all.begin(), fmh_all.end());
                    fmh_all.erase(std::unique(fmh_all.begin(), fmh_all.end()), fmh_all.end());
                    std::lock_guard<std::mutex> lk(fmhr_writer_mx);
                    build_fmhr_w.add(GcovWriter::hash_genus(genus_key), std::move(fmh_all));
                }

                // Score contigs and populate outlier fields.
                const bool genus_ok = (ge.flags & GCOV_FLAG_VALID) != 0;
                for (size_t ii = 0; ii < pending_qrs.size(); ++ii) {
                    auto& qr        = pending_qrs[ii];
                    const auto& cps = all_profiles[ii];
                    qr.qual_flags  |= QualRecord::QUAL_FLAG_GCOV_SCORED;
                    if (!genus_ok || cps.empty()) continue;

                    uint32_t cco_bp=0, spe_bp=0, rho_bp=0, sib_bp=0, scored_bp=0;
                    for (const auto& cp : cps) {
                        scored_bp += cp.bp;
                        float xmu[136];
                        for (int d=0;d<136;++d) xmu[d] = cp.p[d] - ge.mu[d];
                        float mahal, spe_val;
                        gcov_mahalanobis_spe(ge, xmu, &mahal, &spe_val);
                        const float pct  = gcov_percentile(ge, mahal);
                        const float spct = gcov_spe_percentile(ge, spe_val);
                        const bool t2  = pct  >= kGcovOutlierPct;
                        const bool spe = spct >= kGcovOutlierPct;
                        if (t2 || spe) cco_bp += cp.bp;
                        if (spe)       spe_bp += cp.bp;
                        float rhod[16];
                        for (int i=0;i<16;++i) rhod[i] = cp.rho[i] - ge.rho_mean[i];
                        if (gcov_rho_percentile(ge, gcov_rho_distance(ge, rhod)) >= kGcovOutlierPct)
                            rho_bp += cp.bp;
                        if ((t2||spe) && fcov_ok) {
                            float xmuf[136];
                            for (int d=0;d<136;++d) xmuf[d] = cp.p[d] - fe.mu[d];
                            float fmahal, fspe;
                            gcov_mahalanobis_spe(fe, xmuf, &fmahal, &fspe);
                            const float fp  = gcov_percentile(fe, fmahal);
                            const float fsp = gcov_spe_percentile(fe, fspe);
                            if (!(fp >= kGcovOutlierPct || fsp >= kGcovOutlierPct))
                                sib_bp += cp.bp;
                        }
                    }
                    if (scored_bp > 0) {
                        auto enc = [&](uint32_t bp) {
                            return static_cast<uint8_t>(
                                std::min(255u, (bp * 255u) / scored_bp));
                        };
                        qr.contig_outlier_u8  = enc(cco_bp);
                        qr.spe_outlier_u8     = enc(spe_bp);
                        qr.rho_outlier_u8     = enc(rho_bp);
                    }
                }

                // Marker completeness — reuse aamers collected in FMH loop.
                if (do_markers) {
                    const uint64_t gh = GcovWriter::hash_genus(genus_key);
                    auto calib = build_mrk_rd->lookup_lineage(gh);
                    if (calib.valid()) {
                        const bool is_arc = (calib.header->domain == MRKR_DOMAIN_ARC);
                        auto mh  = is_arc ? build_mrk_rd->merged_hashes_arc()
                                          : build_mrk_rd->merged_hashes_bac();
                        auto mid = is_arc ? build_mrk_rd->merged_ids_arc()
                                          : build_mrk_rd->merged_ids_bac();
                        const uint8_t n_markers = calib.header->n_markers;
                        for (size_t ii = 0; ii < pending_qrs.size(); ++ii) {
                            const auto& qv = buf[ii].aamers;
                            if (qv.empty()) continue;
                            uint32_t hits[173] = {};
                            const uint64_t *qp=qv.data(), *qpe=qp+qv.size();
                            const uint64_t *mhp=mh.data(), *mhe=mhp+mh.size();
                            const uint8_t  *mip=mid.data();
                            while (qp!=qpe && mhp!=mhe) {
                                if      (*qp < *mhp) ++qp;
                                else if (*mhp < *qp) { ++mhp; ++mip; }
                                else                 { hits[*mip]++; ++qp; ++mhp; ++mip; }
                            }
                            int n_present=0, n_expected=0;
                            float redun_sum=0.0f; int redun_n=0;
                            for (uint8_t mi=0; mi<n_markers; ++mi) {
                                if (!calib.marker_expected(mi)) continue;
                                ++n_expected;
                                const uint32_t thr = static_cast<uint32_t>(
                                    calib.slots[mi].null_floor_u16)
                                    + static_cast<uint32_t>(cfg.marker_min_hits);
                                if (hits[mi] >= thr) {
                                    ++n_present;
                                    const uint8_t pool_mi = is_arc
                                        ? static_cast<uint8_t>(build_mrk_rd->n_bac()+mi)
                                        : mi;
                                    const uint32_t psz = build_mrk_rd->pool_n_hashes(pool_mi);
                                    if (psz > 0) {
                                        redun_sum += static_cast<float>(hits[mi])
                                                   / static_cast<float>(psz);
                                        ++redun_n;
                                    }
                                }
                            }
                            auto& qr = pending_qrs[ii];
                            if (n_expected > 0) {
                                const float comp = static_cast<float>(n_present)
                                                 / static_cast<float>(n_expected);
                                qr.marker_completeness_u8 = static_cast<uint8_t>(
                                    std::clamp(comp, 0.0f, 1.0f) * 254.0f + 1.0f);
                            }
                        }
                    }
                }

                // Contamination duplication — panel.score is const/thread-safe.
                // v3 panels: store the SCC-restricted core_dup EXCESS (relative to the placed
                // genus's clean ceiling) in the existing u16 slot. It flags same-genus/strain
                // duplication that redundancy_fraction (diagnostic-only) and CheckM/CheckM2 miss.
                // NAN (no SCC set: v2 panel or underpopulated/novel genus) -> encode_dup=0 sentinel
                // -> check reads NA and defers to the observability cap. Clean-but-scored -> 0.
                if (build_contam_panel) {
                    for (size_t ii = 0; ii < pending_qrs.size() && ii < buf.size(); ++ii) {
                        // contam_dup_excess was computed in the worker (FASTA alive). NAN when no
                        // SCC set (v2 panel / underpopulated genus) -> encode_dup(NAN)=0 sentinel
                        // -> check decodes NA -> observability cap.
                        pending_qrs[ii].contamination_duplication_u16 =
                            QualRecord::encode_dup(buf[ii].contam_dup_excess);
                        // Phase-2: non-saturating SCC dup mass stored raw (NAN = not scored;
                        // core_dup_mass is unbounded-ish so a u16 fraction would clip it).
                        pending_qrs[ii].contamination_core_dup_mass = buf[ii].contam_dup_mass;
                    }
                }
            }

            {
                std::lock_guard<std::mutex> lk(qual_writer_mx);
                for (auto& qr : pending_qrs) qual_writer.add(qr);
            }
        }; // end finalize_genus

        auto flush_staging_buf = [&](std::vector<ChunkItem>& buf,
                                     std::unique_ptr<GenusAamerCounter>* accum_holder = nullptr,
                                     std::vector<uint64_t>* accum_member_ids = nullptr) {
            if (buf.empty()) return;
            // Sort by kmer NN chain (if enabled) or oph_fingerprint
            if (cfg.kmer_nn_sort && buf.size() > 1) {
                auto order = greedy_nn_chain(buf);
                std::vector<ChunkItem> sorted;
                sorted.reserve(order.size());
                for (size_t idx : order) sorted.push_back(std::move(buf[idx]));
                buf = std::move(sorted);
            } else {
                std::sort(buf.begin(), buf.end(), [](const ChunkItem& a, const ChunkItem& b) {
                    return a.stats.oph_fingerprint < b.stats.oph_fingerprint;
                });
            }
            open_shard();
            for (auto& item : buf) {
                GenomeMeta meta{};
                meta.genome_id         = item.genome_id;
                meta.genome_type       = static_cast<uint32_t>(item.record.genome_type);
                meta.shard_id          = current_shard_id - 1;
                meta.genome_length     = item.record.genome_length > 0 ? item.record.genome_length : item.stats.genome_length;
                meta.n_contigs         = item.record.n_contigs > 0     ? item.record.n_contigs     : item.stats.n_contigs;
                meta.gc_pct_x100       = item.stats.gc_pct_x100;
                meta.completeness_x10  = static_cast<uint16_t>(item.record.completeness  * 10.0f);
                meta.contamination_x10 = static_cast<uint16_t>(item.record.contamination * 10.0f);
                meta.oph_fingerprint   = item.stats.oph_fingerprint;
                meta.date_added        = date;
                shard_writer->add_genome(item.genome_id, item.stats.oph_fingerprint,
                                         item.fasta.data(), item.fasta.size());
                gidx_infos.push_back({item.genome_id,
                                      static_cast<ShardId>(current_shard_id - 1),
                                      genomes_in_current_shard++});
                source_row_indices.push_back(item.input_row_index);
                catalog_rows.push_back(meta);
                accession_pairs.emplace_back(item.record.accession, item.genome_id);
                kmer_pairs.emplace_back(item.genome_id, item.stats.kmer4_profile);
                if (cfg.build_sketch) {
                    if (skch_writer_mk && !item.sketches_mk.empty()) {
                        std::vector<std::vector<uint16_t>> sigs1_per_k;
                        std::vector<std::vector<uint16_t>> sigs2_per_k;
                        std::vector<uint32_t>              n_real_bins_per_k;
                        std::vector<std::vector<uint64_t>> masks_per_k;
                        for (const auto& sk : item.sketches_mk) {
                            sigs1_per_k.push_back(sk.signature1);
                            sigs2_per_k.push_back(sk.signature2);
                            n_real_bins_per_k.push_back(sk.n_real_bins);
                            masks_per_k.push_back(sk.real_bins_bitmask);
                        }
                        skch_writer_mk->add(item.genome_id, item.stats.genome_length,
                                            sigs1_per_k, sigs2_per_k,
                                            n_real_bins_per_k, masks_per_k);
                    } else if (skch_writer && !item.sketch.signature1.empty()) {
                        const auto& sk = item.sketch;
                        skch_writer->add(item.genome_id, sk.signature1, sk.signature2,
                                         sk.n_real_bins, sk.genome_length,
                                         sk.real_bins_bitmask);
                    }
                }
                for (const auto& [k, v] : item.record.extra_fields)
                    if (k == "taxonomy") { taxonomy_map.emplace(item.record.accession, v); break; }
                parse_fasta_contig_accessions(item.fasta, [&](std::string_view contig_acc) {
                    cidx_writer.add(contig_acc, static_cast<uint32_t>(item.genome_id));
                });
                std::string().swap(item.fasta);  // free ~4 MB FASTA — last use
                meta_out << item.record.accession << "\t" << item.genome_id;
                for (const auto& [k, v] : item.record.extra_fields) meta_out << "\t" << v;
                meta_out << "\n";
            }
            // GSTX/GCOV/QUAL: dispatch off the main thread so worker threads stay busy.
            // Cross-chunk genera share mutable state (accum_holder/accum_member_ids) across
            // consecutive flushes of the same bucket, so run those synchronously to avoid
            // races. Chunk-local genera (accum_holder==nullptr) are fully independent and
            // safe to parallelize.
            if (accum_holder != nullptr) {
                // Cross-chunk: run inline so the counter accumulation is serialized.
                finalize_genus(std::move(buf), accum_holder, accum_member_ids);
            } else {
                // Chunk-local: cap in-flight tasks at 2×n_workers, then dispatch async.
                reap_futs();
                while (genus_futs.size() >= static_cast<size_t>(std::max(4ul, static_cast<size_t>(n_workers) / 4))) {
                    genus_futs.front().wait();
                    reap_futs();
                }
                auto buf_moved = std::move(buf);
                genus_futs.push_back(std::async(std::launch::async,
                    [buf_moved = std::move(buf_moved), &finalize_genus]() mutable {
                        finalize_genus(std::move(buf_moved), nullptr, nullptr);
                    }));
            }


            launch_shard_freeze();
            buf.clear();
            staging_raw_bytes = 0;
        };

        // Free a genus's cross-chunk bucket state once all its members have drained.
        // PCORE is handled by pcore_spiller.finalize() after all genomes drain.
        auto emit_core_pcore = [&](const std::string& /*key*/, TaxonBucket& b) {
            if (b.model_emitted) return;
            b.model_emitted = true;
            b.genus_counter.reset();
            std::vector<uint64_t>().swap(b.member_ids);
        };

        // Drain completion queue until all genomes processed
        for (size_t remaining = total; remaining > 0; --remaining) {
            Done d;
            {
                std::unique_lock lk(done_mx);
                done_cv.wait(lk, [&]{ return !done_q.empty(); });
                d = std::move(done_q.front());
                done_q.pop();
            }
            done_push_cv.notify_one();
            if (d.item) {
                ++n_done;
                if (cfg.taxonomy_group) {
                    // Find taxonomy string from extra_fields
                    std::string_view tax;
                    for (const auto& [k, v] : d.item->record.extra_fields)
                        if (k == "taxonomy") { tax = v; break; }
                    std::string key = extract_taxonomy_bucket(tax, cfg.taxonomy_rank.empty()
                                                              ? 'g' : cfg.taxonomy_rank[0]);
                    auto& bucket = taxon_buckets[key];
                    const uint64_t genome_bytes = d.item->fasta.size();
                    const uint64_t aamer_b      = d.item->aamers.size() * sizeof(uint64_t);
                    bucket.raw_bytes   += genome_bytes;
                    bucket.aamer_bytes += aamer_b;
                    taxon_total_bytes  += genome_bytes;
                    taxon_total_aamer  += aamer_b;
                    bucket.items.push_back(std::move(*d.item));
                    ++bucket.arrived;
                    auto drop = [&](TaxonBucket& b) {            // flush a bucket + clear its byte accounting
                        taxon_total_bytes -= b.raw_bytes; taxon_total_aamer -= b.aamer_bytes;
                        flush_staging_buf(b.items, &b.genus_counter, &b.member_ids);
                        b.raw_bytes = 0; b.aamer_bytes = 0;
                    };
                    const bool genus_done = cfg.build_pcore
                        && genus_total.count(key) && bucket.arrived >= genus_total[key]
                        && genus_total[key] >= cc_micro_threshold;
                    // Micro-genus complete: all expected genomes arrived, no CORE/PCORE model.
                    // Fold immediately into the cross-genus bin-pack buffer to free the memory
                    // rather than letting it accumulate until the end-of-build sweep.
                    const bool micro_done = cc_micro_threshold > 0
                        && genus_total.count(key) && bucket.arrived >= genus_total[key]
                        && genus_total[key] < cc_micro_threshold;
                    if (genus_done) {
                        // All members of a modeled genus have arrived → flush the
                        // remainder (folds it in) and emit one full-density CORE/PCORE.
                        drop(bucket);
                        emit_core_pcore(key, bucket);
                    } else if (micro_done) {
                        taxon_total_bytes -= bucket.raw_bytes;
                        taxon_total_aamer -= bucket.aamer_bytes;
                        if (!micro_combine.empty() && micro_combine_bytes + bucket.raw_bytes
                                >= cfg.shard_cfg.max_shard_size_bytes) {
                            flush_staging_buf(micro_combine);
                            micro_combine_bytes = 0;
                        }
                        micro_combine_bytes += bucket.raw_bytes;
                        for (auto& item : bucket.items) micro_combine.push_back(std::move(item));
                        bucket.raw_bytes = 0; bucket.aamer_bytes = 0;
                        bucket.items.clear(); bucket.items.shrink_to_fit();
                    } else if (bucket.raw_bytes >= cfg.shard_cfg.max_shard_size_bytes) {
                        drop(bucket);
                    } else if (taxon_total_aamer >= taxon_aamer_cap || taxon_total_bytes >= taxon_global_cap) {
                        // Over a RAM cap → flush the bucket holding the most RAM (aamers
                        // dominate, so weight by aamer_bytes) to free it.
                        auto it = std::max_element(taxon_buckets.begin(), taxon_buckets.end(),
                            [](const auto& a, const auto& b){
                                return (a.second.aamer_bytes + a.second.raw_bytes)
                                     < (b.second.aamer_bytes + b.second.raw_bytes); });
                        // When the largest remaining bucket is a micro-genus, all large genera
                        // have already been flushed.  Batch-flush every COMPLETE micro-genus in
                        // sorted-key order (deterministic packing) instead of one-at-a-time,
                        // which can never drain fast enough since each is well below the cap.
                        const bool max_is_micro = cc_micro_threshold > 0
                            && genus_total.count(it->first)
                            && genus_total[it->first] < cc_micro_threshold;
                        if (max_is_micro) {
                            std::vector<std::pair<std::string, uint64_t>> ready;
                            for (auto& [mk, mb] : taxon_buckets) {
                                if (mb.items.empty()) continue;
                                if (!genus_total.count(mk) || genus_total[mk] >= cc_micro_threshold) continue;
                                if (mb.arrived < genus_total[mk]) continue;
                                ready.push_back({mk, mb.raw_bytes});
                            }
                            if (!ready.empty()) {
                                std::sort(ready.begin(), ready.end(),
                                    [](const auto& a, const auto& b){ return a.first < b.first; });
                                std::vector<ChunkItem> mcomb;
                                uint64_t mcomb_bytes = 0;
                                for (auto& [mk2, gbytes] : ready) {
                                    auto& mb2 = taxon_buckets[mk2];
                                    if (!mcomb.empty() && mcomb_bytes + gbytes >= cfg.shard_cfg.max_shard_size_bytes) {
                                        flush_staging_buf(mcomb);
                                        mcomb_bytes = 0;
                                    }
                                    taxon_total_bytes -= mb2.raw_bytes;
                                    taxon_total_aamer -= mb2.aamer_bytes;
                                    for (auto& item : mb2.items) mcomb.push_back(std::move(item));
                                    mcomb_bytes += gbytes;
                                    mb2.raw_bytes = 0; mb2.aamer_bytes = 0;
                                    mb2.items.clear(); mb2.items.shrink_to_fit();
                                }
                                if (!mcomb.empty()) flush_staging_buf(mcomb);
                            } else {
                                // No complete micro-genera yet; flush the largest incomplete one.
                                drop(it->second);
                            }
                        } else {
                            drop(it->second);
                        }
                    }
                } else {
                    staging_raw_bytes += d.item->fasta.size();
                    staging_buffer.push_back(std::move(*d.item));
                    if (staging_raw_bytes >= cfg.shard_cfg.max_shard_size_bytes)
                        flush_staging_buf(staging_buffer);
                }
            } else {
                ++n_failed;
            }
            if (cfg.verbose || n_done % 50000 == 0)
                spdlog::info("Build: {}/{} genomes ({:.1f}%) | {} shards | {} failed",
                             n_done, original_total_records,
                             100.0 * n_done / std::max<size_t>(size_t(1), original_total_records),
                             current_shard_id, n_failed);
        }
        // Flush any micro-genera that completed mid-build but didn't fill the shard.
        if (!micro_combine.empty()) flush_staging_buf(micro_combine);

        // Final flush: remaining items in all buffers.
        // Micro-genera (below threshold) are bin-packed together to avoid
        // thousands of tiny shards from singleton-genus inputs.
        if (cfg.taxonomy_group) {
            const uint32_t micro_genome_threshold =
                resolve_micro_threshold(cfg.micro_genus_threshold, original_total_records);
            if (cfg.micro_genus_threshold == 0)
                spdlog::info("micro-genus-threshold: auto -> {} ({} genomes)",
                             micro_genome_threshold, original_total_records);

            // Collect micro-genus keys + their raw_bytes for bin-packing, without moving
            // all ChunkItems into one vector (which would spike peak RSS by n_micro * shard_size).
            std::vector<std::pair<std::string, uint64_t>> micro_keys; // (key, raw_bytes)
            uint32_t n_micro = 0;
            for (auto& [key, bucket] : taxon_buckets) {
                if (bucket.items.empty()) continue;
                if (static_cast<uint32_t>(bucket.items.size()) < micro_genome_threshold) {
                    micro_keys.push_back({key, bucket.raw_bytes});
                    ++n_micro;
                } else {
                    taxon_total_bytes -= bucket.raw_bytes; taxon_total_aamer -= bucket.aamer_bytes;
                    flush_staging_buf(bucket.items, &bucket.genus_counter, &bucket.member_ids);
                    bucket.raw_bytes = 0; bucket.aamer_bytes = 0;
                    emit_core_pcore(key, bucket);
                }
            }
            // Safety: emit CORE/PCORE for any modeled genus not closed during the drain.
            for (auto& [k2, b2] : taxon_buckets) emit_core_pcore(k2, b2);

            // Deterministic ordering → reproducible packing across builds.
            // Sort only the keys; items stay in taxon_buckets until consumed below.
            std::sort(micro_keys.begin(), micro_keys.end(),
                [](const auto& a, const auto& b) { return a.first < b.first; });

            // Greedy first-fit: pack micro-genera one-at-a-time (bounded combined buffer)
            // rather than materialising all items at once.
            std::vector<ChunkItem> combined;
            uint64_t combined_bytes = 0;
            for (auto& [key, genus_bytes] : micro_keys) {
                auto& bucket = taxon_buckets[key];
                if (!combined.empty() && combined_bytes + genus_bytes >= cfg.shard_cfg.max_shard_size_bytes) {
                    flush_staging_buf(combined);
                    combined_bytes = 0;
                }
                for (auto& item : bucket.items) combined.push_back(std::move(item));
                combined_bytes += genus_bytes;
                bucket.items.clear();
                bucket.items.shrink_to_fit();
            }
            if (!combined.empty()) flush_staging_buf(combined);

            if (n_micro > 0)
                spdlog::info("Packed {} micro-genera into combined shards", n_micro);
        } else {
            if (!staging_buffer.empty()) flush_staging_buf(staging_buffer);
        }

        producer.join();
        for (auto& w : workers) w.join();
        pending.clear();
        pending_input_rows.clear();
        {
            std::unique_lock lk(write_q_mx);
            writer_done = true;
        }
        write_q_cv.notify_one();
        io_writer.join();

        // Drain all in-flight genus-finalization tasks before writing metadata sections.
        for (auto& f : genus_futs) f.wait();
        genus_futs.clear();

        // Sort + scan all PCORE spill partitions; emit one add_sorted() per genus.
        if (cfg.build_pcore && !pcore_spiller.empty()) {
            spdlog::info("build: PCORE — sorting {} spill partitions…",
                         PcoreSpillWriter::N_PARTS);
            pcore_spiller.finalize(std::max(1, static_cast<int>(n_workers) / 8),
                [&](uint64_t gkh, std::vector<uint64_t>& aamers,
                    std::vector<uint32_t>& counts, uint32_t n_mem) {
                    std::lock_guard<std::mutex> lk(core_writer_mx);
                    build_pcore_w.add_sorted(gkh, aamers, counts, n_mem);
                });
        }

        std::fclose(ckpt_meta_file);
        ckpt_meta_file = nullptr;

        spdlog::info("Build: {}/{} genomes ({:.1f}%) | {} shards | {} failed",
                     n_done, original_total_records,
                     100.0 * n_done / std::max<size_t>(size_t(1), original_total_records),
                     current_shard_id, n_failed);
        spdlog::info("Wrote {} shards, {} genomes ({} failed)", current_shard_id,
                     catalog_rows.size(), n_failed);

        // Flush all shard writes to NFS, then write all metadata to a local
        // temp file. NFS write-back + ENOSPC can silently corrupt metadata,
        // leaving a garbage TailLocator that makes the archive unreadable by
        // merge_archives. Writing locally guarantees integrity.
        app_writer.flush();

        const uint64_t meta_base = app_writer.current_offset();
        auto meta_tmp_path = gpk_path_.parent_path() / "gpk_bld_meta_XXXXXX";
        std::string meta_tmp_str = meta_tmp_path.string();
        std::vector<char> bld_meta_tmp(meta_tmp_str.begin(), meta_tmp_str.end());
        bld_meta_tmp.push_back('\0');
        {
            int fd = ::mkstemp(bld_meta_tmp.data());
            if (fd < 0) throw std::runtime_error("Cannot create temp metadata file");
            ::close(fd);
        }
        AppendWriter mw;
        mw.create(bld_meta_tmp.data());
        mw.seek_to(meta_base);  // so section.file_offset values match the NFS file

        // Write CATL section
        uint64_t catalog_root_id = 0;
        {
            CatalogSectionWriter csw;
            for (const auto& m : catalog_rows)
                csw.add(m);
            SectionDesc catl_sd = csw.finalize(mw, next_section_id++);
            catalog_root_id = catl_sd.section_id;
            toc.add_section(catl_sd);
        }

        // Write GIDX section (O(1) genome_id lookup index)
        {
            GidxWriter gw;
            for (uint64_t i = 0; i < gidx_infos.size(); ++i) {
                const auto& gi = gidx_infos[i];
                auto it = shard_id_to_section_id.find(gi.shard_id);
                uint32_t sec_id = (it != shard_id_to_section_id.end())
                    ? static_cast<uint32_t>(it->second) : 0;
                gw.add(gi.genome_id, sec_id, gi.dir_index, i);
            }
            SectionDesc gidx_sd = gw.finalize(mw, next_section_id++);
            toc.add_section(gidx_sd);
        }

        // Write ACCX section
        uint64_t accession_root_id = 0;
        {
            AccessionIndexWriter aiw;
            for (const auto& [accession, genome_id] : accession_pairs)
                aiw.add(accession, genome_id);
            SectionDesc accx_sd = aiw.finalize(mw, next_section_id++);
            accession_root_id = accx_sd.section_id;
            toc.add_section(accx_sd);
        }

        // Write CIDX section (contig accession → genome_id index)
        if (cfg.build_cidx && cidx_writer.size() > 0) {
            spdlog::info("Writing CIDX: {} contig accessions", cidx_writer.size());
            SectionDesc cidx_sd = cidx_writer.finalize(mw, next_section_id++, /*batch_id=*/0);
            toc.add_section(cidx_sd);
        }

        // Write TAXN section (if taxonomy data available)
        if (!taxonomy_map.empty()) {
            TaxonomyIndexWriter tiw;
            for (const auto& [accession, taxonomy] : taxonomy_map)
                tiw.add(accession, taxonomy);
            SectionDesc taxn_sd = tiw.finalize(mw, next_section_id++);
            toc.add_section(taxn_sd);

            std::unordered_map<std::string, std::array<float, 136>> acc_kmer_profiles;
            if (!kmer_pairs.empty()) {
                std::unordered_map<GenomeId, const std::string*> gid_to_acc;
                gid_to_acc.reserve(accession_pairs.size());
                for (const auto& [acc, gid] : accession_pairs)
                    gid_to_acc[gid] = &acc;
                acc_kmer_profiles.reserve(kmer_pairs.size());
                for (const auto& [gid, prof] : kmer_pairs) {
                    auto it = gid_to_acc.find(gid);
                    if (it != gid_to_acc.end())
                        acc_kmer_profiles[*it->second] = prof;
                }
            }

            TxdbWriter txw;
            for (const auto& [accession, taxonomy] : taxonomy_map)
                txw.add(accession, taxonomy);
            if (!acc_kmer_profiles.empty())
                txw.set_kmer_profiles(acc_kmer_profiles);
            SectionDesc txdb_sd = txw.finalize(mw, next_section_id++);
            toc.add_section(txdb_sd);
        }

        // Write KMRX section
        if (!kmer_pairs.empty()) {
            KmrxWriter kw;
            for (auto& [gid, prof] : kmer_pairs)
                kw.add(gid, prof);
            SectionDesc kmrx_sd = kw.finalize(mw, next_section_id++);
            toc.add_section(kmrx_sd);
        }

        // Invariant: every catalogued genome must have a sketch row. A resume that
        // restored catalog entries without re-sketching would otherwise ship silently
        // sketchless genomes (resolvable via ACCX/GIDX but absent from SKCH).
        if (cfg.build_sketch) {
            const size_t skch_rows = skch_writer_mk ? skch_writer_mk->n_added()
                                   : (skch_writer ? skch_writer->n_added() : catalog_rows.size());
            if (skch_rows != catalog_rows.size())
                throw std::runtime_error("SKCH/catalog mismatch: " + std::to_string(skch_rows)
                    + " sketch rows vs " + std::to_string(catalog_rows.size())
                    + " catalog genomes — sketchless genomes would ship");
        }

        // Write SKCH section (OPH sketches)
        if (cfg.build_sketch) {
            if (skch_writer_mk) {
                std::string klist;
                for (size_t i = 0; i < cfg.sketch_kmer_sizes.size(); ++i) {
                    if (i) klist += ',';
                    klist += std::to_string(cfg.sketch_kmer_sizes[i]);
                }
                spdlog::info("Writing SKCH (v2 multi-k): {} sketches (k=[{}], m={})",
                             catalog_rows.size(), klist, cfg.sketch_size);
                SectionDesc skch_sd = skch_writer_mk->finalize(mw, next_section_id++);
                toc.add_section(skch_sd);
                skch_writer_mk.reset();
            } else if (skch_writer) {
                spdlog::info("Writing SKCH: {} sketches (k={}, m={})",
                             catalog_rows.size(), cfg.sketch_kmer_size, cfg.sketch_size);
                SectionDesc skch_sd = skch_writer->finalize(mw, next_section_id++);
                toc.add_section(skch_sd);
                skch_writer.reset();
            }
        }

        // Write GSTX section (genus sketch stats — enables O(1) check queries)
        if (cfg.build_gstx && multi_k_sketch && gstx_writer.n_genera() > 0) {
            spdlog::info("Writing GSTX: {} genera", gstx_writer.n_genera());
            SectionDesc gstx_sd = gstx_writer.finalize(mw, next_section_id++);
            toc.add_section(gstx_sd);
        }

        // Write QCOL section (intrinsic quality, columnar — replaces flat QUAL).
        if (cfg.build_gstx && qual_writer.size() > 0) {
            spdlog::info("Writing QCOL: {} genomes", qual_writer.size());
            SectionDesc qcol_sd = qcol_write(mw, next_section_id++, qual_writer.take());
            toc.add_section(qcol_sd);
        }

        // Write GCOV/FCOV/FMHR into the SAME metadata block so the archive has
        // exactly one footer generation (no post-TailLocator NFS rewrites — P11).
        // Profiles were computed in-memory during flush_staging_buf.
        if (cfg.build_gcov && cfg.build_gstx) {
            if (build_gcov_w.n_genera() > 0) {
                SectionDesc sd = build_gcov_w.finalize(mw, next_section_id++, SEC_GCOV);
                toc.add_section(sd);
                spdlog::info("GCOV: {} genera", build_gcov_w.n_genera());
            }
            if (build_fcov_w.n_genera() > 0) {
                SectionDesc sd = build_fcov_w.finalize(mw, next_section_id++, SEC_FCOV);
                toc.add_section(sd);
                spdlog::info("FCOV: {} genera", build_fcov_w.n_genera());
            }
            if (build_fmhr_w.n_genera() > 0) {
                SectionDesc sd = build_fmhr_w.finalize(mw, next_section_id++,
                                                       cfg.fmh_k, cfg.fmh_c);
                toc.add_section(sd);
                spdlog::info("FMHR: {} genera", build_fmhr_w.n_genera());
            }
            if (cfg.build_pcore && build_pcore_w.n_entries() > 0) {
                SectionDesc sd = build_pcore_w.finalize(mw, next_section_id++, build_pcore_frac_max);
                toc.add_section(sd);
                spdlog::info("PCORE: {} genus references (dense), model_hash={:#018x}",
                             build_pcore_w.n_entries(), sd.aux0);
            }

        }

        // BPRM — self-describing build parameters (mandatory, one per archive).
        // Written into the same metadata block so it is checksummed and covered
        // by the TOC/MetaBundle. Records the concrete param VALUES every
        // param-bearing section was built with (kills hardcoded literals).
        {
            BprmWriter bw(make_bprm_header_from_cfg(cfg));
            SectionDesc bprm_sd = bw.finalize(mw, next_section_id++);
            toc.add_section(bprm_sd);
            spdlog::info("BPRM: params_hash={:#018x}", bw.params_hash());
        }

        // Populate per-section content checksums (P1) for the metadata sections.
        // SHRD sections were hashed at write time; everything else lives in the
        // local temp file and is hashed here (streamed, no size cap), before the
        // TOC / MetaBundle (which cover these descriptors) are finalized.
        mw.flush();
        stamp_section_checksums(bld_meta_tmp.data(), toc.sections());

        // Content-addressed derivation (§5): now that every section's
        // content_xxh128 is populated, fold each non-SOURCE section's params +
        // upstream content hashes into derivation_hash + semantic_schema_hash.
        // Both the MetaBundle and the legacy TOC are written from toc.sections()
        // inside finalize(), so populating here covers both. Fresh build ⇒ no
        // deletions ⇒ deletion_set_hash = 0.
        populate_derivation(toc.sections(),
                            make_bprm_header_from_cfg(cfg).params_hash,
                            /*deletion_set_hash=*/0);

        // Write TOC + TailLocator to local file
        toc.finalize(mw,
                     /*generation=*/1,
                     /*live_count=*/catalog_rows.size(),
                     /*total_count=*/catalog_rows.size(),
                     /*prev_toc_offset=*/0,
                     catalog_root_id,
                     accession_root_id,
                     /*tombstone_root_id=*/0);
        mw.flush();  // local flush — always reliable

        // Copy metadata from local file to NFS using O_SYNC in 64MB chunks
        {
            const uint64_t meta_size = mw.current_offset() - meta_base;
            const size_t   CHUNK     = 64 * 1024 * 1024;
            int local_fd = ::open(bld_meta_tmp.data(), O_RDONLY);
            int nfs_fd   = ::open(partial_path.c_str(), O_WRONLY | O_SYNC);
            if (local_fd < 0 || nfs_fd < 0) {
                if (local_fd >= 0) ::close(local_fd);
                if (nfs_fd   >= 0) ::close(nfs_fd);
                throw std::runtime_error("Cannot open files for metadata copy");
            }
            std::vector<uint8_t> buf(CHUNK);
            uint64_t copied = 0;
            while (copied < meta_size) {
                size_t  chunk = static_cast<size_t>(std::min<uint64_t>(CHUNK, meta_size - copied));
                off_t   soff  = static_cast<off_t>(meta_base + copied);
                ssize_t nr    = ::pread(local_fd, buf.data(), chunk, soff);
                if (nr <= 0) {
                    ::close(local_fd); ::close(nfs_fd);
                    throw std::runtime_error("Local metadata read failed");
                }
                size_t  wr = 0;
                int     retries = 0;
                while (wr < static_cast<size_t>(nr)) {
                    ssize_t nw = ::pwrite(nfs_fd, buf.data() + wr,
                                          static_cast<size_t>(nr) - wr,
                                          soff + static_cast<off_t>(wr));
                    if (nw < 0) {
                        if ((errno == ENOSPC || errno == EIO) && retries < 60) {
                            spdlog::warn("NFS metadata write ENOSPC, retry {}…", ++retries);
                            std::this_thread::sleep_for(std::chrono::seconds(5));
                            continue;
                        }
                        ::close(local_fd); ::close(nfs_fd);
                        throw std::runtime_error("NFS metadata write failed: " +
                                                 std::string(std::strerror(errno)));
                    }
                    wr += static_cast<size_t>(nw);
                }
                copied += static_cast<uint64_t>(nr);
            }
            ::fsync(nfs_fd);
            ::close(local_fd);
            ::close(nfs_fd);
        }
        ::unlink(bld_meta_tmp.data());

        // ── Atomic publish ───────────────────────────────────────────────────
        // The archive was built into <out>.partial (shards + single-generation
        // metadata block). Make it durable, then atomically rename to the final
        // path and fsync the parent directory so the rename itself survives a
        // crash. A crash before the rename leaves <out>.partial + checkpoints
        // for resume; the final .gpk only ever appears fully sealed (P19).
        app_writer.flush();
        app_writer.close();
        std::filesystem::rename(partial_path, gpk_path_);
        {
            int dfd = ::open(gpk_path_.parent_path().c_str(), O_RDONLY | O_DIRECTORY);
            if (dfd >= 0) { ::fsync(dfd); ::close(dfd); }
        }

        // Remove checkpoint files — build completed successfully
        std::error_code ec;
        std::filesystem::remove(ckpt_path, ec);
        std::filesystem::remove(ckpt_meta_path, ec);

        spdlog::info("genopack archive written: {}", gpk_path_.string());
    }
};

ArchiveBuilder::ArchiveBuilder(const std::filesystem::path& dir, Config cfg)
    : impl_(std::make_unique<Impl>(dir, cfg))
{}
ArchiveBuilder::~ArchiveBuilder() = default;

void ArchiveBuilder::add_from_tsv(const std::filesystem::path& tsv_path) {
    impl_->add_from_tsv(tsv_path);
}
void ArchiveBuilder::add_from_gpk(const std::filesystem::path& source) {
    impl_->add_from_gpk(source);
}
void ArchiveBuilder::add(const BuildRecord& rec) { impl_->add(rec); }
void ArchiveBuilder::finalize() { impl_->finalize(); }

} // namespace genopack
