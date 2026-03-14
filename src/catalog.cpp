#include <genopack/catalog.hpp>
#include <algorithm>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace genopack {

// ── CatalogWriter::Impl ───────────────────────────────────────────────────────

struct CatalogWriter::Impl {
    std::filesystem::path path;
    std::vector<GenomeMeta> rows;

    explicit Impl(const std::filesystem::path& p) : path(p) {}

    void finalize() {
        // Sort by taxon_id for binary-search taxonomy lookup
        std::sort(rows.begin(), rows.end(),
                  [](const GenomeMeta& a, const GenomeMeta& b) {
                      return a.taxon_id < b.taxon_id;
                  });

        std::ofstream f(path, std::ios::binary | std::ios::trunc);
        if (!f) throw std::runtime_error("Cannot write catalog: " + path.string());

        const uint32_t n_rows = static_cast<uint32_t>(rows.size());
        const uint32_t rg_size = 4096;
        const uint32_t n_row_groups = (n_rows + rg_size - 1) / rg_size;

        // Header (64 bytes)
        f.write(reinterpret_cast<const char*>(&GPKC_MAGIC), 4);
        uint16_t version = FORMAT_VERSION;
        f.write(reinterpret_cast<const char*>(&version), 2);
        uint16_t pad0 = 0;
        f.write(reinterpret_cast<const char*>(&pad0), 2);
        uint64_t n64 = n_rows;
        f.write(reinterpret_cast<const char*>(&n64), 8);
        f.write(reinterpret_cast<const char*>(&n_row_groups), 4);
        // Placeholder offsets (we'll fix if needed; v1 = sequential write)
        uint64_t placeholder = 0;
        for (int i = 0; i < 5; ++i)
            f.write(reinterpret_cast<const char*>(&placeholder), 8);
        // pad to 64 bytes
        char header_pad[4] = {};
        f.write(header_pad, 4);

        // Column arrays (SoA)
        // Each column: n_rows elements, 8-byte-aligned per column

        auto write_col_u64 = [&](auto field) {
            for (const auto& r : rows) {
                uint64_t v = static_cast<uint64_t>(r.*field);
                f.write(reinterpret_cast<const char*>(&v), 8);
            }
        };
        auto write_col_u32 = [&](auto field) {
            for (const auto& r : rows) {
                uint32_t v = static_cast<uint32_t>(r.*field);
                f.write(reinterpret_cast<const char*>(&v), 4);
            }
            // Pad to 8-byte alignment for next column
            if (n_rows % 2 != 0) {
                uint32_t pad = 0;
                f.write(reinterpret_cast<const char*>(&pad), 4);
            }
        };
        auto write_col_u16 = [&](auto field) {
            for (const auto& r : rows) {
                uint16_t v = static_cast<uint16_t>(r.*field);
                f.write(reinterpret_cast<const char*>(&v), 2);
            }
            // Pad to 8-byte alignment
            size_t rem = (n_rows * 2) % 8;
            if (rem != 0) {
                char p[8] = {};
                f.write(p, 8 - rem);
            }
        };

        write_col_u64(&GenomeMeta::genome_id);
        write_col_u32(&GenomeMeta::taxon_id);
        write_col_u32(&GenomeMeta::shard_id);
        write_col_u64(&GenomeMeta::blob_offset);
        write_col_u32(&GenomeMeta::blob_len_cmp);
        write_col_u32(&GenomeMeta::blob_len_raw);
        write_col_u64(&GenomeMeta::genome_length);
        write_col_u32(&GenomeMeta::n_contigs);
        write_col_u16(&GenomeMeta::gc_pct_x100);
        write_col_u16(&GenomeMeta::completeness_x10);
        write_col_u16(&GenomeMeta::contamination_x10);
        write_col_u64(&GenomeMeta::oph_fingerprint);
        write_col_u32(&GenomeMeta::date_added);
        write_col_u32(&GenomeMeta::flags);

        // Row group stats
        for (uint32_t g = 0; g < n_row_groups; ++g) {
            uint32_t lo = g * rg_size;
            uint32_t hi = std::min(lo + rg_size, n_rows);
            RowGroupStats rgs{};
            rgs.first_row = lo;
            rgs.last_row  = hi - 1;
            rgs.taxon_id_min = rows[lo].taxon_id;
            rgs.taxon_id_max = rows[hi-1].taxon_id;
            rgs.genome_length_min = rows[lo].genome_length;
            rgs.genome_length_max = rows[lo].genome_length;
            rgs.completeness_min = rows[lo].completeness_x10;
            rgs.completeness_max = rows[lo].completeness_x10;
            for (uint32_t i = lo; i < hi; ++i) {
                rgs.genome_length_min = std::min(rgs.genome_length_min, rows[i].genome_length);
                rgs.genome_length_max = std::max(rgs.genome_length_max, rows[i].genome_length);
                rgs.completeness_min = std::min(rgs.completeness_min, rows[i].completeness_x10);
                rgs.completeness_max = std::max(rgs.completeness_max, rows[i].completeness_x10);
                rgs.flags_any |= rows[i].flags;
            }
            f.write(reinterpret_cast<const char*>(&rgs), sizeof(rgs));
        }
    }
};

CatalogWriter::CatalogWriter(const std::filesystem::path& path)
    : impl_(std::make_unique<Impl>(path))
{}
CatalogWriter::~CatalogWriter() = default;
void CatalogWriter::add(const GenomeMeta& meta) { impl_->rows.push_back(meta); }
void CatalogWriter::finalize() { impl_->finalize(); }

// ── CatalogReader (stub for v1 — full mmap implementation is next) ────────────

struct CatalogReader::Impl {
    std::vector<char>    file_data;
    std::vector<GenomeMeta> rows;  // deserialized for v1
    bool open_ = false;

    void open(const std::filesystem::path& path) {
        std::ifstream f(path, std::ios::binary | std::ios::ate);
        if (!f) throw std::runtime_error("Cannot open catalog: " + path.string());
        size_t sz = static_cast<size_t>(f.tellg());
        f.seekg(0);
        file_data.resize(sz);
        f.read(file_data.data(), sz);

        // Parse header
        if (sz < 64) throw std::runtime_error("Catalog too small");
        uint32_t magic;
        memcpy(&magic, file_data.data(), 4);
        if (magic != GPKC_MAGIC) throw std::runtime_error("Invalid catalog magic");

        uint64_t n_rows;
        memcpy(&n_rows, file_data.data() + 8, 8);

        // Read rows from SoA columns (same layout as written above)
        rows.resize(n_rows);
        size_t offset = 64;

        auto read_col_u64 = [&](uint64_t GenomeMeta::* field) {
            for (size_t i = 0; i < n_rows; ++i) {
                uint64_t v; memcpy(&v, file_data.data() + offset + i * 8, 8);
                rows[i].*field = v;
            }
            offset += n_rows * 8;
        };
        auto read_col_u32 = [&](uint32_t GenomeMeta::* field) {
            for (size_t i = 0; i < n_rows; ++i) {
                uint32_t v; memcpy(&v, file_data.data() + offset + i * 4, 4);
                rows[i].*field = v;
            }
            offset += n_rows * 4;
            if (n_rows % 2 != 0) offset += 4;  // alignment pad
        };
        auto read_col_u16 = [&](uint16_t GenomeMeta::* field) {
            for (size_t i = 0; i < n_rows; ++i) {
                uint16_t v; memcpy(&v, file_data.data() + offset + i * 2, 2);
                rows[i].*field = v;
            }
            offset += n_rows * 2;
            size_t rem = (n_rows * 2) % 8;
            if (rem != 0) offset += 8 - rem;
        };

        read_col_u64(&GenomeMeta::genome_id);
        read_col_u32(&GenomeMeta::taxon_id);
        read_col_u32(&GenomeMeta::shard_id);
        read_col_u64(&GenomeMeta::blob_offset);
        read_col_u32(&GenomeMeta::blob_len_cmp);
        read_col_u32(&GenomeMeta::blob_len_raw);
        read_col_u64(&GenomeMeta::genome_length);
        read_col_u32(&GenomeMeta::n_contigs);
        read_col_u16(&GenomeMeta::gc_pct_x100);
        read_col_u16(&GenomeMeta::completeness_x10);
        read_col_u16(&GenomeMeta::contamination_x10);
        read_col_u64(&GenomeMeta::oph_fingerprint);
        read_col_u32(&GenomeMeta::date_added);
        read_col_u32(&GenomeMeta::flags);

        open_ = true;
    }
};

CatalogReader::CatalogReader()  : impl_(std::make_unique<Impl>()) {}
CatalogReader::~CatalogReader() = default;

void CatalogReader::open(const std::filesystem::path& path)  { impl_->open(path); }
void CatalogReader::close() { impl_->file_data.clear(); impl_->rows.clear(); impl_->open_ = false; }
size_t CatalogReader::n_rows() const { return impl_->rows.size(); }
size_t CatalogReader::n_genomes() const {
    return std::count_if(impl_->rows.begin(), impl_->rows.end(),
                         [](const GenomeMeta& m){ return !m.is_deleted(); });
}

const GenomeMeta* CatalogReader::find_genome(GenomeId id) const {
    for (const auto& r : impl_->rows)
        if (r.genome_id == id) return &r;
    return nullptr;
}

std::vector<const GenomeMeta*> CatalogReader::for_taxon(TaxonId taxon_id) const {
    std::vector<const GenomeMeta*> out;
    // Rows are sorted by taxon_id; binary search for range
    auto lo = std::lower_bound(impl_->rows.begin(), impl_->rows.end(), taxon_id,
        [](const GenomeMeta& m, TaxonId v){ return m.taxon_id < v; });
    auto hi = std::upper_bound(lo, impl_->rows.end(), taxon_id,
        [](TaxonId v, const GenomeMeta& m){ return v < m.taxon_id; });
    for (auto it = lo; it != hi; ++it)
        if (!it->is_deleted()) out.push_back(&*it);
    return out;
}

std::vector<const GenomeMeta*> CatalogReader::for_taxon_range(TaxonId lo, TaxonId hi) const {
    std::vector<const GenomeMeta*> out;
    auto it_lo = std::lower_bound(impl_->rows.begin(), impl_->rows.end(), lo,
        [](const GenomeMeta& m, TaxonId v){ return m.taxon_id < v; });
    auto it_hi = std::lower_bound(it_lo, impl_->rows.end(), hi,
        [](const GenomeMeta& m, TaxonId v){ return m.taxon_id < v; });
    for (auto it = it_lo; it != it_hi; ++it)
        if (!it->is_deleted()) out.push_back(&*it);
    return out;
}

std::vector<const GenomeMeta*> CatalogReader::filter(const ExtractQuery& q) const {
    std::vector<const GenomeMeta*> out;
    for (const auto& r : impl_->rows) {
        if (!q.include_deleted && r.is_deleted()) continue;
        if (r.completeness_x10 < static_cast<uint16_t>(q.min_completeness * 10.0f)) continue;
        if (r.contamination_x10 > static_cast<uint16_t>(q.max_contamination * 10.0f)) continue;
        if (r.genome_length < q.min_genome_length) continue;
        if (r.genome_length > q.max_genome_length) continue;
        if (r.n_contigs > q.max_contigs) continue;
        out.push_back(&r);
        if (out.size() >= q.limit) break;
    }
    return out;
}

void CatalogReader::scan(std::function<bool(const GenomeMeta&)> predicate) const {
    for (const auto& r : impl_->rows)
        if (!predicate(r)) break;
}

CatalogReader::Stats CatalogReader::compute_stats(TaxonId taxon_id) const {
    Stats s{};
    double csum = 0, contam_sum = 0, gc_sum = 0;
    for (const auto& r : impl_->rows) {
        if (taxon_id != INVALID_TAXON_ID && r.taxon_id != taxon_id) continue;
        ++s.n_total;
        if (r.is_deleted()) { ++s.n_deleted; continue; }
        s.total_genome_length += r.genome_length;
        csum      += r.completeness_x10 / 10.0;
        contam_sum += r.contamination_x10 / 10.0;
        gc_sum    += r.gc_pct_x100 / 100.0;
    }
    size_t n_live = s.n_total - s.n_deleted;
    if (n_live > 0) {
        s.mean_completeness  = csum / n_live;
        s.mean_contamination = contam_sum / n_live;
        s.mean_gc            = gc_sum / n_live;
    }
    return s;
}

} // namespace genopack
