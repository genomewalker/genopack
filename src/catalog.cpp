#include <genopack/catalog.hpp>
#include <algorithm>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace genopack {

// ── CatalogWriter::Impl ───────────────────────────────────────────────────────

struct CatalogWriter::Impl {
    std::filesystem::path path;
    std::vector<GenomeMeta> rows;

    explicit Impl(const std::filesystem::path& p) : path(p) {}

    void finalize() {
        // Sort by oph_fingerprint for locality-sensitive similarity lookup
        std::sort(rows.begin(), rows.end(),
                  [](const GenomeMeta& a, const GenomeMeta& b) {
                      return a.oph_fingerprint < b.oph_fingerprint;
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
        write_col_u32(&GenomeMeta::_reserved0);
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
            rgs.oph_min = rows[lo].oph_fingerprint;
            rgs.oph_max = rows[lo].oph_fingerprint;
            rgs.genome_length_min = rows[lo].genome_length;
            rgs.genome_length_max = rows[lo].genome_length;
            rgs.completeness_min = rows[lo].completeness_x10;
            rgs.completeness_max = rows[lo].completeness_x10;
            for (uint32_t i = lo; i < hi; ++i) {
                rgs.oph_min = std::min(rgs.oph_min, rows[i].oph_fingerprint);
                rgs.oph_max = std::max(rgs.oph_max, rows[i].oph_fingerprint);
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
        read_col_u32(&GenomeMeta::_reserved0);
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

std::vector<const GenomeMeta*> CatalogReader::for_oph_range(uint64_t lo, uint64_t hi) const {
    std::vector<const GenomeMeta*> out;
    auto it_lo = std::lower_bound(impl_->rows.begin(), impl_->rows.end(), lo,
        [](const GenomeMeta& m, uint64_t v){ return m.oph_fingerprint < v; });
    auto it_hi = std::upper_bound(it_lo, impl_->rows.end(), hi,
        [](uint64_t v, const GenomeMeta& m){ return v < m.oph_fingerprint; });
    for (auto it = it_lo; it != it_hi; ++it)
        if (!it->is_deleted()) out.push_back(&*it);
    return out;
}

std::vector<const GenomeMeta*> CatalogReader::filter(const ExtractQuery& q) const {
    std::vector<const GenomeMeta*> out;
    for (const auto& r : impl_->rows) {
        if (!q.include_deleted && r.is_deleted()) continue;
        if (r.oph_fingerprint < q.min_oph) continue;
        if (r.oph_fingerprint > q.max_oph) continue;
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

CatalogReader::Stats CatalogReader::compute_stats() const {
    Stats s{};
    double csum = 0, contam_sum = 0, gc_sum = 0;
    for (const auto& r : impl_->rows) {
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

// ── CatalogSectionWriter (v2) ─────────────────────────────────────────────────

CatalogSectionWriter::CatalogSectionWriter(size_t row_group_size)
    : row_group_size_(row_group_size)
{}

void CatalogSectionWriter::add(const GenomeMeta& m) {
    rows_.push_back(m);
}

void CatalogSectionWriter::flush_group(size_t first, size_t last) {
    RowGroupStatsV2 g{};
    g.first_row = static_cast<uint32_t>(first);
    g.last_row  = static_cast<uint32_t>(last);

    const GenomeMeta& r0 = rows_[first];
    g.oph_min            = r0.oph_fingerprint;
    g.oph_max            = r0.oph_fingerprint;
    g.genome_length_min  = r0.genome_length;
    g.genome_length_max  = r0.genome_length;
    g.completeness_min   = r0.completeness_x10;
    g.completeness_max   = r0.completeness_x10;
    g.contamination_min  = r0.contamination_x10;
    g.contamination_max  = r0.contamination_x10;
    g.live_count         = 0;
    g.flags_any          = 0;

    for (size_t i = first; i <= last; ++i) {
        const GenomeMeta& r = rows_[i];
        g.oph_min           = std::min(g.oph_min,           r.oph_fingerprint);
        g.oph_max           = std::max(g.oph_max,           r.oph_fingerprint);
        g.genome_length_min = std::min(g.genome_length_min, r.genome_length);
        g.genome_length_max = std::max(g.genome_length_max, r.genome_length);
        g.completeness_min  = std::min(g.completeness_min,  r.completeness_x10);
        g.completeness_max  = std::max(g.completeness_max,  r.completeness_x10);
        g.contamination_min = std::min(g.contamination_min, r.contamination_x10);
        g.contamination_max = std::max(g.contamination_max, r.contamination_x10);
        g.flags_any        |= r.flags;
        if (!r.is_deleted()) ++g.live_count;
    }

    groups_.push_back(g);
}

SectionDesc CatalogSectionWriter::finalize(AppendWriter& writer, uint64_t section_id) {
    const size_t n_rows   = rows_.size();
    const size_t n_groups = (n_rows + row_group_size_ - 1) / row_group_size_;

    // Build row group stats
    groups_.clear();
    groups_.reserve(n_groups);
    for (size_t g = 0; g < n_groups; ++g) {
        size_t first = g * row_group_size_;
        size_t last  = std::min(first + row_group_size_, n_rows) - 1;
        flush_group(first, last);
    }

    // Compute offsets within the section
    const uint64_t stats_offset = sizeof(CatlHeader);
    const uint64_t rows_offset  = stats_offset + groups_.size() * sizeof(RowGroupStatsV2);

    CatlHeader hdr{};
    hdr.magic          = SEC_CATL;
    hdr.n_rows         = static_cast<uint32_t>(n_rows);
    hdr.n_groups       = static_cast<uint32_t>(groups_.size());
    hdr.row_group_size = static_cast<uint32_t>(row_group_size_);
    hdr.stats_offset   = stats_offset;
    hdr.rows_offset    = rows_offset;

    const uint64_t section_start = writer.current_offset();

    writer.append(&hdr, sizeof(hdr));
    if (!groups_.empty())
        writer.append(groups_.data(), groups_.size() * sizeof(RowGroupStatsV2));
    if (!rows_.empty())
        writer.append(rows_.data(), rows_.size() * sizeof(GenomeMeta));

    const uint64_t section_end  = writer.current_offset();
    const uint64_t section_size = section_end - section_start;

    SectionDesc desc{};
    desc.type               = SEC_CATL;
    desc.version            = FORMAT_V2_MAJOR;
    desc.flags              = 0;
    desc.section_id         = section_id;
    desc.file_offset        = section_start;
    desc.compressed_size    = section_size;
    desc.uncompressed_size  = section_size;
    desc.item_count         = static_cast<uint64_t>(n_rows);
    desc.aux0               = static_cast<uint64_t>(groups_.size());
    desc.aux1               = static_cast<uint64_t>(row_group_size_);

    return desc;
}

// ── CatalogSectionReader (v2) ─────────────────────────────────────────────────

void CatalogSectionReader::open(const uint8_t* base, uint64_t section_offset, uint64_t section_size) {
    if (section_size < sizeof(CatlHeader))
        throw std::runtime_error("CatalogSectionReader: section too small");

    const uint8_t* sec = base + section_offset;
    header_ = reinterpret_cast<const CatlHeader*>(sec);

    if (header_->magic != SEC_CATL)
        throw std::runtime_error("CatalogSectionReader: bad magic");

    groups_ = reinterpret_cast<const RowGroupStatsV2*>(sec + header_->stats_offset);
    rows_   = reinterpret_cast<const GenomeMeta*>(sec + header_->rows_offset);
}

size_t CatalogSectionReader::n_rows() const {
    return header_ ? header_->n_rows : 0;
}

size_t CatalogSectionReader::n_groups() const {
    return header_ ? header_->n_groups : 0;
}

bool CatalogSectionReader::group_may_match(const RowGroupStatsV2& g, const ExtractQuery& q) const {
    // Prune if query excludes deleted and all rows in group are deleted
    if (!q.include_deleted && g.live_count == 0)
        return false;

    // oph range
    if (g.oph_max < q.min_oph) return false;
    if (g.oph_min > q.max_oph) return false;

    // completeness
    const uint16_t min_comp = static_cast<uint16_t>(q.min_completeness * 10.0f);
    if (g.completeness_max < min_comp) return false;

    // contamination
    const uint16_t max_contam = static_cast<uint16_t>(q.max_contamination * 10.0f);
    if (g.contamination_min > max_contam) return false;

    // genome length
    if (g.genome_length_max < q.min_genome_length) return false;
    if (g.genome_length_min > q.max_genome_length) return false;

    return true;
}

bool CatalogSectionReader::row_matches(const GenomeMeta& m, const ExtractQuery& q) const {
    if (!q.include_deleted && m.is_deleted()) return false;
    if (m.oph_fingerprint < q.min_oph) return false;
    if (m.oph_fingerprint > q.max_oph) return false;
    if (m.completeness_x10 < static_cast<uint16_t>(q.min_completeness * 10.0f)) return false;
    if (m.contamination_x10 > static_cast<uint16_t>(q.max_contamination * 10.0f)) return false;
    if (m.genome_length < q.min_genome_length) return false;
    if (m.genome_length > q.max_genome_length) return false;
    if (m.n_contigs > q.max_contigs) return false;
    return true;
}

std::vector<const GenomeMeta*> CatalogSectionReader::filter(const ExtractQuery& q) const {
    std::vector<const GenomeMeta*> out;
    if (!header_) return out;

    for (uint32_t gi = 0; gi < header_->n_groups; ++gi) {
        const RowGroupStatsV2& g = groups_[gi];
        if (!group_may_match(g, q)) continue;

        for (uint32_t ri = g.first_row; ri <= g.last_row; ++ri) {
            const GenomeMeta& m = rows_[ri];
            if (!row_matches(m, q)) continue;
            out.push_back(&m);
            if (out.size() >= q.limit) return out;
        }
    }
    return out;
}

const GenomeMeta* CatalogSectionReader::find_genome(GenomeId id) const {
    if (!header_) return nullptr;
    for (uint32_t i = 0; i < header_->n_rows; ++i)
        if (rows_[i].genome_id == id) return &rows_[i];
    return nullptr;
}

void CatalogSectionReader::scan(std::function<bool(const GenomeMeta&)> cb) const {
    if (!header_) return;
    for (uint32_t i = 0; i < header_->n_rows; ++i)
        if (!cb(rows_[i])) return;
}

// ── MergedCatalogReader ───────────────────────────────────────────────────────

void MergedCatalogReader::add_fragment(const uint8_t* base, uint64_t offset, uint64_t size) {
    CatalogSectionReader r;
    r.open(base, offset, size);
    fragments_.push_back(std::move(r));
}

std::vector<const GenomeMeta*> MergedCatalogReader::filter(const ExtractQuery& q) const {
    // Collect results from all fragments; newest fragment wins for duplicate genome_id.
    // fragments_ is ordered oldest-first, so iterate newest-first for dedup.
    std::unordered_map<GenomeId, const GenomeMeta*> seen;

    for (int fi = static_cast<int>(fragments_.size()) - 1; fi >= 0; --fi) {
        auto partial = fragments_[fi].filter(q);
        for (const GenomeMeta* m : partial) {
            // Only insert if not already seen (newer fragment already recorded it)
            seen.emplace(m->genome_id, m);
        }
    }

    std::vector<const GenomeMeta*> out;
    out.reserve(seen.size());
    for (auto& [id, ptr] : seen)
        out.push_back(ptr);

    if (out.size() > q.limit)
        out.resize(q.limit);

    return out;
}

const GenomeMeta* MergedCatalogReader::find_genome(GenomeId id) const {
    // Search newest-first
    for (int fi = static_cast<int>(fragments_.size()) - 1; fi >= 0; --fi) {
        const GenomeMeta* m = fragments_[fi].find_genome(id);
        if (m) return m;
    }
    return nullptr;
}

void MergedCatalogReader::scan(std::function<bool(const GenomeMeta&)> cb) const {
    // Scan all fragments oldest-first; deduplicate by genome_id (newest wins).
    // Two-pass: collect newest pointer per genome_id, then call cb in order.
    std::unordered_map<GenomeId, const GenomeMeta*> seen;

    for (int fi = static_cast<int>(fragments_.size()) - 1; fi >= 0; --fi) {
        fragments_[fi].scan([&](const GenomeMeta& m) -> bool {
            seen.emplace(m.genome_id, &m);
            return true;
        });
    }

    for (auto& [id, ptr] : seen)
        if (!cb(*ptr)) return;
}

size_t MergedCatalogReader::n_rows() const {
    size_t total = 0;
    for (const auto& f : fragments_) total += f.n_rows();
    return total;
}

} // namespace genopack
