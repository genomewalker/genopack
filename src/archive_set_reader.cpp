#include "genopack/archive_set_reader.hpp"
#include "genopack/archive.hpp"
#include "genopack/txdb.hpp"

#include <algorithm>
#include <filesystem>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

namespace genopack {

namespace fs = std::filesystem;

namespace {

bool path_is_single_archive(const fs::path& p) {
    if (fs::exists(p / "toc.bin")) return true;
    if (p.extension() == ".gpk")   return true;
    return false;
}

std::vector<fs::path> discover_parts(const fs::path& dir) {
    std::vector<fs::path> out;
    if (!fs::is_directory(dir)) return out;
    for (auto& e : fs::directory_iterator(dir)) {
        const auto& ep = e.path();
        if (ep.extension() != ".gpk") continue;
        if (!(e.is_directory() || e.is_regular_file())) continue;
        const auto stem = ep.stem().string();
        if (stem.rfind("part_", 0) != 0) continue;
        out.push_back(ep);
    }
    std::sort(out.begin(), out.end(), [](const fs::path& a, const fs::path& b) {
        auto num = [](const fs::path& p) -> long long {
            const auto s = p.stem().string();
            try { return std::stoll(s.substr(5)); } catch (...) { return -1; }
        };
        long long na = num(a), nb = num(b);
        if (na != nb) return na < nb;
        return a < b;
    });
    return out;
}

} // namespace

struct ArchiveSetReader::Impl {
    std::vector<std::unique_ptr<ArchiveReader>> parts;
    std::vector<fs::path>                       part_paths;
    bool                                        multipart = false;

    // Lazy per-part TaxonomyTree cache. Vector is sized to parts.size() at open;
    // entries are populated on demand under tax_mu.
    mutable std::vector<std::optional<TaxonomyTree>> tax_cache;
    mutable std::vector<bool>                        tax_loaded;
    mutable std::mutex                                tax_mu;

    int find_owner_part(std::string_view accession) const {
        for (size_t i = 0; i < parts.size(); ++i) {
            if (parts[i]->genome_meta_by_accession(accession).has_value())
                return static_cast<int>(i);
        }
        return -1;
    }

    const TaxonomyTree* tree_for_part(size_t i) const {
        std::lock_guard<std::mutex> lk(tax_mu);
        if (i >= parts.size()) return nullptr;
        if (!tax_loaded[i]) {
            tax_cache[i]  = parts[i]->taxonomy_tree();
            tax_loaded[i] = true;
        }
        return tax_cache[i].has_value() ? &(*tax_cache[i]) : nullptr;
    }
};

ArchiveSetReader::ArchiveSetReader() : impl_(std::make_unique<Impl>()) {}
ArchiveSetReader::~ArchiveSetReader() = default;
ArchiveSetReader::ArchiveSetReader(ArchiveSetReader&&) noexcept = default;
ArchiveSetReader& ArchiveSetReader::operator=(ArchiveSetReader&&) noexcept = default;

void ArchiveSetReader::open(const fs::path& path) {
    close();
    if (!fs::exists(path))
        throw std::runtime_error("archive path does not exist: " + path.string());

    std::vector<fs::path> part_paths;
    bool multipart = false;

    if (fs::is_directory(path) && path_is_single_archive(path)) {
        part_paths.push_back(path);
    } else if (!fs::is_directory(path)) {
        part_paths.push_back(path);
    } else {
        part_paths = discover_parts(path);
        if (part_paths.empty())
            throw std::runtime_error("no toc.bin and no part_*.gpk entries under: " + path.string());
        multipart = true;
    }

    impl_->parts.reserve(part_paths.size());
    for (const auto& pp : part_paths) {
        auto ar = std::make_unique<ArchiveReader>();
        ar->open(pp);
        impl_->parts.push_back(std::move(ar));
    }
    impl_->part_paths = std::move(part_paths);
    impl_->multipart  = multipart;
    impl_->tax_cache.assign(impl_->parts.size(), std::nullopt);
    impl_->tax_loaded.assign(impl_->parts.size(), false);
}

void ArchiveSetReader::close() {
    impl_->parts.clear();
    impl_->part_paths.clear();
    impl_->tax_cache.clear();
    impl_->tax_loaded.clear();
    impl_->multipart = false;
}

bool ArchiveSetReader::is_open() const { return !impl_->parts.empty(); }
bool ArchiveSetReader::is_multipart() const { return impl_->multipart; }
size_t ArchiveSetReader::part_count() const { return impl_->parts.size(); }
const std::vector<fs::path>& ArchiveSetReader::part_paths() const { return impl_->part_paths; }

ArchiveReader::ArchiveStats ArchiveSetReader::archive_stats() const {
    ArchiveReader::ArchiveStats agg{};
    agg.generation = 0;
    for (const auto& p : impl_->parts) {
        auto s = p->archive_stats();
        agg.generation             = std::max(agg.generation, s.generation);
        agg.n_shards              += s.n_shards;
        agg.n_genomes_total       += s.n_genomes_total;
        agg.n_genomes_live        += s.n_genomes_live;
        agg.total_raw_bp          += s.total_raw_bp;
        agg.total_compressed_bytes += s.total_compressed_bytes;
    }
    agg.compression_ratio = (agg.total_compressed_bytes > 0)
        ? static_cast<double>(agg.total_raw_bp) / static_cast<double>(agg.total_compressed_bytes)
        : 0.0;
    return agg;
}

size_t ArchiveSetReader::count(const ExtractQuery& q) const {
    size_t n = 0;
    for (const auto& p : impl_->parts) n += p->count(q);
    return n;
}

std::vector<LocatedGenomeMeta> ArchiveSetReader::filter_meta(const ExtractQuery& q) const {
    std::vector<LocatedGenomeMeta> out;
    for (size_t i = 0; i < impl_->parts.size(); ++i) {
        auto rows = impl_->parts[i]->filter_meta(q);
        out.reserve(out.size() + rows.size());
        for (auto& m : rows)
            out.push_back({i, impl_->part_paths[i], std::move(m)});
    }
    return out;
}

std::vector<ExtractedGenome> ArchiveSetReader::extract(const ExtractQuery& q) const {
    if (!q.accessions.empty()) {
        auto fetched = batch_fetch_by_accessions(q.accessions);
        std::vector<ExtractedGenome> out;
        out.reserve(fetched.size());
        for (auto& f : fetched) if (f) out.push_back(std::move(*f));
        return out;
    }

    std::vector<ExtractedGenome> out;
    const size_t limit = (q.limit > 0) ? q.limit : 0;
    for (const auto& p : impl_->parts) {
        ExtractQuery sub = q;
        if (limit > 0) {
            if (out.size() >= limit) break;
            sub.limit = limit - out.size();
        }
        auto part_out = p->extract(sub);
        for (auto& g : part_out) {
            out.push_back(std::move(g));
            if (limit > 0 && out.size() >= limit) break;
        }
    }
    return out;
}

std::optional<ExtractedGenome>
ArchiveSetReader::fetch_by_accession(std::string_view accession) const {
    for (const auto& p : impl_->parts) {
        auto eg = p->fetch_by_accession(accession);
        if (eg) return eg;
    }
    return std::nullopt;
}

std::vector<std::optional<ExtractedGenome>>
ArchiveSetReader::batch_fetch_by_accessions(const std::vector<std::string>& accessions) const {
    std::vector<std::optional<ExtractedGenome>> out(accessions.size());

    if (impl_->parts.size() == 1) {
        return impl_->parts[0]->batch_fetch_by_accessions(accessions);
    }

    // Route each accession to its owning part (probe parts in order).
    std::vector<int> owner(accessions.size(), -1);
    std::vector<std::vector<size_t>> by_part(impl_->parts.size());
    for (size_t i = 0; i < accessions.size(); ++i) {
        for (size_t pi = 0; pi < impl_->parts.size(); ++pi) {
            if (impl_->parts[pi]->genome_meta_by_accession(accessions[i]).has_value()) {
                owner[i] = static_cast<int>(pi);
                by_part[pi].push_back(i);
                break;
            }
        }
    }

    for (size_t pi = 0; pi < impl_->parts.size(); ++pi) {
        const auto& idx = by_part[pi];
        if (idx.empty()) continue;
        std::vector<std::string> sub;
        sub.reserve(idx.size());
        for (size_t k : idx) sub.push_back(accessions[k]);
        auto got = impl_->parts[pi]->batch_fetch_by_accessions(sub);
        for (size_t j = 0; j < idx.size(); ++j)
            out[idx[j]] = std::move(got[j]);
    }
    return out;
}

std::optional<std::string>
ArchiveSetReader::fetch_sequence_slice_by_accession(std::string_view accession,
                                                   uint64_t start,
                                                   uint64_t length) const {
    for (const auto& p : impl_->parts) {
        auto s = p->fetch_sequence_slice_by_accession(accession, start, length);
        if (s) return s;
        // genome_meta_by_accession check would be redundant — slice returns
        // nullopt both for "not in this part" and "out of range"; we want to
        // continue probing parts in the former case. Use accession existence:
        if (p->genome_meta_by_accession(accession).has_value()) {
            // owned by this part; nullopt means out-of-range, do not retry others
            return std::nullopt;
        }
    }
    return std::nullopt;
}

std::optional<std::string>
ArchiveSetReader::taxonomy_for_accession(std::string_view accession) const {
    for (const auto& p : impl_->parts) {
        auto t = p->taxonomy_for_accession(accession);
        if (t) return t;
    }
    return std::nullopt;
}

void ArchiveSetReader::scan_taxonomy(
        const std::function<void(std::string_view, std::string_view)>& cb) const {
    for (const auto& p : impl_->parts) p->scan_taxonomy(cb);
}

std::optional<ArchiveSetReader::LocatedTaxonomy>
ArchiveSetReader::taxonomy_tree_for_accession(std::string_view accession) const {
    int owner = impl_->find_owner_part(accession);
    if (owner < 0) return std::nullopt;
    const TaxonomyTree* tr = impl_->tree_for_part(static_cast<size_t>(owner));
    if (!tr) return std::nullopt;
    return LocatedTaxonomy{ static_cast<size_t>(owner), tr };
}

ArchiveSetReader::TaxonomySummary ArchiveSetReader::taxonomy_summary() const {
    TaxonomySummary s{};
    std::unordered_set<uint64_t> taxids;
    for (size_t i = 0; i < impl_->parts.size(); ++i) {
        const TaxonomyTree* tr = impl_->tree_for_part(i);
        if (!tr) continue;
        ++s.n_parts_with_taxonomy;
        s.n_accessions += tr->n_accessions();
        tr->scan_nodes([&](uint64_t taxid, uint64_t /*parent*/,
                           TaxRank /*rank*/, std::string_view /*name*/,
                           bool /*is_synthetic*/) {
            taxids.insert(taxid);
        });
    }
    s.n_nodes_union = taxids.size();
    return s;
}

ArchiveSetReader open_archive_auto(const fs::path& path) {
    ArchiveSetReader r;
    r.open(path);
    return r;
}

} // namespace genopack
