#pragma once
#include <genopack/archive.hpp>
#include <genopack/fmhr.hpp>
#include <genopack/gcov.hpp>
#include <genopack/gstx.hpp>
#include <algorithm>
#include <functional>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace genopack::check {

struct ICheckReader {
    virtual ~ICheckReader() = default;
    virtual void scan_taxonomy_with_id(
        const std::function<void(std::string_view acc,
                                 std::string_view tax,
                                 GenomeId gid)>& cb) const = 0;
    virtual std::optional<GenomeMeta> genome_meta_by_accession(std::string_view acc) const = 0;
    virtual const float* kmer_profile_by_accession(std::string_view acc) const = 0;
    virtual bool has_sketches() const = 0;
    virtual std::vector<uint32_t> available_kmer_sizes() const = 0;
    virtual uint32_t sketch_sketch_size() const = 0;
    virtual const GstxEntry* gstx_for_genus(std::string_view genus) const = 0;
    virtual const GstxReader* gstx_reader() const = 0;
    virtual const GcovEntry* gcov_for_genus(std::string_view genus) const = 0;
    virtual const GcovReader* gcov_reader() const = 0;
    virtual const GcovEntry* fcov_for_family(std::string_view family) const = 0;
    virtual const GcovReader* fcov_reader() const = 0;
    virtual FmhrView fmhr_for_genus(std::string_view genus) const = 0;
    virtual const FmhrReader* fmhr_reader() const = 0;
    virtual bool has_core() const = 0;
    virtual CoreView core_for_genus(std::string_view genus) const = 0;
    virtual const CoreReader* core_reader() const = 0;
    virtual bool has_fcore() const = 0;
    virtual CoreView core_for_family(std::string_view family) const = 0;
    virtual const CoreReader* fcore_reader() const = 0;
    virtual bool has_pcore() const = 0;
    virtual PcoreView pcore_for_genus(std::string_view genus) const = 0;
    virtual PcoreView pcore_for_family(std::string_view family) const = 0;
    virtual const PcoreReader* pcore_reader() const = 0;
    virtual bool has_gami_v2() const { return false; }
    virtual void load_gami_into(GlobalMultiplicityIndex& out) const { (void)out; }
    virtual void visit_sketch_batches_multi_k(
        const std::vector<std::string>& accs,
        const std::vector<uint32_t>& ks, uint32_t sz,
        const std::function<void(size_t idx, uint32_t k,
                                 const SketchResult& sk)>& cb) const = 0;

    // Single-pass sequential scan: sorted_gids must be sorted ascending.
    // cb(idx_in_sorted_gids, ki, result) — one frame decompression per frame, all k values.
    virtual void visit_sketch_gids(
        const std::vector<GenomeId>& sorted_gids,
        uint32_t sz,
        const std::function<void(size_t idx, uint32_t ki,
                                 const SketchResult& sk)>& cb) const = 0;
    virtual void visit_shard_batches(
        const std::vector<std::string>& accs,
        const std::function<void(ArchiveReader::ShardBatch&)>& cb) const = 0;
    virtual void visit_shard_batches_parallel(
        const std::vector<std::string>& accs,
        int n_readers,
        const std::function<void(ArchiveReader::ShardBatch&)>& cb) const = 0;
};

class SingleArchiveCheckReader : public ICheckReader {
public:
    explicit SingleArchiveCheckReader(ArchiveReader& ar) : ar_(ar) {}

    void scan_taxonomy_with_id(
        const std::function<void(std::string_view, std::string_view, GenomeId)>& cb) const override
    {
        std::unordered_map<std::string, GenomeId> acc_to_gid;
        ar_.scan_genome_accessions([&](std::string_view acc, GenomeId gid) {
            acc_to_gid.emplace(std::string(acc), gid);
        });
        ar_.scan_taxonomy([&](std::string_view acc, std::string_view tax) {
            auto it = acc_to_gid.find(std::string(acc));
            cb(acc, tax, it != acc_to_gid.end() ? it->second : GenomeId{0});
        });
    }

    std::optional<GenomeMeta> genome_meta_by_accession(std::string_view acc) const override {
        return ar_.genome_meta_by_accession(acc);
    }
    const float* kmer_profile_by_accession(std::string_view acc) const override {
        return ar_.kmer_profile_by_accession(acc);
    }
    bool has_sketches() const override { return ar_.has_sketches(); }
    std::vector<uint32_t> available_kmer_sizes() const override {
        return ar_.available_sketch_kmer_sizes();
    }
    uint32_t sketch_sketch_size() const override { return ar_.sketch_sketch_size(); }
    const GstxEntry* gstx_for_genus(std::string_view genus) const override {
        return ar_.gstx_for_genus(genus);
    }
    const GstxReader* gstx_reader() const override {
        return ar_.gstx_reader();
    }
    const GcovEntry* gcov_for_genus(std::string_view genus) const override {
        return ar_.gcov_for_genus(genus);
    }
    const GcovReader* gcov_reader() const override {
        return ar_.gcov_reader();
    }
    const GcovEntry* fcov_for_family(std::string_view family) const override {
        return ar_.fcov_for_family(family);
    }
    const GcovReader* fcov_reader() const override {
        return ar_.fcov_reader();
    }
    FmhrView fmhr_for_genus(std::string_view genus) const override {
        return ar_.fmhr_for_genus(genus);
    }
    const FmhrReader* fmhr_reader() const override {
        return ar_.fmhr_reader();
    }
    bool has_core() const override { return ar_.has_core(); }
    bool has_fcore() const override { return ar_.has_fcore(); }
    CoreView core_for_family(std::string_view family) const override {
        return ar_.core_for_family(family);
    }
    const CoreReader* fcore_reader() const override { return ar_.fcore_reader(); }
    bool has_pcore() const override { return ar_.has_pcore(); }
    PcoreView pcore_for_genus(std::string_view genus) const override {
        return ar_.pcore_for_genus(genus);
    }
    PcoreView pcore_for_family(std::string_view family) const override {
        return ar_.pcore_for_family(family);
    }
    const PcoreReader* pcore_reader() const override { return ar_.pcore_reader(); }
    bool has_gami_v2() const override { return ar_.has_gami_v2(); }
    void load_gami_into(GlobalMultiplicityIndex& out) const override { ar_.load_gami_into(out); }
    CoreView core_for_genus(std::string_view genus) const override {
        return ar_.core_for_genus(genus);
    }
    const CoreReader* core_reader() const override { return ar_.core_reader(); }

    void visit_sketch_batches_multi_k(
        const std::vector<std::string>& accs,
        const std::vector<uint32_t>& ks, uint32_t sz,
        const std::function<void(size_t, uint32_t, const SketchResult&)>& cb) const override
    {
        std::vector<std::pair<GenomeId, size_t>> id_idx;
        id_idx.reserve(accs.size());
        for (size_t i = 0; i < accs.size(); ++i) {
            auto meta = ar_.genome_meta_by_accession(accs[i]);
            if (meta) id_idx.emplace_back(meta->genome_id, i);
        }
        std::sort(id_idx.begin(), id_idx.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        std::vector<GenomeId> sorted_ids;
        sorted_ids.reserve(id_idx.size());
        for (const auto& [gid, _] : id_idx) sorted_ids.push_back(gid);
        ar_.sketch_for_ids_multi_k(sorted_ids, sz,
            [&](size_t batch_idx, uint32_t ki, const SketchResult& sk) {
                if (ki >= ks.size()) return;
                cb(id_idx[batch_idx].second, ks[ki], sk);
            }, 1);
    }

    void visit_sketch_gids(
        const std::vector<GenomeId>& sorted_gids,
        uint32_t sz,
        const std::function<void(size_t, uint32_t, const SketchResult&)>& cb) const override
    {
        ar_.sketch_for_ids_multi_k(sorted_gids, sz, cb, 1);
    }

    void visit_shard_batches(
        const std::vector<std::string>& accs,
        const std::function<void(ArchiveReader::ShardBatch&)>& cb) const override
    {
        ar_.visit_shard_batches(accs, cb);
    }
    void visit_shard_batches_parallel(
        const std::vector<std::string>& accs,
        int n_readers,
        const std::function<void(ArchiveReader::ShardBatch&)>& cb) const override
    {
        ar_.visit_shard_batches_parallel(accs, n_readers, cb);
    }

private:
    ArchiveReader& ar_;
};

} // namespace genopack::check
