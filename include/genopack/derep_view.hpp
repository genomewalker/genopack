#pragma once
#include "archive_set_reader.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <memory>
#include <span>
#include <string>
#include <string_view>

namespace genopack {

enum class DerepStaleness {
    Valid,
    LayoutChangedSameLiveSet,
    StaleNewGenomes,
    StaleTombstones,
    Mismatch,
};

struct RepStatus {
    enum class Kind {
        Representative,
        Member,
        Unclustered,
        UnknownSinceGeneration,
        Absent,
        Tombstoned,
    };
    Kind             kind            = Kind::Absent;
    uint32_t         rep_id          = UINT32_MAX;
    std::string_view rep_accession;
};

struct DerepStats {
    uint64_t n_genomes_indexed = 0;
    uint32_t n_reps            = 0;
    uint32_t n_unclustered     = 0;
    uint32_t n_singletons      = 0;
    uint64_t n_members         = 0;
};

class DerepView {
public:
    DerepView();
    ~DerepView();
    DerepView(const DerepView&)            = delete;
    DerepView& operator=(const DerepView&) = delete;
    DerepView(DerepView&&) noexcept;
    DerepView& operator=(DerepView&&) noexcept;

    void open(const std::filesystem::path& gpd_path);
    void close();
    bool is_open() const;

    DerepStaleness check(const ArchiveSetReader& archive) const;

    RepStatus status_for_accession(std::string_view accession) const;

    void scan_representatives(
        const std::function<void(uint32_t rep_id,
                                 std::string_view rep_accession,
                                 uint32_t cluster_size)>& cb) const;

    bool embedding_for_rep(uint32_t rep_id, std::span<float> out) const;

    DerepStats stats() const;
    uint32_t   embedding_dim() const;
    uint64_t   accession_set_hash() const;
    uint32_t   format_major() const;
    uint32_t   format_minor() const;

    std::array<uint8_t, 16> run_id() const;
    uint64_t                created_at_unix() const;
    uint16_t                source_n_parts() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace genopack
