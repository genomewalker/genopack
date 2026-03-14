#include <genopack/archive.hpp>
#include <spdlog/spdlog.h>
#include <fstream>
#include <stdexcept>

namespace genopack {

// ── ArchiveAppender (stub for v1) ─────────────────────────────────────────────
// Full implementation: reads existing manifest, writes new shard, updates catalog,
// increments generation. For now: minimal tombstone support.

struct ArchiveAppender::Impl {
    std::filesystem::path archive_dir;
    std::vector<std::string> tombstones;
    std::vector<BuildRecord> pending;

    explicit Impl(const std::filesystem::path& dir) : archive_dir(dir) {
        if (!std::filesystem::exists(dir / "MANIFEST.bin"))
            throw std::runtime_error("Not a genopack archive: " + dir.string());
    }
};

ArchiveAppender::ArchiveAppender(const std::filesystem::path& dir)
    : impl_(std::make_unique<Impl>(dir))
{}
ArchiveAppender::~ArchiveAppender() = default;

void ArchiveAppender::add_from_tsv(const std::filesystem::path& tsv_path) {
    spdlog::warn("ArchiveAppender::add_from_tsv: not yet implemented");
}
void ArchiveAppender::add(const BuildRecord& rec) {
    impl_->pending.push_back(rec);
}
void ArchiveAppender::remove(GenomeId id) {
    impl_->tombstones.push_back(std::to_string(id));
}
void ArchiveAppender::remove_by_accession(std::string_view accession) {
    impl_->tombstones.emplace_back(accession);
    spdlog::info("Tombstoned: {}", accession);
}
void ArchiveAppender::commit() {
    if (!impl_->tombstones.empty()) {
        // Write tombstones to a simple text file for now
        auto tpath = impl_->archive_dir / "tombstones.txt";
        std::ofstream f(tpath, std::ios::app);
        for (const auto& t : impl_->tombstones)
            f << t << "\n";
        spdlog::info("Wrote {} tombstones", impl_->tombstones.size());
    }
    if (!impl_->pending.empty())
        spdlog::warn("ArchiveAppender::commit: {} pending records not written (not implemented)",
                     impl_->pending.size());
}

} // namespace genopack
