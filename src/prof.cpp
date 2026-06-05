#include <genopack/prof.hpp>
#include <genopack/section_checksum.hpp>
#include <genopack/toc.hpp>
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <xxhash.h>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

namespace genopack {

// ── policy hash ───────────────────────────────────────────────────────────────
uint64_t profile_policy_hash(const Profile& p) {
    std::vector<const ProfileSelection*> sels;
    sels.reserve(p.selections.size());
    for (const auto& s : p.selections) sels.push_back(&s);
    std::sort(sels.begin(), sels.end(),
              [](const ProfileSelection* a, const ProfileSelection* b) { return a->axis < b->axis; });

    std::string buf;
    buf.reserve(p.name.size() + sels.size() * 32 + 16);
    buf.append(p.name);
    buf.push_back('\0');
    auto put = [&](uint64_t v) {
        for (int i = 0; i < 8; ++i) buf.push_back(static_cast<char>((v >> (8 * i)) & 0xFF));
    };
    for (const auto* s : sels) {
        buf.append(s->axis);
        buf.push_back('\0');
        put(s->column_identity);
        put(s->fallback_identity);
    }
    return XXH3_64bits(buf.data(), buf.size());
}

// ── ProfWriter ────────────────────────────────────────────────────────────────
uint32_t ProfWriter::intern(const std::string& s) {
    if (s.empty()) return 0;
    auto off = static_cast<uint32_t>(strtab_.size());
    strtab_.append(s);
    strtab_.push_back('\0');
    return off;
}

SectionDesc ProfWriter::finalize(AppendWriter& w, uint64_t section_id) {
    const auto n_profiles = static_cast<uint32_t>(profiles_.size());
    uint32_t n_selections = 0;
    for (const auto& p : profiles_) n_selections += static_cast<uint32_t>(p.selections.size());

    // Build descriptor + selection arrays first (they intern strings into strtab_).
    std::vector<ProfileDesc>       descs(n_profiles);
    std::vector<ProfSelectionDesc> sels;
    sels.reserve(n_selections);
    for (uint32_t i = 0; i < n_profiles; ++i) {
        const Profile& p = profiles_[i];
        ProfileDesc& d   = descs[i];
        d.policy_hash = profile_policy_hash(p);
        d.name_off    = intern(p.name);
        d.desc_off    = intern(p.description);
        d.sel_start   = static_cast<uint32_t>(sels.size());
        d.sel_count   = static_cast<uint32_t>(p.selections.size());
        for (const auto& s : p.selections) {
            ProfSelectionDesc sd{};
            sd.axis_off          = intern(s.axis);
            sd.sel_flags         = 0;
            sd.column_identity   = s.column_identity;
            sd.fallback_identity = s.fallback_identity;
            sels.push_back(sd);
        }
    }

    w.align(8);
    const uint64_t section_start = w.current_offset();

    ProfHeader hdr{};
    w.append(&hdr, sizeof(hdr));

    if (n_profiles) w.append(descs.data(), static_cast<uint64_t>(n_profiles) * sizeof(ProfileDesc));
    if (n_selections) w.append(sels.data(), static_cast<uint64_t>(n_selections) * sizeof(ProfSelectionDesc));

    const uint64_t strtab_off = w.current_offset() - section_start;
    w.append(strtab_.data(), strtab_.size());

    const uint64_t section_end = w.current_offset();

    hdr.magic          = PROF_CONTAINER_MAGIC;
    hdr.format_version = PROF_FORMAT_VERSION;
    hdr.flags          = 0;
    hdr.n_profiles     = n_profiles;
    hdr.n_selections   = n_selections;
    hdr.strtab_off     = strtab_off;
    hdr.strtab_len     = strtab_.size();
    w.write_at(section_start, &hdr, sizeof(hdr));

    SectionDesc sd{};
    sd.type              = SEC_PROF;
    sd.version           = PROF_FORMAT_VERSION;
    sd.flags             = 0;
    sd.section_id        = section_id;
    sd.file_offset       = section_start;
    sd.compressed_size   = section_end - section_start;
    sd.uncompressed_size = section_end - section_start;
    sd.item_count        = n_profiles;
    sd.aux0              = n_selections;
    sd.header_size       = sizeof(ProfHeader);
    return sd;
}

// ── ProfReader ────────────────────────────────────────────────────────────────
void ProfReader::open(const uint8_t* data, uint64_t offset, uint64_t size) {
    if (size < sizeof(ProfHeader)) throw std::runtime_error("PROF section too small");
    data_   = data + offset;
    header_ = reinterpret_cast<const ProfHeader*>(data_);
    if (header_->magic != PROF_CONTAINER_MAGIC) throw std::runtime_error("PROF: bad magic");
    if (header_->format_version > PROF_FORMAT_VERSION)
        throw std::runtime_error("PROF: format version too new");

    const uint64_t n_profiles   = header_->n_profiles;
    const uint64_t n_selections = header_->n_selections;
    const uint64_t descs_off    = sizeof(ProfHeader);
    const uint64_t sels_off     = descs_off + n_profiles * sizeof(ProfileDesc);
    const uint64_t descs_end    = sels_off;
    const uint64_t sels_end     = sels_off + n_selections * sizeof(ProfSelectionDesc);
    const uint64_t st_end       = header_->strtab_off + header_->strtab_len;
    if (descs_end > size || sels_end > size || st_end > size)
        throw std::runtime_error("PROF section truncated");

    descs_      = reinterpret_cast<const ProfileDesc*>(data_ + descs_off);
    sels_       = reinterpret_cast<const ProfSelectionDesc*>(data_ + sels_off);
    strtab_     = reinterpret_cast<const char*>(data_ + header_->strtab_off);
    strtab_len_ = header_->strtab_len;

    for (uint64_t i = 0; i < n_profiles; ++i) {
        const ProfileDesc& d = descs_[i];
        if (static_cast<uint64_t>(d.sel_start) + d.sel_count > n_selections)
            throw std::runtime_error("PROF: profile selections out of bounds");
    }
}

std::string_view ProfReader::str_at(uint32_t off) const {
    if (off == 0 || off >= strtab_len_) return {};
    const char* p = strtab_ + off;
    size_t maxlen = static_cast<size_t>(strtab_len_ - off);
    return {p, ::strnlen(p, maxlen)};
}

Profile ProfReader::profile(uint32_t i) const {
    const ProfileDesc& d = descs_[i];
    Profile p;
    p.name        = std::string(str_at(d.name_off));
    p.description = std::string(str_at(d.desc_off));
    p.selections.reserve(d.sel_count);
    for (uint32_t s = 0; s < d.sel_count; ++s) {
        const ProfSelectionDesc& sd = sels_[d.sel_start + s];
        ProfileSelection sel;
        sel.axis              = std::string(str_at(sd.axis_off));
        sel.column_identity   = sd.column_identity;
        sel.fallback_identity = sd.fallback_identity;
        p.selections.push_back(std::move(sel));
    }
    return p;
}

uint64_t ProfReader::policy_hash(uint32_t i) const { return descs_[i].policy_hash; }

int ProfReader::find(std::string_view name) const {
    for (uint32_t i = 0; i < n_profiles(); ++i)
        if (str_at(descs_[i].name_off) == name) return static_cast<int>(i);
    return -1;
}

// ── archive append ────────────────────────────────────────────────────────────
void write_prof_to_archive(const std::filesystem::path& gpk,
                           const std::vector<Profile>& profiles) {
    int lock_fd = ::open(gpk.c_str(), O_RDWR);
    if (lock_fd < 0) throw std::runtime_error("profile: cannot open for locking: " + gpk.string());
    struct LockGuard { int fd; explicit LockGuard(int f):fd(f){::flock(fd,LOCK_EX);}
                       ~LockGuard(){::flock(fd,LOCK_UN);::close(fd);} } lock(lock_fd);

    MmapFileReader mmap;
    mmap.open(gpk);
    auto toc = TocReader::read(mmap);
    const auto* tail = mmap.ptr_at<TailLocator>(mmap.size() - sizeof(TailLocator));
    uint64_t prev_toc_offset = tail->toc_offset;
    uint64_t generation      = toc.header.generation + 1;
    mmap.close();

    AppendWriter writer;
    writer.open_append(gpk);
    uint64_t section_id = toc.next_section_id();

    ProfWriter pw;
    for (const auto& p : profiles) pw.add(p);
    SectionDesc prof_sd = pw.finalize(writer, section_id);

    TocWriter new_toc;
    for (const auto& sd : toc.sections)
        if (sd.type != SEC_PROF) new_toc.add_section(sd);
    new_toc.add_section(prof_sd);

    writer.flush();
    stamp_section_checksums(gpk.c_str(), new_toc.sections(), /*only_if_zero=*/true);

    new_toc.finalize(writer, generation,
                     toc.header.live_genome_count, toc.header.total_genome_count,
                     prev_toc_offset, toc.header.catalog_root_section_id,
                     toc.header.accession_root_section_id, toc.header.tombstone_root_section_id);
}

} // namespace genopack
