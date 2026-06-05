#pragma once
// ── SEC_PROF — named reporting/fusion profiles ────────────────────────────────
// A profile is a NAMED, content-hashed resolution policy: for each quality axis
// (completeness, contamination, …) it pins the exact column identity to report
// (plus an optional fallback identity). It makes fusion EXPLICIT and auditable —
// a reported "completeness" is whatever column the active profile selected, and
// the report carries that column's provenance (tool/method). Built-in profiles
// (intrinsic / external / calibrated / best) are resolved by rule at report time
// and need no storage; SEC_PROF persists USER-authored profiles so a site can pin
// and audit an exact resolution that survives across runs.
//
// On-disk layout (within the section, all offsets from section start):
//   [ProfHeader 32B]
//   [ProfileDesc[n_profiles]]        (immediately after header)
//   [ProfSelectionDesc[n_selections]]
//   [string table]                   (NUL-delimited; offset 0 == "")
#include "format.hpp"
#include "mmap_file.hpp"
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace genopack {

static constexpr uint32_t PROF_CONTAINER_MAGIC = 0x31465250u; // "PRF1"
static constexpr uint16_t PROF_FORMAT_VERSION  = 1;

struct ProfHeader {              // 32 bytes
    uint32_t magic;              //  0  PROF_CONTAINER_MAGIC
    uint16_t format_version;     //  4  PROF_FORMAT_VERSION
    uint16_t flags;              //  6
    uint32_t n_profiles;         //  8
    uint32_t n_selections;       // 12  total across all profiles
    uint64_t strtab_off;         // 16
    uint64_t strtab_len;         // 24
};
static_assert(sizeof(ProfHeader) == 32);

struct ProfileDesc {             // 32 bytes
    uint64_t policy_hash;        //  0  content hash of (name + sorted selections)
    uint32_t name_off;           //  8  strtab (>=1)
    uint32_t desc_off;           // 12  strtab (0 == none)
    uint32_t sel_start;          // 16  index into the selection array
    uint32_t sel_count;          // 20
    uint64_t _reserved;          // 24
};
static_assert(sizeof(ProfileDesc) == 32);

struct ProfSelectionDesc {       // 24 bytes
    uint32_t axis_off;           //  0  strtab (axis name, >=1)
    uint32_t sel_flags;          //  4
    uint64_t column_identity;    //  8  chosen column's identity hash
    uint64_t fallback_identity;  // 16  0 == none
};
static_assert(sizeof(ProfSelectionDesc) == 24);

// ── In-memory model ───────────────────────────────────────────────────────────
struct ProfileSelection {
    std::string axis;                 // qual_axis::*
    uint64_t    column_identity   = 0;  // chosen column's column_identity_hash
    uint64_t    fallback_identity = 0;  // 0 == none
};

struct Profile {
    std::string                    name;
    std::string                    description;
    std::vector<ProfileSelection>  selections;
};

// Deterministic content hash of a profile's policy: name + selections sorted by
// axis (axis NUL-delimited, then column/fallback identity in little-endian).
// `description` is cosmetic and excluded, exactly as `label` is for columns.
uint64_t profile_policy_hash(const Profile& p);

// ── Writer ────────────────────────────────────────────────────────────────────
class ProfWriter {
public:
    void add(Profile p) { profiles_.push_back(std::move(p)); }
    SectionDesc finalize(AppendWriter& w, uint64_t section_id);
    size_t n_profiles() const { return profiles_.size(); }

private:
    uint32_t intern(const std::string& s);
    std::vector<Profile> profiles_;
    std::string          strtab_{std::string(1, '\0')};  // offset 0 == ""
};

// ── Reader ────────────────────────────────────────────────────────────────────
class ProfReader {
public:
    void open(const uint8_t* data, uint64_t offset, uint64_t size);
    bool is_open() const { return data_ != nullptr; }

    uint32_t n_profiles() const { return header_ ? header_->n_profiles : 0; }
    // Materialize the i-th profile (resolves all strings into owned std::string).
    Profile  profile(uint32_t i) const;
    uint64_t policy_hash(uint32_t i) const;
    // Find a profile by exact name; -1 if absent.
    int      find(std::string_view name) const;

private:
    std::string_view str_at(uint32_t off) const;

    const uint8_t*           data_       = nullptr;
    const ProfHeader*        header_     = nullptr;
    const ProfileDesc*       descs_      = nullptr;
    const ProfSelectionDesc* sels_       = nullptr;
    const char*              strtab_     = nullptr;
    uint64_t                 strtab_len_ = 0;
};

// Append a fresh SEC_PROF (dropping any prior one) to a single .gpk archive,
// mirroring the XQAL/CQAL archive-append idiom (flock → read TOC → append section
// → rewrite TOC dropping the old SEC_PROF → stamp checksums → finalize TOC).
void write_prof_to_archive(const std::filesystem::path& gpk,
                           const std::vector<Profile>& profiles);

} // namespace genopack
