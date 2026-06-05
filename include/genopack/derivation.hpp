#pragma once
#include "format.hpp"
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <vector>
#include <xxhash.h>

namespace genopack {

// Content-addressed derivation (§5). Every non-SOURCE section stores a
// derivation_hash that folds its params, model versions, semantic schema, and
// the CONTENT hashes of its upstream inputs — so any upstream content / param
// change provably invalidates stale derived bytes (the --from-gpk reuse oracle).

enum class SectionClass { Source, Index, Derived };

// §5.1 classification.
inline SectionClass section_class(uint32_t type) {
    switch (type) {
        case SEC_SHRD: case SEC_CATL:
        case SEC_NTDB: case SEC_TXDB: case SEC_GTAX:
        case SEC_XQAL:                       // externally ingested quality — ground truth
        case SEC_PROF:                       // user-authored reporting/fusion policy
            return SectionClass::Source;     // ground truth / policy, reusable verbatim
        case SEC_ACCX: case SEC_TAXN: case SEC_GIDX: case SEC_CIDX:
            return SectionClass::Index;      // reusable iff rowset + deletions + schema match
        default:
            return SectionClass::Derived;    // SKCH, FMHR, KMRX, GSTX, GCOV, FCOV, QUAL, BPRM…
    }
}

// §5.2 upstream content dependencies (which section TYPES feed each section).
inline std::vector<uint32_t> upstream_types(uint32_t type) {
    switch (type) {
        case SEC_SKCH: return {SEC_SHRD};
        case SEC_FMHR: return {SEC_SHRD};
        case SEC_KMRX: return {SEC_SHRD};
        // CORE per-genus prevalence cores: built from raw sequence (SHRD) keyed
        // by taxonomy (CATL). The intrinsic columnar quality (QCOL) folds CORE.
        case SEC_CORE: return {SEC_SHRD, SEC_CATL};
        // FCORE per-family prevalence cores: same inputs as CORE (raw sequence keyed
        // by taxonomy), aggregated at family rank instead of genus.
        case SEC_FCORE: return {SEC_SHRD, SEC_CATL};
        // PCORE unified dense per-genus aamer reference: same inputs as CORE.
        case SEC_PCORE: return {SEC_SHRD, SEC_CATL};
        // GSTX/GCOV/FCOV bucket and key by taxonomy (CATL) and are computed from
        // raw sequence (SHRD), not from SKCH alone — fold both so a taxonomy
        // relabel with identical sketches/params still changes derivation_hash.
        case SEC_GSTX: return {SEC_SKCH, SEC_SHRD, SEC_CATL};
        case SEC_GCOV: return {SEC_SKCH, SEC_SHRD, SEC_CATL};
        case SEC_FCOV: return {SEC_SKCH, SEC_SHRD, SEC_CATL};
        case SEC_QUAL: return {SEC_FMHR, SEC_GSTX, SEC_GCOV, SEC_SHRD, SEC_CATL};
        // QCOL = intrinsic quality, columnar; supersedes QUAL and additionally folds
        // the CORE model it references (completeness_aamer_core column).
        case SEC_QCOL: return {SEC_FMHR, SEC_GSTX, SEC_GCOV, SEC_CORE, SEC_FCORE, SEC_SHRD, SEC_CATL};
        // XQAL = externally ingested quality (CheckM2/anvi'o): a leaf input with no
        // section-content upstream (classified Source → derivation_hash 0).
        case SEC_XQAL: return {};
        // CQAL = calibrated/fused quality: an explicit function of the intrinsic
        // (QCOL) and external (XQAL) columns plus the calibration model.
        case SEC_CQAL: return {SEC_QCOL, SEC_XQAL, SEC_CATL};
        // PROF = named reporting policy: references columns by content-hash identity,
        // not by section content (classified Source → derivation_hash 0).
        case SEC_PROF: return {};
        // QCONTIG = per-contig quality overlay: same intrinsic inputs as QCOL plus the
        // CORE/FCORE prevalence cores (the per-contig foreign-aamer containment channel).
        case SEC_QCONTIG: return {SEC_FMHR, SEC_GSTX, SEC_GCOV, SEC_CORE, SEC_FCORE, SEC_PCORE, SEC_SHRD, SEC_CATL};
        case SEC_ACCX: case SEC_TAXN:
        case SEC_GIDX: case SEC_CIDX: return {SEC_CATL};
        default: return {};
    }
}

// §4 semantic_schema_hash: identity of ordering/hash/collation for a section.
// v1 folds {type, per-section schema version}. The section `version` is bumped
// whenever sort order / key normalization / hash seed / comparator / partition
// scheme changes, so this value changes in lockstep — O(1) equality/reject.
inline uint64_t semantic_schema_hash_of(uint32_t type, uint16_t version) {
    uint64_t buf[2] = { type, version };
    return XXH3_64bits(buf, sizeof(buf));
}

// Populate semantic_schema_hash + derivation_hash for every section in the
// directory. MUST run after content_xxh128 (checksum) is computed for all of
// them. SOURCE sections get derivation_hash = 0 (always reusable verbatim).
// Each non-SOURCE section folds: type, its semantic_schema_hash,
// build_params_hash, deletion_set_hash (INDEX only), and the SORTED
// content_xxh128 + semantic_schema_hash of every upstream section.
//
// Caveat: a section left unchecksummed by the size cap contributes zero bytes,
// degrading its dependents from content-bound to param-bound for that input.
inline void populate_derivation(std::vector<SectionDesc>& secs,
                                uint64_t build_params_hash,
                                uint64_t deletion_set_hash) {
    for (auto& s : secs)
        s.semantic_schema_hash = semantic_schema_hash_of(s.type, s.version);

    for (auto& s : secs) {
        if (section_class(s.type) == SectionClass::Source) {
            s.derivation_hash = 0;
            continue;
        }
        // Assemble the derivation preimage as a contiguous u64 buffer, hashed
        // once: [type, semantic_schema_hash, build_params_hash, (deletion_set_hash
        // for INDEX), sorted(upstream lo, hi, semantic_schema_hash)...].
        std::vector<uint64_t> pre;
        pre.push_back(static_cast<uint64_t>(s.type));
        pre.push_back(s.semantic_schema_hash);
        pre.push_back(build_params_hash);
        if (section_class(s.type) == SectionClass::Index)
            pre.push_back(deletion_set_hash);

        std::vector<uint64_t> up;
        for (uint32_t ut : upstream_types(s.type))
            for (const auto& o : secs)
                if (o.type == ut) {
                    uint64_t lo, hi;
                    std::memcpy(&lo, o.checksum,     8);
                    std::memcpy(&hi, o.checksum + 8, 8);
                    up.push_back(lo);
                    up.push_back(hi);
                    up.push_back(o.semantic_schema_hash);
                }
        std::sort(up.begin(), up.end());
        pre.insert(pre.end(), up.begin(), up.end());

        s.derivation_hash = XXH3_64bits(pre.data(), pre.size() * sizeof(uint64_t));
    }
}

} // namespace genopack
