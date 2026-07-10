#pragma once
#include "types.hpp"
#include "pass_a.hpp"
#include "pack_iface.hpp"
#include <string>
#include <unordered_map>

namespace genopack::check {

// Lightweight aamer-core completeness pass for fragmented genomes (completeness_cluster_relative=NaN).
// Computes completeness_aamer_core (genus prevalence-core coverage) from FASTA + PCORE only;
// does not load GMI, GCOV, mixture model, or Fiedler.  Run after pass-A, before pass-B.
void run_pass_aamer(ICheckReader& pack,
                    const PassAResult& pass_a,
                    std::unordered_map<std::string, GenomeQuality>& quality,
                    int threads);

} // namespace genopack::check
