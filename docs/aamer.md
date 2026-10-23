# AAMER: The Amino-Acid K-mer Subsystem

> Status as of 2026-10-23: CLI citations point at current `cli/main.cpp` lines. `pcore` bootstraps from `--theta` / `--min-seg-aa` / `--frac-max-hash` when an archive has no `SEC_CORE`; it no longer skips such archives.

AAMER is genopack's marker-gene-**free** completeness/contamination substrate. Its
k=8 amino-acid k-mer engine is shared byte-for-byte between the build side
(protein) and the check side (DNA). It backs three archive sections
(`SEC_CORE`, `SEC_PCORE`, `SEC_GAMI`), the clade-split contamination tier, the
SCG marker panel, and the per-contig foreign-containment scan. This page describes
the implementation across its source files.

## 1. Overview

An **aamer** is an amino-acid k-mer with **k=8** (`AAMER_K`,
`include/genopack/aamer.hpp:18`). Reference genomes almost never ship as protein;
query genomes never do. AAMER's build side ingests **protein** (GTDB-Tk MSA
sequences, `extract_aamers_aa`, `aamer.hpp:687-716`), and its check side ingests
**DNA**, 6-frame translated to protein in software
(`translate_6frame`/`extract_aamers_dna`, `aamer.hpp:357-476`, `620-630`). Both
paths hash the resulting amino-acid octamers with the **same** function
(`aamer_hash`, `aamer.hpp:109-124`), so a protein aamer built from an MSA column
and a DNA-translated aamer from a draft assembly produce the same 64-bit key.
Membership lookup therefore works across the protein↔DNA boundary without
materializing a shared FASTA format. A completeness model built from curated
GTDB-Tk protein alignments can score
(possibly fragmented, error-prone) nucleotide assemblies at check time, with no
external gene caller and no HMM search.

Because it is a raw k-mer signal rather than a curated single-copy-gene set,
AAMER continues to produce signal from partial ORFs and fragmented contigs
(`cladesplit.hpp:14-22`).

## 2. The AAMER engine

### 2.1 Encoding

`AA_ENC[256]` maps ASCII → a 20-letter alphabetical index (A=0 … Y=19, `0xFF`
sentinel for gaps/stop/ambiguous; handles upper/lower case, `-`, `*`, `X`)
(`aamer.hpp:22-46`). `DNA_ENC[256]` maps ASCII → `{A=0,C=1,G=2,T/U=3}`, `0xFF`
for anything else (`aamer.hpp:49-68`). Translation uses genetic code 11
(bacterial/archaeal), `CODON11[64]` (`aamer.hpp:79-92`).

### 2.2 Hashing

```cpp
uint64_t aamer_hash(const uint8_t* aa) noexcept;   // aamer.hpp:109
```

Packs 8 amino acids at 5 bits each (40 bits), then applies a splitmix64 mix
(`aamer.hpp:109-124`). Both sides must use the same AA encoding. They do because
`extract_aamers_aa` and `extract_aamers_dna`/`translate_6frame`
route through the same `AA_ENC`/`CODON11` tables. A Murphy-10 reduced-alphabet
variant (`aamer_hash_murphy10`, `aamer.hpp:128-143`, 10 classes × 4 bits) exists
for the marker panel's reduced-alphabet mode (§4, §6).

### 2.3 Low-complexity masking

```cpp
bool aamer_is_low_complexity(const uint8_t* aa, int k) noexcept;  // aamer.hpp:146
```

Rejects a k-mer if it contains a run of `>= (k+1)/2` identical amino acids. For
k=8, this is a run of ≥4 (integer division). The rule is applied at extraction time on
both build and check sides (`emit_aamers`, `aamer.hpp:156-162`).

### 2.4 6-frame translation

```cpp
template <typename Cb>
void translate_6frame(std::string_view seq, int min_aa_len, Cb&& cb);  // aamer.hpp:357
```

For each inter-stop amino-acid segment `>= min_aa_len` in all 6 reading frames
(0-2 forward, 3-5 reverse-complement), invokes
`cb(frame, aa_ptr, aa_len, nt_start, nt_end)`. Ambiguous bases (`N` etc.) break a
reading frame exactly like a stop codon. An AVX2 path (`aamer_detail::xlate8_dual`,
`aamer.hpp:307-342`) derives the forward and reverse-complement codons from one
shared base-encode/deinterleave pass, folding what would be 6 independent
translation passes into 3; a scalar fallback preserves identical output when
AVX2 is unavailable. `translate_3frame_fwd` (`aamer.hpp:544`) is a forward-only
variant used by marker scoring (§5), which never needs RC frames.

### 2.5 Extraction / FracMinHash

| Function | Side | Notes |
|---|---|---|
| `extract_aamers_aa(protein, k=8)` | build (protein) | MSA gaps (`-`) skipped without breaking the segment; stops break it (`aamer.hpp:687-716`) |
| `extract_aamers_dna(seq, k, min_seg_aa)` | check (DNA) | full 6-frame translate + hash, sorted-unique (`aamer.hpp:620-630`) |
| `extract_aamers_dna(seq, k, min_seg_aa, max_hash)` | check (DNA), subsampled | keeps only `h <= max_hash` using FracMinHash (`aamer.hpp:633-643`) |
| `extract_aamers_dna_into(...)` | check (DNA), in-place | appends into a caller-owned `thread_local` buffer, no per-call allocation (`aamer.hpp:648-653`) |

The FracMinHash rolling variant (`emit_aamers_frac`, `aamer.hpp:169-188`)
maintains the packed 40-bit value across the sliding window with 3 ops instead
of an 8-byte repack, and hash-gates *before* the low-complexity scan so the
expensive check only runs on the ~1/N survivors of the FMH subsample.

### 2.6 Dayhoff-6 variant (marker panel only)

A second, independent k-mer scheme lives in the same file: Dayhoff-6 reduces the
20-letter alphabet to 6 classes (`AA_DAYHOFF6[20]`, `aamer.hpp:721-723`), with
**k=12** (`AAMER_K_D6`) and closed-syncmer selection at inner s-mer length
**s=4** (`AAMER_S_D6`, `aamer.hpp:725-726`). A syncmer criterion selects a
k-mer only if its minimum s-mer sits at position 0, giving ~11% density
(`1/(k-s+1)`, `aamer.hpp:742-757`). This is a separate, coarser k-mer used only
by the `--dayhoff6` marker-panel build mode (§6); it is not part of the
CORE/PCORE/GAMI path, which is fixed at k=8 full-alphabet aamers.

## 3. The three archive sections

All three sections are keyed by a hashed taxon key (genus, or family where
noted) and store per-key aamer collections built at k=8 from the aamer engine
above.

| Section | Density | What it stores | Feeds |
|---|---|---|---|
| `SEC_CORE` (legacy) | sparse; aamers with member-prevalence `>= ceil(θ·N)` only (the "conserved core") | per-genus sorted-unique aamer set + contributing member ids + `core_model_hash` (`core_section.hpp:21-43`) | historically `completeness_aamer_core`; superseded by PCORE (see below) |
| `SEC_PCORE` | dense; **every** aamer seen in ≥1 member, plus a per-aamer member-prevalence | 3-run stratified encoding (singleton / multi / core) + quantized prevalence, optional per-aamer IDF tier byte (v2) (`pcore.hpp:10-16`, `32-89`) | `completeness_aamer_core` (via the θ-threshold core run) **and** small-foreign-contig detection (dense union) |
| `SEC_GAMI` | global; one `(hash, count)` pair per aamer across **all** genera, count capped at 255 | sorted pairs, zstd-compressed on disk, Elias-Fano + 2-bit packed counts at load (`gami.hpp:2-6`, `18-42`) | the specificity gate for foreign-aamer contamination scoring (rare-vs-universal) |

### 3.1 SEC_CORE: genus prevalence core

An entry is the set of aamers present in `>= ceil(theta * N)` of a genus's `N`
reference members (`CoreWriter::add_from_members`, `src/core_section.cpp:17-37`);
default `theta = 0.90` (`core_section.hpp:27` comment, matches `PcoreWriter`'s
default, `pcore.hpp:208`). The section is content-hashed
(`core_model_hash`, XXH3 over `(k, min_seg_aa, theta, sorted per-genus cores)`,
order-independent; `core_section.cpp:65-88`) so the intrinsic completeness
column is a pure function of build inputs, not the raw sibling set
(`core_section.hpp:6-10`).

**This section has no active builder in the current codebase.** `CoreWriter`'s
only call sites are its own defining file and its unit test
(`tests/core_section_test.cpp`); no CLI subcommand constructs one. `CoreReader`
remains available for reading legacy archives. In addition,
`genopack pcore` *sources its own build config* (`k`, `min_seg_aa`,
`theta`, `frac_max_hash`) from an existing archive's `SEC_CORE` header when one
is present. It uses the CoreReader *only* for those four config values, not the
core aamer content. `merge` also **drops** `SEC_CORE` (it is in the non-concatenable set,
`src/merger.cpp:190`), so a merged archive loses the ability to rebuild PCORE.

`pcore` therefore reads its build config from a pre-existing `SEC_CORE` when one
is present. When it is **absent**, `pcore` **bootstraps** from CLI options
`--theta` / `--min-seg-aa` / `--frac-max-hash` (`src/run_pcore.cpp:154-155` logs
`"has no SEC_CORE; bootstrapping build config"`) instead of skipping the archive,
so an archive can build `SEC_PCORE → completeness` without a legacy `SEC_CORE`.
The defaults are taken from the code (`k = AAMER_K = 8`;
`min_seg_aa = k = 8`, `aamer.hpp:621`; `theta = 0.90`, the `PcoreWriter`/unit-test
default; `frac_max_hash = UINT64_MAX`, i.e. dense/no-subsampling, required for the
1-2 kb small-contig detection provided by PCORE). These values affect quality;
override them only to match a specific legacy `SEC_CORE`.

### 3.2 SEC_PCORE: dense per-genus aamer reference

`SEC_PCORE` is the dense reference that supersedes CORE. It stores
**every** aamer seen in ≥1 member together with a per-aamer member prevalence,
from which a consumer derives, at query time (`pcore.hpp:6-9`):

- prevalence ≥ θ → the conserved completeness core (== old CORE, exact)
- prevalence ≥ 1 → the dense foreign/containment reference (small-contig detection)
- raw prevalence → the per-aamer weight for a log-likelihood-ratio contamination score

Two codecs exist. The legacy codec uses an open-addressing hash and a u8 member
count capped at 255. This known bug affects genera with >255 members, so the codec
is kept read-only for backward compatibility. V1 uses three contiguous PFOR-delta
runs per genus (singleton/multi/core), with
quantized log-normalized prevalence for multi+core; the core run is decided from
true u32 counts at build time so completeness stays exact regardless of
quantization (`pcore.hpp:10-16`, `32-36`). v2 adds a per-aamer 4-bit IDF tier
byte (`tier = floor(log2(genus_count))`, capped at 15,
`pcore_tier_from_count`, `pcore.hpp:110-117`), stamped in a separate
`genopack tier stamp` pass (§6) and used to compute a soft IDF weight
`w = 1 - tier/log2(G)` in place of the binary GAMI rare/common call
(`pcore_idf_weight`, `pcore.hpp:121-126`).

`PcoreWriter` builds bounded-memory: only ~64 bytes/genus of metadata stays
resident while per-genus encoded runs spill to a temp file
(`$GENOPACK_SPILL_DIR` or system temp) as they're added (`pcore.hpp:203-206`).

Coverage against the θ-threshold core is a plain sorted-merge intersection:

$$
\mathrm{completeness\_aamer\_core} = \frac{|\text{query aamers} \cap \text{core}|}{|\text{core}|}
$$

(`pcore_core_coverage`, `pcore.hpp:188-200`; re-implemented inline in
`src/check/pass_aamer.cpp:60-68` for the hot path).

### 3.3 SEC_GAMI: global aamer multiplicity index

The dense PCORE union tells a consumer *whether* an aamer is host- or
family-specific; it says nothing about how common that aamer is **outside**
the host's lineage. GAMI answers that: for every aamer seen anywhere in the
corpus, it stores the number of distinct genera that contain it, saturated at
255 (`GlobalMultiplicityIndex::increment`, `gami.hpp:305-310`). V1 was a
Bloom filter over rare aamers (removed); V2 (current) is exact sorted
`(hash, count8)` pairs, zstd-compressed on disk and loaded into an Elias-Fano
structure with 2-bit packed counts (4 states: absent/1/2/≥3) at query time
(`gami.hpp:18-42`). `genopack gami build` decodes PCORE (or CORE, as a
fallback) once and serializes this index, eliminating a 10-30 minute in-memory
GMI rebuild on every `genopack check` run (`gami.hpp:5-6`, `src/run_gami.cpp:69-141`).
If no `SEC_GAMI` section is present, `check` falls back to building the same
index in memory from PCORE/CORE at check-start (`src/check/pass_b.cpp:318-368`).

## 4. Protein k-mer primitives

The clade-split contamination tier (`cladesplit.hpp`, `.csp` panels) is a
second, independent consumer of the aamer engine that supports **three**
interchangeable protein k-mer primitives, selected by `ProteinKmer`
(`cladesplit.hpp:36-40`) and must match between panel build and score:

| Primitive | CLI name | Definition |
|---|---|---|
| `PK_AAMER` (default) | `aamer` | Plain k=8 aamer; reuses `aamer_hash` from the core engine directly (`src/cladesplit.cpp:95-96`) |
| `PK_METAMER` | `metamer` | Joint AA+DNA k-mer: 8 AA (5 bits each, 40 bits) **plus** the 8 third-codon-position bases (2 bits each, 16 bits) packed into 56 bits, then splitmix64; a Metabuli-style metamer (`cladesplit.cpp:97-101`) |
| `PK_STROBE` | `strobemer` | Order-2 protein randstrobe: a 4-AA "strobe" hash (`aa4_hash`) at position `i`, paired with a second 4-AA strobe chosen by minimizing the XOR link within a window `[i+9, i+min(i+24, len-4)]` (`s=4, wmin=5, wmax=20`), combined via splitmix64. Its gapped, longer span tolerates more mutations than a contiguous aamer (`cladesplit.cpp:66-90`) |

All three are FracMinHash-subsampled at `1/frac_c` via `h % c == 0`
(`cladesplit.cpp:88, 103`). This modulo-based subsampling differs from the
max-hash-threshold FracMinHash used in `extract_aamers_dna` (§2.5).
`CladeSplitConfig` defaults: `frac_c = 30`, `min_aa = 8`, `mode = PK_AAMER`
(`cladesplit.hpp:42-46`), matching the CLI defaults (`cli/main.cpp:5212-5218`).

A `.csp` panel maps each **genus-diagnostic** aamer (present in exactly one
genus among clean reference genomes) to that genus; a query genome's aamers
vote, and `minority_fraction` = votes to non-majority genus / total diagnostic
votes (`cladesplit.hpp:14-22`). v3 panels additionally carry a per-genus
**single-copy-core (SCC)** index: an aamer is SCC for a genus if it is
prevalent (present in `>= kSccPrevalence = 0.90` of the genus's genomes,
`cladesplit.hpp:33`) *and* single-copy (occurs exactly once in
`>= 1 - kSccMultiFrac = 0.90` of the genomes where present,
`cladesplit.hpp:34`); genera with fewer than `kSccMinRefs = 10` reference
genomes get no SCC set (`cladesplit.hpp:32`). `core_dup` = fraction of SCC
aamers occurring ≥2× in a query, feeding `contamination_core_dup_mass` /
`duplication_contamination` in the quality score (see
[quality.md, Duplication as a graded term](quality.md#duplication-as-a-graded-term)).

## 5. How AAMER feeds quality

Two independent aamer-derived signals reach `genopack check`
(full model: [quality.md](quality.md)):

**Completeness:** `completeness_aamer_core` is the intrinsic-completeness
**fallback** behind `marker_completeness` in the priority chain
(`completeness_effective`, [quality.md § Completeness](quality.md#completeness)).
`src/check/pass_aamer.cpp` (`run_pass_aamer`) computes it only for accessions
that pass-A left un-scored (`std::isnan(completeness_cluster_relative)`,
`pass_aamer.cpp:20-24`): each candidate's contigs are 6-frame translated and
FracMinHash-subsampled with the *same* `k`/`min_seg_aa`/`frac_max_hash` used to
build the PCORE reader. These values must match because PCORE stores them in its
header (`pass_aamer.cpp:88`). The result is intersected against that genome's genus completeness-core
(`pcore_for_genus(...).core_into(...)`, `pass_aamer.cpp:43-44`).
`marker_completeness` (§ below) is also aamer-engine output, but is scored
against a curated single-copy-gene panel instead of a raw genus core.

**Contamination (small-contig channel):** PCORE's dense union plus GAMI's
global multiplicity together drive `score_contig_foreign_indexed`
(`src/check/foreign_contam.hpp:254-336`), a composition-independent per-contig
signal for contigs down to ~1-2 kb. TNF/GCOV composition needs
5-20 kb of context (`foreign_contam.hpp:5-6`). Each contig's aamers are
classified host-specific (in the host genus PCORE union), foreign-specific (in
some other genus's union, not host-family, and carried by at most
`K` genera according to GAMI, so it is taxon-specific rather than universal), or excluded
(family-conserved or too common). `K` defaults to
`max(2, n_genera/500)` (capped 255), overridable via `GENOPACK_FOREIGN_K`
(`src/check/pass_b.cpp:338-344`). Below 16 classifiable aamers the scorer
abstains, or (if a family reference exists) falls back to a family-hit-density
proxy (`foreign_contam.hpp:314-334`). Unlike `completeness_aamer_core`, this
signal is **per-contig**, not a genome-level D/S/G channel: it's persisted to
`SEC_QCONTIG` (`prot_foreign_frac`, `prot_loglr`, `prot_classifiable`,
`prot_flags`; `include/genopack/qcontig.hpp:44-48`) and surfaced by
`genopack qcontig` / `genopack report` to flag *which* contigs to drop
(the `"foreign"` channel, alongside `"tnf"` and `"leakage"`,
`src/run_report.cpp:508-546`). It is a decontamination signal, not one of the
genome-level D/S/G contamination channels documented in
[quality.md § Contamination](quality.md#contamination).

**SCG markers:** `src/check/marker_score.hpp::score_markers_for_genome` also
runs on the raw aamer engine: each contig is forward-3-frame-translated
(`translate_6frame`, `MRK_MIN_ORF_AA = 50`, `marker_score.hpp:49, 62`), and every
ORF's k=8 aamers (`aamer_hash` or, for reduced-alphabet panels,
`aamer_hash_murphy10`, selected by `mrk_rd.is_murphy10()`) are FracMinHash-gated
against the marker panel's `frac_max` and matched to a marker id via a bloom
prefilter + per-genus flat hit-map, or a binary search fallback
(`marker_score.hpp:65-90`). A contig votes for a marker if any ORF clears
`max(min_hits, orf_len*3%)` hits (`MRK_MIN_FRAC_PC = 3`,
`marker_score.hpp:50, 92-93`); `completeness = present/expected`,
`redundancy = markers with ≥2 votes / expected` (`marker_score.hpp:104-118`).

## 6. Building it

Pipeline order (see [cli.md](cli.md) for the full flag reference; this section
lists only the AAMER-relevant flags):

1. **`genopack pcore <archive>`:** builds `SEC_PCORE` from an archive that
   already carries a `SEC_CORE` header (§3.1). Streams every live genome's
   6-frame-translated, FracMinHash-subsampled aamers, grouped by genus **and**
   family, into `PcoreWriter` (`src/run_pcore.cpp:137-193`).
   `-t/--threads` (default 8, parallel shard readers); `--members <file>`
   restricts which accessions contribute (default: all live genomes)
   (`cli/main.cpp:5110-5134`).

2. **`genopack gami build <archive>`:** builds `SEC_GAMI` v2 from the archive's
   PCORE (or CORE fallback). `-t/--threads` (default 8, decode parallelism)
   (`cli/main.cpp:5137-5151`).

3. **`genopack tier merge` / `genopack tier stamp`** (optional): merges
   per-part `.ptier` side-channel files (emitted during PCORE build via
   `PcoreWriter::set_tier_sidechannel`, `pcore.hpp:219-222`) into a global
   `.gtier` IDF table, then stamps per-aamer tier bytes into PCORE, upgrading
   v1 → v2 in place (`cli/main.cpp:5152-5184`).

4. **`genopack markers build --gtdbtk-db <root> -o markers.mrk`**: builds the
   SCG marker aamer panel from GTDB-Tk MSA files (bac120 + ar53).
   `--threshold` (default 0.30, presence-call threshold); `--scale` (default 1,
   FracMinHash divisor, ignored with `--dayhoff6`); `--dayhoff6` (off by
   default; switches to the k=12 Dayhoff-6 syncmer profile, §2.6);
   `--ic-threshold` (default 0.25, `--dayhoff6` only); `--expected-min-frac`
   (default 0.50; a marker counts as expected for a genus only if detectable
   in at least this fraction of the genus's reference genomes, mirroring
   CheckM2's ~97% single-copy universality criterion) (`cli/main.cpp:5383-5421`).

5. **`genopack cladesplit build -i genomes.tsv -o panel.csp`** (independent of
   1-4): builds a `.csp` contamination panel from clean reference genomes.
   `--frac-c` (default 30), `--min-aa` (default 8), `--primitive`
   (`aamer`\|`metamer`\|`strobemer`, default `aamer`) (§4;
   `cli/main.cpp:5211-5225`). `genopack cladesplit score` and
   `genopack cladesplit aamers` consume/dump the same primitive and must be
   invoked with matching flags.
