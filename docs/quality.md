# Quality Scoring

`genopack check` computes per-genome completeness and contamination signals from
sketch geometry, aggregates them into a continuous quality score and an HQ/MQ/LQ
tier, and writes both a TSV and a 56-byte-per-genome `QUAL` section
(`include/genopack/qual.hpp`). All axes are derived from OPH containment,
aamer coverage, TNF geometry, and single-copy-core duplication — no external tool
is called at runtime. CheckM2 is the calibration yardstick, not a dependency (see
[CheckM2 alignment](#checkm2-alignment)).

The signals below are the columns of the `check` TSV
(`src/check/run_check.cpp:644-673`).

## Terms

| Term | Definition |
|------|------------|
| **aamer** | An amino-acid **k=8** k-mer. Query DNA is 6-frame translated (`translate_6frame`); inter-stop segments of ≥ 8 AA yield the k-mers. Build-side aamers come from protein input (GTDB-Tk MSA). A GUNC-style protein-k-mer signal, **not** a nucleotide k-mer (`include/genopack/aamer.hpp:13-18`, `cladesplit.hpp:15-19`). |
| **OPH** | One-permutation hash — the MinHash sketch stored in SKCH (`include/genopack/skch.hpp`). |
| **SCC** | Single-copy core — aamers prevalent across a genus and single-copy where present (`cladesplit.hpp:27-30`). |
| **TNF** | Tetranucleotide frequency — the 136-dim k=4 composition profile. |
| **CCO** | Contig contamination outlier — per-contig T²/SPE exceedance fraction (requires GCOV). |
| **SPE** | Squared prediction error — the PCA-residual half of the T²/SPE outlier pair. |
| **FMH** (FracMinHash) | Fractional MinHash — the minority-fraction contamination axis. |
| **GSTX** | Per-genus sketch-statistics section: OPH consensus, p90 completeness, containment dispersion (median/MAD). |

---

## Completeness

`completeness_effective` is the single completeness value fed to the score and
tier (`src/check/run_check.cpp:55-76`). It is an **intrinsic** estimate — the
fraction of the genome recovered, not the fraction of a genus pangenome covered.

### Intrinsic estimate priority

$$
\mathrm{intrinsic} =
\begin{cases}
\mathrm{marker\_completeness} & \text{if scored} \\
\mathrm{completeness\_aamer\_core} & \text{else if scored} \\
\mathrm{completeness\_post\_decontam} & \text{otherwise}
\end{cases}
$$

| Estimator | Meaning | Notes |
|-----------|---------|-------|
| `marker_completeness` (QualRecord field; TSV column `completeness_marker`) | present / expected single-copy genus markers (CheckM2-aligned, genus-calibrated) | primary signal; fraction-tracking, declines with fragmentation. Field: `types.hpp:90`; TSV header: `run_check.cpp:673` |
| `completeness_aamer_core` | fraction of the genus prevalence-core aamer set (amino-acid 8-mers) recovered | fallback; presence-**saturating** (small dynamic range near 1.0), never NaN when scored (`run_check.cpp:78-80`, `qual.hpp:44`) |
| `completeness_post_decontam` | bp retained after the contig contamination scan | last resort |

`marker_completeness` is scored on the pass-B FASTA route. When `--markers` is
supplied, **every** genome is routed through pass-B (`pass_b.cpp`,
`score_all_completeness = scan_all || markers given`), so `marker_completeness` is
populated catalog-wide (≈ 99% of GTDB r232) rather than only for pass-A-flagged
genomes; a genome that still lacks it (no marker panel, or served from a stale
QCOL cache) falls through to `aamer_core`.

### Cluster-relative corroborator

`completeness_cluster_relative` is the fraction of the genus **pangenome** (union
aamer content across genus members) a genome covers — accessory-genome breadth,
not intrinsic completeness (`run_check.cpp:62-65`, `657-662`). A finished isolate
in a diverse genus recovers ~100% of its own core but only a small slice of the
genus accessory pangenome, so a low `cluster_relative` there is genus diversity,
not missing sequence. It must never drive `completeness_effective` down on its own.

It is admitted only as a soft corroborator (`run_check.cpp:82-90`):

$$
\mathrm{comp\_eff} =
\begin{cases}
\mathrm{cr} & \text{if intrinsic is NaN} \\[2pt]
\sqrt{\mathrm{intrinsic}\cdot \mathrm{cr}} & \text{if } \mathrm{intrinsic} < 0.50 \ \wedge\ \mathrm{cr} < \mathrm{intrinsic} \ \wedge\ (\mathrm{intrinsic}-\mathrm{cr}) > 0.30 \\[2pt]
\mathrm{intrinsic} & \text{otherwise}
\end{cases}
$$

where $\mathrm{cr} = $ `completeness_cluster_relative`.

The geomean pulls the estimate down only when the intrinsic signal **already**
reads genuinely partial (below the MQ line of 0.50) and `cluster_relative` sits
well below it — genuine agreement that sequence is missing. A supported intrinsic
(≥ 0.50) keeps its value regardless of a low, diverse-genus `cluster_relative`.

---

## Contamination

### Three reported channels (D / S / G)

Contamination is **reported, never a discard gate** — a genome is LQ on
completeness alone (see [Tier chain](#tier-chain)). `contamination_channels`
emits three near-independent channels (`src/check/run_check.cpp:114-163`),
replacing the old NA-safe **max** over seven axes. Four of those axes (the
geometry axes below) co-fire at 20–50× independence — they are one TNF/GCOV
signal — so a max counted that signal up to four times; calibrated on all 9.53 M
genomes with CheckM2 as truth, the resulting gate scored 6.9% PPV / 5.7% recall,
demoting 362,935 genomes to catch ~1,200 real contaminants.

- **D** — `duplication_contamination(...)`, the **only** CheckM2-calibrated
  channel (single-copy-core duplication mass, below).
- **S** — `fmh_minority_fraction`, the FracMinHash k-mer minority mass.
- **G** — the **median** of the present geometry axes
  (`contamination_rho_outlier`, `contamination_spe`,
  `contamination_contig_outlier_adj`, `contamination_tnf_minor`), collapsing the
  one correlated signal to a single vote.

`leakage` is dropped (mathematically dead: max 0.005 ≪ its 0.02 threshold) and
`tnf_excess` is dropped (default-0, untrusted). `contamination_mixture` (a
BIC-selected K=2 minority mass, uncorrelated with CheckM2) and `cross_genus` are
reported TSV diagnostics only, never gates. A `display` value — the noisy-OR
union of the present channels — is emitted for the TSV/human view only;
dereplication consumes the raw D/S/G channels as an explicit reliability-ordered
**D → S → G** tiebreak (`run_check.cpp:104-113`).

### Duplication as a graded term

`contamination_core_dup_mass` is a non-saturating duplication-mass estimator
computed at build time over the single-copy-core (SCC) aamer set of the placed
(majority) genus (`src/cladesplit.cpp:600-610`):

$$
\mathrm{core\_dup\_mass} = \frac{\sum_i (c_i - 1)}{\sum_i c_i}
$$

where $c_i$ is the copy number of SCC aamer $i$ present in the genome (sum over
present SCC aamers). It rises with the mixture
fraction instead of saturating near 1.0, and is NaN when the majority genus has no
SCC set (underpopulated / novel), which defers to the observability cap.

`duplication_contamination` maps it into CheckM2 contamination units by a
spike-panel OLS fit (`run_check.cpp:87-96`):

$$
c_{\mathrm{dup}} =
\begin{cases}
\max\!\left(0,\ 1.448954 \cdot \mathrm{core\_dup\_mass} + 0.009998\right) & \text{if core\_dup\_mass scored} \\[2pt]
0.125 \cdot \mathrm{core\_dup\_excess} & \text{otherwise (legacy, saturating)}
\end{cases}
$$

The legacy `excess ÷ 8` term is saturating; the mass-based term replaces it when
available (constants `kDupMassSlope`/`kDupMassIntercept`/`kDupToContamScale`,
`run_check.cpp:87-89`).

### Channels do not gate the tier

Contamination is decoupled from the LQ tier entirely. The old gate demoted
362,935 genomes at 6.9% PPV; removing it eliminates every one of those wrongful
discards and leaves only the genuine-incompleteness floor (`comp_eff < 0.50`).
Contamination retains **one** narrow role: the single CheckM2-calibrated channel,
duplication **D**, may cap HQ → MQ (`D ≥ 0.05`), so the top tier still means
"clean and complete". This cap never forces LQ.

`contam_channels_fired` (`src/check/run_check.cpp:166-172`) counts channels over
threshold and is reported for the TSV/derep view — it is **not** a tier input.
`cross_genus` is a **ranker, not a gate**: on a 2,036-genome GTDB-Tk stratified
sample its positive predictive value as a demotion rule is 2.8% at ≥ 0.10
(`run_check.cpp:199-213`), so it never demotes. Per-channel report thresholds:

| Channel / axis | Reports at |
|------|----------|
| `duplication_contamination(...)` (D) | ≥ 0.05 |
| `fmh_minority_fraction` (S) | ≥ 0.10 |
| `contamination_contig_outlier_adj` (G) | ≥ 0.05 |
| `contamination_rho_outlier` / `contamination_spe` / `contamination_tnf_minor` (G) | ≥ 0.10 |

### accessory_z

`accessory_z` is a relative-conspecific-containment z-score baked into GSTX at
build time (`src/check/pass_a.cpp:641-669`). Over the genus OPH-containment
distribution `c0_all` (only when the genus is `GenusSaturated` with ≥ 10 members
and MAD > 0; else NaN):

$$
c_{0,\mathrm{query}} = 1 - \mathrm{leakage}[0], \qquad
\mathrm{MAD} = \operatorname{median}\big(\lvert c_{0,\mathrm{all}} - \operatorname{median}(c_{0,\mathrm{all}})\rvert\big)
$$

$$
\mathrm{accessory\_ratio} = \frac{c_{0,\mathrm{query}}}{\operatorname{median}(c_{0,\mathrm{all}})}, \qquad
\mathrm{accessory\_z} = \frac{c_{0,\mathrm{query}} - \operatorname{median}(c_{0,\mathrm{all}})}{1.4826 \cdot \mathrm{MAD}}
$$

The $1.4826$ factor rescales the MAD to a standard-deviation-equivalent for a
normal reference. `accessory_z` is a reported diagnostic (not an input to the
score or tier).

---

## Quality score and tiers

### Continuous score

`compute_quality_score(comp_eff)` is a threshold-free score in [0,1] equal to
`comp_eff` itself (`src/check/run_check.cpp:38-40`; NaN → 0.5 neutral prior). It is
**completeness-only** — contamination no longer folds into the score. The old
multiplicative `cont_factor` was a back-door gate: a 0.10 knee halved the score on
a signal whose PPV vs CheckM2 never exceeds ~15%, silently reordering the catalog
on noise. Contamination is reported separately as the D/S/G channels above and
consumed by dereplication as an explicit D → S → G tiebreak, so the ranking stays
interpretable and the noisy channels can never override the calibrated one.

### Tier chain

The HQ/MQ/LQ tier is **completeness-only** (`src/check/run_check.cpp:181-197`),
shared byte-for-byte between the QCOL (`quality_tier_u8`) and TSV paths so they
cannot drift:

```
ce = completeness_effective(q)
if ce is NaN or ce < 0.50:            LQ
elif ce < 0.90:                       MQ
else:  # ce ≥ 0.90
    if aamer_only:                    MQ   # saturating completeness alone can't support HQ
    elif D ≥ 0.05:                    MQ   # calibrated duplication caps HQ→MQ
    else:                             HQ
```

- Contamination **never** demotes to LQ; the LQ floor is genuine incompleteness
  (`ce < 0.50`) only. The old contamination gate (max-over-axes, cross_genus
  veto, ≥2-axis rule) is gone.
- `aamer_only` means completeness came from the saturating `aamer_core` fallback
  with no `post_decontam`/`cluster_relative` support — HQ is not granted on a
  saturating estimator alone.
- `D` is `duplication_contamination(...)`, the only CheckM2-calibrated channel; it
  is the sole contamination signal that touches the tier, and only to distinguish
  HQ from MQ.

### CheckM2 alignment

The axes are native sketch-geometry signals: OPH conspecific containment, aamer
core/pangenome coverage, TNF Mahalanobis/GMM geometry, and single-copy-core
duplication mass. CheckM2 is used only to calibrate and validate — the
`core_dup_mass → contamination` map is a spike-panel OLS fit to CheckM2 units
(`run_check.cpp:87-96`) and `marker_completeness` is genus-calibrated to be
CheckM2-aligned — but no CheckM/CheckM2 model is loaded or run by `genopack
check`. The one place a contamination signal meets a CheckM2-derived cutoff is the
duplication cap `D ≥ 0.05` on HQ (`run_check.cpp:190-192`), mirroring CheckM2's
5% HQ contamination line; the tier is otherwise completeness-only.
