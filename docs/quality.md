# Quality Scoring

`genopack check` computes per-genome completeness and contamination signals from
sketch geometry, aggregates them into a continuous quality score and an HQ/MQ/LQ
tier, and writes both a TSV and a 104-byte-per-genome `QUAL` section
(`include/genopack/qual.hpp:13-156`). All axes are derived from OPH containment,
conserved-k-mer (aamer) coverage, TNF geometry, and single-copy-core duplication —
no external tool is called at runtime. CheckM2 is the calibration yardstick, not a
dependency (see [CheckM2 alignment](#checkm2-alignment)).

The signals below are the columns of the `check` TSV
(`src/check/run_check.cpp:644-673`).

---

## Completeness

`completeness_effective` is the single completeness value fed to the score and
tier (`src/check/run_check.cpp:70-91`). It is an **intrinsic** estimate — the
fraction of the genome recovered, not the fraction of a genus pangenome covered.

### Intrinsic estimate priority

```
intrinsic = marker_completeness            if scored
          else completeness_aamer_core     if scored
          else completeness_post_decontam
```

| Estimator | Meaning | Notes |
|-----------|---------|-------|
| `marker_completeness` | present / expected single-copy genus markers (CheckM2-aligned, genus-calibrated) | primary signal; fraction-tracking, declines with fragmentation (`run_check.cpp:75-77`) |
| `completeness_aamer_core` | genus CORE conserved-k-mer recovery | fallback; presence-**saturating** (small dynamic range near 1.0), never NaN when scored (`run_check.cpp:78-80`, `qual.hpp:44`) |
| `completeness_post_decontam` | bp retained after the contig contamination scan | last resort |

`marker_completeness` is populated only on the pass-B marker route; QCOL-cache
genomes fall through to `aamer_core` (`run_check.cpp:78-79`).

### Cluster-relative corroborator

`completeness_cluster_relative` is the fraction of the genus **pangenome** (union
aamer content across genus members) a genome covers — accessory-genome breadth,
not intrinsic completeness (`run_check.cpp:62-65`, `657-662`). A finished isolate
in a diverse genus recovers ~100% of its own core but only a small slice of the
genus accessory pangenome, so a low `cluster_relative` there is genus diversity,
not missing sequence. It must never drive `completeness_effective` down on its own.

It is admitted only as a soft corroborator (`run_check.cpp:82-90`):

```
if isnan(intrinsic):                        return cluster_relative   # no intrinsic signal at all
if intrinsic < 0.50
   and cluster_relative < intrinsic
   and (intrinsic - cluster_relative) > 0.30:
                                            return sqrt(intrinsic * cluster_relative)
else:                                       return intrinsic
```

The geomean pulls the estimate down only when the intrinsic signal **already**
reads genuinely partial (below the MQ line of 0.50) and `cluster_relative` sits
well below it — genuine agreement that sequence is missing. A supported intrinsic
(≥ 0.50) keeps its value regardless of a low, diverse-genus `cluster_relative`.

---

## Contamination

### Aggregate determinant

`contamination_aggregate` is an NA-safe **max** over the trusted contamination
axes (`src/check/run_check.cpp:130-157`). The max means no single absent (NaN)
axis can lower the verdict, which is what a lone FMH signal used to do.

Axes that establish observability (NaN when unmeasured — a non-NaN value proves
contamination was actually observed):

- `fmh_minority_fraction` — FracMinHash minority mass
- `contamination_tnf_minor` — TNF-GMM minority mass (near-clade / same-genus merger detector)
- `contamination_contig_outlier_adj` — CCO, contig T²/SPE outlier bp fraction (requires GCOV)
- `contamination_spe` — SPE contig outlier fraction
- `contamination_rho_outlier` — correlation-outlier fraction
- `duplication_contamination(...)` — calibrated single-copy-core duplication (below)

Axes that raise the max but default to 0.0 when unmeasured (they cannot by
themselves establish observability):

- `contamination_leakage` — minimizer-mass leakage outside the expected genus range
- `contamination_tnf_excess` — TNF Mahalanobis distance from the genus centroid

`contamination_mixture` is **excluded** from gating: a BIC-selected K=2 minority
mass with no separation or multi-contig guard, uncorrelated with CheckM2
(`run_check.cpp:146-148`). It remains a reported TSV diagnostic only.

When no observability axis is available, contamination is unconfirmed and HQ is
capped at MQ (`run_check.cpp:733`).

### Duplication as a graded term

`contamination_core_dup_mass` is a non-saturating duplication-mass estimator
computed at build time over the single-copy-core (SCC) aamer set of the placed
(majority) genus (`src/cladesplit.cpp:600-610`):

```
core_dup_mass = Σ(count - 1) / Σ(count)      over SCC aamers present in the genome
```

where `count` is the copy number of each SCC aamer. It rises with the mixture
fraction instead of saturating near 1.0, and is NaN when the majority genus has no
SCC set (underpopulated / novel), which defers to the observability cap.

`duplication_contamination` maps it into CheckM2 contamination units by a
spike-panel OLS fit (`run_check.cpp:99-108`):

```
if core_dup_mass is scored:  c = 1.448954 * core_dup_mass + 0.009998   (clamped ≥ 0)
else (old-stride pack):      c = core_dup_excess * 0.125               (legacy ÷8 fallback)
```

The legacy `excess ÷ 8` term is saturating; the mass-based term replaces it when
available (constants `kDupMassSlope`/`kDupMassIntercept`/`kDupToContamScale`,
`run_check.cpp:101-103`).

### The ≥2-axis rule

A single firing axis only caps the tier at MQ; an LQ contamination demotion
requires at least two trusted axes to fire simultaneously, which kills single-axis
false positives (`trusted_contam_axes`, `src/check/run_check.cpp:163-172`). Per-axis
firing thresholds:

| Axis | Fires at |
|------|----------|
| `fmh_minority_fraction` | ≥ 0.10 |
| `contamination_tnf_minor` | ≥ 0.10 |
| `contamination_contig_outlier_adj` | ≥ 0.05 |
| `contamination_spe` | ≥ 0.10 |
| `contamination_rho_outlier` | ≥ 0.10 |
| `duplication_contamination(...)` | ≥ 0.05 |
| `contamination_leakage` | ≥ 0.02 |

`cross_genus` is handled separately as a hard chimera veto: a genome with
`contamination_cross_genus ≥ 0.10` (fraction of scored bp assigned to a foreign
genus) is demoted to LQ regardless of every other axis
(`cross_genus_chimera`, `run_check.cpp:178-180`).

### accessory_z

`accessory_z` is a relative-conspecific-containment z-score baked into GSTX at
build time (`src/check/pass_a.cpp:641-669`). Over the genus OPH-containment
distribution `c0_all` (only when the genus is `GenusSaturated` with ≥ 10 members
and MAD > 0; else NaN):

```
c0_query   = 1 - leakage[0]
c0_median  = median(c0_all)
c0_mad     = 1.4826 * median(|c0_all - c0_median|)
accessory_ratio = c0_query / c0_median
accessory_z     = (c0_query - c0_median) / c0_mad
```

The `1.4826` factor rescales the MAD to a standard-deviation-equivalent for a
normal reference. `accessory_z` is a reported diagnostic (not an input to the
score or tier).

---

## Quality score and tiers

### Continuous score

`compute_quality_score(comp_eff, contamination, fiedler)` replaces the discrete
cutoffs with a multiplicative score in [0,1] (`src/check/run_check.cpp:44-55`).
Multiplicative because completeness, contamination, and contig coherence are
near-independent failure modes — additive mixing would let high completeness mask
contamination.

```
C           = comp_eff        (0.5 neutral prior if NaN)
X           = contamination   (0 if NaN)
s           = fiedler         (0 if NaN)
cont_factor = 1 / (1 + (X / 0.10)^2)          # knee at 0.10: 5% → ≈0.80, 10% → 0.50
ramp        = clamp((s - 0.5) / 0.4, 0, 1)
coh_penalty = 1 - 0.5 * ramp                  # ≤ 50% cut, engaged only for s > 0.5
score       = C * cont_factor * coh_penalty
```

The contamination knee at 0.10 is a smooth analog of the old HQ/MQ cliffs, aligned
to CheckM2's < 5% HQ cutoff. The coherence penalty is bounded (≤ 50%) rather than a
hard gate because the fiedler threshold is uncalibrated (no test coverage,
`run_check.cpp:41-43`).

### Tier chain

The HQ/MQ/LQ tier is a rule chain over `comp_eff`, the aggregate `cont_val`, the
trusted-axis count, and the fiedler split (`src/check/run_check.cpp:721-733`).
Applied in order, later rules override earlier ones:

```
qtier = LQ
if comp_eff ≥ 0.90 and cont_val < 0.05:                     qtier = HQ
elif comp_eff ≥ 0.50 and cont_val < 0.10:                   qtier = MQ
elif comp_eff ≥ 0.50 and trusted_contam_axes < 2:           qtier = MQ   # single-axis rescue
if fiedler_oph_split < 0.10 and cont_signals ≥ 1:           qtier = LQ
if aamer_only and qtier == HQ:                              qtier = MQ   # saturating completeness only
if cross_genus_chimera and qtier != LQ:                     qtier = LQ   # hard veto
if qtier == HQ and contamination not observed:              qtier = MQ   # unconfirmed → cap at MQ
```

- `cont_val` is `contamination_aggregate` (the NA-safe max above).
- `trusted_contam_axes < 2` is the single-axis MQ rescue: one firing axis never
  demotes to LQ.
- `cont_signals` is a corroboration count (`fmh ≥ 0.10`, `contig_outlier ≥ 0.05`,
  `contig_split ≥ 0.05`, `leakage ≥ 0.02`) required before the uncalibrated
  fiedler-split rule may fire (`run_check.cpp:707-727`).
- `aamer_only` means completeness came from the saturating `aamer_core` fallback
  with no `post_decontam`/`cluster_relative` support — HQ is not granted on a
  saturating estimator alone.

### CheckM2 alignment

The axes are native sketch-geometry signals: OPH conspecific containment, aamer
core/pangenome coverage, TNF Mahalanobis/GMM geometry, and single-copy-core
duplication mass. CheckM2 is used only to calibrate and validate — the
`core_dup_mass → contamination` map is a spike-panel OLS fit to CheckM2 units
(`run_check.cpp:101-102`) and `marker_completeness` is described as CheckM2-aligned
(`run_check.cpp:76`) — but no CheckM/CheckM2 model is loaded or run by `genopack
check`. The score and tier thresholds (contamination knee 0.10, HQ < 0.05, MQ <
0.10) mirror CheckM2's HQ/MQ contamination cutoffs.
