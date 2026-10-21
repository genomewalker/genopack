# Spec: stop inventing species and genus names; apply the self-named species overrides

Status as of 2026-10-21: proposed, not implemented, no Codex stream launched. Companion:
`$REPOS/geodesic/docs/SPEC-2026-10-21-unnamed-species-clustering.md` (stream B). Order: B merged before any
pack built with stream A is used by geodesic (see "Order").

## Why

GTDB-Tk leaves `s__` empty for genomes it cannot place in a species, and `g__` empty when it cannot place a genus.
genopack and the release scripts fill those gaps with invented names, and every consumer then treats them as real:

- `s__<accession>`: 500,846 of 1,313,026 reps in the r232 catalogue (38 %). geodesic groups by the full taxonomy string
  (`geodesic/src/pipeline.cpp:56-80`), so each became its own species and was never compared to the others. A de novo
  GTDB-rule clustering (2026-10-21, `/projects/caeg/scratch/kbd606/tmp/selfnamed_clust/s2/`) puts them in 173,057
  species.
- Family name copied into an empty genus: 20,286 genomes whose source gave no genus (GTDB-Tk empty `g__` for 15,626;
  GlobDB genus replaced for 4,660) carried their family name as genus; 22,417 of the 173,057 clusters have no member
  with a real genus assignment.

Sites that invent names (all to change):

| Site | What it does |
|---|---|
| `cli/main.cpp:2238-2283` `normalize_gtdb_10rank` (used by `taxonomy patch`, `:2424`) | `prop()` copies each empty rank from its parent (`:2264-2271`, incl. `g__` <- `f__`); empty `s__` -> `s__<stem>` (`:2273`); `S__` copies `s__` (`:2274`) |
| `cli/main.cpp:3936-3973` `normalize_gtdb` lambda (`taxonomy normalize`, called `:4010`) | same copy of the same logic; a bare-prefix domain becomes `d__Unclassified` with the stem at every rank (`:3959-3963`) |
| `src/ncbi_taxdb.cpp:391-428` `lineage_to_10rank` | `propagate()` copies down to `g__` (`:407-417`); empty `s__` -> `s__<genus> unclassified` or the accession at every rank (`:419-428`) |
| release repo `scripts/14_globdb_new_taxonomy.py:40-42` (`/maps/projects/caeg/people/kbd606/gtdb_r232_release_1.0`) | a rank GlobDB named after a genome takes its parent's name (`:41`); species = `s__<acc>;S__<acc>` (`:42`) |

## Stream A: genopack `feat/unresolved-ranks`

Repo `$REPOS/genopack`, branch `feat/unresolved-ranks` from `main`. The checkout is on
`feat/taxonomy-partition-contiguous` with uncommitted files: work in a worktree, do not touch that checkout.

A1. One normaliser. Replace both copies (`:2238` and the `:3936` lambda) with one function in `src/` (header in
`include/genopack/`), used by `taxonomy patch`, `taxonomy normalize` and the build path. `ncbi_taxdb.cpp`
`lineage_to_10rank` uses it for the rank-filling decision.

A2. Unresolved ranks stay empty. A rank the source left empty is written as the bare prefix (`g__`, `s__`, `S__`),
GTDB's own convention. No copy from the parent rank, no accession stem, no `<genus> unclassified`. The bare-prefix
domain case writes `d__;l__;...;S__` (all bare) instead of `d__Unclassified` plus stems. Higher ranks (`l__`, `k__`)
that GTDB never supplies keep their current derivation from `d__`, because they are a fixed mapping, not invention.

A3. Compatibility flag `--fill-unresolved` on `taxonomy patch`, `taxonomy normalize` and build reproduces today's output
byte for byte, for rebuilding old releases. Default is off (new behaviour).

A4. Consumers of TAXN / TXDB that assume every rank is named must accept bare prefixes. Audit and list each in the PR
with the change or "no change needed": TXDB build, `taxonomy unify` (`cli/main.cpp:2840` `cmd_taxonomy_unify`),
`taxonomy show/export`, `subset` taxonomy filters, the partition-by-taxonomy code, and the release repo's
`scripts/12_kaiku_taxdump.py`. A bare rank must not become a shared node named `""` that joins unrelated genomes: in the
unified tree an unresolved rank has no node and the genome hangs from its deepest named rank.

A5. Release repo `scripts/14_globdb_new_taxonomy.py`: a rank GlobDB named after a genome ID (`g__TPMCS_05909`) becomes
the bare prefix (`g__`), not the parent's name (today `:40-41`) and not the genome ID; write `s__;S__` (bare) instead
of `s__<acc>;S__<acc>` (`:42`). Keep the graft logic (`:28-35`). Keep the `ORDER` assert.

A gates (build dir, `ctest`, plus):
- unit test with GTDB-Tk strings: full 7 rank; empty `s__`; empty `g__` and `s__`; empty domain; already 10 rank.
  Expected outputs written in the test, for both default and `--fill-unresolved`;
- `--fill-unresolved` output byte-identical to current `main` on the r232 GTDB-Tk summaries
  (`/maps/projects/caeg/scratch/kbd606/geodesic/gtdb/232/gtdbtk_classify/gtdbtk_r232_merged.tsv`);
- default output on the same file: count of rows with bare `s__` equals the count with GTDB-Tk empty `s__`, same for
  `g__`; no row contains its own accession in any rank;
- `genopack verify` passes on a test pack built both ways; `taxonomy unify` on it creates no empty-named node.

## Applying `s4/s4_overrides.tsv` (step C, after A and B land; a script, not a Codex stream)

Input `/projects/caeg/scratch/kbd606/tmp/selfnamed_clust/s4/s4_overrides.tsv` (sha256 prefix 5ce716df, 464,051 rows,
S5 gates G1-G6 all pass, `s5/report.txt`):
- `s2_cluster` rows (463,911): `d__`..`g__` = the GTDB r232 lineage of `new_genus` (always a GTDB r232 genus, gate
  G6), `s__<species name>`, `S__` = same as `s__`. Taking `f__` and above from GTDB rather than the genome's current
  lineage keeps one parent per genus: S5 flags 220 `mixed_family` and 233 `cross_family_merge` genomes whose current
  family differs from their new genus's family. Report how many genomes change family or above.
- `s1_anchor` rows (140): the full 10-rank lineage of the target GTDB r232 species from
  `gtdbtk_db/taxonomy/gtdb_taxonomy.tsv` (normalised with A1).
- `genus_unresolved` clusters (22,417 clusters, 36,795 genomes) are not in the overrides. C1 writes a row for each
  of these genomes too: current `d__`..`f__` kept, bare `g__`, `s__`, `S__` (A2 changes only new normalisation, not
  existing TAXN, so without a row they would keep `s__<accession>`). Geodesic stream B clusters them at run time
  without naming them.

C1. Release repo `scripts/15_species_overrides_patch.py` writes the 2-column patch file (`accession TAB taxonomy`, the
"simple" format `taxonomy patch` already reads, `cli/main.cpp:2386-2427`) from the overrides and the current lineage
(`refdb/v1.0/taxonomy/common_8packs_2026-10-19_columnar`: `acc2taxid.tsv` + `taxonomy.tsv`, as `s4/s4.py` reads it;
genopack has no TAXN-to-TSV export, `genopack taxdump` writes taxdump/columnar). Asserts: every row 10 ranks; `d__`..`f__` unchanged for `s2_cluster` rows; each species
name maps to one genus.

C2. `genopack taxonomy patch` on each archive that holds these genomes (the 8 packs behind
`refdb/v1.0/taxonomy/common_8packs_*` and `reps_pack`), then `taxonomy unify` to a new dated
`common_8packs_<date>_columnar`. Do not overwrite the 2026-10-19 directory.

C gates: patched count = 464,051 + 36,795 = 500,846 per union of archives; no `s__<accession>` left among patched genomes; the 2 anchors
whose species has two nodes in the 2026-10-19 tree (`s__CAYMTR01 sp012719495`, `s__Methanocatella thaueri_B`) end
with one node each after unify, or the duplicate is listed as a known issue; S5 rerun on the new tree passes.

## Order

1. Stream B (geodesic) merged: it must handle bare `s__` before any pack has them, otherwise
   `group_by_taxonomy` would put every unnamed genome of a genus into one taxon and dereplicate across species.
2. Stream A merged; A5 in the release repo.
3. Step C on the existing packs. Rebuilding packs is not required: C patches TAXN in place.

## Not in scope

S6(a) rescoring and S6(b) re-derep (the user's call). Genus-level naming for the 22,417 unresolved clusters.
