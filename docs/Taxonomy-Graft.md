# Taxonomy Graft — `build --ncbi-taxdump`

> **Status as of 2026-08-20. IMPLEMENTED + tested.** Unblocks NCBI-taxid-preserving
> archives (euk/organelle/viral packs). Two entry points:
> - `genopack build -i in.tsv --ncbi-taxdump <dir> -o out.gpk` — build a pack whose
>   taxa keep native NCBI taxids where the scientific name is unambiguous in NCBI.
> - `genopack taxonomy unify a.gpk b.gpk … --ncbi-taxdump <dir> -o <outdir>` — **the
>   unify command**: reconcile the taxonomies of several INDEPENDENT packs into ONE
>   common taxdump (nodes/names/acc2taxid) **without merging sequence data**. Reads only
>   TAXN lineage strings, rebuilds one NCBI-reconciled tree, writes one taxid space
>   covering every accession — so per-pack bowtie2 references and kaiku's merge share
>   consistent taxids. Bacteria/Archaea/Viruses = NCBI 2/2157/10239 across all sources.
>   (`genopack merge … --ncbi-taxdump` also exists but merges the data too — use `unify`
>   when the packs must stay independent.)
>
> Mechanism: an override map (`build_ncbi_taxid_overrides`, `src/txdb.cpp`) resolves
> each lineage token against NCBI's unambiguous name index; `TxdbWriter::assign_taxid`
> uses the native NCBI taxid (bit63=0) when present, else a synthetic concept_id
> (bit63=1). The `taxdump` export preserves NCBI taxids and renumbers only synthetic
> ids into a block above `max_ncbi`. Verified end-to-end: Homo sapiens→9606,
> E. coli→562, Pseudomonadota→1224 (GTDB name, NCBI id), Caudoviricetes→2731619.
> Fixed two latent bugs: `split_dmp` dropped every name (record terminator) and the
> taxdump export flattened all taxids to 1..N.
>
> Files: `src/txdb.cpp`, `include/genopack/txdb.hpp`, `src/ncbi_taxdb.cpp`,
> `src/builder.cpp`, `src/merger.cpp`, `include/genopack/{archive,merger}.hpp`,
> `cli/main.cpp`.

## Goal

Build a single archive (or a mergeable set) whose TXDB is **one tree** where:

- the **eukaryote backbone carries real NCBI taxids** (so `9606` means *Homo
  sapiens* to anyone), and
- our **prokaryote (GTDB) and viral (MetaVR) taxa carry stable concept_ids**,
  grafted onto their real NCBI domain anchors, and
- **organelles inherit their host's NCBI taxid** for free.

Then `LCA` (`TaxonomyTree::lca`, `txdb.hpp:108`) returns a **real NCBI taxid**
on every euk / organelle / anchor node, and a concept_id only inside the
prok/viral subtrees.

## Why this needs new code

The build model today hashes **every** node to a concept_id. Verified:

- `TxdbWriter::add(accession, lineage_string)` is the only ingest — strings only
  (`src/txdb.cpp:342`).
- `assign_taxid(path)` always returns `FNV-64(path) | TAXID_GTDB_BIT`
  (`src/txdb.cpp:231-248`); the root (`taxid=1`) is the sole exception
  (`src/txdb.cpp:251-260`).
- `build_tree` therefore stamps `namespace=NCBI` only on the root
  (`src/txdb.cpp:379-380`).
- `build` has no `--ncbi-taxdump` and no taxid column (`cli/main.cpp:3072-3151`;
  taxonomy arrives only as the TSV `taxonomy` string, `src/builder.cpp:1768`).

The **data model already supports the goal** — `bit63=0` = native NCBI taxid,
`bit63=1` = concept_id, plus a per-node `TaxNamespace{NCBI,GTDB}`
(`txdb.hpp:20-24, 40-43`). The two id spaces are **disjoint by the high bit**,
so NCBI-native and synthetic nodes coexist with **no renumbering and no
collision**. Only the build path fails to populate the NCBI-native half.

## Inputs

```
genopack build \
  --input records.tsv \
  --ncbi-taxdump <dir>   # NCBI new_taxdump: nodes.dmp + names.dmp
  --output pack.gpk
```

`--ncbi-taxdump <dir>` is loaded once via `NcbiTaxdb::load(dir)`
(`ncbi_taxdb.hpp:34`; it already parses nodes.dmp/names.dmp and can auto-fetch).

Each input record declares **how its taxonomy resolves**, via two optional
columns added to the build TSV (both absent ⇒ current behaviour, all-hashed):

| column | meaning |
|---|---|
| `taxid` | an **NCBI taxid** for this genome (euk, organelle host, virus with an NCBI id). Triggers **NCBI-native** ingest. |
| `taxonomy` | a lineage string (GTDB `d__…`, ICTV `r__…`). Triggers **concept-id** ingest, grafted onto an NCBI anchor. |

A record may carry both: `taxid` places the **anchor**, `taxonomy` hangs a custom
subtree below it (used by organelles that want an explicit organelle child).

## Algorithm (one build pass)

Per record, classify and emit into the shared node set:

**A. NCBI-native** (`taxid` present):
1. Resolve the lineage root→`taxid` from `NcbiTaxdb.nodes_`
   (`ncbi_taxdb.hpp:63-68`: `taxid → {parent_taxid, rank}`).
2. Emit every node on that chain **once**, with its **native** taxid & parent
   (`bit63=0`), `namespace=NCBI`, rank mapped NCBI→`TaxRank`, name from
   `names.dmp`. Dedup across records by taxid.
3. The record's deepest taxid is its accession→taxid entry.

**B. Concept-id + graft** (`taxonomy` present, prok/viral):
1. Run the existing token loop (`src/txdb.cpp:279-296`) over the lineage —
   unchanged: `concept_id = FNV-64(path)|bit63`, parent = path-prefix tid.
2. **Seed the parent of the top custom token from a resolved NCBI anchor**
   instead of root `1`: map the lineage's domain/realm to an NCBI taxid
   (Bacteria→`2`, Archaea→`2157`, Viruses→`10239`) via a `--graft-anchor`
   table (name- or token-resolved through `NcbiTaxdb.taxid_for_name`,
   `ncbi_taxdb.hpp:57`). The anchor itself is emitted by path **A**.

The two halves share the root (`taxid=1`) and never collide (disjoint high bit).
Everything else — name pool, accession FNV hash, `finalize`, taxdump export —
is reused unchanged (`src/txdb.cpp:354-419`).

### Code touch list

- `TaxRank rank_from_ncbi(std::string_view)` — NCBI rank string → `TaxRank`
  (new; mirror of the GTDB prefix mapping).
- `TxdbWriter`: add `add_ncbi(accession, uint64_t ncbi_taxid)` and hold a
  `const NcbiTaxdb*` for backbone expansion; `build_internal` grows a
  namespace-aware branch (native taxid, no hashing) alongside the existing
  hashed branch (`src/txdb.cpp:231-338`).
- `assign_taxid`: unchanged for hashed nodes; native nodes bypass it.
- `build_tree`: stamp `namespace` from the node's origin, not from
  `taxid_is_gtdb` alone (`src/txdb.cpp:379-380`) — native nodes are NCBI even
  though the current heuristic already yields NCBI for `bit63=0`.
- `builder.cpp`: parse `taxid` column into a `taxid_map`; call `add_ncbi` vs
  `add` per record (`src/builder.cpp:2095-2125`).
- `cli/main.cpp`: `--ncbi-taxdump`, `--graft-anchor OLD=NCBI` flags on `build`
  (`cli/main.cpp:3072-3151`).

## Corner cases

1. **Anchors retained, descendants replaced** — keep NCBI `2`/`2157`/`10239`;
   NCBI's own bacterial/archaeal/viral descendants are simply not referenced by
   any NCBI-native record, so they never enter the tree. (No excision pass
   needed — unlike taxmerge, which streams a full nodes.dmp and must delete.)
2. **Rank mismatch** — NCBI `superkingdom` vs GTDB `domain`; ICTV `realm (r__)`
   has no `TaxRank`. Map `superkingdom→DOMAIN`; map ICTV `r__` to the anchor
   (Viruses) so the MetaVR lineage hangs below `10239` starting at its kingdom.
3. **Merged/deleted NCBI taxids** — resolve through `merged.dmp` when a `taxid`
   column value is stale (NcbiTaxdb already handles the dump).
4. **Duplicate accession** across packs — deepest-taxid wins per accession;
   prok/viral authoritative in their subtree.
5. **UINT32 / memory** — non-issues here: concept_ids are 64-bit by
   construction, native ids ≤ ~3.5M, and the tree is a sorted array with binary
   search (`txdb.hpp:76-88`), not an `O(max_taxid)` array.

## The three packs

All three are separate archives; merge or query-time union gives one LCA tree
(`genopack merge` already unions TAXN + rebuilds TXDB, `src/merger.cpp:485-496`).

| Pack | Input (verified) | Ingest |
|---|---|---|
| **prok** | canonical release `/maps/projects/caeg/people/kbd606/gtdb_r232_release_1.0/` — `reps_pack/part_0..6.gpk` (dereplicated-rep genopack) + `catalog/part_*/r232_sd*/*_derep_genomes.tsv` (self-describing); rep set `recommended_representatives.txt` (1,291,315). genopack v1.1.0, geodesic v1.2.0. | B — GTDB `d__…` hashed, graft `2/2157` |
| **organelle** | `…/taxDB/20220926/TaxDB_v6/results/organelle/{mitochondria,plastid}.fna.gz` + `organelle_genome_taxonomy.txt`; `genome_type ∈ {MITOCHONDRION,PLASTID}` | A — `taxid` = **host** NCBI taxid (RefSeq organelle acc → host taxid; e.g. `NC_012920→9606`). Optional organelle child = B below the host. |
| **viral** | `…/geodesic/metavirome/METAVR.fna.zst` + `metavr_acc_tax_votu.tsv.gz` (`accession, taxonomy(ICTV), votu`); parquet adds `host_taxonomy`, `genome_type`, checkv quality | B — ICTV hashed, graft `10239`. `host_taxonomy` available as an optional cross-anchor. |

Only **prok** is geodesic-dereplicated; organelle & viral are ingested as-is
(RefSeq organelle set; MetaVR vOTU representatives).

## Prep still required

- **prok**: filter the derep TSV to representatives; confirm `v16` is canonical.
- **organelle**: build the `accession → host NCBI taxid` map (from RefSeq
  organelle acc2taxid) so each record carries a `taxid` column.
- **viral**: choose `taxid` column source — none of the vOTUs are NCBI taxa, so
  viral is pure graft under `10239` (ICTV lineage → concept_id). Realm `r__`
  handling per corner case 2.
