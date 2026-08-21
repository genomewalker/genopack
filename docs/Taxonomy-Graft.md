# Unified Taxonomy

> **Status as of 2026-08-21.** Implemented and tested.

Combine independent reference packs (GTDB prokaryotes, organelles, viruses, and
optionally full NCBI) into **one taxonomy with NCBI taxids where they apply**, so
per-pack mapping references and a downstream classifier share a single taxid space.

Two entry points:

- `genopack build --ncbi-taxdump <dir>`: build a pack whose taxa carry native
  NCBI taxids where the name resolves, synthetic concept_ids otherwise.
- `genopack taxonomy unify <packs…> --ncbi-taxdump <dir> -o <out>`: reconcile the
  taxonomies of several **independent** packs into one taxdump **without merging
  their sequence data**.

## Taxid scheme

Taxids are `uint64_t` with bit 63 as a namespace flag:

- **bit 63 = 0:** native NCBI taxid (32-bit value).
- **bit 63 = 1:** synthetic concept_id, `FNV-64(canonical_path) | (1<<63)`.

The two ranges are disjoint, so NCBI-native and source-specific nodes coexist in
one tree with no renumbering and no collisions. Concept_ids are a deterministic
hash of the lineage path, so they are stable across NCBI releases. On export, the
columnar/taxdump writer keeps native NCBI taxids and renumbers only concept_ids
into a dense block above the maximum NCBI taxid.

## How taxids are assigned

For each lineage token, in order of precedence:

1. **Preserved:** if the accession has an authoritative NCBI taxid
   (`--preserve-taxids`), that taxid's real NCBI lineage is imported and the id is
   kept verbatim, down to a strain/subspecies leaf. Stale ids remap via
   `merged.dmp`.
2. **Name-resolved:** a taxon whose scientific name is **unambiguous** in the
   NCBI taxdump takes its native NCBI taxid. This gives NCBI precedence at every
   rank where names agree (e.g. `Bacteria`→2, `Pseudomonadota`→1224,
   `Escherichia coli`→562, `Viruses`→10239).
3. **Concept_id:** otherwise a synthetic hash id (GTDB-renamed clades, viral
   vOTUs, and other taxa NCBI does not name).

Ambiguous NCBI names are never used, so name resolution cannot cross homonyms.
Entries are processed in accession-sorted order, so the output is deterministic:
identical inputs produce identical taxids.

## `taxonomy unify`

```bash
genopack taxonomy unify prok.gpk organelle.gpk viral.gpk \
    --ncbi-taxdump ./taxdump \
    --preserve-taxids organelle_host_taxids.tsv \
    -o ./common_taxonomy -f columnar
```

Reads each pack's TAXN lineage strings (per-part reads run in parallel), builds
one tree, and writes the taxonomy. GTDB's propagated `l__`/`k__` placeholder ranks
(the domain name repeated at life/kingdom) are dropped so the prokaryote backbone
attaches phyla directly under the domain.

Outputs: `taxdump` format writes `nodes.dmp`, `names.dmp`, `acc2taxid.dmp`;
`columnar` writes `taxnodes.bin`, `acc2taxid.bin`, and `taxonomy.tsv`/`acc2taxid.tsv`.

### Preserving authoritative taxids

`--preserve-taxids` takes a TSV mapping accessions to NCBI taxids (2-column, or a
headered file with `accession` and `taxid`/`host_taxid` columns). Use it where a
pack has an authoritative id. For example, organelle genomes carry their host's
NCBI taxid. On the RefSeq organelle set this raises the fraction of accessions
matching their authoritative host taxid from 98.3% (name resolution) to 99.997%,
recovering strain/subspecies ids that a species name cannot express (e.g.
*S. cerevisiae* S288C 559292, *Canis lupus familiaris* 9615).

### Full NCBI backbone

By default the tree contains only taxa referenced by the packs. To retain taxa
with no reference genome for LCA, include the whole NCBI scaffold by adding
`--ncbi-backbone`. When both full NCBI and GTDB are present, add
`--gtdb-replaces-prok` to excise NCBI's own `Bacteria`(2)/`Archaea`(2157)
descendants (the anchors are kept) so GTDB stays authoritative for prokaryotes and
nothing is double-counted. NCBI prokaryotes that GTDB also has (name-matched) are
retained via the GTDB pack; NCBI-only prokaryotes are dropped.

## Where taxids come from, per source

| source | taxonomy in the pack | id in the unified tree |
|---|---|---|
| GTDB prokaryotes | GTDB lineage strings | native NCBI id where the name matches (domains, and many phyla/genera/species), concept_id where GTDB diverges |
| Organelles | host lineage / `--preserve-taxids` host taxid | the host's native NCBI taxid (species or strain) |
| Viruses (IMG/VR) | ICTV lineage strings under `d__Viruses` | native NCBI id for the ICTV backbone (`Viruses`→10239, kingdoms/phyla/…/families), concept_id for the vOTU leaves |
| Full NCBI (optional) | `--ncbi-backbone` | native NCBI ids for the whole scaffold |

## Guarantees

- **Single-rooted, acyclic.** Each node's parent is a strict prefix of its own
  lineage path, so parent pointers strictly decrease depth. Verified on the
  produced trees (root 1, no cycles, every node reaches the root).
- **Deterministic.** Accession-sorted processing makes concept_id assignment and
  name-resolution order-independent.

## Implementation

`src/txdb.cpp` (tree build + override), `src/ncbi_taxdb.cpp` (name index,
`canonical_path_taxids`, `merged.dmp`, backbone iteration, binary cache),
`src/builder.cpp`, `src/merger.cpp`, `cli/main.cpp` (`taxonomy unify`, taxdump
export). The NCBI taxdump parse is cached to `<dir>/.ncbi_taxdb.cache`
(staleness-checked; any read error falls back to a full parse).
