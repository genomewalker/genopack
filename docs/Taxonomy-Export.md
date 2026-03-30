# Taxonomy Export

`genopack taxdump` exports the taxonomy stored in a `.gpk` archive to two formats: NCBI-compatible taxdump and a high-performance columnar binary.

## Usage

```bash
genopack taxdump <archive.gpk> -f <format> -o <output_dir>
```

| Option | Default | Description |
|--------|---------|-------------|
| `-f, --format` | `columnar` | `taxdump` or `columnar` |
| `-o, --output` | required | Output directory (created if needed) |

## NCBI taxdump format (`-f taxdump`)

Produces files compatible with tools that consume the NCBI taxonomy dump (Kraken2, Kaiju, Bracken, taxonkit, etc.).

```
output_dir/
  names.dmp       — taxid | name | unique_name | name_class |
  nodes.dmp       — taxid | parent_taxid | rank | ... (NCBI pipe format)
  acc2taxid.dmp   — accession | accession.version | taxid | gi=0
  merged.dmp      — empty stub (required by some tools)
  delnodes.dmp    — empty stub (required by some tools)
```

**Rank mapping** — GTDB ranks are mapped to NCBI rank names:

| GTDB rank | NCBI name |
|-----------|-----------|
| domain | superkingdom |
| phylum | phylum |
| class | class |
| order | order |
| family | family |
| genus | genus |
| species | species |
| (intermediate) | no rank |

**Root node** — the root `taxid=1` is self-referential (`parent=1`) per NCBI convention.

**Taxids** — synthetic integer identifiers assigned by the TXDB builder. They are stable within an archive but are not NCBI taxids.

### Example: building a Kraken2 database from genopack taxonomy

```bash
genopack taxdump mydb.gpk -f taxdump -o ./kraken_taxonomy/
kraken2-build --download-taxonomy --taxonomy-dir ./kraken_taxonomy/
```

## Columnar binary format (`-f columnar`)

Designed for applications that need taxonomy without linking the genopack library — a single `mmap` of two small binary files is sufficient.

```
output_dir/
  acc2taxid.bin   — sorted (hash, taxid) binary array
  taxnodes.bin    — sorted node records + name pool
  acc2taxid.tsv   — accession → taxid (tab-separated, human-readable)
  taxonomy.tsv    — taxid, parent_taxid, rank, name, is_synthetic
```

### acc2taxid.bin

Maps accession strings → taxid via FNV-1a-64 hash, sorted ascending. Binary-search for O(log n) lookup.

**Layout:**

```
[Acc2TaxidHeader — 32 bytes]
[Acc2TaxidEntry × n_entries — sorted by acc_hash]
```

```c
struct Acc2TaxidHeader {
    uint32_t magic;       // 0x54415047 ("GPAT")
    uint32_t version;     // 1
    uint64_t n_entries;
    uint8_t  reserved[16];
};  // 32 bytes

struct Acc2TaxidEntry {
    uint64_t acc_hash;    // FNV-1a-64(accession_string)
    uint32_t taxid;
    uint32_t _pad;
};  // 16 bytes, sorted ascending by acc_hash
```

**Hash function** — same FNV-1a-64 used by CIDX, so lookup code is shared:

```c
uint64_t fnv1a64(const char* s, size_t len) {
    uint64_t h = 14695981039346656037ULL;
    for (size_t i = 0; i < len; i++) {
        h ^= (uint8_t)s[i];
        h *= 1099511628211ULL;
    }
    return h;
}
```

**Python lookup example:**

```python
import struct, mmap, os

def load_acc2taxid(path):
    fd = open(path, "rb")
    mm = mmap.mmap(fd.fileno(), 0, access=mmap.ACCESS_READ)
    magic, version, n_entries = struct.unpack_from("<IIQ", mm, 0)
    assert magic == 0x54415047
    return mm, n_entries

def fnv1a64(s: str) -> int:
    h = 14695981039346656037
    for c in s.encode():
        h ^= c
        h = (h * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return h

def lookup(mm, n_entries, accession: str) -> int | None:
    h = fnv1a64(accession)
    lo, hi = 0, n_entries
    base = 32  # sizeof(Acc2TaxidHeader)
    while lo < hi:
        mid = (lo + hi) // 2
        entry_hash = struct.unpack_from("<Q", mm, base + mid * 16)[0]
        if entry_hash < h:   lo = mid + 1
        elif entry_hash > h: hi = mid
        else:
            return struct.unpack_from("<I", mm, base + mid * 16 + 8)[0]
    return None
```

### taxnodes.bin

Stores all taxonomy nodes sorted by taxid. Each node has parent, rank, and a pointer into the name pool appended at the end.

**Layout:**

```
[TaxnodesHeader — 32 bytes]
[TaxnodeEntry × n_nodes — sorted by taxid]
[name_pool — null-terminated name strings]
```

```c
struct TaxnodesHeader {
    uint32_t magic;           // 0x4E545047 ("GPTN")
    uint32_t version;         // 1
    uint32_t n_nodes;
    uint32_t name_pool_size;  // bytes
    uint8_t  reserved[16];
};  // 32 bytes

struct TaxnodeEntry {
    uint32_t taxid;
    uint32_t parent_taxid;    // self-referential for root
    uint32_t name_offset;     // byte offset into name pool
    uint8_t  rank;            // TaxRank enum (0=no_rank … 9=species)
    uint8_t  flags;           // bit 0: synthetic node
    uint16_t name_len;        // bytes (excluding null terminator)
};  // 16 bytes
```

**TaxRank enum:**

| Value | Rank |
|-------|------|
| 0 | no_rank |
| 1 | domain |
| 4 | phylum |
| 5 | class |
| 6 | order |
| 7 | family |
| 8 | genus |
| 9 | species |

**Python lookup example:**

```python
import struct, mmap

def load_taxnodes(path):
    fd = open(path, "rb")
    mm = mmap.mmap(fd.fileno(), 0, access=mmap.ACCESS_READ)
    magic, version, n_nodes, name_pool_size = struct.unpack_from("<IIII", mm, 0)
    assert magic == 0x4E545047
    nodes_base = 32
    pool_base  = nodes_base + n_nodes * 16
    return mm, n_nodes, pool_base

def get_node(mm, n_nodes, pool_base, taxid: int):
    lo, hi = 0, n_nodes
    base = 32
    while lo < hi:
        mid = (lo + hi) // 2
        t = struct.unpack_from("<I", mm, base + mid * 16)[0]
        if t < taxid:   lo = mid + 1
        elif t > taxid: hi = mid
        else:
            t, p, name_off, rank, flags, name_len = struct.unpack_from("<IIIBBh", mm, base + mid * 16)
            name = mm[pool_base + name_off: pool_base + name_off + name_len].decode()
            return {"taxid": t, "parent": p, "rank": rank, "name": name, "synthetic": bool(flags & 1)}
    return None
```

### TSV sidecars

`acc2taxid.tsv` — tab-separated, one genome per line:
```
accession	taxid
GCA_000008085.1	450743105
GCA_000008885.1	1938317629
```

`taxonomy.tsv` — full node table:
```
taxid	parent_taxid	rank	name	is_synthetic
1	1	no_rank	root	0
450743105	81194823	species	Nanoarchaeum equitans	0
```

## Sizes (GTDB r226, 143,614 genomes, 181,960 nodes)

| File | Size |
|------|------|
| `acc2taxid.bin` | 2.2 MB |
| `taxnodes.bin` | 6.3 MB |
| `acc2taxid.tsv` | 3.7 MB |
| `taxonomy.tsv` | 9.0 MB |
| `names.dmp` (taxdump) | 9.8 MB |
| `nodes.dmp` (taxdump) | 13 MB |
| `acc2taxid.dmp` (taxdump) | 6.2 MB |
