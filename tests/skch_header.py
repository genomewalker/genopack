"""Read SKCH section headers straight from a .gpk file (TailLocator -> TOC -> SkchSeekHdr)."""
import pathlib
import struct

SKC4 = struct.pack('<I', 0x34434B53)

def skch_params(pack):
    """(kmer_sizes, syncmer_s) of every SKCH section in the live TOC."""
    b = pathlib.Path(pack).read_bytes()
    toc_off = struct.unpack_from('<Q', b, len(b) - 128 + 8)[0]
    n = struct.unpack_from('<Q', b, toc_off + 24)[0]
    out = []
    for i in range(n):
        off = struct.unpack_from('<Q', b, toc_off + 128 + 112 * i + 16)[0]
        if b[off:off + 4] == SKC4:
            nk = struct.unpack_from('<I', b, off + 20)[0]
            out.append((list(struct.unpack_from(f'<{nk}I', b, off + 24)), struct.unpack_from('<I', b, off + 56)[0]))
    return out
