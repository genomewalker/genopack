"""CLI regression: concurrent builds sharing the default spill dir must not clash.

Six builds with different genera run at once with the same TMPDIR (the default
PCORE spill dir), three rounds. Each must succeed, write the PCOR section it
writes when built alone, and leave no pspill files behind.

Usage: python3 concurrent_builds.py <genopack binary> <tmp parent dir>
"""
import os
import pathlib
import random
import struct
import subprocess
import sys
import tempfile

binary = sys.argv[1]
N = 6

def pcor_checksums(gpk):
    with open(gpk, 'rb') as h:
        h.seek(-128, 2)
        toc, size = struct.unpack_from('<QQ', h.read(128), 8)
        h.seek(toc)
        b = h.read(size)
    n = struct.unpack_from('<Q', b, 24)[0]
    return sorted(b[o + 64:o + 80].hex() for o in (128 + 112 * i for i in range(n)) if b[o:o + 4] == b'PCOR')

def build(tsv, out, tmpdir):
    return subprocess.Popen([binary, 'build', '-i', str(tsv), '-o', str(out)], stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, universal_newlines=True, env={**os.environ, 'TMPDIR': str(tmpdir)})

rng = random.Random(20261023)
with tempfile.TemporaryDirectory(prefix='genopack-concurrent-', dir=sys.argv[2]) as tmp:
    d = pathlib.Path(tmp)
    for j in range(N):
        rows = []
        for i in range(10):  # two genera of 5: PCORE needs >= 4 members (auto micro-genus threshold)
            fa = d / f'f{j}_g{i}.fa'
            fa.write_text(f'>f{j}g{i}\n' + ''.join(rng.choice('ACGT') for _ in range(12000)) + '\n')
            rows.append(f'F{j}G{i}\td__Bacteria;p__P;c__C;o__O;f__F{j};g__G{j}_{i % 2};s__G{j}_{i % 2} sp{i}\t{fa}\n')
        (d / f'in{j}.tsv').write_text('accession\ttaxonomy\tfile_path\n' + ''.join(rows))

    solo = d / 'solo_tmp'
    solo.mkdir()
    want = []
    for j in range(N):
        p = build(d / f'in{j}.tsv', d / f'solo{j}.gpk', solo)
        out = p.communicate()[0]
        assert p.returncode == 0, (j, out)
        want.append(pcor_checksums(d / f'solo{j}.gpk'))
        assert want[-1], f'fixture {j} has no PCOR section'

    shared = d / 'shared_tmp'
    shared.mkdir()
    for r in range(3):
        procs = [build(d / f'in{j}.tsv', d / f'r{r}_{j}.gpk', shared) for j in range(N)]
        outs = [p.communicate()[0] for p in procs]
        for j, (p, out) in enumerate(zip(procs, outs)):
            assert p.returncode == 0, (r, j, p.returncode, out[-600:])
            assert pcor_checksums(d / f'r{r}_{j}.gpk') == want[j], (r, j, 'PCOR differs from the solo build')
        left = [x.name for x in shared.iterdir() if x.name.startswith('pspill')]
        assert not left, (r, 'spill leftovers', left[:5])
print('concurrent_builds: PASS')
