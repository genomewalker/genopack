"""CLI regression: reindex --skch resolves auto syncmer s from the k it actually sketches, including an inherited k.

Usage: python3 reindex_syncmer_inherit.py <genopack binary> <tmp parent dir>
"""
import pathlib
import random
import shutil
import subprocess
import sys
import tempfile

from skch_header import skch_params

binary = sys.argv[1]

def run(*args):
    p = subprocess.run([binary, *map(str, args)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    assert p.returncode == 0, p.stdout + p.stderr
    return p.stdout + p.stderr

rng = random.Random(20261023)
with tempfile.TemporaryDirectory(prefix='genopack-syncmer-', dir=sys.argv[2]) as tmp:
    d = pathlib.Path(tmp)
    rows = []
    for i in range(4):
        (d / f'g{i}.fa').write_text(f'>g{i}\n' + ''.join(rng.choice('ACGT') for _ in range(8000)) + '\n')
        # One shared species: reindex's GSTX step indexes [0] of an empty list when every species is a singleton.
        rows.append(f'G{i}\td__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G sp\t{d}/g{i}.fa\n')
    (d / 'in.tsv').write_text('accession\ttaxonomy\tfile_path\n' + ''.join(rows))
    base = d / 'base.gpk'
    run('build', '-i', d / 'in.tsv', '-o', base, '--sketch-kmers', '', '--sketch-kmer', 21, '--sketch-syncmer', 3)
    assert skch_params(base) == [([21], 3)], ('fixture: expected one k=21 s=3 SKCH', skch_params(base))

    cases = [  # (reindex args, expected SKCH after)
        ((), ([21], 7)),                                   # k inherited (21) -> auto s = 21/3, s=3 not inherited
        (('--sketch-kmer', 24), ([24], 8)),                # explicit k
        (('--sketch-syncmer', 0), ([21], 0)),              # explicit s=0 disables
        (('--sketch-kmers', '15,21'), ([15, 21], 5)),      # multi-k: smallest k
    ]
    for i, (args, want) in enumerate(cases):
        pack = d / f'case{i}.gpk'
        shutil.copyfile(base, pack)
        run('reindex', pack, '--skch', *args)
        got = skch_params(pack)
        assert got == [want], (args, got, want)

print('reindex_syncmer_inherit: PASS')
