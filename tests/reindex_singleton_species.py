"""CLI regression: reindex completes on a pack where every species has one genome (GSTX has no multi-member group).

Usage: python3 reindex_singleton_species.py <genopack binary> <tmp parent dir>
"""
import pathlib
import random
import shutil
import subprocess
import sys
import tempfile

binary = sys.argv[1]

def run(*args):
    p = subprocess.run([binary, *map(str, args)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    assert p.returncode == 0, (p.returncode, p.stdout + p.stderr)
    return p.stdout + p.stderr

rng = random.Random(20261023)
with tempfile.TemporaryDirectory(prefix='genopack-singleton-', dir=sys.argv[2]) as tmp:
    d = pathlib.Path(tmp)
    rows = []
    for i in range(4):
        (d / f'g{i}.fa').write_text(f'>g{i}\n' + ''.join(rng.choice('ACGT') for _ in range(8000)) + '\n')
        rows.append(f'G{i}\td__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G sp{i}\t{d}/g{i}.fa\n')
    (d / 'in.tsv').write_text('accession\ttaxonomy\tfile_path\n' + ''.join(rows))
    base = d / 'base.gpk'
    run('build', '-i', d / 'in.tsv', '-o', base)

    for i, args in enumerate([(), ('--skch',)]):
        pack = d / f'case{i}.gpk'
        shutil.copyfile(base, pack)
        log = run('reindex', pack, *args)
        assert 'GSTX: pass-2 done' in log, log
        assert args or 'QUAL: 4 genomes computed' in log, log   # singleton QUAL records still written
        stat = run('stat', pack)
        assert any(line.split() == ['Genomes', '(live):', '4'] for line in stat.splitlines()), stat

print('reindex_singleton_species: PASS')
