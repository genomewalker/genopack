"""CLI regression: a single --sketch-kmers value sets k for build, enrich and reindex --skch.

Usage: python3 sketch_kmers_single.py <genopack binary> <tmp parent dir>
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
with tempfile.TemporaryDirectory(prefix='genopack-kmers-', dir=sys.argv[2]) as tmp:
    d = pathlib.Path(tmp)
    rows = []
    for i in range(4):
        (d / f'g{i}.fa').write_text(f'>g{i}\n' + ''.join(rng.choice('ACGT') for _ in range(8000)) + '\n')
        rows.append(f'G{i}\td__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G sp\t{d}/g{i}.fa\n')
    tsv = d / 'in.tsv'
    tsv.write_text('accession\ttaxonomy\tfile_path\n' + ''.join(rows))

    def build(name, *args):
        pack = d / name
        run('build', '-i', tsv, '-o', pack, *args)
        return skch_params(pack)

    assert build('k21.gpk', '--sketch-kmers', 21) == [([21], 7)]
    assert build('k21dup.gpk', '--sketch-kmers', '21,21') == [([21], 7)]
    assert build('k16_21.gpk', '--sketch-kmers', '16,21') == [([16, 21], 5)]
    assert build('kmer24.gpk', '--sketch-kmers', '', '--sketch-kmer', 24) == [([24], 8)]

    thin = d / 'thin.gpk'
    run('build', '-i', tsv, '-o', thin, '--thin')
    run('enrich', '-i', thin, '-o', d / 'enriched.gpk', '--sketch-kmers', 18)
    assert skch_params(d / 'enriched.gpk') == [([18], 6)], skch_params(d / 'enriched.gpk')

    re = d / 'reindex.gpk'
    shutil.copyfile(d / 'k21.gpk', re)
    run('reindex', re, '--skch', '--sketch-kmers', 24)
    assert skch_params(re) == [([24], 8)], skch_params(re)

print('sketch_kmers_single: PASS')
