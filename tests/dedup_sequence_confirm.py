"""CLI regression: dedup removes only sequence-identical genomes and its groups TSV says why.

Usage: python3 dedup_sequence_confirm.py <genopack binary> <tmp parent dir>
"""
import csv
import hashlib
import pathlib
import random
import subprocess
import sys
import tempfile

binary = sys.argv[1]

def run(*args):
    p = subprocess.run([binary, *map(str, args)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    assert p.returncode == 0, p.stdout + p.stderr
    return p.stdout + p.stderr

def wrap(seq, width):
    return ''.join(seq[i:i + width] + '\n' for i in range(0, len(seq), width))

rng = random.Random(20261023)
def rand_seq(n):
    return ''.join(rng.choice('ACGT') for _ in range(n))

a = rand_seq(6000)
b1 = rand_seq(6000)
b2 = b1[:3000] + ('A' if b1[3000] != 'A' else 'C') + b1[3001:]   # one SNP, same length
c = rand_seq(6000)
# accession -> (fasta text, completeness). Order fixes genome ids: C2 (one base longer) sits between C1 and C3.
genomes = {
    'A1': ('>A1 first\n' + wrap(a, 80), 90),
    'A2': ('>A2 other header\r\n' + wrap(a, 60).replace('\n', '\r\n'), 95),
    'B1': ('>B1\n' + wrap(b1, 80), 90),
    'B2': ('>B2\n' + wrap(b2, 80), 80),
    'C1': ('>C1\n' + wrap(c, 80), 80),
    'C2': ('>C2\n' + wrap(c + 'A', 80), 99),
    'C3': ('>C3\n' + wrap(c, 70), 70),
}

with tempfile.TemporaryDirectory(prefix='genopack-dedup-', dir=sys.argv[2]) as tmp:
    d = pathlib.Path(tmp)
    rows = []
    for acc, (text, comp) in genomes.items():
        (d / f'{acc}.fa').write_text(text)
        rows.append(f'{acc}\td__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G sp\t{d / acc}.fa\t{comp}\n')
    (d / 'input.tsv').write_text('accession\ttaxonomy\tfile_path\tcompleteness\n' + ''.join(rows))
    pack = d / 'tiny.gpk'
    run('build', '-i', d / 'input.tsv', '-o', pack)

    def groups(name, *extra):
        out = d / name
        run('dedup', pack, '--groups', out, *extra)
        with out.open() as f:
            return {r['accession']: r for r in csv.DictReader(f, delimiter='\t')}

    before = hashlib.sha256(pack.read_bytes()).hexdigest()
    g = groups('dry.tsv', '--dry-run')
    assert hashlib.sha256(pack.read_bytes()).hexdigest() == before, 'dry run modified the archive'

    # Fixture preconditions: each pair really is a fingerprint+length candidate.
    for x, y in (('A1', 'A2'), ('B1', 'B2'), ('C1', 'C3')):
        assert x in g and y in g and g[x]['candidate_group'] == g[y]['candidate_group'], (x, y, g)
    assert 'C2' not in g, g

    action = {acc: r['action'] for acc, r in g.items()}
    assert action == {'A1': 'remove', 'A2': 'keep', 'B1': 'distinct', 'B2': 'distinct',
                      'C1': 'keep', 'C3': 'remove'}, action
    assert g['A1']['sequence_hash'] == g['A2']['sequence_hash'] != g['B1']['sequence_hash'] != g['B2']['sequence_hash']

    groups('real.tsv')
    out = d / 'all.fa'
    run('extract', pack, '-o', out)
    kept = {line[1:].split()[0] for line in out.read_text().splitlines() if line.startswith('>')}
    assert kept == {'A2', 'B1', 'B2', 'C1', 'C2'}, kept

    again = groups('again.tsv', '--dry-run')
    assert {acc: r['action'] for acc, r in again.items()} == {'B1': 'distinct', 'B2': 'distinct'}, again

print('dedup_sequence_confirm: PASS')
