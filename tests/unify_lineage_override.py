"""CLI regression: lineage overrides behave like patch, without changing the pack."""
import csv
import hashlib
import pathlib
import subprocess
import sys
import tempfile

binary = sys.argv[1]

def run(*args, fail=False):
    p = subprocess.run([binary, *map(str, args)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    assert (p.returncode != 0) if fail else (p.returncode == 0), p.stdout + p.stderr
    return p.stdout + p.stderr

def table(path, key):
    with path.open() as f:
        return {r[key]: r for r in csv.DictReader(f, delimiter='\t')}

with tempfile.TemporaryDirectory(prefix='genopack-unify-', dir=sys.argv[2]) as tmp:
    d = pathlib.Path(tmp)
    lineage = 'd__Bacteria;p__P;c__C;o__O;f__F;g__Wrong;s__Right alpha'
    right = lineage.replace('g__Wrong', 'g__Right')
    fasta = d / 'a.fa'
    fasta.write_text('>a\n' + 'ACGT' * 1000 + '\n')
    manifest = d / 'input.tsv'
    manifest.write_text('accession\ttaxonomy\tfile_path\n' + ''.join(
        f'{acc}\t{tax}\t{fasta}\n' for acc, tax in [('A', lineage), ('B', right)]))
    pack = d / 'tiny.gpk'
    run('build', '-i', manifest, '-o', pack)
    original = hashlib.sha256(pack.read_bytes()).digest()
    overrides = d / 'overrides.tsv'
    overrides.write_text('A\t' + right + '\n')
    run('taxonomy', 'unify', pack, '-f', 'columnar', '-o', d / 'before')
    for fmt in ('columnar', 'taxdump'):
        log = run('taxonomy', 'unify', pack, '--lineage-overrides', overrides,
                  '-f', fmt, '-o', d / fmt)
        assert 'applied 1 lineage overrides' in log
    nodes = table(d / 'columnar/taxonomy.tsv', 'taxid')
    accs = table(d / 'columnar/acc2taxid.tsv', 'accession')
    before = table(d / 'before/taxonomy.tsv', 'taxid')
    assert any(n['name'] == 'Wrong' for n in before.values())
    assert not any(n['name'] == 'Wrong' for n in nodes.values())
    assert accs['B'] == table(d / 'before/acc2taxid.tsv', 'accession')['B'], (before, nodes, accs)
    leaf = nodes[accs['A']['taxid']]
    # Unify ignores rank-S tokens for both original and overridden genomes.
    assert leaf['rank'] == 'species'
    species = leaf
    assert species['rank'] == 'species' and species['name'] == 'Right alpha'
    assert nodes[species['parent_taxid']]['name'] == 'Right'
    names = {}
    for line in (d / 'taxdump/names.dmp').read_text().splitlines():
        f = [v.strip() for v in line.split('|')]
        names[f[0]] = f[1]
    ranks = {'domain': 'superkingdom', 'no_rank': 'no rank'}
    taxdump_nodes = {}
    for line in (d / 'taxdump/nodes.dmp').read_text().splitlines():
        f = [v.strip() for v in line.split('|')]
        taxdump_nodes[f[0]] = (f[1], f[2], names[f[0]])
    assert taxdump_nodes == {k: (v['parent_taxid'], ranks.get(v['rank'], v['rank']), v['name']) for k, v in nodes.items()}
    assert {k: v['taxid'] for k, v in accs.items()} == {k: v['taxid'] for k, v in table(d / 'taxdump/acc2taxid.dmp', 'accession').items()}
    # Compare directly with patch's normalization, including explicit ten-rank S__.
    for tax in (right, 'd__Bacteria;l__;k__;p__P;c__C;o__O;f__F;g__Right;s__Right alpha;S__isolate A'):
        overrides.write_text('A\t' + tax + '\n')
        run('taxonomy', 'unify', pack, '--lineage-overrides', overrides, '-f', 'columnar', '-o', d / 'ten')
        patch_input = d / 'patch-input.tsv'
        patch_input.write_text('accession\ttaxonomy\nA\t' + tax + '\n')
        run('taxonomy', 'patch', '--patch', patch_input, '--tsv', manifest, '--tsv-out', d / 'patched.tsv')
        patched = table(d / 'patched.tsv', 'accession')['A']['taxonomy']
        overrides.write_text('A\t' + patched + '\n')
        run('taxonomy', 'unify', pack, '--lineage-overrides', overrides, '-f', 'columnar', '-o', d / 'normalized')
        assert (d / 'ten/taxonomy.tsv').read_bytes() == (d / 'normalized/taxonomy.tsv').read_bytes()
        assert table(d / 'ten/acc2taxid.tsv', 'accession') == accs
    # The override is also consumed before NCBI backbone grafting/excision.
    ncbi = d / 'ncbi'
    ncbi.mkdir()
    ncbi_rows = [(1, 1, 'no rank', 'root'), (2, 1, 'superkingdom', 'Bacteria'),
                 (3, 2, 'genus', 'NCBI genus')]
    (ncbi / 'nodes.dmp').write_text(''.join(
        f'{tid}\t|\t{parent}\t|\t{rank}\t|\n' for tid, parent, rank, name in ncbi_rows))
    (ncbi / 'names.dmp').write_text(''.join(
        f'{tid}\t|\t{name}\t|\t\t|\tscientific name\t|\n' for tid, parent, rank, name in ncbi_rows))
    overrides.write_text('A\t' + right + '\n')
    for replace in (False, True):
        flags = ['--gtdb-replaces-prok'] if replace else []
        for fmt in ('columnar', 'taxdump'):
            dest = d / f'backbone-{replace}-{fmt}'
            run('taxonomy', 'unify', pack, '--lineage-overrides', overrides,
                '--ncbi-taxdump', ncbi, '--ncbi-backbone', *flags, '-f', fmt, '-o', dest)
            accfile = 'acc2taxid.tsv' if fmt == 'columnar' else 'acc2taxid.dmp'
            actual = table(dest / accfile, 'accession')
            baseline = d / f'baseline-{replace}-{fmt}'
            run('taxonomy', 'unify', pack, '--ncbi-taxdump', ncbi, '--ncbi-backbone',
                *flags, '-f', fmt, '-o', baseline)
            assert actual['B'] == table(baseline / accfile, 'accession')['B']
            assert actual['A']['taxid'] == actual['B']['taxid']
            if fmt == 'columnar':
                scaffold = table(dest / 'taxonomy.tsv', 'taxid')
                assert ('3' in scaffold) == (not replace)
                assert not any(n['name'] == 'Wrong' for n in scaffold.values())
    for text, error in [('A\t' + right + '\nA\t' + right + '\n', 'duplicate accession'),
                        ('ABSENT\t' + right + '\n', 'absent from input packs'),
                        ('A\t' + right + '\textra\n', 'expected accession TAB lineage')]:
        overrides.write_text(text)
        assert error in run('taxonomy', 'unify', pack, '--lineage-overrides', overrides, '-o', d / 'bad', fail=True)
    assert hashlib.sha256(pack.read_bytes()).digest() == original
print('lineage override regression passed')
