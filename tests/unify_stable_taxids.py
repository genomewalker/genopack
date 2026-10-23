"""CLI regression: --stable-taxids-from keeps every prior taxid's meaning and never reissues a retired one."""
import csv
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

with tempfile.TemporaryDirectory(prefix='genopack-stable-', dir=sys.argv[2]) as tmp:
    d = pathlib.Path(tmp)
    fasta = d / 'a.fa'
    fasta.write_text('>a\n' + 'ACGT' * 1000 + '\n')
    lineages = {f'G{i}': f'd__Bacteria;p__P;c__C;o__O;f__F{i % 2};g__Gen{i % 3};s__Gen{i % 3} sp{i}' for i in range(8)}
    manifest = d / 'input.tsv'
    manifest.write_text('accession\ttaxonomy\tfile_path\n' + ''.join(f'{a}\t{t}\t{fasta}\n' for a, t in lineages.items()))
    pack = d / 'tiny.gpk'
    run('build', '-i', manifest, '-o', pack)
    # A one-genus NCBI taxdump: Bacteria resolves to its native taxid 2.
    ncbi = d / 'ncbi'
    ncbi.mkdir()
    rows = [(1, 1, 'no rank', 'root'), (2, 1, 'superkingdom', 'Bacteria')]
    (ncbi / 'nodes.dmp').write_text(''.join(f'{t}\t|\t{p}\t|\t{r}\t|\n' for t, p, r, n in rows))
    (ncbi / 'names.dmp').write_text(''.join(f'{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n' for t, p, r, n in rows))

    def unify(out, overrides=None, prior=None, fmt='columnar'):
        args = ['taxonomy', 'unify', pack, '--ncbi-taxdump', ncbi, '-f', fmt, '-o', d / out]
        if overrides is not None:
            (d / f'{out}.ov.tsv').write_text(''.join(f'{a}\t{t}\n' for a, t in sorted(overrides.items())))
            args += ['--lineage-overrides', d / f'{out}.ov.tsv']
        if prior:
            args += ['--stable-taxids-from', d / prior / 'taxonomy.tsv']
        log = run(*args)
        if fmt == 'taxdump':
            return None
        assert not prior or 'stable taxids:' in log, log
        return table(d / out / 'taxonomy.tsv', 'taxid'), table(d / out / 'acc2taxid.tsv', 'accession')

    # T1 baseline: the prior tree.
    prior, prior_accs = unify('prior')
    prior_max = max(map(int, prior))
    def climb(t):
        while prior[t]['parent_taxid'] != t:
            yield t
            t = prior[t]['parent_taxid']

    def check(nodes, accs, retired, moved):
        both = prior.keys() & nodes.keys()
        # R1: no prior taxid names another taxon, and kept nodes keep their parent.
        assert all(prior[t]['rank'] == nodes[t]['rank'] and prior[t]['name'] == nodes[t]['name'] and
                   prior[t]['parent_taxid'] == nodes[t]['parent_taxid'] for t in both), (prior, nodes)
        # Every prior node outside the retired set keeps its taxid; retired ids are not reissued.
        assert set(prior) - retired == both, (sorted(set(prior) - retired - both), sorted(retired & both))
        fresh = nodes.keys() - prior.keys()
        assert fresh and min(map(int, fresh)) > prior_max, (fresh, prior_max)
        assert {a for a in accs if accs[a]['taxid'] != prior_accs[a]['taxid']} == moved
        return fresh

    # T2 rename the node holding the highest prior taxid: a build that reissues retired ids (no floor at the
    # prior maximum) is caught.
    top = prior[str(prior_max)]
    token = {'domain': 'd', 'phylum': 'p', 'class': 'c', 'order': 'o', 'family': 'f', 'genus': 'g',
             'species': 's'}[top['rank']] + '__' + top['name']
    retired = {t for t in prior if str(prior_max) in climb(t)}
    moved = {a for a, r in prior_accs.items() if str(prior_max) in climb(r['taxid'])}
    renamed = {a: ';'.join(x.replace('__', '__Renamed', 1) if x == token else x for x in lineages[a].split(';'))
               for a in moved}
    nodes, accs = unify('rename', renamed, 'prior')
    fresh = check(nodes, accs, retired, moved)
    assert any(nodes[t]['name'] == 'Renamed' + top['name'] for t in fresh)
    unify('rename_dmp', renamed, 'prior', 'taxdump')
    assert {k: v['taxid'] for k, v in accs.items()} == \
        {k: v['taxid'] for k, v in table(d / 'rename_dmp/acc2taxid.dmp', 'accession').items()}

    # T3 move G0 (f__F0, g__Gen0, shared with G6) under g__Gen1 of the same family (held by G4): only the species
    # node is new; the new parent keeps its taxid.
    old = prior_accs['G0']['taxid']
    gen1_f0 = prior[prior_accs['G4']['taxid']]['parent_taxid']
    nodes, accs = unify('move', {'G0': lineages['G0'].replace('g__Gen0', 'g__Gen1')}, 'prior')
    fresh = check(nodes, accs, {old}, {'G0'})
    assert fresh == {accs['G0']['taxid']} and nodes[accs['G0']['taxid']]['parent_taxid'] == gen1_f0

    # T4 unchanged input reproduces the prior byte for byte; T5 the same input twice gives identical bytes.
    unify('same', None, 'prior')
    unify('rename2', renamed, 'prior')
    for a, b in (('same', 'prior'), ('rename2', 'rename')):
        for f in ('taxonomy.tsv', 'acc2taxid.tsv'):
            assert (d / a / f).read_bytes() == (d / b / f).read_bytes(), (a, b, f)

    # T6 the native NCBI taxid stays put.
    assert prior['2']['name'] == nodes['2']['name'] == 'Bacteria'

    assert 'Cannot open --stable-taxids-from' in run('taxonomy', 'unify', pack, '--stable-taxids-from', d / 'absent.tsv',
                                                     '-o', d / 'bad', fail=True)
print('stable taxids regression passed')
