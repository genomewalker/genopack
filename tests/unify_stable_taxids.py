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
    run('taxonomy', 'unify', pack, '-f', 'columnar', '-o', d / 'prior')
    prior = table(d / 'prior/taxonomy.tsv', 'taxid')
    prior_accs = table(d / 'prior/acc2taxid.tsv', 'accession')
    prior_max = max(map(int, prior))
    # Rename the node holding the highest prior taxid: that id and its subtree's retire, so a build that reissues
    # retired ids (no floor at the prior maximum) is caught.
    def climb(t):
        while prior[t]['parent_taxid'] != t:
            yield t
            t = prior[t]['parent_taxid']
    top = prior[str(prior_max)]
    token = {'domain': 'd', 'phylum': 'p', 'class': 'c', 'order': 'o', 'family': 'f', 'genus': 'g',
             'species': 's'}[top['rank']] + '__' + top['name']
    moved = {a for a, r in prior_accs.items() if str(prior_max) in climb(r['taxid'])}
    assert moved
    overrides = d / 'overrides.tsv'
    overrides.write_text(''.join(f'{a}\t' + ';'.join(x.replace('__', '__Renamed', 1) if x == token else x
                                                      for x in lineages[a].split(';')) + '\n' for a in sorted(moved)))
    for fmt in ('columnar', 'taxdump'):
        log = run('taxonomy', 'unify', pack, '--lineage-overrides', overrides,
                  '--stable-taxids-from', d / 'prior/taxonomy.tsv', '-f', fmt, '-o', d / fmt)
        assert 'stable taxids:' in log and f'new from {prior_max + 1}' in log, log
    nodes = table(d / 'columnar/taxonomy.tsv', 'taxid')
    accs = table(d / 'columnar/acc2taxid.tsv', 'accession')
    # R1: no prior taxid names another taxon.
    both = prior.keys() & nodes.keys()
    assert all((prior[t]['rank'], prior[t]['name']) == (nodes[t]['rank'], nodes[t]['name']) for t in both), (prior, nodes)
    # New nodes are numbered above the prior maximum, and the renamed node's old id is retired.
    fresh = nodes.keys() - prior.keys()
    assert fresh and min(map(int, fresh)) > prior_max and str(prior_max) not in nodes, (fresh, prior_max)
    assert any(nodes[t]['name'] == 'Renamed' + top['name'] for t in fresh)
    # Exactly the accessions under the renamed node change taxid.
    assert {a for a in accs if accs[a]['taxid'] != prior_accs[a]['taxid']} == moved
    assert {k: v['taxid'] for k, v in accs.items()} == {k: v['taxid'] for k, v in table(d / 'taxdump/acc2taxid.dmp', 'accession').items()}
    assert 'Cannot open --stable-taxids-from' in run('taxonomy', 'unify', pack, '--stable-taxids-from', d / 'absent.tsv',
                                                     '-o', d / 'bad', fail=True)
print('stable taxids regression passed')
