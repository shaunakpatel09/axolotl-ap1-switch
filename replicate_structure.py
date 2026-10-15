"""
replicate_structure.py -- Analysis 13 (partial): the biological-replicate
structure of GSE162610, and whether the Section 2.9 protocol control is
confounded with animal.

The Jun/Fos ratios reported in Section 2.9 are pseudobulk aggregates over all
cells of a given type at a given timepoint.  GSE162610 is not a single sample
per timepoint: it contains ten libraries.  This module reads the deposited
barcode metadata and reports

  1. how many biological samples contribute to each timepoint,
  2. whether dissociation protocol and sequencing chemistry vary WITHIN a
     timepoint or are confounded with sample identity,
  3. how cell-type composition varies between samples of the same timepoint.

It does not recompute per-sample Jun/Fos ratios: that needs the count matrix,
which is not redistributed here.  What it establishes is whether a per-sample
analysis is possible at all, and what would be confounded with what.

NOTE ON THE FILE: the deposited metadata has a header of 29 names but 30
fields per row, the first being the unnamed barcode index.  Column names are
therefore shifted by one relative to a naive read.  This module corrects for
that explicitly.

Author: Shaunak Patel.  Built with AI assistance (Claude); reviewed and owned
by the author.
"""

import collections
import csv
import gzip
import json
import os
import sys

META = os.environ.get(
    'GSE162610_META',
    '/tmp/claude-0/-home-claude/1e6c9726-040e-5d3b-aa6e-2bb5f4748e2b/'
    'scratchpad/geo/GSE162610_barcode_metadata.tsv.gz')

OUT = 'replicate_structure.json'


def load(path=META):
    with gzip.open(path, 'rt') as f:
        r = csv.reader(f, delimiter='\t')
        hdr = next(r)
        rows = [row for row in r]
    # header is short by one: field 0 is the barcode index
    names = ['barcode'] + hdr
    if len(rows[0]) != len(names):
        raise ValueError(f'expected {len(names)} fields, got {len(rows[0])}')
    idx = {n: i for i, n in enumerate(names)}
    return rows, idx


def report(rows, idx):
    S = idx['orig.ident']     # e.g. 1dpi_sample1
    TP = idx['time']          # 1dpi / 3dpi / 7dpi / Uninjured
    DM = idx['dissociationMethod']
    CH = idx['chemistry']
    CT = idx['celltype']

    out = {}

    # 1. replicate structure
    by_tp = collections.defaultdict(set)
    cells = collections.Counter()
    for row in rows:
        by_tp[row[TP]].add(row[S])
        cells[row[S]] += 1
    out['samples_per_timepoint'] = {k: sorted(v) for k, v in sorted(by_tp.items())}
    out['cells_per_sample'] = dict(sorted(cells.items()))
    print('1. Biological samples per timepoint')
    for tp in sorted(by_tp):
        ss = sorted(by_tp[tp])
        print(f'   {tp:10s} n = {len(ss)}  ' +
              ', '.join(f'{s.split("_")[-1]} ({cells[s]})' for s in ss))

    # 2. is protocol confounded with sample?
    combo = collections.defaultdict(set)
    for row in rows:
        combo[row[S]].add((row[DM], row[CH]))
    print('\n2. Protocol and chemistry by sample')
    conf = {}
    for s in sorted(combo):
        pairs = sorted(combo[s])
        conf[s] = [list(p) for p in pairs]
        print(f'   {s:16s} {"; ".join(f"{a}/{b}" for a, b in pairs)}')
    out['protocol_by_sample'] = conf

    within = {}
    for tp in sorted(by_tp):
        prot = set()
        for s in by_tp[tp]:
            prot |= combo[s]
        within[tp] = sorted({p[0] for p in prot})
    out['protocols_within_timepoint'] = within
    print('\n   dissociation protocols present within each timepoint:')
    for tp, p in within.items():
        flag = '' if len(p) > 1 else '   <-- single protocol: not separable'
        print(f'     {tp:10s} {", ".join(p)}{flag}')

    # does any single sample carry both protocols?
    split = [s for s in combo if len({p[0] for p in combo[s]}) > 1]
    out['samples_spanning_both_protocols'] = sorted(split)
    print(f'\n   samples containing BOTH protocols: '
          f'{len(split)} of {len(combo)}'
          f'{" (" + ", ".join(sorted(split)) + ")" if split else ""}')

    # 3. cell-type composition per sample, for the two types used in 2.9
    targets = [c for c in {row[CT] for row in rows}
               if 'strocyte' in c or 'pendymal' in c]
    comp = collections.defaultdict(collections.Counter)
    for row in rows:
        comp[row[S]][row[CT]] += 1
    print('\n3. Composition of the two cell types analysed in Section 2.9')
    print(f'   (cell-type labels matched: {targets})')
    comp_out = {}
    for s in sorted(comp):
        tot = sum(comp[s].values())
        line = {t: comp[s].get(t, 0) for t in targets}
        line['_total'] = tot
        comp_out[s] = line
        frac = ', '.join(f'{t}: {comp[s].get(t,0)} '
                         f'({100*comp[s].get(t,0)/tot:.1f}%)' for t in targets)
        print(f'   {s:16s} {frac}')
    out['composition_by_sample'] = comp_out
    return out


def main():
    rows, idx = load()
    print(f'{len(rows)} barcodes\n')
    out = report(rows, idx)
    with open(OUT, 'w') as f:
        json.dump(out, f, indent=1)
    print(f'\nwrote {OUT}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
