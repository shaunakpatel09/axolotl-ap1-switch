"""
mouse_ratio.py -- measure the AP-1 production ratio in mouse astrocytes across
a spinal cord injury time course, and place it on the regime map.

THE TEST
--------
regime.py predicts that a cell is switchable only while its baseline miR-200a
lies inside the two-state window, and that the window slides to higher
miR-200a as c-Jun production rises relative to c-Fos. Past a c-Jun/c-Fos
production ratio of roughly 1.3 (at the reference miR-200a of 0.58) the cell
is left below the window: locked in the high-c-Jun scar state.

So: does the mouse astrocyte c-Jun/c-Fos ratio rise across the injury time
course, and does it cross that boundary?

This is a within-dataset comparison -- one platform, one annotation, one
experiment -- so it needs no cross-species normalisation argument.

WHY A PSEUDOBULK RATIO
----------------------
Counts are summed across cells FIRST, then divided:

    a_J/a_F  =  sum(Jun counts) / sum(Fos counts)

Per-cell library size cancels out of the ratio, so no normalisation is needed
and none is applied. Log-normalised values (Seurat's data slot,
AverageExpression, DotPlot) must NOT be used: the mean of logs is not the log
of the mean, and the error compresses exactly the differences being measured.

CAVEAT: the model's a_J is a PROTEIN production rate. mRNA abundance times
translation rate gives protein production, so an mRNA ratio is the production
ratio only if translation efficiency is comparable across these factors. That
assumption is stated, not hidden.

DATA
----
GSE162610 (Milich et al. 2021, J Exp Med). Download these five files into one
folder and do not decompress them:

    GSE162610_sci_mat.mtx.gz
    GSE162610_genes.tsv.gz
    GSE162610_barcodes.tsv.gz
    GSE162610_barcode_metadata.tsv.gz

Then:  inspect(folder)          see what the columns and cell types are called
       counts = extract(folder) stream the matrix, keep 3 genes
       report(counts, meta)     ratios by timepoint, placed on the map

Author: Shaunak Patel. Built with AI assistance (Claude); reviewed and owned
by the author.
"""

import gzip
import os

import numpy as np
import pandas as pd

GENES = ('Jun', 'Fos', 'Junb')          # mouse symbols; inspect() verifies
MTX = 'GSE162610_sci_mat.mtx.gz'
GENE_FILE = 'GSE162610_genes.tsv.gz'
BC_FILE = 'GSE162610_barcodes.tsv.gz'
META_FILE = 'GSE162610_barcode_metadata.tsv.gz'


def _read_tsv(path, **kw):
    with gzip.open(path, 'rt') as f:
        return pd.read_csv(f, sep='\t', **kw)


def _mtx_header(path):
    """Dimensions line of a MatrixMarket file, skipping % comments."""
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('%'):
                continue
            a, b, c = line.split()
            return int(a), int(b), int(c)


def inspect(folder):
    """Look before computing: what are the files, the orientation, the column
    names, and the cell types actually called?"""
    gp = os.path.join(folder, GENE_FILE)
    genes = _read_tsv(gp, header=None).iloc[:, 0].astype(str)
    bcs = _read_tsv(os.path.join(folder, BC_FILE), header=None)
    meta = _read_tsv(os.path.join(folder, META_FILE))
    nr, nc, nnz = _mtx_header(os.path.join(folder, MTX))

    print(f'matrix      {nr} x {nc}, {nnz} non-zero entries')
    print(f'genes.tsv   {len(genes)} rows      barcodes.tsv {len(bcs)} rows')
    orient = ('genes x cells' if nr == len(genes) else
              'cells x genes' if nc == len(genes) else 'UNRECOGNISED')
    print(f'orientation {orient}')
    print(f'\nbarcode_metadata columns:\n  {list(meta.columns)}')
    for col in meta.columns:
        v = meta[col]
        if v.dtype == object and 2 <= v.nunique() <= 40:
            print(f'\n  {col}  ({v.nunique()} values)')
            print('    ' + ', '.join(map(str, sorted(v.unique())[:40])))

    print('\ngene symbols sought:')
    low = {g.lower(): i for i, g in enumerate(genes)}
    for g in GENES:
        i = low.get(g.lower())
        print(f'  {g:6s} ' + (f'found at row {i + 1}' if i is not None
                              else 'NOT FOUND -- check the naming convention'))
    hits = [g for g in genes if g.lower().startswith(('jun', 'fos'))]
    print(f'  all Jun/Fos-family symbols present: {sorted(set(hits))}')
    return genes, meta, (nr, nc, nnz), orient


def extract(folder, genes_wanted=GENES):
    """Stream the matrix and keep only the wanted genes. Returns a DataFrame
    of per-cell counts, indexed in barcode order.

    The matrix is never loaded: it is read one line at a time, so this costs
    almost no memory regardless of file size."""
    gp = os.path.join(folder, GENE_FILE)
    genes = _read_tsv(gp, header=None).iloc[:, 0].astype(str)
    low = {g.lower(): i for i, g in enumerate(genes)}
    rows = {}
    for g in genes_wanted:
        i = low.get(g.lower())
        if i is None:
            raise KeyError(f'{g} not in {GENE_FILE} -- run inspect() first')
        rows[i + 1] = g                       # MatrixMarket is 1-based

    nr, nc, nnz = _mtx_header(os.path.join(folder, MTX))
    n_cells = nc if nr == len(genes) else nr
    gene_axis_is_row = (nr == len(genes))
    out = {g: np.zeros(n_cells, dtype=np.int64) for g in genes_wanted}

    with gzip.open(os.path.join(folder, MTX), 'rt') as f:
        seen_dims = False
        for line in f:
            if line[0] == '%':
                continue
            if not seen_dims:
                seen_dims = True
                continue
            a, b, v = line.split()
            r, c = (int(a), int(b)) if gene_axis_is_row else (int(b), int(a))
            g = rows.get(r)
            if g is not None:
                out[g][c - 1] += int(v)

    df = pd.DataFrame(out)
    bcs = _read_tsv(os.path.join(folder, BC_FILE), header=None)
    df.insert(0, 'barcode', bcs.iloc[:, 0].astype(str).values)
    return df


def join_meta(counts, folder, barcode_col=None):
    """Attach the per-cell metadata, matched on the barcode STRING.

    In GSE162610 the metadata's first column is unnamed and holds the barcode,
    so pandas loads it as the index; barcode_col=None reads it from there.
    The join is on the barcode itself, never on row order -- the two files
    happen to agree here, but relying on that would fail silently on any
    dataset where they don't."""
    if barcode_col is None:
        meta = _read_tsv(os.path.join(folder, META_FILE), index_col=0)
        meta = meta.rename_axis('barcode').reset_index()
    else:
        meta = _read_tsv(os.path.join(folder, META_FILE))
        meta = meta.rename(columns={barcode_col: 'barcode'})
    meta['barcode'] = meta['barcode'].astype(str)
    out = counts.merge(meta, on='barcode', how='inner')
    if len(out) < 0.95 * min(len(counts), len(meta)):
        raise ValueError(f'only {len(out)} of {len(counts)} cells matched the '
                         'metadata by barcode -- check barcode_col')
    return out


def ratios(df, celltype_col, group_col, celltypes=('Astrocyte',),
           order=None, genes=GENES):
    """Pseudobulk AP-1 ratios within the chosen cell type, by group.

    Counts are summed over cells before dividing, so per-cell library size
    cancels and no normalisation is applied."""
    J, F, B = genes
    sub = df[df[celltype_col].astype(str).isin(celltypes)]
    if not len(sub):
        raise ValueError(f'no cells matched {celltypes} in {celltype_col}; '
                         f'available: {sorted(df[celltype_col].unique())[:20]}')
    g = sub.groupby(group_col)
    tab = pd.DataFrame({
        'n_cells': g.size(),
        f'{J}_counts': g[J].sum(),
        f'{F}_counts': g[F].sum(),
        f'{B}_counts': g[B].sum(),
    })
    tab['a_J/a_F'] = tab[f'{J}_counts'] / tab[f'{F}_counts']
    tab['a_B/a_F'] = tab[f'{B}_counts'] / tab[f'{F}_counts']
    if order:
        tab = tab.reindex([o for o in order if o in tab.index])
    return tab


def place(tab, m0=0.58):
    """Classify each row's ratio against the model's regime map."""
    import regime
    print(f'\nregime at baseline miR-200a {m0:.2f}:')
    print(f'  {"group":>12s}{"a_J/a_F":>10s}{"a_B/a_F":>10s}'
          f'{"window":>14s}   regime')
    out = {}
    for name, row in tab.iterrows():
        x, y = row['a_J/a_F'], row['a_B/a_F']
        if not np.isfinite(x) or not np.isfinite(y):
            print(f'  {str(name):>12s}{x:>10.2f}{y:>10.2f}'
                  '          --   not computable (zero denominator)')
            continue
        p = regime.ratio_params(x, y)
        kind, info = regime.classify(p, m0)
        w = info['w']
        span = ('none' if w['lo'] is None else f"{w['lo']:.2f}-{w['hi']:.2f}")
        print(f'  {str(name):>12s}{x:>10.2f}{y:>10.2f}{span:>14s}   {kind}')
        out[name] = (x, y, kind, info)
    return out


# ----------------------------------------------------------------------
# the result
# ----------------------------------------------------------------------

# Measured from GSE162610 by the pipeline above (extract -> join_meta ->
# ratios). Recorded here so the figure is reproducible without re-streaming
# the 458 MB matrix. a_J/a_F = sum(Jun counts) / sum(Fos counts).
TIMES = ('Uninjured', '1dpi', '3dpi', '7dpi')
MEASURED = {
    # label:                  (a_J/a_F by timepoint, a_B/a_F by timepoint)
    'Astrocyte':             ((0.643, 1.410, 0.349, 0.500),
                              (0.499, 0.858, 0.591, 0.828)),
    'Ependymal':             ((0.388, 0.335, 0.267, 0.400),
                              (0.227, 0.288, 0.302, 0.370)),
    'Ependymal (Enriched)':  ((0.401, 0.349, 0.252, 0.413),
                              (0.215, 0.268, 0.309, 0.364)),
    'Ependymal (Standard)':  ((0.263, 0.283, 0.285, 0.363),
                              (0.345, 0.358, 0.294, 0.385)),
}
N_CELLS = {'Astrocyte': (123, 492, 197, 118),
           'Ependymal': (2678, 2869, 1220, 950),
           'Ependymal (Enriched)': (2048, 1987, 650, 567),
           'Ependymal (Standard)': (630, 882, 570, 383)}
COUNTS_J = {'Astrocyte': (446, 1603, 130, 345),
            'Ependymal': (70267, 42122, 8482, 8203)}
COUNTS_F = {'Astrocyte': (694, 1137, 372, 690),
            'Ependymal': (181074, 125870, 31737, 20525)}


def ci(j, f):
    """Approximate 95% interval on a ratio of two Poisson counts. Ignores
    cell-to-cell overdispersion, so it is a FLOOR on the true uncertainty."""
    r = j / f
    s = np.sqrt(1.0 / j + 1.0 / f)
    return r * np.exp(-1.96 * s), r * np.exp(1.96 * s)


def figure10(m0=0.58, xs=None, ys=None, out_dir='figures_round2'):
    """Figure 10: the model's boundary, and the measured data refusing to
    reach it."""
    import matplotlib.pyplot as plt
    import regime

    xs = np.linspace(0.25, 3.0, 16) if xs is None else xs
    ys = np.linspace(0.10, 1.50, 12) if ys is None else ys
    W, CJ, G, LO = regime.load_map(xs, ys)
    j = 3
    lo, hi = LO[:, j], LO[:, j] + W[:, j]
    inside = (lo <= m0) & (m0 <= hi)
    boundary = xs[~inside & (xs > 1.0)].min() if (~inside).any() else np.nan

    fig, ax = plt.subplots(1, 2, figsize=(14.5, 5.2))
    col = {'Astrocyte': '#D85A30', 'Ependymal': '#1D9E75',
           'Ependymal (Enriched)': '#534AB7',
           'Ependymal (Standard)': '#8E86D8'}

    # A -- where the data lands on the model's map
    a = ax[0]
    a.fill_between(xs, lo, hi, color='purple', alpha=0.22,
                   label='two-state window (switchable)')
    a.plot(xs, lo, color='purple', lw=1)
    a.plot(xs, hi, color='purple', lw=1)
    a.axhline(m0, color='k', ls='--', lw=1.4,
              label=f'miR-200a assumed = {m0:.2f} (NOT measured)')
    a.axvline(boundary, color='#D85A30', lw=1.6, ls='-.',
              label=f'locked beyond a_J/a_F = {boundary:.2f}')
    for lab in ('Astrocyte', 'Ependymal'):
        v = MEASURED[lab][0]
        a.plot(v, [m0] * len(v), 'o', color=col[lab], ms=9, mec='k', mew=0.6,
               label=f'{lab}, measured (4 timepoints)', zorder=5)
    a.set_xlim(0.2, 3.05)
    a.set_xlabel('c-Jun production / c-Fos production  (Jun / Fos counts)')
    a.set_ylabel('miR-200a (multiples of its repression constant)')
    a.set_title('A  Seven of eight measured points fall inside the\n'
                'switchable band; only astrocytes at 1 dpi cross',
                loc='left', fontsize=10)
    a.legend(fontsize=7, loc='upper left')

    # B -- the measurement is protocol-dominated, and still nowhere near
    b = ax[1]
    pos = np.arange(len(TIMES))
    for lab in ('Astrocyte', 'Ependymal (Enriched)', 'Ependymal (Standard)'):
        b.plot(pos, MEASURED[lab][0], 'o-', color=col[lab], label=lab, ms=7)
    for k, t in enumerate(TIMES):                      # Poisson floor on error
        for lab in ('Astrocyte',):
            l, u = ci(COUNTS_J[lab][k], COUNTS_F[lab][k])
            b.plot([k, k], [l, u], color=col[lab], lw=1.4, alpha=0.7)
    b.axhline(boundary, color='#D85A30', lw=1.6, ls='-.',
              label=f'locked boundary = {boundary:.2f}')
    top7 = max(MEASURED[l][0][3] for l in
               ('Astrocyte', 'Ependymal (Enriched)', 'Ependymal (Standard)'))
    b.annotate('', xy=(3.18, boundary), xytext=(3.18, top7),
               arrowprops=dict(arrowstyle='<->', color='0.35', lw=1.2))
    b.text(3.24, 0.5 * (boundary + top7), f'{boundary / top7:.1f}x short\nat 7 dpi',
           fontsize=8, color='0.35', va='center')
    b.set_xlim(-0.3, 3.9)
    b.set_xticks(pos)
    b.set_xticklabels(TIMES)
    b.set_ylim(0, 1.55)
    b.set_ylabel('a_J/a_F  (Jun counts / Fos counts)')
    b.set_title('B  Dissociation protocol moves the ratio as much as\n'
                'injury does (bars: Poisson floor on uncertainty)',
                loc='left', fontsize=10)
    b.legend(fontsize=7, loc='upper left')

    fig.suptitle('Figure 10. The regime prediction tested against mouse '
                 'spinal cord injury data (GSE162610) — and refuted',
                 fontsize=12)
    fig.tight_layout()
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, 'fig10_regime_vs_data.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print('saved', path)
    return fig


def report(folder, celltype_col, group_col, celltypes=('Astrocyte',),
           order=None, barcode_col=None, m0=0.58):
    """Whole pipeline: extract, join, tabulate, place."""
    counts = extract(folder)
    df = join_meta(counts, folder, barcode_col)
    tab = ratios(df, celltype_col, group_col, celltypes, order)
    print(tab.round(3).to_string())
    place(tab, m0)
    return tab
