"""
kd_test.py -- Analysis 9: does the model survive the conventional c-Jun
homodimer dissociation constant?

Section 4.1 states that shallow bistability exists only for K_d(Jun-Jun) of
roughly 30-150 nM, and is absent at 10 nM -- the conventional value implied by
the Fos-Jun affinity being ~100x tighter.  That claim was made at the reference
parameter set only.  This module tests it two ways:

  A. SWEEP.  Window width and state separation across K_d, under two
     conventions: varying the c-Jun homodimer alone, and varying every
     Jun-containing weak pairing together (JJ, BB, JB) as the reference
     parameterisation does.  The two are reported separately because the
     distinction matters: if bistability needs all three weakened, the
     assumption is broader than a statement about the homodimer.

  B. COMPENSATION.  At the conventional K_d = 10 nM, can other parameters
     recover shallow bistability?  Latin hypercube over the feedback, basal
     production and degradation terms, scoring each set for a two-state window
     with separation inside the 2-10x band the data require.

K_d = koff / kon, with kon = 10 throughout, so koff = 10 * K_d.

Author: Shaunak Patel.  Built with AI assistance (Claude); reviewed and owned
by the author.
"""

import json
import os
import sys

import numpy as np

from ap1_model_v3 import SPECIES, N, steady_state, totals
from round2 import ref_params, robust_window

SEED = 20261015
OUT = 'kd_test_results.json'

KON = 10.0
KD_CONVENTIONAL = 10.0      # 100x weaker than Fos-Jun (K_d 0.1 nM)
KD_REFERENCE = 100.0        # the value used in this work

KD_GRID = [5, 7.5, 10, 15, 20, 30, 40, 60, 80, 100, 125, 150, 175,
           200, 250, 300, 400]


def set_kd(p, kd, tie=True):
    q = dict(p)
    q['koff_JJ'] = KON * kd
    if tie:
        q['koff_BB'] = KON * kd
        q['koff_JB'] = KON * kd
    return q


def window_and_separation(p):
    """Two-state window width and total-c-Jun separation at its midpoint."""
    try:
        w = robust_window(p)
    except Exception:
        return 0.0, np.nan
    if w is None or not np.isfinite(w.get('width', np.nan)) or w['width'] <= 0:
        return 0.0, np.nan
    q = dict(p)
    q['mir'] = w['m_mid']
    y = np.zeros(N)
    y[SPECIES.index('J')] = 100.0
    hi = steady_state(q, y0=y, t_end=5000.0)
    lo = steady_state(q, y0=np.zeros(N), t_end=5000.0)
    chi, clo = totals(hi)['cJun_total'], totals(lo)['cJun_total']
    if chi < clo:
        chi, clo = clo, chi
    sep = chi / max(clo, 1e-12)
    if sep < 1.25:                      # branches collapsed: not two states
        return 0.0, np.nan
    return float(w['width']), float(sep)


def sweep(verbose=True):
    p = ref_params()
    out = {'tied': [], 'jj_only': []}
    for tie, key in ((True, 'tied'), (False, 'jj_only')):
        for kd in KD_GRID:
            width, sep = window_and_separation(set_kd(p, kd, tie=tie))
            out[key].append({'kd': kd, 'width': width, 'sep': sep})
            if verbose:
                s = f'{sep:.2f}' if np.isfinite(sep) else '--'
                print(f'  [{key:8s}] K_d {kd:6.1f} nM  width {width:5.3f}  '
                      f'separation {s}', flush=True)
    return out


# ------------------------------------------------------------------ B

COMP_SPEC = [
    ('b_J_by_JJ', 4.0, 400.0, 'log'),     # homodimer feedback strength
    ('a_J',       1.0, 60.0,  'log'),     # basal c-Jun production
    ('Km_auto',   0.5, 50.0,  'log'),     # feedback half-maximal threshold
    ('n_auto',    1.0, 6.0,   'lin'),     # feedback cooperativity
    ('g_J',       0.1, 3.0,   'log'),     # c-Jun degradation
    ('a_F',       1.0, 40.0,  'log'),     # c-Fos production
]


def lhs(n, k, rng):
    u = (rng.permuted(np.tile(np.arange(n), (k, 1)), axis=1).T
         + rng.random((n, k))) / n
    return u


def compensation(n_sets=400, kd=KD_CONVENTIONAL, tie=True, seed=SEED,
                 verbose=True):
    """At the conventional K_d, can other parameters recover shallow
    bistability?  Returns one record per sampled set."""
    rng = np.random.default_rng(seed)
    u = lhs(n_sets, len(COMP_SPEC), rng)
    base = set_kd(ref_params(), kd, tie=tie)
    rows = []
    for i in range(n_sets):
        q = dict(base)
        for j, (name, lo, hi, kind) in enumerate(COMP_SPEC):
            q[name] = (lo + u[i, j] * (hi - lo) if kind == 'lin'
                       else 10 ** (np.log10(lo) + u[i, j]
                                   * (np.log10(hi) - np.log10(lo))))
        width, sep = window_and_separation(q)
        rows.append({'width': width, 'sep': sep,
                     **{nm: float(q[nm]) for nm, *_ in COMP_SPEC}})
        if verbose and (i + 1) % 50 == 0:
            nb = sum(1 for r in rows if r['width'] > 0)
            print(f'  compensation {i+1}/{n_sets}: {nb} bistable so far',
                  flush=True)
    return rows


def summarise(rows):
    bist = [r for r in rows if r['width'] > 0 and np.isfinite(r['sep'])]
    shallow = [r for r in bist if 2.0 <= r['sep'] <= 10.0]
    return {
        'n': len(rows),
        'bistable': len(bist),
        'shallow_2_10': len(shallow),
        'sep_median': float(np.median([r['sep'] for r in bist])) if bist else None,
        'sep_min': float(min((r['sep'] for r in bist), default=np.nan)),
        'sep_max': float(max((r['sep'] for r in bist), default=np.nan)),
    }


def main():
    print('A. K_d sweep at the reference parameter set')
    sw = sweep()
    print('\nB. compensation search at the conventional K_d = 10 nM '
          '(all Jun pairings tied)')
    comp_tied = compensation(tie=True)
    s_tied = summarise(comp_tied)
    print('   ', s_tied)
    print('\nB2. compensation search at K_d = 10 nM for the homodimer only')
    comp_jj = compensation(tie=False, seed=SEED + 1)
    s_jj = summarise(comp_jj)
    print('   ', s_jj)

    res = {'seed': SEED, 'kd_grid': KD_GRID, 'sweep': sw,
           'compensation_tied': s_tied, 'compensation_jj_only': s_jj,
           'compensation_spec': [[a, b, c, d] for a, b, c, d in COMP_SPEC],
           'shallow_examples': [r for r in comp_tied
                                if r['width'] > 0 and np.isfinite(r['sep'])
                                and 2.0 <= r['sep'] <= 10.0][:10]}
    with open(OUT, 'w') as f:
        json.dump(res, f, indent=1)
    print(f'\nwrote {OUT}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
