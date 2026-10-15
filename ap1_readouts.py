"""
Generalised readout layer for the v3 AP-1 model.

PROBLEM THIS FIXES
------------------
v3 assumed every scar gene has the same promoter polarity: cFos:cJun activates,
cFos:JunB represses (the polarity Sabin et al. 2019 demonstrated for the GFAP
promoter, Fig 2b). Under that assumption every readout must move in the SAME
direction when the switch flips. But Sabin Fig 6b shows CSPG5 goes DOWN (0.75x)
on miR-200a inhibition while GFAP, vimentin and CSPG4 go UP. The single-polarity
model cannot produce opposite directions -- a structural miss.

THE FIX
-------
Give each readout a polarity s:
  s = +1  cFos:cJun activates, cFos:JunB represses   (GFAP-like, measured)
  s = -1  cFos:JunB activates, cFos:cJun represses   (mirror)
Nothing in the literature requires all promoters to share GFAP's polarity;
AP-1 dimers are known to act oppositely on different promoters.

Supporting observation: in Sabin Fig 6a's heatmap CSPG5 clusters with the
DOWN-regulated group (NEFH, NEFM, CNTN1, CHL1, BCAN...), which is largely
neuronal. CSPG5 (neuroglycan C) is a neuronal/developmental proteoglycan
rather than a classic reactive-scar CSPG, so an opposite polarity -- or
exclusion from the glial scar panel entirely -- is biologically reasonable.

RESULT under the round-2 reference model: THE POLARITY HYPOTHESIS DOES NOT
RESCUE THE READOUTS. The ceiling is dynamic range, not direction.

  largest per-cell ON/OFF ratio any cFos:cJun / cFos:JunB readout can produce
  (phi = 1, i.e. fully AP-1-dependent):   1.246x  (polarity +1)
                                          0.100x  (polarity -1)

  per-cell targets implied by Sabin's bulk ratios and the Fig 5c fractions:
      GFAP 2.14x   vimentin 4.13x   CSPG4 4.71x   CSPG5 -0.04x   Col4a1 1.18x

Only Col4a1 (1.18x) is reachable. GFAP, vimentin and CSPG4 all sit above the
1.246x ceiling, and no choice of phi or polarity reaches them. CSPG5's target
is NEGATIVE, so no readout of any polarity can produce it -- consistent with
dropping CSPG5 from the scar panel. Polarity -1 overshoots in the other
direction (0.100x, a 10-fold decrease) against a measured bulk 0.75.

This strengthens P4 branch (a): scar promoters must integrate c-Jun homodimer
(or otherwise Fos-independent) input, because a c-Fos:c-Jun-only readout tops
out at 1.25x between the two states while the measured per-cell spreads are
2-5x. The reason is the c-Fos ceiling -- FJ is pinned near the limiting c-Fos
pool in both states and simply has no dynamic range to pass on.

WHY THIS IS SAFE
----------------
The readout layer is strictly DOWNSTREAM: no scar gene feeds back into the
AP-1 equations. Verified numerically (verify_downstream below): deliberately
wrecking the GFAP parameters changes the switch species by ~1e-15, i.e.
floating-point noise. Adding, removing or re-polarising readouts therefore
CANNOT alter the bistable window, the state separation, or the ablation result.

Author: Shaunak Patel
Model built with AI assistance (Claude); reviewed and owned by the author.
"""

import numpy as np
import matplotlib.pyplot as plt

from ap1_model_v3 import (SPECIES, N, rhs, hill, repress,
                          steady_state, totals, hysteresis_sweep)
from round2 import ref_params

# Measured bulk fold-changes, Sabin et al. 2019 Fig 6b
# (4 dpi miR-200a inhibitor / 4 dpi control). Read from the published scatter;
# re-check against the figure before using in a write-up.
MEASURED_BULK = {
    'GFAP': 1.23, 'vimentin': 1.55, 'CSPG4': 1.63,
    'CSPG5': 0.75, 'Col4a1': 1.04,
}

# Polarity assignment. +1 is the measured GFAP polarity; -1 is the mirror.
# Only CSPG5 is assigned -1, on the basis of its direction and its clustering
# with neuronal genes in Fig 6a. This is a hypothesis, flagged as such (A17).
POLARITY = {
    'GFAP': +1, 'vimentin': +1, 'CSPG4': +1,
    'CSPG5': -1, 'Col4a1': +1,
}


def readout_level(y, p, phi, polarity=+1, drive=20.0):
    """Steady-state level of a scar gene with AP-1-dependent fraction phi and
    promoter polarity. phi=0 is fully AP-1-independent; phi=1 fully dependent."""
    FJ = y[SPECIES.index('FJ')]
    FB = y[SPECIES.index('FB')]
    if polarity > 0:
        A = hill(FJ, p['Km'], p['n']) * repress(FB, p['Km_rep'], p['n_rep'])
    else:
        A = hill(FB, p['Km_rep'], p['n_rep']) * repress(FJ, p['Km'], p['n'])
    return drive * ((1.0 - phi) + phi * A) / p['g_G']


def window_branch_states(p=None, mirs=None):
    """Full ON- and OFF-branch state vectors at the middle of the bistable
    window, obtained by continuation from each end of the miR axis."""
    p = p or ref_params()
    mirs = mirs if mirs is not None else np.linspace(0, 3, 61)
    up, down = hysteresis_sweep(p, mirs)
    uJ = np.array([t['cJun_total'] for t in up])
    dJ = np.array([t['cJun_total'] for t in down])
    mid = (uJ.max() + uJ.min()) / 2
    win = (dJ < mid) & (uJ > mid)
    if not win.any():
        raise RuntimeError('no bistable window for these parameters')
    idx = np.where(win)[0]
    m_mid = mirs[idx[len(idx) // 2]]

    def branch(start):
        q = dict(p)
        q['mir'] = start
        y = steady_state(q)
        for m in np.linspace(start, m_mid, 7)[1:]:
            q['mir'] = m
            y = steady_state(q, y0=y)
        return y

    return m_mid, branch(0.0), branch(3.0)


def verify_downstream(p=None, tol=1e-9):
    """Prove the readout layer cannot perturb the switch."""
    p = p or ref_params()
    base = steady_state(p)
    q = dict(p)
    q['b_G_by_FJ'] = 999.0
    q['a_G_indep'] = 50.0
    alt = steady_state(q)
    sw = ['F', 'J', 'B', 'FJ', 'FB', 'JJ', 'BB', 'JB']
    diff = max(abs(base[SPECIES.index(s)] - alt[SPECIES.index(s)]) for s in sw)
    print(f'max change in switch species after altering readout params: '
          f'{diff:.2e}  ({"PASS" if diff < tol else "FAIL"})')
    return diff < tol


def fit_phi(target_ratio, polarity, y_on, y_off, p):
    """Solve for the phi reproducing a given per-cell ON/OFF ratio."""
    from scipy.optimize import brentq

    def err(x):
        phi = 1 / (1 + np.exp(-x))
        r = (readout_level(y_on, p, phi, polarity)
             / readout_level(y_off, p, phi, polarity))
        return r - target_ratio
    try:
        x = brentq(err, -12, 12)
        return 1 / (1 + np.exp(-x))
    except ValueError:
        return None


def polarity_demo(cjun_pos_fractions=(0.08, 0.30), plot=True, p=None):
    """Show that with polarity, the model spans BOTH directions, and report the
    phi each measured gene implies. Per-cell targets are obtained from the bulk
    ratios using the Fig 5c cJun+ fractions."""
    p = p or ref_params()
    m_mid, y_on, y_off = window_branch_states(p)
    f0, f1 = cjun_pos_fractions

    def per_cell(bulk):
        return (bulk * (1 - f0) - (1 - f1)) / (f1 - bulk * f0)

    print(f'window midpoint mir = {m_mid:.2f}\n')
    print('gene       bulk   per-cell   polarity   fitted phi   model per-cell')
    rows = []
    for g, bulk in MEASURED_BULK.items():
        tgt = per_cell(bulk)
        s = POLARITY[g]
        phi = fit_phi(tgt, s, y_on, y_off, p)
        if phi is None:
            print(f'{g:9s} {bulk:5.2f}  {tgt:8.2f}   {s:+d}         '
                  f'unreachable   --')
            rows.append((g, tgt, None, None))
            continue
        got = (readout_level(y_on, p, phi, s)
               / readout_level(y_off, p, phi, s))
        print(f'{g:9s} {bulk:5.2f}  {tgt:8.2f}   {s:+d}         '
              f'{phi:.3f}        {got:.2f}x')
        rows.append((g, tgt, phi, got))

    if plot:
        phis = np.linspace(0.01, 0.99, 60)
        up = [readout_level(y_on, p, f, +1) / readout_level(y_off, p, f, +1)
              for f in phis]
        dn = [readout_level(y_on, p, f, -1) / readout_level(y_off, p, f, -1)
              for f in phis]
        plt.figure(figsize=(7, 4.4))
        plt.plot(phis, up, lw=2, label='polarity +1 (up-genes)')
        plt.plot(phis, dn, lw=2, label='polarity -1 (down-genes)')
        plt.axhline(1.0, color='k', lw=1, ls=':')
        for g, tgt, phi, _ in rows:
            if phi is not None:
                plt.scatter([phi], [tgt], zorder=5)
                plt.annotate(g, (phi, tgt), textcoords='offset points',
                             xytext=(6, 4), fontsize=9)
        plt.yscale('log')
        plt.xlabel('AP-1-dependent fraction phi')
        plt.ylabel('per-cell ON/OFF ratio')
        plt.title('Polarity lets one switch drive genes in both directions')
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.show()
    return rows


if __name__ == '__main__':
    verify_downstream()
    print()
    polarity_demo(plot=False)
