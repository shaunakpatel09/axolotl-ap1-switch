"""
gsa.py -- global sensitivity analysis of the round-2 reference AP-1 model.

Method: Latin hypercube sampling with partial rank correlation coefficients
(LHS-PRCC; Marino et al. 2008, J Theor Biol 254:178-196), the approach
Degefu et al. cite for their own parameter sampling.

Eighteen parameters are sampled log-uniformly over a 16-fold range around
the reference (0.25x to 4x) unless noted in SPEC. For each sample:

  width_full  width of the two-state window in miR-200a, in units of the
              miR-200a repression constant Km_mir (0 if monostable)
  width_noJJ  the same with the c-Jun homodimer loop removed, FJ edge kept
  sep_cjun    c-Jun ratio between the two states (bistable samples only)
  sep_gfap    GFAP ratio between the two states (bistable samples only)
  amp_share   with ALL autoregulation removed: effective Hill coefficient of
              the share of c-Fos paired with c-Jun, divided by the miR-200a
              input's Hill coefficient. ~1 means competition adds no steepness.
  amp_freeJ   the same for free c-Jun

Round-2 conclusions tested:
  C1  bistability requires the homodimer loop       -> width_noJJ
  C2  competition alone is graded, not steep        -> amp_share, amp_freeJ
  C3  states differ more in c-Jun than in GFAP      -> sep_cjun vs sep_gfap
  C4  c-Fos abundance must be intermediate          -> binned by a_F
  C5  the FJ edge is tolerated only up to a point   -> binned by r_FJ

PRCC assumes monotone relationships, so band-shaped effects (c-Fos abundance,
basal drive) are also shown as binned bistable fractions. b_G_by_FJ and Km_rep
act only on the GFAP readout, which cannot feed back on the switch; they are
negative controls and should show PRCC ~ 0 on window width.

A candidate window is accepted only if the two branch states stay distinct
after a further 5000 h of integration, to exclude slow convergence near fold
points masquerading as hysteresis.

Run in batches so no single call runs too long; the design is fixed by SEED,
so batches combine into one Latin hypercube:
    import gsa
    gsa.run(0, 200); gsa.run(200, 400)
    res = gsa.load(); gsa.report(res); gsa.figure(res)

Author: Shaunak Patel. Built with AI assistance (Claude); reviewed and owned
by the author.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rankdata, t as tdist

from ap1_model_v3 import N, SPECIES, steady_state, lhs
from round2 import ref_params

SEED = 2026
N_SAMPLES = 400
RESULT_DIR = 'gsa_results'
OUT_DIR = 'figures_round2'

# name, kind, low, high. 'fold' = log-uniform multiplier on the reference;
# 'log' = log-uniform absolute value; 'lin' = uniform absolute value.
SPEC = [
    ('a_F', 'fold', 0.25, 4), ('a_J', 'fold', 0.25, 4),
    ('a_B', 'fold', 0.25, 4), ('b_J_by_JJ', 'fold', 0.25, 4),
    ('r_FJ', 'lin', 0.0, 1.0),
    ('Km_auto', 'fold', 0.25, 4), ('n_auto', 'lin', 1.5, 5.0),
    ('g_F', 'fold', 0.25, 4), ('g_J', 'fold', 0.25, 4),
    ('g_B', 'fold', 0.25, 4),
    ('koff_FJ', 'fold', 0.25, 4), ('koff_JJ', 'fold', 0.25, 4),
    ('junb_weakness', 'log', 1.0, 6.0),
    ('Km_mir', 'fold', 0.25, 4), ('n_mir', 'lin', 1.0, 4.0),
    ('mir_leak', 'lin', 0.02, 0.4),
    ('b_G_by_FJ', 'fold', 0.25, 4), ('Km_rep', 'fold', 0.25, 4),
]
NAMES = [s[0] for s in SPEC]
NEGATIVE_CONTROLS = ('b_G_by_FJ', 'Km_rep')
IDX = {s: SPECIES.index(s) for s in SPECIES}


def design():
    """The full Latin hypercube, as values of each SPEC entry."""
    u = lhs(N_SAMPLES, len(SPEC), np.random.default_rng(SEED))
    X = np.empty_like(u)
    for j, (_, kind, lo, hi) in enumerate(SPEC):
        if kind == 'lin':
            X[:, j] = lo + u[:, j] * (hi - lo)
        else:
            X[:, j] = 10 ** (np.log10(lo) + u[:, j] * (np.log10(hi) - np.log10(lo)))
    return X


def params_from_row(x):
    ref = ref_params()
    p = dict(ref)
    for (name, kind, _, _), v in zip(SPEC, x):
        if name == 'r_FJ':
            continue
        p[name] = ref[name] * v if kind == 'fold' else v
    p['koff_BB'] = p['koff_JB'] = p['koff_JJ']
    p['b_J_by_FJ'] = x[NAMES.index('r_FJ')] * p['b_J_by_JJ']
    return p


def _cjun(Y):
    return Y[:, IDX['J']] + Y[:, IDX['FJ']] + 2 * Y[:, IDX['JJ']] + Y[:, IDX['JB']]


def _sweeps(p, mirs, t_end=800.0):
    """Continuation up (starting high in c-Jun) and down (starting at zero)."""
    q = dict(p)
    y = np.zeros(N)
    y[IDX['J']] = 100.0
    up = []
    for m in mirs:
        q['mir'] = m
        y = steady_state(q, y0=y, t_end=t_end)
        up.append(y.copy())
    q['mir'] = mirs[-1]
    y = steady_state(q, y0=np.zeros(N), t_end=t_end)
    down = []
    for m in mirs[::-1]:
        q['mir'] = m
        y = steady_state(q, y0=y, t_end=t_end)
        down.append(y.copy())
    return np.array(up), np.array(down[::-1])


def window(p, n_pts=81, span=10.0):
    """(width in Km_mir units, c-Jun ratio, GFAP ratio) of the two-state
    window; (0, nan, nan) if monostable."""
    rel_m = np.linspace(0, span, n_pts)
    mirs = rel_m * p['Km_mir']
    up, dn = _sweeps(p, mirs)
    cu, cd = _cjun(up), _cjun(dn)
    hyst = np.abs(cu - cd) / (0.5 * (cu + cd)) > 0.25
    if not hyst.any():
        return 0.0, np.nan, np.nan
    i = np.where(hyst)[0][hyst.sum() // 2]
    q = dict(p)
    q['mir'] = mirs[i]
    yu = steady_state(q, y0=up[i], t_end=5000.0)
    yd = steady_state(q, y0=dn[i], t_end=5000.0)
    ju, jd = _cjun(yu[None])[0], _cjun(yd[None])[0]
    if abs(ju - jd) / (0.5 * (ju + jd)) < 0.1:
        return 0.0, np.nan, np.nan       # slow convergence, not bistability
    hi, lo = (yu, yd) if ju > jd else (yd, yu)
    return (hyst.sum() * (rel_m[1] - rel_m[0]),
            _cjun(hi[None])[0] / _cjun(lo[None])[0],
            hi[IDX['G']] / lo[IDX['G']])


def effective_hill(m, y):
    y = np.asarray(y, float)
    if y.max() - y.min() < 1e-6 * max(1.0, abs(y.max())):
        return np.nan
    r = (y - y.min()) / (y.max() - y.min())
    if r[0] < r[-1]:
        r = 1 - r
    m90 = np.interp(0.9, r[::-1], m[::-1])
    m10 = np.interp(0.1, r[::-1], m[::-1])
    if m90 <= 0 or m10 <= m90:
        return np.nan
    return np.log(81) / np.log(m10 / m90)


def competition_amplification(p):
    q = dict(p)
    q.update(b_J_by_JJ=0.0, b_J_by_FJ=0.0)
    mirs = np.logspace(-2, 1.7, 50) * p['Km_mir']
    up, _ = _sweeps(q, mirs)
    share = up[:, IDX['FJ']] / (up[:, IDX['FJ']] + up[:, IDX['FB']])
    return (effective_hill(mirs, share) / p['n_mir'],
            effective_hill(mirs, up[:, IDX['J']]) / p['n_mir'])


def run(start=0, stop=N_SAMPLES, verbose=True):
    """Evaluate samples start..stop-1 of the fixed design and save them."""
    X = design()
    os.makedirs(RESULT_DIR, exist_ok=True)
    rows = []
    for k in range(start, stop):
        p = params_from_row(X[k])
        w, sj, sg = window(p)
        q = dict(p)
        q['b_J_by_JJ'] = 0.0
        w0, _, _ = window(q, n_pts=61)
        a_sh, a_fj = competition_amplification(p)
        rows.append([w, w0, sj, sg, a_sh, a_fj])
        if verbose and (k + 1 - start) % 25 == 0:
            done = np.array(rows)
            print(f'  {k + 1}/{stop}: {int((done[:, 0] > 0).sum())} bistable '
                  f'so far in this batch', flush=True)
    np.savez(os.path.join(RESULT_DIR, f'batch_{start}_{stop}.npz'),
             idx=np.arange(start, stop), out=np.array(rows))


def load():
    X = design()
    idx, out = [], []
    for f in sorted(os.listdir(RESULT_DIR)):
        if f.endswith('.npz'):
            d = np.load(os.path.join(RESULT_DIR, f))
            idx.append(d['idx'])
            out.append(d['out'])
    idx = np.concatenate(idx)
    out = np.concatenate(out)
    cols = ['width_full', 'width_noJJ', 'sep_cjun', 'sep_gfap',
            'amp_share', 'amp_freeJ']
    res = {c: out[:, i] for i, c in enumerate(cols)}
    res['X'] = X[idx]
    return res


def prcc(X, y):
    """Partial rank correlation of y with each column of X, controlling for
    all other columns. Returns arrays (rho, p)."""
    n, k = X.shape
    R = np.column_stack([rankdata(X[:, j]) for j in range(k)])
    ry = rankdata(y)
    rho, pv = np.empty(k), np.empty(k)
    df = n - 2 - (k - 1)
    for j in range(k):
        A = np.column_stack([np.ones(n), np.delete(R, j, axis=1)])
        rx = R[:, j] - A @ np.linalg.lstsq(A, R[:, j], rcond=None)[0]
        rr = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
        rho[j] = np.corrcoef(rx, rr)[0, 1]
        t = rho[j] * np.sqrt(df / (1 - rho[j] ** 2))
        pv[j] = 2 * tdist.sf(abs(t), df)
    return rho, pv


def report(res):
    X, w = res['X'], res['width_full']
    bist = w > 0
    n = len(w)
    print(f'samples: {n}   bistable: {bist.sum()} ({100 * bist.mean():.1f}%)\n')

    rho, pv = prcc(X, w)
    alpha = 0.05 / len(NAMES)
    print('PRCC on window width (Bonferroni alpha = %.4f)' % alpha)
    for j in np.argsort(-np.abs(rho)):
        tag = ' *' if pv[j] < alpha else ''
        ctrl = '  [negative control]' if NAMES[j] in NEGATIVE_CONTROLS else ''
        print(f'  {NAMES[j]:14s} {rho[j]:+.3f}  p={pv[j]:.1e}{tag}{ctrl}')

    print('\nC1  bistable without the homodimer loop: '
          f'{int((res["width_noJJ"] > 0).sum())} of {n} samples')
    for key, lab in (('amp_share', 'share of c-Fos paired with c-Jun'),
                     ('amp_freeJ', 'free c-Jun')):
        a = res[key][np.isfinite(res[key])]
        print(f'C2  competition-only steepness / input, {lab}: median '
              f'{np.median(a):.2f}, 95th pct {np.percentile(a, 95):.2f}, '
              f'max {a.max():.2f}; >1.5 in {100 * (a > 1.5).mean():.1f}%')
    sj, sg = res['sep_cjun'][bist], res['sep_gfap'][bist]
    print(f'C3  bistable samples with c-Jun ratio > GFAP ratio: '
          f'{int((sj > sg).sum())} of {bist.sum()}; median c-Jun '
          f'{np.median(sj):.2f}x, GFAP {np.median(sg):.2f}x')
    for name in ('a_F', 'a_J', 'r_FJ'):
        v = X[:, NAMES.index(name)]
        edges = np.quantile(v, np.linspace(0, 1, 6))
        fr = [bist[(v >= lo) & (v <= hi)].mean()
              for lo, hi in zip(edges[:-1], edges[1:])]
        print(f'C4/5 bistable fraction across {name} quintiles: '
              + ' '.join(f'{100 * f:.0f}%' for f in fr))
    return rho, pv


def figure(res):
    X, w = res['X'], res['width_full']
    bist = w > 0
    rho, pv = prcc(X, w)
    alpha = 0.05 / len(NAMES)

    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    o = np.argsort(rho)
    colors = ['#888780' if NAMES[j] in NEGATIVE_CONTROLS
              else ('#534AB7' if pv[j] < alpha else '#C9C7EE') for j in o]
    ax[0, 0].barh([NAMES[j] for j in o], rho[o], color=colors)
    ax[0, 0].axvline(0, color='k', lw=1)
    ax[0, 0].set_xlabel('PRCC with width of two-state window')
    ax[0, 0].set_title('A  Which parameters set bistability\n'
                       '(dark: significant; grey: negative controls)',
                       loc='left', fontsize=10)

    for name, color, lab in (('a_F', '#BA7517', 'c-Fos production'),
                             ('a_J', '#D85A30', 'basal c-Jun production')):
        v = X[:, NAMES.index(name)]
        edges = np.quantile(v, np.linspace(0, 1, 7))
        mid = np.sqrt(edges[:-1] * edges[1:])
        fr = np.array([bist[(v >= lo) & (v <= hi)].mean()
                       for lo, hi in zip(edges[:-1], edges[1:])])
        ax[0, 1].plot(mid, 100 * fr, 'o-', color=color, label=lab)
    ax[0, 1].set_xscale('log')
    ax[0, 1].set_xlabel('fold change from reference')
    ax[0, 1].set_ylabel('% of samples bistable')
    ax[0, 1].set_title('B  Abundance effects (C4): scarce c-Fos lowers\nbistability; abundant c-Fos does not prevent it', loc='left', fontsize=10)
    ax[0, 1].legend(fontsize=8)

    v = X[:, NAMES.index('r_FJ')]
    edges = np.linspace(0, 1, 6)
    fr = [bist[(v >= lo) & (v <= hi)].mean()
          for lo, hi in zip(edges[:-1], edges[1:])]
    ax[0, 2].bar(0.5 * (edges[:-1] + edges[1:]), 100 * np.array(fr),
                 width=0.17, color='#534AB7')
    ax[0, 2].set_xlabel('FJ edge strength / homodimer loop strength')
    ax[0, 2].set_ylabel('% of samples bistable')
    ax[0, 2].set_title('C  Tolerance of the FJ edge (C5)', loc='left',
                       fontsize=10)

    for key, color, lab in (('amp_share', '#1D9E75', 'share of c-Fos with c-Jun'),
                            ('amp_freeJ', '#D85A30', 'free c-Jun')):
        a = res[key][np.isfinite(res[key])]
        ax[1, 0].hist(a, bins=30, alpha=0.6, color=color, label=lab)
    ax[1, 0].axvline(1, color='k', ls=':', lw=1)
    ax[1, 0].set_xlabel('steepness with competition only / steepness of '
                        'miR-200a input')
    ax[1, 0].set_ylabel('samples')
    ax[1, 0].set_title('D  Does competition alone sharpen the response? (C2)',
                       loc='left', fontsize=10)
    ax[1, 0].legend(fontsize=8)

    sj, sg = res['sep_cjun'][bist], res['sep_gfap'][bist]
    ax[1, 1].scatter(sj, sg, s=18, color='#534AB7')
    lim = [1, max(sj.max(), sg.max()) * 1.1]
    ax[1, 1].plot(lim, lim, 'k:', lw=1, label='equal')
    ax[1, 1].set_xscale('log')
    ax[1, 1].set_yscale('log')
    ax[1, 1].set_xlabel('c-Jun, high / low state')
    ax[1, 1].set_ylabel('GFAP, high / low state')
    ax[1, 1].set_title('E  Bistable samples: c-Jun vs GFAP separation (C3)',
                       loc='left', fontsize=10)
    ax[1, 1].legend(fontsize=8)

    ref = ref_params()
    ratio = (ref['a_J'] * X[:, NAMES.index('a_J')]) / (ref['a_F'] * X[:, NAMES.index('a_F')])
    edges = np.quantile(ratio, np.linspace(0, 1, 7))
    mid = np.sqrt(edges[:-1] * edges[1:])
    for flag, color, lab in ((bist, '#1D9E75', 'full model'),
                             (res['width_noJJ'] > 0, '#D85A30',
                              'homodimer loop removed')):
        fr = [flag[(ratio >= lo) & (ratio <= hi)].mean()
              for lo, hi in zip(edges[:-1], edges[1:])]
        ax[1, 2].plot(mid, 100 * np.array(fr), 'o-', color=color, label=lab)
    ax[1, 2].axvline(1, color='k', ls=':', lw=1)
    ax[1, 2].set_xscale('log')
    ax[1, 2].set_xlabel('basal c-Jun production / c-Fos production')
    ax[1, 2].set_ylabel('% of samples bistable')
    ax[1, 2].set_title('F  Without the homodimer loop, a switch needs\n'
                       'excess c-Fos (C1)', loc='left', fontsize=10)
    ax[1, 2].legend(fontsize=8)

    fig.suptitle(f'Figure 6. Global sensitivity analysis: {len(w)} Latin '
                 'hypercube samples, 18 parameters, 16-fold ranges',
                 fontsize=12)
    fig.tight_layout()
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, 'fig6_global_sensitivity.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print('saved', path)
    return fig
