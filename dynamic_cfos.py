"""
dynamic_cfos.py -- relaxing the fixed-c-Fos assumption (Degefu, Q5).

In the reference model c-Fos is produced at a constant basal rate. In cells it
is an immediate-early gene: injury induces it transiently (Sabin et al. 2019:
c-Fos and JunB up at 1 dpi, miR-200a peaking at 3 dpi), and c-Fos represses
its own promoter, which terminates the response (Sassone-Corsi, Sisson & Verma
1988, Nature 334:314). Repression requires c-Fos to be dimerised but not to
bind DNA, so here it is driven by dimerised c-Fos (c-Fos:c-Jun + c-Fos:JunB).

Two additions; every other v3 equation is unchanged (rhs_dyn wraps the v3 rhs
and only replaces the c-Fos and JunB production rates):

  1. Negative autoregulation. Production is multiplied by K^n/(K^n + D^n),
     D = dimerised c-Fos. Strength s in [0, 1) is the fractional repression
     at the reference level D_ref, and production is rescaled by 1/(1-s) so
     the reference steady state is unchanged. K and the rescaling are fixed
     at the reference, so changing c-Fos production is resisted by the loop.
  2. Injury pulses. c-Fos production multiplied by 1 + A_F * pulse(t),
     peaking at 24 h; JunB optionally co-induced with the same timing; the
     miR-200a pulse (peak 72 h) comes from the v3 model.

Experiments
  E1  Does regulated c-Fos change the switch, and does it widen the narrow
      range of c-Fos production that supports bistability (Figure 2)?
  E2  A realistic injury, starting in the scar state with baseline miR-200a
      inside the window. Which state does the cell end in?
  E3  How the miR-200a pulse needed for a durable switch depends on when the
      c-Fos pulse peaks.

Author: Shaunak Patel. Built with AI assistance (Claude); reviewed and owned
by the author.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from ap1_model_v3 import N, SPECIES, rhs as rhs_v3
from round2 import ref_params

IDX = {s: SPECIES.index(s) for s in SPECIES}
OUT_DIR = 'figures_round2'


def _pulse(t, peak, sharp=2.0):
    x = max(t, 0.0) / peak
    return (x * np.exp(1.0 - x)) ** sharp


def cfos_production(t, y, p):
    a = p['a_F']
    if p.get('fos_amp', 0.0) > 0:
        a *= 1.0 + p['fos_amp'] * _pulse(t, p.get('fos_peak_t', 24.0))
    s = p.get('fos_autorep', 0.0)
    if s > 0:
        n = p.get('n_fos', 2.0)
        D = max(y[IDX['FJ']] + y[IDX['FB']], 0.0)
        K = p['fos_ref'] / (s / (1.0 - s)) ** (1.0 / n)
        a = a / (1.0 - s) * K ** n / (K ** n + D ** n)
    return a


def junb_production(t, p):
    a = p['a_B']
    if p.get('junb_amp', 0.0) > 0:
        a *= 1.0 + p['junb_amp'] * _pulse(t, p.get('fos_peak_t', 24.0))
    return a


def rhs_dyn(t, y, p):
    q = dict(p)
    q['a_F'] = cfos_production(t, y, p)
    q['a_B'] = junb_production(t, p)
    return rhs_v3(t, y, q)


def integrate(p, y0, t_end, n_pts=None):
    t_eval = None if n_pts is None else np.linspace(0, t_end, n_pts)
    return solve_ivp(rhs_dyn, (0, t_end), y0, args=(p,), method='LSODA',
                     t_eval=t_eval, rtol=1e-8, atol=1e-10)


def steady(p, y0, t_end=800.0):
    return integrate(p, y0, t_end).y[:, -1]


def cjun(y):
    return y[IDX['J']] + y[IDX['FJ']] + 2 * y[IDX['JJ']] + y[IDX['JB']]


def reference_dimerised_cfos(p=None):
    p = p or ref_params()
    q = dict(p)
    q['mir'] = 0.0
    y0 = np.zeros(N)
    y0[IDX['J']] = 100.0
    y = steady(q, y0, 2000.0)
    return y[IDX['FJ']] + y[IDX['FB']]


def branches(p, mirs):
    """Continuation up (from high c-Jun) and down (from zero)."""
    q = dict(p)
    y = np.zeros(N)
    y[IDX['J']] = 100.0
    up = []
    for m in mirs:
        q['mir'] = m
        y = steady(q, y)
        up.append(y.copy())
    q['mir'] = mirs[-1]
    y = steady(q, np.zeros(N))
    down = []
    for m in mirs[::-1]:
        q['mir'] = m
        y = steady(q, y)
        down.append(y.copy())
    return np.array(up), np.array(down[::-1])


def window(p, mirs=None):
    """(width, low edge, high edge, c-Jun ratio) of the two-state window in
    miR-200a; zeros/nan if monostable. Candidates are verified by a further
    5000 h of integration from both branches."""
    mirs = np.linspace(0, 2.5, 51) if mirs is None else mirs
    up, dn = branches(p, mirs)
    cu = np.array([cjun(y) for y in up])
    cd = np.array([cjun(y) for y in dn])
    hyst = np.abs(cu - cd) / (0.5 * (cu + cd)) > 0.25
    if not hyst.any():
        return 0.0, np.nan, np.nan, np.nan
    i = np.where(hyst)[0][hyst.sum() // 2]
    q = dict(p)
    q['mir'] = mirs[i]
    a = steady(q, up[i], 5000.0)
    b = steady(q, dn[i], 5000.0)
    if abs(cjun(a) - cjun(b)) / (0.5 * (cjun(a) + cjun(b))) < 0.1:
        return 0.0, np.nan, np.nan, np.nan
    step = mirs[1] - mirs[0]
    return (hyst.sum() * step, mirs[hyst].min(), mirs[hyst].max(),
            max(cjun(a), cjun(b)) / min(cjun(a), cjun(b)))


# ----------------------------------------------------------------------
# E1 -- regulated vs fixed c-Fos
# ----------------------------------------------------------------------

def e1_regulation(strengths=(0.0, 0.5, 0.8, 0.9),
                  folds=(0.25, 0.375, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0)):
    p0 = ref_params()
    d_ref = reference_dimerised_cfos(p0)
    table = {}
    for s in strengths:
        rows = []
        for f in folds:
            p = dict(p0)
            p.update(a_F=p0['a_F'] * f, fos_autorep=s, fos_ref=d_ref)
            rows.append(window(p)[:3])
        table[s] = rows
        print(f'autorepression {s:.1f}: window [lo-hi] by c-Fos production fold '
              + ' '.join(f'{f:g}:[{lo:.2f}-{hi:.2f}]' for f, (_, lo, hi) in zip(folds, rows)),
              flush=True)
    return folds, table, d_ref


# ----------------------------------------------------------------------
# E2 / E3 -- injury time courses
# ----------------------------------------------------------------------

def start_states(p, m0, n_steps=41):
    """High- and low-c-Jun steady states at baseline miR-200a m0, by fine
    continuation (coarse steps can fall off a narrow high branch). Raises if
    the two are not distinct, so no switching test compares a state with
    itself."""
    q = dict(p)
    y = np.zeros(N)
    y[IDX['J']] = 100.0
    for m in np.linspace(0.0, m0, n_steps):
        q['mir'] = m
        y = steady(q, y)
    hi = y.copy()
    y = np.zeros(N)
    for m in np.linspace(max(3.0, 1.5 * m0), m0, n_steps):
        q['mir'] = m
        y = steady(q, y)
    if cjun(hi) < 1.25 * cjun(y):
        raise ValueError(f'no distinct states at miR-200a {m0:.2f}: '
                         f'c-Jun {cjun(hi):.2f} vs {cjun(y):.2f}')
    return hi, y


def injury(p, m0, y0, A_m, A_F, A_B=0.0, t_fos=24.0, t_mir=72.0,
           t_end=800.0, n_pts=None):
    q = dict(p)
    q.update(mir=m0, mir_pulse=A_m > 0, mir_amp=A_m, mir_peak_t=t_mir,
             mir_sharp=2.0, fos_amp=A_F, fos_peak_t=t_fos, junb_amp=A_B)
    return integrate(q, y0, t_end, n_pts)


def switched(p, m0, A_m, A_F, A_B=0.0, t_fos=24.0, states=None):
    hi, lo = states or start_states(p, m0)
    sol = injury(p, m0, hi, A_m, A_F, A_B, t_fos)
    return cjun(sol.y[:, -1]) < np.sqrt(cjun(hi) * cjun(lo))


def e2_outcome_map(p=None, m0=None, A_F=(0, 1, 2, 4, 8, 16),
                   A_m=(0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0), junb=False):
    p = p or ref_params()
    if m0 is None:
        _, lo_e, hi_e, _ = window(p)
        m0 = 0.5 * (lo_e + hi_e)
    st = start_states(p, m0)
    M = np.zeros((len(A_m), len(A_F)), dtype=int)
    for i, am in enumerate(A_m):
        for j, af in enumerate(A_F):
            M[i, j] = switched(p, m0, am, af, af if junb else 0.0, states=st)
    return m0, np.array(A_F), np.array(A_m), M


def min_mir_pulse(p, m0, A_F, A_B, t_fos, states, hi_cap=3.0, tol=0.02):
    """Smallest miR-200a pulse giving a durable switch (bisection);
    inf if none up to hi_cap."""
    if not switched(p, m0, hi_cap, A_F, A_B, t_fos, states):
        return np.inf
    if switched(p, m0, 0.0, A_F, A_B, t_fos, states):
        return 0.0
    lo, hi = 0.0, hi_cap
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if switched(p, m0, mid, A_F, A_B, t_fos, states):
            hi = mid
        else:
            lo = mid
    return hi


def e3_timing(p=None, m0=None, A_F=4.0,
              t_fos=(0, 12, 24, 48, 72, 96, 120, 144)):
    p = p or ref_params()
    if m0 is None:
        _, lo_e, hi_e, _ = window(p)
        m0 = 0.5 * (lo_e + hi_e)
    st = start_states(p, m0)
    base = min_mir_pulse(p, m0, 0.0, 0.0, 24.0, st)
    out = {}
    for lab, A_B in (('c-Fos alone', 0.0), ('c-Fos + JunB', A_F)):
        out[lab] = [min_mir_pulse(p, m0, A_F, A_B, max(t, 0.5), st)
                    for t in t_fos]
        print(f'{lab:13s} min miR-200a pulse by c-Fos peak time: '
              + ' '.join(f'{t}h:{v:.2f}' for t, v in zip(t_fos, out[lab])),
              flush=True)
    print(f'no c-Fos pulse: {base:.2f}')
    return m0, np.array(t_fos), out, base


# ----------------------------------------------------------------------
# Figure 7
# ----------------------------------------------------------------------

def figure7(e1, e3, thresholds, m0):
    folds, table, _ = e1
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))

    a = ax[0, 0]
    for s, color in ((0.0, '#D85A30'), (0.8, '#1D9E75')):
        lo = np.array([r[1] for r in table[s]], float)
        hi = np.array([r[2] for r in table[s]], float)
        lab = 'c-Fos fixed' if s == 0 else f'c-Fos autorepressed ({s:g})'
        a.fill_between(folds, lo, hi, color=color, alpha=0.18)
        a.plot(folds, hi, 'o-', color=color, label=f'{lab}: upper edge')
    a.set_xscale('log')
    a.set_xlabel('c-Fos production (fold of reference)')
    a.set_ylabel('miR-200a (relative units)')
    a.set_title('A  Autorepression keeps the switch in place as\n'
                'c-Fos production changes (shaded: two-state window)',
                loc='left', fontsize=10)
    a.legend(fontsize=8)

    p = ref_params()
    hi_s, _ = start_states(p, m0)
    b = ax[0, 1]
    for A_F, ls in ((2.0, '-'), (0.5, '--')):
        sol = injury(p, m0, hi_s, 0.0, A_F, n_pts=400, t_end=400.0)
        b.plot(sol.t, [cjun(sol.y[:, i]) for i in range(sol.y.shape[1])],
               ls, color='#534AB7', lw=2, label=f'total c-Jun, c-Fos pulse {A_F:g}')
        if A_F == 2.0:
            b.plot(sol.t, sol.y[IDX['FJ']], '-', color='#D85A30', lw=1.5,
                   label='c-Fos:c-Jun, pulse 2')
            b.plot(sol.t, sol.y[IDX['FB']], '-', color='#1D9E75', lw=1.5,
                   label='c-Fos:JunB, pulse 2')
    b.set_xlabel('time after injury (h)')
    b.set_ylabel('concentration (nM)')
    b.set_title(f'B  A c-Fos pulse alone (peak 24 h, no miR-200a) switches\n'
                f'the cell above a threshold (baseline miR-200a {m0:.2f})',
                loc='left', fontsize=10)
    b.legend(fontsize=7)

    c = ax[1, 0]
    labs = list(thresholds)
    x = np.arange(len(labs))
    for k, (off, color, lab) in enumerate(((-0.2, '#534AB7', 'c-Fos alone'),
                                           (0.2, '#BA7517', 'c-Fos + JunB'))):
        vals = [thresholds[l][1 + k] for l in labs]
        shown = [v if np.isfinite(v) else 0 for v in vals]
        c.bar(x + off, shown, width=0.38, color=color, label=lab)
        for xi, v in zip(x, vals):
            c.text(xi + off, (v if np.isfinite(v) else 0) + 0.05,
                   f'{v:.2f}' if np.isfinite(v) else 'never', ha='center',
                   fontsize=8)
    c.set_xticks(x)
    c.set_xticklabels([l.replace(' (', '\n(') for l in labs], fontsize=8)
    c.set_ylabel('smallest c-Fos pulse that switches the cell')
    c.set_title('C  A stronger c-Fos:c-Jun -> c-Jun edge blocks\n'
                'c-Fos-triggered switching', loc='left', fontsize=10)
    c.legend(fontsize=8)

    _, t_fos, out, base = e3
    d = ax[1, 1]
    d.plot(t_fos, out['c-Fos alone'], 'o-', color='#534AB7',
           label='with c-Fos pulse (half its switching threshold)')
    d.axhline(base, color='#888780', ls='--', label='no c-Fos pulse')
    d.axvline(72, color='#1D9E75', ls=':', label='miR-200a peak (72 h)')
    d.axvline(24, color='#D85A30', ls=':', label='c-Fos peak in Sabin (24 h)')
    d.set_xlabel('time of c-Fos peak (h after injury)')
    d.set_ylabel('smallest miR-200a pulse for a lasting switch')
    d.set_title('D  c-Fos and miR-200a add up only when they overlap',
                loc='left', fontsize=10)
    d.legend(fontsize=7)

    fig.suptitle('Figure 7. Relaxing the fixed-c-Fos assumption: regulated '
                 'and injury-induced c-Fos', fontsize=12)
    fig.tight_layout()
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, 'fig7_dynamic_cfos.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print('saved', path)
    return fig


# ----------------------------------------------------------------------
# Switching thresholds across FJ-edge strengths (Figure 7, panel C)
# ----------------------------------------------------------------------

THRESHOLD_CASES = (
    ('no FJ edge (basal 16)', dict(b_J_by_FJ=0.0, a_J=16.0)),
    ('reference (FJ 9, basal 8)', dict(b_J_by_FJ=9.0, a_J=8.0)),
    ('strong FJ edge (18, basal 4)', dict(b_J_by_FJ=18.0, a_J=4.0)),
    ('FJ = homodimer (36, basal 1)', dict(b_J_by_FJ=36.0, a_J=1.0)),
)


def _threshold(p, m0, states, which, junb=False, cap=32.0, tol=0.02):
    """Smallest pulse that switches the cell: which='F' for a c-Fos pulse
    (miR-200a off), which='m' for a miR-200a pulse (c-Fos off). inf if none
    up to cap; 0 if the cell switches with no pulse at all."""
    def f(x):
        if which == 'F':
            return switched(p, m0, 0.0, x, x if junb else 0.0, states=states)
        return switched(p, m0, x, 0.0, 0.0, states=states)
    if f(0.0):
        return 0.0
    if not f(cap):
        return np.inf
    lo, hi = 0.0, cap
    while hi - lo > tol * max(hi, 0.1):
        mid = 0.5 * (lo + hi)
        if f(mid):
            hi = mid
        else:
            lo = mid
    return hi


def thresholds(cases=THRESHOLD_CASES):
    """For each strength of the c-Fos:c-Jun -> c-Jun edge (v3 parameters
    otherwise): smallest c-Fos pulse alone, with JunB co-induced, and smallest
    miR-200a pulse that switch a cell out of the high-c-Jun state. Cases with
    no baseline holding two distinct states are reported and skipped."""
    from ap1_model_v3 import v3_params
    out = {}
    for lab, upd in cases:
        p = v3_params()
        p.update(upd)
        w = window(p)
        if w[0] == 0:
            print(f'{lab}: not bistable; skipped')
            continue
        st, m0 = None, None
        for m in (0.5 * (w[1] + w[2]), w[1] + 0.25 * (w[2] - w[1]),
                  w[1] + 0.75 * (w[2] - w[1])):
            try:
                st, m0 = start_states(p, m), m
                break
            except ValueError:
                pass
        if st is None:
            print(f'{lab}: no baseline with two distinct states; skipped')
            continue
        r = (m0, _threshold(p, m0, st, 'F'),
             _threshold(p, m0, st, 'F', junb=True), _threshold(p, m0, st, 'm'))
        out[lab] = r
        print(f'{lab:30s} baseline miR {m0:.2f} | c-Fos alone {r[1]:.2f}, '
              f'with JunB {r[2]:.2f} | miR alone {r[3]:.2f}', flush=True)
    return out


def make_all():
    """Run E1, E2, the thresholds and E3, then draw Figure 7."""
    e1 = e1_regulation(strengths=(0.0, 0.5, 0.8))
    m0, AF, AM, M1 = e2_outcome_map(junb=False)
    _, _, _, M2 = e2_outcome_map(m0=m0, junb=True)
    for lab, M in (('c-Fos alone', M1), ('c-Fos + JunB', M2)):
        print(f'E2, {lab}: rows = miR-200a pulse, columns = c-Fos pulse '
              f'{list(AF)}; S = switched')
        for am, row in zip(AM, M):
            print(f'  {am:4.2f}: ' + ' '.join('S' if v else '.' for v in row))
    th = thresholds()
    e3 = e3_timing(m0=m0, A_F=0.5 * th['reference (FJ 9, basal 8)'][1])
    figure7(e1, e3, th, m0)
    return e1, (m0, AF, AM, M1, M2), th, e3
