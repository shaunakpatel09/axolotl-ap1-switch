"""
regime.py -- the species arm: where does a cell's AP-1 production ratio put it
relative to the two-state window?

THE QUESTION THIS REPLACES
--------------------------
"What dose would a human need?" is the wrong question for this model.
Pharmacokinetics describes what a body does to a drug; it does not make a
salamander network into a mammalian one. The axolotl-to-mammal gap lives in
the network, and the sharpest version of it is:

  Do axolotl glia sit INSIDE the two-state window (switchable) while
  mammalian reactive astrocytes sit OUTSIDE it, locked in the high-c-Jun
  state?

If so, that is a one-sentence account of why axolotls regenerate and mammals
scar, and it changes the drug question. In a locked (monostable) regime there
is no second state to latch into, so an intervention suppresses c-Jun only
while it is present. Durability -- the thing Test 1 found inside the window --
is a property of the regime, not of the drug.

WHAT IS SCANNED
---------------
The global sensitivity analysis named the c-Jun/c-Fos production ratio as the
strongest determinant of bistability (Fig 6F), with g_J and a_J next. So the
plane here is the competition plane, with c-Fos production held at its
reference value:

  x = a_J / a_F    c-Jun production relative to c-Fos
  y = a_B / a_F    JunB production relative to c-Fos

Both are ratios of transcription-factor abundances, which is what expression
data reports -- so a measured mammalian c-JUN : c-FOS : JUNB ratio places a
point on this map directly. The axolotl reference sits at (1.0, 0.5).

NO MAMMALIAN NUMBERS ARE HARD-CODED HERE. place() takes a ratio and reports
the regime; the ratio has to come from data, and which data is a choice the
author makes and defends.

Author: Shaunak Patel. Built with AI assistance (Claude); reviewed and owned
by the author.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from ap1_model_v3 import N, SPECIES, steady_state, totals
from round2 import ref_params, robust_window
from pharmacology import (DRUGS, COLORS, apply_drug, max_engagement,
                          high_state, required_engagement, after_withdrawal)

_J = SPECIES.index('J')
OUT_DIR = 'figures_round2'
RESULT_DIR = 'regime_results'

# coarse grid used for the map; verified against the full-resolution
# robust_window() at selected points by check_resolution()
MAP_MIRS = np.linspace(0.0, 2.5, 26)


def ratio_params(x, y, p=None):
    """Reference model with c-Jun and JunB production set to x and y times
    c-Fos production. Reference = (1.0, 0.5)."""
    p = dict(p or ref_params())
    p['a_J'] = x * p['a_F']
    p['a_B'] = y * p['a_F']
    return p


def window_width(p, mirs=MAP_MIRS, t_end=400.0):
    """Width of the two-state window in miR-200a units, 0 if monostable."""
    w = robust_window(p, mirs=mirs)
    return w['width'], w


def state_at(p, m0, t_end=800.0):
    """Steady state reached from HIGH c-Jun at baseline miR-200a m0. In a
    monostable regime this is the only state; inside the window it is the
    scar state."""
    return totals(high_state(p, m0))


def classify(p, m0, mirs=MAP_MIRS):
    """Regime of a parameter set at baseline miR-200a m0.

      'switchable'  m0 lies inside the two-state window
      'locked'      monostable at m0, and c-Jun stays high
      'off'         monostable at m0, and c-Jun is low
      'bistable elsewhere' window exists but not at m0

    High/low is judged against the midpoint of this parameter set's own c-Jun
    range across the miR axis, so it does not depend on absolute scale."""
    width, w = window_width(p, mirs)
    q = dict(p)
    lo_end = totals(steady_state(dict(q, mir=mirs[-1]), y0=np.zeros(N),
                                 t_end=1500.0))['cJun_total']
    y = np.zeros(N)
    y[_J] = 100.0
    hi_end = totals(steady_state(dict(q, mir=0.0), y0=y,
                                 t_end=1500.0))['cJun_total']
    here = state_at(p, m0)['cJun_total']
    mid = 0.5 * (lo_end + hi_end)
    if width > 0 and w['lo'] <= m0 <= w['hi']:
        return 'switchable', dict(width=width, w=w, cJun=here, mid=mid)
    kind = 'locked' if here > mid else 'off'
    if width > 0:
        kind = 'bistable elsewhere'
    return kind, dict(width=width, w=w, cJun=here, mid=mid)


# ----------------------------------------------------------------------
# the map
# ----------------------------------------------------------------------

def map_batch(xs, ys, m0=0.58, start=0, stop=None, tag='map'):
    """Compute the regime map over the grid, in slices so each call stays
    short. Saves to regime_results/<tag>_<start>_<stop>.npz."""
    pts = [(i, j) for i in range(len(xs)) for j in range(len(ys))]
    stop = len(pts) if stop is None else min(stop, len(pts))
    rows = []
    for k in range(start, stop):
        i, j = pts[k]
        p = ratio_params(xs[i], ys[j])
        width, w = window_width(p)
        s = state_at(p, m0)
        rows.append((i, j, width, s['cJun_total'], s['GFAP'],
                     np.nan if w['lo'] is None else w['lo'],
                     np.nan if w['hi'] is None else w['hi']))
    os.makedirs(RESULT_DIR, exist_ok=True)
    path = os.path.join(RESULT_DIR, f'{tag}_{start}_{stop}.npz')
    np.savez(path, rows=np.array(rows), xs=xs, ys=ys, m0=m0)
    print('saved', path, f'({stop - start} points)', flush=True)
    return path


def load_map(xs, ys, tag='map'):
    rows = []
    for f in sorted(os.listdir(RESULT_DIR)):
        if f.startswith(tag + '_') and f.endswith('.npz'):
            rows.append(np.load(os.path.join(RESULT_DIR, f))['rows'])
    r = np.vstack(rows)
    W = np.full((len(xs), len(ys)), np.nan)
    CJ = np.full((len(xs), len(ys)), np.nan)
    G = np.full((len(xs), len(ys)), np.nan)
    LO = np.full((len(xs), len(ys)), np.nan)
    for i, j, width, cj, g, lo, hi in r:
        W[int(i), int(j)] = width
        CJ[int(i), int(j)] = cj
        G[int(i), int(j)] = g
        LO[int(i), int(j)] = lo
    return W, CJ, G, LO


# ----------------------------------------------------------------------
# what a drug does off the island
# ----------------------------------------------------------------------

def drug_in_regime(p, m0, drug, n=161, drop=0.25):
    """Apply a drug to whatever regime p is in. Returns the engagement needed
    to drop c-Jun (by a jump, as in Test 1) or, failing a jump, the engagement
    needed to halve it; plus whether the cell stays down after withdrawal."""
    y = high_state(p, m0)
    grid = np.linspace(0.0, 0.999 * max_engagement(p, drug, m0), n)
    states = []
    for E in grid:
        y = steady_state(apply_drug(p, drug, E, m0), y0=y, t_end=1500.0)
        states.append(y.copy())
    cj = np.array([totals(s)['cJun_total'] for s in states])
    E_jump, k = required_engagement(grid, states, drop)
    half = np.where(cj < 0.5 * cj[0])[0]
    E_half = grid[half[0]] if len(half) else np.nan
    k_ref = k if k is not None else (half[0] if len(half) else None)
    stays = None
    if k_ref is not None:
        back = totals(after_withdrawal(p, m0, states[k_ref]))
        stays = bool(back['cJun_total'] < 0.75 * cj[0])
    return dict(grid=grid, cJun=cj, E_jump=E_jump, E_half=E_half,
                stays=stays, jumped=k is not None)


def compare_regimes(points, m0=0.58, p=None):
    """Run every Test-1 drug in each named regime and report engagement and
    durability side by side.  points = {label: (x, y)}"""
    out = {}
    for label, (x, y) in points.items():
        q = ratio_params(x, y, p)
        kind, info = classify(q, m0)
        print(f'\n{label}:  a_J/a_F = {x:.2f}, a_B/a_F = {y:.2f}'
              f'   regime at miR-200a {m0:.2f}: {kind.upper()}'
              f'   (window width {info["width"]:.2f})', flush=True)
        print(f'  {"drug":40s}{"jump at":>10s}{"halved at":>11s}{"stays down":>12s}')
        res = {}
        for drug, name in DRUGS.items():
            r = drug_in_regime(q, m0, drug)
            res[drug] = r
            j = '  never' if np.isnan(r['E_jump']) else f'{r["E_jump"]:.2f}'
            h = '  never' if np.isnan(r['E_half']) else f'{r["E_half"]:.2f}'
            s = '-' if r['stays'] is None else ('yes' if r['stays'] else 'NO')
            print(f'  {name:40s}{j:>10s}{h:>11s}{s:>12s}', flush=True)
        out[label] = dict(params=q, kind=kind, info=info, drugs=res)
    return out


def junb_axis(x=1.0, ys=(0.1, 0.35, 0.6, 0.85, 1.1, 1.35, 1.5), m0=0.58,
              koff_JB=None):
    """Window width and GFAP along the JunB production axis. With koff_JB
    lowered, the JunB:c-Jun heterodimer becomes a real sink for c-Jun."""
    W, G = [], []
    for y in ys:
        p = ratio_params(x, y)
        if koff_JB is not None:
            p['koff_JB'] = koff_JB
        W.append(robust_window(p)['width'])
        G.append(state_at(p, m0)['GFAP'])
    return np.array(ys, float), np.array(W), np.array(G)


def withdrawal_trace(p, m0, drug, E, t_on=400.0, t_off=1200.0, n=600):
    """Hold the drug at engagement E until t_on, then stop. Returns (t, cJun)."""
    from ap1_model_v3 import simulate
    y = high_state(p, m0)
    s1 = simulate(apply_drug(p, drug, E, m0), y0=y, t_end=t_on, n_pts=n // 3)
    s2 = simulate(dict(p, mir=m0), y0=s1.y[:, -1], t_end=t_off, n_pts=n)
    t = np.concatenate([s1.t, t_on + s2.t])
    cj = np.concatenate([[totals(s1.y[:, k])['cJun_total']
                          for k in range(s1.y.shape[1])],
                         [totals(s2.y[:, k])['cJun_total']
                          for k in range(s2.y.shape[1])]])
    return t, cj, t_on


def figure9(cmp_res, xs=None, ys=None, m0=0.58, drug='jnk'):
    """Figure 9: the species arm."""
    xs = np.linspace(0.25, 3.0, 16) if xs is None else xs
    ys = np.linspace(0.10, 1.50, 12) if ys is None else ys
    W, CJ, G, LO = load_map(xs, ys)
    j = 3                                    # a_B/a_F column used for panel A
    fig, ax = plt.subplots(1, 4, figsize=(22, 4.8))

    # A -- the window slides away from the cell
    a = ax[0]
    lo, hi = LO[:, j], LO[:, j] + W[:, j]
    a.fill_between(xs, lo, hi, color='purple', alpha=0.25,
                   label='two-state window')
    a.plot(xs, lo, color='purple', lw=1)
    a.plot(xs, hi, color='purple', lw=1)
    a.axhline(m0, color='k', ls='--', lw=1.5,
              label=f"cell's baseline miR-200a ({m0:.2f})")
    inside = (lo <= m0) & (m0 <= hi)
    a.plot(xs[inside], np.full(inside.sum(), m0), 'o', color='#1D9E75', ms=6,
           label='switchable')
    a.plot(xs[~inside], np.full((~inside).sum(), m0), 'o', color='#D85A30',
           ms=6, label='locked high (scar)')
    a.axvline(1.0, color='0.5', lw=1, ls=':')
    a.text(1.05, 0.03, 'axolotl\nreference', fontsize=7, color='0.4',
           va='bottom')
    a.set_xlabel('c-Jun production / c-Fos production')
    a.set_ylabel('miR-200a (multiples of its repression constant)')
    a.set_title('A  The switch does not vanish, it slides.\n'
                'More c-Jun needs more miR-200a to hold it down.',
                loc='left', fontsize=10)
    a.legend(fontsize=7, loc='upper left')

    # B -- JunB tunes the output, not the switch
    b = ax[1]
    yj, Wj, Gj = junb_axis(1.0, tuple(ys), m0)
    _, Wk, _ = junb_axis(1.0, tuple(ys), m0, koff_JB=100.0)
    b.plot(yj, Wj, 'o-', color='#534AB7',
           label='window width (c-Fos competition only)')
    b.plot(yj, Wk, 's--', color='#8E86D8',
           label='window width (JunB:c-Jun dimer 10x stronger)')
    b.set_xlabel('JunB production / c-Fos production')
    b.set_ylabel('two-state window width (miR units)')
    b.set_ylim(0, max(Wk.max(), Wj.max()) * 1.25)
    b2 = b.twinx()
    b2.plot(yj, Gj, '^-', color='#1D9E75', label='GFAP in the scar state')
    b2.set_ylabel('GFAP (nM)', color='#1D9E75')
    b2.tick_params(axis='y', colors='#1D9E75')
    h1, l1 = b.get_legend_handles_labels()
    h2, l2 = b2.get_legend_handles_labels()
    b.legend(h1 + h2, l1 + l2, fontsize=7, loc='lower left')
    b.set_title('B  Competing for c-Fos does not move the switch.\n'
                'Sequestering c-Jun directly does.', loc='left', fontsize=10)

    # C -- engagement needed, both regimes
    c = ax[2]
    labels = list(cmp_res)
    width = 0.35
    pos = np.arange(len(DRUGS))
    for k, lab in enumerate(labels):
        v = [cmp_res[lab]['drugs'][d]['E_jump'] for d in DRUGS]
        v = [1.05 if np.isnan(q) else q for q in v]
        bars = c.bar(pos + (k - 0.5) * width, v, width,
                     color=['#1D9E75', '#D85A30'][k], alpha=0.85,
                     label=lab.split(' (')[0])
        for q, bar in zip(v, bars):
            if q == 1.05:
                c.text(bar.get_x() + bar.get_width() / 2, 1.06, 'never',
                       ha='center', fontsize=7, rotation=90)
    short = {'mimic': 'miR-200a\nmimic', 'jnk': 'JNK\n(activity)',
             'jnk_stab': 'JNK\n(+ stability)', 't5224_fos': 'T-5224\n(c-Fos)',
             't5224_fos_jj': 'T-5224\n(+ homodimer)'}
    c.set_xticks(pos)
    c.set_xticklabels([short[d] for d in DRUGS], fontsize=7.5)
    c.set_ylabel('engagement needed (fraction blocked)')
    c.set_ylim(0, 1.25)
    c.legend(fontsize=7)
    c.set_title('C  The same targets cost roughly twice as much\n'
                'engagement in the locked regime', loc='left', fontsize=10)

    # D -- withdrawal
    d = ax[3]
    for k, lab in enumerate(labels):
        p = ratio_params(*[(1.0, 0.5), (2.0, 0.5)][k])
        E = cmp_res[lab]['drugs'][drug]['E_jump']
        t, cj, t_on = withdrawal_trace(p, m0, drug, min(E * 1.05, 0.99))
        d.plot(t, cj / cj[0], color=['#1D9E75', '#D85A30'][k], lw=2,
               label=lab.split(' (')[0])
    d.axvline(400.0, color='k', ls='--', lw=1)
    d.text(410, 1.02, 'drug stopped', fontsize=8)
    d.set_xlabel('time (h)')
    d.set_ylabel('total c-Jun (relative to untreated scar state)')
    d.set_ylim(0, 1.15)
    d.legend(fontsize=7)
    d.set_title(f'D  {DRUGS[drug]} withdrawn.\n'
                'Only the switchable cell stays down.', loc='left', fontsize=10)

    fig.suptitle('Figure 9. The species arm: whether a drug cures or merely '
                 'suppresses is a property of the regime, not the drug',
                 fontsize=12)
    fig.tight_layout()
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, 'fig9_regime.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print('saved', path)
    return fig


def check_resolution(points=((1.0, 0.5), (2.0, 0.5), (0.6, 0.5))):
    """The map uses a 26-point miR grid at t_end=400 for speed. Confirm it
    agrees with the full 51-point, t_end=800 window used everywhere else."""
    print(f'{"a_J/a_F":>8s}{"a_B/a_F":>9s}{"coarse":>9s}{"full":>8s}')
    for x, y in points:
        p = ratio_params(x, y)
        c = robust_window(p, mirs=MAP_MIRS)['width']
        f = robust_window(p)['width']
        print(f'{x:8.2f}{y:9.2f}{c:9.2f}{f:8.2f}', flush=True)
