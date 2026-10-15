"""
round2.py -- analyses requested by Yoni Degefu, following his suggested
sequence for the question:

  Can competition between c-Jun and JunB for limiting c-Fos, together with
  c-Jun positive feedback, generate two stable AP-1 states with distinct GFAP
  output?

  Figure 1  Steps 1-3  Competition alone vs. with c-Jun autoregulation
  Figure 2  Step 4     How c-Fos abundance changes the switch
  Figure 3  Step 5/Q2  Redistribution of c-Fos between c-Jun and JunB
  Figure 4  FJ edge    c-Fos:c-Jun -> c-Jun feedback restored
  Figure 5  Step 6     Are the two states distinct in GFAP?

REFERENCE MODEL. ref_params() is v3 with the cFos:cJun -> cJun edge restored,
because c-Jun autoregulation through Jun/AP-1 is established (Angel et al.
1988, Cell). Its strength (b_J_by_FJ = 9, a quarter of the homodimer loop) is
a CHOICE, not a measurement. Basal c-Jun is lowered from 16 to 8 so total
basal drive is unchanged. Figure 4 shows the trade-off these sit on.

Figures go to ./figures_round2/.
Author: Shaunak Patel. Built with AI assistance (Claude); reviewed and owned
by the author.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from ap1_model_v3 import (N, SPECIES, v3_params, simulate, steady_state,
                          totals, hill, repress)

_J = SPECIES.index('J')

OUT_DIR = 'figures_round2'


def ref_params():
    p = v3_params()
    p.update(b_J_by_FJ=9.0, a_J=8.0)
    return p


def competition_only(p):
    q = dict(p)
    q.update(b_J_by_JJ=0.0, b_J_by_FJ=0.0)
    return q


def branches(p, mirs, t_end=800.0):
    """Continuation sweeps in miR-200a. The rising sweep starts from HIGH
    c-Jun, so it finds the high branch even when the two-state window reaches
    miR-200a = 0. (The first version of this script started from zero c-Jun,
    which lands on the low branch and missed those windows entirely.)"""
    q = dict(p)
    y = np.zeros(N)
    y[_J] = 100.0
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
    return [totals(v) for v in up], [totals(v) for v in down[::-1]], up, down[::-1]


def robust_window(p, mirs=None):
    """Two-state window: miR-200a values where the branches differ by >25% in
    total c-Jun, verified by a further 5000 h of integration from both
    branches. Returns dict(width, lo, hi, m_mid, on, off)."""
    mirs = np.linspace(0, 2.5, 51) if mirs is None else mirs
    tu, td, yu, yd = branches(p, mirs)
    cu = np.array([t['cJun_total'] for t in tu])
    cd = np.array([t['cJun_total'] for t in td])
    hyst = np.abs(cu - cd) / (0.5 * (cu + cd)) > 0.25
    none = dict(width=0.0, lo=None, hi=None, m_mid=None, on=None, off=None)
    if not hyst.any():
        return none
    i = np.where(hyst)[0][hyst.sum() // 2]
    q = dict(p)
    q['mir'] = mirs[i]
    a = totals(steady_state(q, y0=yu[i], t_end=5000.0))
    b = totals(steady_state(q, y0=yd[i], t_end=5000.0))
    if abs(a['cJun_total'] - b['cJun_total']) / (0.5 * (a['cJun_total'] + b['cJun_total'])) < 0.1:
        return none
    on, off = (a, b) if a['cJun_total'] > b['cJun_total'] else (b, a)
    return dict(width=hyst.sum() * (mirs[1] - mirs[0]), lo=mirs[hyst].min(),
                hi=mirs[hyst].max(), m_mid=mirs[i], on=on, off=off)


def _save(fig, name):
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, name)
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print('saved', path)


def _share(t):
    """Share of c-Fos-containing dimers that are paired with c-Jun."""
    return t['FJ'] / (t['FJ'] + t['FB'])


def _window(mirs, up, down):
    cu = np.array([t['cJun_total'] for t in up])
    cd = np.array([t['cJun_total'] for t in down])
    hyst = np.abs(cu - cd) / (0.5 * (cu + cd)) > 0.25
    return (mirs[hyst].min(), mirs[hyst].max()) if hyst.any() else None


def effective_hill(mirs, y):
    """Effective Hill coefficient of a monotone response:
    n = ln(81) / ln(m10 / m90), where the response has fallen 10% at m90 and
    90% at m10. A pure Hill repression with coefficient n returns n."""
    y = np.asarray(y, float)
    r = (y - y.min()) / (y.max() - y.min())
    if r[0] < r[-1]:
        r = 1 - r
    m90 = np.interp(0.9, r[::-1], mirs[::-1])
    m10 = np.interp(0.1, r[::-1], mirs[::-1])
    return np.log(81) / np.log(m10 / m90)


# ----------------------------------------------------------------------
# Figure 1 -- steps 1-3
# ----------------------------------------------------------------------

def fig1_competition_vs_feedback(p=None):
    p = p or ref_params()
    mirs = np.linspace(0, 2.0, 81)

    # step 1: steepness of competition alone, on a log-spaced axis
    logm = np.logspace(-2, 1.3, 140)
    up_c, _, _, _ = branches(competition_only(p), logm)
    n_cjun = effective_hill(logm, [t['cJun_total'] for t in up_c])
    n_share = effective_hill(logm, [_share(t) for t in up_c])
    print(f'competition only: effective Hill coefficient {n_cjun:.2f} in total '
          f'c-Jun, {n_share:.2f} in the share of c-Fos paired with c-Jun '
          f'(miR-200a repression alone has n = {p["n_mir"]:g})')

    cases = (('A  Competition only', competition_only(p)),
             ('B  Competition + c-Jun autoregulation', p))
    rows = ((lambda t: t['cJun_total'], 'total c-Jun (nM)'),
            (_share, 'share of c-Fos paired with c-Jun'))
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True)
    for col, (title, q) in enumerate(cases):
        up, down, _, _ = branches(q, mirs)
        bounds = _window(mirs, up, down)
        for row, (fn, ylab) in enumerate(rows):
            a = ax[row, col]
            a.plot(mirs, [fn(t) for t in up], '-', lw=2,
                   label='miR-200a rising')
            a.plot(mirs, [fn(t) for t in down], '--', lw=2,
                   label='miR-200a falling')
            if bounds:
                a.axvspan(*bounds, color='purple', alpha=0.12,
                          label='two stable states')
            if col == 0:
                a.set_ylabel(ylab)
            if row == 1:
                a.set_xlabel('miR-200a (relative units)')
        ax[0, col].set_title(title, loc='left', fontsize=10)
        print(f'{title}: window = {bounds}')
    ax[1, 0].text(0.97, 0.9, f'graded and single-valued\n'
                  f'effective Hill n = {n_share:.1f}\n'
                  f'(miR-200a input alone: n = {p["n_mir"]:g})',
                  transform=ax[1, 0].transAxes, ha='right', va='top',
                  fontsize=9)
    ax[0, 1].legend(fontsize=8)
    fig.suptitle('Figure 1. Competition alone gives a graded, single-valued '
                 'response; c-Jun autoregulation adds two stable states',
                 fontsize=11)
    fig.tight_layout()
    _save(fig, 'fig1_competition_vs_feedback.png')
    return fig


# ----------------------------------------------------------------------
# Figure 2 -- step 4
# ----------------------------------------------------------------------

def fig2_cfos_abundance(p=None, levels=(2.0, 8.0, 16.0),
                        scan=(1, 2, 3, 4, 6, 8, 12, 16, 24, 32)):
    p = p or ref_params()
    mirs = np.linspace(0, 2.0, 81)
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
    names = ('low', 'medium', 'high')
    for a, aF, name in zip(ax.flat[:3], levels, names):
        q = dict(p)
        q['a_F'] = aF
        up, down, _, _ = branches(q, mirs)
        bounds = _window(mirs, up, down)
        a.plot(mirs, [t['cJun_total'] for t in up], '-', lw=2,
               label='miR-200a rising')
        a.plot(mirs, [t['cJun_total'] for t in down], '--', lw=2,
               label='miR-200a falling')
        if bounds:
            a.axvspan(*bounds, color='purple', alpha=0.12)
        a.set_title(f'{name} c-Fos (production {aF:g} nM/h)', loc='left',
                    fontsize=10)
        a.set_xlabel('miR-200a (relative units)')
        a.set_ylabel('total c-Jun (nM)')
    ax[0, 0].legend(fontsize=8)

    widths, lows, highs = [], [], []
    for aF in scan:
        q = dict(p)
        q['a_F'] = aF
        w = robust_window(q, np.linspace(0, 2.5, 51))
        widths.append(w['width'])
        lows.append(np.nan if w['lo'] is None else w['lo'])
        highs.append(np.nan if w['hi'] is None else w['hi'])
        print(f'c-Fos production {aF:>4g}: window {lows[-1]:.2f} to '
              f'{highs[-1]:.2f} (width {widths[-1]:.2f})')
    ax[1, 1].fill_between(scan, lows, highs, color='purple', alpha=0.2,
                          label='two stable states')
    ax[1, 1].plot(scan, highs, 'o-', color='#534AB7', label='upper edge')
    ax[1, 1].plot(scan, lows, 's-', color='#888780', label='lower edge')
    ax[1, 1].set_xscale('log')
    ax[1, 1].set_xlabel('c-Fos production (nM/h)')
    ax[1, 1].set_ylabel('miR-200a (relative units)')
    ax[1, 1].set_title('D  Where the two-state window sits', loc='left',
                       fontsize=10)
    ax[1, 1].legend(fontsize=8)
    fig.suptitle('Figure 2. c-Fos abundance shifts where the switch sits in '
                 'miR-200a; very high c-Fos removes it',
                 fontsize=11)
    fig.tight_layout()
    _save(fig, 'fig2_cfos_abundance.png')
    return fig, dict(zip(scan, widths))


# ----------------------------------------------------------------------
# Figure 3 -- step 5 / Q2
# ----------------------------------------------------------------------

def fig3_dimer_allocation(p=None, pulse_amp=1.0):
    p = p or ref_params()
    mirs = np.linspace(0, 2.0, 81)
    up, down, _, _ = branches(p, mirs)
    bounds = _window(mirs, up, down)
    m_mid = 0.5 * (bounds[0] + bounds[1])

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    a = ax[0]
    for key, color, lab in (('FJ', '#D85A30', 'c-Fos:c-Jun'),
                            ('FB', '#1D9E75', 'c-Fos:JunB')):
        a.plot(mirs, [t[key] for t in up], '-', color=color, lw=2,
               label=f'{lab}, miR-200a rising')
        a.plot(mirs, [t[key] for t in down], '--', color=color, lw=2,
               label=f'{lab}, miR-200a falling')
    a.axvspan(*bounds, color='purple', alpha=0.12)
    a.set_xlabel('miR-200a (relative units)')
    a.set_ylabel('dimer concentration (nM)')
    a.set_title('A  Steady state across miR-200a', loc='left', fontsize=10)
    a.legend(fontsize=7)

    # time course: start on the c-Jun-high branch at the window midpoint
    q = dict(p)
    q['mir'] = 0.0
    y = steady_state(q)
    for m in np.linspace(0, m_mid, 7)[1:]:
        q['mir'] = m
        y = steady_state(q, y0=y)
    q.update(mir=m_mid, mir_pulse=True, mir_amp=pulse_amp, mir_peak_t=72.0,
             mir_sharp=2.0)
    sol = simulate(q, y0=y, t_end=600, n_pts=300)
    tt = [totals(sol.y[:, i]) for i in range(sol.y.shape[1])]
    b = ax[1]
    b.plot(sol.t, [t['FJ'] for t in tt], color='#D85A30', lw=2,
           label='c-Fos:c-Jun')
    b.plot(sol.t, [t['FB'] for t in tt], color='#1D9E75', lw=2,
           label='c-Fos:JunB')
    b.set_xlabel('time after injury (h)')
    b.set_ylabel('dimer concentration (nM)')
    b.set_title(f'B  After a transient miR-200a rise (baseline {m_mid:.2f})',
                loc='left', fontsize=10)
    b.legend(fontsize=8)
    fig.suptitle('Figure 3. Raising miR-200a moves c-Fos from c-Jun to JunB, '
                 'and the shift persists after the signal fades', fontsize=11)
    fig.tight_layout()
    print(f'window {bounds}; FJ {tt[0]["FJ"]:.2f} -> {tt[-1]["FJ"]:.2f}, '
          f'FB {tt[0]["FB"]:.2f} -> {tt[-1]["FB"]:.2f}')
    _save(fig, 'fig3_dimer_allocation.png')
    return fig


# ----------------------------------------------------------------------
# Figure 4 -- the c-Fos:c-Jun -> c-Jun edge
# ----------------------------------------------------------------------

def fig4_fj_edge(aJ_grid=(1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 24),
                 bfj_grid=(0, 4.5, 9, 13.5, 18, 27, 36)):
    """Across basal c-Jun production and the strength of the c-Fos:c-Jun ->
    c-Jun edge (all else at v3 values; homodimer loop 36): colour = the
    miR-200a level at which the high-c-Jun state disappears (upper edge of the
    two-state window). Grey = no bistability. Crosses = the window reaches
    miR-200a = 0, i.e. both states coexist even without miR-200a."""
    mirs = np.linspace(0, 2.5, 41)
    H = np.full((len(bfj_grid), len(aJ_grid)), np.nan)
    L = np.full_like(H, np.nan)
    for i, b in enumerate(bfj_grid):
        for j, a in enumerate(aJ_grid):
            p = v3_params()
            p.update(a_J=float(a), b_J_by_FJ=float(b))
            w = robust_window(p, mirs)
            if w['width'] > 0:
                H[i, j], L[i, j] = w['hi'], w['lo']
        print(f'FJ edge {b:>5g}: upper edges ' + ' '.join(
            '  -  ' if np.isnan(v) else f'{v:.2f}' for v in H[i]), flush=True)

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    cmap = plt.get_cmap('viridis').copy()
    cmap.set_bad('#D3D1C7')
    im = ax.pcolormesh(aJ_grid, bfj_grid, np.ma.masked_invalid(H),
                       shading='nearest', cmap=cmap)
    fig.colorbar(im, ax=ax, label='miR-200a needed to remove the '
                 'high-c-Jun state')
    zi, zj = np.where(L == 0)
    ax.plot(np.array(aJ_grid)[zj], np.array(bfj_grid)[zi], 'wx', ms=6,
            label='two states coexist even without miR-200a')
    ax.plot(16, 0, 'w^', ms=9, label='v3 (no FJ edge)')
    ax.plot(8, 9, 'wo', ms=9, label='reference with FJ edge')
    ax.set_xlabel('basal c-Jun production (nM/h)')
    ax.set_ylabel('strength of c-Fos:c-Jun -> c-Jun induction (nM/h)')
    ax.legend(fontsize=7, loc='upper right')
    ax.set_title('Figure 4. The c-Fos:c-Jun -> c-Jun edge adds to c-Jun drive: '
                 'a stronger edge\nlowers the basal c-Jun the switch tolerates '
                 '(grey: no bistability)', fontsize=10)
    fig.tight_layout()
    _save(fig, 'fig4_fj_edge.png')
    return fig, H, L


# ----------------------------------------------------------------------
# Figure 5 -- step 6
# ----------------------------------------------------------------------

def fig5_gfap_distinctness(p=None, observed=(1.47, 1.16, 1.86)):
    """Model GFAP in the two states, next to the single-cell estimate
    (JUN+ vs JUN- ependymoglia, 1 week post-injury, Lust et al. 2022).
    The observed value is a lower bound because of dropout."""
    p = p or ref_params()
    mirs = np.linspace(0, 2.0, 81)
    up, down, _, _ = branches(p, mirs)
    bounds = _window(mirs, up, down)
    w = robust_window(p)
    on, off = w['on'], w['off']
    r_g = on['GFAP'] / off['GFAP']
    r_j = on['cJun_total'] / off['cJun_total']

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2),
                           gridspec_kw={'width_ratios': [1.6, 1]})
    a = ax[0]
    a.plot(mirs, [t['GFAP'] for t in up], '-', lw=2, label='miR-200a rising')
    a.plot(mirs, [t['GFAP'] for t in down], '--', lw=2,
           label='miR-200a falling')
    a.axvspan(*bounds, color='purple', alpha=0.12)
    a.set_xlabel('miR-200a (relative units)')
    a.set_ylabel('GFAP (nM)')
    a.set_title('A  GFAP across miR-200a', loc='left', fontsize=10)
    a.legend(fontsize=8)

    b = ax[1]
    vals = [r_j, r_g, observed[0]]
    b.bar([0, 1], vals[:2], color=['#534AB7', '#888780'])
    b.errorbar(2, observed[0], yerr=[[observed[0] - observed[1]],
                                     [observed[2] - observed[0]]],
               fmt='o', color='#1D9E75', capsize=4)
    b.set_xticks([0, 1, 2])
    b.set_xticklabels(['model:\nc-Jun', 'model:\nGFAP',
                       'data:\nGFAP (lower bound)'], fontsize=8)
    b.axhline(1, color='k', ls=':', lw=1)
    b.set_ylabel('high-c-Jun state / low-c-Jun state')
    b.set_title('B  How different are the two states?', loc='left',
                fontsize=10)
    fig.suptitle('Figure 5. The two states differ clearly in c-Jun but only '
                 'modestly in GFAP', fontsize=11)
    fig.tight_layout()
    print(f'two-state ratios: c-Jun {r_j:.2f}x, GFAP {r_g:.2f}x')
    _save(fig, 'fig5_gfap_distinctness.png')
    return fig


def make_all():
    fig1_competition_vs_feedback()
    fig2_cfos_abundance()
    fig3_dimer_allocation()
    fig4_fj_edge()
    fig5_gfap_distinctness()


if __name__ == '__main__':
    make_all()
