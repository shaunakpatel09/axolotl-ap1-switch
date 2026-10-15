"""
v2 downstream experiments, re-run on the v3 architecture.

WHY THIS MODULE EXISTS
----------------------
Several findings (GFAP floor erosion, the developmental cliff, the global
LHS/ablation statistics, and the phi readout calibration) were produced on the
v2 architecture -- the one the Sabin calibration falsified on state
separation. They cannot be reported alongside a v3 model without being
recomputed. This module ports each of them to v3 and reports what survives.

Expect some findings NOT to survive. That is the point of the exercise.

Run:  python ap1_v3_downstream.py      (or import and call individually)

Author: Shaunak Patel
Model built with AI assistance (Claude); reviewed and owned by the author.
"""

import numpy as np
import matplotlib.pyplot as plt

from ap1_model_v3 import (N, SPECIES, FIXED, v3_params, rhs, hill, repress,
                          steady_state, totals, hysteresis_sweep,
                          find_attractors, cluster_states, lhs,
                          sample_initial_conditions)


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------

def window_states(p, mirs=None, t_end=1200.0):
    """Locate the bistable window on total cJun (the switching species) and
    return (mir_at_midpoint, ON_state_totals, OFF_state_totals, gap).
    Returns (None, None, None, 0.0) if there is no window."""
    mirs = mirs if mirs is not None else np.linspace(0, 3, 61)
    up, down = hysteresis_sweep(p, mirs, t_end=t_end)
    uJ = np.array([t['cJun_total'] for t in up])
    dJ = np.array([t['cJun_total'] for t in down])
    if uJ.max() < 1.5 * uJ.min():
        return None, None, None, 0.0
    mid = (uJ.max() + uJ.min()) / 2
    win = (dJ < mid) & (uJ > mid)
    if not win.any():
        return None, None, None, 0.0
    idx = np.where(win)[0]
    i = idx[len(idx) // 2]
    gap = mirs[np.argmax(uJ < mid)] - mirs[np.argmax(dJ < mid)]
    return mirs[i], up[i], down[i], gap


# ----------------------------------------------------------------------
# 1. GFAP floor erosion, re-run on v3
# ----------------------------------------------------------------------

def v3_gfap_floor(floors=None, plot=True):
    """v2 version measured contrast between mir=0 and mir=3 and found 38x
    collapsing to 14x with a small AP-1-independent drive. Under v3 the
    per-cell (within-window) contrast is only ~1.22x, so BOTH quantities are
    reported here and they mean different things:

      switch_contrast : ON vs OFF branch at matched miR INSIDE the window
                        -- the per-cell quantity, comparable to the data
      range_contrast  : GFAP at mir=0 vs mir=3 -- a whole-axis quantity that
                        mixes the switch with the graded decline outside it
    """
    floors = np.linspace(0, 4, 17) if floors is None else np.asarray(floors)
    switch_c, range_c = [], []
    for f in floors:
        p = v3_params()
        p['a_G_indep'] = float(f)
        _, on, off, gap = window_states(p)
        switch_c.append(on['GFAP'] / max(off['GFAP'], 1e-12)
                        if on is not None else np.nan)
        q = dict(p)
        q['mir'] = 0.0
        g_lo = totals(steady_state(q))['GFAP']
        q['mir'] = 3.0
        g_hi = totals(steady_state(q))['GFAP']
        range_c.append(g_lo / max(g_hi, 1e-12))
    switch_c, range_c = np.array(switch_c), np.array(range_c)

    print('floor  switch_contrast(per-cell)  range_contrast(mir 0 vs 3)')
    for f, s, r in zip(floors, switch_c, range_c):
        s_txt = 'no window' if np.isnan(s) else f'{s:.3f}x'
        print(f'{f:5.2f}   {s_txt:>12s}              {r:.2f}x')

    if plot:
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].plot(floors, switch_c, 'o-')
        ax[0].set_xlabel('AP-1-independent GFAP production (nM/h)')
        ax[0].set_ylabel('per-cell GFAP ON/OFF')
        ax[0].set_title('v3: within-window contrast')
        ax[1].plot(floors, range_c, 's-', color='crimson')
        ax[1].set_xlabel('AP-1-independent GFAP production (nM/h)')
        ax[1].set_ylabel('GFAP(mir=0) / GFAP(mir=3)')
        ax[1].set_title('v3: whole-axis contrast')
        plt.tight_layout()
        plt.show()
    return floors, switch_c, range_c


# ----------------------------------------------------------------------
# 2. Developmental axis, re-run on v3
# ----------------------------------------------------------------------

def v3_developmental_axis(fs=None, plot=True):
    """v2 found an abrupt 'birth' of a deep scar attractor at f ~ 0.65.

    Under v3 the reactive programme is parameterised differently (homodimer
    self-activation, high basal cJun). Scaling it down moves a_J toward the
    lower edge of the A14/P6 band, so the expected v3 behaviour is that
    bistability is LOST below some f rather than the scar state becoming
    shallow. Whether that reproduces a discontinuity is the open question.

    Mapping (exploratory, as in v2): f scales basal cJun production, the
    homodimer feedback strength, and the PNN/CSPG-linked GFAP floor together;
    f = 1 is adult-like.
    """
    fs = np.linspace(0.5, 1.4, 10) if fs is None else np.asarray(fs)
    base = v3_params()
    A_J, B_JJ, FLOOR = base['a_J'], base['b_J_by_JJ'], 3.0
    gaps, seps, scar = [], [], []
    for f in fs:
        p = v3_params()
        p['a_J'] = A_J * f
        p['b_J_by_JJ'] = B_JJ * f
        p['a_G_indep'] = FLOOR * f
        _, on, off, gap = window_states(p)
        gaps.append(gap)
        seps.append(on['cJun_total'] / max(off['cJun_total'], 1e-12)
                    if on is not None else np.nan)
        q = dict(p)
        q['mir'] = 0.0
        scar.append(totals(steady_state(q))['GFAP'])
    gaps, seps, scar = np.array(gaps), np.array(seps), np.array(scar)

    print('  f    window_gap   cJun_sep    GFAP at mir=0')
    for f, g, s, sc in zip(fs, gaps, seps, scar):
        s_txt = '   --   ' if np.isnan(s) else f'{s:6.2f}x'
        print(f'{f:5.2f}   {g:8.2f}   {s_txt}    {sc:6.1f}')
    alive = fs[gaps > 0.03]
    if len(alive):
        print(f'bistability present for f in [{alive.min():.2f}, '
              f'{alive.max():.2f}]')
    else:
        print('no bistability anywhere in this range')

    if plot:
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].plot(fs, gaps, 'o-')
        ax[0].set_xlabel('developmental scale f (1 = adult-like)')
        ax[0].set_ylabel('bistable window width')
        ax[0].set_title('v3: does the switch exist?')
        ax[1].plot(fs, scar, 's-', color='crimson')
        ax[1].set_xlabel('developmental scale f (1 = adult-like)')
        ax[1].set_ylabel('GFAP at mir = 0 (nM)')
        ax[1].set_title('v3: scar-state level')
        plt.tight_layout()
        plt.show()
    return fs, gaps, seps, scar


# ----------------------------------------------------------------------
# 3. Global LHS scan + ablation, re-run on v3
# ----------------------------------------------------------------------

V3_SCAN_RANGES = {
    'a_F': (2.0, 20.0),
    'a_J': (4.0, 40.0),
    'a_B': (1.0, 10.0),
    'a_G': (0.1, 2.0),
    'b_J_by_JJ': (8.0, 120.0),
    'Km_auto': (1.0, 20.0),
    'g_F': (0.3, 2.0),
    'g_J': (0.2, 1.5),
    'g_B': (0.3, 2.0),
    'g_G': (0.05, 1.0),
    'junb_weakness': (1.0, 6.0),
    'koff_JJ': (300.0, 2000.0),
}
V3_SCAN_KEYS = list(V3_SCAN_RANGES)


def v3_sample_parameter_sets(n, rng):
    keys = V3_SCAN_KEYS
    u = lhs(n, len(keys), rng)
    out = []
    for i in range(n):
        p = v3_params()
        for d, k in enumerate(keys):
            lo, hi = V3_SCAN_RANGES[k]
            p[k] = 10 ** (np.log10(lo) + u[i, d] *
                          (np.log10(hi) - np.log10(lo)))
        p['koff_BB'] = p['koff_JB'] = p['koff_JJ']
        p['mir'] = 0.0
        out.append(p)
    return out


def v3_global_scan(n_param_sets=250, n_initial=8, overrides=None, seed=42):
    rng = np.random.default_rng(seed)
    psets = v3_sample_parameter_sets(n_param_sets, rng)
    y0s = sample_initial_conditions(n_initial, rng)
    counts = []
    for p in psets:
        if overrides:
            p.update(overrides)
        n_states, _ = find_attractors(p, y0s)
        counts.append(n_states)
    return psets, np.array(counts)


def v3_ablation(n_param_sets=250, n_initial=8, plot=True):
    """v2 gave 18% / 24% / 2% / 0% for full / no-FJ / no-JJ / no-autoreg.
    v3 has the FJ->cJun edge already removed by construction (Z2), so the
    meaningful comparison is full vs homodimer-loop-ablated."""
    conditions = {'v3 full model': {}, 'JJ->cJun loop removed': {'b_J_by_JJ': 0.0}}
    out = {}
    for name, ov in conditions.items():
        _, counts = v3_global_scan(n_param_sets, n_initial, overrides=ov)
        valid = counts[counts > 0]
        pct = 100.0 * float(np.mean(valid > 1)) if len(valid) else 0.0
        out[name] = pct
        print(f'{name:24s} multistable: {pct:5.1f}%  (n={len(valid)})')
    if plot:
        plt.figure(figsize=(6, 4))
        plt.bar(list(out), list(out.values()))
        plt.ylabel('% of parameter sets multistable')
        plt.title('v3: bistability requires the homodimer loop')
        plt.xticks(rotation=15, ha='right')
        plt.tight_layout()
        plt.show()
    return out


# ----------------------------------------------------------------------
# 4. phi recalibration on v3 (expected to FAIL -- that is prediction P4)
# ----------------------------------------------------------------------

MEASURED_PER_CELL = {'GFAP': 2.14, 'vimentin': 4.13, 'CSPG4': 4.71}


def v3_phi_ceiling(verbose=True):
    """Under v3, compute the maximum per-cell fold-change any promoter reading
    cFos:cJun activity can show, and compare with the measured per-cell values.
    phi in [0,1] can only ATTENUATE this ceiling, never exceed it."""
    p = v3_params()
    _, on, off, gap = window_states(p)
    if on is None:
        print('no window; cannot compute ceiling')
        return None
    A_on = (hill(on['FJ'], p['Km'], p['n'])
            * repress(on['FB'], p['Km_rep'], p['n_rep']))
    A_off = (hill(off['FJ'], p['Km'], p['n'])
             * repress(off['FB'], p['Km_rep'], p['n_rep']))
    ceiling = A_on / A_off
    if verbose:
        print(f'v3 AP-1 (cFos:cJun) activity ON/OFF = {ceiling:.3f}x')
        print('measured per-cell targets:')
        for g, v in MEASURED_PER_CELL.items():
            verdict = 'reachable' if v <= ceiling else 'UNREACHABLE'
            print(f'  {g:9s} {v:.2f}x  -> {verdict}')
        print('\n=> phi can only attenuate, so any target above the ceiling '
              'cannot be\n   fitted by an FJ-only readout. This is prediction '
              'P4: the scar\n   promoters must integrate additional '
              '(Jun-homodimer) input.')
    return ceiling


if __name__ == '__main__':
    print('=' * 60)
    print('1. GFAP floor erosion on v3')
    print('=' * 60)
    v3_gfap_floor(plot=False)
    print('\n' + '=' * 60)
    print('2. Developmental axis on v3')
    print('=' * 60)
    v3_developmental_axis(plot=False)
    print('\n' + '=' * 60)
    print('3. Global scan + ablation on v3')
    print('=' * 60)
    v3_ablation(plot=False)
    print('\n' + '=' * 60)
    print('4. phi ceiling on v3')
    print('=' * 60)
    v3_phi_ceiling()
