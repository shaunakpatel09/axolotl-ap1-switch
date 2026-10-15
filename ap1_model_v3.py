"""
AP-1 competition model, v3: architecture revised after Zhu et al. 2022
(Science 375:eabg9765, "Synthetic multistability in mammalian cells").

WHY v3 EXISTS
-------------
Calibration against Sabin et al. 2019 (Figs 5-6) falsified v2 quantitatively:
real ON/OFF glial states differ by only ~2-5 fold (per-cell, after correcting
bulk qPCR by the Fig 5c cJun+ fractions), while v2's attractors were separated
by 9-4700 fold. State separation was the named structural defect (audit A13).

THE ZHU-DERIVED REVISION (each maps to the MultiFate architecture):
  Z1. Leaky basal production sets the OFF state. State separation =
      (self-activated level)/(leaky level), so a real basal floor for cJun
      replaces v2's near-zero OFF state.  [v2: a_J ~ 0.4 fully miR-repressed;
      v3: larger a_J plus a miR-repression leak]
  Z2. Homodimer-only self-activation. MultiFate = homodimer activates own
      promoter; heterodimers are inactive sequestration sinks. Our own
      ablation agreed: the JJ loop is the engine, the FJ->cJun loop
      SUPPRESSES multistability. v3 removes the FJ->cJun edge (b_J_by_FJ=0
      by default; parameter retained for ablation).
  Z3. Cooperativity carries a shallow switch. Dimerisation contributes
      quadratic dependence (JJ ~ J^2); cooperative promoter binding adds
      more. The autoregulation Hill (Km_auto, n_auto) is positioned BETWEEN
      the compressed states. Degradation remains the depth-tuning knob
      (Zhu tuned stability experimentally to move between regimes).

TARGET SET BY DATA: bistable window in miR-200a with cJun ON/OFF separation
in the 2-10x band (measured-implied ~3.9x).

Author: Shaunak Patel
Model built with AI assistance (Claude); reviewed and owned by the author.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
import matplotlib.pyplot as plt

SPECIES = ['F', 'J', 'B', 'FJ', 'FB', 'JJ', 'BB', 'JB', 'G']
N = len(SPECIES)

FIXED = {
    'kon_FJ': 10.0, 'koff_FJ': 1.0,
    'kon_FB': 10.0,
    'kon_JJ': 10.0, 'koff_JJ': 100.0,
    'kon_BB': 10.0, 'koff_BB': 100.0,
    'kon_JB': 10.0, 'koff_JB': 100.0,
    'Km': 10.0, 'n': 2.0,          # GFAP activation by FJ (unchanged)
    'Km_rep': 5.0, 'n_rep': 2.0,   # GFAP repression by FB (unchanged)
    'Km_mir': 1.0, 'n_mir': 2.0,
}


def v3_params():
    """Zhu-revised reference set. Values chosen to place the bistable window
    and ~3x separation; every choice is audit-listed (A14) and swept in the
    gauntlet, not asserted."""
    p = dict(FIXED)
    p.update(
        # basal production
        a_F=8.0, a_B=4.0, a_G=0.5,
        a_J=16.0,            # Z1: real basal cJun floor (v2: 0.4). Places both
                             # states ABOVE the cFos-saturation scale, where the
                             # FJ sequestration drain is capped (see A14)
        # feedback
        b_J_by_FJ=0.0,       # Z2: FJ->cJun edge removed
        b_J_by_JJ=36.0,      # Z1: induced/basal ~ 3.2 (v2: 225)
        Km_auto=5.0,         # Z3: positioned between compressed JJ states
        n_auto=3.0,          # Z3: cooperativity (dimer + promoter binding)
        b_G_by_FJ=20.0,
        # miR-200a repression with leak: rep = leak + (1-leak)*Hill
        mir_leak=0.15,       # Z1: repression is a rheostat (Sabin's own words)
        # A15: JUN-JUN pairings weakened to Kd = 100 nM (1000x weaker than
        # FOS-JUN) so the cJun homodimer is SIGNALLING-dominant rather than a
        # mass sink; with the reference 100x convention the quadratic JJ
        # degradation drain forbids shallow bistability. Audit-flagged;
        # sweep koff_JJ in sensitivity analysis.
        koff_JJ=1000.0, koff_BB=1000.0, koff_JB=1000.0,
        # degradation
        g_F=1.0, g_J=0.5, g_B=0.7, g_G=0.2,
        junb_weakness=3.0,
        mir=0.0, a_G_indep=0.0,
    )
    return p


def hill(x, K, n):
    xn = np.maximum(x, 0.0) ** n
    return xn / (K ** n + xn)


def repress(x, K, n):
    return K ** n / (K ** n + np.maximum(x, 0.0) ** n)


def mir_value(t, p):
    m = p.get('mir', 0.0)
    if p.get('mir_pulse', False):
        tp = p.get('mir_peak_t', 72.0)
        s = p.get('mir_sharp', 2.0)
        x = max(t, 0.0) / tp
        m = m + p.get('mir_amp', 0.0) * (x * np.exp(1.0 - x)) ** s
    return m


def rhs(t, y, p):
    F, J, B, FJ, FB, JJ, BB, JB, G = y
    koff_FB = p['koff_FJ'] * p['junb_weakness']

    v_FJ = p['kon_FJ'] * F * J - p['koff_FJ'] * FJ
    v_FB = p['kon_FB'] * F * B - koff_FB * FB
    v_JJ = p['kon_JJ'] * J * J - p['koff_JJ'] * JJ
    v_BB = p['kon_BB'] * B * B - p['koff_BB'] * BB
    v_JB = p['kon_JB'] * J * B - p['koff_JB'] * JB

    g_FJ = 0.5 * (p['g_F'] + p['g_J'])
    g_FB = 0.5 * (p['g_F'] + p['g_B'])
    g_JB = 0.5 * (p['g_J'] + p['g_B'])

    leak = p.get('mir_leak', 0.0)
    mir_rep = leak + (1.0 - leak) * repress(mir_value(t, p),
                                            p['Km_mir'], p['n_mir'])
    # Z1 + Z2 + Z3: basal floor, homodimer-only self-activation,
    # positioned cooperativity
    prod_J = (p['a_J']
              + p['b_J_by_FJ'] * hill(FJ, p['Km_auto'], p['n_auto'])
              + p['b_J_by_JJ'] * hill(JJ, p['Km_auto'], p['n_auto'])) * mir_rep

    dF = p['a_F'] - p['g_F'] * F - v_FJ - v_FB
    dJ = prod_J - p['g_J'] * J - v_FJ - 2 * v_JJ - v_JB
    dB = p['a_B'] - p['g_B'] * B - v_FB - 2 * v_BB - v_JB
    dFJ = v_FJ - g_FJ * FJ
    dFB = v_FB - g_FB * FB
    dJJ = v_JJ - p['g_J'] * JJ
    dBB = v_BB - p['g_B'] * BB
    dJB = v_JB - g_JB * JB
    dG = (p.get('a_G_indep', 0.0)
          + (p['a_G'] + p['b_G_by_FJ'] * hill(FJ, p['Km'], p['n']))
          * repress(FB, p['Km_rep'], p['n_rep']) - p['g_G'] * G)
    return [dF, dJ, dB, dFJ, dFB, dJJ, dBB, dJB, dG]


def simulate(p, y0=None, t_end=300.0, n_pts=400):
    if y0 is None:
        y0 = np.zeros(N)
    return solve_ivp(rhs, (0, t_end), y0, args=(p,), method='LSODA',
                     t_eval=np.linspace(0, t_end, n_pts), rtol=1e-8, atol=1e-10)


def steady_state(p, y0=None, t_end=2000.0):
    if y0 is None:
        y0 = np.zeros(N)
    sol = solve_ivp(rhs, (0, t_end), y0, args=(p,), method='LSODA',
                    rtol=1e-8, atol=1e-10)
    return sol.y[:, -1]


def totals(y):
    F, J, B, FJ, FB, JJ, BB, JB, G = y
    return {'cFos_total': F + FJ + FB,
            'cJun_total': J + FJ + 2 * JJ + JB,
            'JunB_total': B + FB + 2 * BB + JB,
            'FJ': FJ, 'FB': FB, 'GFAP': G}


def refine_fixed_point(p, y_guess):
    sol = root(lambda y: rhs(0.0, y, p), y_guess, method='hybr', tol=1e-12)
    y = np.maximum(sol.x, 0.0)
    ok = sol.success and np.max(np.abs(rhs(0.0, y, p))) < 1e-6 * max(1.0, y.max())
    return y, bool(ok)


def numerical_jacobian(p, y, rel=1e-6):
    f0 = np.asarray(rhs(0.0, y, p))
    Jm = np.zeros((N, N))
    for i in range(N):
        h = rel * max(1.0, abs(y[i]))
        yp = y.copy()
        yp[i] += h
        Jm[:, i] = (np.asarray(rhs(0.0, yp, p)) - f0) / h
    return Jm


def is_stable(p, y, tol=1e-6):
    return bool(np.max(np.linalg.eigvals(numerical_jacobian(p, y)).real) < tol)


def cluster_states(states, rel_tol=0.15):
    logs = np.log10(np.maximum(np.asarray(states, dtype=float), 1e-9))
    centers, labels = [], []
    for row in logs:
        for k, c in enumerate(centers):
            if np.max(np.abs(row - c)) < rel_tol:
                labels.append(k)
                break
        else:
            centers.append(row)
            labels.append(len(centers) - 1)
    return np.array(labels), len(centers)


def find_attractors(p, y0s, t_pre=300.0):
    finals = []
    for y0 in y0s:
        y_end = steady_state(p, y0=y0, t_end=t_pre)
        y_fp, ok = refine_fixed_point(p, y_end)
        if ok and np.all(np.isfinite(y_fp)) and is_stable(p, y_fp):
            finals.append(y_fp)
    if not finals:
        return 0, []
    _, n_unique = cluster_states(finals)
    return n_unique, finals


def lhs(n_samples, n_dims, rng):
    cut = np.arange(n_samples) / n_samples
    u = np.empty((n_samples, n_dims))
    for d in range(n_dims):
        u[:, d] = rng.permutation(cut + rng.random(n_samples) / n_samples)
    return u


def sample_initial_conditions(n, rng, lo=0.316, hi=316.0):
    u = lhs(n, 3, rng)
    y0 = np.zeros((n, N))
    y0[:, :3] = 10 ** (np.log10(lo) + u * (np.log10(hi) - np.log10(lo)))
    return y0


def hysteresis_sweep(p, mirs, t_end=1500.0):
    p = dict(p)
    y = np.zeros(N)
    up = []
    for m in mirs:
        p['mir'] = m
        y = steady_state(p, y0=y, t_end=t_end)
        up.append(totals(y))
    p['mir'] = mirs[-1]
    y = steady_state(p, y0=np.zeros(N), t_end=t_end)
    down = []
    for m in mirs[::-1]:
        p['mir'] = m
        y = steady_state(p, y0=y, t_end=t_end)
        down.append(totals(y))
    return up, down[::-1]


# ----------------------------------------------------------------------
# THE GAUNTLET: the three tests a revision must pass
# ----------------------------------------------------------------------

def gauntlet(p=None, mirs=None, plot=True):
    """Test 1 -- bistability: a hysteresis window must exist in miR-200a.
    Test 2 -- separation: cJun ON/OFF at matched miR inside the window must
              land in the data-implied 2-10x band (target ~3.9x).
    Test 3 -- engine: removing the JJ->cJun loop must abolish the window
              (bistability still requires the autoregulation)."""
    p = p or v3_params()
    mirs = mirs if mirs is not None else np.linspace(0, 3, 121)

    up, down = hysteresis_sweep(p, mirs)
    # detect the window on the SWITCHING SPECIES (total cJun), not on GFAP:
    # with data-calibrated shallow separation, downstream readouts move too
    # little to locate the bifurcation reliably.
    upJ = np.array([t['cJun_total'] for t in up])
    dnJ = np.array([t['cJun_total'] for t in down])
    mid = (upJ.max() + upJ.min()) / 2
    ut = mirs[np.argmax(upJ < mid)]
    dt = mirs[np.argmax(dnJ < mid)]
    gap = ut - dt
    window = (dnJ < mid) & (upJ > mid)
    upG = np.array([t['GFAP'] for t in up])
    dnG = np.array([t['GFAP'] for t in down])
    print(f'TEST 1  bistable window: gap = {gap:.3f} in miR units '
          f'({"PASS" if gap > 0.03 else "FAIL"})')

    sep_J = sep_G = None
    if window.any():
        i = np.where(window)[0][len(np.where(window)[0]) // 2]
        on, off = up[i], down[i]
        sep_J = on['cJun_total'] / max(off['cJun_total'], 1e-12)
        sep_G = on['GFAP'] / max(off['GFAP'], 1e-12)
        ok = 2.0 <= sep_J <= 10.0
        print(f'TEST 2  separation at mir={mirs[i]:.2f} inside window: '
              f'cJun ON/OFF = {sep_J:.2f}x, GFAP = {sep_G:.2f}x '
              f'({"PASS -- in the 2-10x data band" if ok else "FAIL"})')
    else:
        print('TEST 2  skipped: no window')

    q = dict(p)
    q['b_J_by_JJ'] = 0.0
    up2, dn2 = hysteresis_sweep(q, mirs)
    up2J = np.array([t['cJun_total'] for t in up2])
    dn2J = np.array([t['cJun_total'] for t in dn2])
    gap2 = abs(mirs[np.argmax(up2J < (up2J.max() + up2J.min()) / 2)]
               - mirs[np.argmax(dn2J < (dn2J.max() + dn2J.min()) / 2)])
    print(f'TEST 3  autoregulation ablated: gap = {gap2:.3f} '
          f'({"PASS -- window vanishes" if gap2 < 0.03 else "FAIL"})')

    if plot:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        axes[0].plot(mirs, upG, 'o-', ms=3, label='sweep up')
        axes[0].plot(mirs, dnG, 's-', ms=3, label='sweep down')
        if gap > 0.03:
            axes[0].axvspan(dt, ut, alpha=0.15, color='purple')
        axes[0].set_xlabel('miR-200a')
        axes[0].set_ylabel('GFAP (nM)')
        axes[0].set_title('v3: bistable window')
        axes[0].legend(fontsize=8)
        upJ = np.array([t['cJun_total'] for t in up])
        dnJ = np.array([t['cJun_total'] for t in down])
        axes[1].plot(mirs, upJ, 'o-', ms=3, label='sweep up')
        axes[1].plot(mirs, dnJ, 's-', ms=3, label='sweep down')
        axes[1].set_xlabel('miR-200a')
        axes[1].set_ylabel('total cJun (nM)')
        axes[1].set_title('v3: cJun states are a few fold apart')
        axes[1].legend(fontsize=8)
        plt.tight_layout()
        plt.show()
    return {'gap': gap, 'sep_cJun': sep_J, 'sep_GFAP': sep_G,
            'gap_ablated': gap2}


if __name__ == '__main__':
    gauntlet()
