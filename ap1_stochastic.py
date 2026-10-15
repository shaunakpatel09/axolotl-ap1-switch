"""
Stochastic ensemble for the v3 AP-1 model: cell-to-cell variability via slow
extrinsic noise on cJun production (Ornstein-Uhlenbeck, tau = 10 h), so a
population of bistable cells partitions into cJun+ / cJun- fractions -- the
quantity Sabin et al. Fig 5c actually measures (8% control -> 30% inhibitor).

Design (fit-one-predict-one):
  - sigma (noise amplitude) is fitted to the CONTROL fraction (8%) ONLY
  - the INHIBITOR fraction is then an out-of-sample prediction
  - miR-axis placement of the two conditions is fixed A PRIORI from window
    geometry (control just above the window, inhibitor inside it) and must
    not be tuned to the outcome

RESULT under the round-2 reference model, over four seeds (use fit_and_predict
for a single run, seed_check for the honest version):

  control   (miR-200a 0.95, 120 cells)  0.046 +- 0.035   target   0.08  [fitted]
  inhibitor (miR-200a 0.58, 150 cells)  0.482 +- 0.039   measured 0.30  [predicted]

Direction and scale correct; magnitude ~1.6x over. Essentially unchanged by
the v3 -> round-2 architecture revision, which gave 0.52.

TWO HONEST CAVEATS, both found by running more than one seed:
  1. The between-seed spread on the control fraction (sd 0.035) is comparable
     to the target itself (0.08), so a single 120-cell run cannot pin sigma
     down. Control fraction by sigma, seed 11: 0.000 at 0.28, 0.008 at 0.30,
     0.017 at 0.32, 0.083 at 0.35. sigma is constrained only to roughly
     0.33-0.40.
  2. The single-seed numbers first reported here flattered the result in both
     directions: the fit used the highest control fraction of four seeds and
     the prediction the lowest inhibitor fraction of four.

Dominant remaining uncertainty: the a-priori miR mapping. Do NOT retune the
mapping to improve this number -- report it as-is or redesign the mapping on
independent grounds first.

Author: Shaunak Patel
Model built with AI assistance (Claude); reviewed and owned by the author.
"""
import numpy as np
from scipy.integrate import solve_ivp
from ap1_model_v3 import rhs, steady_state, totals, N
from round2 import ref_params, robust_window

TAU, SEG, T_TOTAL = 10.0, 5.0, 250.0


def geometry(p=None):
    """Threshold and the two miR-200a placements, derived from the window
    rather than hard-coded.

      threshold     geometric midpoint of the ON and OFF c-Jun levels
      control       just above the window's upper edge (+0.15), so most cells
                    sit in the low-c-Jun state
      inhibitor     the window midpoint, where both states are available

    Both placements follow the SAME a-priori rule used before the round-2
    revision, and must not be retuned to improve the prediction."""
    p = p or ref_params()
    w = robust_window(p)
    thr = float(np.sqrt(w['on']['cJun_total'] * w['off']['cJun_total']))
    return thr, w['hi'] + 0.15, 0.5 * (w['lo'] + w['hi']), w


def ensemble_fraction(mir, sigma, n_cells=120, seed=0, p_base=None,
                      threshold=None):
    """Fraction of cells above the cJun+ threshold at steady state, under OU
    extrinsic noise applied multiplicatively to cJun production."""
    r = np.random.default_rng(seed)
    p0 = p_base or ref_params()
    THRESHOLD_CJUN = threshold if threshold is not None else geometry(p0)[0]
    nseg = int(T_TOTAL / SEG)
    hits = 0
    for _ in range(n_cells):
        p = dict(p0)
        p['mir'] = mir
        y = steady_state(p, y0=np.zeros(N), t_end=800)
        eta = 0.0
        a0, b0 = p['a_J'], p['b_J_by_JJ']
        decay = np.exp(-SEG / TAU)
        kick = sigma * np.sqrt(1 - decay ** 2)
        for _ in range(nseg):
            eta = eta * decay + kick * r.standard_normal()
            p['a_J'], p['b_J_by_JJ'] = a0 * np.exp(eta), b0 * np.exp(eta)
            sol = solve_ivp(rhs, (0, SEG), y, args=(p,), method='LSODA',
                            rtol=1e-6, atol=1e-8)
            y = sol.y[:, -1]
        hits += totals(y)['cJun_total'] > THRESHOLD_CJUN
    return hits / n_cells


def fit_and_predict(sigmas=(0.15, 0.2, 0.25, 0.35), n_cells=120, p=None):
    """Fit sigma on the control fraction, then predict the inhibitor fraction.

    Only the control fraction (0.08) is fitted. The inhibitor fraction (0.30)
    is never shown to the fit -- it is the out-of-sample test."""
    p = p or ref_params()
    thr, m_ctrl, m_inh, w = geometry(p)
    print(f'window {w["lo"]:.2f}-{w["hi"]:.2f}; ON/OFF c-Jun '
          f'{w["on"]["cJun_total"]:.1f}/{w["off"]["cJun_total"]:.1f}; '
          f'threshold {thr:.1f}')
    print(f'a-priori placements: control miR-200a {m_ctrl:.2f}, '
          f'inhibitor {m_inh:.2f}\n')
    fit = None
    for s in sigmas:
        f = ensemble_fraction(m_ctrl, s, n_cells=n_cells, seed=11, p_base=p,
                              threshold=thr)
        print(f'sigma={s:.2f} -> control fraction {f:.2f} (target 0.08)',
              flush=True)
        if fit is None and abs(f - 0.08) <= 0.05:
            fit = s
    fit = fit if fit is not None else sigmas[0]
    f_kd = ensemble_fraction(m_inh, fit, n_cells=n_cells + 30, seed=29,
                             p_base=p, threshold=thr)
    print(f'\nfitted sigma={fit:.2f}; PREDICTED inhibitor fraction = '
          f'{f_kd:.2f} (measured 0.30)')
    return fit, f_kd


def seed_check(sigma=0.35, seeds=(11, 29, 77, 101), p=None):
    """The reportable version: the same two fractions over several seeds, so
    the between-seed spread is visible rather than hidden in one run."""
    p = p or ref_params()
    thr, m_ctrl, m_inh, w = geometry(p)
    out = {}
    for lab, m, n, target in (('control (fitted)', m_ctrl, 120, 0.08),
                              ('inhibitor (predicted)', m_inh, 150, 0.30)):
        v = np.array([ensemble_fraction(m, sigma, n_cells=n, seed=s,
                                        p_base=p, threshold=thr)
                      for s in seeds])
        out[lab] = v
        print(f'{lab:24s} {v.mean():.3f} +- {v.std(ddof=1):.3f}   '
              f'(target {target:.2f})   seeds: '
              + ' '.join(f'{x:.3f}' for x in v), flush=True)
    return out


if __name__ == '__main__':
    seed_check()
