"""
test_model_properties.py -- property tests for the four invariants the paper
argues from.

Each test corresponds to a claim made in the manuscript, not to an
implementation detail:

  T1  dimer stoichiometry      homodimerisation consumes and releases TWO
                               monomers.  Section 5.1 asserts the factor-of-two
                               bookkeeping; every concentration in the paper
                               depends on it.
  T2  equilibrium stability    states reported as steady and stable actually
                               satisfy the residual and Jacobian criteria of
                               Section 5.3, and the branches are separated by a
                               genuine unstable equilibrium (saddle-node), not
                               by a sweep discontinuity.
  T3  readout independence     changing GFAP-only parameters cannot move the
                               upstream AP-1 species.  Section 2.4 reports
                               1.8e-15 for this; the test re-derives it.
  T4  intervention equivalence JNK-activity and T-5224-with-homodimer modify
                               the same three c-Jun parameters and differ only
                               in a term entering the GFAP equation, so they
                               are REQUIRED to give identical c-Jun.  Section
                               2.8 reports 0.00e+00; the test re-derives it.

Run:  python3 test_model_properties.py
Exit code 0 = all properties hold.

Author: Shaunak Patel.  Built with AI assistance (Claude); reviewed and owned
by the author.
"""

import sys

import numpy as np

from ap1_model_v3 import (SPECIES, N, v3_params, rhs, steady_state,
                          refine_fixed_point, numerical_jacobian, totals)
from round2 import ref_params, robust_window

SEED = 20261015
UPSTREAM = [i for i, s in enumerate(SPECIES) if s != 'G']

# GFAP-only parameters: these appear in dG and nowhere else.
READOUT_ONLY = ['a_G', 'b_G_by_FJ', 'Km', 'n', 'Km_rep', 'n_rep', 'g_G',
                'a_G_indep']


def _report(name, ok, detail):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    for line in detail.splitlines():
        print(f"         {line}")
    return ok


# ----------------------------------------------------------------- T1

def test_dimer_stoichiometry(n_draws=500):
    """With production and degradation switched off, dimerisation alone must
    conserve each monomer's total pool:

        cFos_total = F + FJ + FB
        cJun_total = J + FJ + 2*JJ + JB
        JunB_total = B + FB + 2*BB + JB

    The factors of two are the property under test.  If dJ carried -v_JJ
    instead of -2*v_JJ, cJun_total would drift and this test would fail.
    """
    p = dict(ref_params())
    for k in ('a_F', 'a_J', 'a_B', 'a_G', 'a_G_indep',
              'b_J_by_FJ', 'b_J_by_JJ', 'b_G_by_FJ',
              'g_F', 'g_J', 'g_B', 'g_G'):
        p[k] = 0.0

    rng = np.random.default_rng(SEED)
    worst = {'cFos_total': 0.0, 'cJun_total': 0.0, 'JunB_total': 0.0}
    for _ in range(n_draws):
        y = rng.uniform(0.0, 50.0, N)
        dy = np.asarray(rhs(0.0, y, p))
        d = {
            'cFos_total': dy[SPECIES.index('F')] + dy[SPECIES.index('FJ')]
                          + dy[SPECIES.index('FB')],
            'cJun_total': dy[SPECIES.index('J')] + dy[SPECIES.index('FJ')]
                          + 2 * dy[SPECIES.index('JJ')]
                          + dy[SPECIES.index('JB')],
            'JunB_total': dy[SPECIES.index('B')] + dy[SPECIES.index('FB')]
                          + 2 * dy[SPECIES.index('BB')]
                          + dy[SPECIES.index('JB')],
        }
        scale = max(1.0, float(np.max(np.abs(dy))))
        for k, v in d.items():
            worst[k] = max(worst[k], abs(v) / scale)

    tol = 1e-12
    ok = all(v < tol for v in worst.values())
    detail = (f"{n_draws} random states, seed {SEED}\n"
              + "\n".join(f"max |d({k})/dt| / scale = {v:.3e}"
                          for k, v in worst.items())
              + f"\ntolerance {tol:.0e}")

    # negative control: a deliberately wrong factor must be caught
    y = rng.uniform(0.0, 50.0, N)
    dy = np.asarray(rhs(0.0, y, p))
    wrong = (dy[SPECIES.index('J')] + dy[SPECIES.index('FJ')]
             + 1 * dy[SPECIES.index('JJ')] + dy[SPECIES.index('JB')])
    detail += (f"\nnegative control (factor 1 instead of 2): "
               f"|drift| = {abs(wrong):.3e}  (must be non-zero)")
    ok = ok and abs(wrong) > 1e-9
    return _report("T1 dimer stoichiometry", ok, detail)


# ----------------------------------------------------------------- T2

def test_equilibrium_stability():
    """Both reported states satisfy the Section 5.3 acceptance criteria, and
    an unstable equilibrium separates them."""
    p = dict(ref_params())
    w = robust_window(p)
    p['mir'] = w['m_mid']

    y_hi = np.zeros(N)
    y_hi[SPECIES.index('J')] = 100.0
    hi = steady_state(p, y0=y_hi, t_end=5000.0)
    lo = steady_state(p, y0=np.zeros(N), t_end=5000.0)
    if totals(hi)['cJun_total'] < totals(lo)['cJun_total']:
        hi, lo = lo, hi

    lines, ok = [], True
    fps = {}
    for name, y0 in (('high', hi), ('low', lo)):
        y, refined = refine_fixed_point(p, y0)
        res = float(np.max(np.abs(rhs(0.0, y, p))))
        scale = max(1.0, float(np.max(np.abs(y))))
        ev = float(np.max(np.linalg.eigvals(numerical_jacobian(p, y)).real))
        good = refined and (res < 1e-6 * scale) and ev < 0.0
        ok = ok and good
        fps[name] = y
        lines.append(f"{name:5s} cJun_total={totals(y)['cJun_total']:8.3f}  "
                     f"residual={res:.3e} (<{1e-6*scale:.1e})  "
                     f"max Re(lambda)={ev:+.4f}")

    # a saddle must sit between the two stable branches
    saddle, found = None, False
    for frac in np.linspace(0.05, 0.95, 37):
        guess = fps['low'] + frac * (fps['high'] - fps['low'])
        y, refined = refine_fixed_point(p, guess)
        if not refined or np.any(y < -1e-9):
            continue
        cj = totals(y)['cJun_total']
        if not (totals(fps['low'])['cJun_total'] * 1.02 < cj
                < totals(fps['high'])['cJun_total'] * 0.98):
            continue
        ev = float(np.max(np.linalg.eigvals(numerical_jacobian(p, y)).real))
        if ev > 0:
            saddle, found = (cj, ev), True
            break
    if found:
        lines.append(f"saddle cJun_total={saddle[0]:8.3f}  "
                     f"max Re(lambda)={saddle[1]:+.4f} (unstable, as required)")
    else:
        lines.append("saddle NOT FOUND between the branches")
    ok = ok and found
    return _report("T2 equilibrium stability", ok, "\n".join(lines))


# ----------------------------------------------------------------- T3

def test_readout_independence(n_draws=40):
    """Perturbing GFAP-only parameters must leave all eight upstream species
    untouched."""
    p = dict(ref_params())
    w = robust_window(p)
    p['mir'] = w['m_mid']
    y_hi = np.zeros(N)
    y_hi[SPECIES.index('J')] = 100.0
    base = steady_state(p, y0=y_hi, t_end=5000.0)

    rng = np.random.default_rng(SEED + 1)
    worst = 0.0
    for _ in range(n_draws):
        q = dict(p)
        for k in READOUT_ONLY:
            if k in ('n', 'n_rep'):
                q[k] = float(rng.uniform(1.0, 4.0))
            else:
                q[k] = float(p.get(k, 1.0)) * float(rng.uniform(0.1, 10.0)) \
                       + float(rng.uniform(0.0, 5.0))
        y = steady_state(q, y0=y_hi, t_end=5000.0)
        rel = np.abs(y[UPSTREAM] - base[UPSTREAM]) / np.maximum(
            np.abs(base[UPSTREAM]), 1e-12)
        worst = max(worst, float(np.max(rel)))

    tol = 1e-8
    ok = worst < tol
    detail = (f"{n_draws} random readout parameterisations, seed {SEED + 1}\n"
              f"parameters varied: {', '.join(READOUT_ONLY)}\n"
              f"max relative change in the 8 upstream species = {worst:.3e}\n"
              f"tolerance {tol:.0e}")
    return _report("T3 readout independence", ok, detail)


# ----------------------------------------------------------------- T4

def _apply(p, drug, E):
    """Mirror of pharmacology.apply_drug for the two interventions under test,
    restated here so the test does not depend on that module."""
    q = dict(p)
    if drug == 'jnk':
        q['b_J_by_JJ'] = p['b_J_by_JJ'] * (1 - E)
        q['b_J_by_FJ'] = p['b_J_by_FJ'] * (1 - E)
        q['b_G_by_FJ'] = p['b_G_by_FJ'] * (1 - E)
    elif drug == 't5224_fos_jj':
        q['b_J_by_FJ'] = p['b_J_by_FJ'] * (1 - E)
        q['b_G_by_FJ'] = p['b_G_by_FJ'] * (1 - E)
        q['Km_rep'] = p['Km_rep'] / (1 - E)
        q['b_J_by_JJ'] = p['b_J_by_JJ'] * (1 - E)
    return q


def test_intervention_equivalence(n_steps=41):
    """The two interventions differ only in Km_rep, which enters dG alone.
    Total c-Jun must therefore be identical at every engagement."""
    p = dict(ref_params())
    w = robust_window(p)
    p['mir'] = w['m_mid']
    y_hi = np.zeros(N)
    y_hi[SPECIES.index('J')] = 100.0

    worst_cjun, worst_up, gfap_differs = 0.0, 0.0, False
    for E in np.linspace(0.0, 0.8, n_steps):
        a = steady_state(_apply(p, 'jnk', E), y0=y_hi, t_end=5000.0)
        b = steady_state(_apply(p, 't5224_fos_jj', E), y0=y_hi, t_end=5000.0)
        ca, cb = totals(a)['cJun_total'], totals(b)['cJun_total']
        worst_cjun = max(worst_cjun, abs(ca - cb))
        worst_up = max(worst_up, float(np.max(np.abs(a[UPSTREAM]
                                                     - b[UPSTREAM]))))
        if abs(a[SPECIES.index('G')] - b[SPECIES.index('G')]) > 1e-6:
            gfap_differs = True

    tol = 1e-9
    ok = worst_cjun < tol and worst_up < tol and gfap_differs
    detail = (f"engagement 0.00 to 0.80 in {n_steps} steps\n"
              f"max |cJun_total(JNK) - cJun_total(T-5224+JJ)| = "
              f"{worst_cjun:.2e}\n"
              f"max |difference| over all 8 upstream species = {worst_up:.2e}\n"
              f"tolerance {tol:.0e}\n"
              f"GFAP differs somewhere (it must, Km_rep differs): "
              f"{gfap_differs}")
    return _report("T4 intervention equivalence", ok, detail)


def main():
    print("Property tests for the AP-1 glial switch model")
    print(f"master seed {SEED}\n")
    results = [
        test_dimer_stoichiometry(),
        test_equilibrium_stability(),
        test_readout_independence(),
        test_intervention_equivalence(),
    ]
    print(f"\n{sum(results)} of {len(results)} properties hold")
    return 0 if all(results) else 1


if __name__ == '__main__':
    sys.exit(main())
