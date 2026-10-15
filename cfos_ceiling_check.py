"""
Does the c-Fos ceiling explanation generalise?

The corrected mechanism (see PREDICTIONS_v3.md, "Corrections") says the
cFos:cJun heterodimer cannot carry the autoregulatory feedback because c-Fos is
limiting: FJ is pinned near the c-Fos ceiling in BOTH stable states, so it
barely differs between them and cannot sense the cell's state. The c-Jun
homodimer is not capped by c-Fos and does carry the signal.

That was checked at the reference parameter set only. This script tests it
across every bistable parameter set found by an LHS scan.

For each bistable set, at the two stable states (ordered by total c-Jun):
  fj_ratio      FJ(ON) / FJ(OFF)                  -- state information in FJ
  jj_ratio      JJ(ON) / JJ(OFF)                  -- state information in JJ
  fj_frac       FJ / total c-Fos, in each state   -- how close to the ceiling
  hill_on/off   the feedback term driven by FJ    -- saturation

Mechanism holds if fj_ratio stays near 1 while jj_ratio is much larger, and
fj_frac is high in both states.
"""
import numpy as np
import pandas as pd

from ap1_model_v3 import (SPECIES, hill, totals, find_attractors,
                          cluster_states, sample_initial_conditions)
from ap1_v3_downstream import v3_sample_parameter_sets


def ceiling_check(n_sets=150, n_ic=10, seed=7, verbose=True):
    rng = np.random.default_rng(seed)
    psets = v3_sample_parameter_sets(n_sets, rng)
    y0s = sample_initial_conditions(n_ic, rng)
    rows = []
    for i, p in enumerate(psets):
        n, states = find_attractors(p, y0s)
        if n < 2:
            continue
        labels, k = cluster_states(states)
        reps = [states[np.where(labels == c)[0][0]] for c in range(k)]
        reps.sort(key=lambda y: totals(y)['cJun_total'])
        off, on = reps[0], reps[-1]
        t_on, t_off = totals(on), totals(off)
        iJJ = SPECIES.index('JJ')
        rows.append(dict(
            cjun_ratio=t_on['cJun_total'] / max(t_off['cJun_total'], 1e-12),
            fj_ratio=t_on['FJ'] / max(t_off['FJ'], 1e-12),
            jj_ratio=on[iJJ] / max(off[iJJ], 1e-12),
            fj_frac_on=t_on['FJ'] / max(t_on['cFos_total'], 1e-12),
            fj_frac_off=t_off['FJ'] / max(t_off['cFos_total'], 1e-12),
            hill_on=hill(t_on['FJ'], p['Km_auto'], p['n_auto']),
            hill_off=hill(t_off['FJ'], p['Km_auto'], p['n_auto']),
        ))
        if verbose and len(rows) % 10 == 0:
            print(f'  {len(rows)} bistable sets found ({i+1}/{n_sets} scanned)',
                  flush=True)
    return pd.DataFrame(rows)


def report(df):
    if df.empty:
        print('no bistable parameter sets found')
        return
    print(f'\nbistable parameter sets: {len(df)}\n')
    q = df.quantile([0.25, 0.5, 0.75]).T
    q.columns = ['25%', 'median', '75%']
    print(q.round(3).to_string())
    print('\ninformation ratio (jj_ratio / fj_ratio), median: '
          f'{(df.jj_ratio / df.fj_ratio).median():.1f}x')
    print(f'sets where FJ differs < 1.5x between states: '
          f'{100*(df.fj_ratio < 1.5).mean():.0f}%')
    print(f'sets where JJ differs > 3x between states:   '
          f'{100*(df.jj_ratio > 3).mean():.0f}%')
    print(f'sets where FJ holds > 70% of c-Fos in BOTH states: '
          f'{100*((df.fj_frac_on > 0.7) & (df.fj_frac_off > 0.7)).mean():.0f}%')


if __name__ == '__main__':
    report(ceiling_check())
