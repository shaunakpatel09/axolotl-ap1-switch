import warnings; warnings.filterwarnings('ignore')
import numpy as np, gsa
from ap1_model_v3 import steady_state, refine_fixed_point, is_stable, totals

res = gsa.load(); X = res['X']; ref = gsa.ref_params()

def branch_states(p, n_pts=81, span=10.0):
    rel = np.linspace(0, span, n_pts); mirs = rel * p['Km_mir']
    up, dn = gsa._sweeps(p, mirs)
    cu, cd = gsa._cjun(up), gsa._cjun(dn)
    hyst = np.abs(cu - cd) / (0.5 * (cu + cd)) > 0.25
    if not hyst.any(): return None
    i = np.where(hyst)[0][hyst.sum() // 2]
    q = dict(p); q['mir'] = mirs[i]
    a = steady_state(q, y0=up[i], t_end=20000.0); b = steady_state(q, y0=dn[i], t_end=20000.0)
    hi, lo = (a, b) if gsa._cjun(a[None])[0] > gsa._cjun(b[None])[0] else (b, a)
    return q, hi, lo

def stable_pair(q, hi, lo):
    ok = True
    for y in (hi, lo):
        yf, conv = refine_fixed_point(q, y)
        ok &= conv and is_stable(q, yf)
    sep = gsa._cjun(hi[None])[0] / gsa._cjun(lo[None])[0]
    return ok, sep

print('SPOT-CHECK: 8 random bistable samples of the full model')
rng = np.random.default_rng(0)
bist = np.where(res['width_full'] > 0)[0]
for k in rng.choice(bist, 8, replace=False):
    out = branch_states(gsa.params_from_row(X[k]))
    if out is None: print(f'  sample {k}: window not re-found'); continue
    ok, sep = stable_pair(*out)
    print(f'  sample {k}: two stable fixed points after 20000 h: {ok}   c-Jun {sep:.2f}x')

print('\nC1 FAILURES: bistable with the homodimer loop removed')
idx = np.where(res['width_noJJ'] > 0)[0]
rows = []
for k in idx:
    p = gsa.params_from_row(X[k]); q = dict(p); q['b_J_by_JJ'] = 0.0
    out = branch_states(q, n_pts=61)
    if out is None: continue
    qq, hi, lo = out
    ok, sep = stable_pair(qq, hi, lo)
    th, tl = totals(hi), totals(lo)
    rows.append(dict(ok=ok, sep=sep, fj_ratio=th['FJ']/tl['FJ'],
        jun_fos_lo=tl['cJun_total']/tl['cFos_total'], jun_fos_hi=th['cJun_total']/th['cFos_total'],
        r_FJ=X[k, gsa.NAMES.index('r_FJ')], a_J=X[k, gsa.NAMES.index('a_J')],
        a_F=X[k, gsa.NAMES.index('a_F')], n_auto=X[k, gsa.NAMES.index('n_auto')]))
import pandas as pd
d = pd.DataFrame(rows)
print(f'  re-found {len(d)} of {len(idx)}; verified stable pairs: {int(d.ok.sum())}')
print(d.drop(columns='ok').describe().loc[['25%','50%','75%']].round(2).to_string())
allX = pd.DataFrame(X, columns=gsa.NAMES)
print('\n  population medians: r_FJ %.2f  a_J fold %.2f  a_F fold %.2f  n_auto %.2f' %
      (allX.r_FJ.median(), allX.a_J.median(), allX.a_F.median(), allX.n_auto.median()))
d.to_csv('gsa_results/c1_failures.csv', index=False)
