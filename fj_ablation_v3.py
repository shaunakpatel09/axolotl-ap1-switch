"""
Re-test of the lead finding under v3: does routing c-Jun autoregulation through
cFos:cJun (FJ) SUPPRESS multistability, as the v2 decomposition found?

Matched design: for every LHS-sampled v3 parameter set, run four conditions
using the SAME feedback strength b (the sampled b_J_by_JJ), routed differently:
  JJ only   : b_J_by_JJ = b, b_J_by_FJ = 0   (v3 default)
  JJ + FJ   : b_J_by_JJ = b, b_J_by_FJ = b
  FJ only   : b_J_by_JJ = 0, b_J_by_FJ = b
  neither   : both 0
Caveat: FJ and JJ sit at different concentrations but share Km_auto, so equal b
is equal maximal strength, not equal effective activation.
"""
import sys, numpy as np
from ap1_model_v3 import find_attractors, sample_initial_conditions
from ap1_v3_downstream import v3_sample_parameter_sets

def run(conds, n_sets=150, n_ic=8, seed=42):
    rng = np.random.default_rng(seed)
    psets = v3_sample_parameter_sets(n_sets, rng)
    y0s = sample_initial_conditions(n_ic, rng)
    out = {c: [] for c in conds}
    for p0 in psets:
        b = p0['b_J_by_JJ']
        spec = {'JJ only': (b, 0.0), 'JJ + FJ': (b, b),
                'FJ only': (0.0, b), 'neither': (0.0, 0.0)}
        for c in conds:
            p = dict(p0); p['b_J_by_JJ'], p['b_J_by_FJ'] = spec[c]
            n, _ = find_attractors(p, y0s)
            out[c].append(n)
    for c in conds:
        a = np.array(out[c]); v = a[a > 0]
        print(f'{c:9s} multistable: {100*np.mean(v>1):5.1f}%  (n={len(v)})', flush=True)
    return out

if __name__ == '__main__':
    run(sys.argv[1].split(','))
