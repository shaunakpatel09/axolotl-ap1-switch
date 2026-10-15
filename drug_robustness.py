"""
drug_robustness.py -- does the Test-1 target ranking depend on the parameters,
or on the wiring?

WHY THIS MATTERS MORE THAN PHARMACOKINETICS
-------------------------------------------
Test 1 ranked five interventions at ONE reference parameter set -- a set found
by grid search, not measurement. A ranking that only holds there is a fact
about that guess. A ranking that holds across every bistable parameter set the
global sensitivity analysis found is a consequence of the network's wiring,
and wiring is the part that is conserved between species. That is the only
honest route from this model to a claim about mammalian cells, since the
model's miR-200a axis is expressed in multiples of a repression constant that
no expression dataset can supply.

So: take the bistable parameter sets from gsa.py, put every drug through each
one, and ask whether the ordering survives.

Reference ordering (round-2 reference set, baseline miR-200a 0.58):
    JNK (activity + stability)  0.11
    miR-200a mimic              0.16
    JNK (activity)              0.19   ) identical in c-Jun by construction:
    T-5224 (c-Fos + homodimer)  0.19   ) they differ only in the GFAP terms
    T-5224 (c-Fos dimers)       0.82

The two tied entries are a built-in check. jnk and t5224_fos_jj scale the same
three c-Jun parameters and differ only in Km_rep, which enters the GFAP
equation alone. If they ever disagree on c-Jun, something is wrong.

Each sample is evaluated at the MIDPOINT of its own two-state window, since
every sample has a different Km_mir and therefore a different miR-200a scale.

Author: Shaunak Patel. Built with AI assistance (Claude); reviewed and owned
by the author.
"""

import os

import numpy as np

import gsa
from ap1_model_v3 import N, SPECIES, steady_state, totals
from pharmacology import DRUGS, COLORS, apply_drug, max_engagement
from regime import drug_in_regime

_J = SPECIES.index('J')
RESULT_DIR = 'drug_robustness_results'
OUT_DIR = 'figures_round2'

# the ordering found at the reference set, cheapest first
REF_ORDER = ['jnk_stab', 'mimic', 'jnk', 't5224_fos_jj', 't5224_fos']


def bistable_from_load():
    """(row, params) for every GSA sample with a two-state window. `row` is
    the position in gsa.load()'s arrays, not the original sample number."""
    res = gsa.load()
    keep = np.where(res['width_full'] > 0)[0]
    return [(int(k), gsa.params_from_row(res['X'][k])) for k in keep]


def window_mid(p, n_pts=61, span=10.0):
    """Midpoint of this parameter set's two-state window, in ABSOLUTE
    miR-200a units. None if monostable. Mirrors gsa.window()'s detector."""
    rel = np.linspace(0, span, n_pts)
    mirs = rel * p['Km_mir']
    up, dn = gsa._sweeps(p, mirs)
    cu, cd = gsa._cjun(up), gsa._cjun(dn)
    hyst = np.abs(cu - cd) / (0.5 * (cu + cd)) > 0.25
    if not hyst.any():
        return None
    return float(mirs[np.where(hyst)[0][hyst.sum() // 2]])


def sweep_one(p, m0, n=121):
    """Required engagement for each drug at this parameter set."""
    row = {}
    for drug in DRUGS:
        r = drug_in_regime(p, m0, drug, n=n)
        row[drug] = (r['E_jump'], r['E_half'], r['stays'])
    return row


def run(start=0, stop=None, tag='dr'):
    """Evaluate bistable samples start..stop-1 and save."""
    sets = bistable_from_load()
    stop = len(sets) if stop is None else min(stop, len(sets))
    os.makedirs(RESULT_DIR, exist_ok=True)
    rows = []
    for k in range(start, stop):
        gi, p = sets[k]
        m0 = window_mid(p)
        if m0 is None:
            continue
        r = sweep_one(p, m0)
        rows.append([gi, m0] + [v for d in DRUGS for v in
                               (r[d][0], r[d][1],
                                np.nan if r[d][2] is None else float(r[d][2]))])
        print(f'  {k + 1}/{stop}  sample {gi}', flush=True)
    path = os.path.join(RESULT_DIR, f'{tag}_{start}_{stop}.npz')
    np.savez(path, rows=np.array(rows, float), drugs=np.array(list(DRUGS)))
    print('saved', path, flush=True)
    return path


def load(tag='dr'):
    rows = [np.load(os.path.join(RESULT_DIR, f))['rows']
            for f in sorted(os.listdir(RESULT_DIR))
            if f.startswith(tag + '_') and f.endswith('.npz')]
    r = np.vstack([x for x in rows if len(x)])
    cols = list(DRUGS)
    jump = {d: r[:, 2 + 3 * i] for i, d in enumerate(cols)}
    half = {d: r[:, 3 + 3 * i] for i, d in enumerate(cols)}
    stays = {d: r[:, 4 + 3 * i] for i, d in enumerate(cols)}
    return dict(idx=r[:, 0].astype(int), m0=r[:, 1], jump=jump, half=half,
                stays=stays, n=len(r))


def report(res=None, metric='jump'):
    """metric='jump': engagement at which the high state ceases to exist (the
    bifurcation; nan if it never does). metric='half': engagement at which
    c-Jun is halved -- defined more often, but conflates losing the state with
    merely suppressing it. The ranking should be read off 'jump'; 'half' is
    the robustness check."""
    res = res or load()
    n = res['n']
    half = res[metric]
    print(f'{n} bistable parameter sets, each at the midpoint of its own '
          f'two-state window   [metric: {metric}]\n')

    print('CONSISTENCY CHECK  jnk and t5224_fos_jj must agree on c-Jun:')
    d = np.abs(np.nan_to_num(half['jnk'], nan=-1)
               - np.nan_to_num(half['t5224_fos_jj'], nan=-1))
    print(f'  max |difference| = {d.max():.2e}'
          f'   ({"PASS" if d.max() < 1e-9 else "FAIL"})\n')

    print(f'{"drug":42s}{"clears":>8s}{"median E":>10s}{"IQR":>16s}'
          f'{"cheapest":>10s}')
    M = np.vstack([half[d] for d in DRUGS])
    finite = np.isfinite(M)
    best = np.full(n, -1)
    for k in range(n):
        col = M[:, k]
        if np.isfinite(col).any():
            best[k] = int(np.nanargmin(np.where(np.isfinite(col), col, np.inf)))
    for i, (drug, name) in enumerate(DRUGS.items()):
        v = half[drug][np.isfinite(half[drug])]
        q1, q3 = (np.percentile(v, [25, 75]) if len(v) else (np.nan, np.nan))
        print(f'  {name:40s}{finite[i].sum():>8d}'
              f'{(np.median(v) if len(v) else np.nan):>10.2f}'
              f'{f"{q1:.2f}-{q3:.2f}":>16s}'
              f'{(best == i).sum():>10d}')

    print('\nPAIRWISE ORDERINGS (fraction of sets where the reference '
          'ordering holds;\n  only sets where BOTH drugs clear the state '
          'are counted):')
    pairs = [('jnk_stab', 'jnk'), ('jnk_stab', 'mimic'), ('mimic', 'jnk'),
             ('jnk', 't5224_fos'), ('mimic', 't5224_fos'),
             ('jnk_stab', 't5224_fos')]
    for a, b in pairs:
        ok = np.isfinite(half[a]) & np.isfinite(half[b])
        frac = (half[a][ok] < half[b][ok]).mean() if ok.sum() else np.nan
        print(f'  {DRUGS[a]:38s} < {DRUGS[b]:34s}'
              f'{frac:7.2f}   (n={ok.sum()})')

    print('\nT-5224 SPARING THE HOMODIMER:')
    f_fos = np.isfinite(half['t5224_fos'])
    print(f'  clears the high state in {f_fos.sum()}/{n} sets '
          f'({100 * f_fos.mean():.0f}%)')
    both = f_fos & np.isfinite(half['t5224_fos_jj'])
    if both.sum():
        ratio = half['t5224_fos'][both] / half['t5224_fos_jj'][both]
        print(f'  where both work, it costs {np.median(ratio):.1f}x the '
              'engagement of the homodimer-including version '
              f'(median, n={both.sum()})')
    return res


def figure11(res=None, metric='jump'):
    import matplotlib.pyplot as plt
    res = res or load()
    half = res[metric]
    fig, ax = plt.subplots(1, 2, figsize=(14.5, 5.2))

    a = ax[0]
    data = [half[d][np.isfinite(half[d])] for d in DRUGS]
    bp = a.boxplot(data, vert=True, patch_artist=True, widths=0.6,
                   medianprops=dict(color='k', lw=1.6))
    for patch, d in zip(bp['boxes'], DRUGS):
        patch.set_facecolor(COLORS[d])
        patch.set_alpha(0.75)
    short = ['miR-200a\nmimic', 'JNK\n(activity)', 'JNK\n(+ stability)',
             'T-5224\n(c-Fos)', 'T-5224\n(+ homodimer)']
    a.set_xticklabels(short, fontsize=8)
    a.set_ylabel('engagement needed (fraction blocked)')
    a.set_title('A  Engagement across every bistable parameter set\n'
                '(box = IQR; only sets where the drug works)',
                loc='left', fontsize=10)
    for i, d in enumerate(DRUGS):
        frac = np.isfinite(half[d]).mean()
        a.text(i + 1, 1.04, f'{100 * frac:.0f}%', ha='center', fontsize=7.5,
               color=COLORS[d])
    a.text(0.5, 1.09, 'percentage of sets in which the drug clears the state',
           fontsize=7.5, color='0.35', transform=a.get_yaxis_transform(),
           ha='left')
    a.set_ylim(0, 1.13)

    b = ax[1]
    pairs = [('jnk_stab', 'jnk'), ('jnk_stab', 'mimic'), ('mimic', 'jnk'),
             ('jnk_stab', 't5224_fos'), ('mimic', 't5224_fos'),
             ('jnk', 't5224_fos')]
    tag = {'mimic': 'miR-200a mimic', 'jnk': 'JNK (activity)',
           'jnk_stab': 'JNK (+ stability)', 't5224_fos': 'T-5224 (c-Fos only)',
           't5224_fos_jj': 'T-5224 (+ homodimer)'}
    labs, fracs = [], []
    for x, y in pairs:
        ok = np.isfinite(half[x]) & np.isfinite(half[y])
        fracs.append((half[x][ok] < half[y][ok]).mean() if ok.sum() else np.nan)
        labs.append(f'{tag[x]}  <\n{tag[y]}      (n={ok.sum()})')
    pos = np.arange(len(pairs))
    b.barh(pos, fracs, color='#534AB7', alpha=0.8)
    b.axvline(0.5, color='0.4', ls='--', lw=1, label='no better than a coin')
    b.set_yticks(pos)
    b.set_yticklabels(labs, fontsize=7.5)
    b.set_xlim(0, 1)
    b.set_xlabel('fraction of parameter sets where the ordering holds')
    b.set_title('B  Does the reference ranking survive the\n'
                'parameter space?', loc='left', fontsize=10)
    b.legend(fontsize=7.5, loc='lower right')

    fig.suptitle('Figure 11. The target ranking across the bistable parameter '
                 'space: wiring or guesswork?', fontsize=12)
    fig.tight_layout()
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, 'fig11_target_robustness.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print('saved', path)
    return fig
