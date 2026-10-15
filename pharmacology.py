"""
pharmacology.py -- drug interventions on the round-2 reference AP-1 model.

TEST 1 (this file so far): which target needs the smallest push?
No pharmacokinetics yet. Each drug is applied as a SUSTAINED effect, and we
find the smallest effect that removes the high-c-Jun state. Then the drug is
stopped, to see whether the cell stays switched.

Drugs are represented by what they do to existing model parameters. Nothing in
the v3 equations changes.

  mimic            miR-200a mimic: adds miR-200a. Engagement is the fraction of
                   c-Jun production silenced relative to the cell's baseline,
                   so it shares an axis with the inhibitors. The model's
                   repression leak (0.15) caps how much a mimic can silence.
  jnk              JNK inhibitor, activity only. JNK phosphorylation of c-Jun
                   Ser63/73 strengthens its transactivation and feeds c-Jun
                   autoregulation (Hibi et al. 1993; Angel & Karin 1991).
                   Engagement E scales c-Jun self-activation (homodimer and
                   c-Fos:c-Jun edges) and c-Fos:c-Jun activation of GFAP by
                   (1 - E).
  jnk_stab         JNK inhibitor, activity + stability. As above, plus c-Jun
                   degradation multiplied by (1 + E): Ser73 phosphorylation
                   protects c-Jun from ubiquitination (Fuchs et al. 1996;
                   Musti et al. 1997). JNK also helps TARGET c-Jun for
                   ubiquitination, so the true size of this effect is
                   uncertain; the two JNK variants bracket it. The factor
                   (1 + E), i.e. at most a doubling, is a placeholder.
  t5224_fos        T-5224 blocking DNA binding of c-Fos-containing dimers
                   (c-Fos:c-Jun and c-Fos:JunB): scales the c-Fos:c-Jun ->
                   c-Jun edge and c-Fos:c-Jun activation of GFAP by (1 - E),
                   and weakens c-Fos:JunB repression of GFAP by the same factor.
  t5224_fos_jj     As above, plus the c-Jun homodimer loop. The literature
                   describes T-5224 as targeting c-Fos/c-Jun but does not rule
                   homodimers in or out, so both mappings are modelled.

Engagement E is the fraction of the target's activity blocked (0 to 1). Equal
engagement does NOT mean equal dose: converting engagement to dose needs
potency and pharmacokinetics, which come in the next test.

Author: Shaunak Patel. Built with AI assistance (Claude); reviewed and owned
by the author.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from ap1_model_v3 import N, SPECIES, steady_state, totals, repress
from round2 import ref_params, robust_window

_J = SPECIES.index('J')
OUT_DIR = 'figures_round2'

DRUGS = {
    'mimic': 'miR-200a mimic',
    'jnk': 'JNK inhibitor (activity)',
    'jnk_stab': 'JNK inhibitor (activity + stability)',
    't5224_fos': 'T-5224 (c-Fos dimers)',
    't5224_fos_jj': 'T-5224 (c-Fos dimers + homodimer)',
}
COLORS = {'mimic': '#1D9E75', 'jnk': '#534AB7', 'jnk_stab': '#8E86D8',
          't5224_fos': '#D85A30', 't5224_fos_jj': '#BA7517'}


# ----------------------------------------------------------------------
# How each drug changes the model
# ----------------------------------------------------------------------

def mir_rep(p, m):
    """Fraction of maximal c-Jun production left at miR-200a level m
    (same expression as in rhs)."""
    leak = p.get('mir_leak', 0.0)
    return leak + (1.0 - leak) * repress(m, p['Km_mir'], p['n_mir'])


def max_engagement(p, drug, m0):
    """Largest achievable engagement. For the mimic this is set by the leak:
    production can never fall below leak x maximal."""
    if drug == 'mimic':
        return 1.0 - p['mir_leak'] / mir_rep(p, m0)
    return 1.0


def mir_for_silencing(p, m0, s):
    """miR-200a level that silences a fraction s of c-Jun production relative
    to baseline m0. Inverts mir_rep: needs s below max_engagement."""
    leak = p['mir_leak']
    r = (mir_rep(p, m0) * (1.0 - s) - leak) / (1.0 - leak)   # Hill term
    return p['Km_mir'] * (1.0 / r - 1.0) ** (1.0 / p['n_mir'])


def apply_drug(p, drug, E, m0):
    """Parameters for a cell at baseline miR-200a m0 under sustained drug
    engagement E."""
    q = dict(p)
    q['mir'] = m0
    if drug == 'mimic':
        q['mir'] = mir_for_silencing(p, m0, E)
        return q
    E = min(max(E, 0.0), 0.999)
    if drug in ('jnk', 'jnk_stab'):
        q['b_J_by_JJ'] = p['b_J_by_JJ'] * (1 - E)
        q['b_J_by_FJ'] = p['b_J_by_FJ'] * (1 - E)
        q['b_G_by_FJ'] = p['b_G_by_FJ'] * (1 - E)
        if drug == 'jnk_stab':
            q['g_J'] = p['g_J'] * (1 + E)
    elif drug in ('t5224_fos', 't5224_fos_jj'):
        q['b_J_by_FJ'] = p['b_J_by_FJ'] * (1 - E)
        q['b_G_by_FJ'] = p['b_G_by_FJ'] * (1 - E)
        # blocking a fraction E of c-Fos:JunB means repression by FB*(1-E),
        # which is the same as raising the repression constant to Km_rep/(1-E)
        q['Km_rep'] = p['Km_rep'] / (1 - E)
        if drug == 't5224_fos_jj':
            q['b_J_by_JJ'] = p['b_J_by_JJ'] * (1 - E)
    else:
        raise ValueError(drug)
    return q


# ----------------------------------------------------------------------
# Test 1
# ----------------------------------------------------------------------

def high_state(p, m0):
    """Untreated high-c-Jun steady state at baseline miR-200a m0, reached by
    continuation from miR-200a = 0 starting at high c-Jun."""
    q = dict(p)
    y = np.zeros(N)
    y[_J] = 100.0
    for m in np.linspace(0.0, m0, 21):
        q['mir'] = m
        y = steady_state(q, y0=y, t_end=800.0)
    return y


def low_state(p, m0):
    """Untreated low-c-Jun steady state at m0 (continuation from high miR-200a
    starting at zero). Equals the high state outside the two-state window."""
    q = dict(p)
    y = np.zeros(N)
    for m in np.linspace(2.5, m0, 41):
        q['mir'] = m
        y = steady_state(q, y0=y, t_end=800.0)
    return y


def engagement_scan(p, m0, drug, n=201):
    """Continuation in drug engagement, starting from the untreated high-c-Jun
    state at m0. Returns (engagement grid, list of states)."""
    y = high_state(p, m0)
    grid = np.linspace(0.0, 0.999 * max_engagement(p, drug, m0), n)
    states = []
    for E in grid:
        y = steady_state(apply_drug(p, drug, E, m0), y0=y, t_end=1500.0)
        states.append(y.copy())
    return grid, states


def required_engagement(grid, states, drop=0.25):
    """Smallest engagement at which the high-c-Jun state disappears: the first
    step where total c-Jun falls by more than `drop` in one step (a jump, not a
    graded decline). Returns (engagement, index after the jump), or
    (nan, None) if no jump occurs."""
    cj = np.array([totals(y)['cJun_total'] for y in states])
    jumps = np.where((cj[:-1] - cj[1:]) / cj[:-1] > drop)[0]
    if not len(jumps):
        return np.nan, None
    return grid[jumps[0] + 1], jumps[0] + 1


def after_withdrawal(p, m0, y, t_end=5000.0):
    """Stop the drug: integrate the untreated model at m0 from state y."""
    q = dict(p)
    q['mir'] = m0
    return steady_state(q, y0=y, t_end=t_end)


def test1(baselines=(0.0, 0.15, 0.3, 0.45, 0.58, 0.7), p=None, n=201):
    """For every drug at every baseline miR-200a: the sustained engagement
    needed to remove the high-c-Jun state, GFAP under the drug just past that
    point, and whether the cell stays low once the drug is stopped."""
    p = p or ref_params()
    w = robust_window(p)
    print(f"reference two-state window: miR-200a {w['lo']:.2f} to {w['hi']:.2f}"
          " (inside it, a switch made by the drug can persist after it stops)\n")

    hi = {m0: totals(high_state(p, m0)) for m0 in baselines}
    lo = {m0: totals(low_state(p, m0)) for m0 in baselines}
    scans, need, gfap_on, stays = {}, {}, {}, {}
    for drug in DRUGS:
        need[drug], gfap_on[drug], stays[drug] = [], [], []
        for m0 in baselines:
            g, states = engagement_scan(p, m0, drug, n)
            scans[(drug, m0)] = (g, states)
            E, k = required_engagement(g, states)
            need[drug].append(E)
            if k is None:
                gfap_on[drug].append(np.nan)
                stays[drug].append(None)
                continue
            gfap_on[drug].append(totals(states[k])['GFAP'] / hi[m0]['GFAP'])
            back = totals(after_withdrawal(p, m0, states[k]))
            stays[drug].append(back['cJun_total'] < 0.75 * hi[m0]['cJun_total'])

    head = 'baseline miR-200a'
    cols = ''.join(f'{m0:>8.2f}' for m0 in baselines)

    def show(title, table, fmt):
        print(title)
        print(f'  {head:40s}{cols}')
        for drug, name in DRUGS.items():
            print(f'  {name:40s}' + ''.join(fmt(v) for v in table[drug]))
        print()

    show('ENGAGEMENT NEEDED to remove the high-c-Jun state (fraction blocked)',
         need, lambda v: '   never' if np.isnan(v) else f'{v:8.2f}')
    print('  mimic ceiling set by the leak: ' + ' '.join(
        f'{max_engagement(p, "mimic", m0):.2f}' for m0 in baselines) + '\n')
    show('GFAP UNDER THE DRUG just past that point (x untreated high state)',
         gfap_on, lambda v: '       -' if np.isnan(v) else f'{v:8.2f}')
    show('STAYS LOW AFTER THE DRUG IS STOPPED?',
         stays, lambda v: '       -' if v is None else ('     yes' if v else '      no'))
    print('untreated GFAP, low / high state: ' + ' '.join(
        f'{lo[m0]["GFAP"] / hi[m0]["GFAP"]:.2f}' for m0 in baselines))
    return dict(p=p, window=w, baselines=baselines, need=need,
                gfap_on=gfap_on, stays=stays, scans=scans, hi=hi, lo=lo)


def figure8(res, m_detail=0.58):
    p, w, baselines = res['p'], res['window'], res['baselines']
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))

    a = ax[0]
    a.axvspan(w['lo'], w['hi'], color='purple', alpha=0.12,
              label='two-state window')
    for drug, name in DRUGS.items():
        v = np.array(res['need'][drug], float)
        a.plot(baselines, v, 'o-', color=COLORS[drug], label=name)
        never = np.isnan(v)
        a.plot(np.array(baselines)[never], np.full(never.sum(), 1.03), 'x',
               color=COLORS[drug], ms=9, mew=2)
    a.set_ylim(0, 1.08)
    a.set_xlabel('baseline miR-200a (multiples of its repression constant)')
    a.set_ylabel('target engagement needed (fraction blocked)')
    a.set_title('A  Sustained engagement needed to remove the\n'
                'high-c-Jun state (x at top = never)', loc='left', fontsize=10)
    a.legend(fontsize=7)

    if m_detail in baselines:
        scans = {d: res['scans'][(d, m_detail)] for d in DRUGS}
    else:
        scans = {d: engagement_scan(p, m_detail, d) for d in DRUGS}
    ref_hi = totals(high_state(p, m_detail))
    ref_lo = totals(low_state(p, m_detail))
    for k, (key, ylab, title) in enumerate((
            ('cJun_total', 'total c-Jun (nM)', 'B  c-Jun as engagement rises'),
            ('GFAP', 'GFAP (nM)', 'C  GFAP as engagement rises'))):
        b = ax[k + 1]
        for drug, name in DRUGS.items():
            g, states = scans[drug]
            b.plot(g, [totals(y)[key] for y in states], color=COLORS[drug],
                   lw=2, label=name)
        b.axhline(ref_hi[key], color='0.4', ls=':', lw=1)
        b.axhline(ref_lo[key], color='0.4', ls='--', lw=1)
        b.set_xlabel('target engagement (fraction blocked)')
        b.set_ylabel(ylab)
        b.set_title(f'{title}\n(baseline miR-200a {m_detail:.2f}; dotted = '
                    'untreated high, dashed = untreated low)', loc='left',
                    fontsize=10)
    ax[2].legend(fontsize=7)
    fig.suptitle('Figure 8. Test 1: which drug target needs the smallest '
                 'sustained push to remove the high-c-Jun state', fontsize=12)
    fig.tight_layout()
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, 'fig8_target_ranking.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    print('saved', path)
    return fig


if __name__ == '__main__':
    figure8(test1())
