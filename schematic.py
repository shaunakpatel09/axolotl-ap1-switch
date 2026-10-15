"""
schematic.py -- Figure 1: the network, and which parts of it carry the state.

WHY TWO PANELS

A pathway cartoon drawn as arrows cannot distinguish the two mechanisms this
paper separates. "JunB competes with c-Jun for c-Fos" and "JunB sequesters
c-Jun directly" are the same picture and different equations. Panel A is
therefore the topology with every edge labelled by the parameter that governs
it; panel B reports, from the model rather than by assertion, how much each
species actually differs between the two stable states.

Panel B is the point. Most of the network is static between states. Reading it
beside panel A is what makes the c-Fos ceiling argument visible rather than
verbal.

Author: Shaunak Patel. Built with AI assistance (Claude); reviewed and owned
by the author.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from ap1_model_v3 import SPECIES, steady_state, totals
from round2 import ref_params, robust_window

OUT_DIR = 'figures_round2'

INK = '#16302B'
MUTED = '#5A6B72'
GREEN = '#1D7A5F'
RUST = '#C25A32'
VIOLET = '#534AB7'
PAPER = '#FFFFFF'

# node: (x, y, label, facecolour, edgecolour)
NODES = {
    'mir':  (1.1, 8.6, 'miR-200a',  '#EDE6F2', VIOLET),
    'J':    (3.2, 6.6, 'c-Jun',     '#E3F0E9', GREEN),
    'F':    (6.0, 8.6, 'c-Fos',     '#F7E6DE', RUST),
    'B':    (8.9, 6.6, 'JunB',      '#EDE6F2', VIOLET),
    'JJ':   (1.6, 3.9, 'c-Jun:c-Jun', '#E3F0E9', GREEN),
    'FJ':   (4.7, 3.9, 'c-Fos:c-Jun', '#F7E6DE', RUST),
    'FB':   (7.9, 3.9, 'c-Fos:JunB',  '#F7E6DE', RUST),
    'G':    (4.7, 1.2, 'GFAP',      '#F1EDE4', '#8A7F6A'),
}
NODE_W = {'mir': 1.5, 'J': 1.25, 'F': 1.25, 'B': 1.2,
          'JJ': 2.0, 'FJ': 2.0, 'FB': 2.0, 'G': 1.3}


def _node(ax, key, fontsize=12.5):
    x, y, label, fc, ec = NODES[key]
    w, h = NODE_W[key], 0.62
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle='round,pad=0.10,rounding_size=0.12',
                                facecolor=fc, edgecolor=ec, linewidth=1.5,
                                zorder=3))
    ax.text(x, y, label, ha='center', va='center', fontsize=fontsize,
            color=INK, zorder=4)


def _edge(ax, p0, p1, color, style='-|>', rad=0.0, lw=1.6, ls='-', z=2):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style,
                                 connectionstyle=f'arc3,rad={rad}',
                                 color=color, linewidth=lw, linestyle=ls,
                                 mutation_scale=16, shrinkA=6, shrinkB=6,
                                 zorder=z))


def _bar_end(ax, p0, p1, color, rad=0.0, lw=1.6):
    """Repression: a line ending in a perpendicular bar."""
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle='-[,widthB=0.30,lengthB=0',
                                 connectionstyle=f'arc3,rad={rad}',
                                 color=color, linewidth=lw,
                                 mutation_scale=16, shrinkA=6, shrinkB=8,
                                 zorder=2))


def _plabel(ax, x, y, text, color=MUTED, weight='normal', size=10.5):
    ax.text(x, y, text, ha='center', va='center', fontsize=size, color=color,
            style='italic', fontweight=weight, zorder=5,
            bbox=dict(boxstyle='round,pad=0.18', facecolor=PAPER,
                      edgecolor='none', alpha=0.92))


def panel_a(ax):
    for k in NODES:
        _node(ax, k)

    xm, ym, *_ = NODES['mir']
    xJ, yJ, *_ = NODES['J']
    xF, yF, *_ = NODES['F']
    xB, yB, *_ = NODES['B']
    xJJ, yJJ, *_ = NODES['JJ']
    xFJ, yFJ, *_ = NODES['FJ']
    xFB, yFB, *_ = NODES['FB']
    xG, yG, *_ = NODES['G']

    # miR-200a represses c-Jun production
    _bar_end(ax, (xm, ym - 0.35), (xJ - 0.25, yJ + 0.35), VIOLET, rad=-0.15)
    _plabel(ax, 1.55, 7.45, 'K$_m$, n$_m$, leak', VIOLET)

    # dimerisation
    _edge(ax, (xJ - 0.45, yJ - 0.35), (xJJ + 0.25, yJJ + 0.35), GREEN,
          style='<|-|>', rad=0.12, lw=2.2)
    _plabel(ax, 1.95, 5.35, 'k$_{on}$ / k$_{off}$  (JJ)', GREEN, 'bold')

    _edge(ax, (xJ + 0.35, yJ - 0.30), (xFJ - 0.45, yFJ + 0.35), MUTED,
          style='<|-|>', rad=-0.10)
    _edge(ax, (xF - 0.40, yF - 0.35), (xFJ + 0.10, yFJ + 0.38), MUTED,
          style='<|-|>', rad=0.16)
    _plabel(ax, 4.25, 5.45, 'k$_{on}$ / k$_{off}$')

    _edge(ax, (xB - 0.25, yB - 0.35), (xFB + 0.35, yFB + 0.35), MUTED,
          style='<|-|>', rad=-0.12)
    _edge(ax, (xF + 0.40, yF - 0.35), (xFB - 0.20, yFB + 0.38), MUTED,
          style='<|-|>', rad=-0.16)
    _plabel(ax, 8.55, 5.45, 'k$_{on}$ / 3k$_{off}$')

    # autoregulatory feedback
    _edge(ax, (xJJ - 0.95, yJJ + 0.35), (xJ - 0.70, yJ - 0.35), GREEN,
          rad=0.55, lw=2.6)
    _plabel(ax, 0.70, 5.35, 'b$_{JJ}$ = 36', GREEN, 'bold', 11.5)

    _edge(ax, (xFJ - 0.55, yFJ + 0.40), (xJ + 0.15, yJ - 0.40), MUTED,
          rad=0.40, lw=1.4, ls=(0, (5, 2)))
    _plabel(ax, 3.95, 5.00, 'b$_{FJ}$ = 9')

    # GFAP
    _edge(ax, (xFJ, yFJ - 0.38), (xG, yG + 0.38), RUST, lw=2.0)
    _plabel(ax, 5.25, 2.55, 'b$_G$, K, n', RUST)
    _bar_end(ax, (xFB - 0.30, yFB - 0.38), (xG + 0.55, yG + 0.30), VIOLET,
             rad=0.18)
    _plabel(ax, 8.15, 2.95, 'K$_{rep}$, n$_{rep}$', VIOLET)

    # production and degradation markers
    for key, txt in (('J', 'a$_J$'), ('F', 'a$_F$'), ('B', 'a$_B$')):
        x, y, *_ = NODES[key]
        ax.text(x, y + 0.60, txt, ha='center', va='bottom', fontsize=10.5,
                color=MUTED, style='italic')
    for key, txt in (('J', r'$\gamma_J$'), ('F', r'$\gamma_F$'),
                     ('B', r'$\gamma_B$')):
        x, y, *_ = NODES[key]
        ax.text(x + NODE_W[key] / 2 + 0.22, y - 0.05, txt, ha='left',
                va='center', fontsize=10.5, color=MUTED, style='italic')

    # c-Fos scarcity note
    ax.annotate('limiting pool', xy=(xF + 0.72, yF + 0.10),
                xytext=(xF + 1.75, yF + 0.95), fontsize=11, color=RUST,
                ha='left', va='center',
                arrowprops=dict(arrowstyle='->', color=RUST, lw=1.3))

    ax.text(6.55, 1.05,
            'solid arrow = activation\n'
            'bar = repression\n'
            'double arrow = reversible binding\n'
            'bold green = the loop that carries the state',
            fontsize=10.5, color=MUTED, ha='left', va='center',
            linespacing=1.6)

    ax.set_xlim(-0.2, 11.0)
    ax.set_ylim(-0.3, 9.9)
    ax.axis('off')
    ax.set_title('A  Network topology, with the parameter governing each edge',
                 loc='left', fontsize=13, color=INK, pad=12)


def state_separations(p=None):
    """Fold difference between the two stable states, per species, computed at
    the midpoint of the two-state window."""
    p = p or ref_params()
    w = robust_window(p)
    q = dict(p)
    q['mir'] = w['m_mid']
    y_hi = np.zeros(len(SPECIES))
    y_hi[SPECIES.index('J')] = 100.0
    hi = steady_state(q, y0=y_hi, t_end=5000.0)
    lo = steady_state(q, y0=np.zeros(len(SPECIES)), t_end=5000.0)
    if totals(hi)['cJun_total'] < totals(lo)['cJun_total']:
        hi, lo = lo, hi
    out = {}
    for i, s in enumerate(SPECIES):
        out[s] = hi[i] / max(lo[i], 1e-12)
    out['cJun_total'] = (totals(hi)['cJun_total']
                         / max(totals(lo)['cJun_total'], 1e-12))
    return out, w


def panel_b(ax, sep):
    order = [('JJ', 'c-Jun:c-Jun', GREEN),
             ('J', 'c-Jun (free)', GREEN),
             ('cJun_total', 'c-Jun (total)', GREEN),
             ('JB', 'c-Jun:JunB', MUTED),
             ('FJ', 'c-Fos:c-Jun', RUST),
             ('G', 'GFAP', '#8A7F6A'),
             ('F', 'c-Fos (free)', RUST),
             ('FB', 'c-Fos:JunB', RUST),
             ('B', 'JunB (free)', VIOLET)]
    labels = [lab for _, lab, _ in order]
    vals = [sep[k] for k, _, _ in order]
    cols = [c for _, _, c in order]
    pos = np.arange(len(order))

    ax.barh(pos, vals, color=cols, alpha=0.85, height=0.62)
    ax.axvline(1.0, color=INK, lw=1.2)
    ax.set_yticks(pos)
    ax.set_yticklabels(labels, fontsize=11.5)
    ax.invert_yaxis()
    ax.set_xscale('log')
    ax.set_xlabel('fold difference between the two stable states  '
                  '(log scale; 1 = identical)', fontsize=11.5)
    for i, v in enumerate(vals):
        ax.text(v * 1.25 if v >= 1 else v * 0.78, i,
                f'{v:,.2f}×' if v < 100 else f'{v:,.0f}×',
                va='center', ha='left' if v >= 1 else 'right',
                fontsize=11, color=INK)
    ax.set_xlim(min(vals) * 0.35, max(vals) * 6)
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis='y', length=0)
    ax.set_title('B  How much each species differs between the two states',
                 loc='left', fontsize=13, color=INK, pad=12)


def figure1(p=None, show_window=True):
    sep, w = state_separations(p)
    fig, ax = plt.subplots(1, 2, figsize=(19, 8.2),
                           gridspec_kw=dict(width_ratios=[1.32, 1.0]))
    panel_a(ax[0])
    panel_b(ax[1], sep)

    ratio = sep['JJ'] / sep['FJ']
    note = (f'Most of the network is unchanged between states. The c-Jun '
            f'homodimer separates {sep["JJ"]:.1f}-fold while c-Fos:c-Jun, '
            f'pinned near the limiting c-Fos pool, separates only '
            f'{sep["FJ"]:.2f}-fold — {ratio:.0f}x less state information.\n'
            f'That gap is why the heterodimer cannot carry the switch, and '
            f'why any gene reading only c-Fos:c-Jun inherits a narrow dynamic '
            f'range. Free c-Fos falls to {sep["F"]:.2f}x in the high-c-Jun '
            f'state: it is being consumed, not accumulating.')
    if show_window:
        note += (f'\nComputed at the reference parameter set, at '
                 f'miR-200a = {w["m_mid"]:.2f} (the midpoint of its two-state '
                 f'window). The larger separations quoted in the text '
                 f'(FJ 1.46x, JJ 736x, medians) come from the 27 bistable '
                 f'sets of the cfos_ceiling_check.py scan, evaluated at '
                 f'miR-200a = 0 where both states lie further apart and '
                 f'where this reference set is itself monostable. Those are '
                 f'not the 400-sample sensitivity analysis, which did not '
                 f'record per-dimer separations. The ordering is the same in '
                 f'both.')
    fig.text(0.012, 0.012, note, fontsize=11.5, color=MUTED, ha='left',
             va='bottom', wrap=True)

    fig.suptitle('Figure 1. The AP-1 glial switch: what is connected, and '
                 'what actually changes', fontsize=15, color=INK, y=0.985)
    fig.tight_layout(rect=(0, 0.155, 1, 0.96))
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, 'fig1_schematic.png')
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    print('saved', path)
    for k, v in sep.items():
        print(f'  {k:12s} {v:12.3f}x')
    return fig, sep


if __name__ == '__main__':
    figure1()
