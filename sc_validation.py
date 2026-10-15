"""
Single-cell consistency check of the v3 model against Lust et al. 2022
(ependymoglia, axolotl pallium; Zenodo 6390083, file ep_wpi_srat.rds, extracted
with extract_ap1_genes.R).

Tests
-----
P2: per-cell GFAP differs only modestly (~1.22x) between c-Jun-high and
    c-Jun-low glia.
P6: c-Jun production exceeds c-Fos (~1.2-3x).

Findings (first run)
--------------------
- FOS has exactly zero counts in all 7,343 nuclei, while FOSB/FOSL1/FOSL2 are
  detected. Likely a gene-model/annotation issue; either way P6 is UNTESTABLE
  here. Do NOT read "JUN detected, FOS not" as support for P6.
- P2 at 1 wpi: GFAP 1.47x [1.16-1.86] in JUN+ vs JUN- nuclei (depth-adjusted).
  Direction matches the model. Dropout biases this toward 1, so it is a lower
  bound; the model's 1.22x is probably too small.
- Uninjured: no JUN-associated difference in GFAP, vimentin or CSPG4.
- CSPG4 does not track JUN even at 1 wpi, contrary to the ~4.7x per-cell value
  back-calculated from Sabin 2019.

Caveats
-------
- Earliest injured sample is 1 week; Sabin's switch occurs at 1-4 days.
- Pallium, not spinal cord. mRNA, not protein.
- Every injured nucleus is Div-seq (EdU-selected); every uninjured nucleus is
  10x v3.1 or multiome. Injured-vs-uninjured comparisons are confounded with
  technology, so all tests here are WITHIN condition.
- VIM and CSPG5 contain non-integer values (likely merged gene models); they are
  rounded for the count model.

Author: Shaunak Patel
Analysis built with AI assistance (Claude); reviewed and owned by the author.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm

CSV_PATH = 'ep_wpi_srat_ap1_genes.csv'
ORDER = ['SS', '1_wpi', '2_wpi', '4_wpi', '6_wpi', '8_wpi', '12_wpi']


def load(path=CSV_PATH):
    return pd.read_csv(path)


def data_overview(df):
    """Technology confound, sequencing depth, and detection rates."""
    print('technology by condition (SS = uninjured):')
    print(pd.crosstab(df.sample_simp, df.sample_dat).reindex(ORDER).to_string())
    print('\nmedian UMIs per nucleus:')
    print(df.groupby('sample_simp').total_counts.median()
            .reindex(ORDER).round(0).to_string())
    genes = ['JUN', 'JUNB', 'JUND', 'FOS', 'FOSL2', 'GFAP', 'VIM', 'CSPG4']
    print('\nfraction of nuclei with >= 1 count:')
    det = df.groupby('sample_simp')[genes].apply(lambda d: (d > 0).mean())
    print(det.reindex(ORDER).round(3).to_string())


def p6_check(df):
    """c-Jun vs c-Fos. Reports raw totals; refuses to compute a ratio on zeros."""
    for g in ['FOS', 'FOSB', 'FOSL1', 'FOSL2', 'JUN', 'JUNB', 'JUND']:
        print(f'{g:6s} total counts {int(df[g].sum()):5d}   '
              f'nuclei > 0: {int((df[g] > 0).sum()):4d}')
    if df.FOS.sum() == 0:
        print('\nFOS is zero in every nucleus -> P6 cannot be tested with this '
              'dataset.')


def nb_rate_ratio(sub, gene):
    """Depth-adjusted fold difference in a gene, JUN+ vs JUN- nuclei.
    Negative binomial GLM with log(total UMIs) as offset."""
    y = np.round(sub[gene].values).astype(int)
    X = sm.add_constant((sub.JUN > 0).astype(float).values)
    fit = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=1.0),
                 offset=np.log(sub.total_counts.values)).fit()
    b, se = fit.params[1], fit.bse[1]
    return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se)


def p2_check(df, conditions=('SS', '1_wpi', '2_wpi'),
             genes=('GFAP', 'VIM', 'CSPG4')):
    """P2: per-cell readout difference between JUN+ and JUN- nuclei, within
    each condition, for all nuclei and for the better-sequenced half.
    Dropout makes these LOWER BOUNDS on the true difference."""
    rows = []
    for cond in conditions:
        s = df[df.sample_simp == cond]
        for label, d in (('all', s),
                         ('top-half depth',
                          s[s.total_counts >= s.total_counts.median()])):
            row = {'condition': cond, 'subset': label,
                   'JUN+': int((d.JUN > 0).sum())}
            for g in genes:
                rr, lo, hi = nb_rate_ratio(d, g)
                row[g] = f'{rr:.2f} [{lo:.2f}-{hi:.2f}]'
            rows.append(row)
    out = pd.DataFrame(rows)
    print(out.to_string(index=False))
    print('\nmodel prediction P2: GFAP ~1.22x per cell')
    return out


def dropout_sensitivity(observed_rr, contamination):
    """If a fraction c of 'JUN-' nuclei truly express c-Jun, the true fold
    difference R satisfies observed = R / ((1 - c) + c * R)."""
    c = contamination
    return observed_rr * (1 - c) / (1 - observed_rr * c)


if __name__ == '__main__':
    df = load()
    data_overview(df)
    print('\n' + '=' * 60 + '\nP6\n' + '=' * 60)
    p6_check(df)
    print('\n' + '=' * 60 + '\nP2\n' + '=' * 60)
    p2_check(df)
    print('\nIf the 1 wpi GFAP estimate (1.47x) is diluted by misclassified '
          'nuclei, the true value would be:')
    for c in (0.1, 0.2, 0.3, 0.5):
        print(f'  {int(c * 100)}% of JUN- nuclei truly positive -> '
              f'{dropout_sensitivity(1.47, c):.2f}x')
