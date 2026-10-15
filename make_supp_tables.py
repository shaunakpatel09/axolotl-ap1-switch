"""
make_supp_tables.py -- builds Supplementary Tables S1-S3 as one workbook.

  S1  complete reference parameter set: value, units, role, provenance
  S2  sampling ranges and distributions for the 18 GSA parameters
  S3  per-parameter-set engagement thresholds for all five interventions

Author: Shaunak Patel.  Built with AI assistance (Claude); reviewed and owned
by the author.
"""

import os
import sys

import numpy as np
import pandas as pd

from round2 import ref_params
import gsa
import drug_robustness as dr

OUT = 'Supplementary_Tables_S1_S3.xlsx'

# name -> (units, role, provenance)
META = {
    'a_F': ('nM/h', 'c-Fos basal production', 'Modelling choice; swept in GSA'),
    'a_J': ('nM/h', 'c-Jun basal production (leaky floor)', 'Selected by parameter search (Table 1)'),
    'a_B': ('nM/h', 'JunB basal production', 'Modelling choice; swept in GSA'),
    'a_G': ('nM/h', 'GFAP basal production', 'Modelling choice (readout layer)'),
    'a_G_indep': ('nM/h', 'AP-1-independent GFAP production', 'Set to zero in the reference model (Section 4.1)'),
    'b_J_by_JJ': ('nM/h', 'c-Jun production driven by the c-Jun:c-Jun homodimer', 'Selected by parameter search (Table 1)'),
    'b_J_by_FJ': ('nM/h', 'c-Jun production driven by c-Fos:c-Jun', 'Modelling choice: one quarter of the homodimer loop (Table 1)'),
    'b_G_by_FJ': ('nM/h', 'GFAP production driven by c-Fos:c-Jun', 'Modelling choice (readout layer); GSA negative control'),
    'Km_auto': ('nM', 'Half-maximal threshold of autoregulatory feedback', 'Selected by parameter search (Table 1)'),
    'n_auto': ('dimensionless', 'Hill coefficient of autoregulatory feedback', 'Selected by parameter search (Table 1)'),
    'Km': ('nM', 'Half-maximal threshold, c-Fos:c-Jun activation of GFAP', 'Modelling choice (readout layer)'),
    'n': ('dimensionless', 'Hill coefficient, c-Fos:c-Jun activation of GFAP', 'Modelling choice (readout layer)'),
    'Km_rep': ('nM', 'Half-maximal threshold, c-Fos:JunB repression of GFAP', 'Modelling choice (readout layer); GSA negative control'),
    'n_rep': ('dimensionless', 'Hill coefficient, c-Fos:JunB repression of GFAP', 'Modelling choice (readout layer)'),
    'Km_mir': ('relative units', 'miR-200a concentration giving half-maximal repression of c-Jun', 'Normalisation constant; defines the control axis (Section 4.3)'),
    'n_mir': ('dimensionless', 'Hill coefficient of miR-200a repression', 'Modelling choice'),
    'mir_leak': ('fraction', 'Residual c-Jun production under maximal miR-200a repression', 'Modelling assumption; motivated by Sabin et al. (2019)'),
    'kon_FJ': ('nM^-1 h^-1', 'Association, c-Fos + c-Jun', 'Adopted convention (Degefu et al.)'),
    'koff_FJ': ('h^-1', 'Dissociation, c-Fos:c-Jun (K_d = 0.1 nM)', 'Adopted convention (Degefu et al.)'),
    'kon_FB': ('nM^-1 h^-1', 'Association, c-Fos + JunB', 'Adopted convention'),
    'kon_JJ': ('nM^-1 h^-1', 'Association, c-Jun + c-Jun', 'Adopted convention'),
    'koff_JJ': ('h^-1', 'Dissociation, c-Jun:c-Jun (K_d = 100 nM)', 'Selected in this work; tested across 5-400 nM (Section 4.1)'),
    'kon_BB': ('nM^-1 h^-1', 'Association, JunB + JunB', 'Adopted convention'),
    'koff_BB': ('h^-1', 'Dissociation, JunB:JunB', 'Tied to koff_JJ'),
    'kon_JB': ('nM^-1 h^-1', 'Association, c-Jun + JunB', 'Adopted convention'),
    'koff_JB': ('h^-1', 'Dissociation, c-Jun:JunB', 'Tied to koff_JJ'),
    'junb_weakness': ('dimensionless', 'Factor by which JunB binds c-Fos more weakly than c-Jun', 'Provisional assumption; author correspondence (Table 1)'),
    'g_F': ('h^-1', 'c-Fos degradation', 'Modelling choice; swept in GSA'),
    'g_J': ('h^-1', 'c-Jun degradation', 'Modelling choice; swept in GSA'),
    'g_B': ('h^-1', 'JunB degradation', 'Modelling choice; swept in GSA'),
    'g_G': ('h^-1', 'GFAP degradation', 'Modelling choice (readout layer)'),
    'mir': ('relative units', 'miR-200a concentration (control parameter)', 'Varied throughout; not a fixed parameter'),
}

ORDER = ['a_F', 'a_J', 'a_B', 'a_G', 'a_G_indep',
         'b_J_by_JJ', 'b_J_by_FJ', 'b_G_by_FJ',
         'Km_auto', 'n_auto', 'Km', 'n', 'Km_rep', 'n_rep',
         'Km_mir', 'n_mir', 'mir_leak',
         'kon_FJ', 'koff_FJ', 'kon_FB', 'kon_JJ', 'koff_JJ',
         'kon_BB', 'koff_BB', 'kon_JB', 'koff_JB', 'junb_weakness',
         'g_F', 'g_J', 'g_B', 'g_G', 'mir']


def table_s1():
    p = ref_params()
    rows = []
    seen = set()
    for k in ORDER:
        if k not in p:
            continue
        seen.add(k)
        u, role, prov = META.get(k, ('', '', ''))
        rows.append({'Parameter': k, 'Reference value': p[k], 'Units': u,
                     'Role': role, 'Provenance': prov})
    for k in sorted(set(p) - seen):
        u, role, prov = META.get(k, ('', '(not used in the reference model)', ''))
        rows.append({'Parameter': k, 'Reference value': p[k], 'Units': u,
                     'Role': role, 'Provenance': prov})
    return pd.DataFrame(rows)


def table_s2():
    p = ref_params()
    rows = []
    for name, kind, lo, hi in gsa.SPEC:
        if kind == 'fold':
            ref = p.get(name, np.nan)
            dist = 'log-uniform'
            low, high = lo * ref, hi * ref
            rng = f'{lo}x to {hi}x the reference'
        elif kind == 'log':
            dist, low, high = 'log-uniform', lo, hi
            rng = f'{lo} to {hi} (absolute)'
        else:
            dist, low, high = 'uniform', lo, hi
            rng = f'{lo} to {hi} (absolute)'
        rows.append({
            'Parameter': name,
            'Distribution': dist,
            'Lower bound': low,
            'Upper bound': high,
            'Range specification': rng,
            'Negative control': 'yes' if name in gsa.NEGATIVE_CONTROLS else '',
            'Units': META.get(name, ('', '', ''))[0],
        })
    note = pd.DataFrame([{
        'Parameter': 'r_FJ',
        'Distribution': 'uniform',
        'Lower bound': 0.0, 'Upper bound': 1.0,
        'Range specification': 'ratio of c-Fos:c-Jun feedback to homodimer feedback; b_J_by_FJ = r_FJ * b_J_by_JJ',
        'Negative control': '', 'Units': 'dimensionless'}])
    df = pd.DataFrame(rows)
    df = df[df.Parameter != 'r_FJ']
    return pd.concat([df, note], ignore_index=True)


def table_s3():
    """One row per sampled parameter set; engagement thresholds under both
    criteria, plus durability after withdrawal, for all five interventions."""
    res = dr.load()
    out = {'parameter_set': res['idx'], 'window_midpoint_miR200a': res['m0']}
    for d in res['jump']:
        out[f'{d} | E at bifurcation'] = res['jump'][d]
        out[f'{d} | E at 50% c-Jun'] = res['half'][d]
        out[f'{d} | low state persists after withdrawal'] = res['stays'][d]
    df = pd.DataFrame(out)
    return df.sort_values('parameter_set').reset_index(drop=True)


def main():
    s1 = table_s1()
    s2 = table_s2()
    print(f'S1: {len(s1)} parameters')
    print(f'S2: {len(s2)} sampled parameters')
    s3 = table_s3()
    with pd.ExcelWriter(OUT, engine='openpyxl') as xw:
        s1.to_excel(xw, sheet_name='S1 reference parameters', index=False)
        s2.to_excel(xw, sheet_name='S2 GSA sampling ranges', index=False)
        s3.to_excel(xw, sheet_name='S3 engagement thresholds', index=False)
    print(f'S3: {s3.shape[0]} rows x {s3.shape[1]} columns')
    print('wrote', OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
