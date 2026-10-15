# Random seeds and sampling designs

Every stochastic or sampled result in the manuscript, with the seed and design
needed to regenerate it. Deterministic analyses are listed too, marked as such,
so that the absence of a seed is explicit rather than an omission.

Generated 2026-10-15. Verify against the modules before posting.

## Sampled ensembles

| Analysis | Module | Seed | Design | Feeds |
|---|---|---|---|---|
| Global sensitivity analysis | `gsa.py` | `SEED = 2026` | Latin hypercube, `N_SAMPLES = 400`, 18 parameters, 16-fold range about the reference; log-uniform for rates and production terms, uniform for Hill exponents and the repression leak | Figure 7; the 193 bistable sets; the 177 with a locatable window midpoint |
| Intervention ranking across parameter space | `drug_robustness.py` | none — consumes `gsa.load()` | deterministic given the GSA ensemble; engagement raised by continuation | Figure 12, Table S3 |
| c-Fos ceiling across parameter space | `cfos_ceiling_check.py` | `seed = 7` | `n_sets = 150` LHS draws, `n_ic = 10` initial conditions; 27 returned two stable states | Section 2.4 ensemble medians (1.46-fold, 736-fold, 97.5%, 72.4%) |
| FJ-vs-JJ routing ablation | `fj_ablation_v3.py` | `seed = 42` | `n_sets = 150`, `n_ic = 8` | Section 2.4 ablation |
| v3 global scan (ablation support) | `ap1_v3_downstream.py` | `seed = 42` | `n_param_sets = 250`, `n_initial = 8` | ablation percentages |
| Stochastic population fraction | `ap1_stochastic.py` | `seeds = (11, 29, 77, 101)` | Ornstein–Uhlenbeck extrinsic noise, tau = 10 h, 250 h total in 5 h segments; 120 cells control, 150 cells inhibitor | Section 2.7 (control 0.046 +/- 0.035; inhibitor 0.482 +/- 0.039) |
| Property tests | `test_model_properties.py` | `SEED = 20261015` | 500 random states (T1), 40 random readout parameterisations (T3) | Sections 2.4, 2.8, 5.3 |

Note on the stochastic ensemble: seeds 11 and 29 were the single seeds used in
the superseded report (Section 4.2). The four-seed means are the reported
values. `seed_check()` regenerates all four.

## Deterministic analyses (no seed)

| Analysis | Module | Grid |
|---|---|---|
| Reference continuation and window detection | `round2.py` | continuation in miR-200a from both ends of the axis |
| Regime map | `regime.py` | `xs = linspace(0.25, 3.0, 16)` (c-Jun/c-Fos production ratio), `ys = linspace(0.10, 1.50, 12)` (JunB/c-Fos), `MAP_MIRS = linspace(0.0, 2.5, 26)` |
| Engagement scans | `pharmacology.py` | engagement raised in 201 steps by continuation, each step initialised from the preceding steady state |
| Dynamic c-Fos | `dynamic_cfos.py` | threshold bisection and pulse-timing sweeps |
| GSE162610 extraction | `mouse_ratio.py` | pseudobulk counts summed per cell type and timepoint; no sampling |
| Network schematic and state separations | `schematic.py` | evaluated at the reference window midpoint, miR-200a = 0.60 |

## Reference parameter set

`round2.ref_params()` returns the reference set. It is deterministic and
carries no seed. Table S1 should list its numerical values in full.
