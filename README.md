# The glial AP-1 switch is carried by autoregulatory feedback, not by competition for c-Fos

Computational model and analysis code for an independent study of the AP-1 /
miR-200a regulatory switch in axolotl spinal cord glia, and its relationship to
persistent mammalian glial scarring.

Shaunak Patel, North Carolina School of Science and Mathematics.

---

## What this is

A nine-species mass-action model of c-Fos, c-Jun and JunB dimerisation under
miR-200a-dependent repression of c-Jun, with GFAP as a downstream readout. The
model is calibrated against published axolotl injury data and used to ask which
part of the network sustains a bistable glial state, what that implies for
transcriptional readouts, and whether the mammalian–urodele difference can be
explained by AP-1 production ratios.

The headline results, with the sections of the report that carry them:

- Competition for limiting c-Fos generates neither ultrasensitivity nor
  bistability; autoregulatory feedback does (§2.3).
- Because c-Fos is near-saturated in both states, c-Fos-containing dimers are
  poor state sensors, which bounds any readout driven by them alone (§2.4).
- The intervention ranking is largely but not completely robust across
  parameter space, and makes T-5224's efficacy contingent on whether it blocks
  c-Jun homodimer DNA binding (§2.8).
- The model's principal interspecies prediction is **not supported** by murine
  spinal cord injury transcriptomes (§2.9). This is the strongest evidential
  result in the project.

A complete record of predictions, corrections and superseded claims is in
[`docs/PREDICTIONS_v3.md`](docs/PREDICTIONS_v3.md). It is not an appendix —
several conclusions in the report were retracted during development, and that
file is where they are documented.

## Requirements

```
python >= 3.11
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Pinned versions are those used to produce the published figures. `statsmodels`
is required only by `sc_validation.py` (negative binomial regression against
the Lust et al. single-nucleus atlas); every other module runs without it.

## Layout

| Path | Contents |
|---|---|
| `ap1_model_v3.py` | Core nine-species model. Every other module imports it. |
| `round2.py` | Reference parameterisation and principal analyses (Figures 2–6) |
| `schematic.py` | Network topology and per-species state separations (Figure 1) |
| `gsa.py`, `gsa_diag.py` | Global sensitivity analysis and verification (Figure 7) |
| `dynamic_cfos.py` | Injury-induced c-Fos and autorepression (Figure 8) |
| `pharmacology.py` | Intervention representations and engagement scans (Figure 9) |
| `regime.py` | Regime mapping and the interspecies prediction (Figure 10) |
| `mouse_ratio.py` | GSE162610 extraction and comparison (Figure 11) |
| `drug_robustness.py` | Intervention ranking across parameter space (Figure 12) |
| `sc_validation.py` | Single-nucleus comparisons (Lust et al., 2022) |
| `ap1_stochastic.py` | Population-fraction prediction under extrinsic noise |
| `ap1_readouts.py` | Promoter polarity layer and the dynamic-range bound |
| `cfos_ceiling_check.py`, `fj_ablation_v3.py` | c-Fos saturation and routing ablation controls |
| `test_model_properties.py` | Property tests for four model invariants |
| `kd_test.py` | Dissociation-constant sweep and compensation search |
| `replicate_structure.py` | Replicate structure and protocol confounding in GSE162610 |
| `make_supp_tables.py` | Generates Supplementary Tables S1–S3 |
| `gsa_results/`, `drug_robustness_results/`, `regime_results/` | Cached outputs of the long-running scans |
| `figures_round2/` | Published figures |
| `docs/` | Corrections record, seed provenance, supplementary tables |

## Reproducing the results

Long scans write to `results/` and are cached, so figures regenerate in seconds
once the scans have been run. Run everything from ``:

```bash
cd src
```

**Property tests** — verify the four invariants the report argues from. Run this
first; it exits non-zero if any property fails.

```bash
python test_model_properties.py
```

**Figure 1** — network topology and state separations.

```bash
python schematic.py
```

**Figures 2–6** — reference-model analyses.

```bash
python round2.py
```

**Figure 7** — global sensitivity analysis. The scan itself is long; run it in
batches, then draw.

```bash
python run_gsa.py          # populates results/gsa_results/
python -c "import gsa; gsa.figure(gsa.load())"
```

**Figure 8** — dynamic c-Fos.

```bash
python -c "import dynamic_cfos; dynamic_cfos.make_all()"
```

**Figure 9** — intervention engagement.

```bash
python pharmacology.py
```

**Figures 10–12** and the murine comparison require either cached results or the
GSE162610 count matrix; see "Data" below.

**Dissociation-constant analysis** (§4.1):

```bash
python kd_test.py
```

**Supplementary tables:**

```bash
python make_supp_tables.py
```

## Data

No new data were generated. Three published datasets are used:

- **GSE162610** (Milich et al., 2021) — murine spinal cord contusion
  transcriptomes, from the Gene Expression Omnibus. The count matrix is not
  redistributed here; download it from GEO and point `mouse_ratio.py` at the
  directory. `replicate_structure.py` needs only the barcode metadata.
- **Lust et al. (2022)** — axolotl single-nucleus atlas, from the authors'
  deposited data.
- **Sabin et al. (2019)** — published figures and tables.

Note that `FOS` registers zero counts across all 7,343 nuclei in the Lust
atlas while `FOSB`, `FOSL1` and `FOSL2` are detected. This is most plausibly an
annotation gap, and it is why the axolotl arm of the interspecies comparison
could not be evaluated.

## Seeds

Every stochastic analysis and its seed is listed in
[`docs/SEEDS_AND_SAMPLES.md`](docs/SEEDS_AND_SAMPLES.md). In brief: global
sensitivity analysis seed 2026 (400 samples); c-Fos ceiling scan seed 7 (150
draws, 27 bistable); stochastic ensemble seeds 11, 29, 77, 101; ablations seed
42; property tests and dissociation-constant analysis seed 20261015. The regime
map, engagement scans and GSE162610 extraction are deterministic.

## Status and attribution

This work has not been peer reviewed. The author welcomes correction.

Claude (Anthropic) was used for coding.

Yonatan Degefu (Fallahi-Sichani Laboratory, University of Virginia) reviewed the
model and advised on the sequence of analyses presented in Figures 2–6.

Conducted as an independent study under the supervision of Dr. Benton at NCSSM.

## Licence

Code is released under the MIT Licence; see [`LICENSE`](LICENSE).
