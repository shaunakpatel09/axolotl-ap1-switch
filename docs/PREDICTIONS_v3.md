# Falsifiable Predictions of the Calibrated AP-1 Competition Model (v3)

Status key: each item states what the model *predicts*, what was *fitted* (and
therefore cannot count as evidence), and how the prediction could be tested.
Fitted quantities: state separation band (targeted to the Fig 5c-corrected
per-cell values); stochastic noise amplitude sigma (fitted to the control 8%
cJun+ fraction only); the miR-axis placement of control/inhibitor conditions
(chosen a priori from window geometry, not tuned to outcomes).

**Scorecard against data the model was not fitted to.** Three comparisons
exist, and none of them is a clean win:

| | result | verdict |
|---|---|---|
| P2/P4 per-cell GFAP, Lust 2022 | 1.47x measured vs 1.22x model | direction right, dropout-limited |
| P3 cJun+ fraction, Sabin Fig 5c | 0.48 ± 0.04 predicted vs 0.30 | ~1.6x over, sigma poorly constrained |
| P10 production ratio, GSE162610 | 0.27-0.64 measured vs >1.35 predicted | **refuted** |

The project's strongest evidential claim is therefore a refutation of one of
its own predictions, not a confirmation. That is stated here deliberately.

---

## P1 — The durability law for transient miR-200a therapy
A transient miR-200a elevation produces a *permanent* conversion of glia to
the low-cJun (pro-regenerative) state **only if the post-treatment baseline
lies inside the bistable window**; there, even small pulses flip cells
durably (v3: amplitude >= 0.5 relative units). Outside the window, cells
revert regardless of pulse size.
**Test:** compare sustained vs transient miR-200a mimic delivery (cf. the
pluronic-gel delivery of Diaz Quiroz & Echeverri 2014) with a washout arm;
score cJun+ fraction weeks after washout.
*(Survived the v2 -> v3 architecture revision unchanged.)*

## P2 — GFAP barely differs between flipped and unflipped cells
Under the calibrated switch, per-cell GFAP in cJun+ vs cJun- glia should
differ by only ~1.2x — far less than vimentin or CSPG4. The dramatic
"GFAP-high scar state" is a population/tissue-level property, not a
single-cell one.
**Test:** quantify GFAP immunofluorescence intensity in cJun+ vs cJun- GFAP+
cells in Sabin-style sections (data of exactly this kind already exist in
Fig 5b-type images).
**Status (single-cell check, Lust et al. 2022 ependymoglia):** at 1 week
post-injury, JUN+ nuclei carry 1.47x more GFAP than JUN- nuclei [95% CI
1.16-1.86], depth-adjusted. Direction matches. Dropout biases this toward 1,
so it is a lower bound; plausible true range ~1.5-2.8x. The model's 1.22x is
probably too small -- GFAP may be over-buffered. No JUN-associated difference
in uninjured tissue. Caveats: 1 week is after Sabin's 1-4 day switching window;
pallium not spinal cord; mRNA not protein. See sc_validation.py.

## P3 — The responding variable is the fraction of cells, with hysteresis
Graded miR-200a inhibition moves the *percentage* of cJun+ glia along a
steep, history-dependent (hysteretic) dose-response; single-cell intensities
shift far less than the population fraction. Qualitatively matches Fig 5c
(8% -> 30%).
**Test:** titrate inhibitor dose; score cJun+ fraction; look for
history-dependence with sequential dosing.
**Status (stochastic ensemble re-run under the round-2 reference,
ap1_stochastic.py).** Extrinsic OU noise on c-Jun production; sigma fitted to
the CONTROL fraction only, both miR-200a placements fixed a priori from window
geometry (control = upper edge + 0.15; inhibitor = window midpoint). Over four
seeds:

| | model | target |
|---|---|---|
| control, 120 cells *(fitted)* | 0.046 ± 0.035 | 0.08 |
| inhibitor, 150 cells *(predicted)* | **0.482 ± 0.039** | **0.30** |

Direction and scale correct, magnitude ~1.6x over. Essentially unchanged by
the v3 -> round-2 architecture revision (which gave 0.52), so the result does
not depend on the FJ edge.
**Two caveats found by running more than one seed, both of which weaken the
original single-seed report:** (i) the between-seed sd on the control fraction
(0.035) is comparable to the target itself (0.08), so a single 120-cell run
cannot constrain sigma — the control fraction runs 0.000 / 0.008 / 0.017 /
0.083 at sigma 0.28 / 0.30 / 0.32 / 0.35, pinning sigma only to roughly
0.33-0.40; (ii) the previously reported single-seed numbers used the *highest*
of four control seeds for the fit and the *lowest* of four inhibitor seeds for
the prediction. Reported here as the four-seed mean instead.

## P4 — Scar promoters must read Jun-homodimer input (a fork)
Under the calibrated switch, cFos:cJun activity changes only ~1.23x per cell
between states — so **no promoter reading cFos:cJun alone can reproduce the
measured per-cell spreads** (vimentin ~4.1x, CSPG4 ~4.7x). Either:
(a) these promoters additionally integrate Jun-homodimer (or Fos-independent
Jun) input — supported by PROMO finding c-Jun motifs in GFAP/vimentin/CSPG
promoters (Sabin 2019); or
(b) the per-cell back-calculation (applying cJun+ fractions to all genes) is
wrong, which is itself testable by P2/P3-style per-cell measurements.
**Test:** promoter-reporter assays comparing response to cFos:cJun vs cJun
homodimer; or mutate the cJun-only sites in the vimentin/CSPG4 promoters.
**Status (single-cell check):** directly measured per-nucleus differences at
1 week are vimentin 1.38x [1.15-1.65] and CSPG4 0.94x [0.65-1.36] -- CSPG4
does not track JUN at all. Even under heavy dropout, a true ~4.7x difference
would appear larger than this. Leans toward branch (b) for CSPG4, but the
tissue and timepoint differ from Sabin, so it cannot resolve the fork.
**Status (readout-polarity layer, ap1_readouts.py, re-run under the round-2
reference).** The polarity hypothesis — that some scar genes read the mirror
promoter logic — was built to explain CSPG5's opposite direction. It does not
rescue the readouts, because the ceiling is dynamic range, not direction. The
largest per-cell ON/OFF ratio ANY cFos:cJun / cFos:JunB readout can produce,
at phi = 1, is **1.246x** (polarity +1) or 0.100x (polarity -1), against
per-cell targets of GFAP 2.14x, vimentin 4.13x, CSPG4 4.71x, Col4a1 1.18x and
CSPG5 -0.04x. Only Col4a1 is reachable; CSPG5's target is negative and
therefore unreachable at any polarity, consistent with dropping it from the
panel. This is direct support for branch (a): the cause is the c-Fos ceiling —
FJ is pinned near the limiting c-Fos pool in both states and has no dynamic
range to pass downstream — so scar promoters must integrate homodimer or
otherwise Fos-independent c-Jun input. The readout layer is verified strictly
downstream (wrecking the GFAP parameters moves the switch species by 1.8e-15).

## P5 — A bound on Jun-Jun dimer affinity
Shallow bistability (states 2-10x apart, as the data require) exists only if
Jun-Jun dimerisation is substantially weaker than the standard
"100x-weaker-than-Fos-Jun" convention: robust for Kd ~30-150 nM
(koff/koff_FosJun ~ 300-1500 in model units), absent at the conventional
Kd = 10 nM and again at Kd >= ~250 nM.
**Test:** against direct kinetic measurements (Kohler & Schepartz 2001 and
successors); falsifiable by biophysics alone, no animal needed.

## P6 — A stoichiometric requirement: cJun production must exceed cFos production
Bistability exists only when basal cJun production sits in a bounded band
relative to cFos production, so that both states lie above the cFos
sequestration scale. Expressed as a ratio, the model requires
**a_J / a_F ~ 1.2 - 3** (measured across cFos pool sizes: a_F=8 -> ratio
2.0-3.0; a_F=16 -> 1.4-2.1; a_F=32 -> ~1.2, window collapsing; a_F=64 ->
no bistability at any cJun level).
**Test:** quantitative measurement of total cJun vs total cFos protein in
axolotl glia after injury -- the kind of single-cell multiplexed measurement
the Fallahi-Sichani lab performs routinely. If cJun does not exceed cFos in
roughly this band, the bistable interpretation fails.
**Status: untestable with Lust et al. 2022.** FOS has exactly zero counts in
all 7,343 nuclei, while FOSB, FOSL1 and FOSL2 are detected -- more likely a
gene-model/annotation issue than biology. Not counted as support: a ratio
with zero in the denominator is not evidence.
**Status in mouse (GSE162610): measurable, and it does NOT hold.** See the
regime test below. Mouse spinal cord ependymal cells and astrocytes carry
Jun/Fos count ratios of 0.27-0.64 at every timepoint, with c-Fos transcript
running ~2.5x c-Jun. If the mRNA ratio tracks the production ratio, c-Fos is
in EXCESS rather than limiting in these cells, which is the opposite of the
regime the model's reference parameters assume. mRNA is a weak proxy here
(c-Fos protein is short-lived and heavily post-transcriptionally controlled),
so this is evidence against the premise rather than a refutation of it -- but
it is evidence against, and it should be stated as such.

---

## Corrections and scope changes since the first v3 run

- **FJ-loop mechanism corrected.** Adding cFos:cJun-driven autoregulation on
  top of the homodimer loop drops multistability from 21.3% to 1.3% (150
  matched parameter sets; FJ-only and neither both 0.0%). The earlier
  explanation -- that FJ feedback consumes cFos -- is WRONG. At the reference
  point the FJ Hill term reads 0.891 (ON) vs 0.846 (OFF): it is saturated, and
  replacing it with a constant of equal size abolishes the window identically.
  Correct mechanism: cFos is limiting, so FJ is pinned near the cFos ceiling in
  both states (10.06 vs 8.83, 1.14x) and cannot sense the cell's state; the
  homodimer escapes that ceiling and carries the state signal.
  **Generalised across 27 bistable LHS parameter sets** (cfos_ceiling_check.py):
  FJ differs between states by only 1.46x (IQR 1.27-1.77) while JJ differs by
  736x (IQR 80-2585) -- a ~535x difference in state information, with JJ > 3x
  in 100% of sets. Refinement: c-Fos is essentially fully bound as FJ in the ON
  state (median 97.5%) but 72% in the OFF state, so FJ's dynamic range is
  compressed rather than identically pinned in both. With the Hill exponent
  fixed at 3, even an optimally placed threshold could amplify FJ's 1.46x to
  only ~3x. Measured at attractors at mir = 0, so separations run larger than
  the 3.48x at the calibrated reference. See fj_ablation_v3.py.
- **v2 phi calibration void.** Under v3, cFos:cJun activity changes only
  1.23x between states, so the per-cell targets are unreachable by any
  FJ-only readout. Superseded by P4.
- **CSPG5 dropped from the glial scar panel.** It decreases while the
  cJun+ fraction increases, which yields a negative per-cell target -- it
  cannot be a readout of this switch. It clusters with neuronal genes in
  Sabin Fig 6a.
- **Developmental arm unsupported.** Under v3 the v2 cliff reverses and
  depends on an invented development-to-parameter mapping. Rebuild from the
  literature (glial vs neuronal locus of the age effect; McKeon 1995, Laabs
  2005) or cut.
- **Model scope: JunD and Fra-2.** In the Lust atlas JUND is expressed about
  as much as JUN, and FOSL2 (Fra-2) is the most-detected FOS-family member.
  Neither is in the model. JunD would be a third Jun competitor; Fra-2 is a
  candidate FOS-family partner at later timepoints.

---

## Round 2 revisions (after Degefu's feedback)

- **Reference model now includes the c-Fos:c-Jun -> c-Jun edge.** c-Jun
  autoregulation through Jun/AP-1 is established (Angel et al. 1988, Cell),
  so removing it in v3 was not defensible. Reference: b_J_by_FJ = 9 (a quarter
  of the homodimer loop) with basal c-Jun lowered 16 -> 8. Both are choices,
  not measurements. Behaviour is nearly unchanged: window 0.35-0.80, c-Jun
  separation 3.6x, GFAP 1.24x. See round2.py, Figure 4.
- **The FJ edge acts as extra basal drive** (second refinement of this
  finding). Bistability follows a band of constant total drive, basal +
  ~0.87 x FJ edge ~ 16. It survives an FJ edge up to about half the homodimer
  loop and collapses by three-quarters (b_J_by_FJ >= 27). So the earlier
  statement that FJ feedback "suppresses" bistability was a drive-balance
  effect at fixed basal c-Jun, not an intrinsic property of the edge. The
  switching feedback must still come from the homodimer, because FJ sits near
  the c-Fos ceiling and cannot sense the cell's state.
- **Competition alone is graded, not steep.** Without autoregulation the
  response is single-valued with effective Hill coefficient 2.0-2.7 depending
  on the species read out, against 2.0 for the miR-200a input. This revises
  hypothesis H1, which expected competition alone to be sharp. Both the
  steepness and the two states come from the feedback. (Reference regime
  only; the global sensitivity analysis should test other regimes.)
- **P6 revised.** Two stable states exist only at intermediate c-Fos
  abundance: no window at c-Fos production <= 2 or >= 12 nM/h, a window at
  3-10 (reference parameters with the FJ edge). The earlier claim that c-Jun
  must exceed c-Fos by 1.2-3x was too narrow; c-Fos abundance is bounded on
  both sides.
- **JunB autoregulation.** The literature describes c-Jun as positively
  autoregulated and JunB as an NF-kB target rather than self-induced (Mathas
  et al. 2002, EMBO J), consistent with the model's asymmetry. Shown in
  Hodgkin lymphoma cells, not glia.
- **Parameter provenance, stated plainly.** Two conventions come from Degefu
  et al. (Fos-Jun binding ~100x tighter than Jun-Jun; dimer degradation as the
  mean of its monomers). JunB weakness (3x) is a placeholder from Degefu's
  email. The bistable regime (basal c-Jun, homodimer feedback 36, Km 5,
  Hill 3, Jun-Jun off-rate 1000) was found by grid search. The parameters are
  NOT literature-derived values.
- **Global sensitivity analysis (gsa.py, Figure 6).** LHS-PRCC (Marino et
  al. 2008), 400 samples, 18 parameters each over a 16-fold range around the
  reference. Bistable somewhere on the miR-200a axis in 193 of 400 samples;
  8 random bistable samples verified as two Jacobian-stable fixed points after
  20,000 h. Results by conclusion:
  - Controls of bistability (PRCC on window width, Bonferroni-significant):
    homodimer feedback strength +0.50, c-Jun degradation -0.30, basal c-Jun
    -0.28, feedback cooperativity +0.27, feedback threshold Km -0.19. Binding
    constants (c-Fos:c-Jun and Jun-Jun off-rates, JunB weakness) and c-Fos /
    JunB production were not significant. Both negative controls (GFAP-readout
    parameters) non-significant, as required.
  - C2 HOLDS: competition alone never sharpens the response meaningfully.
    Steepness relative to the miR-200a input: median 1.00, max 1.35 for the
    c-Fos allocation; median 1.03, max 1.59 for free c-Jun.
  - C1 REFINED: bistability without the homodimer loop occurred in 42 of 400
    samples (all 42 verified). In those, c-Jun is far below c-Fos (c-Jun/c-Fos
    median 0.12 in the low state), FJ differs 4.4x between states, and they
    concentrate at low c-Jun/c-Fos production ratios (30% bistable below 0.29,
    0-1% above 1). So the homodimer is required when c-Jun exceeds c-Fos; when
    c-Fos is in excess, FJ is no longer capped and the c-Fos:c-Jun edge can
    carry the switch itself. This is the c-Fos ceiling mechanism's prediction,
    and it makes the c-Jun:c-Fos abundance ratio in axolotl glia (P6) the
    measurement that decides which loop carries the switch.
  - C3 MOSTLY HOLDS: c-Jun separates more than GFAP in 164 of 193 bistable
    samples (median 4.15x vs 1.52x). The global GFAP median (1.52x) is closer
    to the single-cell estimate (1.47x, lower bound) than the reference value
    (1.24x), so the reference GFAP readout sits at the low end.
  - C4 ONLY PARTLY HOLDS globally: the lowest c-Fos quintile is less often
    bistable (35% vs ~50%), but abundant c-Fos does not prevent bistability
    once other parameters vary. The window loss at high c-Fos in Figure 2 is
    specific to the reference parameters.
  - C5 WEAKENED globally: a strong FJ edge lowers the bistable fraction only
    modestly (~52% -> ~43%); other parameters compensate.
  - Not tested: the A15 deviation itself. The Jun-Jun off-rate range
    (250-4000) sits around the reference, so it does not reach the
    conventional value of 100.
- **c-Fos as a regulated, injury-induced input** (Degefu Q5; dynamic_cfos.py,
  Figure 7). Adds c-Fos autorepression driven by dimerised c-Fos (Sassone-Corsi,
  Sisson & Verma 1988, Nature) and an injury pulse peaking at 24 h, with JunB
  optionally co-induced as in Sabin. All other equations unchanged; the wrapper
  reproduces v3 exactly with the additions off.
  - Fixing c-Fos was a safe simplification at the reference: with autorepression
    at matched baseline the window is identical (0.35-0.80, 3.60x).
  - Autorepression buffers the switch. With c-Fos fixed, raising production 4x
    removes it (upper edge 0.80 -> none); with autorepression (0.8), a 16-fold
    range moves the upper edge only 0.85 -> 0.60.
  - A c-Fos pulse alone, with no miR-200a, switches a cell sitting in the window
    to the low-c-Jun state: the surge pairs up c-Jun and starves the homodimer
    loop. Smallest pulse (peak extra production): 0.73 with no FJ edge, 0.91 at
    the reference, never with the FJ edge equal to the homodimer loop. So how
    strongly c-Fos:c-Jun induces c-Jun decides whether injury-induced c-Fos
    promotes or blocks the switch. (FJ edge at half the loop omitted: no
    baseline with two distinct states could be located.)
  - JunB co-induction barely changes these thresholds (0.74, 0.94).
  - Timing: with c-Fos at half its threshold, the miR-200a pulse needed falls
    from 0.25 to 0.12 when c-Fos peaks with or after miR-200a (48-120 h); a peak
    before 24 h does almost nothing. At Sabin's c-Fos peak (24 h): 0.25 -> 0.20.
  - Caveats: the pulses are broad (days) where real c-Fos induction is sharper;
    lasting switches still require the baseline to sit inside the window.

---

## Round 3: the drug layer, and the species question tested against data

### P7 — Durability is a property of the regime, not of the drug
*(pharmacology.py, Figure 8; regime.py, Figure 9)*

Five interventions are represented by what they do to existing parameters;
no equation changes. Engagement E is the fraction of the target's activity
blocked. Equal engagement is NOT equal dose — converting engagement to dose
needs potency and pharmacokinetics, which this model deliberately does not
attempt (see "Why there is no pharmacokinetics arm" below).

  * **miR-200a mimic** — adds miR-200a; engagement reported as the fraction
    of c-Jun production silenced, capped by the repression leak (0.15).
  * **JNK inhibitor (activity)** — scales c-Jun self-activation (homodimer
    and FJ edges) and FJ activation of GFAP by (1-E). JNK phosphorylation of
    c-Jun Ser63/73 strengthens transactivation (Hibi 1993; Angel & Karin).
  * **JNK inhibitor (activity + stability)** — as above plus c-Jun
    degradation x (1+E); Ser73 phosphorylation protects c-Jun from
    ubiquitination (Fuchs 1996; Musti 1997), but JNK also helps target it, so
    the two JNK variants bracket an unknown. The (1+E) factor is a
    placeholder, not a measurement.
  * **T-5224 (c-Fos dimers)** and **(c-Fos dimers + homodimer)** — the
    literature does not settle whether T-5224 touches the c-Jun homodimer, so
    both mappings are carried throughout.

At the reference set, the sustained engagement needed to remove the
high-c-Jun state falls as baseline miR-200a rises (mimic 0.34 -> 0.08 across
miR-200a 0 -> 0.70), and **the switch persists after withdrawal only when the
baseline lies inside the two-state window** — at miR-200a 0.45 and above, not
below. The drug does not hold the cell down; it walks it across a boundary.
This is P1's durability law restated for sustained interventions, and it
means the regime a cell sits in decides whether a drug is curative or merely
suppressive. In a simulated locked regime every drug costs roughly twice the
engagement and nothing persists (Figure 9C-D).

### P8 — T-5224's effect depends entirely on whether it blocks the c-Jun homodimer
*(pharmacology.py; drug_robustness.py, Figure 11)*

If T-5224 spares the c-Jun homodimer, blocking **every** c-Fos-containing
dimer leaves the high-c-Jun state intact below miR-200a ~0.5 at the reference
set. Across the 177 bistable GSA parameter sets at which this could be
evaluated, T-5224 sparing the homodimer **fails to clear the state in 39% of
sets**, and where it does work it costs **2.4x** the engagement of the
mapping that also hits the homodimer (median).
**Test:** a DNA-binding or reporter assay asking whether T-5224 blocks c-Jun
homodimer binding, in any mammalian cell line. T-5224 has been through human
phase II trials in rheumatoid arthritis, so this is a question about an
existing compound, not a hypothetical one. This is the sharpest testable claim
the model produces.

### P9 — The target ranking is a property of the wiring, not of the parameters
*(drug_robustness.py, Figure 11)*

Each of the 193 bistable GSA samples was evaluated at the midpoint of its own
two-state window (177 had a locatable midpoint). Median engagement, and the
fraction of sets in which each intervention clears the state at all:

| intervention | clears | median E | IQR | cheapest in |
|---|---|---|---|---|
| JNK (activity + stability) | 100% | 0.11 | 0.05-0.18 | 150/177 |
| miR-200a mimic | 83% | 0.12 | 0.04-0.22 | 24/177 |
| JNK (activity) | 100% | 0.17 | 0.09-0.28 | 3/177 |
| T-5224 (+ homodimer) | 100% | 0.17 | 0.09-0.28 | 0 |
| T-5224 (c-Fos only) | 61% | 0.27 | 0.11-0.48 | 0 |

Pairwise orderings from the reference set survive the parameter space:
mimic < JNK-activity in 98% of sets, JNK-with-stability < JNK-activity in 94%,
< mimic in 84%, < T-5224-c-Fos-only in 96%; mimic < T-5224-c-Fos-only in 100%.
The weakest is JNK-activity < T-5224-c-Fos-only at 74% — clear of chance, but
not assertable without the qualifier. The same conclusions hold under both
definitions of "removes the state" (the bifurcation, and c-Jun halved).
**Internal check:** JNK-activity and T-5224-with-homodimer scale the same
three c-Jun parameters and differ only in Km_rep, which enters the GFAP
equation alone, so they MUST agree on c-Jun in every sample. Maximum
difference across 177 sets: 0.00e+00.

**Why this is the route to human relevance.** Every structural element of the
model comes from mammalian literature — c-Jun autoregulation (Angel 1988,
human/mouse), AP-1 dimer competition (Comandante-Lou/Degefu, mammalian lines),
the GFAP promoter polarity (Sabin), the JNK phosphorylation biology. The
axolotl supplied the motivating question, not a single parameter. A ranking
that holds across the parameter space is therefore a claim about any cell with
this architecture, mammalian included.

### P10 — REFUTED: the production ratio does not explain the species difference
*(regime.py, Figure 9; mouse_ratio.py, Figure 10)*

**Prediction made before looking:** axolotl glia sit inside the two-state
window and mammalian reactive astrocytes sit outside it, locked in the
high-c-Jun state, because a higher c-Jun/c-Fos production ratio slides the
window to miR-200a levels the cell does not reach. The model places that
boundary at a_J/a_F = 1.35 (at baseline miR-200a 0.58).

**Tested against GSE162610** (Milich et al. 2021, mouse spinal cord contusion,
uninjured + 1/3/7 dpi), as a pseudobulk ratio sum(Jun)/sum(Fos) within each
cell type — counts summed before dividing, so per-cell library size cancels
and no normalisation is applied:

| | Uninjured | 1 dpi | 3 dpi | 7 dpi |
|---|---|---|---|---|
| Astrocyte (n = 123/492/197/118) | 0.64 | 1.41 | 0.35 | 0.50 |
| Ependymal (n = 2678/2869/1220/950) | 0.388 | 0.335 | 0.267 | 0.400 |

Approximate 95% intervals (Poisson on the counts, a floor on the real
uncertainty): astrocyte 0.57-0.72, 1.31-1.52, 0.29-0.43, 0.44-0.57;
ependymal 0.385-0.391, 0.331-0.339, 0.261-0.274, 0.390-0.410.

**The prediction fails.** Seven of eight points lie inside the switchable
band. The only crossing is astrocytes at 1 dpi (1.41), which is the acute
immediate-early response — Jun, Fos and Junb all rise together — and it is
gone by 3 dpi, when the ratio sits *below* uninjured. Ependymal cells, the
homolog of the axolotl ependymoglia the model is about, do not move at all
across the time course, and their intervals are tight enough that this is not
a power problem. At 7 dpi, with the scar established, the highest value across
all groups is 0.50 — 2.7x short of the boundary.

**Technical control (dissociation).** Fos, Jun and Junb are immediate-early
genes induced by enzymatic dissociation itself (van den Brink et al. 2017).
The dataset carries two dissociation protocols, and in ependymal cells the
protocol moves the ratio as much as injury does, with inconsistent sign:
Enriched minus Standard is +0.138 (uninjured), +0.066 (1 dpi), **-0.033**
(3 dpi), +0.050 (7 dpi), against a within-protocol biological range of 0.161
(Enriched) and 0.100 (Standard). Chemistry adds more (v2 0.283, v3 0.395), and
a_B/a_F swings 60% on protocol alone. So the fine movements in the table are
not trustworthy — but the refutation does not rest on them, because protocol
effects of ±0.14 cannot close a gap of 0.9.

**What this test did NOT interrogate.** The map has two coordinates and only
one was measured. Baseline miR-200a is expressed in multiples of its
repression constant — a property of the miR-200a-to-c-Jun interaction that no
expression dataset contains, and that miRNA-seq could not supply either
without a dose-response calibration. A miR-200a-based explanation of the
species difference remains untested, and is now a conjecture rather than a
supported claim. Figure 10A labels the assumed coordinate explicitly.

### Structural finding: competition for c-Fos does not move the switch
*(regime.py, Figure 9B)*

Across a fifteen-fold change in JunB production (a_B/a_F 0.1 -> 1.5) the
two-state window changes from 0.50 to 0.45 and total c-Jun by 1.8%, while
GFAP falls 17% (54.8 -> 45.5 nM). JunB sets the OUTPUT, not the SWITCH.
Sweeping JunB's affinity for c-Fos from 1x to 6x changes nothing (window
0.45-0.50 throughout). But strengthening the JunB:c-Jun heterodimer tenfold
(koff_JB 1000 -> 100) widens the window to 0.75.

So in this model JunB influences bistability only by sequestering c-Jun
directly, never by competing for c-Fos — because the engine is the c-Jun
homodimer, which contains no c-Fos, so removing c-Fos from c-Jun leaves the
engine untouched. **This qualifies the project's own framing**, which is built
on competition for limiting c-Fos, and it is aimed squarely at the framework
the collaborator's own model is built on. Caveat: koff_JB = 1000 is the
audit-flagged A15 choice, and this conclusion is sensitive to exactly that
parameter.

### Why there is no pharmacokinetics arm
A PK layer — half-lives, dosing intervals, clearance — was considered and
rejected twice, the second time on a harder ground than the first. PK
describes what a body does to a drug; it does not convert a model of one
species' network into another's. And after P10 the obstacle is structural: a
cell cannot be placed on the model's miR-200a axis from any expression
dataset, so a dosing schedule would be computed for a cell whose regime is
unknown. The ranking's robustness across parameter space (P9) is the honest
substitute, and the T-5224 question (P8) is the honest deliverable.

---

## CORRECTIONS: window-detection bug (found during the dynamic c-Fos work)

The detector behind round-2 Figures 2 and 4, the earlier c-Fos abundance test,
and the A14 sweep started its rising miR-200a sweep from zero c-Jun. When the
two-state window extends down to miR-200a = 0, that sweep lands on the low
branch and reports no window. round2.py now uses a detector that starts from
high c-Jun and verifies both states. These earlier statements are superseded:

- "Two stable states exist only at intermediate c-Fos abundance: no window at
  c-Fos production <= 2 or >= 12" (round-2 P6 revision): WRONG. Windows exist
  from 1 to 24. c-Fos abundance shifts where the window sits: its upper edge is
  0.85 at low c-Fos, 0.80 at the reference, 0.65 at 12, 0.50 at 16, 0.10 at 24,
  with no high-c-Jun state at 32. More c-Fos means less miR-200a is needed to
  break the high-c-Jun state.
- "The window loss at high c-Fos in Figure 2 is specific to the reference
  parameters" (GSA C4 note): it was a detection artifact. The GSA itself used a
  correct detector, so its results stand.
- "Bistability follows a band of constant total drive ... survives an FJ edge up
  to about half the homodimer loop and collapses by three-quarters": WRONG.
  Bistability fills the region below an upper limit of total c-Jun drive. The
  FJ edge adds to that drive, so a stronger edge lowers the basal c-Jun the
  switch tolerates; with the edge equal to the homodimer loop, bistability
  survives at low basal c-Jun as a deep switch. At strong FJ edges the
  two-state region can be fragmented.
- A14 "bistability requires basal c-Jun ~16-24, gone at <= 13": WRONG. With v3
  parameters bistability spans basal c-Jun 8 to 28 (windows reach miR-200a = 0
  for 8-14), with c-Jun separations of 2.1-6.5x. A14 is far less fragile than
  reported.
- Unaffected: the GSA, all ablations and the c-Fos ceiling scan (these start
  from many initial conditions), the reference window, and Figures 1, 3 and 5.

---

## Resolution of A14 (previously listed as the model's fragile assumption)
A14 was: "bistability requires basal cJun in a narrow band (a_J ~ 16-24)."
Swept across that entire band, every conclusion is INVARIANT:

| a_J | window gap | cJun ON/OFF | GFAP ON/OFF | ablation removes window |
|-----|-----------|-------------|-------------|--------------------------|
| 15  | 0.55      | 3.68x       | 1.22x       | yes |
| 16  | 0.45      | 3.48x       | 1.22x       | yes |
| 18  | 0.35      | 3.13x       | 1.20x       | yes |
| 20  | 0.20      | 2.85x       | 1.18x       | yes |
| 22  | 0.15      | 2.61x       | 1.15x       | yes |
| 24  | 0.05      | 2.37x       | 1.12x       | yes |

Separation stays inside the 2-10x data band at every value; the per-cell GFAP
claim (P2) holds throughout; the autoregulation-is-the-engine result holds at
every value. A14 therefore governs WHETHER the switch exists, not WHAT is
concluded from it. Restated as P6 above, it is a testable stoichiometric
prediction rather than a free parameter. Note also that a ~1.5-fold parameter
window is unremarkable for bistability -- the reference framework found only
~5% of its sampled parameter space multistable.

---

## Sensitive assumptions the predictions depend on
- A15: the Jun-Jun weakening of P5 (a deviation from the Comandante-Lou
  convention; P5 converts the assumption into the falsifiable claim).
- A16: the GFAP equation. Sabin Fig 2b establishes the DIRECTIONS (cFos:cJun
  activates, cFos:JunB represses); the magnitudes, the multiplicative form of
  repression, and the Hill exponents are unmeasured choices. The basal
  AP-1-dependent term is inside the repression, the AP-1-independent floor is
  not. The STAT3/NF-kB floor (a_G_indep) is 0 in the reference parameters,
  so default runs overstate GFAP's AP-1 dependence.
- A5 (restated): every result depends on cJun homodimer autoregulation in
  axolotl glia, which has not been shown.
- The Fig 5c fractions (8%/30%) were read from a published scatter plot and
  applied to bulk qPCR ratios; re-derive both from the figure before any
  submission.
- A17: the drug mappings of P7-P9. Each intervention is a scaling of existing
  parameters chosen to match a described mechanism, not a fitted
  pharmacodynamic model. The JNK stability factor (1+E) is the weakest of
  them — a bare placeholder for an effect whose sign is contested, since JNK
  both protects c-Jun via Ser73 and helps target it for degradation. The two
  JNK variants are reported precisely so the unknown is visible rather than
  hidden in one curve.
- A18: the identifiability limit exposed by P10. Baseline miR-200a is in
  multiples of an unmeasured repression constant, so no expression dataset can
  place a real cell on the model's miR-200a axis. Any claim about which regime
  a real cell occupies — and therefore any claim about whether a drug would be
  curative or suppressive in it — is beyond what this model can currently
  support without a dose-response calibration.

---

## Files

| file | what it holds |
|---|---|
| `ap1_model_v3.py` | the core model; everything imports it |
| `schematic.py` | network topology and per-species state separations, Figure 1 |
| `round2.py` | `ref_params()` and Figures 2-6 |
| `gsa.py`, `gsa_diag.py` | global sensitivity analysis, Figure 7 |
| `dynamic_cfos.py` | injury-induced c-Fos, Figure 8 |
| `pharmacology.py` | drug layer, Test 1, Figure 9 |
| `regime.py` | regime map and the species arm, Figure 10 |
| `mouse_ratio.py` | GSE162610 extraction and Figure 11 |
| `drug_robustness.py` | the ranking across parameter space, Figure 12 |
| `sc_validation.py` | Lust 2022 single-cell checks (P2, P4, P6) |
| `ap1_stochastic.py` | the Fig 5c fraction prediction (P3) |
| `ap1_readouts.py` | readout polarity layer and the 1.25x ceiling (P4) |
| `ap1_v3_downstream.py` | dependency of the two below; v2 experiments ported to v3 |
| `cfos_ceiling_check.py` | the c-Fos ceiling across bistable parameter sets |
| `fj_ablation_v3.py` | the FJ-vs-JJ routing ablation under v3 |

Figure numbers above are those of the manuscript as of 2026-10-12, after
`schematic.py` was inserted as Figure 1 and everything downstream shifted by
one. The deck prepared for Dr. Benton predates that shift and still uses the
old numbering (its "Figure 8" is the manuscript's Figure 9, its "Figure 10"
is Figure 11, its "Figure 11" is Figure 12).

## CORRECTION (2026-10-12): the c-Fos ceiling numbers were misattributed in the manuscript

This log had the provenance right; the manuscript did not. Section 2.4 of the
paper quoted the ensemble medians -- FJ 1.46x, JJ 736x, and "97.5% of total
c-Fos bound as FJ in the high state, 72% in the low state" -- and attributed
them to "bistable parameter sets sampled by the sensitivity analysis". Two
things were wrong with that.

1. **Wrong source.** Those numbers come from `cfos_ceiling_check.py`: 27
   bistable sets out of 150 drawn by its own LHS scan. They are not from the
   400-sample GSA in `gsa.py`, which records `sep_cjun` and `sep_gfap` but no
   per-dimer separations at all. Re-ran `cfos_ceiling_check.py` (n_sets=150,
   seed=7) on 2026-10-12 to confirm: 27 bistable sets, fj_ratio median 1.464
   (IQR 1.267-1.765), jj_ratio median 735.909 (IQR 79.894-2585.164), JJ > 3x
   in 100% of sets, fj_frac_on median 0.975, fj_frac_off median 0.724. Saved
   to `ceiling_check_results.csv`.

2. **Ensemble medians presented as reference-set values.** The 97.5% / 72%
   occupancies read as though measured at the reference parameter set. They
   are not. At the reference set, at the midpoint of its own two-state window
   (miR-200a = 0.60), the figures are **95.0% and 84.2%**, and the separations
   are **FJ 1.15x, JJ 14.8x** -- an order of magnitude below the ensemble
   values for JJ. The ensemble is evaluated at miR-200a = 0, where the two
   attractors lie furthest apart and where the reference set is itself
   monostable (both initial conditions converge to cJun_total = 97.6), so the
   two sets of numbers are not comparable without saying so.

Both the manuscript (Section 2.4) and the Figure 1 caption now state the
reference-set values, the ensemble values, the scan they came from, and the
input at which each was measured. The ordering -- the claim the section
actually makes -- is identical either way; only the magnitudes move. Flagging
this because the larger numbers are the more impressive ones and they were the
ones on display.
