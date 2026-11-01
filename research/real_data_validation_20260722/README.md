# Real-data validation of linear and nonlinear TOF beta (2026-07-22)

Historical validation using AMS real-data ntuples. Research sources, this note, and selected numerical summaries are versioned; event-level outputs and generated plots remain in the work directory.

## Existing PDF audit

`macro/results/beta_comparison/real_data_beta_comparison.pdf` uses `betaRigidity` as the reference. In `src/Util.cc`, that quantity and the nonlinear fit both assume `mass = 2*Z*0.9314941 GeV` for every real event. This is invalid for mixed isotope data (especially Z=1, which is treated as deuterium, and Z=4, which is treated as A=8). The PDF also mixes charges unless an external restriction was supplied, and `results/10.root` lacks event IDs and all RICH/quality variables.

The merged `results/10.root` cannot be friended to the 194 input files by entry order: although both contain 11,112,424 entries, 107,854 of 111,013 sampled entries had `linearBeta != tof_betah`. Future outputs must retain `run`, `event`, source-file ID, and local entry.

## Exact RICH cross-check

The repository's own `build/bin/betaFitter` was run on:

`/lustre02/user/yachen/public/amsd69n_B1130_Data_TofCalibN2/1306976244_0.root`

The resulting 57,860-entry beta tree was compared event-by-event with the RICH variables in the same source file. Selection was exploratory: `rich_good`, `rich_clean`, at least four used ring hits, ring probability >0.01, good TOF geometry, positive rigidity, and valid beta values.

For 25,644 selected events:

- All RICH: robust 68% width in delta(1/beta) was 0.02262 linear and 0.02269 nonlinear.
- NaF: 0.02249 linear and 0.02293 nonlinear.
- NaF beta 0.70-0.80: 0.1303 linear and 0.1252 nonlinear, a ~4% core improvement.
- The fraction beyond 0.5 in delta(1/beta) increased from 0.031% to 0.094% overall.

Direct TOF-RICH means are not absolute bias measurements until both beta values are propagated to a common reference plane. Narrow-bin widths and tails remain useful comparisons because RICH resolution is much better than TOF resolution.

## TOF leave-one-out prediction

For 7,008 real events with 0.4 < built-in TOF beta < 0.9 and good geometry, the repository's exact `BetaNL` was fitted using three TOF times and used to predict a withheld time. S1 was not withheld because current `BetaNL` initializes its time offset from `_hitTime[0]`; setting S1 to -1 makes that test invalid.

Using S2-S4 as holdouts:

- Overall robust 68% width: 0.09854 linear, 0.09855 nonlinear.
- beta 0.4-0.6: 0.07572 linear, 0.07334 nonlinear (~3% improvement).
- beta 0.6-0.7: 0.07007 linear, 0.06961 nonlinear (<1%).
- beta 0.7-0.8: identical within 0.001%.
- beta 0.8-0.9: nonlinear slightly worse.

Thus real-data out-of-sample timing prediction does not reproduce the ~50% MC resolution improvement.

## Four-layer blind model validation

The main real-data validation now keeps all four TOF measurements. Linear and nonlinear models both fit exactly two event parameters: inverse beta and a common time offset. The energy-loss scale is fixed before validation, so the models have equal event-level parameter counts.

The 194 real-data files were split chronologically:

- first 65 files: training only;
- last 129 files: blind validation;
- training: 3,160,781 input entries, 238,034 selected events, deterministic 1/5 prescale (47,599 events) for the zeta scan;
- validation: 7,951,643 input entries, 600,204 selected events.

Selection: `tof_qs==1111`, four valid TOF times/deposits/path lengths, Z in {2,6,8}, and `0.4 < tof_betah < 0.9`. A diagnostic core sample requires the four-layer linear chi-square to be below 10; 594,405 validation events pass.

Before the low-beta test, median layer residuals were estimated from an independent high-beta control sample (`0.97 < tof_betah < 1.03`) separately in the training and validation periods. The layer patterns agree between periods at approximately the 1 ps level. An otherwise identical uncalibrated run was also made.

For fixed zeta=1.9 in the blind core sample:

- all Z: mean four-layer timing SSE decreased by 0.270%; 59.1% of events improved;
- Z2: 0.105% SSE decrease; 55.9% improved;
- Z6: 0.569% SSE decrease; 61.9% improved;
- Z8: 0.849% SSE decrease; 63.5% improved.

The effect is concentrated where energy loss should matter. For zeta=1.9:

| Z | beta range | events | SSE change | fraction improved |
|---:|:----------:|-------:|-----------:|------------------:|
| 2 | 0.40-0.55 | 2,741 | -3.86% | 66.9% |
| 6 | 0.40-0.55 | 1,096 | -14.37% | 78.3% |
| 8 | 0.40-0.55 | 656 | -18.27% | 82.6% |
| 2 | 0.55-0.65 | 15,496 | -0.93% | 61.6% |
| 6 | 0.55-0.65 | 7,034 | -3.65% | 70.3% |
| 8 | 0.55-0.65 | 4,881 | -6.15% | 74.1% |

The same low-beta/high-Z hierarchy persists without subtracting the high-beta layer baseline: in beta 0.40-0.55, the SSE reductions are 2.93% (Z2), 12.10% (Z6), and 16.64% (Z8).

The improvement is stable in time. Among validation files with at least 20 core events, mean delta-chi-square is negative in 126/127 files overall, 126/127 for Z2, 126/127 for Z6, and 127/127 for Z8.

A curvature projection removes the two nuisance directions corresponding to a common time offset and an overall constant-beta slope. It shows an important identifiability limitation: after refitting beta, the zeta=1 nonlinear curvature template has only a few-ps RMS at beta 0.40-0.55 (about 2.1 ps for Z2, 4.7 ps for Z6, and 5.5 ps for Z8), and becomes sub-ps at higher beta. Consequently, residual time-walk and layer calibration effects can dominate an internally fitted zeta. The low-beta training scan prefers an effective scale around 5-7, but this should not be interpreted as the physical energy-loss scale or used to replace the production value: larger scales improve a minority of events very strongly while worsening tails.

For zeta=1.9, the nonlinear fit also shifts the reconstructed beta substantially relative to the exact four-layer linear fit. In beta 0.40-0.55 the median beta shifts are +0.0199 (Z2), +0.0419 (Z6), and +0.0478 (Z8). Therefore the four-layer timing likelihood demonstrates a real nonlinear timing structure, but cannot by itself prove that the shifted initial beta is more accurate. That still requires an independent common-plane reference such as RICH or known-isotope mass peaks.

Reproducible source:

`src/four_layer_model_validation.cpp`

Primary outputs:

`results/four_layer_blind_final/`

Uncalibrated robustness outputs:

`results/four_layer_blind_uncalibrated/`

## AMS-center beta with MC global zeta (2026-07-28)

The blind four-layer validation was repeated with the current
`BetaNL::propagate(AMSCenter)` timing definition.  The fixed zeta values came
from full-statistics MC measured-time global profiles, not from the real-data
training or validation split:

| Z | MC stable events | measured-time global zeta | checkpoint-truth global zeta |
|---:|-----------------:|--------------------------:|-----------------------------:|
| 2 | 1,217,863 | 3.25872 | 4.45434 |
| 6 | 453,395 | 2.83809 | 4.40118 |
| 8 | 332,367 | 2.60123 | 4.29041 |

The real-data sample, chronological 65/129-file blind split, high-beta layer
calibration, and event selections are unchanged.  In the 594,405-event blind
core sample, the center/global-zeta model reduces the mean four-layer timing
SSE by 0.403%, versus 0.270% for zeta=1.9.  The effect remains concentrated at
low beta and high charge:

| Z | beta range | events | SSE reduction | fraction improved | median center-beta shift |
|---:|:----------:|-------:|--------------:|------------------:|-------------------------:|
| 2 | 0.40-0.55 | 2,741 | 5.97% | 65.6% | +0.000062 |
| 6 | 0.40-0.55 | 1,096 | 21.27% | 75.2% | +0.000231 |
| 8 | 0.40-0.55 | 656 | 25.82% | 80.5% | +0.000310 |
| 2 | 0.55-0.65 | 15,496 | 1.59% | 61.2% | +0.000012 |
| 6 | 0.55-0.65 | 7,034 | 5.22% | 69.0% | +0.000060 |
| 8 | 0.55-0.65 | 4,881 | 8.11% | 73.0% | +0.000082 |

The mean improvement is larger than at zeta=1.9, while the fraction of events
improved is slightly smaller.  Thus the larger scale increases the gain for
the favorable events rather than improving more events.  The all-core mean
delta-chi-square is negative in 126/127 validation files with at least 20
events.  Unlike the old before-TOF beta parameter, the fitted center beta moves
only by O(1e-4) in the low-beta bins; most of the nonlinear effect is absorbed
as the physically expected change of reference plane rather than a large beta
shift.

Independent checks do not show a beta-resolution gain.  On the exact single-
file RICH comparison restricted to Z={2,6,8}, the robust 68% width of
delta(1/beta) changes from 0.018644 (linear) to 0.018674 (nonlinear), while the
RMS changes from 0.04342 to 0.06761 and the fraction beyond 0.5 increases from
0.025% to 0.076%.  Direct means remain non-absolute because RICH and TOF beta
are not propagated to one common plane.  In the S2/S3 leave-one-out timing
test, the overall robust width changes only from 0.08257 ns to 0.08231 ns; the
0.4-0.6 bin is slightly worse (0.06647 ns to 0.06756 ns).  S4 is not a valid
holdout in center mode because the current Minuit time-offset initialization
uses both S1 and S4.

These results support the existence of the nonlinear four-layer timing
curvature, but not an independent improvement in reconstructed beta.  The
large separation between measured-time and checkpoint-truth global zeta also
shows that the fitted center/global scale is an effective model parameter, not
a universal physical energy-loss multiplier.

Primary outputs:

- `results/four_layer_blind_center_global/`
- `results/rich_crosscheck_center_global.txt`
- `results/tof_leave_one_out_center_global.txt`

## Blind disjoint-pair TOF closure (2026-07-28)

The tracker rigidity is not an independent low-beta reference: the rigidity
used by `betaFitter` is a GBL inner-track state at the first inner measurement
plane, and its reconstruction has its own material and resolution model.  A
new data-only test therefore uses two disjoint sets of TOF times rather than a
tracker-derived mass peak or a leave-one-out fit.

For every event, center inverse beta is solved independently from S1--S3 and
S2--S4.  These estimators share no timing hit.  The primary closure variable
is

`delta_13_24 = 1/beta(S1,S3) - 1/beta(S2,S4)`.

The complementary comparison uses S1--S4 and S2--S3.  Linear, center with
zeta=1.9, center with the fixed per-charge MC global zeta, and all null models
use exactly the same events.  The first 65 files supply the high-beta layer
calibration; that calibration is frozen before evaluating the last 129 files.
The validation-period high-beta residuals are written only as a drift
diagnostic and are not applied.  The primary sample has no four-layer SSE cut;
the old linear-chi2<10 sample is retained only as a labelled diagnostic.

All 600,204 selected blind events have valid solutions for both disjoint-pair
comparisons under every tested model.  In the most sensitive built-in-beta
range, 0.40--0.55, the primary closure is:

| Z | events | mean: linear | mean: center zeta=1.9 | mean: center global zeta | median: linear -> global | sigma68: linear -> global | median abs(delta) reduction | global-zeta events with smaller abs(delta) |
|---:|-------:|-------------:|----------------------:|-------------------------:|--------------------------:|----------------------------:|----------------------------:|--------------------------------------------:|
| 2 | 2,771 | -0.01496 | -0.01248 | -0.01063 | -0.01103 -> -0.00697 | 0.02546 -> 0.02508 (-1.48%) | 7.29% | 64.7% |
| 6 | 1,115 | -0.02453 | -0.01905 | -0.01603 | -0.01609 -> -0.00826 | 0.02053 -> 0.01937 (-5.64%) | 27.16% | 75.7% |
| 8 | 662 | -0.02332 | -0.01681 | -0.01404 | -0.01937 -> -0.01142 | 0.02363 -> 0.02191 (-7.28%) | 29.40% | 80.7% |

The 95% absolute tail decreases by 5.55%, 12.44%, and 17.12% for Z2, Z6,
and Z8.  Center/global zeta also outperforms center/zeta=1.9 in aggregate
centering, core width, and median absolute closure, particularly for Z6 and
Z8.  However, zeta=1.9 has a slightly larger per-event win fraction.  The
larger global zeta therefore helps favorable events more strongly rather than
helping more events, matching the pattern in the four-layer SSE test.

The correction has the required charge, beta, energy-deposit, and direction
structure:

- From beta 0.40--0.55 to beta 0.85--0.90, the mean center/global correction
  shrinks from 0.00433/0.00850/0.00928 to
  0.000070/0.000146/0.000166 for Z2/Z6/Z8; the high-beta width change is at
  or below 0.02%.
- Within beta 0.40--0.55, the Z6 mean correction grows from 0.00713 at
  upper edep/Z^2=5.5--6.5 MeV to 0.01169 above 6.5 MeV.  For Z8 it grows
  from 0.00765 at 5.0--5.5 MeV to 0.01069 at 5.5--6.5 MeV.
- Reversing the energy-loss sign moves the low-beta mean in the wrong
  direction: to -0.01914, -0.03226, and -0.03136, and broadens the core.
  In every validation file with at least five low-beta events, the physical
  model's shift has the expected sign (119/119 Z2, 112/112 Z6, 80/80 Z8).
- Jointly shuffling S2/S3 deposits among events matched in Z and 0.01-wide
  beta bins changes the low-beta mean by at most 6.3e-5.  Swapping S2 and S3
  deposits also changes little.  Thus the test validates mainly the average
  charge/beta-dependent curvature and its sign; it does not show that
  event-by-event deposit fluctuations improve beta.

The S1--S4 versus S2--S3 comparison is less discriminating because the
central pair is short and the upper/lower energy changes partly cancel.  Its
low-beta robust median moves closer to zero for all three charges and its
sigma68 improves by 1.4%, 3.7%, and 6.6%, but its ordinary mean is not
uniformly improved for Z2 and Z6.  It is supporting, not primary, evidence.

This is stronger real-data evidence than the four-layer SSE alone that the
physical-sign nonlinear propagation model is preferred at low beta, and it
supports the current global zeta over 1.9 as a model ranking.  It still does
not establish absolute center-beta accuracy or the full four-layer beta
resolution.  Residual layer calibration, common-mode timing effects, the
A=2Z mass assumption, and the use of built-in TOF beta for binning remain.
An independent low-beta velocity or momentum reference would still be needed
for an absolute claim.

Reproducible source and outputs:

- `src/disjoint_pair_validation.cpp`
- `results/disjoint_pair_blind_center_global/pair_closure_summary.csv`
- `results/disjoint_pair_blind_center_global/pair_closure_file_stability.csv`
- `results/disjoint_pair_blind_center_global/high_beta_layer_calibration.csv`

## Recommended references

1. RICH beta in NaF/Agl overlap, after common-plane propagation.
2. Three-layer fit predicting the fourth TOF time (data-only, no mass or MC).
3. Known-isotope mass peaks from tracker rigidity and TOF beta, with isotope tags independent of the tested TOF beta.
4. High-rigidity beta approximately 1 and upward/downward symmetry for velocity-scale calibration.
5. Geomagnetic-cutoff selected near-pure isotope samples, following the AMS lithium-isotope analysis.
