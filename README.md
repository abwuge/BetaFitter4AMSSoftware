# betaEstimator

Archived on 2026-11-01 as a code module of the velocity-reconstruction work,
together with `../betaFitter`. Further velocity-reconstruction research is stopped.
The unified source snapshot and both histories are published in the existing
BetaFitter4AMSSoftware repository on branch `archive/velocity-reconstruction`.

AMS TOF beta reconstruction from flight measurements. RICH is an independent
validation reference; RICH beta and MC truth do not enter the current training.

The current carbon study covers 2.5–200 GV with three independently trained
mass-constrained models. Training and validation use disjoint actual runs.

- Train/evaluate: `analysis/flight_beta_wide.py`.
- Frozen models and paired predictions: `results/flight_beta_wide/`.
- Rigidity/beta curves: `analysis/plot_flight_beta_wide.py`.
- Per-bin residual distributions: `analysis/plot_flight_beta_distributions.py`.
- Larger NEWTOF validation: `analysis/prepare_tofcalib_validation.py` and
  `analysis/tofcalib_validation_worker.py`; evaluate with
  `analysis/evaluate_tofcalib_validation.py`, verify with
  `analysis/verify_tofcalib_validation.py`; outputs in `results/tofcalib_validation/`.

The NEWTOF validation preserves existing TOF measurements and supplements five
missing native measurement branches. The completed 4096-file holdout contains
907830 independent RICH reference events, with 48 rigidity bins and 29 beta bins.
The frozen model gives an overall beta-residual W68 reduction of 2.9%; the
2.5–6 GV region deteriorates by 3.0%, so the improvement is not uniform.
Numerical outputs, models and necessary logs are retained; routine runs do not
require generated reports.

Earlier approaches remain under `results/`; their method descriptions are in
`docs/`. `betaFitter` remains the comparison baseline.
