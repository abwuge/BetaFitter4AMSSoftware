# Detector-only velocity reconstruction

The training target is now measured four-layer TOF self-consistency. The
MC-supervised inverse-beta error regressor is removed from the active pipeline.
Real flight data from `amsd69n_getrigidity_dt_p000-009` are used for training and
held-out validation. No mass, rigidity, MC charge, MC beta, or MC momentum enters
training. Reconstructed charge chooses the Z=1,2,6,8 model.

## Model and identifiability

Let `x = tof_leng/c` in ns, and `t = tof_tl` in ns. For each event, construct a
2-by-4 orthonormal matrix Q satisfying `Q 1 = Q x = 0`. Then `y = Q t` contains
neither the unknown common time zero nor the event's inverse velocity.

Each layer has 20 measured-input basis functions: a constant, log(deposit/Z^2),
its square, transverse positions, their squares and product, deposit-position
products, and ten bar indicators. Continuous features are winsorized at the
training-only 1st/99th percentiles and standardized with training-only moments.
No fitted speed or measured time is a conditioning feature.

The model predicts four layer mean offsets b and four positive layer noise
scales v from these features. In the two-dimensional contrast space:

```
e = Q t - Q b(X)
S = Q diag(v(X)) Q^T
loss = 0.5 log det S + (nu+2)/2 log(1 + e^T S^-1 e / nu) + ridge
```

A Student-t likelihood with fixed `nu=6` limits sensitivity to timing tails.
The log determinant prevents the trivial solution of inflating all errors.
An L2 coefficient penalty of 0.01 regularizes weakly identified directions;
log variances use a bounded smooth transform. Parameters, regularization, and
sampling rules are fixed before MC evaluation. All 160 parameters per charge
are learned by L-BFGS with an analytic gradient. Nonconverged fits are rejected.

The reconstruction subtracts **only** `Q^T Q b`, then performs a weighted
straight-line fit using `1/v`. It does not apply the unconstrained part of b
parallel to the time-zero or velocity directions. Consequently an injected
`t -> t + a + delta*x` gives `u_reco -> u_reco + delta`, with unchanged features,
on the clean-event reconstruction path.

This is an effective TOF timing-slope estimator. It does not learn an absolute
center-plane energy-loss correction or remove every beta-dependent bias. Such
bias is invisible to these contrasts. Noise correlations parallel to the
slope are likewise unobservable: reconstructing velocity from the fitted
contrast covariance requires the stated diagonal layer-noise model. Better
contrast likelihood is therefore a detector-consistency result, not proof of
an absolute beta-resolution gain.

## Data boundary and splitting

- Extraction uses an explicit allowlist of measured branches and stores no truth.
- Training rejects caches containing extra fields, including historical truth
  arrays. Truth-based training cuts and true-charge matching are removed.
- Measured charge, valid Tracker L2-L8 deposits/positions, and center-crossing
  geometry define the extracted sample. Clean training additionally requires
  four TOF layers (`4444`), finite paths/times, positive measured slope, valid
  positions/bar IDs, and positive deposits.
- All 439 existing flight files are assigned by SHA256(run) modulo 5: 341
  training, 98 validation. Splitting is independent of measured or true speed.
  Duplicate `(run,event)` identities and run overlap are fatal errors.
- At most 120,000 clean events per charge are sampled uniformly from training
  with a fixed seed. All available clean heavy-nucleus events are used if fewer.
- Held-out metrics include contrast negative log likelihood and leave-one-layer-out
  timing prediction widths. No held-out timing or MC truth is used to fit parameters.
- Incomplete or invalid inference events fall back to the existing producer
  linear beta. They are retained in evaluation; fallback counts are reported.

The extraction audit establishes that all required branches can be read, not
that every original production job completed. This source snapshot is the
existing-file subset, not a full-production completeness claim.

## Final evaluation

`write_comparison.py` verifies model hashes against the training summary before
opening MC truth. For each MC source it makes detector-only predictions first,
then reads truth for evaluation. Historical Hybrid rows are aligned using raw
source-entry selection, reconstructed charge, true charge, producer linear beta,
and true beta. The same final common-event cuts and true-beta bins are used for
both methods. No new model fit or hyperparameter search occurs at this stage.

The existing ROOT macro measures the Gaussian core of
`1/beta_reco - 1/beta_MC` in eight bins from 0.60 to 1.00 per charge. MC truth
uses the historical checkpoint-10 center reference of the Hybrid baseline.
Data-trained results are written to `results/data_trained_comparison`; the old
supervised results in `results/comparison` remain historical only.

MC evaluation after flight training tests transfer to that MC detector model.
A worse MC score can reflect a model limitation or a data/MC response mismatch;
it does not establish which cause dominates. Conversely a good score does not
prove unbiased isotope reconstruction on flight data. An independent physical
reference is needed for the absolute velocity scale.

## Verification

Automated checks cover the analytic likelihood gradient, slope/time-zero
invariance, projection of the mean correction, model serialization, invalid-hit
fallback, batched inference, rejection of truth-bearing caches, and reconstruction on an independent
synthetic heteroscedastic sample. Synthetic latent velocities generate test
measurements and are used only to score the held-out reconstruction.

## Flight-data result (2026-09-19)

All 439 sources were read successfully (134,666,562 raw rows), with no duplicate
selected events or overlap of training and validation runs. The optimizer
converged for all four charges. The following detector-only holdout results
were recorded before opening MC truth:

| Z | Clean training pool | Used for fit | Clean holdout | NLL decrease [nats/event] | Leave-one-layer-out width decrease across S1-S4 |
|---:|---:|---:|---:|---:|---:|
| 1 | 22,595,003 | 120,000 | 6,328,322 | 0.07991 | 1.13% to 2.59% |
| 2 | 5,311,577 | 120,000 | 1,529,784 | 0.30396 | 7.42% to 12.61% |
| 6 | 149,227 | 120,000 | 43,105 | 0.71208 | 21.31% to 23.72% |
| 8 | 113,984 | 113,984 | 32,593 | 0.82104 | 25.95% to 26.74% |

These widths are central-68% half-widths of **time prediction errors**, not
velocity residual widths. Each layer's time is predicted from the other three
measured times, with model parameters fixed from independent training runs.
The model features contain no time, so the omitted time does not leak into its
prediction through the learned mean or weights. All four layers improve at the
point-estimate level for every charge.

The full extracted holdout, including fallback events, has 8,663,639 / 1,632,795 /
44,824 / 33,904 events for Z=1/2/6/8. Fallback counts are 2,335,368 / 103,052 /
1,719 / 1,311. The exported flight ROOT trees retain all these events and mark
whether the model was applied; the timing metrics above concern clean input
events only. The original calibration's absolute velocity scale remains the
reference, and these checks cannot measure an unobservable common slope bias.

## Frozen-model MC result (2026-09-19)

All 96 historical Hybrid validation partitions aligned. Common samples reproduce
193,655 / 118,554 / 45,894 / 24,747 events for Z=1/2/6/8, including fallbacks.
All 32 Gaussian-core fits produced finite positive widths for both methods.

The data-trained model **does not give a consistent MC velocity-resolution gain**.
In the 0.95-1.00 true-beta bin:

| Z | Hybrid inverse-beta core sigma | Data-trained sigma | Narrowing vs Hybrid |
|---:|---:|---:|---:|
| 1 | 0.036259 | 0.036098 | +0.44% |
| 2 | 0.017718 | 0.017593 | +0.71% |
| 6 | 0.011117 | 0.011458 | -3.06% |
| 8 | 0.010923 | 0.011458 | -4.90% |

Negative narrowing means a wider distribution. Z1/Z2 gains here are small; Z6/Z8
are about 3.1%/4.9% wider. Across the eight beta bins, Z1 improves at the point
estimate in five bins, Z2 in four, and Z6/Z8 in none. Low-statistics bin errors
are large; no per-bin significance claim is made. The old supervised model's
8%-13% high-beta gains must not be attributed to this estimator.

Absolute bias is also not uniformly removed. At beta 0.95-1.00 the data-trained
core means in inverse beta are +0.005177 / +0.002294 / -0.000074 / +0.000393 for
Z1/Z2/Z6/Z8. Heavy-nucleus mean bias improves in several bins, whereas proton
bias is generally larger than Hybrid. The full per-bin means and fit errors are
in `results/data_trained_comparison/per_bin_summary.csv`.

The demonstrated improvement is flight-data timing self-consistency. Neither
that improvement nor the MC transfer result establishes an improved absolute
center-beta resolution on flight data. The present diagonal noise model,
unobservable slope-error components, physical energy loss, and differences in
flight/MC response are possible limitations; this comparison does not separate
their contributions. Model parameters were not retuned after opening MC truth.
An independent measured velocity reference or independently constrained detector
response is needed to test the missing absolute-scale information.

The ROOT output contains 44 pages: eight residual-fit pages per charge and four
pages each for core sigma, core mean and relative sigma change. Fit bands show
fit-parameter errors, not training or detector-model systematic uncertainties.
