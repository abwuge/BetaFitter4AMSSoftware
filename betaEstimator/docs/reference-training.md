# Independent-reference TOF calibration study

This study did not establish a method that improves AMS-center inverse-beta
Gaussian-core resolution for all four charges and all eight beta bins.
A near-light-speed RICH calibration improves the independent flight
TOF-minus-RICH distribution, but fails to transfer uniformly to the final MC.
It is a research implementation, not a replacement for the active estimator.

Training, feature selection, regularization and correction strength used flight
measurements only. All models and the selected candidate were frozen at
2026-09-19 08:41:57 UTC. Final MC evaluation started at 08:42:26 UTC.
No model parameter or model choice changed after truth evaluation.

## Why the contrast objective is insufficient

With x=tof_leng/c, u=1/beta, and t=t0*1+x*u+b+epsilon, the current two-row Q obeys
Q*1=Q*x=0. Its likelihood observes Q*b and Q*V*Q^T. Projecting the fitted mean
back with Q^T*Q leaves the ordinary equal-weight slope exactly unchanged.
The original estimate changes through unequal layer weights and their
interaction with that projected mean.

At fixed geometry, even independent layer variances have one unidentifiable
direction. Let n_i=product_{j!=i}(x_i-x_j). Then
Q*diag(n)*Q^T=0. Four variances cannot generally be recovered from the three
independent entries of a two-dimensional covariance. Shared parameters and
varying geometry can reduce this degeneracy, but regularization is not proof
of identification. LOO layer-time errors contain the same two residual degrees
of freedom and do not restore the missing slope information.

The numerical test uses x=(-1,-0.92,0.92,1),
v_A=(1,1,1,1), v_B=(0.1,1.828,0.172,1.9). These produce the same contrast
covariance. GLS using v_B increases the actual slope width under v_A by 32.2%.
This is an algebraic counterexample, not a MC tuning experiment.

## Implemented candidates

1. Independent measured reference:
   u_hat=u_OLS-f(geometry, bars, measured deposits/Z^2, orthogonal time residuals).
   A bounded polynomial/interaction basis and Huber-weighted ridge fit target
   u_OLS-1/beta_RICH. Inference never uses RICH, rigidity, a mass hypothesis,
   fitted TOF beta, or absolute event time as correction features.
2. Identifiable noise correction:
   a0=x_centered/sum(x_centered^2), g_n=Q*diag(n)*a0, and h perpendicular to g_n.
   With S=Q*diag(v)*Q^T and g=Q*diag(v)*a0, use
   k=h*(h^T*g)/(h^T*S*h) and u_hat=u_OLS-k^T*(Q*t-Q*b).
   Its conditional variance improvement is independent of the missing variance
   direction, under the diagonal-noise and correct conditional-mean assumptions.
   This guarantee does not cover correlated noise, energy-loss curvature or
   an incorrect noise/mean model.
3. A measured-rigidity R>100 GV, u approximately 1 calibration was kept as an
   explicitly approximate cross-check. It uses the producer default rigidity,
   without a separately aligned fit-quality cut; it was not selected.
   R and the approximate mass-to-charge argument are never inference inputs.

The RICH fit and inner selection require measured beta_RICH>=0.99. This reduces,
but does not exactly remove, the physical speed difference from center to RICH.
No low-beta center-to-radiator transport model was fitted. Applying the learned
correction at lower beta is an extrapolation test and failed for heavy nuclei.

Modes geometry/full, ridge 0.01/0.1/1, and shrinkage 0/0.25/0.5/0.75/1 were
compared on independent inner flight runs. All charges selected full, ridge
0.01, shrinkage 1 before MC access. Parameters were not refitted on holdout runs.

## Reference quality and source audit

All 439 flight sources in config/data_sources.conf were read. Its existing
341/98 train/holdout whole-run split was retained; an additional hash of the
training run defines the inner selection split. Per-source light-charge
sampling caps the broad reference/anchor pool at 3000 events; all eligible
heavy-charge events are retained. These are sampled reference measurements,
not the entire flight event population.

The RICH reference must satisfy irich>=0, itrtrack>=0 and
rich_itrtrack==itrtrack. Without these checks, stale good/clean/error fields
can survive from the preceding event in the producer.
Quality uses only RICH reconstruction/PE/geometry checks and measured detector
cleanliness. The pre-existing SelectRich TOF-RICH agreement cut is not used.
BetaConsistency is the difference between two RICH algorithms; its -1 missing
sentinel is rejected. Full quality thresholds are implemented in
analysis/train_reference.py::quality and traced in reference-audit.md.

| Z | RICH fit events | Inner selection events | Outer RICH holdout |
|---:|---:|---:|---:|
| 1 | 100000 | 83602 | 123487 |
| 2 | 100000 | 73860 | 154973 |
| 6 | 14038 | 4985 | 10727 |
| 8 | 11424 | 4099 | 8316 |

RICH is downstream of lower TOF. Audited flight radiator points are near
z=-74.705 cm (aerogel) and -75.552 cm (NaF); they are not L5.
The high-beta restriction is a calibration approximation, not exact
identification of AMS-center beta. RICH systematics and shared Tracker
trajectory errors remain possible correlated label errors.

## Independent flight result

For 0.99<=measured beta_RICH<=1.02, use the central-68% half-width W68 of
u_TOF-u_RICH. Values above measured beta=1 are retained as measurement
fluctuations. W68 is not a Gaussian sigma and is not pure TOF resolution.

| Z | Holdout N | OLS W68 | New W68 | Narrowing | Paired source-bootstrap 95% |
|---:|---:|---:|---:|---:|---:|
| 1 | 88671 | 0.0436424 | 0.0424076 | 2.83% | 2.40% to 3.16% |
| 2 | 80066 | 0.0208149 | 0.0190772 | 8.35% | 7.79% to 9.00% |
| 6 | 5310 | 0.0145709 | 0.0106841 | 26.68% | 24.36% to 28.99% |
| 8 | 4231 | 0.0144745 | 0.0096719 | 33.18% | 30.93% to 36.26% |

Relative to the old NLL estimator, the corresponding narrowing is
0.60%, 7.00%, 27.43% and 35.04%.
Intervals hold the fitted models fixed and exclude retraining/model-selection
uncertainty, reference systematics and detector transfer uncertainty.
RMS/tails and all broader reference bins remain in the JSON outputs;
large light-charge tails are not removed to make the W68 numbers look better.

## Frozen final MC evaluation

Historical common counts exactly reproduce 193655/118554/45894/24747.
All 96 Hybrid partitions were aligned. A further 24 source files per charge,
excluded from the previous MC training/validation manifest, were selected by
filename hash before truth access. Their common counts are
194772/117794/46594/24010, with no duplicate selected events or overlap with
the historical evaluated events.

Historical truth uses the existing checkpoint-10 convention for direct Hybrid
comparability. Fresh truth is the primary L5 MC cluster:
mtrpri[4]==1, mtrmom[4]>0, abs(mtrcoo[4,2]-1.698)<0.03 cm, combined with the
MC particle mass. This is a layer-center reference cut, not a material interval.

Fresh-MC Gaussian-core result, 0.95<=beta_L5<=1:

| Z | Linear sigma | New sigma | Narrowing | New core mean |
|---:|---:|---:|---:|---:|
| 1 | 0.0363460 | 0.0355533 | +2.18% | +0.003066 |
| 2 | 0.0177108 | 0.0190544 | -7.59% | +0.002241 |
| 6 | 0.0110577 | 0.0114667 | -3.70% | +0.003712 |
| 8 | 0.0107231 | 0.0114174 | -6.48% | +0.003561 |

Negative narrowing means degradation. Over eight 0.05-wide bins from 0.60 to
1.00, improvement point estimates occur in 5/0/0/1 bins for Z1/Z2/Z6/Z8,
on both historical and fresh MC. All 64 method fits on each sample are finite
and positive. Heavy-charge maximum absolute new core means across the fresh
bins reach 0.04725 (Z6) and 0.03322 (Z8). Low-beta extrapolation is unacceptable.

Historical high-beta narrowing versus Hybrid is +2.48%/-7.94%/-3.94%/-3.66%.
The independent result therefore does not rescue the old all-charge claim.

The variance diagnostic checks
Var(e-c)=Var(e)+Var(c)-2*Cov(e,c) on the full fixed high-beta sample.
For fresh Z2 at beta>=0.99, correction variance is 1.27494e-4 while twice the
error-correction covariance is only 8.91770e-5, so variance necessarily rises.
For Z6/Z8 full variance falls substantially while the Gaussian core widens:
tail improvement must not be mislabeled as improved core resolution.
These checks diagnose outcomes; they were not used to alter any model.

## Reproduction and retained artifacts

Run in /nas01/home/hxwu/AMSSoft/betaEstimator:

~~~bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 analysis/extract_reference.py --output /tmp/betaEstimator-reference
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 analysis/augment_reference.py /tmp/betaEstimator-reference
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 analysis/train_reference.py
# Only after the preceding step creates the frozen summary:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 analysis/evaluate_reference.py
root -l -b -q 'macro/plotHybridComparison.C("results/reference_comparison/historical","results/reference_comparison/historical/plots","Hybrid","RICH-calibrated")'
root -l -b -q 'macro/plotHybridComparison.C("results/reference_comparison/fresh","results/reference_comparison/fresh/plots","Linear","RICH-calibrated")'
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m unittest discover -s tests -v
~~~

The reference trainer retains every frozen candidate and the inner selection
metrics in results/reference_training. Final MC ROOT/CSV/plots are in
results/reference_comparison. The existing estimator and betaFitter physics
baseline remain available; no new result is promoted as uniformly superior.
The ROOT macro now accepts optional series labels and branch names; its original
defaults remain unchanged.

Seventeen tests passed across the existing likelihood checks, independent
reference model, affine response, forbidden feature access, label-independent
feature scaling, serialization, partial-layer producer fallback, and variance
nullspace/risk algebra. Synthetic latent values are used only to score frozen
test models and are not detector MC training labels.

The next scientifically justified work is independent center-speed/transport
constraints and data-MC response closure. Merely selecting another frozen
candidate after reading this MC result would violate this study's contract.
