# betaEstimator

Standalone study of AMS TOF beta estimators that is intentionally independent of `betaFitter` nonlinear propagation.

The current estimator starts from the four-layer affine timing slope

\[
u_0 = \frac{\sum_i (x_i-\bar{x})(t_i-\bar{t})}{\sum_i (x_i-\bar{x})^2},\qquad x_i=L_i/c,
\]

and studies event information orthogonal to the fitted common time and slope.  The first diagnostic is the per-layer timing residual

\[
r_i=t_i-\bar{t}-u_0(x_i-\bar{x}).
\]

`betaFitter` is used only as a source of the existing train/validation file split and as a physics-performance baseline; this repository does not link or call its nonlinear fitter.

## Layout

- `config/sources.conf`: snapshot of the existing independent train/validation MC split.
- `macro/`: ROOT 6 analysis and plotting macros.
- `results/`: durable generated results.
- `docs/`: method notes.

## Layer residual plots

Run from a clean shell on `ssh:sdiat` with ROOT 6:

```bash
root -l -b -q 'macro/plotLayerResiduals.C()'
```

The macro uses validation files only, requires clean four-layer timing (`tof_bpatt == 4444`), reconstructed charge equal to MC charge, and `0.6 < beta_MC <= 1`.
