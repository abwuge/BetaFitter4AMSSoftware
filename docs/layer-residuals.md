# Layer timing residual diagnostics

The standalone estimator defines `x_i = L_i/c` and first fits the four TOF times with an affine line. The plotted layer residual is

\[
r_i=t_i-\bar t-u_0(x_i-\bar x).
\]

The plots use the existing independent validation files only, with reconstructed charge equal to MC charge, complete Tracker L2-L8 information, clean four-layer timing (`tof_bpatt=4444`), finite positive fitted slope, and `0.6 < beta_MC <= 1`.

## Current Gaussian core widths

Values below are in ns.

| Z | S1 | S2 | S3 | S4 |
|---:|---:|---:|---:|---:|
| 1 | 0.10093 | 0.10342 | 0.10253 | 0.10007 |
| 2 | 0.04523 | 0.04639 | 0.04691 | 0.04572 |
| 6 | 0.01963 | 0.02013 | 0.02117 | 0.02069 |
| 8 | 0.01736 | 0.01783 | 0.01943 | 0.01890 |

The residual distributions already show a strong charge dependence. S3/S4 are modestly broader than S1/S2 for the heavy nuclei.

The profile plots also show structured dependence on true beta, measured TOF deposit, and along-bar hit position. These structures are diagnostics of predictable event-level error; they are not yet assigned uniquely to time-walk, geometry, path reconstruction, or energy-loss fluctuations.

Generated outputs are under `results/layer_residuals/` and are intentionally ignored by Git so regenerated plots do not clutter source history.
