# New estimator vs physical Hybrid

**Historical supervised result (2026-09-13).** This document and
`results/comparison/` describe the former MC-truth-trained regressor. That
training method was replaced on 2026-09-19. These gains do not apply to the
current flight-trained model; see [detector-only training](detector-only-training.md).

## Scope

This is the main betaEstimator performance comparison. It does not plot the total beta spectrum and it does not use scatter plots.

For every charge Z=1,2,6,8 and every 0.05-wide true-beta bin from 0.60 to 1.00, the same validation events are used for both methods. The plotted residual is

\[
\Delta(1/\beta)=1/\beta_{reco}-1/\beta_{MC}.
\]

Each bin is fitted independently. The Gaussian core window is derived from the residual median and IQR, using `median ± 3 * 0.741301109252801 * IQR`.

## Sample closure

The standalone training counts exactly reproduce the previous validated study:

- Z1: 191592
- Z2: 117906
- Z6: 45311
- Z8: 24133

The common validation sample also reproduces the previous counts:

- Z1: 193655 total, 185283 with beta >= 0.6
- Z2: 118554 total, 112598 with beta >= 0.6
- Z6: 45894 total, 44729 with beta >= 0.6
- Z8: 24747 total, 24299 with beta >= 0.6

All 96 Hybrid validation partitions were aligned against their raw source rows before training or comparison.

## Resolution result

The new estimator has a smaller fitted core sigma than Hybrid in all 32 Z-by-beta bins at the point-estimate level.

The improvement `(1 - sigma_new/sigma_hybrid) * 100%` spans:

- Z1: 1.46% to 10.33%
- Z2: 2.08% to 8.37%
- Z6: 2.87% to 10.66%
- Z8: 0.60% to 14.55%

The lowest-statistics heavy-nucleus bins have correspondingly large fit uncertainties, so individual low-statistics points should not be over-interpreted. The stable high-beta region is clearer. In the 0.95-1.00 bin the improvement is:

| Z | Hybrid sigma | New sigma | Improvement |
|---:|---:|---:|---:|
| 1 | 0.036259 | 0.032512 | 10.33% |
| 2 | 0.017718 | 0.016242 | 8.33% |
| 6 | 0.011117 | 0.009932 | 10.66% |
| 8 | 0.010923 | 0.009495 | 13.07% |

The full values are in `results/comparison/per_bin_summary.csv`.

## Bias result

For Z1 and Z2, the new estimator also moves the fitted residual mean substantially closer to zero across most of beta > 0.6.

For Z6 and Z8 the picture is more nuanced: the resolution improves, but some intermediate-beta bins become more negative than Hybrid. Examples:

- Z6, 0.75-0.80: Hybrid mean about -0.00007, New about -0.00069.
- Z8, 0.75-0.80: Hybrid mean about -0.00017, New about -0.00094.

At beta 0.95-1.00 the New means are small for all four charges: approximately -0.00017, -0.00004, -0.00019, -0.00015 for Z1/Z2/Z6/Z8.

Therefore the current conclusion is: **the new estimator gives a consistent resolution gain, but heavy-nucleus bias still needs bin-by-bin control rather than relying on a single beta>0.6 average.**

## Main figures

- `results/comparison/plots/residual_fits_Z*.pdf`: eight pages per charge, one page per beta bin, Hybrid and New histograms plus fitted core curves.
- `results/comparison/plots/core_sigma_vs_beta.pdf`: direct resolution comparison as line graphs.
- `results/comparison/plots/core_mean_vs_beta.pdf`: direct bias comparison as line graphs.
- `results/comparison/plots/relative_improvement_vs_beta.pdf`: direct relative resolution gain as a line graph.

The shaded bands on summary plots are the Gaussian fit parameter errors only. They do not include retraining, model-selection, MC-production, isotope, or data/MC systematic uncertainties.
