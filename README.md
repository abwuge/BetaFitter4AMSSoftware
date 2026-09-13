# betaEstimator

Standalone AMS TOF beta-estimator study, intentionally independent of `betaFitter` nonlinear propagation.

The main analysis compares the new estimator directly with the current physical-boundary Hybrid on the **same validation events and the same true-beta bins**. `betaFitter` is read only: its existing train/validation source split and Hybrid `betaTree` outputs provide the baseline.

## Current analysis

The primary outputs are under `results/comparison/plots/`:

- `residual_fits_Z1.pdf`, `residual_fits_Z2.pdf`, `residual_fits_Z6.pdf`, `residual_fits_Z8.pdf`: one page per 0.05-wide beta bin, overlaying Hybrid and New residual distributions plus Gaussian core fits.
- `core_sigma_vs_beta.pdf`: Hybrid and New Gaussian core width as line graphs.
- `core_mean_vs_beta.pdf`: Hybrid and New Gaussian core mean as line graphs.
- `relative_improvement_vs_beta.pdf`: `(1 - sigma_new/sigma_hybrid) * 100%` as a line graph.
- `results/comparison/per_bin_summary.csv`: numerical values used by the plots.

No total-beta spectrum or scatter plot is part of the main result.

See `docs/hybrid-comparison.md` for the current findings.

## Reproduce on ssh:sdiat

Extraction and ROOT writing use the existing system Python/uproot environment. The estimator training uses a temporary Python 3.6 environment because the host glibc/compiler cannot use current NumPy wheels; nothing is installed into the system Python.

```bash
python3 analysis/extract_inputs.py

python3 -m venv /tmp/betaEstimator-venv
/tmp/betaEstimator-venv/bin/pip install --only-binary=:all: -r requirements-sdiat.txt
/tmp/betaEstimator-venv/bin/python analysis/train_estimator.py

python3 analysis/write_comparison.py
root -l -b -q 'macro/plotHybridComparison.C()'
```

The canonical source split is the snapshot in `config/sources.conf`. Generated analysis products are ignored by Git.
