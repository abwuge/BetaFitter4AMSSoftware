# B1308 nonlinear 1/beta bias decomposition (2026-08-02)

This study diagnoses why the center-reference nonlinear reconstruction in
`results/nonlinear_B1308_20260802` is not flat and improves little over the
linear result. It does **not** modify production reconstruction code.

## Data and method

The first 20 B1308 files of He4, C12, and O16 were processed event-by-event.
All comparisons use the same events and the same reconstructed TOF path
lengths.

The diagnostic program compares:

- production `linearBeta`, `nonlinearBetaLegacy`, and `nonlinearBeta`;
- a robust one-dimensional reimplementation of the current center/tracker
  model, with the common time offset profiled analytically;
- checkpoint-derived truth flight times constructed from `mevmom1/mevcoo1`;
- an oracle propagation model that preserves each event's checkpoint energy
  change profile while varying only the center beta hypothesis;
- standard reconstructed `tof_tl`, raw two-side `tof_rtime`, and
  `tof_betahmc = BetaHR::GetMCBeta(1)`;
- tracker active deposits versus checkpoint total energy changes;
- a blind bar-pair test: files 0--9 determine the high-beta S2/S3 bar-pair
  baseline, and files 10--19 test the low-beta residual.

The checkpoint time is an interpolation diagnostic, not a direct Geant global
clock. Its use is independently checked by `tof_betahmc`, which is formed from
the earliest TofMCCluster times plus detector smearing.

## Main result: two large opposite biases cancel

The small production residual is not evidence that the nonlinear model is
nearly correct. It is the partial cancellation of two larger effects.

### 1. Propagation model underestimates the S2--S3 flight time

At representative low-beta bins, the current propagation model evaluated at
the true center beta differs from checkpoint truth almost entirely in the long
S2--S3 interval:

| species | beta | S1--S2 | S2--S3 | S3--S4 |
|---|---:|---:|---:|---:|
| He4 | 0.375 | +5.0 ps | -79.9 ps | -5.9 ps |
| C12 | 0.525 | +3.7 ps | -39.2 ps | -3.6 ps |
| O16 | 0.525 | +4.1 ps | -56.0 ps | -4.7 ps |

The sign means that the current model predicts too short a flight time. When
it is fitted to checkpoint truth times, the resulting median 1/beta biases are
+0.0180, +0.0090, and +0.0132 respectively.

The direct energy accounting explains the scale of this error. The current
center model uses unscaled tracker active deposits for the central region, but
the checkpoint total energy changes are much larger:

| species | beta | upper truth / active | lower truth / active |
|---|---:|---:|---:|
| He4 | 0.375 | 60.5 / 8.9 MeV = 6.8 | 60.6 / 7.0 MeV = 8.8 |
| C12 | 0.525 | 326.9 / 38.7 MeV = 8.6 | 319.1 / 30.0 MeV = 11.2 |
| O16 | 0.525 | 561.7 / 68.2 MeV = 8.5 | 592.5 / 53.6 MeV = 11.6 |

A single multiplier on `tk_edep` reduces the error but does not flatten every
species and beta bin, because the missing loss is not distributed at the
active silicon locations alone.

Changing the TOF zeta from 1.9 to about 4.3 has little leverage on this error:
in the center implementation zeta acts mainly on the short outer intervals,
not on the long central S2--S3 region.

### 2. Standard TOF cluster times make S2--S3 too short at low beta

Using the oracle propagation profile with standard `tof_tl` produces the
opposite low-beta 1/beta bias:

| species | beta | oracle model fitted to `tof_tl` |
|---|---:|---:|
| He4 | 0.375 | -0.0178 |
| C12 | 0.525 | -0.0111 |
| O16 | 0.525 | -0.0166 |

At the same bins, measured S2--S3 times are shorter than checkpoint truth by
about -72.7, -42.1, and -64.5 ps.

This is not caused by the optional charge-specific `DoMCtune`: the B1308
binary was built with `Makefile534_1SL6`, which does not define
`_DOFINALTUNE_`, so `tof_tl` is the ordinary `TofClusterHR::Time`.

The standard hit-time reconstruction is

```text
tm = -0.5 * (sdtm0 + sdtm1)
     - slew0 / A0^powindex
     - slew1 / A1^powindex
     - delay
```

using `TofAlign2` parameters. The raw two-side S2--S3 interval is typically
wrong by 0.1--0.3 ns; the slew/delay calibration removes most of it, but leaves
a smooth low-beta residual of tens of ps.

`tof_betahmc` is close to flat (typically within about 0.002 in the relevant
C12/O16 bins), which independently supports the checkpoint reference and
locates the negative component in the standard cluster-time reconstruction.

### 3. The remaining timing trend is not a fixed bar delay

A blind control learned one constant S2/S3 bar-pair baseline from high-beta
(0.85--0.95) events in files 0--9 and applied it to files 10--19.

| species | low-beta range | raw median | after blind bar-pair baseline | negative populated pairs |
|---|---|---:|---:|---:|
| He4 | 0.35--0.50 | -31.7 ps | -38.1 ps | 54 / 59 |
| C12 | 0.48--0.63 | -32.8 ps | -37.1 ps | 33 / 39 |
| O16 | 0.48--0.63 | -37.1 ps | -39.1 ps | 4 / 6 |

After the bar-pair baseline is frozen, the residual still becomes more
negative toward lower beta and retains correlations with S2/S3 ADC or
energy-deposit asymmetry. Therefore the trend is not an artifact of mixing
bars with different constant delays; it is consistent with a residual
high-ionization/time-walk response outside the calibration model's effective
range. Exact refitting must remain bar- and side-specific because the current
slew coefficients differ by counter.

## Excluded causes

- **Minuit failure is not the trend source.** A global-scan one-dimensional
  fit with the time offset profiled analytically differs from production
  `nonlinearBeta` by only about 1e-7--1e-5 in median, depending on bin.
- **Using MC `mtofdep` instead of reconstructed TOF edep changes little.**
- **There is no established checkpoint-10 versus z=0 reference shift.**
  Checkpoint 10 is the post-L5 state and its momentum persists through the
  material-free central gap.  The earlier linear interpolation of momentum
  between checkpoints 10 and 11 was invalid.

## Consequence for optimization

Fixing only one side is unsafe:

- correcting propagation while retaining the current `tof_tl` exposes the
  negative timing bias;
- correcting timing while retaining the current central energy-loss model
  exposes the positive propagation bias.

The next model should therefore be tested as four frozen variants on held-out
files:

1. current timing + current propagation (baseline);
2. corrected timing + current propagation;
3. current timing + corrected propagation;
4. corrected timing + corrected propagation.

For propagation, fit a segment-level average material-loss model from
checkpoint energy changes using inference-time observables such as charge,
seed beta/rigidity, direction/path, and deposits. Do not use event MC truth at
inference, and do not represent all missing material by one `tk_edep` scale.

For timing, fit the residual after the existing `TofAlign2` correction with a
train/validation split, retaining layer, bar, side ADC, charge, and amplitude
information. A constant high-beta bar baseline must be frozen separately from
the low-beta amplitude-dependent component.

The target metric remains event-level `nonlinearBeta` versus center `mcBeta`:
per-beta-bin core bias, width, tails, finite coverage, and the corresponding
linear comparison. Zeta fits are diagnostics only.

## Reproduction

Use the local ROOT 6 environment:

```bash
export PATH=/home/ams/hxwu/.local/bin:$PATH
g++ -O3 -std=c++17 src/bias_decomposition.cpp \
  $(root-config --cflags --libs) -o /tmp/bias_decomposition

/tmp/bias_decomposition \
  ../../results/nonlinear_B1308_20260802/input_C12.list \
  ../../results/nonlinear_B1308_20260802/event_beta/C12 \
  6 20 \
  results/C12_residual_components.csv \
  results/C12_timing_components.csv \
  results/C12_run_summary.txt \
  results/C12_energy_loss_components.csv
```

The stored result CSVs contain the complete 20-file binned summaries and are versioned together with the research sources. Event-level plotting inputs remain generated artifacts.

## ROOT diagnostic plots (2026-08-03)

The event-level plotting inputs are regenerated from the same first 20 B1308
files used by this decomposition and stored in `results/plot_events/`.
Run the ROOT macro with:

```bash
cd ~/AMSSoft/betaFitter/research/nonlinear_bias_decomposition_20260802
root -l -b -q 'src/plot_low_beta_diagnostics.C+()'
```

The ROOT file, PNGs, and PDFs are written to `plots/root_diagnostics/` and
copied to the main result tree at
`results/nonlinear_B1308_20260802/root_diagnostics/`.

For each species, `tof_timewalk_vs_beta_*.pdf` shows the raw S2--S3 timing
residual, the applied standard slew/delay correction, and the remaining
`tof_tl` residual relative to checkpoint truth. The points are medians and the
vertical bars span the 16--84% event interval.

`central_energy_loss_vs_beta_*.pdf` compares checkpoint total energy change
between S2 and S3 with the tracker active-silicon `tk_edep`, then shows the
missing energy change, their ratio (total and each half), and the event-level
2D comparison. All ROOT graphs, histograms, and canvases are retained in
`low_beta_diagnostics.root`.
