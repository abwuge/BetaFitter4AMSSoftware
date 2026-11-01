# ζ≈2 study (2026-07-22)

Independent validation of the energy-loss scale used by `BetaNL`.
Research sources and compact numerical summaries accompany this historical study.

## Main results (all valid four-hit events, β_MC < 0.9)

| Species | checkpoint-truth global ζ | measured-time global ζ |
|---|---:|---:|
| Z=2 | 1.9631001 | 1.8927570 |
| Z=6 | 2.0689139 | 1.9556252 |
| Z=8 (O-16 mass 14.899169 GeV) | 2.0286581 | 1.9112594 |

The global fit minimizes the sum of event χ² after profiling out each event's common time offset. It does not average eventwise Minuit estimates.

## Z=2 first-order timing-kernel decomposition

For β_MC < 0.9:

- S1 active bracket: 0.293458
- passive gap S1→S2: 0.420786
- S2 active bracket: 0.307646
- passive S2→Tracker L2: 0.498616
- inner tracker/magnet: 0.404099
- S3 and downstream combined: 0.007381
- total: 1.931986

Thus S1+S2 active brackets contribute about 0.601, while early passive material contributes about 1.324. The exact nonlinear checkpoint fit gives 1.963, agreeing within about 1.6%.

## Reproduction

From the repository root with ROOT available:

```bash
g++ -O3 -std=c++17 research/zeta_study_20260722/src/global_zeta_profile.cpp \
  $(root-config --cflags --libs) -o /tmp/global_zeta_profile
/tmp/global_zeta_profile input_Z2.list 100

g++ -O3 -std=c++17 research/zeta_study_20260722/src/timing_kernel_decomposition.cpp \
  $(root-config --cflags --libs) -o /tmp/timing_kernel_decomposition
/tmp/timing_kernel_decomposition input_Z2.list 100
```

## Important limitations

- The old input ntuples do not contain the newer `mtofdep[4]` Geant truth-deposit branch, so reconstructed `GetEdep()` cannot yet be separated cleanly into primary, secondary, and detector-response components.
- Checkpoint truth time is reconstructed by interpolating momentum between acceptance planes and using reconstructed path lengths; it is not a direct Geant global-time branch.
- `tof_qs` is based on `IsGoodQGeom`; the source plotting macro labels `tof_qs/100==11` as unoverlapped, but it is more precisely a good-path geometry selection excluding problematic edge/overlap-like regions.
