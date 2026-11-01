# Nonlinear particle velocity reconstruction for AMS TOF

Reconstruct particle velocity at the AMS center or before TOF from four timing measurements. The original nonlinear model converts measured TOF and Tracker deposits into propagation losses. Hybrid combines a charge-aware material mean with a correction based on centered measured deposits.

## Build and run

On sdiat, from the project directory:

```bash
./build.sh
./run.sh input.root output.root 0 1.9 1.0 all center
```

`build.sh` and `run.sh` load the AMSVX/ROOT5 environment through `/home/ams/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/amsroot_TOFBetaFit.sh`. Use ROOT6 from a clean shell for merging and plotting. `./build.sh debug` and `./build.sh clean` select the other build modes.

External dependencies are AMSVX, the AMS data directory, and its `LAPP/dEdxPDF/g4mscan.root` material map. Inputs are `amstreea` ntuples. `input_Z*.list` and `config/*.sources` retain source paths; raw ROOT data and the AMS software are stored separately from the code archive.

```text
./run.sh input output [option=0] [tofScale=1] [trackerScale=1]
  [scaleMode=all] [reference=center] [scaleConfig=none]
  [hybridConfig=none] [materialStepCm=0]
```

- Options: `0` reconstruct, `1` fit a common loss scale at fixed MC beta, `2` benchmark, `3` save beta differences, `4` export Hybrid training samples, `5` train Hybrid parameters, `-1/-2` export field/loss information.
- `tofScale` and `trackerScale` control the original nonlinear model. Direct executable invocation defaults to TOF scale `2`; the wrapper defaults to `1`.
- `scaleMode`: `all` scales S1/S2/S3, `s1s2` scales S1/S2, `s2` scales S2. Other deposits retain coefficient 1.
- `reference`: `center` or `before_tof`.
- `scaleConfig`: three columns, `charge tofScale trackerScale`. Charges missing from the file are skipped.
- `hybridConfig`: calibration path, or `none` to disable Hybrid.
- `materialStepCm`: `0` for coarse propagation; positive values set the maximum z step. The formal fine configuration uses `1` cm.

The output `betaTree` contains `linearBeta`, Tracker-aware `nonlinearBeta`, old TOF-only `nonlinearBetaLegacy`, `hybridBeta`, Hybrid status/chi-square/step count, `mcBeta`, rigidity fields, and reconstructed/MC charges. Disabled Hybrid produces NaN beta and status -1.

## Physical-boundary Hybrid

```bash
./run.sh input.root beta.root 0 1 1 all center \
  config/energy_loss_scales_center_common_z1268.conf \
  config/hybrid_material_fine.conf 1
```

The regional model is:

```text
loss = materialScale * materialMean
     + fluctuationScale * (measuredDeposit - meanDeposit)
```

`materialMean` includes projectile-charge Z². For downward particles, the four calibration regions are S1 entrance to S2 exit; S2 exit to L5 at z=1.698 cm; L5 to S3 entrance; and S3 entrance to the post-S4 acceptance plane at z=-69.975 cm. TOF surfaces follow the event's hit bars. The last endpoint includes material downstream of S4.

The deposit variables are S1+S2; L2+L3+L4+half L5; half L5+L6+L7+L8; and S3+S4. Material integration follows native map cells and detector nodes. Each event fits inverse beta and one common time offset. The fine calibration has 128 rows: four charges, four regions, and eight local-beta bins.

```bash
./run.sh input.root training_part.root 4 1 1 all center none none 1
find hybrid_training_parts -name '*.root' | sort > training.list
./run.sh training.list config/hybrid_material_fine.conf 5
```

The training tree retains boundary audit fields and primary L5 truth. Historical comparison output `mcBeta` uses acceptance checkpoint 10. See [the full method and results](docs/hybrid-material-fine.md). The coarse configuration with the historical `20260904` filename has also been retrained with physical boundaries.

## Results and research

Current physical-boundary plots: `macro/results/beta_comparison/hybrid_material_physical*`. Validation data: `results/hybrid_physical_20260913/merged/`. Old `hybrid_material_fine*` plots describe the earlier fixed-boundary prototype.

Common validation counts for Z1/Z2/Z6/Z8 are 193,655 / 118,554 / 45,894 / 24,747. For residual `1/beta_reco - 1/beta_MC`, Hybrid improves low-beta core means for Z1/Z2 and produces small negative means for Z6/Z8. Core widths remain close to the original nonlinear model.

Use `macro/plotChargeBetaComparison.C` with the four charge-specific merged files under ROOT6. The complete command and numerical tables are in [the Chinese README](README.md) and [the Hybrid note](docs/hybrid-material-fine.md).

- `docs/energy-loss-scale-scope.md`: cumulative-loss explanation and paired scale-scope tests.
- `research/zeta_study_20260722/`: global scale fits and timing-kernel decomposition.
- `research/nonlinear_bias_decomposition_20260802/`: propagation and cluster-time bias decomposition.
- `research/real_data_validation_20260722/`: blind four-layer, disjoint-pair, RICH, and leave-one-out tests.

Research sources, notes, and compact numerical summaries are versioned. Generated ROOT files, plots, binaries, and batch logs remain outside the code snapshot. `ref/LxBetalhd.*` is an external implementation reference; production uses the AMSVX interfaces.

`run_local.sh` retains its existing eight-argument interface for Linear/Nonlinear. For Hybrid batches, call `run.sh` per file. The later independent estimator is maintained in the neighboring `betaEstimator` project.
