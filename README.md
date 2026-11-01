# Velocity reconstruction — archived 2026-11-01

One closed research effort containing two code modules:

- `betaFitter/`: Linear, original Nonlinear and material-based Hybrid TOF propagation.
- `betaEstimator/`: flight-trained reconstruction, mass constraints, detector-response studies and final MC audits.

Velocity-reconstruction research has stopped. Existing algorithms, failed studies and final diagnostics are preserved as research records.

## Restore

Clone this repository on branch `archive/velocity-reconstruction`, or clone the enclosed `history.bundle` with that branch. The two module directories must remain siblings. Run commands from the relevant module directory; its README identifies the original entry points. The live sdiat paths remain `/nas01/home/hxwu/AMSSoft/betaFitter` and `/nas01/home/hxwu/AMSSoft/betaEstimator`.

The combined tar archive also includes `report/nonlinear_report.pptx`, selected scientific evidence/frozen model files, an inventory of data retained on sdiat and this complete Git history. Native AMS/ROOT dependencies and original ROOT/NPZ event datasets remain in their original storage. The report covers the nonlinear/Hybrid propagation work; later flight-trained and terminal-audit evidence is retained with the estimator code.

## Historical evidence

- Hybrid physical-region comparison: `betaFitter/research/` and its configuration snapshots.
- Frozen flight validation: `betaEstimator/results/tofcalib_validation/`, 907,830 independent RICH references; overall beta-residual W68 reduction 2.9%, with 3.0% degradation at 2.5–6 GV.
- Stopping evidence: `betaEstimator/results/terminal_ceiling_20261015/` and `results/mass_prior_snap_audit_20261015/`.

The latter audits distinguish beta error, mass-peak compression and measurable PID performance. Their numerical comparisons apply to the recorded samples and inputs.
