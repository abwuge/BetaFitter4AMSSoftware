<div align="center">
  <h1>Research on Non-linear Reconstruction Method for Particle Velocity in AMS Time-of-Flight Detector</h1>
  
  [![English](https://badgen.net/badge/Language/English/blue?icon=github)](README_EN.md) [![简体中文](https://badgen.net/badge/语言/简体中文/red?icon=github)](README.md)
</div>

## Project Introduction

This project aims to study the non-linear reconstruction method for particle velocity in the Time-of-Flight (TOF) detector of the Alpha Magnetic Spectrometer (AMS). AMS is a particle physics experiment operating on the International Space Station, and the TOF detector is used to measure the flight direction and velocity of charged particles. However, due to ionization energy loss of charged particles in TOF materials, particles decelerate. Using a linear function to fit the time-space relationship of particles leads to certain biases in velocity reconstruction, particularly noticeable at low energies.

This research introduces particle energy loss terms in the velocity fitting process, implementing non-linear fitting to reduce reconstruction bias and improve precision. This will help AMS more accurately identify nuclear isotopes.

## Features

- Particle motion simulation and energy loss calculation
- Linear and non-linear particle velocity reconstruction methods
- Magnetic field data visualization and analysis
- Energy loss correction models
- Beta (velocity/speed of light) comparison analysis tools

## Software Architecture

This project is developed based on the ROOT framework and includes the following components:

- `ParticleData`: Particle data structure definitions
- `ParticlePropagator`: Particle propagation algorithm implementation
- `BetaFitter`: Core class for particle velocity fitting
- `macro`: Data analysis and visualization scripts

## Installation and Usage

### System Requirements
- Linux operating system
- ROOT 5.34+ (ROOT 6.x supported)
- Compiler with C++11 standard support (gcc or icpx)

### Compilation Steps

1. Standard compilation:
```bash
./build.csh
```

2. Debug mode compilation:
```bash
./build.csh debug
```

3. Clean compilation files:
```bash
./build.csh clean
```

### Running Examples

```bash
./run.csh [input_file] [output_file] [parameters]
```

Or use the local processing script:
```bash
./run_local.sh
```

The final argument selects the beta reconstruction reference point: `center`
means the AMS center between S2 and S3 (default), while `before_tof` means the
point immediately before the particle enters TOF:

```bash
./run.sh input.root output.root 0 1.9 1.0 all center
./run.sh input.root output.root 0 1.9 1.0 all before_tof
```

The TOF and tracker energy-loss scales are independent. The tracker scale is
the fifth argument and defaults to `1.0`:

```bash
./run.sh input.root output.root 0 1.9 1.2 all center
```

### Joint TOF/tracker scale fit

Option `4` jointly fits one TOF zeta and one tracker energy-loss scale shared by all selected MC events. Each event's common time offset is profiled out analytically; the implementation does not fit per-event scales and then average them. The input may be one ROOT file or a `.list` file containing one ROOT path per line:

```bash
./run.sh input_Z2.list global_scale_Z2.root 4 0.9 all center 0 6 0 20
```

Here `0.9` is the MC beta upper limit, `0 6` is the zeta range, and `0 20` is the tracker-scale range. The two scales are minimized together on one event sample that remains valid across the full parameter box. The fit reports measured-time and checkpoint-truth results in `globalScaleTree`, plus 41-point zeta and tracker-scale chi-square slices through the joint optima.

This entry point promotes the global-profile method from [`research/zeta_study_20260722`](research/zeta_study_20260722/README.md) into the main program and extends it to the selectable AMS-center reference point.

## Data Analysis

The project provides various data analysis and visualization tools:

Plots are stored under `macro/results/<analysis-category>/` by default. For
example, beta comparisons go to `macro/results/beta_comparison/`, while energy
loss plots go to `macro/results/energy_loss/`. The `outputName` argument still
accepts an explicit output path when needed.

### Beta Comparison Analysis

```bash
root -l 'macro/plotBetaComparison.C("test.root")'
```

### Energy Loss Analysis

```bash
root -l 'macro/plotEnergyLoss.C("test_el.root")'
```

### Magnetic Field Analysis

```bash
root -l 'macro/plotMagneticField.C(0.0, 0.0, "test.root")'
```

## Results Analysis

The project demonstrates significant precision improvement with the non-linear method at lower β values through residual analysis comparing linear and non-linear velocity reconstruction methods. Research shows that energy loss of charged particles in detector materials produces systematic biases in velocity measurements, which can be significantly improved by introducing energy loss correction terms.
