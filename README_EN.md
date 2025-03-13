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

## Data Analysis

The project provides various data analysis and visualization tools:

### Beta Comparison Analysis

```bash
root -l 'macro/plotBetaComparison.C("test.root", "beta_comparison.pdf")'
```

### Energy Loss Analysis

```bash
root -l 'macro/plotEnergyLoss.C("test_el.root", "energy_loss.pdf")'
```

### Magnetic Field Analysis

```bash
root -l 'macro/plotMagneticField.C(0.0, 0.0, "test.root", "magnetic_field.pdf")'
```

## Results Analysis

The project demonstrates significant precision improvement with the non-linear method at lower β values through residual analysis comparing linear and non-linear velocity reconstruction methods. Research shows that energy loss of charged particles in detector materials produces systematic biases in velocity measurements, which can be significantly improved by introducing energy loss correction terms.
