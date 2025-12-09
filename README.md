
# 2D Schrödinger Equation & Slit Experiments — FYS3150 Project 5

## 📖 Overview

This repository contains the implementation for **Project 5** in *FYS3150 — Computational Physics*.

We solve the **two-dimensional time-dependent Schrödinger equation** on a unit-square box using a **Crank–Nicolson** scheme in space and time. The code simulates the evolution of a **Gaussian wave packet** interacting with slit potentials, allowing us to study:

- Free propagation in a box
- Single-slit, double-slit, and multi-slit barrier setups
- Diffraction and interference patterns
- Numerical conservation of total probability

Simulation data are written to file and can be analysed and visualised using the Python scripts in `scripts/`.

https://github.uio.no/user-attachments/assets/2c172ae4-0f34-4fea-994d-63e4b7717327

---

## 🎯 Project Goals

- Implement a Crank–Nicolson discretisation of the 2D time-dependent Schrödinger equation.
- Work with **complex-valued wave functions** and sparse matrices.
- Verify **probability conservation** and numerical stability.
- Construct time-independent slit potentials (single, double, triple).
- Study:
  - Free wave packet propagation
  - Transmission and reflection at barriers
  - Interference patterns for different slit configurations
- Provide a fully reproducible workflow (C++ code + analysis scripts).

---

## 📂 Repository Structure

```text
├── Makefile                           # Build File
├── README.md                          # Project Documentation (this file)
├── build/                             # Compiled executables
├── include/                           # C++ headers (interfaces)
│   └── utils.hpp                      # Utility function declerations
├── output/                            # Generated Results and Plots go here
├── scripts/                           # Scripts for Animation and Plotting
│   ├── animate_double_slit.py         # Creates animations of the slit experiment time evolution
│   ├── plot_double_slit.py            # Generates wavefunction plots for the (single/double/triple) slit simulation.
│   ├── plot_detection_probability.py  # Computes and plots detection probability at a screen
│   └── plot_probability_deviation.py  # Plots Probability Deviation over time
└── src/                               # Source code for the simulation
    ├── main.cpp                       # Main program: argument parsing, setup, time stepping and data generation
    └── utils.cpp                      # Utility function implementations
```

---

## ⚙️ Compilation & Dependencies

### Requirements

- C++17 compiler
- GNU Make
- Armadillo
- Python 3 (`numpy`, `matplotlib`)
- ffmpeg for animations

Build:

```bash
make main
```

Executable:

```
build/main
```

---

## ▶️ Running Simulations

Usage:

```
build/main <file_name_prefix> <h> <dt> <T>   <x_c> <sigma_x> <p_x>   <y_c> <sigma_y> <p_y>   <v0> <n_slits> [options]
```

### Arguments

| Name | Description |
|------|-------------|
| `file_name_prefix` | Prefix for output files (base name for all data written to `output/`). |
| `h` | Spatial step size (domain is `[0,1] × [0,1]`). |
| `dt` | Time step. |
| `T` | Total simulation time. |
| `x_c` | Initial wave packet center in x-direction. |
| `sigma_x` | Wave packet width in x-direction. |
| `p_x` | Initial momentum in x-direction. |
| `y_c` | Initial wave packet center in y-direction. |
| `sigma_y` | Wave packet width in y-direction. |
| `p_y` | Initial momentum in y-direction. |
| `v0` | Barrier height for the slit potential. |
| `n_slits` | Number of slits (`1` = single, `2` = double, `3` = triple). |

### Options

| Flag | Description |
|------|-------------|
| `--M <int>` | Override grid size `M` (grid points per spatial dimension). |

### Example

```
# Case 1: Free propagation (barrier off)
build/main free_propagation_barrier_off 0.005 2.5e-5 0.008 0.25 0.05 200.0 0.5 0.05 0.0 0.0 2

# Case 2: Double-slit (barrier on)
build/main double_slit_barrier_on 0.005 2.5e-5 0.008 0.25 0.05 200.0 0.5 0.10 0.0 1.0e10 2

# Single Slit Experiment
build/main single_slit_experiment 0.005 2.5e-5 0.002 0.25 0.05 200 0.5 0.20 0.0 1e10 1

# Double Slit Experiment
build/main double_slit_experiment 0.005 2.5e-5 0.002 0.25 0.05 200 0.5 0.20 0.0 1e10 2

# Tripple Slit Experiment
build/main triple_slit_experiment 0.005 2.5e-5 0.002 0.25 0.05 200 0.5 0.20 0.0 1e10 3
```

---

## 📊 Plotting

Example:

```bash
python scripts/plot_probability_deviation.py
python scripts/plot_double_slit.py
python scripts/plot_detection_probability.py
python scripts/animate_double_slit.py
```

---

## 👥 Authors

- Solveig Juntunen
- Tom Jussiaume
- Erling Sandbekk

---

## References
- [Project description and background](https://anderkve.github.io/FYS3150/book/projects/project5.html)
- Course: *FYS3150 – Computational Physics* (University of Oslo)

---

## TODO

### Refactoring / OOP
- Encapsulate solver logic into classes
- Improve naming conventions and file organization

### Python Script Input Handling
- Auto-detect available simulation files in `output/`
- Add command-line arguments (e.g., `--prefix`, `--times`, `--save-dir`)
- Improve error messages when files or folders are missing

### Performance Improvements
- Precompute LU factorization of A  
  *(Armadillo >= 12.2 provides `spsolve_factoriser` — integrate for major speedup)*

### Testing and Validation
- Add unit tests for:
  - Wave packet initialization
  - Potential generation (single/double/multi-slit)
  - CN update step and boundary conditions
