# Charged Polymer Simulation and Analysis

This project simulates charged polymer chains and provides tools for analyzing the resulting configurations and observables. The code is split between a C++ simulation engine and Python scripts for data analysis and plotting.

## Features

- **Polymer Simulation:** Uses Monte Carlo techniques to simulate the behavior of charged polymers.
- **Observable Measurements:** Computes quantities such as the structure factor, radius of gyration, and tangent-tangent correlation function.
- **Data Analysis & Plotting:** Python modules generate plots for polymer conformations, structure factors, persistence length fitting, and more.
- **Modular Analysis:** Separate scripts for different analyses (mechanical response, two-length scale analysis, machine learning based inversion, etc.).

## File Structure

- **/code/**
  Contains the C++ source code:
  - `charged_polymer.cpp` and `charged_polymer.h`: Implementation of the polymer simulation.
  - `main.cpp`: Main entry point for running simulations.

- **/plot/**
  Python scripts for generating figures:
  - `main_plot.py`: Central script to produce various plots (configuration, Sq, ML-based plots, etc.).
  - `conformation_plot.py`: Plots polymer conformations and related scales.
  - `ML_Sq_plot.py`: Generates plots related to the machine learning analysis of structure factor data.

- **/analyze/**
  Python modules for detailed data analysis:
  - `main_analyze.py`: Runs a full analysis of simulation data.
  - `plot_analyze.py`: Contains functions to plot observables and persistence lengths.
  - `ML_analyze.py`: Scripts for Gaussian Process prediction and feature analysis.

- **/data/** (assumed)
  Simulation output (CSV files with polymer configurations and observables).

## Requirements

- **C++ Compiler:** The simulation code requires a C++ compiler supporting C++17.
- **Python:** Python 3 with the following packages:
  - NumPy
  - Matplotlib
  - SciPy
  - Other standard libraries (os, time, etc.)


## Requirements

- **C++ Compiler:**
  A C++ compiler supporting C++17 is required to build and run the simulation.

- **Python:**
  Python 3 is required along with the following packages:
  - [NumPy](https://numpy.org/)
  - [Matplotlib](https://matplotlib.org/)
  - [SciPy](https://www.scipy.org/)
  - Other standard libraries (os, time, etc.)


## Compilation and Running the Simulation

### Build the Simulation

A simple way to build the C++ simulation is to use a command like:

```sh
g++ -std=c++17 code/main.cpp code/charged_polymer.cpp -o charged_polymer
```

Make sure that all required source files are included in the build.

### Run the Simulation

The simulation expects command-line arguments. For example:

```sh
./charged_polymer L kappa A invK save_more_config folder
```

Where:
- `L` is the polymer chain length,
- `kappa` is the bending stiffness,
- `A` is an energy amplitude parameter,
- `invK` is the inverse of some characteristic constant,
- `save_more_config` (0 or 1) indicates whether to save additional configurations,
- `folder` is the directory where output files (configuration and observables CSVs) will be saved.

Example:

```sh
./charged_polymer 100 30 5 3 1 ./output
```

## Data Analysis and Plotting

After running a simulation, you can analyze and visualize the results using the provided Python scripts.

### Plotting Figures

To generate the default set of plots, run:

```sh
python3 plot/main_plot.py
```

This script will call various plotting functions (e.g., `plot_configuration()`, `plot_two_length_scale()`, etc.) and save figures in the appropriate directory.

### Running Analysis

For a more detailed analysis of your simulation data, execute:

```sh
python3 analyze/main_analyze.py
```

This script processes observable data, runs persistence length fittings, and generates plots for mechanical response and tangent-correlation functions.

## Customization

- **Parameter Settings:** Many of the functions (both in C++ and Python) allow setting parameters such as polymer length, number of Monte Carlo sweeps, and plotting styles.
- **Plot Appearance:** The Python scripts use LaTeX formatting for text in plots. Customize these settings (e.g., line widths, fonts) in the respective modules if needed.

