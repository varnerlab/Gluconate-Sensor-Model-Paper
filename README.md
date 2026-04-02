# Resource Competition Modeling of a Cell-Free Gluconate Biosensor

This repository contains the model code, experimental data, parameter estimation scripts, and manuscript source for:

> **Resource Competition Modeling of a Cell-Free Gluconate Biosensor**
> A. Adhikari, A. Mukherjee, A. Nguyen, H. Lim, and J.D. Varner
> *Frontiers in Bioengineering and Biotechnology* (submitted 2026)

## Overview

We developed a mechanistic two-layer resource competition model for a D-gluconate biosensor circuit operating in the reconstituted cell-free system PURExpress. The model describes:

1. **Machinery allocation** --- competitive allocation of RNA polymerase and ribosome between genes using Michaelis-Menten-like expressions derived from elementary transcription/translation kinetics.
2. **Consumable resource depletion** --- irreversible exhaustion of nucleotides, energy cofactors, and amino acids through ODEs driven by cumulative biosynthetic flux.

The model predicted that 76% of transcriptional resources were consumed by 12 hours compared to 2% of translational resources, identifying transcription as the primary bottleneck. Without re-estimation, the model predicted the dose-response relationship of Venus protein over four orders of magnitude of D-gluconate concentration.

## Repository Structure

```
Gluconate-Sensor-Model-Paper/
|-- code/
|   |-- src/                    # Model source code (Julia)
|   |   |-- GluconateBiosensor.jl   # Module entry point
|   |   |-- Model.jl                # ODE system (balances!, transcription_control, etc.)
|   |   |-- Solve.jl                # ODE solver wrapper
|   |   |-- Types.jl                # Data structures
|   |   |-- Parameters.jl           # Parameter bounds, mapping, defaults
|   |   |-- Objective.jl            # 6 training objectives
|   |   |-- Biophysical.jl          # Load biophysical constants
|   |   |-- Data.jl                 # Load experimental data from CSV
|   |   |-- CellFree.json           # Biophysical constants (RNAP/ribosome conc, rates)
|   |-- scripts/                # Runnable scripts
|   |   |-- estimate_parameters.jl  # MO/SO cycling parameter estimation
|   |   |-- run_dynamics.jl         # Simulate time courses for ensemble
|   |   |-- run_dose_response.jl    # Simulate dose-response curves
|   |   |-- run_sensitivity.jl      # Morris global sensitivity analysis
|   |   |-- run_fedbatch_sweep.jl   # Synthetic fed-batch supplementation
|   |   |-- make_figures.jl         # Generate all publication figures
|   |   |-- analyze_ensemble.jl     # Parameter distributions and correlations
|   |-- data/                   # Experimental measurements
|   |   |-- mRNA_data.csv           # qPCR mRNA measurements (Venus, GntR)
|   |   |-- protein_data.csv        # Venus protein fluorescence time courses
|   |   |-- dose_response.csv       # Venus protein at 12h vs gluconate
|   |-- results/                # Generated outputs (not tracked; see below)
|   |-- Project.toml            # Julia package dependencies
|-- paper/
|   |-- Paper.tex               # Main manuscript (LaTeX)
|   |-- Supplement.tex          # Supplementary material (derivations)
|   |-- References.bib          # Bibliography
|   |-- Makefile                # Build: make all / make clean / make distclean
|   |-- figs/                   # Publication figures (PDF)
|-- LICENSE                     # MIT License
```

## Requirements

- **Julia** >= 1.10 (tested with 1.12.5)
- Required packages are specified in `code/Project.toml` and include:
  - `DifferentialEquations.jl` --- ODE solver
  - `ParetoEnsembles.jl` --- multi-objective optimization
  - `CairoMakie.jl` --- figure generation
  - `GlobalSensitivity.jl` --- Morris sensitivity analysis
  - `CSV.jl`, `DataFrames.jl`, `DelimitedFiles.jl`, `JSON.jl`, `Statistics.jl`

## Quick Start

### 1. Install dependencies

```bash
cd code
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### 2. Estimate parameters

This runs 15 cycles of multi-objective/single-objective estimation (~30-45 min with 8+ threads):

```bash
julia -t 8 --project=. scripts/estimate_parameters.jl
```

Output: `results/PC_final.dat` (25 x N parameter matrix), `results/EC_final.dat` (6 x N error matrix), and per-cycle archives `results/EC_cycle*.dat`, `results/PC_cycle*.dat`.

### 3. Simulate ensemble dynamics

```bash
julia -t 8 --project=. scripts/run_dynamics.jl
```

Output: ensemble trajectories for mRNA, protein, and resource fractions saved to `results/`.

### 4. Run dose-response prediction

```bash
julia -t 8 --project=. scripts/run_dose_response.jl
```

### 5. Run sensitivity analysis

```bash
julia --project=. scripts/run_sensitivity.jl
```

### 6. Run synthetic fed-batch analysis

```bash
julia -t 8 --project=. scripts/run_fedbatch_sweep.jl
```

### 7. Generate figures

After running steps 2-6, generate all publication figures:

```bash
julia --project=. scripts/make_figures.jl
```

Figures are saved to both `results/figures/` and `paper/figs/`.

### 8. Build the paper

```bash
cd paper
make          # runs pdflatex + bibtex
make clean    # remove auxiliary files
make distclean  # also remove PDF
```

## Experimental Data

All experimental measurements are in `code/data/` in CSV format:

| File | Description |
|------|-------------|
| `mRNA_data.csv` | Venus and GntR mRNA concentrations (qPCR) at 0, 2, 4, 6, 12 hours for three conditions: unrepressed (-GntR), de-repressed (+GntR +gluconate), repressed (+GntR -gluconate) |
| `protein_data.csv` | Venus protein fluorescence time courses (5-min intervals, 12 hours) at 10 gluconate concentrations (0-20 mM) plus no-GntR control |
| `dose_response.csv` | Venus protein endpoint (12h) vs gluconate concentration |

## Model Description

The model consists of 11 ODEs tracking 3 gene concentrations (constant), 3 mRNA concentrations, 3 protein concentrations, and 2 resource fractions. Key features:

- **Transcription control**: Boltzmann partition function formalism with pseudo-energy parameters for each promoter configuration (see Supplement S3)
- **Kinetic rates**: derived from 4-step McClure mechanism for transcription initiation (see Supplement S1)
- **Machinery allocation**: competitive allocation of free RNAP and ribosome (see Supplement S2)
- **Resource depletion**: flux-driven ODE model for consumable resources (see Supplement S4)
- **25 estimated parameters**: pseudo-energies, Hill coefficients, time constants, degradation modifiers, resource consumption rates

## Parameter Estimation

Parameters are estimated using ParetoEnsembles.jl with 6 training objectives:

1. Venus mRNA SSE at 10 mM gluconate
2. Venus protein SSE at 10 mM gluconate
3. GntR mRNA SSE at 10 mM gluconate
4. Venus protein SSE at 0 mM gluconate (repression floor)
5. GntR protein regularization
6. Venus protein SSE without GntR (expression ceiling)

The MO/SO cycling strategy alternates between multi-objective Pareto exploration and single-objective refinement of the worst-performing objective. The final ensemble is post-filtered for physiological consistency and prediction quality.

## Citation

If you use this code or data, please cite:

```bibtex
@article{adhikari2026gluconate,
  title={Resource Competition Modeling of a Cell-Free Gluconate Biosensor},
  author={Adhikari, Abhinav and Mukherjee, Aniruddha and Nguyen, Andrew and Lim, Ha Eun and Varner, Jeffrey D},
  journal={Frontiers in Bioengineering and Biotechnology},
  year={2026}
}
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
