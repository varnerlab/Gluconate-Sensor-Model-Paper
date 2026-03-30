# Gluconate Biosensor Model Code

Computational code for the paper:

**"Modeling and Analysis of a Cell-Free Gluconate Responsive Biosensor"**
Adhikari, Murti, Narayanan, Lim, and Varner

## Model

The model describes a D-gluconate responsive biosensor circuit operating in a reconstituted cell-free (PURExpress) system.
The circuit consists of two genes (P70-GntR and mP70-Venus) where the transcription factor GntR represses Venus expression,
and D-gluconate relieves this repression by binding to GntR.

The model captures:
- **Transcription and translation kinetics** using effective biophysical rate expressions
- **Transcriptional regulation** via a Boltzmann partition function formalism
- **Resource competition** through two mechanistic layers:
  - *Machinery allocation*: RNAP and ribosomes are shared across genes/mRNAs (algebraic competition)
  - *Consumable depletion*: NTP and amino acid pools deplete proportionally to actual TX/TL flux (ODE)

### State variables (11)
| Index | Species | Units |
|-------|---------|-------|
| 1-3 | Gene concentrations (GntR, Venus, σ70) | μM |
| 4-6 | mRNA (GntR, Venus, σ70) | μM |
| 7-9 | Protein (GntR, Venus, σ70) | μM |
| 10 | εX — TX consumable resource fraction | dimensionless |
| 11 | εL — TL consumable resource fraction | dimensionless |

### Parameters
25 estimated parameters (see `src/Parameters.jl` for names, bounds, and defaults).

## Requirements

- Julia ≥ 1.10
- 5+ threads recommended for parallel parameter estimation

## Reproducing Results

```bash
# 1. Install dependencies
cd code
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# 2. Estimate parameters (parallel, ~hours)
julia -t 5 --project=. scripts/estimate_parameters.jl

# 3. Simulate dynamics at training condition
julia --project=. scripts/run_dynamics.jl

# 4. Predict dose-response (validation)
julia --project=. scripts/run_dose_response.jl

# 5. Global sensitivity analysis
julia --project=. scripts/run_sensitivity.jl

# 6. Generate figures
julia --project=. scripts/make_figures.jl
```

## Directory Structure

```
code/
├── Project.toml              # Julia dependencies
├── data/                     # Experimental data (CSV)
│   ├── protein_data.csv      # Venus protein time courses
│   ├── mRNA_data.csv         # Venus + GntR mRNA measurements
│   └── dose_response.csv     # Dose-response at 12h
├── src/                      # Model source code
│   ├── GluconateBiosensor.jl # Main module
│   ├── Types.jl              # Data structures
│   ├── Biophysical.jl        # Biophysical constants
│   ├── CellFree.json         # Cell-free system parameters
│   ├── Parameters.jl         # Parameter mapping and bounds
│   ├── Model.jl              # ODE system (balances, kinetics, control)
│   ├── Solve.jl              # ODE solver wrapper
│   ├── Data.jl               # Experimental data loading
│   └── Objective.jl          # Objective functions for estimation
├── scripts/                  # Analysis scripts
│   ├── estimate_parameters.jl
│   ├── run_dynamics.jl
│   ├── run_dose_response.jl
│   ├── run_sensitivity.jl
│   └── make_figures.jl
└── results/                  # Output (generated, not tracked in git)
```

## License

MIT
