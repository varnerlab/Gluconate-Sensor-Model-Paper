# GluconateBiosensor.jl - Main module for the gluconate biosensor model
#
# Adhikari, Murti, Narayanan, Lim, Varner (2023)
# "Modeling and Analysis of a Cell-Free Gluconate Responsive Biosensor"

module GluconateBiosensor

using DifferentialEquations
using CSV, DataFrames, JSON
using Interpolations, LinearAlgebra, Statistics

# Source files
include("Types.jl")
include("Biophysical.jl")
include("Parameters.jl")
include("Model.jl")
include("Solve.jl")
include("Data.jl")
include("Objective.jl")

# Export types
export BiophysicalConstants, GeneInfo, ModelParameters, ExperimentalData

# Export model functions
export balances!, transcription_control, machinery_allocation, kinetic_rates
export simulate

# Export parameter utilities
export vector_to_parameters, parameters_to_vector, clamp_to_bounds, default_initial_guess
export PARAMETER_NAMES, PARAMETER_LOWER, PARAMETER_UPPER, PARAMETER_BOUNDS, N_PARAMETERS

# Export data and objective functions
export load_biophysical_constants, build_gene_info, load_experimental_data
export evaluate_objectives, evaluate_single_objective

# Export constants
export N_STATES, N_GENES, IDX_GENE, IDX_MRNA, IDX_PROTEIN, IDX_EPSILON_X, IDX_EPSILON_L
export I_GNTR, I_VENUS, I_SIGMA70

end # module
