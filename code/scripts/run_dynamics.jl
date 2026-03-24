# run_dynamics.jl - Simulate time courses for the parameter ensemble
#
# Usage: julia --project=code scripts/run_dynamics.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "GluconateBiosensor.jl"))
using .GluconateBiosensor
using DelimitedFiles
using Statistics

const RESULTS_DIR = joinpath(@__DIR__, "..", "results")

# --- Load ensemble ---
PC = readdlm(joinpath(RESULTS_DIR, "PC_final.dat"))
n_ensemble = size(PC, 2)
println("Loaded ensemble: $(n_ensemble) parameter sets")

# --- Setup ---
bio = load_biophysical_constants(joinpath(@__DIR__, "..", "src", "CellFree.json"))
genes = build_gene_info()

# --- Simulate at training condition (10 mM gluconate) ---
gluconate_conc = 10.0
dt = 0.1  # output time step (hr)
t_eval = collect(0.0:dt:12.0)
n_times = length(t_eval)

# Storage for ensemble trajectories
mRNA_GntR_ensemble = zeros(n_times, n_ensemble)
mRNA_Venus_ensemble = zeros(n_times, n_ensemble)
protein_GntR_ensemble = zeros(n_times, n_ensemble)
protein_Venus_ensemble = zeros(n_times, n_ensemble)
epsilon_X_ensemble = zeros(n_times, n_ensemble)
epsilon_L_ensemble = zeros(n_times, n_ensemble)

println("Simulating $(n_ensemble) parameter sets at $(gluconate_conc) mM gluconate...")
n_success = 0
for k in 1:n_ensemble
    pvec = PC[:, k]
    params = vector_to_parameters(pvec)
    sol = try
        simulate(params, bio, genes, gluconate_conc)
    catch
        continue
    end

    if sol.t[end] < 11.5
        continue
    end

    global n_success += 1
    for (i, t) in enumerate(t_eval)
        state = sol(t)
        mRNA_GntR_ensemble[i, n_success] = state[4]
        mRNA_Venus_ensemble[i, n_success] = state[5]
        protein_GntR_ensemble[i, n_success] = state[7]
        protein_Venus_ensemble[i, n_success] = state[8]
        epsilon_X_ensemble[i, n_success] = state[10]
        epsilon_L_ensemble[i, n_success] = state[11]
    end
end

println("Successful simulations: $(n_success)/$(n_ensemble)")

# Trim to successful simulations
mRNA_GntR_ensemble = mRNA_GntR_ensemble[:, 1:n_success]
mRNA_Venus_ensemble = mRNA_Venus_ensemble[:, 1:n_success]
protein_GntR_ensemble = protein_GntR_ensemble[:, 1:n_success]
protein_Venus_ensemble = protein_Venus_ensemble[:, 1:n_success]
epsilon_X_ensemble = epsilon_X_ensemble[:, 1:n_success]
epsilon_L_ensemble = epsilon_L_ensemble[:, 1:n_success]

# --- Save results ---
writedlm(joinpath(RESULTS_DIR, "time_eval.dat"), t_eval)
writedlm(joinpath(RESULTS_DIR, "mRNA_GntR_ensemble.dat"), mRNA_GntR_ensemble)
writedlm(joinpath(RESULTS_DIR, "mRNA_Venus_ensemble.dat"), mRNA_Venus_ensemble)
writedlm(joinpath(RESULTS_DIR, "protein_GntR_ensemble.dat"), protein_GntR_ensemble)
writedlm(joinpath(RESULTS_DIR, "protein_Venus_ensemble.dat"), protein_Venus_ensemble)
writedlm(joinpath(RESULTS_DIR, "epsilon_X_ensemble.dat"), epsilon_X_ensemble)
writedlm(joinpath(RESULTS_DIR, "epsilon_L_ensemble.dat"), epsilon_L_ensemble)

# --- Print summary statistics ---
println("\nAt t = 12 hr:")
println("  Venus mRNA:    mean=$(round(mean(mRNA_Venus_ensemble[end,:]), sigdigits=4)) ± $(round(std(mRNA_Venus_ensemble[end,:]), sigdigits=4))")
println("  Venus protein: mean=$(round(mean(protein_Venus_ensemble[end,:]), sigdigits=4)) ± $(round(std(protein_Venus_ensemble[end,:]), sigdigits=4))")
println("  GntR mRNA:     mean=$(round(mean(mRNA_GntR_ensemble[end,:]), sigdigits=4)) ± $(round(std(mRNA_GntR_ensemble[end,:]), sigdigits=4))")
println("  GntR protein:  mean=$(round(mean(protein_GntR_ensemble[end,:]), sigdigits=4)) ± $(round(std(protein_GntR_ensemble[end,:]), sigdigits=4))")
println("  ε_X:           mean=$(round(mean(epsilon_X_ensemble[end,:]), sigdigits=4))")
println("  ε_L:           mean=$(round(mean(epsilon_L_ensemble[end,:]), sigdigits=4))")

println("\nDone. Trajectories saved to $(RESULTS_DIR)")
