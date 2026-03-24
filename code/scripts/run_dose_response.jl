# run_dose_response.jl - Dose-response prediction (VALIDATION, not training)
#
# Usage: julia --project=code scripts/run_dose_response.jl

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
exp_data = load_experimental_data()

# Gluconate concentrations to sweep (same as experimental dose-response)
gluconate_doses = [0.0001, 0.001, 0.01, 0.1, 0.5, 1.0, 5.0, 10.0, 20.0]
n_doses = length(gluconate_doses)

# Storage: Venus protein at t=12h for each dose and parameter set
venus_12h = zeros(n_doses, n_ensemble)

println("Simulating dose-response ($(n_doses) doses × $(n_ensemble) parameter sets)...")
for (d, gluc) in enumerate(gluconate_doses)
    n_ok = 0
    for k in 1:n_ensemble
        pvec = PC[:, k]
        params = vector_to_parameters(pvec)
        sol = try
            simulate(params, bio, genes, gluc)
        catch
            continue
        end
        if sol.t[end] < 11.5
            continue
        end
        n_ok += 1
        venus_12h[d, k] = sol(12.0)[8]  # protein_Venus at t=12h
    end
    println("  [$(gluc) mM] $(n_ok)/$(n_ensemble) successful, Venus = $(round(mean(venus_12h[d,:]), sigdigits=4)) ± $(round(std(venus_12h[d,:]), sigdigits=4)) μM")
end

# --- Save results ---
writedlm(joinpath(RESULTS_DIR, "dose_gluconate.dat"), gluconate_doses)
writedlm(joinpath(RESULTS_DIR, "dose_venus_12h.dat"), venus_12h)

# --- Compare to experimental data ---
println("\n--- Dose-Response Comparison ---")
println("$(rpad("Gluconate (mM)", 18)) $(rpad("Exp Mean (μM)", 16)) $(rpad("Sim Mean (μM)", 16)) $(rpad("Sim Std", 12))")
for (i, gluc) in enumerate(exp_data.dose_gluconate)
    d_idx = findfirst(x -> isapprox(x, gluc, rtol=0.1), gluconate_doses)
    if d_idx !== nothing
        sim_mean = mean(venus_12h[d_idx, :])
        sim_std = std(venus_12h[d_idx, :])
        println("  $(rpad(gluc, 16)) $(rpad(round(exp_data.dose_venus_mean[i], sigdigits=4), 14)) $(rpad(round(sim_mean, sigdigits=4), 14)) $(round(sim_std, sigdigits=4))")
    end
end

println("\nDone. Dose-response results saved to $(RESULTS_DIR)")
