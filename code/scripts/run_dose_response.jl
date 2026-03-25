# run_dose_response.jl - Dense dose-response prediction (VALIDATION)
#
# Simulates at 150 gluconate concentrations across the full ensemble.
# Saves raw data so make_figures.jl can plot without re-simulating.
#
# Usage: julia -t 15 --project=code code/scripts/run_dose_response.jl

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

# Dense gluconate sweep: 150 points on log10 scale
log_doses = collect(range(-4.2, 1.5, length = 150))
dense_doses = 10.0 .^ log_doses
n_doses = length(dense_doses)

# --- Simulate (parallel over parameter sets) ---
println("Simulating dense dose-response ($(n_doses) doses × $(n_ensemble) param sets) on $(Threads.nthreads()) threads...")
venus_dense = zeros(n_doses, n_ensemble)

Threads.@threads for k in 1:n_ensemble
    params = vector_to_parameters(PC[:, k])
    for (d, gluc) in enumerate(dense_doses)
        sol = try
            simulate(params, bio, genes, gluc)
        catch
            continue
        end
        if sol.t[end] >= 11.5
            venus_dense[d, k] = sol(12.0)[8]
        end
    end
end
println("Done simulating.")

# --- Compute statistics ---
med = vec(mapslices(median, venus_dense, dims = 2))
q025 = vec(mapslices(x -> quantile(x, 0.025), venus_dense, dims = 2))
q975 = vec(mapslices(x -> quantile(x, 0.975), venus_dense, dims = 2))
μ = vec(mean(venus_dense, dims = 2))
σ = vec(std(venus_dense, dims = 2))

# --- Save everything ---
writedlm(joinpath(RESULTS_DIR, "dose_dense_log_gluconate.dat"), log_doses)
writedlm(joinpath(RESULTS_DIR, "dose_dense_venus_median.dat"), med)
writedlm(joinpath(RESULTS_DIR, "dose_dense_venus_q025.dat"), q025)
writedlm(joinpath(RESULTS_DIR, "dose_dense_venus_q975.dat"), q975)
writedlm(joinpath(RESULTS_DIR, "dose_dense_venus_mean.dat"), μ)
writedlm(joinpath(RESULTS_DIR, "dose_dense_venus_std.dat"), σ)

# --- Print comparison ---
println("\n--- Dose-Response Comparison ---")
exp_doses = exp_data.dose_gluconate
for (i, gluc) in enumerate(exp_doses)
    d_idx = argmin(abs.(dense_doses .- gluc))
    println("  $(gluc) mM: median=$(round(med[d_idx], sigdigits=4)), 95%CI=[$(round(q025[d_idx], sigdigits=3)), $(round(q975[d_idx], sigdigits=3))], exp=$(round(exp_data.dose_venus_mean[i], sigdigits=4))")
end

println("\nDone. Data saved to $(RESULTS_DIR)")
