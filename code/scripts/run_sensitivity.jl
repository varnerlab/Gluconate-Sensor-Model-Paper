# run_sensitivity.jl - Morris global sensitivity analysis
#
# Usage: julia --project=code scripts/run_sensitivity.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "GluconateBiosensor.jl"))
using .GluconateBiosensor
using GlobalSensitivity
using DelimitedFiles
using Statistics
using Random

const RESULTS_DIR = joinpath(@__DIR__, "..", "results")

# --- Setup ---
bio = load_biophysical_constants(joinpath(@__DIR__, "..", "src", "CellFree.json"))
genes = build_gene_info()
gluconate_conc = 10.0  # training condition

# --- Model performance function ---
# Returns a vector of 4 outputs: [Venus_mRNA_12h, Venus_protein_12h, GntR_mRNA_6h, GntR_protein_12h]
function model_performance(pvec::Vector{Float64})
    params = vector_to_parameters(pvec)
    sol = try
        simulate(params, bio, genes, gluconate_conc)
    catch
        return [0.0, 0.0, 0.0, 0.0]
    end

    if sol.t[end] < 11.5
        return [0.0, 0.0, 0.0, 0.0]
    end

    return [
        sol(6.0)[5],   # Venus mRNA at 6h (near peak)
        sol(12.0)[8],  # Venus protein at 12h
        sol(6.0)[4],   # GntR mRNA at 6h (near peak)
        sol(12.0)[7],  # GntR protein at 12h
    ]
end

# --- Morris sensitivity analysis ---
println("Running Morris sensitivity analysis ($(N_PARAMETERS) parameters, 4 outputs)...")
println("Parameter bounds:")
for i in 1:N_PARAMETERS
    println("  $(rpad(PARAMETER_NAMES[i], 22)) [$(PARAMETER_LOWER[i]), $(PARAMETER_UPPER[i])]")
end

# Use ensemble-derived bounds for Morris sampling (mean ± 3*std, clamped to original bounds)
# This ensures Morris samples from the region where the model is well-behaved
PC = readdlm(joinpath(RESULTS_DIR, "PC_final.dat"))
ensemble_bounds = Vector{Vector{Float64}}()
for i in 1:N_PARAMETERS
    vals = PC[i, :]
    μ = mean(vals)
    σ = std(vals)
    lo = max(PARAMETER_LOWER[i], μ - 3.0 * σ)
    hi = min(PARAMETER_UPPER[i], μ + 3.0 * σ)
    if lo >= hi
        lo = PARAMETER_LOWER[i]
        hi = PARAMETER_UPPER[i]
    end
    push!(ensemble_bounds, [lo, hi])
end

# Morris method settings
n_trajectories = 200  # number of Morris trajectories

# GlobalSensitivity.jl Morris calls f(x::Vector) -> Vector per sample
morris_result = gsa(model_performance, Morris(num_trajectory = n_trajectories),
                    ensemble_bounds)

# --- Extract results ---
output_names = ["Venus_mRNA_6h", "Venus_protein_12h", "GntR_mRNA_6h", "GntR_protein_12h"]

# Morris means (μ*) and variances (σ)
# morris_result.means is n_params × n_outputs
# morris_result.variances is n_params × n_outputs
means = morris_result.means      # absolute mean elementary effects
variances = morris_result.variances  # variance of elementary effects

# --- Save results ---
writedlm(joinpath(RESULTS_DIR, "sensitivity_means.dat"), means)
writedlm(joinpath(RESULTS_DIR, "sensitivity_variances.dat"), variances)

# --- Print summary ---
println("\n--- Morris Sensitivity Results (|μ*|) ---")
for (j, oname) in enumerate(output_names)
    println("\n$(oname):")
    # Sort parameters by influence
    sorted_idx = sortperm(abs.(means[:, j]), rev = true)
    for rank in 1:min(10, N_PARAMETERS)
        i = sorted_idx[rank]
        println("  $(rank). $(rpad(PARAMETER_NAMES[i], 22)) μ*=$(round(means[i,j], sigdigits=4))  σ=$(round(sqrt(abs(variances[i,j])), sigdigits=4))")
    end
end

println("\nDone. Sensitivity results saved to $(RESULTS_DIR)")
