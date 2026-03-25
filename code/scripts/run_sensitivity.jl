# run_sensitivity.jl - Morris global sensitivity analysis
#
# Usage: julia -t 15 --project=code code/scripts/run_sensitivity.jl

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
gluconate_conc = 10.0

# --- Model performance function ---
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
        sol(6.0)[5],   # Venus mRNA at 6h
        sol(12.0)[8],  # Venus protein at 12h
        sol(6.0)[4],   # GntR mRNA at 6h
        sol(12.0)[7],  # GntR protein at 12h
    ]
end

# --- Ensemble-derived bounds ---
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

# --- Morris sensitivity analysis ---
n_trajectories = 200
println("Running Morris sensitivity analysis ($(N_PARAMETERS) params, 4 outputs, $(n_trajectories) trajectories)...")

morris_result = gsa(model_performance, Morris(num_trajectory = n_trajectories), ensemble_bounds)

# --- Save results ---
# morris_result.means: 4 x 25 (outputs x params)
# morris_result.variances: 4 x 25
writedlm(joinpath(RESULTS_DIR, "sensitivity_means.dat"), morris_result.means)
writedlm(joinpath(RESULTS_DIR, "sensitivity_variances.dat"), morris_result.variances)

# --- Print top parameters per output ---
output_names = ["Venus_mRNA_6h", "Venus_protein_12h", "GntR_mRNA_6h", "GntR_protein_12h"]
means_t = collect(morris_result.means')  # 25 x 4

for (j, oname) in enumerate(output_names)
    println("\n$(oname) (top 10):")
    sorted = sortperm(abs.(means_t[:, j]), rev = true)
    for rank in 1:min(10, N_PARAMETERS)
        i = sorted[rank]
        println("  $(rank). $(rpad(PARAMETER_NAMES[i], 22)) μ*=$(round(means_t[i,j], sigdigits=4))")
    end
end

println("\nDone. Sensitivity data saved to $(RESULTS_DIR)")
