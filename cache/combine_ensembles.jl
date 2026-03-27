# combine_ensembles.jl - Merge ensembles from parallel runs and re-rank
#
# Reads PC_all.dat and EC_all.dat from each run directory,
# stacks them, computes global Pareto rank, and filters.
#
# Usage: julia --project=code code/scripts/combine_ensembles.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "GluconateBiosensor.jl"))
using .GluconateBiosensor
using ParetoEnsembles
using Statistics
using DelimitedFiles

const RESULTS_DIR = joinpath(@__DIR__, "..", "results")

# --- Find and load all run directories ---
run_dirs = filter(d -> startswith(basename(d), "run_") && isdir(d),
                  readdir(RESULTS_DIR, join=true))
sort!(run_dirs)

println("Found $(length(run_dirs)) run directories:")
ec_combined = zeros(6, 0)
pc_combined = zeros(N_PARAMETERS, 0)

for rd in run_dirs
    ec_file = joinpath(rd, "EC_all.dat")
    pc_file = joinpath(rd, "PC_all.dat")
    if isfile(ec_file) && isfile(pc_file)
        ec = readdlm(ec_file)
        pc = readdlm(pc_file)
        println("  $(basename(rd)): $(size(pc, 2)) solutions")
        global ec_combined = hcat(ec_combined, ec)
        global pc_combined = hcat(pc_combined, pc)
    else
        println("  $(basename(rd)): MISSING data files, skipping")
    end
end

println("\nTotal solutions: $(size(pc_combined, 2))")

# --- Global Pareto ranking ---
println("Computing global Pareto ranking...")
global_ranks = rank_function(ec_combined)

# Filter: rank ≤ 1
ensemble_idx = findall(global_ranks .<= 1)
println("Solutions with rank ≤ 1: $(length(ensemble_idx))")
println("Pareto-optimal (rank 0): $(count(global_ranks .== 0))")

# --- Save combined results ---
writedlm(joinpath(RESULTS_DIR, "EC_final.dat"), ec_combined[:, ensemble_idx])
writedlm(joinpath(RESULTS_DIR, "PC_final.dat"), pc_combined[:, ensemble_idx])
writedlm(joinpath(RESULTS_DIR, "RA_final.dat"), global_ranks[ensemble_idx])

# --- Parameter statistics ---
println("\nEnsemble parameter statistics (N = $(length(ensemble_idx))):")
for (i, name) in enumerate(PARAMETER_NAMES)
    vals = pc_combined[i, ensemble_idx]
    μ = mean(vals); σ = std(vals)
    println("  $(rpad(name, 22)) $(rpad(round(μ, sigdigits=4), 14)) ± $(round(σ, sigdigits=4))")
end

# --- Quick validation: simulate ensemble mean ---
println("\nValidating ensemble mean...")
bio = load_biophysical_constants(joinpath(@__DIR__, "..", "src", "CellFree.json"))
genes = build_gene_info()
pvec_mean = vec(mean(pc_combined[:, ensemble_idx], dims=2))
params = vector_to_parameters(pvec_mean)
sol = try
    simulate(params, bio, genes, 10.0)
catch e
    println("  WARNING: mean parameter simulation failed: $e")
    nothing
end
if sol !== nothing && sol.t[end] >= 11.5
    println("  Venus protein at 12h: $(round(sol(12.0)[8], digits=3)) μM (exp ~1.7)")
    println("  GntR protein at 12h: $(round(sol(12.0)[7], digits=3)) μM")
    println("  ε_X at 12h: $(round(sol(12.0)[10], digits=4))")
    println("  ε_L at 12h: $(round(sol(12.0)[11], digits=4))")
end

println("\nDone. Combined ensemble saved to $(RESULTS_DIR)")
