# estimate_parameters.jl - MO/SO cycling parameter estimation
#
# Usage: julia -t 8 --project=code code/scripts/estimate_parameters.jl
#
# 6 objectives:
#   1. Venus mRNA SSE (10 mM gluconate)
#   2. Venus protein SSE (10 mM gluconate)
#   3. GntR mRNA SSE (10 mM gluconate)
#   4. Venus protein SSE (0 mM gluconate — full repression floor)
#   5. GntR protein regularization (keep in [5, 15] μM)
#   6. Venus protein SSE (no-GntR control — unrepressed ceiling)

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "GluconateBiosensor.jl"))
using .GluconateBiosensor
using ParetoEnsembles
using Statistics
using DelimitedFiles
using Random

# --- Configuration ---
const N_CYCLES = 15
const N_CHAINS = 5
const MO_ITERATIONS = 50
const SO_ITERATIONS = 30
const RANK_CUTOFF = 4.0
const PERTURB_MO = 0.10
const PERTURB_RESEED = 0.05
const COOLING_RATE = 0.9
const RESULTS_DIR = joinpath(@__DIR__, "..", "results")
const N_OBJ = 6
const OBJ_NAMES = ["Venus_mRNA", "Venus_protein", "GntR_mRNA", "Venus_prot_0mM", "GntR_prot_reg", "Venus_noGntR"]

# --- Setup ---
bio = load_biophysical_constants(joinpath(@__DIR__, "..", "src", "CellFree.json"))
genes = build_gene_info()
exp_data = load_experimental_data()

# Objective function (returns N_OBJ×1 Matrix for ParetoEnsembles)
function OF(pvec)
    errs = evaluate_objectives(pvec, bio, genes, exp_data)
    return reshape(errs, :, 1)
end

function OF_single(pvec, obj_idx)
    errs = evaluate_objectives(pvec, bio, genes, exp_data)
    return reshape([errs[obj_idx]], :, 1)
end

function NF(pvec)
    new_pvec = pvec .* (1.0 .+ 0.05 .* randn(N_PARAMETERS))
    return clamp_to_bounds(new_pvec)
end

accept_fn(rank_array, T) = exp(-rank_array[end] / T)
cool_fn(T) = COOLING_RATE * T

# --- Cycling loop ---
κ_best = default_initial_guess()
ec_global = zeros(N_OBJ, 0)
pc_global = zeros(N_PARAMETERS, 0)

println("Starting MO/SO cycling estimation with $(N_CYCLES) cycles, $(N_CHAINS) chains, $(N_OBJ) objectives...")
println("Running with $(Threads.nthreads()) threads")

for cycle in 1:N_CYCLES
    global κ_best, ec_global, pc_global

    println("\n" * "="^60)
    println("CYCLE $(cycle)/$(N_CYCLES)")
    println("="^60)

    # ---- Phase 1: Multi-objective (parallel chains) ----
    println("\n--- MO Phase ($(N_CHAINS) chains, $(MO_ITERATIONS) iterations) ---")
    initial_states = [clamp_to_bounds(κ_best .* (1.0 .+ PERTURB_MO .* randn(N_PARAMETERS))) for _ in 1:N_CHAINS]

    (EC_mo, PC_mo, RA_mo) = estimate_ensemble_parallel(
        OF, NF, accept_fn, cool_fn, initial_states;
        maximum_number_of_iterations = MO_ITERATIONS,
        rank_cutoff = RANK_CUTOFF,
        rng_seed = cycle,
        show_trace = false,
    )

    println("  MO archive size: $(size(PC_mo, 2))")
    println("  MO Pareto front size: $(count(RA_mo .== 0))")

    ec_global = hcat(ec_global, EC_mo)
    pc_global = hcat(pc_global, PC_mo)

    # ---- Identify worst objective ----
    global_ranks = rank_function(ec_global)
    front_idx = findall(global_ranks .== 0)
    if isempty(front_idx)
        front_idx = findall(global_ranks .<= 2)
    end
    mean_errors = vec(mean(ec_global[:, front_idx], dims = 2))
    worst_obj = argmax(mean_errors)
    println("  Worst objective: $(worst_obj) ($(OBJ_NAMES[worst_obj]))")
    println("  Mean errors on front: $(round.(mean_errors, digits=4))")

    total_err = vec(sum(ec_global[:, front_idx], dims = 1))
    best_mo_global_idx = front_idx[argmin(total_err)]
    κ_seed = pc_global[:, best_mo_global_idx]

    # ---- Phase 2: Single-objective (drill down on worst) ----
    println("\n--- SO Phase (objective $(worst_obj): $(OBJ_NAMES[worst_obj]), $(SO_ITERATIONS) iterations) ---")
    OF_so(pvec) = OF_single(pvec, worst_obj)
    κ_so_start = clamp_to_bounds(κ_seed .* (1.0 .+ 0.05 .* randn(N_PARAMETERS)))

    (EC_so, PC_so, RA_so) = estimate_ensemble(
        OF_so, NF, accept_fn, cool_fn, κ_so_start;
        maximum_number_of_iterations = SO_ITERATIONS,
        rank_cutoff = RANK_CUTOFF,
        show_trace = false,
    )

    println("  SO archive size: $(size(PC_so, 2))")

    best_so_idx = argmin(vec(EC_so))
    κ_so_best = PC_so[:, best_so_idx]
    full_err = evaluate_objectives(κ_so_best, bio, genes, exp_data)
    println("  SO best full errors: $(round.(full_err, digits=4))")

    ec_global = hcat(ec_global, reshape(full_err, :, 1))
    pc_global = hcat(pc_global, reshape(κ_so_best, :, 1))

    κ_best = clamp_to_bounds(κ_so_best .* (1.0 .+ PERTURB_RESEED .* randn(N_PARAMETERS)))

    writedlm(joinpath(RESULTS_DIR, "EC_cycle$(cycle).dat"), ec_global)
    writedlm(joinpath(RESULTS_DIR, "PC_cycle$(cycle).dat"), pc_global)
    println("  Saved cycle $(cycle) results (global archive: $(size(pc_global, 2)) solutions)")
end

# --- Final ensemble ---
println("\n" * "="^60)
println("FINAL ENSEMBLE")
println("="^60)

final_ranks = rank_function(ec_global)
ensemble_idx = findall(final_ranks .<= 1)
println("Solutions with rank ≤ 1: $(length(ensemble_idx))")
println("Pareto-optimal (rank 0): $(count(final_ranks .== 0))")

writedlm(joinpath(RESULTS_DIR, "EC_final.dat"), ec_global[:, ensemble_idx])
writedlm(joinpath(RESULTS_DIR, "PC_final.dat"), pc_global[:, ensemble_idx])
writedlm(joinpath(RESULTS_DIR, "RA_final.dat"), final_ranks[ensemble_idx])

println("\nEnsemble parameter statistics:")
for (i, name) in enumerate(PARAMETER_NAMES)
    vals = pc_global[i, ensemble_idx]
    μ = mean(vals); σ = std(vals)
    cv = σ / (abs(μ) + eps()) * 100
    println("  $(rpad(name, 22)) mean=$(rpad(round(μ, sigdigits=4), 12)) std=$(rpad(round(σ, sigdigits=4), 12)) CV=$(round(cv, digits=1))%")
end

println("\nDone. Results saved to $(RESULTS_DIR)")
