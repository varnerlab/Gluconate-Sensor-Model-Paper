# estimate_parameters_run.jl - Single estimation run with per-cycle saves
#
# Usage: julia -t N --project=code code/scripts/estimate_parameters_run.jl <run_id>

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "GluconateBiosensor.jl"))
using .GluconateBiosensor
using ParetoEnsembles
using Statistics
using DelimitedFiles
using Random

# --- Get run ID from command line ---
run_id = parse(Int, ARGS[1])
RESULTS_DIR = joinpath(@__DIR__, "..", "results", "run_$(run_id)")
mkpath(RESULTS_DIR)

# --- Configuration ---
const N_CYCLES = 15
const N_CHAINS = 8
const MO_ITERATIONS = 50
const SO_ITERATIONS = 30
const RANK_CUTOFF = 4.0
const PERTURB_MO = 0.10
const PERTURB_RESEED = 0.05
const COOLING_RATE = 0.9
const N_OBJ = 6
const OBJ_NAMES = ["Venus_mRNA", "Venus_protein", "GntR_mRNA", "Venus_prot_0mM", "GntR_prot_reg", "Venus_noGntR"]

Random.seed!(run_id * 100 + 42)

# --- Setup ---
bio = load_biophysical_constants(joinpath(@__DIR__, "..", "src", "CellFree.json"))
genes = build_gene_info()
exp_data = load_experimental_data()

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

# --- Diverse initial guess ---
κ_best = clamp_to_bounds(default_initial_guess() .* (1.0 .+ 0.20 .* randn(N_PARAMETERS)))
ec_global = zeros(N_OBJ, 0)
pc_global = zeros(N_PARAMETERS, 0)

t_start = time()
println("[Run $(run_id)] Starting ($(N_CYCLES) cycles, $(N_CHAINS) chains, $(Threads.nthreads()) threads)")
flush(stdout)

for cycle in 1:N_CYCLES
    global κ_best, ec_global, pc_global
    t_cycle = time()

    # MO Phase
    initial_states = [clamp_to_bounds(κ_best .* (1.0 .+ PERTURB_MO .* randn(N_PARAMETERS))) for _ in 1:N_CHAINS]
    (EC_mo, PC_mo, RA_mo) = estimate_ensemble_parallel(
        OF, NF, accept_fn, cool_fn, initial_states;
        maximum_number_of_iterations = MO_ITERATIONS,
        rank_cutoff = RANK_CUTOFF,
        rng_seed = run_id * 1000 + cycle,
        show_trace = false,
    )
    ec_global = hcat(ec_global, EC_mo)
    pc_global = hcat(pc_global, PC_mo)

    # Identify worst objective
    global_ranks = rank_function(ec_global)
    front_idx = findall(global_ranks .== 0)
    if isempty(front_idx)
        front_idx = findall(global_ranks .<= 2)
    end
    mean_errors = vec(mean(ec_global[:, front_idx], dims = 2))
    worst_obj = argmax(mean_errors)

    total_err = vec(sum(ec_global[:, front_idx], dims = 1))
    best_mo_global_idx = front_idx[argmin(total_err)]
    κ_seed = pc_global[:, best_mo_global_idx]

    # SO Phase
    OF_so(pvec) = OF_single(pvec, worst_obj)
    κ_so_start = clamp_to_bounds(κ_seed .* (1.0 .+ 0.05 .* randn(N_PARAMETERS)))
    (EC_so, PC_so, RA_so) = estimate_ensemble(
        OF_so, NF, accept_fn, cool_fn, κ_so_start;
        maximum_number_of_iterations = SO_ITERATIONS,
        rank_cutoff = RANK_CUTOFF,
        show_trace = false,
    )
    best_so_idx = argmin(vec(EC_so))
    κ_so_best = PC_so[:, best_so_idx]
    full_err = evaluate_objectives(κ_so_best, bio, genes, exp_data)
    ec_global = hcat(ec_global, reshape(full_err, :, 1))
    pc_global = hcat(pc_global, reshape(κ_so_best, :, 1))
    κ_best = clamp_to_bounds(κ_so_best .* (1.0 .+ PERTURB_RESEED .* randn(N_PARAMETERS)))

    dt = round(time() - t_cycle, digits=1)
    elapsed = round((time() - t_start) / 60, digits=1)
    println("[Run $(run_id)] Cycle $(cycle)/$(N_CYCLES): archive=$(size(pc_global,2)), front=$(length(front_idx)), worst=$(OBJ_NAMES[worst_obj]), $(dt)s ($(elapsed) min total)")
    flush(stdout)

    # Per-cycle save
    writedlm(joinpath(RESULTS_DIR, "EC_all.dat"), ec_global)
    writedlm(joinpath(RESULTS_DIR, "PC_all.dat"), pc_global)
end

# Final filtered ensemble
final_ranks = rank_function(ec_global)
ensemble_idx = findall(final_ranks .<= 1)
writedlm(joinpath(RESULTS_DIR, "EC_final.dat"), ec_global[:, ensemble_idx])
writedlm(joinpath(RESULTS_DIR, "PC_final.dat"), pc_global[:, ensemble_idx])

elapsed = round((time() - t_start) / 60, digits=1)
println("[Run $(run_id)] Done in $(elapsed) min. Ensemble: $(length(ensemble_idx)), total archive: $(size(pc_global,2))")
flush(stdout)
