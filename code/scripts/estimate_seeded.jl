# estimate_seeded.jl - Warm-start estimation from seed solutions
#
# Usage: julia -t N --project=code code/scripts/estimate_seeded.jl <run_id>

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "GluconateBiosensor.jl"))
using .GluconateBiosensor
using ParetoEnsembles
using Statistics
using DelimitedFiles
using Random

run_id = parse(Int, ARGS[1])
RESULTS_DIR = joinpath(@__DIR__, "..", "results", "seed_run_$(run_id)")
mkpath(RESULTS_DIR)

const N_CYCLES = 15
const N_CHAINS = 8
const MO_ITERATIONS = 50
const SO_ITERATIONS = 30
const RANK_CUTOFF = 4.0
const COOLING_RATE = 0.9
const N_OBJ = 6
const OBJ_NAMES = ["Venus_mRNA", "Venus_protein", "GntR_mRNA", "Venus_prot_0mM", "GntR_prot_reg", "Venus_noGntR"]

Random.seed!(run_id * 100 + 7)

bio = load_biophysical_constants(joinpath(@__DIR__, "..", "src", "CellFree.json"))
genes = build_gene_info()
exp_data = load_experimental_data()

OF(pvec) = reshape(evaluate_objectives(pvec, bio, genes, exp_data), :, 1)
OF_single(pvec, obj_idx) = reshape([evaluate_objectives(pvec, bio, genes, exp_data)[obj_idx]], :, 1)
NF(pvec) = clamp_to_bounds(pvec .* (1.0 .+ 0.05 .* randn(N_PARAMETERS)))
accept_fn(ra, T) = exp(-ra[end] / T)
cool_fn(T) = COOLING_RATE * T

# Load seeds and pick a random one for this run
seeds = readdlm(joinpath(@__DIR__, "..", "results", "seeds.dat"))
n_seeds = size(seeds, 2)
seed_idx = rand(1:n_seeds)
κ_best = clamp_to_bounds(seeds[:, seed_idx] .* (1.0 .+ 0.02 .* randn(N_PARAMETERS)))

ec_global = zeros(N_OBJ, 0)
pc_global = zeros(N_PARAMETERS, 0)

t_start = time()
println("[Seed run $(run_id)] Starting from seed $seed_idx ($(N_CYCLES) cycles, $(N_CHAINS) chains, $(Threads.nthreads()) threads)")
flush(stdout)

for cycle in 1:N_CYCLES
    global κ_best, ec_global, pc_global
    t_cycle = time()

    # MO: seed chains from perturbed best + random seeds
    initial_states = Vector{Float64}[]
    push!(initial_states, clamp_to_bounds(κ_best .* (1.0 .+ 0.05 .* randn(N_PARAMETERS))))
    for _ in 2:N_CHAINS
        si = rand(1:n_seeds)
        push!(initial_states, clamp_to_bounds(seeds[:, si] .* (1.0 .+ 0.05 .* randn(N_PARAMETERS))))
    end

    (EC_mo, PC_mo, RA_mo) = estimate_ensemble_parallel(
        OF, NF, accept_fn, cool_fn, initial_states;
        maximum_number_of_iterations = MO_ITERATIONS,
        rank_cutoff = RANK_CUTOFF, rng_seed = run_id * 1000 + cycle, show_trace = false)
    ec_global = hcat(ec_global, EC_mo)
    pc_global = hcat(pc_global, PC_mo)

    global_ranks = rank_function(ec_global)
    front_idx = findall(global_ranks .== 0)
    isempty(front_idx) && (front_idx = findall(global_ranks .<= 2))
    mean_errors = vec(mean(ec_global[:, front_idx], dims = 2))
    worst_obj = argmax(mean_errors)
    total_err = vec(sum(ec_global[:, front_idx], dims = 1))
    κ_seed = pc_global[:, front_idx[argmin(total_err)]]

    OF_so(pvec) = OF_single(pvec, worst_obj)
    κ_so_start = clamp_to_bounds(κ_seed .* (1.0 .+ 0.05 .* randn(N_PARAMETERS)))
    (EC_so, PC_so, RA_so) = estimate_ensemble(
        OF_so, NF, accept_fn, cool_fn, κ_so_start;
        maximum_number_of_iterations = SO_ITERATIONS, rank_cutoff = RANK_CUTOFF, show_trace = false)
    best_so_idx = argmin(vec(EC_so))
    κ_so_best = PC_so[:, best_so_idx]
    full_err = evaluate_objectives(κ_so_best, bio, genes, exp_data)
    ec_global = hcat(ec_global, reshape(full_err, :, 1))
    pc_global = hcat(pc_global, reshape(κ_so_best, :, 1))
    κ_best = clamp_to_bounds(κ_so_best .* (1.0 .+ 0.05 .* randn(N_PARAMETERS)))

    dt = round(time() - t_cycle, digits=1)
    elapsed = round((time() - t_start) / 60, digits=1)
    println("[Seed run $(run_id)] Cycle $(cycle)/$(N_CYCLES): archive=$(size(pc_global,2)), worst=$(OBJ_NAMES[worst_obj]), $(dt)s ($(elapsed) min)")
    flush(stdout)

    writedlm(joinpath(RESULTS_DIR, "EC_all.dat"), ec_global)
    writedlm(joinpath(RESULTS_DIR, "PC_all.dat"), pc_global)
end

final_ranks = rank_function(ec_global)
ensemble_idx = findall(final_ranks .<= 1)
writedlm(joinpath(RESULTS_DIR, "EC_final.dat"), ec_global[:, ensemble_idx])
writedlm(joinpath(RESULTS_DIR, "PC_final.dat"), pc_global[:, ensemble_idx])

elapsed = round((time() - t_start) / 60, digits=1)
println("[Seed run $(run_id)] Done in $(elapsed) min. Ensemble: $(length(ensemble_idx)), archive: $(size(pc_global,2))")
flush(stdout)
