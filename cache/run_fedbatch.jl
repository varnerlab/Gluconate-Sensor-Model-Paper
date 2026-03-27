# run_fedbatch.jl - Synthetic fed-batch experiment
#
# Simulates 4 scenarios to test which resource is limiting:
#   1. Control (no boost)
#   2. Boost ε_X only at t_boost (replenish TX resources — mimics CP/energy addition)
#   3. Boost ε_L only at t_boost (replenish TL resources — mimics AA/tRNA addition)
#   4. Boost both ε_X and ε_L at t_boost (mimics full supplementation)
#
# Qualitative comparison with Li et al. (2017) fed-batch experiments in PURE system.
#
# Usage: julia --project=code scripts/run_fedbatch.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "GluconateBiosensor.jl"))
using .GluconateBiosensor
using DifferentialEquations
using DelimitedFiles
using Statistics

const RESULTS_DIR = joinpath(@__DIR__, "..", "results")

# --- Fed-batch parameters ---
const T_BOOST = 3.0       # hours — when to add fresh resources (ε_X is ~0.5 by this point)
const BOOST_FRAC = 0.5    # replenish 50% of depleted resource (partial fed-batch)

# --- Load ensemble ---
PC = readdlm(joinpath(RESULTS_DIR, "PC_final.dat"))
n_ensemble = size(PC, 2)
println("Loaded ensemble: $(n_ensemble) parameter sets")

# --- Setup ---
bio = load_biophysical_constants(joinpath(@__DIR__, "..", "src", "CellFree.json"))
genes = build_gene_info()
gluconate_conc = 10.0  # de-repressed condition

# --- Time grid ---
dt = 0.1
t_eval = collect(0.0:dt:12.0)
n_times = length(t_eval)

# --- Simulation with optional ε boost via callback ---
function simulate_fedbatch(params, bio, genes, gluconate_conc, t_boost;
                           boost_eX::Bool=false, boost_eL::Bool=false,
                           boost_frac::Float64=0.5,
                           tspan::Tuple{Float64,Float64}=(0.0, 12.0))

    x0 = zeros(N_STATES)
    x0[1:3] .= genes.initial_abundance
    x0[9] = 0.035
    x0[IDX_EPSILON_X] = 1.0
    x0[IDX_EPSILON_L] = 1.0

    p = (params, bio, genes, gluconate_conc)

    # Callback: at t_boost, replenish specified resource(s)
    condition(u, t, integrator) = t == t_boost
    function affect!(integrator)
        if boost_eX
            current = integrator.u[IDX_EPSILON_X]
            integrator.u[IDX_EPSILON_X] = current + boost_frac * (1.0 - current)
        end
        if boost_eL
            current = integrator.u[IDX_EPSILON_L]
            integrator.u[IDX_EPSILON_L] = current + boost_frac * (1.0 - current)
        end
    end

    cb = DiscreteCallback(condition, affect!)

    prob = ODEProblem(balances!, x0, tspan, p)
    sol = solve(prob, AutoTsit5(Rosenbrock23(autodiff=false));
                reltol=1e-8, abstol=1e-8, maxiters=1_000_000,
                save_everystep=true,
                callback=cb, tstops=[t_boost])

    return sol
end

# --- Scenario definitions ---
scenarios = [
    (name="control",   boost_eX=false, boost_eL=false),
    (name="boost_eX",  boost_eX=true,  boost_eL=false),
    (name="boost_eL",  boost_eX=false, boost_eL=true),
    (name="boost_both", boost_eX=true,  boost_eL=true),
]

# --- Storage ---
venus_protein = Dict{String, Matrix{Float64}}()
epsilon_X_traj = Dict{String, Matrix{Float64}}()
epsilon_L_traj = Dict{String, Matrix{Float64}}()
for s in scenarios
    venus_protein[s.name] = zeros(n_times, n_ensemble)
    epsilon_X_traj[s.name] = zeros(n_times, n_ensemble)
    epsilon_L_traj[s.name] = zeros(n_times, n_ensemble)
end

# --- Run ensemble for each scenario ---
n_success = 0
for k in 1:n_ensemble
    pvec = PC[:, k]
    params = vector_to_parameters(pvec)

    # Run all 4 scenarios for this parameter set
    all_ok = true
    sols = Dict{String, Any}()
    for s in scenarios
        sol = try
            simulate_fedbatch(params, bio, genes, gluconate_conc, T_BOOST;
                              boost_eX=s.boost_eX, boost_eL=s.boost_eL,
                              boost_frac=BOOST_FRAC)
        catch
            all_ok = false
            break
        end
        if sol.t[end] < 11.5
            all_ok = false
            break
        end
        sols[s.name] = sol
    end

    if !all_ok
        continue
    end

    global n_success += 1
    for s in scenarios
        sol = sols[s.name]
        for (i, t) in enumerate(t_eval)
            state = sol(t)
            venus_protein[s.name][i, n_success] = state[8]
            epsilon_X_traj[s.name][i, n_success] = state[10]
            epsilon_L_traj[s.name][i, n_success] = state[11]
        end
    end

    if k % 500 == 0
        println("  Processed $k / $n_ensemble parameter sets ($n_success successful)")
    end
end

println("Successful simulations: $(n_success)/$(n_ensemble)")

# --- Trim and save ---
for s in scenarios
    venus_protein[s.name] = venus_protein[s.name][:, 1:n_success]
    epsilon_X_traj[s.name] = epsilon_X_traj[s.name][:, 1:n_success]
    epsilon_L_traj[s.name] = epsilon_L_traj[s.name][:, 1:n_success]

    writedlm(joinpath(RESULTS_DIR, "fedbatch_venus_$(s.name).dat"), venus_protein[s.name])
    writedlm(joinpath(RESULTS_DIR, "fedbatch_eX_$(s.name).dat"), epsilon_X_traj[s.name])
    writedlm(joinpath(RESULTS_DIR, "fedbatch_eL_$(s.name).dat"), epsilon_L_traj[s.name])
end
writedlm(joinpath(RESULTS_DIR, "fedbatch_time.dat"), t_eval)

# --- Summary ---
println("\n=== Fed-batch results at t = 12 hr ===")
println("Boost at t = $(T_BOOST) hr, replenishment fraction = $(BOOST_FRAC)")
println()
for s in scenarios
    v_mean = mean(venus_protein[s.name][end, :])
    v_std = std(venus_protein[s.name][end, :])
    eX_mean = mean(epsilon_X_traj[s.name][end, :])
    eL_mean = mean(epsilon_L_traj[s.name][end, :])
    println("$(rpad(s.name, 12)) Venus protein: $(round(v_mean, digits=3)) ± $(round(v_std, digits=3)) μM  |  ε_X=$(round(eX_mean, digits=3))  ε_L=$(round(eL_mean, digits=3))")
end

# Fold-changes relative to control
ctrl_mean = mean(venus_protein["control"][end, :])
println("\n=== Fold-change vs control ===")
for s in scenarios
    if s.name == "control"
        continue
    end
    fc = mean(venus_protein[s.name][end, :]) / ctrl_mean
    println("$(rpad(s.name, 12)) $(round(fc, digits=3))×")
end

println("\nDone. Results saved to $(RESULTS_DIR)/fedbatch_*.dat")
