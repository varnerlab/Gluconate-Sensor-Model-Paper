# run_fedbatch_sweep.jl - Sweep supplementation levels for ε_X-only and ε_L-only
#
# Usage: julia --project=code scripts/run_fedbatch_sweep.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "GluconateBiosensor.jl"))
using .GluconateBiosensor
using DifferentialEquations
using DelimitedFiles
using Statistics

const RESULTS_DIR = joinpath(@__DIR__, "..", "results")
const T_BOOST = 3.0

# --- Load ensemble ---
PC = readdlm(joinpath(RESULTS_DIR, "PC_final.dat"))
n_ensemble = size(PC, 2)
bio = load_biophysical_constants(joinpath(@__DIR__, "..", "src", "CellFree.json"))
genes = build_gene_info()
gluconate_conc = 10.0

dt = 0.1
t_eval = collect(0.0:dt:12.0)
n_times = length(t_eval)

# Sweep levels
boost_fracs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

function simulate_fedbatch(params, bio, genes, gluconate_conc, t_boost;
                           boost_eX::Bool=false, boost_eL::Bool=false,
                           boost_frac::Float64=0.5)

    x0 = zeros(N_STATES)
    x0[1:3] .= genes.initial_abundance
    x0[9] = 0.035
    x0[IDX_EPSILON_X] = 1.0
    x0[IDX_EPSILON_L] = 1.0
    p = (params, bio, genes, gluconate_conc)

    condition(u, t, integrator) = t == t_boost
    function affect!(integrator)
        if boost_eX
            c = integrator.u[IDX_EPSILON_X]
            integrator.u[IDX_EPSILON_X] = c + boost_frac * (1.0 - c)
        end
        if boost_eL
            c = integrator.u[IDX_EPSILON_L]
            integrator.u[IDX_EPSILON_L] = c + boost_frac * (1.0 - c)
        end
    end
    cb = DiscreteCallback(condition, affect!)

    prob = ODEProblem(balances!, x0, (0.0, 12.0), p)
    sol = solve(prob, AutoTsit5(Rosenbrock23(autodiff=false));
                reltol=1e-8, abstol=1e-8, maxiters=1_000_000,
                save_everystep=true, callback=cb, tstops=[t_boost])
    return sol
end

# Pre-convert parameters
all_params = [vector_to_parameters(PC[:, k]) for k in 1:n_ensemble]

# Storage: rows = boost_fracs, cols = [mean_eX, std_eX, mean_eL, std_eL]
results_eX = zeros(length(boost_fracs), 2)
results_eL = zeros(length(boost_fracs), 2)

for (bi, bf) in enumerate(boost_fracs)
    venus_eX = zeros(n_ensemble)
    venus_eL = zeros(n_ensemble)
    n_ok = Threads.Atomic{Int}(0)

    Threads.@threads for k in 1:n_ensemble
        params = all_params[k]
        
        sol_eX = try
            simulate_fedbatch(params, bio, genes, gluconate_conc, T_BOOST;
                              boost_eX=true, boost_eL=false, boost_frac=bf)
        catch; nothing; end

        sol_eL = try
            simulate_fedbatch(params, bio, genes, gluconate_conc, T_BOOST;
                              boost_eX=false, boost_eL=true, boost_frac=bf)
        catch; nothing; end

        if sol_eX !== nothing && sol_eX.t[end] >= 11.5
            venus_eX[k] = sol_eX(12.0)[8]
        end
        if sol_eL !== nothing && sol_eL.t[end] >= 11.5
            venus_eL[k] = sol_eL(12.0)[8]
        end
    end

    # Filter zeros (failed sims)
    vX = filter(x -> x > 0, venus_eX)
    vL = filter(x -> x > 0, venus_eL)
    
    results_eX[bi, 1] = mean(vX)
    results_eX[bi, 2] = std(vX)
    results_eL[bi, 1] = mean(vL)
    results_eL[bi, 2] = std(vL)

    println("boost_frac=$(bf): eX_boost → $(round(mean(vX), digits=3))±$(round(std(vX), digits=3))  |  eL_boost → $(round(mean(vL), digits=3))±$(round(std(vL), digits=3))  (n=$(length(vX)))")
end

# Save
writedlm(joinpath(RESULTS_DIR, "fedbatch_sweep_fracs.dat"), boost_fracs)
writedlm(joinpath(RESULTS_DIR, "fedbatch_sweep_eX.dat"), results_eX)
writedlm(joinpath(RESULTS_DIR, "fedbatch_sweep_eL.dat"), results_eL)

println("\nDone. Sweep results saved.")
