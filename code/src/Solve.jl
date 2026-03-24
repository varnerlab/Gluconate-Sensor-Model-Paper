# Solve.jl - ODE solver wrapper

using DifferentialEquations

"""
    simulate(params, bio, genes, gluconate_conc; tspan=(0.0, 12.0))

Simulate the gluconate biosensor circuit model.

Returns the DifferentialEquations.jl solution object, which is callable:
  `sol(t)` returns the state vector at any time t in the span.

State indices: 1:3 = genes, 4:6 = mRNA, 7:9 = protein, 10 = ε_X, 11 = ε_L
"""
function simulate(params::ModelParameters, bio::BiophysicalConstants,
                  genes::GeneInfo, gluconate_conc::Float64;
                  tspan::Tuple{Float64,Float64} = (0.0, 12.0))

    # Initial conditions
    x0 = zeros(N_STATES)
    x0[1:3] .= genes.initial_abundance       # gene concentrations (constant)
    # mRNA and protein start at 0, except sigma_70 protein
    x0[9] = 0.035                             # sigma_70 protein initial (μM) — added exogenously
    x0[IDX_EPSILON_X] = 1.0                   # TX resources fully available
    x0[IDX_EPSILON_L] = 1.0                   # TL resources fully available

    # Pack parameters
    p = (params, bio, genes, gluconate_conc)

    # Solve
    prob = ODEProblem(balances!, x0, tspan, p)
    sol = solve(prob, AutoTsit5(Rosenbrock23(autodiff = false));
                reltol = 1e-8, abstol = 1e-8, maxiters = 1_000_000,
                save_everystep = true)

    return sol
end
