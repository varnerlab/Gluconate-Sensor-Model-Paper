# Objective.jl - Objective functions for parameter estimation

"""
    evaluate_objectives(pvec, bio, genes, exp_data; gluconate_training=10.0) → Vector{Float64}(5)

Evaluate the 5 training objectives:
  1. Venus mRNA SSE at 10 mM gluconate (5 time points)
  2. Venus protein SSE at 10 mM gluconate (dense time course)
  3. GntR mRNA SSE at 10 mM gluconate (5 time points)
  4. Venus protein SSE at 0 mM gluconate — full repression (dense time course)
  5. GntR protein regularization — penalizes if GntR protein at 12h
     is outside a plausible range, enforcing ensemble consistency

Returns a vector of 5 error values. Returns fill(1e12, 5) if simulation fails.
"""
function evaluate_objectives(pvec::Vector{Float64}, bio::BiophysicalConstants,
                              genes::GeneInfo, exp_data::ExperimentalData;
                              gluconate_training::Float64 = 10.0)::Vector{Float64}

    BIG = 1e12
    N_OBJ = 5

    # Convert parameter vector to struct
    params = vector_to_parameters(pvec)

    # --- Simulation 1: 10 mM gluconate (training condition) ---
    sol_10 = try
        simulate(params, bio, genes, gluconate_training)
    catch
        return fill(BIG, N_OBJ)
    end
    if string(sol_10.retcode) ∉ ("Success", "Terminated") || sol_10.t[end] < 11.5
        return fill(BIG, N_OBJ)
    end

    # --- Simulation 2: 0 mM gluconate (full repression condition) ---
    sol_0 = try
        simulate(params, bio, genes, 0.0)
    catch
        return fill(BIG, N_OBJ)
    end
    if string(sol_0.retcode) ∉ ("Success", "Terminated") || sol_0.t[end] < 11.5
        return fill(BIG, N_OBJ)
    end

    # --- Objective 1: Venus mRNA SSE (10 mM) ---
    err_venus_mRNA = 0.0
    for (i, t) in enumerate(exp_data.mRNA_time)
        sim_val = sol_10(t)[5]  # mRNA_Venus
        exp_val = exp_data.mRNA_venus_mean[i]
        err_venus_mRNA += (sim_val - exp_val)^2
    end

    # --- Objective 2: Venus protein SSE (10 mM) ---
    err_venus_protein = 0.0
    for (i, t) in enumerate(exp_data.protein_time)
        if t > sol_10.t[end]
            break
        end
        sim_val = sol_10(t)[8]  # protein_Venus
        exp_val = exp_data.protein_venus_mean[i]
        err_venus_protein += (sim_val - exp_val)^2
    end

    # --- Objective 3: GntR mRNA SSE (10 mM) ---
    err_gntr_mRNA = 0.0
    for (i, t) in enumerate(exp_data.mRNA_time)
        sim_val = sol_10(t)[4]  # mRNA_GntR
        exp_val = exp_data.mRNA_gntr_mean[i]
        err_gntr_mRNA += (sim_val - exp_val)^2
    end

    # --- Objective 4: Venus protein SSE (0 mM — full repression) ---
    err_venus_0mM = 0.0
    for (i, t) in enumerate(exp_data.protein_time)
        if t > sol_0.t[end]
            break
        end
        sim_val = sol_0(t)[8]  # protein_Venus at 0 mM gluconate
        exp_val = exp_data.protein_venus_0mM_mean[i]
        err_venus_0mM += (sim_val - exp_val)^2
    end

    # --- Objective 5: GntR protein regularization ---
    # Penalize GntR protein at 12h if outside a plausible range.
    # Venus protein at 12h is ~1.7 μM. GntR has higher mRNA (~30x) but shares
    # resources. A reasonable range for GntR protein: 1-20 μM.
    # Use a soft quadratic penalty outside [1, 20] μM.
    gntr_protein_12h = sol_10(12.0)[7]
    gntr_target_low = 1.0    # μM
    gntr_target_high = 20.0  # μM
    err_gntr_reg = 0.0
    if gntr_protein_12h < gntr_target_low
        err_gntr_reg = (gntr_protein_12h - gntr_target_low)^2
    elseif gntr_protein_12h > gntr_target_high
        err_gntr_reg = (gntr_protein_12h - gntr_target_high)^2
    end

    return [err_venus_mRNA, err_venus_protein, err_gntr_mRNA, err_venus_0mM, err_gntr_reg]
end

"""
    evaluate_single_objective(pvec, bio, genes, exp_data, obj_index; kwargs...) → Vector{Float64}(1)

Evaluate a single objective (for SO phase of cycling).
Returns a 1-element vector for compatibility with ParetoEnsembles.jl.
"""
function evaluate_single_objective(pvec::Vector{Float64}, bio::BiophysicalConstants,
                                    genes::GeneInfo, exp_data::ExperimentalData,
                                    obj_index::Int; kwargs...)::Vector{Float64}
    errs = evaluate_objectives(pvec, bio, genes, exp_data; kwargs...)
    return [errs[obj_index]]
end
