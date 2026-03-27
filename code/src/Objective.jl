# Objective.jl - Objective functions for parameter estimation

"""
    evaluate_objectives(pvec, bio, genes, exp_data; gluconate_training=10.0) → Vector{Float64}(6)

Evaluate the 6 training objectives:
  1. Venus mRNA SSE at 10 mM gluconate (5 time points)
  2. Venus protein SSE at 10 mM gluconate (dense time course)
  3. GntR mRNA SSE at 10 mM gluconate (5 time points)
  4. Venus protein SSE at 0 mM gluconate — full repression (dense time course)
  5. GntR protein regularization — penalizes if outside [1, 20] μM
  6. Venus protein SSE for no-GntR control — unrepressed ceiling (dense time course)

Returns a vector of 6 error values. Returns fill(1e12, 6) if simulation fails.
"""
function evaluate_objectives(pvec::Vector{Float64}, bio::BiophysicalConstants,
                              genes::GeneInfo, exp_data::ExperimentalData;
                              gluconate_training::Float64 = 10.0)::Vector{Float64}

    BIG = 1e12
    N_OBJ = 6

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

    # --- Simulation 3: no-GntR control (unrepressed ceiling) ---
    # Simulate with GntR gene at 0 nM DNA — only Venus gene present
    genes_noGntR = GeneInfo(
        genes.names,
        genes.coding_length_nt,
        genes.protein_length_aa,
        [0.0, genes.initial_abundance[2], genes.initial_abundance[3]],  # GntR DNA = 0
    )
    sol_noGntR = try
        simulate(params, bio, genes_noGntR, 0.0)  # no gluconate needed (no GntR)
    catch
        return fill(BIG, N_OBJ)
    end
    if string(sol_noGntR.retcode) ∉ ("Success", "Terminated") || sol_noGntR.t[end] < 11.5
        return fill(BIG, N_OBJ)
    end

    # --- Objective 1: Venus mRNA weighted SSE (10 mM) ---
    # Weighted by 1/σ² so each point contributes proportionally to its precision
    err_venus_mRNA = 0.0
    for (i, t) in enumerate(exp_data.mRNA_time)
        σ = exp_data.mRNA_venus_std[i]
        σ < 1e-9 && continue  # skip t=0
        sim_val = sol_10(t)[5]
        exp_val = exp_data.mRNA_venus_mean[i]
        err_venus_mRNA += ((sim_val - exp_val) / σ)^2
    end

    # --- Objective 2: Venus protein SSE (10 mM) ---
    err_venus_protein = 0.0
    for (i, t) in enumerate(exp_data.protein_time)
        if t > sol_10.t[end]; break; end
        sim_val = sol_10(t)[8]
        exp_val = exp_data.protein_venus_mean[i]
        err_venus_protein += (sim_val - exp_val)^2
    end

    # --- Objective 3: GntR mRNA weighted SSE (10 mM) ---
    err_gntr_mRNA = 0.0
    for (i, t) in enumerate(exp_data.mRNA_time)
        σ = exp_data.mRNA_gntr_std[i]
        σ < 1e-9 && continue  # skip t=0
        sim_val = sol_10(t)[4]
        exp_val = exp_data.mRNA_gntr_mean[i]
        err_gntr_mRNA += ((sim_val - exp_val) / σ)^2
    end

    # --- Objective 4: Venus protein SSE (0 mM — full repression) ---
    err_venus_0mM = 0.0
    for (i, t) in enumerate(exp_data.protein_time)
        if t > sol_0.t[end]; break; end
        sim_val = sol_0(t)[8]
        exp_val = exp_data.protein_venus_0mM_mean[i]
        err_venus_0mM += (sim_val - exp_val)^2
    end

    # --- Objective 5: GntR protein regularization ---
    gntr_protein_12h = sol_10(12.0)[7]
    err_gntr_reg = 0.0
    if gntr_protein_12h < 1.0
        err_gntr_reg = (gntr_protein_12h - 1.0)^2
    elseif gntr_protein_12h > 20.0
        err_gntr_reg = (gntr_protein_12h - 20.0)^2
    end

    # --- Objective 6: Venus protein SSE (no-GntR — unrepressed ceiling) ---
    err_venus_noGntR = 0.0
    for (i, t) in enumerate(exp_data.protein_time)
        if t > sol_noGntR.t[end]; break; end
        sim_val = sol_noGntR(t)[8]
        exp_val = exp_data.protein_venus_noGntR_mean[i]
        err_venus_noGntR += (sim_val - exp_val)^2
    end

    return [err_venus_mRNA, err_venus_protein, err_gntr_mRNA, err_venus_0mM, err_gntr_reg, err_venus_noGntR]
end

"""
    evaluate_single_objective(pvec, bio, genes, exp_data, obj_index; kwargs...) → Vector{Float64}(1)

Evaluate a single objective (for SO phase of cycling).
"""
function evaluate_single_objective(pvec::Vector{Float64}, bio::BiophysicalConstants,
                                    genes::GeneInfo, exp_data::ExperimentalData,
                                    obj_index::Int; kwargs...)::Vector{Float64}
    errs = evaluate_objectives(pvec, bio, genes, exp_data; kwargs...)
    return [errs[obj_index]]
end
