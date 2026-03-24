# Model.jl - ODE system for the gluconate biosensor circuit
#
# State vector (11 states):
#   x[1:3]  = gene concentrations [GntR, Venus, sigma_70] (constant)
#   x[4:6]  = mRNA concentrations [GntR, Venus, sigma_70] (nM → μM internally)
#   x[7:9]  = protein concentrations [GntR, Venus, sigma_70] (μM)
#   x[10]   = ε_X: fraction of transcription consumable resources remaining
#   x[11]   = ε_L: fraction of translation consumable resources remaining

const N_STATES = 11
const N_GENES = 3

# Gene indices
const IDX_GENE = 1:3
const IDX_MRNA = 4:6
const IDX_PROTEIN = 7:9
const IDX_EPSILON_X = 10
const IDX_EPSILON_L = 11

# Species sub-indices (within their category)
const I_GNTR = 1
const I_VENUS = 2
const I_SIGMA70 = 3

"""
    transcription_control(x, params, gluconate_conc) → Vector{Float64}(3)

Compute the transcription control array ū for each gene using
the Boltzmann partition function formalism.

Returns values in [0, 1] for [GntR, Venus, sigma_70].
"""
function transcription_control(x::AbstractVector, params::ModelParameters, gluconate_conc::Float64)

    # Alias protein concentrations
    protein_GntR = x[7]
    protein_sigma70 = x[9]

    # Gluconate binding: fraction of GntR bound to gluconate (inactive as repressor)
    n_gluc = params.n_gluconate
    K_gluc = params.K_gluconate
    f_bound = (gluconate_conc^n_gluc) / (gluconate_conc^n_gluc + K_gluc^n_gluc + eps())
    effective_GntR = (1.0 - f_bound) * protein_GntR

    # W values from pseudo-energies: W = exp(-dG)
    W_GntR_RNAP = exp(-params.dG_GntR_RNAP)
    W_GntR_sigma70 = exp(-params.dG_GntR_sigma70)
    W_Venus_RNAP = exp(-params.dG_Venus_RNAP)
    W_Venus_sigma70 = exp(-params.dG_Venus_sigma70)
    W_Venus_GntR = exp(-params.dG_Venus_GntR)

    # sigma_70 has a constitutive promoter with fixed W
    W_sigma70_RNAP = exp(-0.5)  # fixed, not estimated (sigma_70 is exogenous)

    # Hill functions for sigma_70 activation
    b_GntR_sigma70 = (protein_sigma70^params.n_GntR_sigma70) /
        (params.K_GntR_sigma70^params.n_GntR_sigma70 + protein_sigma70^params.n_GntR_sigma70 + eps())

    b_Venus_sigma70 = (protein_sigma70^params.n_Venus_sigma70) /
        (params.K_Venus_sigma70^params.n_Venus_sigma70 + protein_sigma70^params.n_Venus_sigma70 + eps())

    # Hill function for GntR repression of Venus
    b_Venus_GntR = (effective_GntR^params.n_Venus_GntR) /
        (params.K_Venus_GntR^params.n_Venus_GntR + effective_GntR^params.n_Venus_GntR + eps())

    # Control functions (Boltzmann partition function)
    # GntR: activated by sigma_70
    u_GntR = (W_GntR_RNAP + W_GntR_sigma70 * b_GntR_sigma70) /
             (1.0 + W_GntR_RNAP + W_GntR_sigma70 * b_GntR_sigma70)

    # Venus: activated by sigma_70, repressed by GntR
    u_Venus = (W_Venus_RNAP + W_Venus_sigma70 * b_Venus_sigma70) /
              (1.0 + W_Venus_RNAP + W_Venus_sigma70 * b_Venus_sigma70 + W_Venus_GntR * b_Venus_GntR)

    # sigma_70: constitutive (RNAP only, no regulation)
    u_sigma70 = W_sigma70_RNAP / (1.0 + W_sigma70_RNAP)

    return [u_GntR, u_Venus, u_sigma70]
end

"""
    machinery_allocation(x, params, bio, genes) → (R_X_free, R_L_free, f_X, f_L)

Compute free machinery concentrations and allocation fractions.

Layer 2 of the resource model: RNAP and ribosomes are shared across genes/mRNAs.
R_X_free is constant (genes are constant); R_L_free is dynamic (mRNA changes).
"""
function machinery_allocation(x::AbstractVector, params::ModelParameters,
                               bio::BiophysicalConstants, genes::GeneInfo)

    v_X = bio.transcription_elongation_rate  # nt/s
    v_L = bio.translation_elongation_rate    # aa/s
    KX = bio.KX
    KL = params.KL  # estimated parameter, not the default from bio

    # Transcription allocation fractions (constant since genes are constant)
    # tau_X_j = (kE/kI) * time_constant_modifier
    # f_X_j = G_j / (tau_j * KX + (1 + tau_j) * G_j)
    tau_mRNA = [params.tau_mRNA_GntR, params.tau_mRNA_Venus, 1.0]  # sigma_70 tau fixed at 1
    f_X = zeros(N_GENES)
    for j in 1:N_GENES
        kE = v_X / genes.coding_length_nt[j]
        kI = 1.0 / bio.characteristic_initiation_time_transcription
        tau_j = (kE / kI) * tau_mRNA[j]
        G_j = x[j]  # gene concentration
        f_X[j] = G_j / (tau_j * KX + (1.0 + tau_j) * G_j + eps())
    end

    # Translation allocation fractions (dynamic — depends on mRNA levels)
    tau_protein = [params.tau_protein_GntR, params.tau_protein_Venus, 1.0]
    f_L = zeros(N_GENES)
    for j in 1:N_GENES
        kE = v_L / genes.protein_length_aa[j]
        kI = 1.0 / bio.characteristic_initiation_time_translation
        tau_j = (kE / kI) * tau_protein[j]
        m_j = x[3 + j]  # mRNA concentration
        f_L[j] = m_j / (tau_j * KL + (1.0 + tau_j) * m_j + eps())
    end

    # Free machinery (closed-form competitive allocation)
    R_X_free = bio.RNAP_concentration / (1.0 + sum(f_X))
    R_L_free = bio.ribosome_concentration / (1.0 + sum(f_L))

    return R_X_free, R_L_free, f_X, f_L
end

"""
    kinetic_rates(f_X, f_L, R_X_free, R_L_free, bio, genes) → (r_X, r_L)

Compute transcription and translation kinetic rates for each gene.
Units: μM/hr
"""
function kinetic_rates(f_X::Vector{Float64}, f_L::Vector{Float64},
                       R_X_free::Float64, R_L_free::Float64,
                       bio::BiophysicalConstants, genes::GeneInfo)

    r_X = zeros(N_GENES)
    r_L = zeros(N_GENES)

    for j in 1:N_GENES
        # Transcription rate: r_X_j = R_X_free * (v_X / l_G_j) * f_X_j
        kE_X = bio.transcription_elongation_rate / genes.coding_length_nt[j]
        r_X[j] = R_X_free * kE_X * f_X[j] * 3600.0  # convert 1/s → 1/hr, result: μM/hr

        # Translation rate: r_L_j = R_L_free * K_P * (v_L / l_P_j) * f_L_j
        # K_P is a polysome factor (set to 1.0 for now, consistent with old code)
        kE_L = bio.translation_elongation_rate / genes.protein_length_aa[j]
        r_L[j] = R_L_free * kE_L * f_L[j] * 3600.0  # μM/hr
    end

    return r_X, r_L
end

"""
    balances!(dx, x, p, t)

In-place ODE right-hand side for DifferentialEquations.jl.

Parameter tuple p = (params::ModelParameters, bio::BiophysicalConstants,
                     genes::GeneInfo, gluconate_conc::Float64)
"""
function balances!(dx::AbstractVector, x::AbstractVector, p::Tuple, t::Float64)

    params, bio, genes, gluconate_conc = p

    # --- Transcription control ---
    u_bar = transcription_control(x, params, gluconate_conc)

    # --- Machinery allocation and kinetic rates ---
    R_X_free, R_L_free, f_X, f_L = machinery_allocation(x, params, bio, genes)
    r_X, r_L = kinetic_rates(f_X, f_L, R_X_free, R_L_free, bio, genes)

    # --- Resource fractions ---
    epsilon_X = max(x[IDX_EPSILON_X], 0.0)
    epsilon_L = max(x[IDX_EPSILON_L], 0.0)

    # --- Degradation rate constants ---
    # θ_m = deg_modifier * ln(2) / mRNA_half_life
    # θ_p = deg_modifier * ln(2) / protein_half_life
    ln2 = log(2.0)
    deg_m = [params.deg_mRNA_GntR, params.deg_mRNA_Venus, 1.0]
    deg_p = [params.deg_protein_GntR, params.deg_protein_Venus, params.deg_protein_sigma70]
    theta_m = deg_m .* (ln2 / bio.mRNA_half_life_hr)
    theta_p = deg_p .* (ln2 / bio.protein_half_life_hr)

    # --- Gene balances (constant) ---
    dx[1] = 0.0
    dx[2] = 0.0
    dx[3] = 0.0

    # --- mRNA balances ---
    # dm_j/dt = r_{X,j} * ū_j * ε_X − θ_{m,j} * m_j
    for j in 1:N_GENES
        dx[3 + j] = r_X[j] * u_bar[j] * epsilon_X - theta_m[j] * x[3 + j]
    end

    # --- Protein balances ---
    # dp_j/dt = r_{L,j} * w̄_j * ε_L − θ_{p,j} * p_j
    # w̄_j = 1.0 (no translational regulation)
    for j in 1:N_GENES
        dx[6 + j] = r_L[j] * 1.0 * epsilon_L - theta_p[j] * x[6 + j]
    end

    # --- Consumable resource balances ---
    # dε_X/dt = −α_X * Σ_j(l_{G,j} * r_{X,j} * ū_j) * ε_X
    tx_consumption = 0.0
    for j in 1:N_GENES
        tx_consumption += genes.coding_length_nt[j] * r_X[j] * u_bar[j]
    end
    dx[IDX_EPSILON_X] = -params.alpha_X * tx_consumption * epsilon_X

    # dε_L/dt = −α_L * Σ_j(l_{P,j} * r_{L,j} * w̄_j) * ε_L
    tl_consumption = 0.0
    for j in 1:N_GENES
        tl_consumption += genes.protein_length_aa[j] * r_L[j] * 1.0
    end
    dx[IDX_EPSILON_L] = -params.alpha_L * tl_consumption * epsilon_L

    return nothing
end
