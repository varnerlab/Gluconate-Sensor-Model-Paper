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
    machinery_allocation(x, params, bio, genes) → (R_X_free, R_L_free, tau_X, tau_L)

Compute free machinery concentrations from the exact multi-gene RNAP/ribosome balance.

Derived from the 4-step transcription mechanism (McClure 1980, Adhikari et al. 2020):
  G_j + R_X ⇌ (G_j:R_X)_C → (G_j:R_X)_O → m_j + R_X + G_j

At steady state, the total RNAP balance gives:
  R_X,T = R_X [1 + Σ_j (1 + τ_j⁻¹) G_j / K_X]

Solving for free RNAP:
  R_X = R_X,T K_X / [K_X + Σ_j (τ_j + 1)/τ_j · G_j]

The transcription rate for gene j is then:
  r_X,j = k_{E,j} · R_X,T · G_j / {τ_j · [K_X + Σ_i (τ_i + 1)/τ_i · G_i]}

An analogous derivation holds for translation with ribosome replacing RNAP
and mRNA replacing gene concentration.
"""
function machinery_allocation(x::AbstractVector, params::ModelParameters,
                               bio::BiophysicalConstants, genes::GeneInfo)

    v_X = bio.transcription_elongation_rate  # nt/s
    v_L = bio.translation_elongation_rate    # aa/s
    KX = bio.KX
    KL = params.KL  # estimated parameter, not the default from bio

    # Compute tau values for transcription
    tau_mRNA = [params.tau_mRNA_GntR, params.tau_mRNA_Venus, 1.0]
    tau_X = zeros(N_GENES)
    for j in 1:N_GENES
        kE = v_X / genes.coding_length_nt[j]
        kI = 1.0 / bio.characteristic_initiation_time_transcription
        tau_X[j] = (kE / kI) * tau_mRNA[j]
    end

    # Compute tau values for translation
    tau_protein = [params.tau_protein_GntR, params.tau_protein_Venus, 1.0]
    tau_L = zeros(N_GENES)
    for j in 1:N_GENES
        kE = v_L / genes.protein_length_aa[j]
        kI = 1.0 / bio.characteristic_initiation_time_translation
        tau_L[j] = (kE / kI) * tau_protein[j]
    end

    # Free RNAP: R_X = R_{X,T} K_X / D_X where D_X = K_X + Σ_j (τ_j+1)/τ_j · G_j
    D_X = KX
    for j in 1:N_GENES
        G_j = x[j]
        if G_j > 0
            D_X += (tau_X[j] + 1.0) / (tau_X[j] + eps()) * G_j
        end
    end
    R_X_free = bio.RNAP_concentration * KX / D_X

    # Free ribosome: R_L = R_{L,T} K_L / D_L where D_L = K_L + Σ_j (τ_j+1)/τ_j · m_j
    D_L = KL
    for j in 1:N_GENES
        m_j = x[3 + j]
        if m_j > 0
            D_L += (tau_L[j] + 1.0) / (tau_L[j] + eps()) * m_j
        end
    end
    R_L_free = bio.ribosome_concentration * KL / D_L

    return R_X_free, R_L_free, tau_X, tau_L
end

"""
    kinetic_rates(x, tau_X, tau_L, R_X_free, R_L_free, bio, genes) → (r_X, r_L)

Compute transcription and translation kinetic rates for each gene using the exact
multi-gene expressions derived from the 4-step mechanism.

  r_{X,j} = k_{E,j} · R_X · G_j / (τ_j · K_X)    [μM/hr]
  r_{L,j} = k_{E,j} · R_L · m_j / (τ_j · K_L)    [μM/hr]

where R_X and R_L are the free machinery concentrations that already account
for multi-gene competition through the shared denominator.
"""
function kinetic_rates(x::AbstractVector, tau_X::Vector{Float64}, tau_L::Vector{Float64},
                       R_X_free::Float64, R_L_free::Float64,
                       bio::BiophysicalConstants, genes::GeneInfo,
                       params::ModelParameters)

    KX = bio.KX
    KL = params.KL

    r_X = zeros(N_GENES)
    r_L = zeros(N_GENES)

    for j in 1:N_GENES
        G_j = x[j]
        m_j = x[3 + j]

        kE_X = bio.transcription_elongation_rate / genes.coding_length_nt[j]
        r_X[j] = kE_X * R_X_free * G_j / (tau_X[j] * KX + eps()) * 3600.0

        kE_L = bio.translation_elongation_rate / genes.protein_length_aa[j]
        r_L[j] = kE_L * R_L_free * m_j / (tau_L[j] * KL + eps()) * 3600.0
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
    R_X_free, R_L_free, tau_X, tau_L = machinery_allocation(x, params, bio, genes)
    r_X, r_L = kinetic_rates(x, tau_X, tau_L, R_X_free, R_L_free, bio, genes, params)

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
