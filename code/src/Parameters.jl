# Parameters.jl - Single source of truth for parameter vector mapping and bounds

const N_PARAMETERS = 25

const PARAMETER_NAMES = [
    "dG_GntR_RNAP",        # 1
    "dG_GntR_sigma70",     # 2
    "dG_Venus_RNAP",       # 3
    "dG_Venus_sigma70",    # 4
    "dG_Venus_GntR",       # 5
    "n_GntR_sigma70",      # 6
    "K_GntR_sigma70",      # 7
    "n_Venus_sigma70",     # 8
    "K_Venus_sigma70",     # 9
    "n_Venus_GntR",        # 10
    "K_Venus_GntR",        # 11
    "tau_mRNA_GntR",       # 12
    "tau_mRNA_Venus",      # 13
    "tau_protein_GntR",    # 14
    "tau_protein_Venus",   # 15
    "deg_mRNA_GntR",       # 16
    "deg_mRNA_Venus",      # 17
    "deg_protein_GntR",    # 18
    "deg_protein_Venus",   # 19
    "deg_protein_sigma70", # 20
    "KL",                  # 21
    "n_gluconate",         # 22
    "K_gluconate",         # 23
    "alpha_X",             # 24
    "alpha_L",             # 25
]

# Parameter bounds: [lower upper] for each parameter
const PARAMETER_BOUNDS = [
    #  lower    upper
     0.01      5.0    ;  # 1  dG_GntR_RNAP
    -5.0      -0.01   ;  # 2  dG_GntR_sigma70
     0.01      5.0    ;  # 3  dG_Venus_RNAP
    -5.0      -0.01   ;  # 4  dG_Venus_sigma70
    -5.0      -0.1    ;  # 5  dG_Venus_GntR
     0.5      10.0    ;  # 6  n_GntR_sigma70
     0.001   100.0    ;  # 7  K_GntR_sigma70
     0.5      10.0    ;  # 8  n_Venus_sigma70
     0.001   100.0    ;  # 9  K_Venus_sigma70
     0.5      10.0    ;  # 10 n_Venus_GntR
     0.001   100.0    ;  # 11 K_Venus_GntR
     0.001   100.0    ;  # 12 tau_mRNA_GntR
     0.001   100.0    ;  # 13 tau_mRNA_Venus
     0.001   100.0    ;  # 14 tau_protein_GntR
     0.001   100.0    ;  # 15 tau_protein_Venus
     0.001   100.0    ;  # 16 deg_mRNA_GntR
     0.001   100.0    ;  # 17 deg_mRNA_Venus
     0.001     1.0    ;  # 18 deg_protein_GntR
     0.001     1.0    ;  # 19 deg_protein_Venus
     0.001     1.0    ;  # 20 deg_protein_sigma70
    10.0    1000.0    ;  # 21 KL (μM)
     1.0       5.0    ;  # 22 n_gluconate
     0.1     100.0    ;  # 23 K_gluconate (mM)
     1e-6      1e-1   ;  # 24 alpha_X (TX resource depletion)
     1e-6      1e-1   ;  # 25 alpha_L (TL resource depletion)
]

const PARAMETER_LOWER = PARAMETER_BOUNDS[:, 1]
const PARAMETER_UPPER = PARAMETER_BOUNDS[:, 2]

"""
    vector_to_parameters(pvec::Vector{Float64})::ModelParameters

Convert a flat parameter vector (length 25) to a ModelParameters struct.
"""
function vector_to_parameters(pvec::Vector{Float64})::ModelParameters
    return ModelParameters(
        dG_GntR_RNAP       = pvec[1],
        dG_GntR_sigma70    = pvec[2],
        dG_Venus_RNAP      = pvec[3],
        dG_Venus_sigma70   = pvec[4],
        dG_Venus_GntR      = pvec[5],
        n_GntR_sigma70     = pvec[6],
        K_GntR_sigma70     = pvec[7],
        n_Venus_sigma70    = pvec[8],
        K_Venus_sigma70    = pvec[9],
        n_Venus_GntR       = pvec[10],
        K_Venus_GntR       = pvec[11],
        tau_mRNA_GntR      = pvec[12],
        tau_mRNA_Venus     = pvec[13],
        tau_protein_GntR   = pvec[14],
        tau_protein_Venus  = pvec[15],
        deg_mRNA_GntR      = pvec[16],
        deg_mRNA_Venus     = pvec[17],
        deg_protein_GntR   = pvec[18],
        deg_protein_Venus  = pvec[19],
        deg_protein_sigma70 = pvec[20],
        KL                 = pvec[21],
        n_gluconate        = pvec[22],
        K_gluconate        = pvec[23],
        alpha_X            = pvec[24],
        alpha_L            = pvec[25],
    )
end

"""
    parameters_to_vector(p::ModelParameters)::Vector{Float64}

Convert a ModelParameters struct to a flat parameter vector.
"""
function parameters_to_vector(p::ModelParameters)::Vector{Float64}
    return [
        p.dG_GntR_RNAP, p.dG_GntR_sigma70, p.dG_Venus_RNAP, p.dG_Venus_sigma70, p.dG_Venus_GntR,
        p.n_GntR_sigma70, p.K_GntR_sigma70, p.n_Venus_sigma70, p.K_Venus_sigma70,
        p.n_Venus_GntR, p.K_Venus_GntR,
        p.tau_mRNA_GntR, p.tau_mRNA_Venus, p.tau_protein_GntR, p.tau_protein_Venus,
        p.deg_mRNA_GntR, p.deg_mRNA_Venus, p.deg_protein_GntR, p.deg_protein_Venus, p.deg_protein_sigma70,
        p.KL, p.n_gluconate, p.K_gluconate, p.alpha_X, p.alpha_L,
    ]
end

"""
    clamp_to_bounds(pvec::Vector{Float64})::Vector{Float64}

Clamp each parameter to its bounds.
"""
function clamp_to_bounds(pvec::Vector{Float64})::Vector{Float64}
    return clamp.(pvec, PARAMETER_LOWER, PARAMETER_UPPER)
end

"""
    default_initial_guess()::Vector{Float64}

Return a reasonable initial guess for the parameter vector.
"""
function default_initial_guess()::Vector{Float64}
    return [
        1.0,      # 1  dG_GntR_RNAP
       -2.5,      # 2  dG_GntR_sigma70
        1.0,      # 3  dG_Venus_RNAP
       -1.0,      # 4  dG_Venus_sigma70
       -5.0,      # 5  dG_Venus_GntR
        2.5,      # 6  n_GntR_sigma70
       30.0,      # 7  K_GntR_sigma70
        2.5,      # 8  n_Venus_sigma70
      100.0,      # 9  K_Venus_sigma70
        2.5,      # 10 n_Venus_GntR
        1.0,      # 11 K_Venus_GntR
        1.0,      # 12 tau_mRNA_GntR
        1.0,      # 13 tau_mRNA_Venus
        8.0,      # 14 tau_protein_GntR
        1.0,      # 15 tau_protein_Venus
        1.0,      # 16 deg_mRNA_GntR
        1.0,      # 17 deg_mRNA_Venus
        1.0,      # 18 deg_protein_GntR
        1.0,      # 19 deg_protein_Venus
        1.0,      # 20 deg_protein_sigma70
      100.0,      # 21 KL
        1.0,      # 22 n_gluconate
        1.0,      # 23 K_gluconate
        1e-4,     # 24 alpha_X
        1e-4,     # 25 alpha_L
    ]
end
