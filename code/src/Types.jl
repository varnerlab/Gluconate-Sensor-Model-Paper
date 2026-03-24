# Types.jl - Data structures for the gluconate biosensor model

"""
    BiophysicalConstants

Biophysical constants for the cell-free TX/TL system (from CellFree.json).
All concentrations in μM, rates in native units (nt/s or aa/s), times in seconds.
"""
struct BiophysicalConstants
    RNAP_concentration::Float64          # μM (RNA polymerase holoenzyme)
    ribosome_concentration::Float64      # μM
    transcription_elongation_rate::Float64  # nt/s
    translation_elongation_rate::Float64    # aa/s
    KX::Float64                          # transcription saturation constant (μM)
    KL_default::Float64                  # translation saturation constant (μM, default from literature)
    characteristic_initiation_time_transcription::Float64  # s
    characteristic_initiation_time_translation::Float64    # s
    mRNA_half_life_hr::Float64           # hr
    protein_half_life_hr::Float64        # hr
end

"""
    GeneInfo

Known properties of the genes in the circuit.
Order: [GntR, Venus, sigma_70]
"""
struct GeneInfo
    names::Vector{Symbol}
    coding_length_nt::Vector{Float64}    # gene length in nucleotides
    protein_length_aa::Vector{Float64}   # protein length in amino acids
    initial_abundance::Vector{Float64}   # gene concentration in μM
end

"""
    ModelParameters

All 25 estimated model parameters. Mutable to allow parameter vector unpacking.
"""
Base.@kwdef mutable struct ModelParameters
    # Pseudo-energy parameters (W = exp(-dG))
    dG_GntR_RNAP::Float64 = 1.0
    dG_GntR_sigma70::Float64 = -2.5
    dG_Venus_RNAP::Float64 = 1.0
    dG_Venus_sigma70::Float64 = -1.0
    dG_Venus_GntR::Float64 = -5.0

    # Binding parameters (Hill coefficient n, dissociation constant K)
    n_GntR_sigma70::Float64 = 2.5
    K_GntR_sigma70::Float64 = 30.0
    n_Venus_sigma70::Float64 = 2.5
    K_Venus_sigma70::Float64 = 100.0
    n_Venus_GntR::Float64 = 2.5
    K_Venus_GntR::Float64 = 1.0

    # Time constant modifiers (scale initiation time for TX/TL)
    tau_mRNA_GntR::Float64 = 1.0
    tau_mRNA_Venus::Float64 = 1.0
    tau_protein_GntR::Float64 = 8.0
    tau_protein_Venus::Float64 = 1.0

    # Degradation modifiers (scale base degradation rate)
    deg_mRNA_GntR::Float64 = 1.0
    deg_mRNA_Venus::Float64 = 1.0
    deg_protein_GntR::Float64 = 1.0
    deg_protein_Venus::Float64 = 1.0
    deg_protein_sigma70::Float64 = 1.0

    # Translation saturation constant (μM)
    KL::Float64 = 100.0

    # Gluconate-GntR binding
    n_gluconate::Float64 = 1.0
    K_gluconate::Float64 = 1.0          # mM

    # Consumable resource depletion rates (NEW)
    alpha_X::Float64 = 1e-4             # TX resource consumption rate
    alpha_L::Float64 = 1e-4             # TL resource consumption rate
end

"""
    ExperimentalData

Parsed experimental measurements for model training and validation.
"""
struct ExperimentalData
    # mRNA data (training, 10 mM gluconate + GntR condition)
    mRNA_time::Vector{Float64}           # hr
    mRNA_venus_mean::Vector{Float64}     # μM (converted from nM during loading)
    mRNA_venus_std::Vector{Float64}      # μM
    mRNA_gntr_mean::Vector{Float64}      # μM (converted from nM during loading)
    mRNA_gntr_std::Vector{Float64}       # μM

    # Protein data (training, 10 mM gluconate condition)
    protein_time::Vector{Float64}        # hr
    protein_venus_mean::Vector{Float64}  # μM (10 mM gluconate)
    protein_venus_std::Vector{Float64}   # μM

    # Protein data (training, 0 mM gluconate — full repression condition)
    protein_venus_0mM_mean::Vector{Float64}  # μM
    protein_venus_0mM_std::Vector{Float64}   # μM

    # Dose-response data (validation)
    dose_gluconate::Vector{Float64}      # mM
    dose_venus_mean::Vector{Float64}     # μM
    dose_venus_std::Vector{Float64}      # μM
end
