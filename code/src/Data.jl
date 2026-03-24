# Data.jl - Load experimental data from CSV files

using CSV, DataFrames

"""
    load_experimental_data(data_dir::String)::ExperimentalData

Load and parse all experimental data files.
Training data: 10 mM gluconate condition (mRNA + protein time courses).
Validation data: dose-response at 12h.
"""
function load_experimental_data(data_dir::String = joinpath(@__DIR__, "..", "data"))::ExperimentalData

    # --- mRNA data ---
    # Columns: time, Venus(+GntR+Gluc), Venus(+GntR-Gluc), Venus(-GntR-Gluc),
    #          SD1, SD2, SD3, GntR1, GntR2, GntR3, SD_1, SD_2, SD_3
    # We use the +GntR+Gluconate condition (training): cols 2 (Venus), 5 (SD), 8 (GntR), 11 (SD)
    mRNA_df = CSV.read(joinpath(data_dir, "mRNA_data.csv"), DataFrame)
    mRNA_time = Float64.(mRNA_df[!, 1])
    # Convert mRNA data from nM to μM to match model units
    mRNA_venus_mean = Float64.(mRNA_df[!, 2]) ./ 1000.0   # nM → μM, +GntR +Gluconate
    mRNA_venus_std = Float64.(mRNA_df[!, 5]) ./ 1000.0     # SD
    mRNA_gntr_mean = Float64.(mRNA_df[!, 8]) ./ 1000.0     # nM → μM, +GntR +Gluconate
    mRNA_gntr_std = Float64.(mRNA_df[!, 11]) ./ 1000.0     # SD

    # --- Protein data ---
    # Column 1: time(h)
    # Columns 5,6,7: Mean_10mM, StDev_10mM, StdErr_10mM (training condition)
    prot_df = CSV.read(joinpath(data_dir, "protein_data.csv"), DataFrame)
    protein_time = Float64.(prot_df[!, 1])
    protein_venus_mean = Float64.(prot_df[!, 5])   # μM, 10 mM gluconate
    protein_venus_std = Float64.(prot_df[!, 6])     # StDev
    protein_venus_0mM_mean = Float64.(prot_df[!, 29])  # μM, 0 mM gluconate (full repression)
    protein_venus_0mM_std = Float64.(prot_df[!, 30])   # StDev

    # --- Dose-response data ---
    dose_df = CSV.read(joinpath(data_dir, "dose_response.csv"), DataFrame)
    dose_gluconate = Float64.(dose_df[!, 1])      # mM
    dose_venus_mean = Float64.(dose_df[!, 2])     # μM
    dose_venus_std = Float64.(dose_df[!, 3])      # StdErr

    return ExperimentalData(
        mRNA_time, mRNA_venus_mean, mRNA_venus_std, mRNA_gntr_mean, mRNA_gntr_std,
        protein_time, protein_venus_mean, protein_venus_std,
        protein_venus_0mM_mean, protein_venus_0mM_std,
        dose_gluconate, dose_venus_mean, dose_venus_std,
    )
end
