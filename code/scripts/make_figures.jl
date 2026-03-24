# make_figures.jl - Generate all publication figures
#
# Usage: julia --project=code scripts/make_figures.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "GluconateBiosensor.jl"))
using .GluconateBiosensor
using CairoMakie
using DelimitedFiles
using Statistics

const RESULTS_DIR = joinpath(@__DIR__, "..", "results")
const FIGURES_DIR = joinpath(@__DIR__, "..", "results", "figures")
const PAPER_FIGS_DIR = joinpath(@__DIR__, "..", "..", "paper", "figs")
mkpath(FIGURES_DIR)
mkpath(PAPER_FIGS_DIR)

# --- Load data ---
exp_data = load_experimental_data()
t_eval = vec(readdlm(joinpath(RESULTS_DIR, "time_eval.dat")))

mRNA_Venus = readdlm(joinpath(RESULTS_DIR, "mRNA_Venus_ensemble.dat"))
mRNA_GntR = readdlm(joinpath(RESULTS_DIR, "mRNA_GntR_ensemble.dat"))
protein_Venus = readdlm(joinpath(RESULTS_DIR, "protein_Venus_ensemble.dat"))
protein_GntR = readdlm(joinpath(RESULTS_DIR, "protein_GntR_ensemble.dat"))
epsilon_X = readdlm(joinpath(RESULTS_DIR, "epsilon_X_ensemble.dat"))
epsilon_L = readdlm(joinpath(RESULTS_DIR, "epsilon_L_ensemble.dat"))

# --- Helper: ensemble statistics ---
function ensemble_stats(M::Matrix)
    μ = vec(mean(M, dims = 2))
    σ = vec(std(M, dims = 2))
    return μ, σ
end

# --- Figure 1: Model training (4-panel: Venus mRNA, Venus protein, GntR mRNA, GntR protein) ---
begin
    fig = Figure(size = (900, 700), fontsize = 12)

    # Venus mRNA
    ax1 = Axis(fig[1, 1], xlabel = "Time (hr)", ylabel = "Venus mRNA (μM)", title = "Venus mRNA")
    μ, σ = ensemble_stats(mRNA_Venus)
    band!(ax1, t_eval, μ .- 1.96 .* σ, μ .+ 1.96 .* σ, color = (:dodgerblue, 0.3))
    lines!(ax1, t_eval, μ, color = :dodgerblue, linewidth = 2)
    errorbars!(ax1, exp_data.mRNA_time, exp_data.mRNA_venus_mean, exp_data.mRNA_venus_std,
               color = :black, whiskerwidth = 6)
    scatter!(ax1, exp_data.mRNA_time, exp_data.mRNA_venus_mean, color = :black, markersize = 8)

    # Venus protein
    ax2 = Axis(fig[1, 2], xlabel = "Time (hr)", ylabel = "Venus protein (μM)", title = "Venus protein")
    μ, σ = ensemble_stats(protein_Venus)
    band!(ax2, t_eval, μ .- 1.96 .* σ, μ .+ 1.96 .* σ, color = (:orange, 0.3))
    lines!(ax2, t_eval, μ, color = :orange, linewidth = 2)
    scatter!(ax2, exp_data.protein_time, exp_data.protein_venus_mean,
             color = :black, markersize = 3)

    # GntR mRNA
    ax3 = Axis(fig[2, 1], xlabel = "Time (hr)", ylabel = "GntR mRNA (μM)", title = "GntR mRNA")
    μ, σ = ensemble_stats(mRNA_GntR)
    band!(ax3, t_eval, μ .- 1.96 .* σ, μ .+ 1.96 .* σ, color = (:green, 0.3))
    lines!(ax3, t_eval, μ, color = :green, linewidth = 2)
    errorbars!(ax3, exp_data.mRNA_time, exp_data.mRNA_gntr_mean, exp_data.mRNA_gntr_std,
               color = :black, whiskerwidth = 6)
    scatter!(ax3, exp_data.mRNA_time, exp_data.mRNA_gntr_mean, color = :black, markersize = 8)

    # GntR protein (no experimental data)
    ax4 = Axis(fig[2, 2], xlabel = "Time (hr)", ylabel = "GntR protein (μM)", title = "GntR protein (predicted)")
    μ, σ = ensemble_stats(protein_GntR)
    band!(ax4, t_eval, μ .- 1.96 .* σ, μ .+ 1.96 .* σ, color = (:red, 0.3))
    lines!(ax4, t_eval, μ, color = :red, linewidth = 2)

    save(joinpath(FIGURES_DIR, "fig_model_training.pdf"), fig)
    save(joinpath(PAPER_FIGS_DIR, "fig_model_training.pdf"), fig)
    println("Saved: fig_model_training.pdf")
end

# --- Figure 2: Resource depletion ---
begin
    fig = Figure(size = (500, 350), fontsize = 12)
    ax = Axis(fig[1, 1], xlabel = "Time (hr)", ylabel = "Resource fraction remaining",
              title = "Consumable Resource Depletion")

    μX, σX = ensemble_stats(epsilon_X)
    μL, σL = ensemble_stats(epsilon_L)

    band!(ax, t_eval, μX .- 1.96 .* σX, μX .+ 1.96 .* σX, color = (:blue, 0.2))
    lines!(ax, t_eval, μX, color = :blue, linewidth = 2, label = "εX (NTPs)")

    band!(ax, t_eval, μL .- 1.96 .* σL, μL .+ 1.96 .* σL, color = (:red, 0.2))
    lines!(ax, t_eval, μL, color = :red, linewidth = 2, label = "εL (amino acids)")

    axislegend(ax, position = :rt)
    ylims!(ax, 0, 1.05)

    save(joinpath(FIGURES_DIR, "fig_resource_depletion.pdf"), fig)
    save(joinpath(PAPER_FIGS_DIR, "fig_resource_depletion.pdf"), fig)
    println("Saved: fig_resource_depletion.pdf")
end

# --- Figure 3: Dose-response ---
begin
    dose_gluc = vec(readdlm(joinpath(RESULTS_DIR, "dose_gluconate.dat")))
    dose_venus = readdlm(joinpath(RESULTS_DIR, "dose_venus_12h.dat"))

    fig = Figure(size = (500, 400), fontsize = 12)
    ax = Axis(fig[1, 1], xlabel = "D-Gluconate (mM)", ylabel = "Venus protein at 12h (μM)",
              title = "Dose-Response (Validation)", xscale = log10)

    # Simulated ensemble
    μ_dose = vec(mean(dose_venus, dims = 2))
    σ_dose = vec(std(dose_venus, dims = 2))
    band!(ax, dose_gluc, μ_dose .- 1.96 .* σ_dose, μ_dose .+ 1.96 .* σ_dose,
          color = (:dodgerblue, 0.3))
    lines!(ax, dose_gluc, μ_dose, color = :dodgerblue, linewidth = 2, label = "Model")

    # Experimental data
    errorbars!(ax, exp_data.dose_gluconate, exp_data.dose_venus_mean, exp_data.dose_venus_std,
               color = :black, whiskerwidth = 6)
    scatter!(ax, exp_data.dose_gluconate, exp_data.dose_venus_mean,
             color = :black, markersize = 8, label = "Experiment")

    axislegend(ax, position = :lb)

    save(joinpath(FIGURES_DIR, "fig_dose_response.pdf"), fig)
    save(joinpath(PAPER_FIGS_DIR, "fig_dose_response.pdf"), fig)
    println("Saved: fig_dose_response.pdf")
end

# --- Figure 4: Sensitivity heatmap ---
if isfile(joinpath(RESULTS_DIR, "sensitivity_means.dat"))
    means_raw = readdlm(joinpath(RESULTS_DIR, "sensitivity_means.dat"))  # 4 x 25
    output_names = ["Venus mRNA", "Venus protein", "GntR mRNA", "GntR protein"]

    # Transpose to n_params x n_outputs (25 x 4)
    means_t = collect(means_raw')

    fig = Figure(size = (700, 600), fontsize = 10)
    ax = Axis(fig[1, 1], xlabel = "Output", ylabel = "Parameter",
              title = "Morris Sensitivity (|μ*|, normalized per output)",
              xticks = (1:4, output_names),
              yticks = (1:N_PARAMETERS, PARAMETER_NAMES))

    # Normalize each column (output) to [0, 1] for visualization
    means_abs = abs.(means_t)
    for j in 1:4
        col_max = maximum(means_abs[:, j])
        if col_max > 0
            means_abs[:, j] ./= col_max
        end
    end

    # CairoMakie heatmap: data matrix is indexed as [x_idx, y_idx]
    # We want x=outputs (4), y=parameters (25), so pass a 4x25 matrix
    heatmap!(ax, 1:4, 1:N_PARAMETERS, means_abs', colormap = :viridis)

    save(joinpath(FIGURES_DIR, "fig_sensitivity.pdf"), fig)
    save(joinpath(PAPER_FIGS_DIR, "fig_sensitivity.pdf"), fig)
    println("Saved: fig_sensitivity.pdf")
else
    println("Skipping sensitivity figure (no sensitivity_means.dat)")
end

println("\nAll figures saved to $(FIGURES_DIR)")
