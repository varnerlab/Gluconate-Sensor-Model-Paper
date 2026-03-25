# make_figures.jl - Generate all publication figures from SAVED data
#
# This script does NO simulation -- it only reads data files and plots.
# Run run_dynamics.jl, run_dose_response.jl, and run_sensitivity.jl first.
#
# Usage: julia --project=code code/scripts/make_figures.jl

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

# --- Style constants (preprint matching) ---
const BAND_COLOR = (:lightblue, 0.5)
const LINE_COLOR = :black
const DATA_COLOR = :black
const LINE_WIDTH = 2.0

# --- Load experimental data ---
exp_data = load_experimental_data()

# ==========================================================================
# Figure 3: Model training (4-panel)
# Data from: run_dynamics.jl
# ==========================================================================
if isfile(joinpath(RESULTS_DIR, "time_eval.dat"))
    t_eval = vec(readdlm(joinpath(RESULTS_DIR, "time_eval.dat")))
    mRNA_Venus = readdlm(joinpath(RESULTS_DIR, "mRNA_Venus_ensemble.dat"))
    mRNA_GntR = readdlm(joinpath(RESULTS_DIR, "mRNA_GntR_ensemble.dat"))
    protein_Venus = readdlm(joinpath(RESULTS_DIR, "protein_Venus_ensemble.dat"))
    protein_GntR = readdlm(joinpath(RESULTS_DIR, "protein_GntR_ensemble.dat"))

    es(M) = (vec(mean(M, dims=2)), vec(std(M, dims=2)))

    fig = Figure(size = (900, 700), fontsize = 14, figure_padding = 20)

    # A: Venus mRNA (nM)
    ax1 = Axis(fig[1, 1], xlabel = "Time (hr)", ylabel = "[Venus mRNA] (nM)")
    μ, σ = es(mRNA_Venus .* 1000.0)
    band!(ax1, t_eval, max.(μ .- 1.96σ, 0), μ .+ 1.96σ, color = BAND_COLOR)
    lines!(ax1, t_eval, μ, color = LINE_COLOR, linewidth = LINE_WIDTH)
    errorbars!(ax1, exp_data.mRNA_time, exp_data.mRNA_venus_mean .* 1000, exp_data.mRNA_venus_std .* 1000,
               color = DATA_COLOR, whiskerwidth = 6, linewidth = 1.5)
    scatter!(ax1, exp_data.mRNA_time, exp_data.mRNA_venus_mean .* 1000, color = DATA_COLOR, markersize = 10)
    text!(ax1, 0.05, 0.95, text = "A", font = :bold, fontsize = 18, space = :relative, align = (:left, :top))

    # B: Venus protein (μM)
    ax2 = Axis(fig[1, 2], xlabel = "Time (hr)", ylabel = "[Venus] (μM)")
    μ, σ = es(protein_Venus)
    band!(ax2, t_eval, max.(μ .- 1.96σ, 0), μ .+ 1.96σ, color = BAND_COLOR)
    lines!(ax2, t_eval, μ, color = LINE_COLOR, linewidth = LINE_WIDTH)
    errorbars!(ax2, exp_data.protein_time, exp_data.protein_venus_mean, exp_data.protein_venus_std,
               color = DATA_COLOR, whiskerwidth = 4, linewidth = 1)
    scatter!(ax2, exp_data.protein_time, exp_data.protein_venus_mean, color = DATA_COLOR, markersize = 6)
    text!(ax2, 0.05, 0.95, text = "B", font = :bold, fontsize = 18, space = :relative, align = (:left, :top))

    # C: GntR mRNA (nM)
    ax3 = Axis(fig[2, 1], xlabel = "Time (hr)", ylabel = "[GntR mRNA] (nM)")
    μ, σ = es(mRNA_GntR .* 1000.0)
    band!(ax3, t_eval, max.(μ .- 1.96σ, 0), μ .+ 1.96σ, color = BAND_COLOR)
    lines!(ax3, t_eval, μ, color = LINE_COLOR, linewidth = LINE_WIDTH)
    errorbars!(ax3, exp_data.mRNA_time, exp_data.mRNA_gntr_mean .* 1000, exp_data.mRNA_gntr_std .* 1000,
               color = DATA_COLOR, whiskerwidth = 6, linewidth = 1.5)
    scatter!(ax3, exp_data.mRNA_time, exp_data.mRNA_gntr_mean .* 1000, color = DATA_COLOR, markersize = 10)
    text!(ax3, 0.05, 0.95, text = "C", font = :bold, fontsize = 18, space = :relative, align = (:left, :top))

    # D: GntR protein (μM)
    ax4 = Axis(fig[2, 2], xlabel = "Time (hr)", ylabel = "[GntR] (μM)")
    μ, σ = es(protein_GntR)
    band!(ax4, t_eval, max.(μ .- 1.96σ, 0), μ .+ 1.96σ, color = BAND_COLOR)
    lines!(ax4, t_eval, μ, color = LINE_COLOR, linewidth = LINE_WIDTH)
    text!(ax4, 0.05, 0.95, text = "D", font = :bold, fontsize = 18, space = :relative, align = (:left, :top))

    save(joinpath(FIGURES_DIR, "fig_model_training.pdf"), fig)
    save(joinpath(PAPER_FIGS_DIR, "fig_model_training.pdf"), fig)
    println("Saved: fig_model_training.pdf")
else
    println("SKIP fig_model_training (run run_dynamics.jl first)")
end

# ==========================================================================
# Figure 4: Resource depletion
# Data from: run_dynamics.jl
# ==========================================================================
if isfile(joinpath(RESULTS_DIR, "epsilon_X_ensemble.dat"))
    t_eval = vec(readdlm(joinpath(RESULTS_DIR, "time_eval.dat")))
    epsilon_X = readdlm(joinpath(RESULTS_DIR, "epsilon_X_ensemble.dat"))
    epsilon_L = readdlm(joinpath(RESULTS_DIR, "epsilon_L_ensemble.dat"))

    es(M) = (vec(mean(M, dims=2)), vec(std(M, dims=2)))

    fig = Figure(size = (500, 400), fontsize = 14, figure_padding = 20)
    ax = Axis(fig[1, 1], xlabel = "Time (hr)", ylabel = "Resource fraction remaining")

    μX, σX = es(epsilon_X)
    μL, σL = es(epsilon_L)

    band!(ax, t_eval, max.(μX .- 1.96σX, 0), min.(μX .+ 1.96σX, 1.05), color = (:lightblue, 0.4))
    lines!(ax, t_eval, μX, color = :black, linewidth = LINE_WIDTH, label = "εX (NTPs)")

    band!(ax, t_eval, max.(μL .- 1.96σL, 0), min.(μL .+ 1.96σL, 1.05), color = (:lightsalmon, 0.4))
    lines!(ax, t_eval, μL, color = :black, linewidth = LINE_WIDTH, linestyle = :dash, label = "εL (amino acids)")

    axislegend(ax, position = :lb, framevisible = false, labelsize = 12)
    ylims!(ax, 0, 1.05)

    save(joinpath(FIGURES_DIR, "fig_resource_depletion.pdf"), fig)
    save(joinpath(PAPER_FIGS_DIR, "fig_resource_depletion.pdf"), fig)
    println("Saved: fig_resource_depletion.pdf")
else
    println("SKIP fig_resource_depletion (run run_dynamics.jl first)")
end

# ==========================================================================
# Figure 5: Dose-response (dense, from saved data)
# Data from: run_dose_response.jl
# ==========================================================================
if isfile(joinpath(RESULTS_DIR, "dose_dense_log_gluconate.dat"))
    log_doses = vec(readdlm(joinpath(RESULTS_DIR, "dose_dense_log_gluconate.dat")))
    med = vec(readdlm(joinpath(RESULTS_DIR, "dose_dense_venus_median.dat")))
    q025 = vec(readdlm(joinpath(RESULTS_DIR, "dose_dense_venus_q025.dat")))
    q975 = vec(readdlm(joinpath(RESULTS_DIR, "dose_dense_venus_q975.dat")))
    log_dose_exp = log10.(exp_data.dose_gluconate)

    fig = Figure(size = (550, 450), fontsize = 14, figure_padding = 20)
    ax = Axis(fig[1, 1], xlabel = "log [Gluconate] (mM)", ylabel = "[Venus] (μM)",
              xlabelsize = 13, ylabelsize = 13, xticks = -4:1:1)

    band!(ax, log_doses, q025, q975, color = BAND_COLOR)
    lines!(ax, log_doses, med, color = LINE_COLOR, linewidth = LINE_WIDTH)
    errorbars!(ax, log_dose_exp, exp_data.dose_venus_mean, exp_data.dose_venus_std,
               color = DATA_COLOR, whiskerwidth = 6, linewidth = 1.5)
    scatter!(ax, log_dose_exp, exp_data.dose_venus_mean, color = DATA_COLOR, markersize = 10)

    ylims!(ax, 0.3, 2.5)
    xlims!(ax, -4.5, 1.5)

    save(joinpath(FIGURES_DIR, "fig_dose_response.pdf"), fig)
    save(joinpath(PAPER_FIGS_DIR, "fig_dose_response.pdf"), fig)
    println("Saved: fig_dose_response.pdf")
else
    println("SKIP fig_dose_response (run run_dose_response.jl first)")
end

# ==========================================================================
# Figure 6: Sensitivity heatmap (preprint-style binned)
# Data from: run_sensitivity.jl
# ==========================================================================
if isfile(joinpath(RESULTS_DIR, "sensitivity_means.dat"))
    means_raw = readdlm(joinpath(RESULTS_DIR, "sensitivity_means.dat"))   # 4 x 25
    vars_raw = readdlm(joinpath(RESULTS_DIR, "sensitivity_variances.dat")) # 4 x 25
    means = collect(means_raw')  # 25 x 4
    vars = collect(vars_raw')

    # 8-column matrix: Venus(μ-mRNA, σ-mRNA, μ-prot, σ-prot) + GntR(μ-mRNA, σ-mRNA, μ-prot, σ-prot)
    sensitivity_matrix = hcat(
        means[:, 1], sqrt.(abs.(vars[:, 1])),
        means[:, 2], sqrt.(abs.(vars[:, 2])),
        means[:, 3], sqrt.(abs.(vars[:, 3])),
        means[:, 4], sqrt.(abs.(vars[:, 4])),
    )

    # Bin by percentile of absolute value
    function bin_sensitivity(val, thresholds)
        av = abs(val)
        for (i, t) in enumerate(thresholds)
            if av < t; return i - 1; end
        end
        return length(thresholds)
    end

    all_abs = filter(x -> x > 0, abs.(sensitivity_matrix[:]))
    thresholds = [quantile(all_abs, p) for p in [0.25, 0.50, 0.75, 0.90]]
    binned = [bin_sensitivity(sensitivity_matrix[i, j], thresholds) for i in 1:25, j in 1:8]

    param_labels = [
        "dG GntR RNAP", "dG GntR σ70", "dG Venus RNAP", "dG Venus σ70", "dG Venus GntR",
        "n GntR σ70", "K GntR σ70", "n Venus σ70", "K Venus σ70", "n Venus GntR", "K Venus GntR",
        "τ mRNA GntR", "τ mRNA Venus", "τ prot GntR", "τ prot Venus",
        "δ mRNA GntR", "δ mRNA Venus", "δ prot GntR", "δ prot Venus", "δ prot σ70",
        "KL", "n gluconate", "K gluconate", "αX", "αL",
    ]
    col_labels = ["μ-mRNA", "σ-mRNA", "μ-prot", "σ-prot", "μ-mRNA", "σ-mRNA", "μ-prot", "σ-prot"]

    colors = [
        RGBf(1.0, 1.0, 1.0), RGBf(0.82, 0.85, 0.90), RGBf(0.60, 0.65, 0.72),
        RGBf(0.35, 0.40, 0.50), RGBf(0.10, 0.15, 0.22),
    ]

    fig = Figure(size = (500, 900), fontsize = 10)
    # Columns 1-4 at positions 1-4, columns 5-8 shifted to 5.5-8.5
    xtick_positions = [1, 2, 3, 4, 5.5, 6.5, 7.5, 8.5]
    ax = Axis(fig[1, 1], yticks = (1:25, param_labels),
              xticks = (xtick_positions, col_labels),
              xticklabelrotation = π/4, yreversed = false, xaxisposition = :top,
              aspect = DataAspect())

    # Draw cells with a gap at column 4/5 boundary (Venus | GntR separator)
    for i in 1:25, j in 1:8
        # Shift columns 5-8 right by 0.5 to create buffer around divider
        jpos = j <= 4 ? Float64(j) : j + 0.5
        poly!(ax, Rect(jpos - 0.45, i - 0.45, 0.9, 0.9),
              color = colors[binned[i, j] + 1], strokecolor = :black, strokewidth = 0.5)
    end

    text!(ax, 2.5, 26.3, text = "Venus", font = :bold, fontsize = 13, align = (:center, :bottom))
    text!(ax, 7.0, 26.3, text = "GntR", font = :bold, fontsize = 13, align = (:center, :bottom))
    vlines!(ax, [4.75], color = :black, linewidth = 1.5)
    text!(ax, 1.5, 25.8, text = "mRNA", fontsize = 10, font = :italic, align = (:center, :bottom))
    text!(ax, 3.5, 25.8, text = "protein", fontsize = 10, font = :italic, align = (:center, :bottom))
    text!(ax, 6.0, 25.8, text = "mRNA", fontsize = 10, font = :italic, align = (:center, :bottom))
    text!(ax, 8.0, 25.8, text = "protein", fontsize = 10, font = :italic, align = (:center, :bottom))

    xlims!(ax, 0.3, 9.2)
    ylims!(ax, 0.3, 27)
    hidespines!(ax)

    # Vertical discrete colorbar on right side -- squares matching the heatmap cells
    level_labels = ["none", "low", "medium", "high", "very high"]
    ax_cb = Axis(fig[1, 2],
                 yticks = ([], []), xticks = ([], []),
                 aspect = DataAspect())
    hidespines!(ax_cb)
    hideydecorations!(ax_cb)
    hidexdecorations!(ax_cb)

    # Draw 5 squares vertically, spaced to match heatmap row spacing
    for (k, (c, lab)) in enumerate(zip(colors, level_labels))
        ypos = Float64(k)
        poly!(ax_cb, Rect(-0.45, ypos - 0.45, 0.9, 0.9), color = c, strokecolor = :black, strokewidth = 0.5)
        text!(ax_cb, 0.7, ypos, text = lab, fontsize = 9, align = (:left, :center))
    end
    text!(ax_cb, 0.5, 6.0, text = "Influence", fontsize = 10, font = :bold, align = (:center, :bottom))
    xlims!(ax_cb, -0.6, 3.5)
    ylims!(ax_cb, 0, 6.5)

    colsize!(fig.layout, 2, 100)

    save(joinpath(FIGURES_DIR, "fig_sensitivity.pdf"), fig)
    save(joinpath(PAPER_FIGS_DIR, "fig_sensitivity.pdf"), fig)
    println("Saved: fig_sensitivity.pdf")
else
    println("SKIP fig_sensitivity (run run_sensitivity.jl first)")
end

println("\nAll figures saved to $(FIGURES_DIR) and $(PAPER_FIGS_DIR)")
