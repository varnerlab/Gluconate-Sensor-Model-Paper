# analyze_ensemble.jl - Characterize parameter distributions and correlations
#
# Usage: julia --project=code code/scripts/analyze_ensemble.jl
#
# Generates:
#   - Parameter distribution summary (mean, std, CV, median, IQR, skewness)
#   - Pairwise Spearman correlation matrix
#   - Identification of strongly correlated parameter pairs
#   - Figures: parameter histograms, correlation heatmap, scatter plots of top correlated pairs

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "GluconateBiosensor.jl"))
using .GluconateBiosensor
using CairoMakie
using DelimitedFiles
using Statistics
using LinearAlgebra

const RESULTS_DIR = joinpath(@__DIR__, "..", "results")
const FIGURES_DIR = joinpath(@__DIR__, "..", "results", "figures")
const PAPER_FIGS_DIR = joinpath(@__DIR__, "..", "..", "paper", "figs")
mkpath(FIGURES_DIR)
mkpath(PAPER_FIGS_DIR)

# --- Load ensemble ---
PC = readdlm(joinpath(RESULTS_DIR, "PC_final.dat"))
n_params, n_ensemble = size(PC)
println("Loaded ensemble: $(n_ensemble) parameter sets, $(n_params) parameters")

# ============================================================
# 1. Parameter Distribution Summary
# ============================================================
println("\n" * "="^80)
println("PARAMETER DISTRIBUTION SUMMARY")
println("="^80)
println(rpad("Parameter", 24), rpad("Mean", 14), rpad("Std", 14), rpad("CV%", 10),
        rpad("Median", 14), rpad("IQR", 14), rpad("Skew", 10))
println("-"^100)

function skewness(x)
    n = length(x)
    m = mean(x)
    s = std(x)
    s == 0 && return 0.0
    return (n / ((n-1)*(n-2))) * sum(((xi - m) / s)^3 for xi in x)
end

for i in 1:n_params
    vals = PC[i, :]
    μ = mean(vals)
    σ = std(vals)
    cv = σ / (abs(μ) + eps()) * 100
    med = median(vals)
    q25 = quantile(vals, 0.25)
    q75 = quantile(vals, 0.75)
    iqr = q75 - q25
    sk = skewness(vals)
    println(rpad(PARAMETER_NAMES[i], 24),
            rpad(round(μ, sigdigits=4), 14),
            rpad(round(σ, sigdigits=4), 14),
            rpad(round(cv, digits=1), 10),
            rpad(round(med, sigdigits=4), 14),
            rpad(round(iqr, sigdigits=4), 14),
            round(sk, digits=2))
end

# ============================================================
# 2. Spearman Rank Correlation Matrix
# ============================================================
println("\n" * "="^80)
println("SPEARMAN RANK CORRELATIONS")
println("="^80)

function spearman_corr(x, y)
    rx = sortperm(sortperm(x))
    ry = sortperm(sortperm(y))
    n = length(x)
    d = Float64.(rx .- ry)
    return 1.0 - 6.0 * sum(d .^ 2) / (n * (n^2 - 1))
end

corr_matrix = zeros(n_params, n_params)
for i in 1:n_params
    for j in 1:n_params
        if i == j
            corr_matrix[i, j] = 1.0
        else
            corr_matrix[i, j] = spearman_corr(PC[i, :], PC[j, :])
        end
    end
end

# Find strongly correlated pairs (|rho| > 0.4)
println("\nStrongly correlated parameter pairs (|ρ| > 0.4):")
println(rpad("Parameter 1", 24), rpad("Parameter 2", 24), "ρ")
println("-"^60)

pairs_printed = Set{Tuple{Int,Int}}()
for threshold in [0.7, 0.5, 0.4]
    for i in 1:n_params
        for j in (i+1):n_params
            if abs(corr_matrix[i, j]) > threshold && !((i,j) in pairs_printed)
                push!(pairs_printed, (i, j))
                println(rpad(PARAMETER_NAMES[i], 24),
                        rpad(PARAMETER_NAMES[j], 24),
                        round(corr_matrix[i, j], digits=3))
            end
        end
    end
end

if isempty(pairs_printed)
    println("  No pairs with |ρ| > 0.4 found.")
end

# Save correlation matrix
writedlm(joinpath(RESULTS_DIR, "parameter_correlation_matrix.dat"), corr_matrix)

# ============================================================
# 3. Figures
# ============================================================

# --- Figure: Parameter histograms (5x5 grid) ---
begin
    fig = Figure(size = (1200, 1000), fontsize = 9)
    for i in 1:n_params
        row = div(i - 1, 5) + 1
        col = mod(i - 1, 5) + 1
        ax = Axis(fig[row, col], title = PARAMETER_NAMES[i],
                  titlesize = 9)
        vals = PC[i, :]
        hist!(ax, vals, bins = 30, color = (:steelblue, 0.7))
        # Mark mean
        vlines!(ax, [mean(vals)], color = :red, linewidth = 1.5)
        hideydecorations!(ax)
    end
    save(joinpath(FIGURES_DIR, "fig_parameter_distributions.pdf"), fig)
    save(joinpath(PAPER_FIGS_DIR, "fig_parameter_distributions.pdf"), fig)
    println("\nSaved: fig_parameter_distributions.pdf")
end

# --- Figure: Correlation heatmap ---
begin
    fig = Figure(size = (800, 700), fontsize = 8)
    ax = Axis(fig[1, 1],
              xticks = (1:n_params, PARAMETER_NAMES),
              yticks = (1:n_params, PARAMETER_NAMES),
              xticklabelrotation = π/3,
              title = "Spearman Rank Correlation")

    hm = heatmap!(ax, 1:n_params, 1:n_params, corr_matrix,
                   colormap = :RdBu, colorrange = (-1, 1))
    Colorbar(fig[1, 2], hm, label = "ρ")

    save(joinpath(FIGURES_DIR, "fig_correlation_heatmap.pdf"), fig)
    save(joinpath(PAPER_FIGS_DIR, "fig_correlation_heatmap.pdf"), fig)
    println("Saved: fig_correlation_heatmap.pdf")
end

# --- Figure: Scatter plots of top correlated pairs ---
begin
    # Collect top pairs
    top_pairs = Tuple{Int,Int,Float64}[]
    for i in 1:n_params
        for j in (i+1):n_params
            push!(top_pairs, (i, j, abs(corr_matrix[i, j])))
        end
    end
    sort!(top_pairs, by = x -> x[3], rev = true)

    n_plots = min(9, length(top_pairs))
    if n_plots > 0
        fig = Figure(size = (900, 900), fontsize = 9)
        for k in 1:n_plots
            i, j, rho = top_pairs[k]
            row = div(k - 1, 3) + 1
            col = mod(k - 1, 3) + 1
            ax = Axis(fig[row, col],
                      xlabel = PARAMETER_NAMES[i],
                      ylabel = PARAMETER_NAMES[j],
                      title = "ρ = $(round(corr_matrix[i,j], digits=3))",
                      titlesize = 10)
            scatter!(ax, PC[i, :], PC[j, :], markersize = 2, color = (:black, 0.3))
        end
        save(joinpath(FIGURES_DIR, "fig_parameter_scatter.pdf"), fig)
        save(joinpath(PAPER_FIGS_DIR, "fig_parameter_scatter.pdf"), fig)
        println("Saved: fig_parameter_scatter.pdf")
    end
end

println("\nDone. All ensemble analysis saved.")
