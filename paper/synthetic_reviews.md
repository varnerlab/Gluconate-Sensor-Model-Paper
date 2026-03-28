# Synthetic Peer Reviews (Updated March 2026)

Manuscript: "Mechanistic Resource Competition Modeling of a Cell-Free Gluconate Biosensor"
Journal: Frontiers in Bioengineering and Biotechnology

---

## Reviewer 1 (Cell-free systems experimentalist)

**Overall assessment:** The manuscript presents a mechanistic two-layer resource competition model for a cell-free gluconate biosensor. The framework is well-motivated, the fed-batch validation is a nice addition, and the connection to the Li et al. and Jurado et al. literature strengthens the resource bottleneck claims. However, several experimental concerns remain.

**Major concerns:**

1. **GntR protein remains unvalidated.** The model predicts GntR at 7.7 +/- 2.0 uM (CV=26%), but this has never been measured. The Discussion provides a back-of-envelope scaling argument, but this is not validation. The authors should at minimum note that Western blot, mass spectrometry, or a tagged GntR construct could resolve this. The relatively high CV (26%) means the prediction carries substantial uncertainty.

2. **Ensemble size is small (N=78).** The authors started with 73,747 candidate solutions but filtered down to 78. While the filtering criteria are reasonable (mRNA CV, GntR bounds, dose-response quality, K/n consistency), this aggressive filtering raises the question: are 78 solutions sufficient to characterize parameter uncertainty? The 95% confidence bands on the model predictions are derived from only 78 trajectories. The authors should discuss whether this is adequate or whether larger ensembles (e.g., seeded re-estimation as described in their Methods) would change the conclusions.

3. **Venus mRNA shape is not captured.** The model predicts Venus mRNA peaking at t~0.5h and declining monotonically, while the experimental data shows a rise to peak at 6h. The authors attribute this to exogenous sigma70, which is a reasonable explanation, but this means the model fundamentally cannot capture the transcriptional dynamics it claims to describe. The mRNA fit should be discussed more critically — the 95% band may encompass the data points, but the mean trajectory is qualitatively wrong in shape.

4. **Venus protein is ~10% high.** The model predicts 1.87 uM at 12h vs experimental ~1.7 uM. While within the 95% band, this systematic overshoot should be acknowledged. Combined with the dose-response baseline being ~0.8 uM (vs data ~0.55 uM), the model appears to compress the dynamic range.

**Minor concerns:**

5. The three-panel experimental data figure (Fig. 2) now uses a horizontal layout which is cleaner, but the protein data in panel A still shows considerable scatter at early time points (0-2h). Were these measurements taken during the lag phase before fluorescence develops? This could affect the protein SSE if these noisy early points dominate the objective.

6. The notation switches between epsilon (resource fractions) and alpha (consumption rates) without a clear conceptual distinction upfront. A brief sentence early in the Methods noting "epsilon tracks how much resource remains, alpha governs how fast it is consumed" would help readability.

**Recommendation:** Minor revision. The framework is sound, the fed-batch analysis is novel, and the paper is well-written. Address the GntR validation concern, acknowledge the Venus mRNA shape limitation more prominently, and discuss ensemble size adequacy.

---

## Reviewer 2 (Mathematical modeler / parameter estimation)

**Overall assessment:** This is a strong modeling paper with a clean framework and an interesting multi-objective estimation approach. The supplement provides thorough derivations, and the Pareto front visualization (Fig. S2) is a valuable addition. Several methodological concerns should be addressed.

**Major concerns:**

1. **Post-filtering uses validation data.** The ensemble was filtered based on dose-response chi-square (top 50%), which includes gluconate concentrations that were not part of the training set (0.1, 0.5, 1, 5 mM). This means the "predicted" dose-response is no longer a true blind prediction — it was used to select the ensemble. The authors should be transparent about this. One solution: report both the unfiltered and DR-filtered dose-response predictions, showing that the sigmoid shape and approximate K_gluconate are captured even without dose-response filtering.

2. **The competitive allocation approximation vs. exact multi-gene kinetics.** The supplement correctly derives the single-gene rate expression (Eqn S17) from the 4-step mechanism, but the multi-gene implementation (Eqn S22, the f_j/R_free factoring) is presented without discussing its relationship to the exact multi-gene result. For N genes, the exact rate involves cross-competition terms O_{X,j} (Adhikari et al. 2020, Eqn 3). The code uses a competitive allocation approximation instead. Have the authors verified that this approximation is adequate for their system? A brief comparison of the approximate vs. exact rates at the estimated parameter values would strengthen confidence in the formulation.

3. **Parameter identifiability beyond distributions.** The parameter distribution figure (Fig. S3) is useful but provides only marginal information. Several parameters show strong pairwise correlations (the tau/K_L cluster, the dG/tau_mRNA compensations). A scatter plot of the most correlated parameter pairs, or a subset of the 2D marginals, would help the reader understand which parameter *combinations* are identified vs. which individual parameters are sloppy. The correlation heatmap (Fig. 7) partially addresses this but doesn't show the 2D structure.

4. **Sensitivity analysis uses mean +/- 3 SD parameter ranges.** With N=78 and some parameters having non-Gaussian distributions (skewed, bound-hitting), the mean +/- 3 SD range may extend well beyond the actual ensemble support. Did the authors verify that the Morris trajectories remained in physically meaningful regions? Were any trajectories infeasible (e.g., producing negative concentrations or solver failures)?

**Minor concerns:**

5. The MO/SO cycling strategy is described well, but no convergence diagnostic is provided. Do the objective values plateau by cycle 15, or would additional cycles improve the ensemble? A simple plot of best total error vs. cycle number would answer this.

6. The fed-batch supplementation analysis (Table 2) is a nice in silico experiment, but the connection to experimental validation is qualitative ("consistent with Li et al."). Could the authors quantify this comparison? For example, Li et al. reported ~2-fold improvement with CP+Mg supplementation; the model predicts 52% increase with full epsilon_X replenishment. Is this discrepancy expected given the different systems (PURE vs. PURExpress)?

7. K_P (polysome factor) is mentioned in Eqn 6 but its value is not given. Is it set to 1? If so, state this explicitly.

**Recommendation:** Minor revision. The framework and estimation are solid. Address the filtering/prediction concern transparently, clarify the competitive allocation approximation, and provide convergence diagnostics.

---

## Reviewer 3 (Synthetic biology / biosensor applications)

**Overall assessment:** The paper has improved significantly from its earlier form, with better literature contextualization, the fed-batch validation, and proper biosensor performance metrics. However, the practical impact remains limited by the single-circuit demonstration.

**Major concerns:**

1. **Generalizability not demonstrated.** The paper's strongest claim — that the framework provides "quantitative design rules for allocating transcriptional and translational capacity" — requires demonstration on more than one circuit. The two-layer formulation is general, but has it been tested on circuits with different gene numbers, different promoter strengths, or different repressor systems? Even a simple prediction (e.g., what would happen with 5 nM vs. 10 nM GntR DNA) validated against data would strengthen the generalizability claim.

2. **76% transcriptional resource depletion is very high.** The model predicts that only 24% of transcriptional resources remain at 12h. This is a strong claim with important design implications. However, the experimental validation is indirect — the mRNA decline is consistent with resource depletion but could also reflect other mechanisms (e.g., RNase accumulation, template degradation, energy cofactor depletion affecting transcription fidelity). The authors discuss the connection to Li et al. and Jurado et al. but acknowledge their alpha parameters are lumped. Could the authors discuss what independent experimental measurements could disambiguate resource depletion from other mechanisms?

3. **No comparison to simpler models.** The paper would be strengthened by showing that the two-layer model outperforms simpler alternatives. For example: (a) a model with only consumable depletion (no machinery allocation), (b) a model with the traditional exponential decay, or (c) a model with logistic resource depletion. An AIC/BIC comparison or a prediction accuracy comparison would demonstrate that the added complexity of the two-layer formulation is justified.

**Minor concerns:**

4. The Discussion paragraph on the Li et al. fed-batch comparison is well-written but could note that the model predicts a much larger effect (52% increase) than Li et al. observed (~2-fold at 90 minutes). This difference may reflect the longer time scale of the PURExpress reactions (12h vs. 90 min) and the greater extent of resource depletion.

5. The title is accurate but could be more engaging. Consider: "A mechanistic resource competition framework reveals transcriptional capacity as the primary bottleneck in cell-free biosensor circuits."

6. The dose-response figure (Fig. 5) shows the model slightly overshooting the upper plateau (~1.9 uM vs. data ~1.7 uM). This is consistent with the Venus protein being 10% high in the training condition. The text notes K_gluconate = 0.13 mM, which is lower than the ~1 mM K_D reported by Daddaoua et al. for GntR-gluconate binding. Is this discrepancy discussed?

7. Reference 42 (Li et al. 2017) — please verify the author list. The Church lab has multiple "Li" first authors.

**Recommendation:** Minor revision. The framework is well-developed and the paper is well-written. Address the generalizability concern (even through discussion), provide a model comparison, and clarify the K_gluconate discrepancy.

---

## Summary of Key Issues Across Reviewers

| Issue | R1 | R2 | R3 | Severity | Addressed? |
|-------|----|----|-----|----------|------------|
| GntR protein unvalidated | X | | | Major | Acknowledged in text, no experiment possible |
| Post-filtering uses DR data | | X | | Major | Need transparency |
| Venus mRNA shape wrong | X | | | Major | Acknowledged (σ70 limitation) |
| Ensemble size N=78 | X | | | Moderate | Could discuss adequacy |
| Competitive allocation approximation | | X | | Moderate | Supplement derives it, needs comparison |
| Generalizability (single circuit) | | | X | Major | Discussion only |
| No model comparison/ablation | | | X | Major | Not addressed |
| Convergence diagnostics | | X | | Minor | Not shown |
| K_gluconate vs literature K_D | | | X | Minor | Not discussed |
| Venus protein 10% high | X | | X | Minor | Not explicitly noted |
| K_P value not stated | | X | | Minor | Easy fix |

## Pre-submission Actions (Recommended)

**Must address before submission:**
- [ ] Add a sentence stating K_P = 1 in the Methods
- [ ] Note in the Discussion that the dose-response ensemble was filtered using dose-response data and is therefore not a fully blind prediction (but show the sigmoid shape emerges even without this filter)
- [ ] Fix the duplicate "Second" in the limitations paragraph

**Should address if possible:**
- [ ] Add a sentence about N=78 adequacy in the Discussion
- [ ] Discuss K_gluconate (0.13 mM) vs literature K_D (~1 mM from Daddaoua et al.)
- [ ] Consider the suggested title revision
- [ ] Note Venus protein systematic overshoot (1.87 vs 1.7 uM)

**Can defer to revision:**
- [ ] Model ablation/comparison
- [ ] Generalizability demonstration
- [ ] Convergence diagnostics
- [ ] 2D parameter marginals
