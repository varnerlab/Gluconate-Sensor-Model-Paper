# Synthetic Peer Reviews

Manuscript: "Mechanistic Resource Competition Modeling of a Cell-Free Gluconate Biosensor"
Journal: Frontiers in Bioengineering and Biotechnology

---

## Reviewer 1 (Cell-free systems experimentalist)

**Overall assessment:** The manuscript presents a mathematical model for a gluconate-responsive biosensor in a reconstituted cell-free system. The resource competition framework is an interesting contribution. However, several experimental and presentation issues should be addressed before publication.

**Major concerns:**

1. **Lack of GntR protein validation.** The model predicts GntR protein at ~10 uM, but this was never measured. The authors acknowledge this limitation, but it significantly weakens the claim that the resource competition model "naturally bounds" GntR protein. Without experimental validation, the predicted GntR level is an assertion, not a result. Could the authors at minimum compare their predicted GntR concentration to what would be expected based on the known expression levels from P70 in PURExpress from the literature? Is 10 uM reasonable given 10 nM DNA template and 12 hours of expression?

2. **The 0 mM gluconate condition was used for both training and validation.** The authors state they trained on 10 mM and 0 mM conditions and then validated on the dose-response. However, the 0 mM condition is effectively one endpoint of the dose-response curve. This means the model was anchored at both extremes (0 mM and 10 mM), and the "prediction" is really interpolation between two known points. The authors should be more transparent about this. True blind prediction would involve training on only the 10 mM condition (or some intermediate condition) and predicting both endpoints.

3. **Missing time-course data at intermediate gluconate concentrations.** The protein time courses at individual gluconate doses (0.0001 to 20 mM) were collected but are not shown. The old preprint included these as supplementary figures (S2-S10). The current manuscript only shows the 12h endpoint. Including the full time courses would demonstrate whether the model captures the dynamics at intermediate doses, not just the endpoints.

4. **Experimental details from the preprint appear incomplete.** The mRNA data figure (Fig. 2) references panel C for GntR mRNA, but in the three-panel figure, panel C appears to show GntR mRNA only for the repressed and de-repressed cases. The caption mentions protein levels in panel A, but panel A appears to show protein (not mRNA). Please verify the panel labels match the caption.

**Minor concerns:**

5. The sensitivity heatmap (Fig. 6) is dominated by alpha_X and alpha_L to the point that all other parameters appear as zero. A log-scale or ranking-based visualization would be more informative for comparing the relative importance of non-resource parameters.

6. Line 260: "the four empirical resource depletion parameters (onset, slope, T_{L,1/2}, and translation capacity initial condition) were removed." The previous model had these parameters, but the manuscript does not clearly enumerate what the old parameters were. A brief comparison table (old vs. new model) would help readers who are not familiar with the prior work.

7. The number of biological replicates for the mRNA measurements is not stated (only "at least three" for protein). How many replicates were used for qPCR?

8. The manuscript should cite more recent cell-free biosensor modeling work (2023-2025). The reference list is sparse for a modeling paper in this active field.

**Recommendation:** Major revision.

---

## Reviewer 2 (Mathematical modeling / systems biology)

**Overall assessment:** This is a solid modeling contribution. The two-layer resource competition framework (machinery allocation + consumable depletion) is well-motivated, mathematically clean, and represents a genuine advance over the ad hoc decay functions used previously. The parameter estimation approach using MO/SO cycling is novel and interesting. I have several concerns about model identifiability and the interpretation of results.

**Major concerns:**

1. **Parameter identifiability is a significant concern.** Table 1 shows that many parameters have coefficients of variation exceeding 100% (e.g., K_{P70,hRNAP} at 630%, K_L at 136%, tau_{L,GntR} at 78%). The correlation analysis (Fig. 7) reveals strong correlations among translation parameters (rho > 0.82 for the tau_protein/K_L cluster). This suggests that these parameters are structurally non-identifiable from the available data. The authors discuss this briefly but should perform a more rigorous identifiability analysis. At minimum, they should: (a) compute the Fisher Information Matrix or profile likelihoods for the key parameters, (b) report which parameters are practically identifiable given the data, and (c) discuss whether the poorly identified parameters affect the model's predictive capability (they may not, if the correlated combinations are well-determined).

2. **The ensemble filtering procedure introduces bias.** The authors generated 52,542 solutions through 15 MO/SO cycles and then post-filtered to 5,464 based on GntR protein being in [7, 13] uM, Venus protein in [1.3, 2.2] uM, Venus at 0.5 mM in [0.5, 1.1] uM, and Venus at 0 mM in [0.3, 0.9] uM. This is effectively using the dose-response data (specifically the 0.5 mM and 0 mM points) as additional training constraints. The filtering criteria include outcomes that are supposed to be validation targets. This compromises the claim of blind prediction. The authors should either: (a) filter based only on training-condition criteria (10 mM and 0 mM time courses, GntR regularization) and report the resulting dose-response prediction, or (b) be transparent that the ensemble was selected to match the full dose-response and acknowledge that the dose-response is therefore not a true prediction.

3. **The resource consumption parameters alpha_X and alpha_L dominate the sensitivity analysis to a degree that is concerning.** If two parameters control essentially all model outputs, the remaining 23 parameters may be underdetermined. The authors should discuss whether the model is effectively a two-parameter model (alpha_X, alpha_L) with 23 nuisance parameters that provide fine-tuning. This would not invalidate the contribution but would change the interpretation.

4. **The machinery allocation layer (Eqs. 10-13) may not contribute meaningfully.** Since gene concentrations are constant, R_{X,free} is also constant throughout the simulation. The authors should demonstrate that the machinery allocation actually improves the model fit compared to a simpler model with only the consumable depletion layer. An ablation study (consumable-only vs. full two-layer model) would strengthen the paper.

**Minor concerns:**

5. The MO/SO cycling strategy is interesting but not well-characterized. How many cycles are needed for convergence? Does the SO phase always target Venus protein (as the output suggests)? A convergence plot showing objective values across cycles would be valuable.

6. Equation 14 includes u_bar_j in the epsilon_X depletion term, meaning that the transcription control function affects resource depletion. Is this correct? If a gene is repressed (u_bar_j near 0), it should not consume NTPs. But the kinetic rate r_{X,j} already includes machinery allocation -- does double-counting occur between r_{X,j} and u_bar_j?

7. The polysome factor K_P in Eq. 6 is mentioned as a constant but its value is not given in Table 1 or the text. What value was used?

8. The units of alpha_X and alpha_L in Table 1 are (uM * nt * hr)^-1 and (uM * aa * hr)^-1. These are unusual units. Could the authors provide a physical interpretation? For example, what NTP concentration does alpha_X correspond to if the total NTP pool in PURExpress is known?

**Recommendation:** Minor revision. The modeling framework is sound, but the identifiability and filtering concerns should be addressed transparently.

---

## Reviewer 3 (Biosensor design / synthetic biology applications)

**Overall assessment:** The paper presents a modeling framework for a cell-free gluconate biosensor. While the mathematical framework is rigorous, the practical utility of the work for the biosensor community is limited. The paper reads more as a modeling exercise than a biosensor study.

**Major concerns:**

1. **Limited practical relevance of gluconate sensing.** The introduction claims gluconate detection is "relevant to food safety and industrial fermentation monitoring" but provides no citations or context for this claim. What are the actual applications? What concentration ranges matter in practice? How does the biosensor's dynamic range (0.5-1.8 uM Venus over 0.1-10 mM gluconate) compare to existing detection methods? Without this context, the choice of gluconate as the target analyte appears arbitrary.

2. **No comparison to other cell-free biosensor models.** The authors compare their resource model to their own prior work (Adhikari et al. 2020) but do not discuss how their approach relates to other resource competition models in the cell-free literature. Several groups have proposed resource-aware models for cell-free systems (e.g., Gyorgy and Murray 2016 on resource competition in genetic circuits; Weisse et al. 2015 on host-circuit interactions; Stogbauer et al. 2012 on cell-free TX/TL resource limitations). The authors should position their two-layer framework relative to these existing approaches.

3. **The biosensor performance metrics are not adequately characterized.** For a biosensor paper, the authors should report: (a) limit of detection (LOD), (b) dynamic range, (c) response time, (d) selectivity against common interferents, and (e) comparison to existing gluconate detection methods. The dose-response curve provides some of this information implicitly, but the authors should extract and report these standard metrics.

4. **Single circuit, single analyte.** The paper demonstrates the resource model on one circuit with one analyte. The claim that the framework provides "quantitative design rules for allocating transcriptional and translational capacity in synthetic gene circuits" (Conclusions) is not supported without testing on additional circuits. Can the authors predict what would happen with a different gene dosage ratio? A different promoter? A different repressor?

**Minor concerns:**

5. The model predicts 81% transcription resource depletion by 12 hours. This is an interesting finding, but the practical implication for biosensor operation is not discussed. If the biosensor is limited by NTP depletion, would supplementing NTPs extend the operational window? This is a straightforward experiment that would validate the model prediction.

6. The title emphasizes "mechanistic resource competition" but the biosensor application is equally important. Consider revising the title to better reflect both contributions.

7. The paper would benefit from a schematic showing the two-layer resource model architecture (machinery allocation + consumable depletion) as a separate figure. Currently, the reader must reconstruct this from the equations.

8. The bibliography has only 27 references. For a journal paper combining biosensing, modeling, and cell-free systems, this is thin. Recent reviews and primary research in cell-free biosensors (2022-2025) should be cited.

**Recommendation:** Major revision. The modeling is solid but the biosensor framing needs significant strengthening, and the generalizability of the approach should be demonstrated or at least discussed more thoroughly.

---

## Summary of Key Issues Across Reviewers

| Issue | R1 | R2 | R3 | Severity |
|-------|----|----|-----|----------|
| GntR protein not validated | X | | | Major |
| Ensemble filtering uses validation data | X | X | | Major |
| Parameter identifiability | | X | | Major |
| Limited biosensor context/metrics | | | X | Major |
| No comparison to other resource models | | | X | Major |
| Sensitivity heatmap dominated by alpha | X | X | | Minor-Major |
| Missing intermediate dose time courses | X | | | Minor |
| Ablation study (two-layer vs. one-layer) | | X | | Minor-Major |
| Thin bibliography | X | | X | Minor |
| Machinery allocation may not contribute | | X | | Minor-Major |
| Generalizability not demonstrated | | | X | Major |
