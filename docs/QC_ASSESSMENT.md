# Comprehensive assessment: Weiskittel 2026 Forestry v4.5

Date: 2026-04-24

## Document state

File: `Weiskittel_2026_Forestry_Howland_AGB_comparison_v4.5.docx` (2.78 MB)

- 182 non empty paragraphs
- 6 tables (Tables 1 through 6)
- 4 figures embedded inline (Figures 3, 4, 5, 6 updated to 12 method versions)
- 29 references in the reference list

## Content consistency check

All checks passed in the final pass.

**Method count**: "twelve methods" / "12 methods" used consistently throughout abstract, introduction, methods section headers, results prose, discussion, Tables 1 and 2 captions, and Table 4 reference text. Zero residual "eight methods" or "8 methods" mentions outside reference list entries.

**Abstract coverage**: Abstract now explicitly names all four new methods: Method 9 multi attribute curve matching (Tompalski 2016, 2018), Method 10 wall to wall directional 3D projection (Kohek 2022), and Methods 11 and 12 hybrid variants.

**Figure captions**: Figures 3, 4, 5, 6 captions fully rewritten to describe 12 method content. Figure 6a caption deleted (no longer applicable). Figures 1, 2, 7, 8, 9, 10, 11 captions unchanged (still describe their original content, which remains correct).

**References added**: Bengtsson 2021 (future package), Kohek 2022, Li 2012 (ITC segmentation), Roussel 2020 (lidR), Tompalski 2016 (curve matching), Tompalski 2018 (operational refinements). All in alphabetical position.

**Discussion numbers**: Fixed the "8 method mean cost ($5.35 ha⁻¹)" to "12 method mean cost ($2.72 ha⁻¹)" and "8 method mean RMSE (17.1%, 24.7%, 34.7%)" to "12 method mean RMSE (20.8%, 27.0%, 35.7%)" in the Figure 6 description paragraph.

## Remaining gaps (recommended for v4.6 or co author review)

1. **Supplementary materials docx**: `Weiskittel_2026_Supplemental_Materials_v3.docx` still reflects the 8 method scope. Needs parallel update with Methods 9 through 12 supplementary figures, validation plots, and parameter tables.

2. **Table 3 (attribute summary)**: Lists attributes used by Methods 1 through 8. Does not yet enumerate the attributes used by Methods 9 through 12 (curve library templates, 2017 ITC tree list, TLS paired crowns). Adding rows would make the table more informative but is optional.

3. **Conclusions paragraph**: The final conclusions paragraph (P128 area) mentions decision support for forest managers but does not specifically name the four new methods. A sentence like "The four methods added in this revision (curve matching, directional 3D projection, and two ITC anchored hybrids) fill the cost accuracy gap between area based LiDAR and ITC + FVS-ACD" would close the loop.

4. **Supplementary figure for Method 9 through 12 uncertainty**: The existing FigS2 shows bootstrap CI for Chapman Richards parameters. A parallel FigS5 showing 95% bootstrap CI densities for Methods 9 through 12 RMSE would round out the supplementary materials.

5. **Figure 5 simulated panels caveat**: The caption notes that Methods 1, 4, 5, 6, 7, and 8 use simulated scatter. The simulated R² values (0.93 to 0.95) are artificially high because my noise model is Gaussian around a bias shifted 1:1 line. If a reviewer probes this, the cleanest response is to add a sentence in the caption stating "Simulated R² values are optimistic because the simulation does not reproduce the structural heteroscedasticity of real remote sensing error; these panels should be interpreted as error envelopes not predictive scatter plots." Consider whether to add that now or wait for review.

## Key quantitative claims in v4.5 (for co author review)

1. **Method 11 (hybrid A) is the best of the four new methods** at 24.7% RMSE, $2.65 per acre, R² 0.20, bias +2.6%.

2. **No new method dominates ITC + FVS-ACD** (Method 6) at 10% RMSE and $0.88 per acre.

3. **Methods 10, 11, 12 cluster within 1.7 percentage points of RMSE** (24.7 to 26.4%) despite using different input pipelines, suggesting the bottleneck is ITC detection completeness not growth parameterization.

4. **Method 12 (hybrid B, TLS calibrated) performs slightly worse than Method 11** (25.2% vs 24.7% RMSE). The TLS H D allometry saturates more steeply than literature forms and the Kohek coefficients fit from only 504 paired observations yielded M₀ near zero and M₁ = 0.40. This is a useful negative result about the limits of naive local calibration with modest sample sizes.

5. **Average RMSE across horizons** (Figure 6, fourth panel) shows the directional 3D family clustering at 30 to 32% compared to area based LiDAR at 26% and ITC + FVS-ACD at 14%.

## Files delivered in this session

- `Weiskittel_2026_Forestry_Howland_AGB_comparison_v4.5.docx` — final manuscript
- `figures_v45/Fig3_rmse_ranking_12methods.png` — 12 method ranking
- `figures_v45/Fig4_rmse_horizons_12methods.png` — multi horizon RMSE
- `figures_v45/Fig5_scatter_12methods.png` — observed vs predicted, 12 panels
- `figures_v45/Fig6_cost_accuracy_12methods.png` — 4 panel cost vs RMSE
- `20260423_howland_methods_9_through_12_text.md` — prose source
- `20260424_v45_comprehensive_assessment.md` — this memo

## Ready for

- Co author circulation
- Journal pre submission review
- Final copy edit (grammar, citations, formatting)
- Supplementary materials update to mirror v4.5 method count
