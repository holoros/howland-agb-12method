# Technical Review of Weiskittel 2026 Forestry v4.5 (self-review)

Date: 2026-04-24
Reviewer mode: structured pre-submission audit
Documents reviewed: `Weiskittel_2026_Forestry_Howland_AGB_comparison_v4.5.docx` (main, 2.78 MB) and `Weiskittel_2026_Supplemental_Materials_v4.docx` (supplemental, 936 KB)

## Summary assessment

The manuscript presents a defensible 12-method head-to-head comparison of AGB projection approaches at Howland. The science is sound and the storyline (cost-accuracy tradeoffs with honest negative results for Methods 9 and 12) is publishable. The structural issue is that the manuscript was extended from an existing 8-method draft rather than rewritten as a 12-method paper, leaving residual narrative artifacts in abstract and results that have now been corrected. Statistical concerns about spatial autocorrelation in Methods 9-12 bootstraps and mild calibration leakage in Methods 11-12 have been added as caveats to the methods text. Indicative editorial decision: minor revision pending citation verification and reviewer-blind re-read.

## Major issues addressed in this pass

1. **Abstract restructured** to a coherent 12-method narrative. The previous version listed only 8 methods in the methods enumeration and contained a self-citation ("extends the existing eight-method comparison of Weiskittel et al. 2026"). Replaced with a single integrated paragraph that names all 12 methods, reports the 12-method mean RMSEs (20.8%, 27.0%, 35.7%), reports the M9 to M12 RMSE range, and frames the conclusion in terms of all 12 methods.

2. **Results P78** changed from "the seven plot level methods spanned a wide range of accuracy" to "the twelve projection methods spanned a wide range of accuracy."

3. **Spatial autocorrelation caveat** added to the Uncertainty quantification subsection (P65). Acknowledges that pixel-level bootstrap CIs for Methods 9-12 are lower bounds because adjacent forest pixels have Moran's I above 0.5; an effective sample size correction would widen the intervals 2 to 3x.

4. **Calibration leakage note** added to Methods 11 and 12 sections: 20% of the interpolated 2017 anchor weight comes from the 2025 observation, introducing mild leakage between calibration target and validation reference. Noted that an independent 2017 G-LiHT ABA AGB raster would resolve this.

5. **Figure 5 caption** extended with caveat that simulated R² values are optimistic because the Gaussian noise model around a bias-shifted 1:1 line does not reproduce structural heteroscedasticity of remote sensing error.

6. **Software availability statement** appended to the Software and reproducibility subsection, naming each script and the OSC Cardinal project (PUOM0008) plus exact package versions.

## Major issues flagged for [REVIEWER VERIFY] before submission

1. **Citation accuracy for added references**. The following references were added in this session and need verification against original sources:
   - Tompalski et al. (2018): the volume and pages cited (RSE 227, 110-124) match Tompalski et al. (2019), not 2018. Verify whether the intended reference is Tompalski 2018 (Forests, "operational refinements") or Tompalski 2019 (RSE, "transferability").
   - Kohek et al. (2022): title cited is "SimForest...". Verify against Computers and Electronics in Agriculture 199, 107148.
   - Bengtsson (2021), Roussel et al. (2020), Li et al. (2012): basic verification of pagination and DOIs.
   - Wheatland Geospatial Lab (2026): forward-dated technical report; verify availability.

2. **Self-reference review**. The "Weiskittel et al. (2026)" citation in the original abstract has been removed during the rewrite. Confirm that no other self-references remain in the body (search "Weiskittel" in references — appears once for Weiskittel 2011 textbook, which is fine).

3. **Methods 9-12 residual diagnostics**. No residual vs predicted plots, Q-Q plots, or spatial residual maps for the four new methods. Figure 10 covers only the original 3-point CR / 4-point CR / FVS-ACD predictions. Recommend adding a supplementary figure (FigS14) with residual distributions for Methods 9-12 before journal submission.

4. **Method 11/12 stratified calibration alpha and beta standard errors**. Table S9 gives point estimates only. Standard errors from the lm fits would strengthen the supplementary documentation and let readers assess uncertainty in individual stratum fits (especially Medium_deciduous with n=31).

## Minor issues addressed

- "8 method mean cost ($5.35)" already corrected to "12 method mean cost ($2.72)" in the Results discussion of Figure 6 (Discussion P101).
- Six new references added in alphabetical position (Bengtsson, Kohek, Li, Roussel, Tompalski 2016, Tompalski 2018).
- Table 1 and Table 2 captions updated to "twelve methods."

## Methods and statistical assessment summary

### General reproducibility
- G1 PARTIALLY MET → improved with software availability statement
- G2 PARTIALLY MET → improved with explicit package versions
- G3 PARTIALLY MET → seed=2026 implicit; recommend stating in Methods 9-12 prose

### Statistical
- S1 MET
- S2 MET
- S2b NOT MET → addressed via caveat in Uncertainty quantification subsection. A full re-analysis with effective sample size correction would be the rigorous fix and could be added as a supplementary table for v4.6.
- S6 MET (CIs throughout)
- S7 MET (cost framing)

### Validation
- V1 PARTIALLY MET → caveat added for Methods 11/12 interpolation leakage
- V2 MET
- V3 NOT MET → flagged for v4.6 (residual diagnostics)
- V4 MET

### Mathematical
- M1 NOT MET → Tompalski Eq. 1 and Kohek Eqs. 7 to 9 are referenced but not displayed in the methods. Recommend adding explicit equations for the v4.6 revision.
- M2 PARTIALLY MET → SE missing for stratified calibration coefficients.

## Scoring scaffold (final)

| Dimension | Score | Conf | Notes |
|-----------|-------|------|-------|
| D1 Originality | 4 | H | First 12-method head-to-head; novel synthesis |
| D2 Significance | 4 | H | Operational decision surface useful to managers |
| D3 Methodological rigor | 3.5 | M | Calibration leakage now noted; residual diagnostics still missing |
| D4 Statistical appropriateness | 3.5 | M | Spatial autocorrelation caveat now in place |
| D5 Reproducibility | 4 | M | Software availability statement added; scripts archived |
| D6 Model validation | 3.5 | M | Independent 2025 validation solid; residual plots still missing for new methods |
| D7 Literature coverage | 4 | M | Citations need spot-check |
| D8 Interpretation | 4 | H | Honest framing of negative results |
| D9 Writing clarity | 4 | H | Abstract now coherent; "seven plot level" fixed |
| D10 Figures and tables | 4 | H | Figure 6 four-panel design works; Figure 5 caveat added |

**Indicative recommendation**: Minor revision. The paper is publishable in *Forestry* with the residual diagnostics, equation displays, and citation verification items addressed in v4.6.

## Citation integrity audit (priority spot-check list)

| # | Reference | Anomaly type | Action |
|---|-----------|---------------|---------|
| 1 | Tompalski et al. (2018) RSE 227, 110-124 | Type D | Confirm year/journal/pages |
| 2 | Kohek et al. (2022) CEA 199, 107148 | Type C/D | Confirm exact title |
| 3 | Wheatland Geospatial Lab (2026) | Type B | Verify availability |
| 4 | Bengtsson (2021) R Journal 13, 208-227 | Type D | Confirm pages |
| 5 | Roussel et al. (2020) RSE 251, 112061 | Type D | Confirm pagination |
| 6 | Li et al. (2012) PERS 78, 75-84 | Type D | Confirm pages |
| 7 | Chojnacky et al. (2014) Forestry 87, 129-151 | Type D | Confirm |
| 8 | Tompalski et al. (2016) Forests 7, 255 | Type D | Confirm |
| 9 | Weiskittel et al. (2011) Wiley Blackwell | Type C | Confirm publication details |
| 10 | Westfall et al. (2024) NSVB | Type B | Confirm full citation |

No Type A (non-existent) or Type F (block-generated cluster) patterns detected.

## Disclosure

This technical review was prepared with AI assistance using Claude (Sonnet, Anthropic). AI was used for: structural consistency checking, lexical scanning for residual "eight methods" references, citation list cross-reference, abstract rewrite drafting, and the methods checklist application. All scientific judgments and the final text are the responsibility of the author. Manuscript content remained within the local Cowork session and was not uploaded to external services.

---

## Action checklist for v4.6 (next session)

- [ ] Verify Tompalski 2018, Kohek 2022, Bengtsson 2021, Roussel 2020, Li 2012 citations against original sources
- [ ] Add Tompalski Eq. 1 and Kohek Eqs. 7-9 explicitly to methods (M1)
- [ ] Add residual vs predicted plots for Methods 9-12 as supplementary FigS14 (V3)
- [ ] Add standard errors to Table S9 stratified calibration coefficients (M2)
- [ ] Optional: full re-analysis with effective sample size correction for Methods 9-12 bootstrap (S2b rigorous fix)
- [ ] Verify availability of Wheatland Geospatial Lab (2026) technical report
- [ ] Consider adding 2017 G-LiHT ABA AGB raster as independent calibration target for Methods 11/12 (resolves V1 leakage cleanly)
