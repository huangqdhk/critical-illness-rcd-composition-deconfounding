# GATE.md — Judgment-Gate Language (判定门语言文档化)

This file documents the **frozen judgment-gate language** implemented by the tool.
It is the generalization of the prospectively registered M10 Test A gate
(`M10_M13_M14_pre_registration_20260824.md`, OSF osf.io/3y2rb, §1.2), whose results
in the study are reported in Results §13.1. The gates below must be read **before**
running the tool on a new cohort (prospective discipline: gates precede results).

## 1. Definitions (frozen)

- `UCS` = within-cohort z of the upstream-collapse-arm mean (log2 space).
- `EIS` = within-cohort z of the execution-induction-arm mean (log2 space).
- `MDI` = EIS − UCS (`score_version = mdi_v1.0`).
- z standardization uses **all samples of the cohort, ungrouped** (no leakage).
- `MDI_comp` = the same computation applied to the composition-only synthetic profile
  X̂_comp(s) = Σ_c π_c(s)·μ_c (NNLS proportions against a purified-cell reference
  spectrum; OLS-CLS is the method cross-check).
- `MDI_resid` = MDI − (b0 + b1·MDI_comp) from a per-cohort OLS.
- Arm membership has a **single source**: `Mitoxyperilysis_Gene_Manifest_v1.0.csv`
  (upstream 30 / execution 33; HGNC-alias-normalized).

## 2. Reporting rules (mandatory)

1. Always report the **three readouts jointly**: unadjusted / composition-adjusted /
   composition-residual — never the residual alone.
2. Whole-blood deconvolution must use a **granulocyte-containing** reference
   (Monaco 29-type, ABIS RNA-seq 17-type, or an equivalent declared at run time).
   Granulocyte-free references (e.g., LM22) are invalid for whole blood.
3. Report ≥2 proportion-estimation methods (NNLS primary; OLS-CLS or CIBERSORTx
   cross-check). CIBERSORTx is infeasible offline — the OLS-CLS substitution is an
   accepted, disclosed deviation (registered).
4. Over-correction disclosure: S100A8/A9-class genes carry both composition and state
   information; composition adjustment may absorb true signal (bias direction:
   conservative). State this in any interpretation.
5. Statistical unit: samples are scored individually, but any inferential pooling
   must be at the **donor/dataset level**; cells/FOVs/spots are never independent
   units.

## 3. Decision gate for a new whole-blood cohort (M10 Test A, generalized)

Inputs: per-cohort R² (observed MDI vs MDI_comp, primary reference|method) and the
disease contrast of MDI_resid (Hedges' g with CI, and its direction consistency
across strata if multiple strata exist).

| Condition | Reading |
|---|---|
| **COMPOSITION** (Track-A direction) | R² ≥ 0.5 **and** the residual disease effect is not significant |
| **INTRINSIC** (Track-B direction) | R² < 0.5 **and** the residual effect is significant and direction-consistent (≥70% of strata when strata exist) |
| **INTERMEDIATE** | any other combination — report honestly, weighting by the supported direction; **do not force a dichotomy** |

The study itself adjudicated INTERMEDIATE: per-cohort R² ≥ 0.5 in 5/6 cohorts
(pooled-sample R² = 0.419), yet the residual MDI disease effect survived across
k = 16 source-non-overlapping layers (g = +0.588, one-sided p = 0.014, 75%
direction-consistent, LOSO 16/16) → "composition is the principal but not the sole
explanation."

## 4. Interpretation discipline (red lines)

1. No dataset enters an analysis before verification against official metadata
   (GEO/ArrayExpress records: title, design, platform, file availability).
2. No judgment on other papers' conclusions — only objective quantities
   (arm scores, composition predictions, residuals, overlaps) are stated.
3. Negative or intermediate results are reported as-is; no endpoint switching.
4. The words "law"/"disruptive" are not used as self-description.

## 5. Adoption-demo adjudication

For each third-party demo cohort the tool reports: arm coverage, per-family disease
effects (MW p, Cliff's δ, Hedges' g) for observed / composition-predicted /
composition-residual readouts, the R² of the composition prediction, the attenuation
g_obs → g_resid, and the group-wise neutrophil/monocyte proportion shifts. The
Gate §3 language is then applied verbatim; the demo results are **descriptive
recomputations**, not new claims about the datasets' biology.
