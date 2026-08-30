# mitoxdi — Two-Arm Decomposition + MDI (mdi_v1.0) + Composition-Predicted / Residual MDI

Open tool released with the study *"Myeloid Cell Composition Confounds Cell-Death
Transcriptomic Signatures in Critical Illness: Pre-Registered Deconfounding Reveals a
Myeloid-Intrinsic Mitochondrial Suppression Axis"*.

**What it computes** (all frozen under `score_version = mdi_v1.0`,
`gene_set_version = Mitoxy-80_v1.0`):

1. **Two-arm decomposition** — arm membership is read from the single source of truth
   `data/Mitoxyperilysis_Gene_Manifest_v1.0.csv` (upstream-collapse arm, 30 genes;
   execution-induction arm, 33 genes; HGNC-alias-normalized).
2. **mdi_v1.0 scoring** — per sample, UCS = within-cohort z of the upstream-arm mean
   (log2 space), EIS = within-cohort z of the execution-arm mean,
   MDI = EIS − UCS. z is computed over **all samples in the cohort (ungrouped)** to
   prevent leakage. Sensitivities: `MDI_nomt` (6 MT-* genes dropped) and running-sum
   ssGSEA (β = 0.25) arm scores.
3. **Composition-predicted MDI** — NNLS deconvolution against purified-cell reference
   spectra (Monaco GSE107011, 29 immune types incl. neutrophils, log2(TPM+1);
   ABIS RNA-seq 17 types / ABIS Micro 11 types, log2(+1); OLS-CLS method cross-check),
   synthetic profile X̂_comp(s) = Σ_c π_c(s)·μ_c scored under the identical mdi_v1.0
   definition → UCS_comp / EIS_comp / MDI_comp.
4. **Composition-residual MDI** — per-cohort OLS of observed MDI on MDI_comp:
   MDI_resid = MDI − (b0 + b1·MDI_comp). Reporting discipline: the three readouts
   (unadjusted / composition-adjusted / composition-residual) are always co-reported.
5. **ISED / confusability readout** — per-cohort R² (observed vs composition-predicted),
   disease-effect attenuation g_obs → g_resid, and group-wise neutrophil/monocyte
   proportion shifts.

**Judgment-gate language:** see [`GATE.md`](GATE.md) — the generalization of the
prospectively registered M10 Test A gate (osf.io/3y2rb) that the tool implements.

## Installation (zero-download for the project run; zero new dependencies)

```bash
pip install -r requirements.txt   # numpy, pandas, scipy — already satisfied by any
                                  # standard scientific environment
```

No downloads are needed for the project-local run: the Monaco reference cache
(`_intermediate/M10A_monaco_ct_means.csv`), the ABIS signature matrices
(`00_RAW_DATA/ABIS/`), the 80-gene manifest, and the two demo input matrices are all
resolved from the local project tree (paths are configurable via the
`MITOXDI_PROJECT_ROOT` environment variable). For the public release, the reference
spectra and demo inputs are distributed via Zenodo (DOI filled at release; the GPL
probe→symbol maps are exported from Bioconductor annotation packages, also archived
with the DOI).

## Quick start

```python
import sys
sys.path.insert(0, "M15_tool")
from mitoxdi.mdi_core import load_manifest, arm_scores
from mitoxdi.composition import load_monaco_ref, score_and_compose

_, arms, _ = load_manifest()
monaco = load_monaco_ref()
out, prop_monaco, _ = score_and_compose(gene_matrix, arms, monaco)
# out: UCS/EIS/MDI, MDI_nomt, ssGSEA sensitivities, UCS_comp/EIS_comp/MDI_comp,
#      neutrophil/monocyte proportions, MDI_resid, b0, b1
```

One-click adoption demos (third-party datasets, not analyzed in the study):

```bash
python demos/demo_GSE157103.py   # COVID-19 whole-blood RNA-seq (Overmyer 2020, TPM)
python demos/demo_GSE66099.py    # pediatric SIRS/sepsis/septic shock (GPSSSI, GPL570)
```

Both demos were verified against official GEO metadata before enlistment
(red line #1: verify-before-citing; record:
`03_LOGS/M15_demo_GEO_verification_20260827.md`).

## Internal regression test (11 cohorts)

The tool core is regression-tested sample-by-sample against the study's frozen
per-sample score tables (produced under the prospective registrations):

| Frozen table | Cohorts | Checked quantities |
|---|---|---|
| `_intermediate/M2_per_sample_scores.csv` | 6 (GSE185263, GSE32707, GSE212865, GSE188309, GSE310929, GSE148871) | UCS/EIS/MDI, nomt, ssGSEA |
| `_intermediate/M10A_synthetic_scores.csv` + `M10A_residual_per_sample.csv` | 6 (same) | MDI_comp (Monaco\|NNLS), MDI_resid, b0, b1 |
| `_intermediate/M11M12_step1_per_sample.csv` | 5 (GSE215865, GSE54514, GSE148871, GSE106878, GSE212865) | scores + composition + residual + proportions |

Run: `python tests/test_regression_11cohorts.py` (all-PASS ⇒ exit code 0).
Result (2026-08-27): **254 family-level checks across the 11 cohorts, 254/254 PASS
(cohort-level 17/17); largest family max |Δ| = 3.77×10⁻¹³ (machine precision)**;
the packaged scoring core is bit-identical to the project's shared `mdi_lib` on the
same matrix (max |Δ| = 0).
Summary tables: `_intermediate/M15_regression_summary.csv` (per family max |Δ|) and
`M15_regression_cohort_summary.csv`.

## Repository layout

```
M15_tool/
├── mitoxdi/            # the tool package
│   ├── mdi_core.py     #   two-arm decomposition + mdi_v1.0 scoring (frozen)
│   ├── composition.py  #   reference spectra, NNLS/OLS-CLS, synthetic + residual MDI
│   ├── report.py       #   contrast stats, R², attenuation, ISED report writer
│   └── loaders.py      #   per-cohort loaders (regression + demos)
├── demos/              # one-click adoption demos (2 third-party datasets)
├── tests/              # 11-cohort internal regression test
├── data/               # bundled 80-gene manifest (SHA-verified copy)
├── GATE.md             # judgment-gate language
├── LICENSE             # MIT
└── requirements.txt
```

## Input conventions

- **Gene-level matrix** (RNA-seq): samples × genes in log2 space, symbols uppercased;
  arm genes are resolved through the manifest alias table (legacy names such as
  DFNA5→GSDME handled), duplicates averaged per sample.
- **Series matrix** (microarray): the official GEO `!series_matrix`; probes are
  collapsed to genes with the frozen max-mean rule (per gene, ≤3 probes with the
  highest overall mean, per-sample mean), symbol mapping from the GEO GPL annotation
  or a Bioconductor probe→symbol export.
- **Groups** for contrast reporting are passed as a per-sample series; scoring itself
  never uses groups (leakage prevention).

## Reproducibility

- Version stamps: `mitoxdi v1.0.0`, `score_version = mdi_v1.0`,
  `gene_set_version = Mitoxy-80_v1.0`.
- Frozen reference caches are identical to those used in the study (SHA256 noted in
  the audit governance directory); the bundled manifest is SHA256-identical to the
  canonical governance copy (`346537be…e9965`).
- All randomness is fixed (`seed = 0` in bootstrap CIs); the regression test is
  deterministic.

## Citation

See `CITATION.cff`. If you use the tool, cite the study and the tool DOI (Zenodo,
assigned at release).
