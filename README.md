# Composition Deconfounding of Cell-Death Transcriptomic Signatures in Critical Illness

Analysis code, open tool, and governance artifacts for:

> **Deconfounding cell-death signatures in sepsis and ARDS reveals monocyte-intrinsic suppression of mitochondrial infrastructure** (manuscript under review)

The study decomposes an 80-gene regulated-cell-death framework into an upstream-collapse arm
(UCS, 30 genes) and an execution-induction arm (EIS, 33 genes), scores their dissociation
(MDI = EIS − UCS, fixed definition `mdi_v1.0`), and shows that (i) blood compositional drift
(myeloid fraction) is the dominant confounder of published RCD transcriptomic readouts, while
(ii) an upstream myeloid-intrinsic mitochondrial suppression axis survives composition
correction across 4 independent single-cell cohorts.

## Contents

| Directory | Contents |
|---|---|
| `analysis_scripts/` | 105 analysis/QC scripts (P0 frozen-plan core re-analysis, M1 spatial, M2 clinical MDI, M3 causal four-layer pQTL-MR/SMR/TWAS, M4 biochemical, composition-hub/temporal/intervention modules, package lint) |
| `mitoxdi_tool/` | **mitoxdi v1.0.0** — standalone MIT-licensed tool: dual-arm decomposition + `mdi_v1.0` scoring + composition-predicted / composition-residual MDI (own README, LICENSE, demos, 11-cohort regression tests, gate document `GATE.md`) |
| `governance/` | RESULTS_MANIFEST v2.0 (SHA256-audited inventory), gene manifest v1.0 (single source of the 80-gene/arm definitions), sample manifest, frozen analysis plan, judgment-gate reports, audit reports, novelty-search logs, software snapshot |
| `01_FIGURE_DATA_CSV/` | Per-panel CSV data behind the main and supplementary figures (incl. the wet-lab Figure 8 panels, `Main/Figure_8w*`; see that folder's README for the ms/pkg numbering bridge) |

## Data availability

All data are public: GEO series (e.g., GSE145926, GSE158055, GSE185263, GSE212865, GSE32707,
GSE165659, GSE67530, GSE271370, GSE253474, GSE310929, GSE188309, GSE148871, GSE216009,
GSE180578, GSE215865, GSE54514, GSE106878, GSE235046, GSE221321, GSE157103, GSE66099),
MetaboLights MTBLS6844, and genetics reference resources (IEU OpenGWAS, FinnGen R10,
GTEx v8, eQTLGen, deCODE pQTL, UKB-PPP). The authoritative per-file source mapping is the
`location` column of `governance/RESULTS_MANIFEST_v2.0.csv`. Supplementary tables accompany
the article. No controlled-access data were used.

## Reproducing

Scripts were executed on Windows in the frozen environment recorded in
`governance/P0_Software_Snapshot_20260821.md` (Python 3.14 / R 4.6.0; pip freeze ×2 and
`sessionInfo` snapshots included). Three notes:

1. **Project root**: scripts reference the project root by absolute path (Windows). Set the
   environment variable `MITOXDI_PROJECT_ROOT` to your unpacked project directory for the
   mitoxdi package; for the full pipeline, re-download public data into `00_RAW_DATA/`
   per the manifest and run scripts by module order (P0 → M1–M4).
2. **Intermediate objects (~12 GB)** are not redistributed; they are rebuilt by the scripts.
   Random seeds are fixed throughout.
3. **IEU OpenGWAS credentials**: two M3 scripts read a local `.opengwas_token` file
   (path reference only; no credentials are stored in this repository). Provide your own
   token from https://ieu-open-gwas-api.readthedocs.io if re-running M3.

The pre-registered analysis plans are archived openly (OSF):
10.17605/OSF.IO/C7RYD (frozen plan), 10.17605/OSF.IO/ETVMJ (M10/M13/M14),
10.17605/OSF.IO/98CM3 (M11/M12).

## License

- Code (including mitoxdi): **MIT** — see `LICENSE` and `mitoxdi_tool/LICENSE`.
- Governance documents and figure data: **CC BY 4.0**.

## Citation

See `CITATION.cff` (and `mitoxdi_tool/CITATION.cff` for the tool). A Zenodo DOI is minted
for each release; please cite the release used.
