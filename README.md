# Composition Deconfounding of Cell-Death Transcriptomic Signatures in Critical Illness

Analysis code, open tool, and governance artifacts for:

> **Cell composition and cell-death signatures in sepsis and acute respiratory distress syndrome** (2026)

Yurong Huang, Huamei Li, Huangwei Lin, Xiang Wang, Yuanshui Liu

Most human evidence that regulated cell death drives sepsis and acute respiratory distress
syndrome comes from death-pathway gene sets scored in bulk blood or lavage transcriptomes.
Critical illness markedly reshapes leukocyte proportions, and death-execution genes are
constitutively myeloid-enriched, so a bulk cell-death signature can rise without any change
inside a single cell. This project supplies the control that this literature lacks.

The study decomposes an 80-gene regulated-cell-death framework into a 30-gene
mitochondrial / ER–mitochondria contact (MAM) infrastructure arm (UCS) and a 33-gene
death-execution arm (EIS), and scores their dissociation (MDI = EIS − UCS, fixed definition
`mdi_v1.0`). A composition-only null built from deconvolved cell proportions alone explains
53.7–81.7% of the variance in that composite index across critical-illness cohorts, and
across 249 canonical gene sets disease effects are predicted by myeloid gene content rather
than by execution-machinery content. What survives composition control is a reproducible
suppression of mitochondrial and ER–mitochondria contact transcripts inside monocytes and
macrophages (pooled g = −1.00, I² = 0%), with a null or opposite pattern in neutrophils. The
same repository ships `mitoxdi`, which re-scores any published signature against its own
composition null.

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

## Related records

| Record | Identifier |
|---|---|
| Analysis code and the standalone `mitoxdi` software (Zenodo) | [10.5281/zenodo.22655285](https://doi.org/10.5281/zenodo.22655285) |
| Large supplementary matrices (Zenodo) | [10.5281/zenodo.23018056](https://doi.org/10.5281/zenodo.23018056) |
| Frozen analysis plan (OSF) | [10.17605/OSF.IO/C7RYD](https://doi.org/10.17605/OSF.IO/C7RYD) |

## License

- Code (including mitoxdi): **MIT** — see `LICENSE` and `mitoxdi_tool/LICENSE`.
- Governance documents and figure data: **CC BY 4.0**.

## Citation

See `CITATION.cff` (and `mitoxdi_tool/CITATION.cff` for the tool). A Zenodo DOI is minted
for each release; please cite the release used.

If you cite the study:

```bibtex
@article{huang2026composition,
  title  = {Cell composition and cell-death signatures in sepsis and acute respiratory distress syndrome},
  author = {Huang, Yurong and Li, Huamei and Lin, Huangwei and Wang, Xiang and Liu, Yuanshui},
  year   = {2026},
  note   = {Analysis code and software: https://doi.org/10.5281/zenodo.22655285}
}
```
