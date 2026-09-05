# -*- coding: utf-8 -*-
"""
P0_build_release_dir.py — §三-6 GitHub/Zenodo 发布目录组装（2026-08-30）
输出 05_RELEASE_GITHUB/：
  repo/            = 上传到 GitHub 的全部内容（白名单复制 + 生成元文件）
  发布操作说明.md   = 作者操作指南（不上传）
"""
import os, re, shutil, sys
from pathlib import Path
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = Path(r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS")
REL = ROOT / "05_RELEASE_GITHUB"
REPO = REL / "repo"
if REPO.exists():
    shutil.rmtree(REPO)
for d in ("analysis_scripts", "mitoxdi_tool", "governance", "figure_data"):
    (REPO / d).mkdir(parents=True)

def copy_all(patterns, dest):
    n = 0
    for pat in patterns:
        for p in ROOT.glob(pat):
            if p.is_file():
                shutil.copy2(p, dest / p.name); n += 1
    return n

def copy_tree(src, dest):
    # 排除运行产物：__pycache__/*.pyc 不进发布仓（网页拖拽上传不受 .gitignore 保护）
    ignore = shutil.ignore_patterns("__pycache__", "*.pyc")
    shutil.copytree(src, dest, dirs_exist_ok=True, ignore=ignore)
    return sum(1 for p in Path(dest).rglob("*") if p.is_file() and p.suffix != ".pyc")

n_scripts = copy_all(["*.py", "*.R"], REPO / "analysis_scripts")
n_tool = copy_tree(ROOT / "M15_tool", REPO / "mitoxdi_tool")
n_gov = copy_tree(ROOT / "04_AUDIT_GOVERNANCE", REPO / "governance")
n_gov += copy_tree(ROOT / "M_judgment_reports", REPO / "governance" / "M_judgment_reports")
n_fig = copy_tree(ROOT / "01_FIGURE_DATA_CSV", REPO / "figure_data")
print(f"复制：脚本 {n_scripts} | M15_tool {n_tool} | 治理 {n_gov} | 图数据 {n_fig}")

# ---------- 元文件 ----------
README = f"""# Composition Deconfounding of Cell-Death Transcriptomic Signatures in Critical Illness

Analysis code, open tool, and governance artifacts for:

> **Myeloid Cell Composition Confounds Cell-Death Transcriptomic Signatures in Critical Illness: Composition Deconfounding Reveals a Myeloid-Intrinsic Mitochondrial Suppression Axis** (manuscript under review)

The study decomposes an 80-gene regulated-cell-death framework into an upstream-collapse arm
(UCS, 30 genes) and an execution-induction arm (EIS, 33 genes), scores their dissociation
(MDI = EIS − UCS, fixed definition `mdi_v1.0`), and shows that (i) blood compositional drift
(myeloid fraction) is the dominant confounder of published RCD transcriptomic readouts, while
(ii) an upstream myeloid-intrinsic mitochondrial suppression axis survives composition
correction across 4 independent single-cell cohorts.

## Contents

| Directory | Contents |
|---|---|
| `analysis_scripts/` | {n_scripts} analysis/QC scripts (P0 frozen-plan core re-analysis, M1 spatial, M2 clinical MDI, M3 causal four-layer pQTL-MR/SMR/TWAS, M4 biochemical, composition-hub/temporal/intervention modules, package lint) |
| `mitoxdi_tool/` | **mitoxdi v1.0.0** — standalone MIT-licensed tool: dual-arm decomposition + `mdi_v1.0` scoring + composition-predicted / composition-residual MDI (own README, LICENSE, demos, 11-cohort regression tests, gate document `GATE.md`) |
| `governance/` | RESULTS_MANIFEST v2.0 (SHA256-audited inventory), gene manifest v1.0 (single source of the 80-gene/arm definitions), sample manifest, frozen analysis plan, judgment-gate reports, audit reports, novelty-search logs, software snapshot |
| `figure_data/` | Per-panel CSV data behind all main figures |

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
"""
(REPO / "README.md").write_text(README, encoding="utf-8", newline="")

LICENSE = """MIT License

Copyright (c) 2026 Yurong Huang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
(REPO / "LICENSE").write_text(LICENSE, encoding="utf-8", newline="")

CITATION = """cff-version: 1.2.0
message: "If you use this repository, please cite it as below."
title: "Composition deconfounding of cell-death transcriptomic signatures in critical illness: analysis code and mitoxdi"
authors:
  - family: Huang
    given: Yurong
    affiliation: Hainan General Hospital; Hainan Affiliated Hospital of Hainan Medical University
version: 1.0.0
date-released: "2026-08-30"
license: MIT
keywords:
  - ARDS
  - regulated cell death
  - cellular composition
  - deconvolution
  - mitochondria
"""
(REPO / "CITATION.cff").write_text(CITATION, encoding="utf-8", newline="")

ZENODO = {
    "title": "Composition deconfounding of cell-death transcriptomic signatures in critical illness: analysis code and mitoxdi v1.0.0",
    "description": "Analysis code (P0-M4 modules), the mitoxdi v1.0.0 open tool (dual-arm decomposition + mdi_v1.0 scoring + composition-predicted/residual MDI), and SHA256-audited governance artifacts (results manifest, gene manifest, frozen analysis plan, judgment-gate reports) for a manuscript under review.",
    "upload_type": "software",
    "creators": [{"name": "Huang, Yurong",
                  "affiliation": "Hainan General Hospital; Hainan Affiliated Hospital of Hainan Medical University"}],
    "license": "MIT",
    "version": "v1.0.0",
    "keywords": ["ARDS", "sepsis", "regulated cell death", "composition deconvolution", "mitochondria", "transcriptomics"],
    "related_identifiers": [
        {"identifier": "10.17605/OSF.IO/C7RYD", "relation": "documents"},
        {"identifier": "10.17605/OSF.IO/ETVMJ", "relation": "documents"},
        {"identifier": "10.17605/OSF.IO/98CM3", "relation": "documents"},
    ],
}
import json
(REPO / ".zenodo.json").write_text(json.dumps(ZENODO, ensure_ascii=False, indent=2), encoding="utf-8", newline="")

GITIGNORE = """__pycache__/
*.pyc
.venv*/
venv/
.opengwas_token
opengwas_token
.DS_Store
*.h5ad
*.mtx
*.mtx.gz
_intermediate/
00_RAW_DATA/
R_libs/
"""
(REPO / ".gitignore").write_text(GITIGNORE, encoding="utf-8", newline="")
print("元文件已生成：README/LICENSE/CITATION.cff/.zenodo.json/.gitignore")

# ---------- 发布前复查 ----------
files = [p for p in REPO.rglob("*") if p.is_file()]
total = sum(p.stat().st_size for p in files)
big = [(p.stat().st_size, p) for p in files if p.stat().st_size > 50 * 1024 * 1024]
pats = {"JWT": r"eyJ[A-Za-z0-9_-]{10,}", "密钥赋值": r"(api_key|secret)\s*=\s*[\"\x27][^\"\x27]{8,}"}
secrets = []
for p in files:
    if p.suffix.lower() in {".py", ".r", ".md", ".json", ".txt", ".cff", ".csv"} and p.stat().st_size < 5_000_000:
        try: t = p.read_text(encoding="utf-8", errors="ignore")
        except Exception: continue
        for k, pat in pats.items():
            if re.search(pat, t): secrets.append((p.name, k))
print(f"\n发布物：{len(files)} 文件 / {total/1024/1024:.1f} MB | >50MB 文件: {len(big)} | 密钥命中: {len(secrets)}")
assert not big and not secrets, "存在大文件或疑似密钥，中止"
print("复查通过（无 >50MB、无密钥）")
