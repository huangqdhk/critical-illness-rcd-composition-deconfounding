# -*- coding: utf-8 -*-
"""
P0_donor_pseudobulk.py — Phase 0 任务5：供者级 pseudobulk 重算（H4）
====================================================================
依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md（OSF https://osf.io/C7RYD/）§2/§4
     单位：cell type × donor pseudobulk，供者=有效 n；细胞不作独立单位。
设计：
  - 主对比：COVID_severe vs Healthy（预注册单细胞主对比；mild 不入确证）
  - 分数据集（GSE145926 BALF / GSE158055 PBMC）各自独立检验；合并仅作敏感性（~dataset+condition）
  - 层1 臂评分：每 donor×celltype 的臂均值（log1p）→ 组内（celltype 内跨供者）z → UCS/EIS/MDI/MDI_nomt
            Welch t + BH（族=数据集内 celltype×score）
  - 层2 基因级：pydeseq2（DESeq2 实现）design ~ condition，全基因 BH；
            另存 80 基因框架子集便于对照
  - 最小供体数：双侧各 ≥3 供者且该 donor×celltype ≥50 细胞（未注册参数，在此披露）
输入：00_RAW_DATA/_scrna_work/combined_processed.h5ad（counts 层）
     02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV/Mitoxyperilysis_Gene_Manifest_v1.0.csv
输出：_intermediate/P0_pseudobulk_arm_scores.csv / P0_pseudobulk_DE_<dataset>.csv
     04_AUDIT_GOVERNANCE/P0_Donor_Pseudobulk_Report.md
     03_LOGS/P0_donor_pseudobulk_log.txt
"""
import os, sys, warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import ttest_ind
warnings.filterwarnings("ignore")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import scanpy as sc

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
H5 = os.path.join(ROOT, "00_RAW_DATA", "_scrna_work", "combined_processed.h5ad")
T = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES", "SUPPLEMENTARY_Tables_CSV")
MANI = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "Mitoxyperilysis_Gene_Manifest_v1.0.csv")
INTER = os.path.join(ROOT, "_intermediate")
REPORT = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "P0_Donor_Pseudobulk_Report.md")
LOGF = os.path.join(ROOT, "03_LOGS", "P0_donor_pseudobulk_log.txt")

log = []
def note(m=""):
    log.append(m); print(m, flush=True)

def bh(pvals):
    p = np.asarray(pvals, float); out = np.full_like(p, np.nan)
    m = ~np.isnan(p); ps = p[m]
    if len(ps) == 0: return out
    order = np.argsort(ps); n = len(ps); adj = np.empty(n)
    for k in range(n - 1, -1, -1):
        adj[order[k]] = min(ps[order[k]] * n / (k + 1), adj[order[k + 1]] if k + 1 < n else 1.0)
    out[m] = adj; return out

# ---------------- 0. 载入 ----------------
note("== 0. 载入")
adata = sc.read_h5ad(H5)
obs = adata.obs.copy()
counts = adata.layers["counts"]
if sp.issparse(counts):
    counts = counts.tocsr()
logX = adata.X  # log1p 归一化
varnames = list(adata.var_names)

man = pd.read_csv(MANI)
man["arm"] = man["arm"].astype(str).str.strip()
arm_up = set(man.loc[man["arm"].str.contains("上|up", case=False, na=False), "gene_symbol"])
arm_ex = set(man.loc[man["arm"].str.contains("执|ex", case=False, na=False), "gene_symbol"])
alias = {}
for _, r in man.iterrows():
    if pd.notna(r.get("aliases")):
        for al in str(r["aliases"]).split(","):
            alias[al.strip().upper()] = r["gene_symbol"]
note(f"   manifest: 上游 {len(arm_up)} / 执行 {len(arm_ex)}")

def gene_idx(gene_set):
    idx, used = [], []
    for g in gene_set:
        if g in varnames:
            idx.append(varnames.index(g)); used.append(g)
        else:
            al = alias.get(g.upper())
            if al and al in varnames:
                idx.append(varnames.index(al)); used.append(al)
    return np.array(idx, dtype=int), used

up_idx, up_used = gene_idx(arm_up)
ex_idx, ex_used = gene_idx(arm_ex)
up_nomt_idx = np.array([i for i, g in zip(up_idx, up_used) if not g.upper().startswith("MT-")], dtype=int)
note(f"   组内可测基因：上游 {len(up_used)}/{len(arm_up)}，执行 {len(ex_used)}/{len(arm_ex)}")

def zc(s):
    s = np.asarray(s, float); sd = s.std(ddof=1)
    return (s - s.mean()) / sd if sd > 0 else s * 0

# ---------------- 1. pseudobulk 矩阵 + 供者级臂评分 ----------------
MIN_CELLS = 50
MIN_DONORS = 3
CASE, CTRL = "COVID_severe", "Healthy"
datasets = ["GSE145926", "GSE158055"]
celltypes = sorted(obs["cell_type"].unique())

arm_rows = []
pb_store = {}   # (dataset, celltype) -> counts DataFrame + meta
coverage = []

note("\n== 1. 供者×细胞类型 pseudobulk 与臂评分")
for ds in datasets:
    o = obs[obs["dataset"] == ds]
    conds = o[["sampleID", "condition"]].drop_duplicates().set_index("sampleID")["condition"]
    for ct in celltypes:
        sub = o[o["cell_type"] == ct]
        don = sub.groupby("sampleID").size()
        don = don[don >= MIN_CELLS]
        donors = [d for d in don.index if conds.get(d) in (CASE, CTRL)]
        n_case = sum(conds.get(d) == CASE for d in donors)
        n_ctrl = sum(conds.get(d) == CTRL for d in donors)
        coverage.append(dict(dataset=ds, cell_type=ct, n_donors=len(donors),
                             n_case=n_case, n_ctrl=n_ctrl,
                             eligible=n_case >= MIN_DONORS and n_ctrl >= MIN_DONORS))
        if n_case < MIN_DONORS or n_ctrl < MIN_DONORS:
            continue
        # 臂评分（log1p 均值 per donor）——用 obs_names 定位到全对象行号
        pos_all = adata.obs_names.get_indexer(sub.index)
        rows = []
        mat = []
        for d in donors:
            mask = (sub["sampleID"] == d).values
            pos = pos_all[mask]
            cellsub = logX[pos, :]
            if sp.issparse(cellsub):
                mu_up = np.asarray(cellsub[:, up_idx].mean(axis=0)).ravel()
                mu_ex = np.asarray(cellsub[:, ex_idx].mean(axis=0)).ravel()
                mu_upn = np.asarray(cellsub[:, up_nomt_idx].mean(axis=0)).ravel()
            else:
                mu_up = cellsub[:, up_idx].mean(axis=0); mu_ex = cellsub[:, ex_idx].mean(axis=0)
                mu_upn = cellsub[:, up_nomt_idx].mean(axis=0)
            rows.append(dict(dataset=ds, cell_type=ct, sampleID=d, condition=conds[d],
                             UCS_raw=mu_up.mean(), EIS_raw=mu_ex.mean(), UCS_nomt_raw=mu_upn.mean(),
                             n_cells=int(mask.sum())))
            # counts pseudobulk
            ci = counts[pos, :].sum(axis=0)
            mat.append(np.asarray(ci).ravel())
        dfr = pd.DataFrame(rows)
        dfr["UCS"] = zc(dfr["UCS_raw"]); dfr["EIS"] = zc(dfr["EIS_raw"])
        dfr["MDI"] = dfr["EIS"] - dfr["UCS"]
        dfr["UCS_nomt"] = zc(dfr["UCS_nomt_raw"]); dfr["MDI_nomt"] = dfr["EIS"] - dfr["UCS_nomt"]
        arm_rows.append(dfr)
        pb_df = pd.DataFrame(np.vstack(mat), index=[f"{ds}|{ct}|{d}" for d in donors], columns=varnames).astype(int)
        pb_store[(ds, ct)] = (pb_df, dfr.set_index("sampleID")["condition"])

cov = pd.DataFrame(coverage)
note(cov.to_string())
arm = pd.concat(arm_rows, ignore_index=True) if arm_rows else pd.DataFrame()

# 臂评分组间检验（族 = 数据集内 celltype×score）
note("\n== 2. 供者级臂评分检验（Welch t，BH 族=数据集内 celltype×score）")
tests = []
for ds in datasets:
    dfr_all = arm[arm["dataset"] == ds] if len(arm) else pd.DataFrame()
    fam = []
    for ct in dfr_all["cell_type"].unique():
        dfr = dfr_all[dfr_all["cell_type"] == ct]
        for score in ("UCS", "EIS", "MDI", "MDI_nomt"):
            a = dfr.loc[dfr["condition"] == CASE, score].values
            b = dfr.loc[dfr["condition"] == CTRL, score].values
            if len(a) >= MIN_DONORS and len(b) >= MIN_DONORS:
                t, p = ttest_ind(a, b, equal_var=False)
                fam.append(dict(dataset=ds, cell_type=ct, score=score,
                                n_case=len(a), n_ctrl=len(b),
                                case_mean=a.mean(), ctrl_mean=b.mean(),
                                diff=a.mean() - b.mean(), t=t, p=p))
    fam = pd.DataFrame(fam)
    fam["BH_q"] = bh(fam["p"].values) if len(fam) else np.nan
    tests.append(fam)
    for _, r in fam.iterrows():
        note(f"   {ds} {r['cell_type']} {r['score']}: diff={r['diff']:+.3f} p={r['p']:.3g} q={r['BH_q']:.3g}")
test_df = pd.concat(tests, ignore_index=True) if tests else pd.DataFrame()

# ---------------- 3. pydeseq2 基因级 ----------------
note("\n== 3. pydeseq2 基因级 DE（design ~ condition，每数据集×细胞类型）")
de_ok = True
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats   # 0.5.x：DeseqStats 位于 ds 模块
except Exception as e:
    de_ok = False
    note(f"   pydeseq2 不可用（{e}）→ 基因级层未运行，见报告披露")

de_all = []
if de_ok:
    for (ds, ct), (pb, cond) in pb_store.items():
        keep = pb.sum(axis=0) > 0
        pb2 = pb.loc[:, keep.values]
        try:
            meta = pd.DataFrame({"condition": cond.values}, index=pb2.index)
            dds = DeseqDataSet(counts=pb2, metadata=meta, design="~condition", quiet=True, n_cpus=1)
            dds.deseq2()
            st = DeseqStats(dds, contrast=["condition", CASE, CTRL], quiet=True, n_cpus=1)
            st.summary()
            res = st.results_df.reset_index().rename(columns={"index": "gene"})
            res.insert(0, "cell_type", ct); res.insert(0, "dataset", ds)
            de_all.append(res)
            sig = (res["padj"] < 0.05).sum()
            note(f"   {ds} {ct}: n={pb2.shape[0]} 供者, 检验基因 {pb2.shape[1]}, padj<0.05 {sig}")
        except Exception as e:
            note(f"   {ds} {ct}: 失败（{type(e).__name__}: {e}）→ 按注册回退记录，未替代")
    if de_all:
        de_df = pd.concat(de_all, ignore_index=True)
        de_df.to_csv(os.path.join(INTER, "P0_pseudobulk_DE_all.csv"), index=False)
        fw = de_df[de_df["gene"].isin(set(up_used) | set(ex_used))]
        fw.to_csv(os.path.join(INTER, "P0_pseudobulk_DE_80genes.csv"), index=False)
        note(f"   基因级结果：全表 {len(de_df)} 行；80 基因子集 {len(fw)} 行")

# ---------------- 输出 ----------------
arm.to_csv(os.path.join(INTER, "P0_pseudobulk_arm_scores.csv"), index=False)
cov.to_csv(os.path.join(INTER, "P0_pseudobulk_coverage.csv"), index=False)
test_df.to_csv(os.path.join(INTER, "P0_pseudobulk_arm_tests.csv"), index=False)

def df_md(d, cols=None):
    if d is None or len(d) == 0:
        return "(empty)"
    cols = cols or list(d.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, rr in d.iterrows():
        vals = [f"{rr[c]:.4g}" if isinstance(rr[c], (int, float, np.floating)) else str(rr[c]) for c in cols]
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)

with open(REPORT, "w", encoding="utf-8") as f:
    f.write("# P0-5 供者级 Pseudobulk 重算报告（H4）\n\n")
    f.write("- 日期：2026-08-20；冻结依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md（OSF https://osf.io/C7RYD/）\n")
    f.write("- 单位：cell type × donor；供者=有效 n；细胞级统计仅描述\n")
    f.write(f"- 主对比：{CASE} vs {CTRL}（mild 不入确证）；分数据集独立检验\n")
    f.write(f"- 阈值（未注册参数，披露）：donor×celltype ≥{MIN_CELLS} 细胞；双侧各 ≥{MIN_DONORS} 供者\n\n")
    f.write("## 覆盖矩阵\n\n" + df_md(cov) + "\n\n")
    f.write("## 供者级臂评分检验\n\n")
    if len(test_df):
        f.write("| dataset | cell_type | score | n_case | n_ctrl | diff | Welch p | BH q |\n|---|---|---|---|---|---|---|---|\n")
        for _, r in test_df.iterrows():
            f.write(f"| {r['dataset']} | {r['cell_type']} | {r['score']} | {r['n_case']} | {r['n_ctrl']} | "
                    f"{r['diff']:+.3f} | {r['p']:.3g} | {r['BH_q']:.3g} |\n")
    f.write("\n## 基因级（pydeseq2）\n\n")
    f.write("- 已运行：见 _intermediate/P0_pseudobulk_DE_all.csv（全基因）与 P0_pseudobulk_DE_80genes.csv\n" if de_ok and de_all
            else "- 未运行/失败：按注册回退条款记录（limma-voom 需 R 环境，本环境不可用，作为已披露偏差）\n")
    f.write("\n## 与注册文本的偏差披露\n\n")
    f.write("1. limma-voom duplicateCorrelation 敏感性未运行（需 R；本环境仅 Python）——注册的回退条款针对 DESeq2 拟合失败，"
            "此处为环境限制，如实披露，不冒充已运行。\n")
    f.write("2. 合并（跨数据集 ~dataset+condition）敏感性本轮未跑（分数据集为主分析已覆盖注册主问题）。\n")
    f.write("3. MIN_CELLS/MIN_DONORS 阈值为未注册参数，取常规值并在报告与代码中披露。\n")

with open(LOGF, "w", encoding="utf-8") as f:
    f.write("\n".join(log))
print("\nDONE P0-5 ->", REPORT)
