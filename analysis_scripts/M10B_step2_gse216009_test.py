# -*- coding: utf-8 -*-
"""
M10B_step2_gse216009_test.py — GSE216009 全血 scRNA 供者×细胞型 arm 评分（M10 检验 B 队列）
========================================================================================
数据：R 提取件（meta.csv + data_sub.csv；Seurat data 槽 = log1p(CPM*1e4) 逐细胞）
口径（沿用 P0 供者级 pseudobulk）：
  - donor x celltype >= 50 细胞、双侧各 >= 3 供者
  - arm 评分 = 供者×细胞型内 arm 基因 log1p 均值 -> 细胞型内 z -> Welch t（BH 族 = celltype x score）
  - 主对比：Sepsis（急性，26 例）vs CS+HV（13 例）；Sepsis_conv 复采样本不入主对比
输出：_intermediate/M10B_gse216009_arm_scores.csv / M10B_gse216009_arm_tests.csv
"""
import os, sys, time, warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
INTER = ROOT + r"\_intermediate"
MANIFEST = ROOT + r"\04_AUDIT_GOVERNANCE\Mitoxyperilysis_Gene_Manifest_v1.0.csv"
EXTRACT = INTER + r"\gse216009_out"   # R 提取件（M10B_step0_extract_gse216009.R 产物，经 E:\SCI\_rtmp 中转后归档）
MT_GENES = {"MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB",
            "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6"}

sys.path.insert(0, ROOT)

t0 = time.time()
log = []
def note(m=""):
    log.append(m); print(m, flush=True)

def zc(s):
    s = np.asarray(s, float); sd = s.std(ddof=1)
    return (s - s.mean()) / sd if sd > 0 else s * 0

def bh(pvals):
    p = np.asarray(pvals, float); out = np.full_like(p, np.nan)
    m = ~np.isnan(p); ps = p[m]
    if len(ps) == 0: return out
    order = np.argsort(ps); n = len(ps); adj = np.empty(n)
    for k in range(n - 1, -1, -1):
        adj[order[k]] = min(ps[order[k]] * n / (k + 1), adj[order[k + 1]] if k + 1 < n else 1.0)
    out[m] = adj; return out

note("== 1. 读入提取件")
meta = pd.read_csv(EXTRACT + r"\meta.csv", index_col=0)
meta.index = [str(i) for i in meta.index]
man = pd.read_csv(MANIFEST)
up_genes = set(man.loc[man["arm"] == "upstream_collapse", "hgnc_symbol"].astype(str).str.upper())
ex_genes = set(man.loc[man["arm"] == "execution_induction", "hgnc_symbol"].astype(str).str.upper())
note(f"   meta {meta.shape}；arm 上游 {len(up_genes)} / 执行 {len(ex_genes)}")

data = pd.read_csv(EXTRACT + r"\data_sub.csv", index_col=0).T    # cells x genes
data.index = [str(i) for i in data.index]
data.columns = [str(c).upper() for c in data.columns]
note(f"   data_sub {data.shape}")
common = sorted(set(data.index) & set(meta.index))
data = data.loc[common]; meta = meta.loc[common]
up_in = [g for g in data.columns if g in up_genes]
ex_in = [g for g in data.columns if g in ex_genes]
up_nomt = [g for g in up_in if g not in MT_GENES]
note(f"   可测 arm 基因：上游 {len(up_in)}/{len(up_genes)}（去MT {len(up_nomt)}），执行 {len(ex_in)}/{len(ex_genes)}")

note("\n== 2. 供者×细胞型聚合")
MIN_CELLS, MIN_DONORS = 50, 3
meta["source"] = meta["source"].astype(str)
# 主对比：Sepsis（急性）vs CS+HV；Sepsis_conv 排除
keep = meta["source"].isin(["Sepsis", "CS", "HV"])
m2 = meta[keep]
d2 = data.loc[m2.index]
cond_map = {"Sepsis": "Sepsis", "CS": "Control", "HV": "Control"}
d2["cond"] = [cond_map[s] for s in m2["source"]]
d2["sample"] = m2["sample_id"].astype(str)
d2["ct"] = m2["fine_annot"].astype(str)
recs = []
for (sm, ct), sub in d2.groupby(["sample", "ct"]):
    n = len(sub)
    if n < MIN_CELLS:
        continue
    mu = sub[up_in + ex_in].mean(axis=0)
    recs.append(dict(cohort="GSE216009", sampleID=sm, cell_type=ct,
                     condition=sub["cond"].iloc[0], n_cells=n,
                     UCS_raw=float(mu[up_in].mean()),
                     EIS_raw=float(mu[ex_in].mean()),
                     UCS_nomt_raw=float(mu[up_nomt].mean())))
df = pd.DataFrame(recs)
note(f"   供者×细胞型记录 {len(df)}")

note("\n== 3. 细胞型内 z + Welch t")
df["UCS"] = df["EIS"] = df["MDI"] = df["UCS_nomt"] = df["MDI_nomt"] = np.nan
for ct in df["cell_type"].unique():
    m = df["cell_type"] == ct
    df.loc[m, "UCS"] = zc(df.loc[m, "UCS_raw"])
    df.loc[m, "EIS"] = zc(df.loc[m, "EIS_raw"])
    df.loc[m, "MDI"] = df.loc[m, "EIS"] - df.loc[m, "UCS"]
    df.loc[m, "UCS_nomt"] = zc(df.loc[m, "UCS_nomt_raw"])
    df.loc[m, "MDI_nomt"] = df.loc[m, "EIS"] - df.loc[m, "UCS_nomt"]
from scipy.stats import ttest_ind
tests = []
for ct in df["cell_type"].unique():
    sub = df[df["cell_type"] == ct]
    for score in ("UCS", "EIS", "MDI", "MDI_nomt"):
        a_ = sub.loc[sub["condition"] == "Sepsis", score].values
        b_ = sub.loc[sub["condition"] == "Control", score].values
        if len(a_) >= MIN_DONORS and len(b_) >= MIN_DONORS:
            t, p = ttest_ind(a_, b_, equal_var=False)
            tests.append(dict(cohort="GSE216009", cell_type=ct, score=score,
                              n_case=len(a_), n_ctrl=len(b_),
                              diff=float(a_.mean() - b_.mean()), p=float(p)))
tests = pd.DataFrame(tests)
tests["BH_q"] = bh(tests["p"].values) if len(tests) else np.nan
note("   检验结果：")
for _, r in tests.iterrows():
    note(f"     {r['cell_type']:<36} {r['score']:<9} n={r['n_case']}v{r['n_ctrl']} "
         f"diff={r['diff']:+.3f} p={r['p']:.3g} q={r['BH_q']:.3g}")

df.to_csv(INTER + r"\M10B_gse216009_arm_scores.csv", index=False)
tests.to_csv(INTER + r"\M10B_gse216009_arm_tests.csv", index=False)
with open(ROOT + r"\03_LOGS\M10B_gse216009_log.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(log))
note(f"\nDONE in {time.time()-t0:.1f}s")
