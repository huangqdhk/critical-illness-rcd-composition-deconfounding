#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N3 衍生修复（2026-08-15）：真·单细胞 DEG 重算
============================================================================
背景：Table_S2_DEGs_Analysis.csv 的 ARDS_vs_Control / Sepsis_vs_Control 列经
数值指纹核对（ACTB≈13.0、EPCAM≈0.5、样本设计 82/44、348/44）实为
GSE185263 全血 Bulk DESeq2 结果，此前被主稿与 Bridge 收敛表误标为
"单细胞BALF"（sc_log2FC）。5.1 管线从未产出过单细胞 DEG 表。

本脚本在修正标签后的 combined_processed.h5ad 上重算 80 通路基因的
单细胞差异表达（COVID_severe vs Healthy）：
  - 分数据集（GSE145926 BALF / GSE158055 PBMC）+ 合并（Combined）
  - 输出 log1p 归一化表达均值、Δmean（自然对数尺度）、表达比例、
    Mann-Whitney p、BH padj（80 基因内校正）
输出：02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV/Table_S2b_scRNA_DEG_80genes.csv
============================================================================
"""
import os, sys, warnings
import numpy as np, pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu
warnings.filterwarnings("ignore")
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass

import scanpy as sc

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
T = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV")
H5 = os.path.join(ROOT, "00_RAW_DATA", "_scrna_work", "combined_processed.h5ad")

adata = sc.read_h5ad(H5)
obs = adata.obs
print(f"cells: {adata.shape[0]}, cond: {dict(obs['condition'].value_counts())}")

# 80 基因清单（权威符号）+ manifest 别名映射
auth = pd.read_csv(os.path.join(T, "Mitoxyperilysis_Pathway_Gene_List.csv"))
manifest = pd.read_csv(os.path.join(T, "Mitoxyperilysis_Gene_Manifest_v1.0.csv"))
alias = {}
for _, r in manifest.iterrows():
    if pd.notna(r["aliases"]):
        for al in str(r["aliases"]).split(","):
            alias[al.strip()] = r["gene_symbol"]

varidx = list(adata.var_names)
sym2idx = {g: i for i, g in enumerate(varidx)}

def bh(pvals):
    """Benjamini-Hochberg FDR（保留 NaN）。"""
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    m = ~np.isnan(p)
    ps = p[m]
    order = np.argsort(ps)
    n = len(ps)
    adj = np.empty(n)
    for k in range(n - 1, -1, -1):
        adj[order[k]] = min(ps[order[k]] * n / (k + 1), adj[order[k + 1]] if k + 1 < n else 1.0)
    out[m] = adj
    return out

def deg_block(X, sev_mask, syms, hits):
    rows = []
    pvals = []
    for sym, hit in zip(syms, hits):
        if hit is None:
            rows.append(None); pvals.append(np.nan); continue
        j = sym2idx[hit]
        col = X[:, j]
        v_sev = np.asarray(col[sev_mask].toarray()).ravel() if sp.issparse(col) else np.asarray(col[sev_mask]).ravel()
        v_hlt = np.asarray(col[~sev_mask].toarray()).ravel() if sp.issparse(col) else np.asarray(col[~sev_mask]).ravel()
        u, p = mannwhitneyu(v_sev, v_hlt, alternative="two-sided")
        rows.append(dict(sev_mean=v_sev.mean(), hlt_mean=v_hlt.mean(),
                         sev_pct=(v_sev > 0).mean(), hlt_pct=(v_hlt > 0).mean(), mw_p=p))
        pvals.append(p)
    padj = bh(pvals)
    out = []
    for r, pa in zip(rows, padj):
        if r is None:
            out.append(None)
        else:
            r["bh_padj"] = pa
            out.append(r)
    return out

syms = auth["gene_symbol"].tolist()
mods = auth["Mitoxyperilysis_module"].tolist()
hits = []
for s in syms:
    cands = [s]
    al = str(manifest.loc[manifest.gene_symbol == s, "aliases"].values[0]) \
         if (manifest.gene_symbol == s).any() else ""
    cands += [x.strip() for x in al.split(",") if x.strip()]
    hits.append(next((c for c in cands if c in sym2idx), None))
print(f"命中 {sum(h is not None for h in hits)}/{len(syms)}")

# 每组用布尔掩码
cond = obs["condition"].values
ds = obs["dataset"].values
sev_all = (cond == "COVID_severe")
hlt_all = (cond == "Healthy")

base = pd.DataFrame(dict(gene_symbol=syms, Mitoxyperilysis_module=mods, hit=hits))

results = {}
for label, m_sev, m_hlt in [
    ("GSE145926_BALF", sev_all & (ds == "GSE145926"), hlt_all & (ds == "GSE145926")),
    ("GSE158055_PBMC", sev_all & (ds == "GSE158055"), hlt_all & (ds == "GSE158055")),
    ("Combined", sev_all, hlt_all),
]:
    print(f"\n[{label}] sev={m_sev.sum()} hlt={m_hlt.sum()}")
    blk = deg_block(adata.X, m_sev, syms, hits)
    results[label] = blk
    # 模块级摘要
    df = pd.DataFrame([dict(sym=s, mod=m, **b) for s, m, b in zip(syms, mods, blk) if b])
    df["d_mean"] = df["sev_mean"] - df["hlt_mean"]
    print("  模块级 Δmean 均值：")
    for mod, g in df.groupby("mod"):
        print(f"    {mod:38s} n={len(g):2d}  mean d={g.d_mean.mean():+.3f}  up率={ (g.d_mean>0).mean():.2f}")
    # 显著数（BH<0.05）
    nsig = (df.bh_padj < 0.05).sum()
    print(f"    BH<0.05: {nsig}/80")

# 汇总输出
cols = []
for s in syms:
    row = dict(gene_symbol=s, Mitoxyperilysis_module=mods[syms.index(s)], hit=hits[syms.index(s)])
    for label, blk in results.items():
        r = blk[syms.index(s)]
        pfx = label
        if r is None:
            row.update({f"{pfx}_sev_mean": np.nan, f"{pfx}_hlt_mean": np.nan,
                        f"{pfx}_d_mean": np.nan, f"{pfx}_sev_pct": np.nan,
                        f"{pfx}_hlt_pct": np.nan, f"{pfx}_mw_p": np.nan, f"{pfx}_bh_padj": np.nan})
        else:
            row.update({f"{pfx}_sev_mean": r["sev_mean"], f"{pfx}_hlt_mean": r["hlt_mean"],
                        f"{pfx}_d_mean": r["sev_mean"] - r["hlt_mean"],
                        f"{pfx}_sev_pct": r["sev_pct"], f"{pfx}_hlt_pct": r["hlt_pct"],
                        f"{pfx}_mw_p": r["mw_p"], f"{pfx}_bh_padj": r["bh_padj"]})
    cols.append(row)
out = pd.DataFrame(cols)
out.to_csv(os.path.join(T, "Table_S2b_scRNA_DEG_80genes.csv"), index=False, encoding="utf-8-sig")
print(f"\n[已保存] Table_S2b_scRNA_DEG_80genes.csv ({len(out)} 行)")

# 关键基因打印
print("\n关键基因（d = sev_mean − hlt_mean, log1p 尺度）：")
print(f"{'gene':10s} | {'BALF d':>8s} {'padj':>9s} | {'PBMC d':>8s} {'padj':>9s} | {'Comb d':>8s} {'padj':>9s}")
for g in ["MT-ATP8", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ATP6", "MT-CO1", "TOMM20", "TOMM40",
          "IP3R1", "HSPA9", "VDAC1", "CANX", "CALR", "PACS2", "SIGMAR1", "VDAC2", "VDAC3",
          "SLC40A1", "FTL", "FTH1", "TFRC", "SLC11A2", "ACO1", "SLC25A37", "IREB2",
          "CASP5", "GSDMD", "CASP1", "NLRP3", "IL1B", "GPX4", "SOD2", "TXN", "CAT",
          "NFE2L2", "HIF1A", "ATF4", "DDIT3", "SQSTM1", "MAP1LC3B"]:
    r = out[out.gene_symbol == g]
    if len(r) == 0 or pd.isna(r.iloc[0]["Combined_d_mean"]):
        continue
    r = r.iloc[0]
    print(f"{g:10s} | {r.GSE145926_BALF_d_mean:+8.3f} {r.GSE145926_BALF_bh_padj:9.1e} | "
          f"{r.GSE158055_PBMC_d_mean:+8.3f} {r.GSE158055_PBMC_bh_padj:9.1e} | "
          f"{r.Combined_d_mean:+8.3f} {r.Combined_bh_padj:9.1e}")
