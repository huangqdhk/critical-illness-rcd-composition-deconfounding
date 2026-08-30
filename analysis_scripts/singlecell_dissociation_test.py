#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
singlecell_dissociation_test.py
============================================================================
风险2 的核心检验：解离签名是否发生在【同一批细胞】里？

上游塌陷臂 = MAM_integrity + mitochondrial_function + ferroptosis_cuproptosis
执行诱导臂 = iron_metabolism + cell_death + oxidative_stress
（模块划分与 bulk Table_S5b 一致；neutral 模块 autophagy/transcription_factors 不入两臂）

判据（用户框架）：
  - 若【同细胞类型内】上游↔执行呈负相关（尤其 Mono/Macrophage）→ 真正的单细胞布线证据；
  - 若仅【全体混合】负相关、各细胞类型内不成立 → 分室化（composition artifact）。
无论哪种都是可发表结论，但写法不同。本脚本输出两者，让数据自陈。

三种读出（同一批细胞、同一基因集），互为敏感性：
  (A) score_genes 评分的 Spearman/Pearson 相关 —— 主结果（已减全局水平，去技术正相关偏倚）；
  (B) score_genes 残差化（逐细胞类型内回归掉 pct_counts_mt + log n_counts）后的相关 ——
      回答 mt-content 伪影问题；
  (C) 上游臂去掉 6 个 MT-* 基因（核编码子集）后的相关 —— 上游塌陷是否依赖 mt-gene。

输出（02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV、02_SUPPLEMENTARY_TABLES/Supplementary_Notes、01_FIGURE_DATA_CSV/Supplementary）：
  Table_S48_Dissociation_singlecell_correlation.csv   主表（pooled + per cell_type × condition，三读出）
  Table_S48_Dissociation_Method_Note.txt                            方法 + 结论
  Figure_S48_Dissociation_hexbin.png                各主要细胞类型 上游↔执行 密度散点
============================================================================
"""
import os, sys, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass

import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
GOV  = os.path.join(ROOT, "04_AUDIT_GOVERNANCE")
TAB  = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV")
NOTES = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES/Supplementary_Notes")
FIG  = os.path.join(ROOT, "01_FIGURE_DATA_CSV/Supplementary")
H5   = os.path.join(ROOT, "00_RAW_DATA", "_scrna_work", "combined_processed.h5ad")
GENELIST = os.path.join(GOV, "Mitoxyperilysis_Pathway_Gene_List.csv")
SEED = 0

# ---------- 两臂模块定义（与 bulk Table_S5b 一致）----------
UPSTREAM_MODS = ["mitoxy_MAM_integrity",
                 "mitoxy_mitochondrial_function",
                 "mitoxy_ferroptosis_cuproptosis"]
EXEC_MODS = ["mitoxy_iron_metabolism",
             "mitoxy_cell_death",
             "mitoxy_oxidative_stress"]
MT_PREFIX = "MT-"

print("载入数据 ...")
adata = ad.read_h5ad(H5)
# 确保 X 是 log 归一化（管道约定）；保留 raw 以防万一
print(f"  {adata.shape[0]} 细胞 × {adata.shape[1]} 基因")
print(f"  条件: {dict(adata.obs['condition'].value_counts())}")
print(f"  细胞类型: {dict(adata.obs['cell_type'].value_counts())}")

# ---------- 基因符号 → var 命中 ----------
auth = pd.read_csv(GENELIST)
varidx = set(adata.var_names)

# HGNC 别名映射（N3 修复，2026-08-15）：canonical 清单符号 IP3R1/HSP60 在
# h5ad 表达矩阵中以别名 ITPR1/HSPD1 存在。按 canonical manifest 的 aliases
# 列做矩阵符号映射，使上游臂 30/30 全部命中。
MANIFEST = os.path.join(GOV, "Mitoxyperilysis_Gene_Manifest_v1.0.csv")
if os.path.exists(MANIFEST):
    mf = pd.read_csv(MANIFEST)
    ALIAS_MAP = {r["gene_symbol"]: r["aliases"]
                 for _, r in mf.iterrows()
                 if pd.notna(r["aliases"]) and str(r["aliases"]).strip()}
else:
    ALIAS_MAP = {"IP3R1": "ITPR1", "HSP60": "HSPD1"}

def genes_of(modules):
    syms = auth.loc[auth["Mitoxyperilysis_module"].isin(modules), "gene_symbol"].tolist()
    present = []
    for s in syms:
        cands = [s] + [a.strip() for a in str(ALIAS_MAP.get(s, "")).split(",") if a.strip()]
        hit = next((c for c in cands if c in varidx), None)
        if hit is not None:
            present.append(hit)
    return syms, present

up_all_syms,  up_all  = genes_of(UPSTREAM_MODS)
ex_syms,      ex      = genes_of(EXEC_MODS)
up_nuc = [g for g in up_all if not g.startswith(MT_PREFIX)]          # 核编码上游（去 6 个 MT-*）
print(f"\n上游臂: {len(up_all)}/{len(up_all_syms)} 命中 | 核编码子集 {len(up_nuc)}")
print(f"执行臂: {len(ex)}/{len(ex_syms)} 命中")
print(f"  上游: {up_all}")
print(f"  执行: {ex}")

# ---------- 评分：score_genes（减全局水平）----------
print("\nscore_genes 评分 ...")
sc.tl.score_genes(adata, gene_list=up_all, score_name="upstream_score",
                  random_state=SEED, use_raw=False)
sc.tl.score_genes(adata, gene_list=ex,   score_name="execution_score",
                  random_state=SEED, use_raw=False)
sc.tl.score_genes(adata, gene_list=up_nuc, score_name="upstream_nuc_score",
                  random_state=SEED, use_raw=False)

obs = adata.obs[["cell_type", "condition", "pct_counts_mt", "n_counts",
                 "upstream_score", "execution_score", "upstream_nuc_score"]].copy()
obs["log_ncounts"] = np.log10(obs["n_counts"].clip(lower=1))

# ---------- 相关性辅助 ----------
def corr_block(df, a, b):
    """返回 n, pearson_r/p, spearman_r/p，跳过 NA。"""
    d = df[[a, b]].dropna()
    n = len(d)
    if n < 30:
        return dict(n=n, pear_r=np.nan, pear_p=np.nan, spear_r=np.nan, spear_p=np.nan)
    pr, pp = stats.pearsonr(d[a], d[b])
    sr, sp = stats.spearmanr(d[a], d[b])
    # bootstrap CI on spearman
    rng = np.random.default_rng(SEED)
    x, y = d[a].values, d[b].values
    bs = []
    for _ in range(500):
        idx = rng.integers(0, n, n)
        if np.unique(x[idx]).size > 1 and np.unique(y[idx]).size > 1:
            bs.append(stats.spearmanr(x[idx], y[idx])[0])
    ci = (float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))) if bs else (np.nan, np.nan)
    return dict(n=n, pear_r=float(pr), pear_p=float(pp),
                spear_r=float(sr), spear_p=float(sp), spear_ci_lo=ci[0], spear_ci_hi=ci[1])

def residualize(df, ycol, covs):
    """在 df 子集内，对 ycol 用 covs 线性回归取残差。"""
    d = df.copy()
    keep = [ycol] + covs
    dd = d[keep].dropna()
    X = np.column_stack([np.ones(len(dd))] + [dd[c].values for c in covs])
    coef, *_ = np.linalg.lstsq(X, dd[ycol].values, rcond=None)
    resid = dd[ycol].values - X @ coef
    out = d[ycol].copy()
    out.loc[dd.index] = resid
    return out

rows = []
def add(label, df, a, b, readout, note=""):
    r = corr_block(df, a, b)
    rows.append(dict(stratum=label, readout=readout, gene_a=a, gene_b=b,
                     n=r["n"], pearson_r=r["pear_r"], pearson_p=r["pear_p"],
                     spearman_r=r["spear_r"], spearman_p=r["spear_p"],
                     spear_ci_lo=r.get("spear_ci_lo", np.nan),
                     spear_ci_hi=r.get("spear_ci_hi", np.nan), note=note))

CELLTYPES = ["Mono_c14", "T_cell", "B_cell", "NK", "Epithelial", "DC", "Club", "Macrophage"]
CONDS = ["Healthy", "COVID_severe", "COVID_mild"]

# ============== 读出 (A) score_genes 原始相关 ==============
print("\n[A] score_genes 相关 ...")
add("ALL_pooled", obs, "upstream_score", "execution_score", "A_scoregenes")
for ct in CELLTYPES:
    add(f"{ct}__all", obs[obs.cell_type == ct], "upstream_score", "execution_score", "A_scoregenes")
    for cond in CONDS:
        sub = obs[(obs.cell_type == ct) & (obs.condition == cond)]
        if len(sub) >= 30:
            add(f"{ct}__{cond}", sub, "upstream_score", "execution_score", "A_scoregenes")

# ============== 读出 (B) 回归掉 pct_counts_mt + log n_counts ==============
print("[B] 回归 pct_counts_mt + log n_counts 后相关 ...")
obsB = obs.copy()
# 全局残差化（全体一次回归，保留可比性；再做 per-ct 版本）
obsB["up_resid_global"]   = residualize(obs, "upstream_score", ["pct_counts_mt", "log_ncounts"])
obsB["ex_resid_global"]   = residualize(obs, "execution_score", ["pct_counts_mt", "log_ncounts"])
add("ALL_pooled", obsB, "up_resid_global", "ex_resid_global", "B_mtregressed_global")
# per cell-type 残差化（更严格：每细胞类型内各自去 mt/confound）
obsB["up_resid_ct"]  = np.nan
obsB["ex_resid_ct"]  = np.nan
for ct in CELLTYPES:
    m = obs.cell_type == ct
    if m.sum() < 30:
        continue
    obsB.loc[m, "up_resid_ct"] = residualize(obs[m], "upstream_score", ["pct_counts_mt", "log_ncounts"]).values
    obsB.loc[m, "ex_resid_ct"] = residualize(obs[m], "execution_score", ["pct_counts_mt", "log_ncounts"]).values
    add(f"{ct}__all", obsB[m], "up_resid_ct", "ex_resid_ct", "B_mtregressed_perCT")
    for cond in CONDS:
        sub = obsB[(obsB.cell_type == ct) & (obsB.condition == cond)]
        if len(sub) >= 30:
            add(f"{ct}__{cond}", sub, "up_resid_ct", "ex_resid_ct", "B_mtregressed_perCT")

# ============== 读出 (C) 核编码上游子集（去 6 MT-*）==============
print("[C] 核编码上游（去 MT-*）相关 ...")
add("ALL_pooled", obs, "upstream_nuc_score", "execution_score", "C_nuclear_upstream")
for ct in CELLTYPES:
    sub = obs[obs.cell_type == ct]
    if len(sub) >= 30:
        add(f"{ct}__all", sub, "upstream_nuc_score", "execution_score", "C_nuclear_upstream")
    for cond in CONDS:
        s2 = obs[(obs.cell_type == ct) & (obs.condition == cond)]
        if len(s2) >= 30:
            add(f"{ct}__{cond}", s2, "upstream_nuc_score", "execution_score", "C_nuclear_upstream")

res = pd.DataFrame(rows)
res.to_csv(os.path.join(TAB, "Table_S48_Dissociation_singlecell_correlation.csv"),
           index=False, encoding="utf-8-sig")
print("\n[已保存] Table_S48_Dissociation_singlecell_correlation.csv")

# ============== 控制台摘要 ==============
def fmt(r):
    if np.isnan(r): return "   —   "
    return f"{r:+.3f}"

print("\n" + "=" * 96)
print("主结果：上游 ↔ 执行 同细胞相关（Spearman ρ）")
print("=" * 96)
for rd, title in [("A_scoregenes", "[A] score_genes"),
                  ("C_nuclear_upstream", "[C] 核编码上游(去MT-*)"),
                  ("B_mtregressed_perCT", "[B] 回归mt+counts(逐CT)")]:
    sub = res[res.readout == rd]
    print(f"\n{title}")
    print(f"  {'stratum':28s} {'n':>8s} {'ρ':>8s} {'95%CI':>20s} {'pearson':>9s}")
    pool = sub[sub.stratum == "ALL_pooled"]
    if len(pool):
        p = pool.iloc[0]
        print(f"  {'ALL_pooled':28s} {int(p.n):8d} {fmt(p.spearman_r):>8s} "
              f"[{fmt(p.spear_ci_lo)},{fmt(p.spear_ci_hi)}] {fmt(p.pearson_r):>9s}")
    for ct in CELLTYPES:
        row = sub[sub.stratum == f"{ct}__all"]
        if len(row):
            r = row.iloc[0]
            print(f"  {ct+'__all':28s} {int(r.n):8d} {fmt(r.spearman_r):>8s} "
                  f"[{fmt(r.spear_ci_lo)},{fmt(r.spear_ci_hi)}] {fmt(r.pearson_r):>9s}")

# 关键细胞类型 × 条件（Mono / Macrophage 在 Healthy vs severe）
print("\n关键细胞类型 × 条件（Mono_c14 / Macrophage / Epithelial）：")
for ct in ["Mono_c14", "Macrophage", "Epithelial"]:
    print(f"  [{ct}] (score_genes ρ)")
    for cond in CONDS:
        row = res[(res.readout == "A_scoregenes") & (res.stratum == f"{ct}__{cond}")]
        if len(row):
            r = row.iloc[0]
            print(f"     {cond:14s} n={int(r.n):7d}  ρ={fmt(r.spearman_r)}  "
                  f"CI[{fmt(r.spear_ci_lo)},{fmt(r.spear_ci_hi)}]")

# ============== 图：主要细胞类型 上游↔执行 密度散点 ==============
print("\n绘图 ...")
plot_cts = ["Mono_c14", "T_cell", "NK", "Epithelial", "DC", "Macrophage"]
fig, axes = plt.subplots(2, 3, figsize=(13, 8))
for ax, ct in zip(axes.ravel(), plot_cts):
    d = obs[obs.cell_type == ct][["upstream_score", "execution_score", "condition"]].dropna()
    if len(d) < 50:
        ax.text(0.5, 0.5, f"{ct}\nn={len(d)}", ha="center", va="center", transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        continue
    sev = d[d.condition == "COVID_severe"]
    hlth = d[d.condition == "Healthy"]
    ax.hexbin(sev["upstream_score"], sev["execution_score"], gridsize=40,
              cmap="Reds", mincnt=1, alpha=0.85, label="severe")
    if len(hlth):
        ax.scatter(hlth["upstream_score"], hlth["execution_score"], s=2,
                   c="steelblue", alpha=0.10, label="healthy", rasterized=True)
    rho_s = stats.spearmanr(d["upstream_score"], d["execution_score"])[0]
    rho_sev = stats.spearmanr(sev["upstream_score"], sev["execution_score"])[0] if len(sev) > 30 else np.nan
    ax.set_title(f"{ct} (n={len(d)})\nρ_all={rho_s:+.2f}  ρ_severe={rho_sev:+.2f}", fontsize=10)
    ax.set_xlabel("upstream score"); ax.set_ylabel("execution score")
    ax.axhline(0, color="grey", lw=0.5); ax.axvline(0, color="grey", lw=0.5)
fig.suptitle("Within-cell upstream vs execution (score_genes)\nRed = COVID_severe density, blue = Healthy",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(os.path.join(FIG, "Figure_S48_Dissociation_hexbin.png"), dpi=150)
print("[已保存] Figure_S48_Dissociation_hexbin.png")

# ============== 方法/结论 note ==============
def pick(rd, strat):
    r = res[(res.readout == rd) & (res.stratum == strat)]
    return r.iloc[0] if len(r) else None

note = []
note.append("Table S48 单细胞解离检验 — 方法与结论")
note.append("=" * 70)
note.append("数据: combined_processed.h5ad, 138,941 细胞 (GSE145926 + GSE158055)")
note.append("条件: Healthy 44,215 / COVID_mild 29,754 / COVID_severe 64,972")
note.append("评分: scanpy score_genes（基因集均值 − 表达匹配的对照基因集均值），")
note.append("      use_raw=False（已 log 归一化 X）。score_genes 减去全局水平，")
note.append("      消除『高 counts 细胞→两臂都高→假性正相关』的技术偏倚。")
note.append("")
note.append("两臂（模块划分与 bulk Table_S5b 一致）：")
note.append(f"  上游塌陷臂 ({len(up_all)} 基因): MAM_integrity + mitochondrial_function + ferroptosis_cuproptosis")
note.append(f"  执行诱导臂 ({len(ex)} 基因): iron_metabolism + cell_death + oxidative_stress")
note.append(f"  核编码上游子集 ({len(up_nuc)} 基因): 上游臂去掉 6 个 MT-* 基因")
note.append("")
note.append("判据: 同细胞类型内 ρ<0 → 真单细胞布线（解离成立）；")
note.append("      仅全体 ρ<0 而各 CT 内 ≈0 → 分室化（组成假象）。")
note.append("")
note.append("-" * 70)
p = pick("A_scoregenes", "ALL_pooled")
note.append(f"全体混合 ρ = {p.spearman_r:+.3f}" if p is not None else "全体混合 ρ = NA")
note.append("-" * 70)
note.append("各细胞类型内 ρ（score_genes，全条件）：")
for ct in CELLTYPES:
    r = pick("A_scoregenes", f"{ct}__all")
    if r is not None:
        note.append(f"  {ct:12s} n={int(r.n):7d}  ρ={r.spearman_r:+.3f} "
                    f"[{r.spear_ci_lo:+.3f},{r.spear_ci_hi:+.3f}]  pearson={r.pearson_r:+.3f}")
note.append("")
note.append("敏感性（核编码上游 / mt-回归）见 Table_S48 ..._correlation.csv 的 readout B、C。")
with open(os.path.join(NOTES, "Table_S48_Dissociation_Method_Note.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(note))
print("\n[已保存] Table_S48_Dissociation_Method_Note.txt")
print("\n[完成]")
