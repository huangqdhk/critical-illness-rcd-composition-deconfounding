#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
S48_donor_cluster_bootstrap.py
============================================================================
§3.1 问题表 #2 的补算：单细胞解离相关的【供者 cluster bootstrap】重算。

背景：Table S48 三读出的 ρ（−0.116 / +0.01 / +0.079）以细胞为统计单位
（138,941 细胞），细胞间不独立（同供者细胞相关），naive p 严重低估不确定性。
本脚本把推断单位升到供者：有放回重抽供者簇（GSE145926 按样本=供者 12 例；
GSE158055 按 subject 40 例，重复采样 -1/-2 合并到供者），每次重抽后合并
其全部细胞重算 Spearman ρ，取 2.5/97.5 分位为供者级 95% CI，双侧 p 由
bootstrap 分布零点两侧比例给出。

输出：
  02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV/Table_S48b_Dissociation_donor_bootstrap.csv
  （并追加方法说明至 Table_S48_Dissociation_Method_Note.txt）
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

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
GOV  = os.path.join(ROOT, "04_AUDIT_GOVERNANCE")
TAB  = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV")
NOTES = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES/Supplementary_Notes")
LOG  = os.path.join(ROOT, "03_LOGS", "S48_donor_bootstrap_log.txt")
H5   = os.path.join(ROOT, "00_RAW_DATA", "_scrna_work", "combined_processed.h5ad")
GENELIST = os.path.join(GOV, "Mitoxyperilysis_Pathway_Gene_List.csv")
MANIFEST = os.path.join(GOV, "Mitoxyperilysis_Gene_Manifest_v1.0.csv")
SEED = 0
B    = 2000   # bootstrap 次数

UPSTREAM_MODS = ["mitoxy_MAM_integrity",
                 "mitoxy_mitochondrial_function",
                 "mitoxy_ferroptosis_cuproptosis"]
EXEC_MODS = ["mitoxy_iron_metabolism",
             "mitoxy_cell_death",
             "mitoxy_oxidative_stress"]
MT_PREFIX = "MT-"
CELLTYPES = ["Mono_c14", "T_cell", "B_cell", "NK", "Epithelial", "DC", "Club", "Macrophage"]

_log = []
def say(msg):
    print(msg, flush=True)
    _log.append(msg)

say("载入数据 ...")
adata = ad.read_h5ad(H5)
say(f"  {adata.shape[0]} 细胞 × {adata.shape[1]} 基因")

# ---------- 基因命中（与 singlecell_dissociation_test.py 完全一致）----------
auth = pd.read_csv(GENELIST)
varidx = set(adata.var_names)
mf = pd.read_csv(MANIFEST)
ALIAS_MAP = {r["gene_symbol"]: r["aliases"] for _, r in mf.iterrows()
             if pd.notna(r["aliases"]) and str(r["aliases"]).strip()}

def genes_of(modules):
    syms = auth.loc[auth["Mitoxyperilysis_module"].isin(modules), "gene_symbol"].tolist()
    present = []
    for s in syms:
        cands = [s] + [a.strip() for a in str(ALIAS_MAP.get(s, "")).split(",") if a.strip()]
        hit = next((c for c in cands if c in varidx), None)
        if hit is not None:
            present.append(hit)
    return syms, present

up_syms, up = genes_of(UPSTREAM_MODS)
ex_syms, ex = genes_of(EXEC_MODS)
up_nuc = [g for g in up if not g.startswith(MT_PREFIX)]
say(f"上游臂 {len(up)}/{len(up_syms)} 命中（核编码 {len(up_nuc)}）| 执行臂 {len(ex)}/{len(ex_syms)} 命中")

# ---------- 评分（与原脚本一致：score_genes, seed=0）----------
say("score_genes 评分 ...")
sc.tl.score_genes(adata, gene_list=up,     score_name="upstream_score",    random_state=SEED, use_raw=False)
sc.tl.score_genes(adata, gene_list=ex,     score_name="execution_score",   random_state=SEED, use_raw=False)
sc.tl.score_genes(adata, gene_list=up_nuc, score_name="upstream_nuc_score", random_state=SEED, use_raw=False)

obs = adata.obs[["cell_type", "condition", "dataset", "sampleID", "subject",
                 "pct_counts_mt", "n_counts",
                 "upstream_score", "execution_score", "upstream_nuc_score"]].copy()
obs["log_ncounts"] = np.log10(obs["n_counts"].clip(lower=1))

# ---------- 供者簇键 ----------
# GSE145926: sampleID（12 例，每样本一供者）；GSE158055: subject（40 例，-1/-2 重复采样并入供者）
donor = []
for ds, sid, su in zip(obs["dataset"].astype(str), obs["sampleID"].astype(str), obs["subject"].astype(str)):
    if ds == "GSE145926":
        donor.append("G1_" + sid)
    else:
        donor.append("G2_" + su)
obs["donor"] = donor
donors = sorted(obs["donor"].unique())
say(f"供者簇 {len(donors)} 个（GSE145926 {sum(d.startswith('G1_') for d in donors)} + "
    f"GSE158055 {sum(d.startswith('G2_') for d in donors)}）")

# ---------- 残差化（读出 B，全局口径与原脚本一致）----------
def residualize(df, ycol, covs):
    d = df.copy()
    keep = [ycol] + covs
    dd = d[keep].dropna()
    X = np.column_stack([np.ones(len(dd))] + [dd[c].values for c in covs])
    coef, *_ = np.linalg.lstsq(X, dd[ycol].values, rcond=None)
    resid = dd[ycol].values - X @ coef
    out = d[ycol].copy()
    out.loc[dd.index] = resid
    return out

obs["up_resid"] = residualize(obs, "upstream_score", ["pct_counts_mt", "log_ncounts"])
obs["ex_resid"] = residualize(obs, "execution_score", ["pct_counts_mt", "log_ncounts"])

# ---------- 供者 cluster bootstrap ----------
def donor_bootstrap(df, acol, bcol, min_donors=8):
    """有放回重抽供者簇 → 合并细胞重算 Spearman ρ。返回 (点估计, naive_p, ci_lo, ci_hi, p_bs, n_donors, n_cells)。"""
    d = df[[acol, bcol, "donor"]].dropna()
    groups = [g[[acol, bcol]].to_numpy() for _, g in d.groupby("donor", observed=True)]
    k = len(groups)
    if k < min_donors:
        return (np.nan,) * 5 + (k, len(d))
    x = d[acol].to_numpy(); y = d[bcol].to_numpy()
    rho0 = stats.spearmanr(x, y)[0]
    naive_p = stats.spearmanr(x, y)[1]
    rng = np.random.default_rng(SEED)
    boots = np.empty(B)
    for i in range(B):
        idx = rng.integers(0, k, k)
        xy = np.concatenate([groups[j] for j in idx])
        boots[i] = stats.spearmanr(xy[:, 0], xy[:, 1])[0]
    ci = (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5)))
    # 双侧 bootstrap p（含 +1 校正）
    n_ge = int(np.sum(boots >= 0)); n_le = int(np.sum(boots <= 0))
    p_bs = min(1.0, 2.0 * min((n_ge + 1) / (B + 1), (n_le + 1) / (B + 1)))
    return (float(rho0), float(naive_p), ci[0], ci[1], p_bs, k, len(d))

rows = []
def add(stratum, readout, df, acol, bcol, note=""):
    r = donor_bootstrap(df, acol, bcol)
    rows.append(dict(stratum=stratum, readout=readout,
                     rho_point=r[0], naive_cell_p=r[1],
                     donor_bs_ci_lo=r[2], donor_bs_ci_hi=r[3], donor_bs_p=r[4],
                     n_donors=r[5], n_cells=r[6], B=B, note=note))
    say(f"  {stratum:28s} [{readout:16s}] ρ={r[0]:+.4f}  供者CI[{r[2]:+.4f},{r[3]:+.4f}]  "
        f"p_bs={r[4]:.4f}  供者={r[5]}  细胞={r[6]}")

say("\n== 供者 cluster bootstrap（主读出）==")
say("[A] score_genes 上游 vs 执行：")
add("ALL_pooled", "A_scoregenes", obs, "upstream_score", "execution_score")
for ct in CELLTYPES:
    sub = obs[obs.cell_type == ct]
    if len(sub) >= 30:
        add(f"{ct}__all", "A_scoregenes", sub, "upstream_score", "execution_score")

say("[B] 全局残差化（mt+counts）上游 vs 执行：")
add("ALL_pooled", "B_mtregressed_global", obs, "up_resid", "ex_resid")

say("[C] 核编码上游（去 MT-*）vs 执行：")
add("ALL_pooled", "C_nuclear_upstream", obs, "upstream_nuc_score", "execution_score")

say("[A×cond] CD14+ 单核细胞 Healthy vs COVID_severe：")
for cond in ["Healthy", "COVID_severe"]:
    sub = obs[(obs.cell_type == "Mono_c14") & (obs.condition == cond)]
    if len(sub) >= 30:
        add(f"Mono_c14__{cond}", "A_scoregenes", sub, "upstream_score", "execution_score")

res = pd.DataFrame(rows)
out_csv = os.path.join(TAB, "Table_S48b_Dissociation_donor_bootstrap.csv")
res.to_csv(out_csv, index=False, encoding="utf-8-sig")
say(f"\n[已保存] {out_csv}")

# ---------- 与 S48 冻结值交叉核对 ----------
say("\n== 与 Table S48 冻结值核对（点估计）==")
frozen = {"ALL_pooled|A_scoregenes": -0.116, "ALL_pooled|B_mtregressed_global": 0.01,
          "ALL_pooled|C_nuclear_upstream": 0.079}
for key, val in frozen.items():
    st, rd = key.split("|")
    row = res[(res.stratum == st) & (res.readout == rd)]
    if len(row):
        r0 = row.iloc[0]["rho_point"]
        flag = "OK" if abs(r0 - val) < 0.005 else "偏离>0.005，需查"
        say(f"  {key:44s} 冻结={val:+.3f}  重算={r0:+.4f}  {flag}")

# ---------- 追加方法说明 ----------
note_path = os.path.join(NOTES, "Table_S48_Dissociation_Method_Note.txt")
append = [
    "",
    "=" * 70,
    "供者 cluster bootstrap（Table S48b，2026-08-23 补算）",
    "=" * 70,
    "动机: 上述 ρ 以细胞为单位（138,941 细胞），同供者细胞不独立，",
    "      naive p 高估置信度。本补算把推断单位升到供者：有放回重抽",
    "      供者簇（B=2000, seed=0），合并所抽供者全部细胞重算 Spearman ρ。",
    "供者簇: 52 个 = GSE145926 按样本 12 例 + GSE158055 按 subject 40 例",
    "      （-1/-2 重复采样并入同一供者）。",
    "读出: 点估计与 naive p 复刻 Table S48；新增 供者级 95% CI 与双侧 p。",
    "",
    "主结果（供者级推断）:",
]
for _, r in res.iterrows():
    append.append(f"  {r['stratum']:28s} [{r['readout']:20s}] ρ={r['rho_point']:+.4f}  "
                  f"供者CI[{r['donor_bs_ci_lo']:+.4f},{r['donor_bs_ci_hi']:+.4f}]  "
                  f"p_donor={r['donor_bs_p']:.4f}  供者={int(r['n_donors'])}")
append.append("")
append.append("结论: 细胞级 ρ 仅作描述性呈现；显著性一律以供者级 bootstrap 为准。")
with open(note_path, "a", encoding="utf-8") as f:
    f.write("\n" + "\n".join(append) + "\n")
say(f"[已追加] {note_path}")

with open(LOG, "w", encoding="utf-8") as f:
    f.write("\n".join(_log) + "\n")
say(f"[已保存] {LOG}")
say("\n[完成]")
