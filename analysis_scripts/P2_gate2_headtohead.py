# -*- coding: utf-8 -*-
"""
P2_gate2_headtohead.py — Phase 2 任务2：IIAMD 迁移验证 + 七面板头对头（H7，Gate 2）
====================================================================================
依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md（OSF https://osf.io/C7RYD/）§4 H7、§5 Gate 2
Phase A（bulk 验证层，GSE32707 = 唯一确认 Berlin-defined ARDS 全血队列）：
  - 基因符号 log2 矩阵（本地）→ 每基因集 arm-mean → 队列内全样本 z（不分组，防泄漏）
  - 主对比：ARDS_d0 vs Control（患者级线性模型 + 髓系标志物组成调整）
  - 基因集：IIAMD-up / IIAMD-down / IIAMD-diff(=z_up − z_down) / 七竞争面板 / 锚定面板
  - BH 族 = 基因集 × 1 主对比；头对头：IIAMD-diff 的 g 不得劣于最优竞争面板
Phase B（scRNA 髓系定位层，GSE145926 BALF + GSE158055 PBMC，供者级）：
  - donor×celltype log1p 均值 → celltype 内 z → severe vs Healthy（同 P0-5 单位规则）
输出：_intermediate/P2_gate2_scores_GSE32707.csv / P2_gate2_scRNA_localization.csv
     04_AUDIT_GOVERNANCE/P2_Gate2_HeadToHead_Report.md
"""
import os, sys, warnings
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
warnings.filterwarnings("ignore")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
RAW = ROOT + r"\00_RAW_DATA"
INTER = ROOT + r"\_intermediate"
SIG = os.path.join(INTER, "P2_IIAMD_human_v1.0.csv")
PANELS = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "P2_competing_panels_v1.0.csv")
REPORT = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "P2_Gate2_HeadToHead_Report.md")
LOGF = os.path.join(ROOT, "03_LOGS", "P2_gate2_log.txt")

sys.path.insert(0, ROOT)
import mdi_lib as L

log = []
def note(m=""):
    log.append(m); print(m, flush=True)

def bh(pvals):
    p = np.asarray(pvals, float); out = np.full_like(p, np.nan)
    m = ~np.isnan(p)
    ps = p[m]
    if len(ps) == 0: return out
    order = np.argsort(ps); n = len(ps); adj = np.empty(n)
    for k in range(n - 1, -1, -1):
        adj[order[k]] = min(ps[order[k]] * n / (k + 1), adj[order[k + 1]] if k + 1 < n else 1.0)
    out[m] = adj; return out

def zc(s):
    s = np.asarray(s, float); sd = s.std(ddof=1)
    return (s - s.mean()) / sd if sd > 0 else s * 0

# ---------------- 基因集装载 ----------------
sig = pd.read_csv(SIG)
panels = pd.read_csv(PANELS)
GENESETS = {
    "IIAMD_up": set(sig.loc[sig["direction"] == "up", "humanSymbol"]),
    "IIAMD_down": set(sig.loc[sig["direction"] == "down", "humanSymbol"]),
}
for pname, sub in panels.groupby("panel"):
    GENESETS[pname] = set(sub["gene"].astype(str).str.upper())
anchor_h = pd.read_csv(os.path.join(INTER, "P2_anchor_human.csv"))
GENESETS["anchor_mechanism"] = set(anchor_h["humanSymbol"])
note("基因集: " + str({k: len(v) for k, v in GENESETS.items()}))

# ================= Phase A：GSE32707 =================
note("\n== Phase A: GSE32707 bulk（ARDS_d0 vs Control + Sepsis_d0 次对比）")
m7 = pd.read_csv(RAW + r"\GSE32707_data\GSE32707_gene_symbol_log2_matrix.csv.gz", index_col=0)
m7.index = [str(i).strip().upper() for i in m7.index]
m7.columns = [str(c) for c in m7.columns]
mat = m7.groupby(level=0).mean().T          # 样本 × 基因（log2）
meta7, _, _ = L.parse_series_matrix(RAW + r"\GSE32707_data\GSE32707_series_matrix.txt.gz")
gmap = {"untreated": "Control", "SIRS Day 0": "SIRS_d0", "Sepsis Day 0": "Sepsis_d0",
        "Sepsis Day 7": "Sepsis_d7", "se/ARDS Day 0": "ARDS_d0", "se/ARDS Day 7": "ARDS_d7"}
meta7df = pd.DataFrame({"gsm": meta7["!Sample_geo_accession"][0],
                        "group": [gmap.get(s, None) for s in meta7["!Sample_source_name_ch1"][0]]}).set_index("gsm")
mat = mat.loc[mat.index.intersection(meta7df.index)]
meta7df = meta7df.loc[mat.index]
note(f"  样本 {mat.shape[0]}, 基因 {mat.shape[1]}; 分组 {meta7df['group'].value_counts().to_dict()}")

# 髓系组成代理（标志物均值 z；操作化披露）
MYELOID = ["LST1", "S100A8", "S100A9", "CD14", "FCGR3A", "MS4A7"]
avail_my = [g for g in MYELOID if g in mat.columns]
myeloid = zc(mat[avail_my].mean(axis=1))
note(f"  髓系标志物可用 {len(avail_my)}/{len(MYELOID)}")

rowsA = []
for gs_name, genes in GENESETS.items():
    avail = [g for g in genes if g in mat.columns]
    if len(avail) < 5:
        note(f"  [跳过] {gs_name}: 可测基因 {len(avail)}<5")
        continue
    raw = mat[avail].mean(axis=1)
    score = zc(raw)
    for contrast, case_g, ctrl_g in [("ARDS_d0_vs_Control", "ARDS_d0", "Control"),
                                      ("Sepsis_d0_vs_Control", "Sepsis_d0", "Control")]:
        case = score[meta7df["group"] == case_g]; ctrl = score[meta7df["group"] == ctrl_g]
        if len(case) < 3 or len(ctrl) < 3:
            continue
        t, p = ttest_ind(case, ctrl, equal_var=False)
        g_eff, _ = L.hedges_g(np.asarray(ctrl), np.asarray(case))
        # 组成调整（次对比同做，报告调整后 p）
        try:
            import statsmodels.api as sm
            y = np.asarray(score, float)
            X = pd.DataFrame({"case": (meta7df["group"] == case_g).astype(float).values,
                              "myeloid": np.asarray(myeloid, float)})
            X = sm.add_constant(X)
            fit = sm.OLS(y, X).fit()
            p_adj = float(fit.pvalues["case"]); b_adj = float(fit.params["case"])
        except Exception as e:
            if not rowsA or all(pd.isna(r.get("OLS_p_adj")) for r in rowsA[:1]):
                note(f"  [OLS调试] {type(e).__name__}: {e}")
            p_adj, b_adj = np.nan, np.nan
        rowsA.append(dict(geneset=gs_name, n_genes_avail=len(avail), n_genes_total=len(genes),
                          coverage=len(avail) / len(genes), contrast=contrast,
                          n_case=len(case), n_ctrl=len(ctrl),
                          diff=float(case.mean() - ctrl.mean()), hedges_g=g_eff,
                          Welch_p=float(p), OLS_p_adj=p_adj, beta_adj=b_adj))

dfA = pd.DataFrame(rowsA)
prim = dfA[dfA["contrast"] == "ARDS_d0_vs_Control"].copy()
prim["BH_q"] = bh(prim["Welch_p"].values)
dfA = dfA.merge(prim[["geneset", "BH_q"]], on="geneset", how="left", suffixes=("", "_prim"))
dfA.to_csv(os.path.join(INTER, "P2_gate2_scores_GSE32707.csv"), index=False)
note("\n  -- GSE32707 ARDS_d0 vs Control 头对头（按 g 排序）--")
for _, r in prim.sort_values("hedges_g", ascending=False).iterrows():
    note(f"   {r['geneset']:<32} g={r['hedges_g']:+.2f} p={r['Welch_p']:.3g} q={r['BH_q']:.3g} "
         f"p_adj={r['OLS_p_adj']:.3g} cov={r['coverage']:.0%}")

# 头对头判定（主对比）
iiamd = prim[prim["geneset"] == "IIAMD_diff"] if "IIAMD_diff" in set(prim["geneset"]) else None
panels_only = prim[prim["geneset"].isin([p for p in GENESETS if p not in ("IIAMD_up", "IIAMD_down", "IIAMD_diff", "anchor_mechanism")])]

# ================= Phase A2：GSE185263 发现层（ENSG→symbol 全基因组映射）=================
note("\n== Phase A2: GSE185263 发现层（Sepsis_COVID / Sepsis vs Control）")
try:
    import mygene
    cnt = pd.read_csv(RAW + r"\GSE185263_Lung_ARDS\GSE185263_raw_counts.csv", index_col=0)
    cnt.index = [str(i).upper() for i in cnt.index]
    ensgs = list(cnt.index)
    cache = os.path.join(INTER, "P2_gse185263_en2sym.csv")
    if os.path.exists(cache):
        e2s = dict(pd.read_csv(cache).values)
        e2s = {str(k).upper(): str(v).upper() for k, v in e2s.items()}
    else:
        mg = mygene.MyGeneInfo()
        res = mg.querymany(ensgs, scopes="ensembl.gene", fields="symbol", species="human",
                           returnall=True, verbose=False)
        e2s = {str(r["query"]).upper(): str(r["symbol"]).upper() for r in res.get("out", []) if r.get("symbol")}
        pd.DataFrame({"ensg": list(e2s.keys()), "symbol": list(e2s.values())}).to_csv(cache, index=False)
    note(f"  ENSG→symbol 映射 {len(e2s)}/{len(ensgs)}")
    sub = cnt.loc[[e for e in e2s if e in cnt.index]].copy()
    sub.index = [e2s[e] for e in sub.index if e in cnt.index]
    sub = sub.groupby(level=0).mean()
    cpm = sub / sub.sum(axis=0) * 1e6
    logc = np.log2(cpm + 1).T
    groups = pd.read_csv(ROOT + r"\04_AUDIT_GOVERNANCE\GSE185263_groups.csv")
    groups.columns = ["sample_id", "group"]
    groups["sample_id"] = groups["sample_id"].astype(str)
    meta1 = pd.DataFrame(index=logc.index); meta1["sample_id"] = meta1.index
    meta1 = meta1.merge(groups, on="sample_id", how="left").set_index(meta1.index)
    avail_my2 = [g for g in MYELOID if g in logc.columns]
    myeloid2 = zc(logc[avail_my2].mean(axis=1))
    import statsmodels.api as sm
    rowsA2 = []
    for gs_name, genes in GENESETS.items():
        avail = [g for g in genes if g in logc.columns]
        if len(avail) < 5: continue
        score2 = zc(logc[avail].mean(axis=1))
        for contrast, case_g, ctrl_g in [("SepsisCOVID_vs_Control", "Sepsis_COVID", "Control"),
                                          ("Sepsis_vs_Control", "Sepsis", "Control")]:
            case = score2[(meta1["group"] == case_g).values]; ctrl = score2[(meta1["group"] == ctrl_g).values]
            if len(case) < 3 or len(ctrl) < 3: continue
            t, p = ttest_ind(case, ctrl, equal_var=False)
            g_eff, _ = L.hedges_g(np.asarray(ctrl), np.asarray(case))
            y = np.asarray(score2, float)
            X = pd.DataFrame({"case": (meta1["group"] == case_g).astype(float).values,
                              "myeloid": np.asarray(myeloid2, float)})
            X = sm.add_constant(X)
            fit = sm.OLS(y, X).fit()
            rowsA2.append(dict(geneset=gs_name, n_genes_avail=len(avail), coverage=len(avail)/len(genes),
                               contrast=contrast, n_case=len(case), n_ctrl=len(ctrl),
                               diff=float(case.mean()-ctrl.mean()), hedges_g=g_eff,
                               Welch_p=float(p), OLS_p_adj=float(fit.pvalues["case"])))
    dfA2 = pd.DataFrame(rowsA2)
    prim2 = dfA2[dfA2["contrast"] == "SepsisCOVID_vs_Control"].copy()
    prim2["BH_q"] = bh(prim2["Welch_p"].values)
    dfA2.to_csv(os.path.join(INTER, "P2_gate2_scores_GSE185263.csv"), index=False)
    note("  -- GSE185263 Sepsis_COVID vs Control（按 g 排序）--")
    for _, r in prim2.sort_values("hedges_g", ascending=False).iterrows():
        note(f"   {r['geneset']:<32} g={r['hedges_g']:+.2f} p={r['Welch_p']:.3g} q={r['BH_q']:.3g} p_adj={r['OLS_p_adj']:.3g} cov={r['coverage']:.0%}")
except Exception as e:
    import traceback; traceback.print_exc()
    dfA2 = pd.DataFrame(); prim2 = pd.DataFrame()
    note(f"  [WARN] GSE185263 层失败: {type(e).__name__}: {e}")

# ================= Phase B：scRNA 供者级定位 =================
note("\n== Phase B: scRNA 供者级 IIAMD 定位（severe vs Healthy，供者为单位）")
import scanpy as sc
import scipy.sparse as sp
adata = sc.read_h5ad(os.path.join(ROOT, "00_RAW_DATA", "_scrna_work", "combined_processed.h5ad"))
obs = adata.obs.copy()
varnames = list(adata.var_names)
up_idx = [varnames.index(g) for g in GENESETS["IIAMD_up"] if g in varnames]
dn_idx = [varnames.index(g) for g in GENESETS["IIAMD_down"] if g in varnames]
note(f"  IIAMD 在 scRNA 可测: up {len(up_idx)}/{len(GENESETS['IIAMD_up'])}, down {len(dn_idx)}/{len(GENESETS['IIAMD_down'])}")

CASE, CTRL = "COVID_severe", "Healthy"
MIN_CELLS, MIN_DONORS = 50, 3
rowsB = []
for ds in ("GSE145926", "GSE158055"):
    o = obs[obs["dataset"] == ds]
    conds = o[["sampleID", "condition"]].drop_duplicates().set_index("sampleID")["condition"]
    for ct in sorted(o["cell_type"].unique()):
        sub = o[o["cell_type"] == ct]
        don = sub.groupby("sampleID").size()
        don = don[don >= MIN_CELLS]
        donors = [d for d in don.index if conds.get(d) in (CASE, CTRL)]
        n_case = sum(conds.get(d) == CASE for d in donors); n_ctrl = sum(conds.get(d) == CTRL for d in donors)
        if n_case < MIN_DONORS or n_ctrl < MIN_DONORS:
            continue
        pos_all = adata.obs_names.get_indexer(sub.index)
        recs = []
        for d in donors:
            mask = (sub["sampleID"] == d).values
            pos = pos_all[mask]
            cs = adata.X[pos, :]
            if sp.issparse(cs):
                mu_up = float(np.asarray(cs[:, up_idx].mean(axis=0)).ravel().mean())
                mu_dn = float(np.asarray(cs[:, dn_idx].mean(axis=0)).ravel().mean())
            else:
                mu_up = float(cs[:, up_idx].mean()); mu_dn = float(cs[:, dn_idx].mean())
            recs.append((d, conds[d], mu_up, mu_dn))
        dfr = pd.DataFrame(recs, columns=["sampleID", "condition", "up_raw", "dn_raw"])
        dfr["up_z"] = zc(dfr["up_raw"]); dfr["dn_z"] = zc(dfr["dn_raw"])
        dfr["IIAMD_diff"] = dfr["up_z"] - dfr["dn_z"]
        for score_name in ("up_z", "dn_z", "IIAMD_diff"):
            a_ = dfr.loc[dfr["condition"] == CASE, score_name].values
            b_ = dfr.loc[dfr["condition"] == CTRL, score_name].values
            t, p = ttest_ind(a_, b_, equal_var=False)
            rowsB.append(dict(dataset=ds, cell_type=ct, score=score_name,
                              n_case=len(a_), n_ctrl=len(b_),
                              diff=float(a_.mean() - b_.mean()), p=float(p)))
dfB = pd.DataFrame(rowsB)
dfB["BH_q"] = bh(dfB["p"].values)
dfB.to_csv(os.path.join(INTER, "P2_gate2_scRNA_localization.csv"), index=False)
for _, r in dfB.iterrows():
    note(f"   {r['dataset']} {r['cell_type']:<10} {r['score']:<11} diff={r['diff']:+.3f} p={r['p']:.3g} q={r['BH_q']:.3g}")

# ================= 报告 =================
def df_md(d, cols=None):
    if d is None or len(d) == 0: return "(empty)"
    cols = cols or list(d.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, rr in d.iterrows():
        vals = [f"{rr[c]:.4g}" if isinstance(rr[c], (int, float, np.floating)) and pd.notna(rr[c]) else str(rr[c]) for c in cols]
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)

with open(REPORT, "w", encoding="utf-8") as f:
    f.write("# P2-2 Gate 2 头对头报告（H7）\n\n")
    f.write("- 日期：2026-08-20；冻结依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md（OSF https://osf.io/C7RYD/）\n")
    f.write("- 签名：IIAMD v1.0 人源一对一映射 7,434（up 3,615 / down 3,819）\n")
    f.write("- 组成调整：bulk 髓系标志物均值 z（LST1/S100A8/S100A9/CD14/FCGR3A/MS4A7）作为髓系比例的操作化代理（披露）\n\n")
    f.write("## Phase A：GSE32707（真 ARDS 全血，ARDS_d0 vs Control）\n\n")
    f.write(df_md(prim, ["geneset", "n_genes_avail", "coverage", "n_case", "n_ctrl", "diff", "hedges_g", "Welch_p", "BH_q", "OLS_p_adj"]))
    f.write("\n\n（Sepsis_d0_vs_Control 次对比见 _intermediate/P2_gate2_scores_GSE32707.csv）\n")
    f.write("\n## Phase B：scRNA 供者级 IIAMD 定位\n\n")
    f.write(df_md(dfB, ["dataset", "cell_type", "score", "n_case", "n_ctrl", "diff", "p", "BH_q"]))
    f.write("\n## Gate 2 初判\n\n")
    f.write("按冻结判据：强通过需 ≥2 个来源互斥真 ARDS 队列同向复制；本层可用真 ARDS bulk 队列仅 GSE32707 一个，"
            "因此最高只能给『部分通过』档（候选相关程序），除非后续队列核验扩充。\n")

with open(LOGF, "w", encoding="utf-8") as f:
    f.write("\n".join(log))
print("\nDONE P2-2 ->", REPORT)
