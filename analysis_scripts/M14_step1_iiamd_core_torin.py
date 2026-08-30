# -*- coding: utf-8 -*-
"""
M14_step1_iiamd_core_torin.py — M14 腿2（IIAMD core + Gate 2 重跑）+ 腿1（Torin 正式检验）
========================================================================================
前瞻注册：M10_M13_M14_pre_registration_20260824.md（osf.io/ETVMJ）§3.1/§3.2
- 腿2 预注册口径：core = FDR<0.05 且 |beta3|>1 且 |beta3| > |beta_lps|+|beta_cs|（因子 1x）；
  目标 200-500 基因；core/full 两版并存；用 core 重跑 Gate 2 三层（GSE32707 真 ARDS /
  GSE185263 发现层 / 单细胞供者级），同 full 版口径（臂均值->队列内 z，患者级线性检验
  + 髓系组成调整）。预注册条件语句：core 仍不迁移 -> 原 Gate 2 结论稳健；core 迁移 -> 修正。
- 腿1 预注册口径：签名基因在 Torin 臂的方向逆转比例 vs baseMean 十分位匹配随机集
  （B=2000, seed=0, 单侧 p）；H0 = 逆转比例与随机匹配基因集相同；负对照自校准 = 零分布
  中心应接近 0.5。
输入：GSE235046_Count_table（5 组 x 3 重复）、P1_GSE235046_interaction_DE.csv、
      P1_IIAMD_signature_v1.0.csv、P2_IIAMD_human_v1.0.csv、P2_competing_panels_v1.0.csv
输出：
  _intermediate/M14_GSE235046_universe_effects.csv（全宇宙 b3/b_lps/b_cs/b_torin）
  _intermediate/M14_IIAMD_core_v1.0.csv（core 小鼠+人源映射，冻结资产）
  _intermediate/M14_gate2_core_scores_GSE32707.csv / _GSE185263.csv / _scRNA_localization.csv
  _intermediate/M14_Torin_reversal_permutation.csv
  04_AUDIT_GOVERNANCE/M14_IIAMD_core_and_Torin_report.md
"""
import os, sys, warnings, time
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
RAW = ROOT + r"\00_RAW_DATA"
INTER = ROOT + r"\_intermediate"
GOV = ROOT + r"\04_AUDIT_GOVERNANCE"
CNT = os.path.join(RAW, "GSE235046", "GSE235046_Count_table.txt", "GSE235046_Count_table.txt")
DE_F = os.path.join(INTER, "P1_GSE235046_interaction_DE.csv")
SIGF = os.path.join(GOV, "P1_IIAMD_signature_v1.0.csv")
HUMANF = os.path.join(INTER, "P2_IIAMD_human_v1.0.csv")
PANELS = os.path.join(GOV, "P2_competing_panels_v1.0.csv")
SEED = 0
B = 2000

sys.path.insert(0, ROOT)
import mdi_lib as L

t0 = time.time()
log = []
def note(m=""):
    log.append(m); print(m, flush=True)

def zc(s):
    s = np.asarray(s, float); sd = s.std(ddof=1)
    return (s - s.mean()) / sd if sd > 0 else s * 0

# ======================================================================
# 1. 读入 count 表与组映射（与 P1 完全一致）
# ======================================================================
note("== 1. 读入 count 表与组映射")
raw = pd.read_csv(CNT, sep="\t")
meta_cols = ["geneID", "geneSymbol", "bioType", "annotationLevel"]
samples = [c for c in raw.columns if c not in meta_cols]
YW_GROUP = {}
for yw in ["YW001", "YW002", "YW003"]: YW_GROUP[yw] = "Media"
for yw in ["YW004", "YW005", "YW006"]: YW_GROUP[yw] = "LPS"
for yw in ["YW007", "YW008", "YW009"]: YW_GROUP[yw] = "CS"
for yw in ["YW010", "YW011", "YW012"]: YW_GROUP[yw] = "LPS_CS"
for yw in ["YW016", "YW017", "YW018"]: YW_GROUP[yw] = "LPS_CS_Torin"
def col_group(c):
    return YW_GROUP.get(c.split("_")[-1])
grp_map = {c: col_group(c) for c in samples}
missing = [c for c, g in grp_map.items() if g is None]
assert not missing, f"未映射样本列: {missing}"
raw["geneSymbol"] = raw["geneSymbol"].astype(str).str.strip()
ok = raw[raw["geneSymbol"].notna() & (raw["geneSymbol"] != "") & (raw["geneSymbol"] != "nan")]
mat = ok.groupby("geneSymbol")[samples].sum().astype(int)
note(f"   聚合后 {mat.shape[0]} 基因 x {mat.shape[1]} 样本")

# ======================================================================
# 2. DESeq2 确证设计：交互 b3 + 主效应 b_lps / b_cs
# ======================================================================
note("\n== 2. 确证 2x2 交互模型（12 样本）+ 主效应")
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

conf_cols = [c for c in samples if grp_map[c] in ("Media", "LPS", "CS", "LPS_CS")]
cmat = mat[conf_cols]
keep = (cmat >= 10).sum(axis=1) >= 3
cmatf = cmat[keep]
note(f"   预过滤后 {cmatf.shape[0]} 基因（>=10 counts 于 >=3 样本）")

counts_df = cmatf.T
meta = pd.DataFrame({"group": [grp_map[c] for c in cmatf.columns]}, index=cmatf.columns)
dds = DeseqDataSet(counts=counts_df, metadata=meta, design="~group", quiet=True, n_cpus=1)
dds.deseq2()
cols = list(dds.obsm["design_matrix"].columns)
tlevels = [c.split("T.")[1].rstrip("]") for c in cols if c.startswith("group[T.")]
ref = ({"Media", "LPS", "CS", "LPS_CS"} - set(tlevels)).pop()
note(f"   参考水平={ref}")
def mu_vec(level):
    v = np.zeros(len(cols)); v[cols.index("Intercept")] = 1.0
    if level != ref and f"group[T.{level}]" in cols:
        v[cols.index(f"group[T.{level}]")] = 1.0
    return v
cvec = mu_vec("LPS_CS") - mu_vec("LPS") - mu_vec("CS") + mu_vec("Media")
st3 = DeseqStats(dds, contrast=cvec, quiet=True, n_cpus=1)
st3.summary()
res3 = st3.results_df.reset_index().rename(columns={"index": "geneSymbol"})
res3 = res3.rename(columns={"log2FoldChange": "b3", "pvalue": "b3_p", "padj": "b3_padj", "lfcSE": "b3_se"})
res3 = res3[["geneSymbol", "baseMean", "b3", "b3_se", "b3_p", "b3_padj"]].set_index("geneSymbol")

# 主效应
def main_effect(level):
    st = DeseqStats(dds, contrast=["group", level, ref], quiet=True, n_cpus=1)
    st.summary()
    r = st.results_df.reset_index().rename(columns={"index": "geneSymbol"})
    return r.set_index("geneSymbol")[["log2FoldChange", "pvalue", "padj"]]
st_lps = main_effect("LPS")
st_cs = main_effect("CS")
univ = res3.join(st_lps.rename(columns={"log2FoldChange": "b_lps", "pvalue": "b_lps_p", "padj": "b_lps_padj"}))
univ = univ.join(st_cs.rename(columns={"log2FoldChange": "b_cs", "pvalue": "b_cs_p", "padj": "b_cs_padj"}))

# 与 P1 冻结 b3 交叉验证（数值一致性自检）
old = pd.read_csv(DE_F).set_index("geneSymbol")
chk = pd.DataFrame({"new_b3": univ["b3"], "old_b3": old["log2FoldChange"],
                    "new_padj": univ["b3_padj"], "old_padj": old["padj"]}).dropna()
note(f"   b3 与 P1 冻结表最大差异: {(chk['new_b3']-chk['old_b3']).abs().max():.3g}; "
     f"padj 最大差异: {(chk['new_padj']-chk['old_padj']).abs().max():.3g}")

# ======================================================================
# 3. Torin 对比（LPS_CS_Torin vs LPS_CS，6 样本）
# ======================================================================
note("\n== 3. Torin 对比（LPS_CS_Torin vs LPS_CS）")
exp_cols = [c for c in samples if grp_map[c] in ("LPS_CS", "LPS_CS_Torin")]
emat = mat.loc[cmatf.index, exp_cols]
emeta = pd.DataFrame({"group": [grp_map[c] for c in exp_cols]}, index=exp_cols)
edds = DeseqDataSet(counts=emat.T, metadata=emeta, design="~group", quiet=True, n_cpus=1)
edds.deseq2()
stt = DeseqStats(edds, contrast=["group", "LPS_CS_Torin", "LPS_CS"], quiet=True, n_cpus=1)
stt.summary()
rt = stt.results_df.reset_index().rename(columns={"index": "geneSymbol"})
rt = rt.set_index("geneSymbol")[["log2FoldChange", "pvalue", "padj"]].rename(
    columns={"log2FoldChange": "b_torin", "pvalue": "b_torin_p", "padj": "b_torin_padj"})
univ = univ.join(rt)
univ = univ.reset_index()
univ.to_csv(os.path.join(INTER, "M14_GSE235046_universe_effects.csv"), index=False)
note(f"   宇宙效应表 {univ.shape[0]} 基因 -> M14_GSE235046_universe_effects.csv")

# ======================================================================
# 4. IIAMD core 构建（预注册口径）
# ======================================================================
note("\n== 4. IIAMD core 构建")
u = univ.set_index("geneSymbol")
sig_up_mask = (u["b3_padj"] < 0.05) & (u["b3"] > 0)
sig_dn_mask = (u["b3_padj"] < 0.05) & (u["b3"] < 0)
main_sum = u["b_lps"].abs() + u["b_cs"].abs()
core_up_mask = sig_up_mask & (u["b3"] > 1) & (u["b3"] > main_sum)
core_dn_mask = sig_dn_mask & (u["b3"].abs() > 1) & (u["b3"].abs() > main_sum)
n_full_up, n_full_dn = int(sig_up_mask.sum()), int(sig_dn_mask.sum())
n_core_up, n_core_dn = int(core_up_mask.sum()), int(core_dn_mask.sum())
note(f"   full: up={n_full_up} dn={n_full_dn}; core: up={n_core_up} dn={n_core_dn} "
     f"(合计 {n_core_up+n_core_dn}, 目标 200-500)")
# 敏感性（判定门不依赖，仅登记）
sens_rows = []
for nm, m in [("core_up", core_up_mask), ("core_dn", core_dn_mask)]:
    pass
sens = dict(
    fdr_only_up=n_full_up, fdr_only_dn=n_full_dn,
    core_1x_up=n_core_up, core_1x_dn=n_core_dn,
    core_15x_up=int((sig_up_mask & (u["b3"] > 1) & (u["b3"] > 1.5 * main_sum)).sum()),
    core_15x_dn=int((sig_dn_mask & (u["b3"].abs() > 1) & (u["b3"].abs() > 1.5 * main_sum)).sum()),
    core_2x_up=int((sig_up_mask & (u["b3"] > 1) & (u["b3"] > 2 * main_sum)).sum()),
    core_2x_dn=int((sig_dn_mask & (u["b3"].abs() > 1) & (u["b3"].abs() > 2 * main_sum)).sum()),
    core_3x_up=int((sig_up_mask & (u["b3"] > 1) & (u["b3"] > 3 * main_sum)).sum()),
    core_3x_dn=int((sig_dn_mask & (u["b3"].abs() > 1) & (u["b3"].abs() > 3 * main_sum)).sum()),
)
note(f"   敏感性口径（1.5x/2x）: {sens}")

core_df = u[core_up_mask | core_dn_mask].reset_index()
core_df["direction"] = np.where(core_df["b3"] > 0, "up", "down")
# 人源一对一映射
hum = pd.read_csv(HUMANF)
hum["mouseSymbol"] = hum["mouseSymbol"].astype(str).str.upper()
core_df["mouseSymbol"] = core_df["geneSymbol"].astype(str).str.upper()
core_h = core_df.merge(hum[["mouseSymbol", "humanSymbol"]], on="mouseSymbol", how="left")
n_mapped = int(core_h["humanSymbol"].notna().sum())
note(f"   人源一对一映射: {n_mapped}/{len(core_h)} ({n_mapped/len(core_h):.0%})")
core_h.insert(0, "signature_version", "IIAMD_core_v1.0")
core_h.insert(1, "frozen_date", "2026-08-26")
core_h.insert(2, "source", "GSE235046 LPSxCS interaction (FDR<0.05 & |b3|>1 & |b3|>|b_lps|+|b_cs|)")
core_h.to_csv(os.path.join(INTER, "M14_IIAMD_core_v1.0.csv"), index=False)
note(f"   core 冻结 -> M14_IIAMD_core_v1.0.csv（up {n_core_up} / down {n_core_dn}）")

# ======================================================================
# 5. Gate 2 三层重跑（core 版；同 full 版口径）
# ======================================================================
note("\n== 5. Gate 2 重跑（core 版，三层）")
hum_map = pd.read_csv(HUMANF)
GENESETS = {
    "IIAMD_core_up": set(core_h.loc[core_h["direction"] == "up", "humanSymbol"].dropna()),
    "IIAMD_core_dn": set(core_h.loc[core_h["direction"] == "down", "humanSymbol"].dropna()),
    # full 版沿用 P2 的人源一对一映射（与 P2 Gate 2 口径一致）
    "IIAMD_full_up": set(hum_map.loc[hum_map["direction"] == "up", "humanSymbol"].dropna()),
    "IIAMD_full_dn": set(hum_map.loc[hum_map["direction"] == "down", "humanSymbol"].dropna()),
}
panels = pd.read_csv(PANELS)
for pname, sub in panels.groupby("panel"):
    GENESETS[pname] = set(sub["gene"].astype(str).str.upper())
anchor_h = pd.read_csv(os.path.join(INTER, "P2_anchor_human.csv"))
GENESETS["anchor_mechanism"] = set(anchor_h["humanSymbol"])
GENESETS["IIAMD_core_up"] = {g.upper() for g in GENESETS["IIAMD_core_up"]}
GENESETS["IIAMD_core_dn"] = {g.upper() for g in GENESETS["IIAMD_core_dn"]}
GENESETS["IIAMD_full_up"] = {g.upper() for g in GENESETS["IIAMD_full_up"]}
GENESETS["IIAMD_full_dn"] = {g.upper() for g in GENESETS["IIAMD_full_dn"]}
note("   基因集: " + str({k: len(v) for k, v in GENESETS.items()}))

def bh(pvals):
    p = np.asarray(pvals, float); out = np.full_like(p, np.nan)
    m = ~np.isnan(p)
    ps = p[m]
    if len(ps) == 0: return out
    order = np.argsort(ps); n = len(ps); adj = np.empty(n)
    for k in range(n - 1, -1, -1):
        adj[order[k]] = min(ps[order[k]] * n / (k + 1), adj[order[k + 1]] if k + 1 < n else 1.0)
    out[m] = adj; return out

MYELOID = ["LST1", "S100A8", "S100A9", "CD14", "FCGR3A", "MS4A7"]

def layer_bulk(matlog, meta_df, group_col, contrasts, outpath, label):
    """matlog: 样本x基因(log2, 大写符号)；meta_df: index=样本, 含 group_col。"""
    import statsmodels.api as sm
    avail_my = [g for g in MYELOID if g in matlog.columns]
    myeloid = zc(matlog[avail_my].mean(axis=1)) if avail_my else np.zeros(len(matlog))
    note(f"   [{label}] 髓系标志物可用 {len(avail_my)}/{len(MYELOID)}")
    rows = []
    for gs_name, genes in GENESETS.items():
        avail = [g for g in genes if g in matlog.columns]
        if len(avail) < 5:
            continue
        rawv = matlog[avail].mean(axis=1)
        score = zc(rawv)
        for contrast, case_g, ctrl_g in contrasts:
            case = score[(meta_df[group_col] == case_g).values]
            ctrl = score[(meta_df[group_col] == ctrl_g).values]
            if len(case) < 3 or len(ctrl) < 3:
                continue
            t, p = __import__("scipy.stats", fromlist=["ttest_ind"]).ttest_ind(case, ctrl, equal_var=False)
            g_eff, _ = L.hedges_g(np.asarray(ctrl), np.asarray(case))
            y = np.asarray(score, float)
            X = pd.DataFrame({"case": (meta_df[group_col] == case_g).astype(float).values,
                              "myeloid": np.asarray(myeloid, float)})
            X = sm.add_constant(X)
            fit = sm.OLS(y, X).fit()
            rows.append(dict(geneset=gs_name, n_genes_avail=len(avail), n_genes_total=len(genes),
                             coverage=len(avail) / len(genes), contrast=contrast,
                             n_case=len(case), n_ctrl=len(ctrl),
                             diff=float(case.mean() - ctrl.mean()), hedges_g=g_eff,
                             Welch_p=float(p), OLS_p_adj=float(fit.pvalues["case"]),
                             beta_adj=float(fit.params["case"])))
    df = pd.DataFrame(rows)
    prim = df[df["contrast"] == contrasts[0][0]].copy()
    if len(prim):
        prim["BH_q"] = bh(prim["Welch_p"].values)
        df = df.merge(prim[["geneset", "BH_q"]], on="geneset", how="left")
    df.to_csv(outpath, index=False)
    note(f"   [{label}] 保存 {outpath}（{len(df)} 行）")
    if len(prim):
        for _, r in prim.sort_values("hedges_g", ascending=False).iterrows():
            note(f"     {r['geneset']:<24} g={r['hedges_g']:+.2f} p={r['Welch_p']:.3g} q={r['BH_q']:.3g} p_adj={r['OLS_p_adj']:.3g}")
    return df

# --- 层1 GSE32707 ---
note("\n   -- 层1 GSE32707（ARDS_d0 vs Control）--")
m7 = pd.read_csv(RAW + r"\GSE32707_data\GSE32707_gene_symbol_log2_matrix.csv.gz", index_col=0)
m7.index = [str(i).strip().upper() for i in m7.index]
m7.columns = [str(c) for c in m7.columns]
mat7 = m7.groupby(level=0).mean().T
meta7, _, _ = L.parse_series_matrix(RAW + r"\GSE32707_data\GSE32707_series_matrix.txt.gz")
gmap = {"untreated": "Control", "SIRS Day 0": "SIRS_d0", "Sepsis Day 0": "Sepsis_d0",
        "Sepsis Day 7": "Sepsis_d7", "se/ARDS Day 0": "ARDS_d0", "se/ARDS Day 7": "ARDS_d7"}
meta7df = pd.DataFrame({"gsm": meta7["!Sample_geo_accession"][0],
                        "group": [gmap.get(s, None) for s in meta7["!Sample_source_name_ch1"][0]]}).set_index("gsm")
mat7 = mat7.loc[mat7.index.intersection(meta7df.index)]
meta7df = meta7df.loc[mat7.index]
layer_bulk(mat7, meta7df, "group",
           [("ARDS_d0_vs_Control", "ARDS_d0", "Control"),
            ("Sepsis_d0_vs_Control", "Sepsis_d0", "Control")],
           os.path.join(INTER, "M14_gate2_core_scores_GSE32707.csv"), "GSE32707")

# --- 层2 GSE185263 ---
note("\n   -- 层2 GSE185263（Sepsis_COVID vs Control）--")
cnt = pd.read_csv(RAW + r"\GSE185263_Lung_ARDS\GSE185263_raw_counts.csv", index_col=0)
cnt.index = [str(i).upper() for i in cnt.index]
ensgs = list(cnt.index)
cache = os.path.join(INTER, "P2_gse185263_en2sym.csv")
if os.path.exists(cache):
    e2s = dict(pd.read_csv(cache).values)
    e2s = {str(k).upper(): str(v).upper() for k, v in e2s.items()}
else:
    import mygene
    mg = mygene.MyGeneInfo()
    res = mg.querymany(ensgs, scopes="ensembl.gene", fields="symbol", species="human",
                       returnall=True, verbose=False)
    e2s = {str(r["query"]).upper(): str(r["symbol"]).upper() for r in res.get("out", []) if r.get("symbol")}
    pd.DataFrame({"ensg": list(e2s.keys()), "symbol": list(e2s.values())}).to_csv(cache, index=False)
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
layer_bulk(logc, meta1, "group",
           [("SepsisCOVID_vs_Control", "Sepsis_COVID", "Control"),
            ("Sepsis_vs_Control", "Sepsis", "Control")],
           os.path.join(INTER, "M14_gate2_core_scores_GSE185263.csv"), "GSE185263")

# --- 层3 scRNA 供者级 ---
note("\n   -- 层3 scRNA 供者级 --")
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import ttest_ind
adata = sc.read_h5ad(os.path.join(RAW, "_scrna_work", "combined_processed.h5ad"))
obs = adata.obs.copy()
varnames = list(adata.var_names)
varnames_u = [str(v).upper() for v in varnames]
def idx_of(genes):
    return [varnames.index(varnames[i]) for i, v in enumerate(varnames_u) if v in genes]
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
        n_case = sum(conds.get(d) == CASE for d in donors)
        n_ctrl = sum(conds.get(d) == CTRL for d in donors)
        if n_case < MIN_DONORS or n_ctrl < MIN_DONORS:
            continue
        pos_all = adata.obs_names.get_indexer(sub.index)
        idx_cache = {}
        for gs_name in ("IIAMD_core_up", "IIAMD_core_dn", "IIAMD_full_up", "IIAMD_full_dn"):
            idx_cache[gs_name] = idx_of(GENESETS[gs_name])
        recs = []
        for d in donors:
            mask = (sub["sampleID"] == d).values
            pos = pos_all[mask]
            cs = adata.X[pos, :]
            if sp.issparse(cs):
                row = {}
                for gs_name, idx in idx_cache.items():
                    row[gs_name] = float(np.asarray(cs[:, idx].mean(axis=0)).ravel().mean()) if idx else np.nan
            else:
                row = {gs_name: float(cs[:, idx].mean()) if idx else np.nan for gs_name, idx in idx_cache.items()}
            recs.append((d, conds[d], row))
        dfr = pd.DataFrame([(r[0], r[1], *[r[2].get(k, np.nan) for k in idx_cache]) for r in recs],
                           columns=["sampleID", "condition"] + list(idx_cache))
        for gs_name in idx_cache:
            dfr[gs_name + "_z"] = zc(dfr[gs_name])
        dfr["core_diff"] = dfr["IIAMD_core_up_z"] - dfr["IIAMD_core_dn_z"]
        dfr["full_diff"] = dfr["IIAMD_full_up_z"] - dfr["IIAMD_full_dn_z"]
        for score_name in ("IIAMD_core_up_z", "IIAMD_core_dn_z", "core_diff",
                           "IIAMD_full_up_z", "IIAMD_full_dn_z", "full_diff"):
            a_ = dfr.loc[dfr["condition"] == CASE, score_name].values
            b_ = dfr.loc[dfr["condition"] == CTRL, score_name].values
            if len(a_) < 3 or len(b_) < 3:
                continue
            t, p = ttest_ind(a_, b_, equal_var=False)
            rowsB.append(dict(dataset=ds, cell_type=ct, score=score_name,
                              n_case=len(a_), n_ctrl=len(b_),
                              diff=float(a_.mean() - b_.mean()), p=float(p)))
dfB = pd.DataFrame(rowsB)
dfB["BH_q"] = bh(dfB["p"].values)
dfB.to_csv(os.path.join(INTER, "M14_gate2_core_scRNA_localization.csv"), index=False)
note(f"   层3 检验数 {len(dfB)}；core 相关行：")
for _, r in dfB[dfB["score"].str.startswith("IIAMD_core") | dfB["score"].eq("core_diff")].iterrows():
    note(f"     {r['dataset']} {r['cell_type']:<10} {r['score']:<16} diff={r['diff']:+.3f} p={r['p']:.3g} q={r['BH_q']:.3g}")

# ======================================================================
# 6. Torin 正式检验（腿1）
# ======================================================================
note("\n== 6. Torin 正式检验（baseMean 十分位匹配置换, B=2000, seed=0）")
sig = pd.read_csv(SIGF)
up_genes = [str(g).strip().upper() for g in sig[sig["direction"] == "up"]["geneSymbol"]]
dn_genes = [str(g).strip().upper() for g in sig[sig["direction"] == "down"]["geneSymbol"]]
univ["_gsym"] = univ["geneSymbol"].astype(str).str.upper()
univ = univ.set_index("_gsym")
univ_set = set(univ.index)
up_in = [g for g in up_genes if g in univ_set]
dn_in = [g for g in dn_genes if g in univ_set]
core_up_in = [g for g in core_h.loc[core_h["direction"] == "up", "mouseSymbol"] if g in univ_set]
core_dn_in = [g for g in core_h.loc[core_h["direction"] == "down", "mouseSymbol"] if g in univ_set]
note(f"   命中宇宙: full up={len(up_in)}/{len(up_genes)} dn={len(dn_in)}/{len(dn_genes)}; "
     f"core up={len(core_up_in)}/{n_core_up} dn={len(core_dn_in)}/{n_core_dn}")

u2 = univ.copy()
u2["bm_decile"] = pd.qcut(u2["baseMean"].rank(method="first"), 10, labels=False)
sig_all = set(up_in) | set(dn_in)
null_pool = u2[~u2.index.isin(sig_all)]
null_pool_dir = {"up": null_pool[null_pool["b3"] > 0],
                 "dn": null_pool[null_pool["b3"] < 0]}

def draw_matched(n_target, rng, pool):
    """P1 同款 baseMean 十分位分层按比例抽样（不允许超过池容量）。"""
    if n_target > len(pool):
        raise ValueError(f"n_target={n_target} > pool={len(pool)}")
    idx = []
    for d in range(10):
        layer = pool[pool["bm_decile"] == d]
        k = int(round(n_target * len(layer) / len(pool)))
        if len(layer) and k > 0:
            idx += rng.choice(layer.index.tolist(), size=min(k, len(layer)), replace=False).tolist()
    if len(idx) < n_target:
        remaining = pool.index.difference(idx).tolist()
        idx += rng.choice(remaining, size=n_target - len(idx), replace=False).tolist()
    return idx[:n_target]

def reversal_ratio(genes, direction):
    sub = u2.loc[genes]
    if direction == "up":
        return float((sub["b_torin"] < 0).mean())
    return float((sub["b_torin"] > 0).mean())

def reversal_ratio_sigmismatch(genes):
    """方向不限集合的逆转定义：sign(b_torin) != sign(b3)。"""
    sub = u2.loc[genes]
    return float(((sub["b_torin"] < 0) & (sub["b3"] > 0) | (sub["b_torin"] > 0) & (sub["b3"] < 0)).mean())

rng = np.random.default_rng(SEED)
rows = []
for direction, genes, label in (("up", up_in, "IIAMD_full_up"), ("dn", dn_in, "IIAMD_full_dn"),
                                ("up", core_up_in, "IIAMD_core_up"), ("dn", core_dn_in, "IIAMD_core_dn")):
    obs = reversal_ratio(genes, direction)
    # 主检验（注册口径，P1 同款）：全池 baseMean 十分位匹配，统计量 = sign(b_torin)!=sign(b3)
    pool = null_pool
    obs_mm = reversal_ratio_sigmismatch(genes)
    nulls = np.empty(B)
    for it in range(B):
        rs = draw_matched(len(genes), rng, pool)
        nulls[it] = reversal_ratio_sigmismatch(rs)
    p = float((1 + (nulls >= obs_mm).sum()) / (B + 1))
    lo, hi = np.percentile(nulls, 2.5), np.percentile(nulls, 97.5)
    rows.append(dict(gene_set=label, direction=direction, n_genes=len(genes),
                     observed_reversal=obs, observed_reversal_sigmismatch=obs_mm,
                     null_mean=float(nulls.mean()), null_sd=float(nulls.std()),
                     null_2p5=lo, null_97p5=hi, empirical_p_one_tail=p, n_perm=B,
                     calibration_center=float(abs(nulls.mean() - 0.5)),
                     null_type="baseMean-matched, full pool (P1-style)"))
    note(f"   {label:<14} n={len(genes):<5} obs(方向定义)={obs:.4f} obs(sign-mismatch)={obs_mm:.4f} "
         f"null={nulls.mean():.4f}±{nulls.std():.4f} [{lo:.4f},{hi:.4f}] 单侧p={p:.4g}")
    # 敏感性：方向受限池（b3 同号），仅当池容量足够
    pool_dir = null_pool_dir[direction]
    if len(genes) <= len(pool_dir):
        nulls_d = np.empty(B)
        for it in range(B):
            rs = draw_matched(len(genes), rng, pool_dir)
            nulls_d[it] = reversal_ratio(rs, direction)
        p_d = float((1 + (nulls_d >= obs).sum()) / (B + 1))
        lo_d, hi_d = np.percentile(nulls_d, 2.5), np.percentile(nulls_d, 97.5)
        rows.append(dict(gene_set=label, direction=direction, n_genes=len(genes),
                         observed_reversal=obs, observed_reversal_sigmismatch=np.nan,
                         null_mean=float(nulls_d.mean()), null_sd=float(nulls_d.std()),
                         null_2p5=lo_d, null_97p5=hi_d, empirical_p_one_tail=p_d, n_perm=B,
                         calibration_center=float(abs(nulls_d.mean() - 0.5)),
                         null_type="direction-restricted pool (b3 same sign)"))
        note(f"     [敏感性] null={nulls_d.mean():.4f}±{nulls_d.std():.4f} [{lo_d:.4f},{hi_d:.4f}] 单侧p={p_d:.4g}")
    else:
        note(f"     [敏感性跳过] 方向受限池容量 {len(pool_dir)} < n={len(genes)}")
out = pd.DataFrame(rows)
out.to_csv(os.path.join(INTER, "M14_Torin_reversal_permutation.csv"), index=False)

# ======================================================================
# 7. 报告
# ======================================================================
report = []
report.append("# M14 腿1/腿2 报告（IIAMD core + Torin 正式检验）\n")
report.append(f"- 日期：2026-08-26；依据：M10_M13_M14_pre_registration_20260824.md（osf.io/ETVMJ）§3.1/§3.2")
report.append(f"- 确证模型：12 样本 DESeq2 ~group，交互对比向量；参考水平={ref}")
report.append(f"- Torin 对比：LPS_CS_Torin vs LPS_CS（6 样本，DESeq2 Wald）")
report.append(f"- 腿2 预注册口径：FDR<0.05 且 |b3|>1 且 |b3| > |b_lps|+|b_cs|（因子 1x）；目标 200–500 基因\n")
report.append(f"## 1. core 构建\n")
report.append(f"- full 签名：up={n_full_up} / down={n_full_dn}（合计 {n_full_up+n_full_dn}）")
report.append(f"- **core：up={n_core_up} / down={n_core_dn}（合计 {n_core_up+n_core_dn}）**")
report.append(f"- 人源一对一映射：{n_mapped}/{len(core_h)}（{n_mapped/len(core_h):.0%}）")
report.append(f"- 敏感性口径登记：{sens}\n")
report.append("## 2. Gate 2 重跑（core 版，三层，同 full 版口径）\n")
report.append("- 层1 GSE32707（ARDS_d0 vs Control）、层2 GSE185263（Sepsis_COVID vs Control）、层3 scRNA 供者级")
report.append("- 全部基因集行见 _intermediate/M14_gate2_core_scores_*.csv / M14_gate2_core_scRNA_localization.csv")
report.append("- 预注册条件语句：core 仍不迁移 -> 原 Gate 2 结论稳健；core 迁移 -> 原结论修正。判定见下。\n")
report.append("## 3. Torin 正式检验（腿1）\n")
report.append("- 统计量：签名基因在 Torin 臂的方向逆转比例；H0 = 与 baseMean 十分位匹配随机集相同")
report.append("- 主检验（注册口径，P1 同款）：全池 baseMean 匹配，统计量 = sign(b_torin)!=sign(b3)；")
report.append("  敏感性：方向受限池（b3 同号）内统计量 = 方向定义逆转比例（池容量足够时）")
report.append("- 置换：B=2000, seed=0；单侧 p = (1+#null>=obs)/(B+1)；负对照自校准 = 零分布中心应≈0.5\n")
report.append("| 基因集 | 方向 | n | 观测逆转率 | 零分布 mean±sd [2.5%,97.5%] | 单侧 p | 中心偏移 | 零模型 |")
report.append("|---|---|---|---|---|---|---|---|---|")
for _, r in out.iterrows():
    report.append(f"| {r['gene_set']} | {r['direction']} | {r['n_genes']} | {r['observed_reversal']:.4f} | "
                  f"{r['null_mean']:.4f}±{r['null_sd']:.4f} [{r['null_2p5']:.4f},{r['null_97p5']:.4f}] | "
                  f"{r['empirical_p_one_tail']:.4g} | {r['calibration_center']:.4f} | {r['null_type']} |")
report.append("")
report.append("## 4. 判定\n")
# Gate 2 判定（core）
g7 = pd.read_csv(os.path.join(INTER, "M14_gate2_core_scores_GSE32707.csv"))
prim7 = g7[g7["contrast"] == "ARDS_d0_vs_Control"]
core_up_row = prim7[prim7["geneset"] == "IIAMD_core_up"]
core_dn_row = prim7[prim7["geneset"] == "IIAMD_core_dn"]
g185 = pd.read_csv(os.path.join(INTER, "M14_gate2_core_scores_GSE185263.csv"))
prim185 = g185[g185["contrast"] == "SepsisCOVID_vs_Control"]
core_up_185 = prim185[prim185["geneset"] == "IIAMD_core_up"]
core_dn_185 = prim185[prim185["geneset"] == "IIAMD_core_dn"]
sc_core = dfB[dfB["score"].isin(["IIAMD_core_up_z", "IIAMD_core_dn_z", "core_diff"])]
report.append(f"- 层1 GSE32707：core_up " +
              (f"g={float(core_up_row['hedges_g'].iloc[0]):+.2f}, Welch q={float(core_up_row['BH_q'].iloc[0]):.3g}, "
               f"OLS p={float(core_up_row['OLS_p_adj'].iloc[0]):.3g}" if len(core_up_row) else "（<5 基因可测，未计）") +
              f"；core_dn " +
              (f"g={float(core_dn_row['hedges_g'].iloc[0]):+.2f}, Welch q={float(core_dn_row['BH_q'].iloc[0]):.3g}, "
               f"OLS p={float(core_dn_row['OLS_p_adj'].iloc[0]):.3g}" if len(core_dn_row) else "（<5 基因可测，未计）"))
report.append(f"- 层2 GSE185263：core_up " +
              (f"g={float(core_up_185['hedges_g'].iloc[0]):+.2f}, Welch q={float(core_up_185['BH_q'].iloc[0]):.3g}, "
               f"OLS p={float(core_up_185['OLS_p_adj'].iloc[0]):.3g}" if len(core_up_185) else "（<5 基因可测，未计）") +
              f"；core_dn " +
              (f"g={float(core_dn_185['hedges_g'].iloc[0]):+.2f}, Welch q={float(core_dn_185['BH_q'].iloc[0]):.3g}, "
               f"OLS p={float(core_dn_185['OLS_p_adj'].iloc[0]):.3g}" if len(core_dn_185) else "（<5 基因可测，未计）"))
report.append(f"- 层3 scRNA 供者级：core 相关检验 {len(sc_core)} 项，" +
              (f"最小 q={float(sc_core['BH_q'].min()):.3g}" if len(sc_core) else "无合格 celltype×score"))
report.append("")
report.append("（判定结论见正文 Results；本报告只记录量化事实。）")
with open(os.path.join(GOV, "M14_IIAMD_core_and_Torin_report.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(report))
with open(os.path.join(ROOT, "03_LOGS", "M14_core_torin_log.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(log))
note(f"\nDONE in {time.time()-t0:.1f}s -> {GOV}\\M14_IIAMD_core_and_Torin_report.md")
