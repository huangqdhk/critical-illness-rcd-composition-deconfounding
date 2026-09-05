# -*- coding: utf-8 -*-
"""
M16_step4_bridge.py — M16B：人体桥接层（mech_v1.0 × UCS/EIS/MDI 共变）
=====================================================================
预注册：M16_pre_registration_20260831.md（判据 D，描述性为主；空间层探索性、家族外）
数据层：
  (a) GSE185263 全血 bulk（本地原始 counts）→ log2(CPM+1) → 队列内基因 z → mech 分数；
      与 Table_S56 逐样本 UCS/EIS/MDI 按样本合并；Spearman ρ + 偏 Spearman
      （组成协变量=Monaco|NNLS 中性粒+单核比例，Table_S64 同口径）。
  (b) GSE212865 外周血芯片（本地 series matrix + GPL23159_family.soft 重建探针映射）同口径。
  (c) GSE158055 scRNA（本地 Mono h5ad）：供体×细胞类型 mech 均值，与 Table_S67 供体级
      UCS/EIS/MDI 合并，单核/巨噬细胞内跨供体 Spearman（≥50 细胞/供体×类型沿用 S67 口径）。
  (d) Visium（GSE271370 归档评分 h5ad 复用，来源注册）：每切片 per-spot mech 与
      up_score_z/ex_score_z/mdi_z 的 Spearman + mech×arm 双变量 Moran's I
      （kNN k=8 行标准化，999 置换，种子 20260831，Stouffer 合并）——探索性。
BH 家族 = (a)(b)(c) 主相关（mech×UCS/EIS/MDI）；(d) 家族外。
输出：Table_S91_M16_Bridge_Correlations.csv、_intermediate/M16_visium_mech_spatial.csv、
      03_LOGS/M16_step4_log.txt
"""
import gzip
import os
import re
import sys
import time
import json
import urllib.request
import urllib.parse

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.abspath(__file__))
RAW = r"D:\!!!Research\!!!课题组\ARDS核心基因筛选策略方案\00_RAW_DATA"
DSET = r"D:\!!!Research\数据集-生物信息学分析"
TAB = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES", "SUPPLEMENTARY_Tables_CSV")
INTER = os.path.join(ROOT, "_intermediate")
LOG = os.path.join(ROOT, "03_LOGS", "M16_step4_log.txt")
VISIUM = r"D:\!!!Research\!!!课题组\ARDS核心基因筛选策略方案\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\归档\SCI论文1黄裕荣_Mitoxyperilysis_ARDS-20260819\_intermediate\M1_visium_scored.h5ad"

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
log = []
def note(m=""):
    log.append(m); print(m, flush=True)

def spearman(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 8:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return r, p, int(m.sum())

def partial_spearman(x, y, covs):
    """以秩转换残差化实现偏 Spearman（covs: 2D 数组）。"""
    m = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(covs), axis=1)
    x, y, covs = x[m], y[m], covs[m]
    if len(x) < 12:
        return np.nan, np.nan, len(x)
    rx = stats.rankdata(x); ry = stats.rankdata(y)
    rc = np.vstack([stats.rankdata(covs[:, j]) for j in range(covs.shape[1])]).T
    X = np.column_stack([np.ones(len(rx)), rc])
    bx = np.linalg.lstsq(X, rx, rcond=None)[0]; by = np.linalg.lstsq(X, ry, rcond=None)[0]
    ex = rx - X @ bx; ey = ry - X @ by
    r = np.corrcoef(ex, ey)[0, 1]
    n = len(rx); k = covs.shape[1]
    t = r * np.sqrt((n - 2 - k) / (1 - r ** 2))
    p = 2 * stats.t.sf(abs(t), n - 2 - k)
    return r, p, n

def symbol2ensembl(symbols):
    cache = os.path.join(INTER, "M16_symbol2ensembl.json")
    if os.path.exists(cache):
        return json.load(open(cache, encoding="utf-8"))
    q = ",".join(symbols)
    data = urllib.parse.urlencode({"q": q, "scopes": "symbol", "fields": "symbol,ensembl.gene",
                                   "species": "human", "size": str(len(symbols) * 3)}).encode()
    for i in range(4):
        try:
            req = urllib.request.Request("https://mygene.info/v3/query", data=data,
                                         headers={"User-Agent": "M16/1.0"})
            res = json.loads(urllib.request.urlopen(req, timeout=60).read())
            break
        except Exception as e:
            note(f"  [retry {i+1}] mygene: {e}"); time.sleep(3)
    else:
        raise RuntimeError("mygene symbol2ensembl failed")
    out = {}
    for h in res:
        if h.get("notfound"):
            continue
        e = h.get("ensembl", {})
        gid = e.get("gene") if isinstance(e, dict) else (e[0].get("gene") if e else None)
        out[str(h.get("symbol", "")).upper()] = gid
    json.dump(out, open(cache, "w", encoding="utf-8"))
    return out

def load_mech():
    mech = pd.read_csv(os.path.join(TAB, "Table_S88_M16_Mechanosensing_Module_v10.csv"), encoding="utf-8-sig")
    return sorted(mech["hgnc_symbol"])

def mech_score_from_z(Zdf, genes_upper):
    idx = [g for g in Zdf.index if g.upper() in genes_upper]
    return Zdf.loc[idx].mean(), len(idx)

# ---------------- (a) GSE185263 ----------------
def bridge_gse185263(mech_genes, s56):
    note("[a] GSE185263 全血 bulk 桥接")
    cnt = pd.read_csv(os.path.join(DSET, "GSE185263_Lung_ARDS", "GSE185263_raw_counts.csv"),
                      index_col=0)
    cnt.index = cnt.index.astype(str).str.split(".").str[0]
    s2e = symbol2ensembl(mech_genes)
    ens2sym = {v: k for k, v in s2e.items() if v}
    hit = [e for e in cnt.index if e in ens2sym]
    sub = cnt.loc[hit]
    sub.index = [ens2sym[e] for e in hit]
    lib = cnt.sum(axis=0)
    cpm = sub.divide(lib, axis=1) * 1e4  # CPM(10,000) 与项目口径一致
    lg = np.log2(cpm + 1)
    Z = lg.sub(lg.mean(axis=1), axis=0).div(lg.std(axis=1, ddof=1), axis=0)
    mech, nh = mech_score_from_z(Z, set(s2e.keys()))
    note(f"  mech 基因命中 {nh}/23")
    df = s56[s56["cohort"] == "GSE185263"].copy()
    df["mech"] = mech.reindex(df["sample"]).values
    # 组成协变量（Monaco|NNLS）
    s64 = pd.read_csv(os.path.join(TAB, "Table_S64_M10A_Deconvolution_Proportions.csv"))
    s64 = s64[(s64["cohort"] == "GSE185263") & (s64["reference"] == "Monaco") & (s64["method"] == "NNLS")]
    s64["nm"] = s64["Neutrophils"] + s64["C_mono"] + s64["I_mono"] + s64["NC_mono"]
    df = df.merge(s64[["sample", "nm"]], on="sample", how="left")
    rows = []
    for target in ["UCS", "EIS", "MDI"]:
        r, p, n = spearman(df["mech"].values, df[target].values)
        rp, pp, np_ = partial_spearman(df["mech"].values, df[target].values,
                                       df[["nm"]].fillna(0).values)
        rows.append({"layer": "a_bulk_blood", "cohort": "GSE185263", "target": target,
                     "spearman_rho": r, "p_two": p, "n": n,
                     "partial_rho_adj_composition": rp, "partial_p_two": pp, "n_partial": np_,
                     "composition_covariate": "Monaco|NNLS neutrophil+monocyte (Table_S64)"})
        note(f"  mech×{target}: ρ={r:+.3f} (p={p:.3g}, n={n}); 校正后 ρ={rp:+.3f} (p={pp:.3g})")
    return rows

# ---------------- (b) GSE212865 ----------------
def bridge_gse212865(mech_genes, s56):
    note("[b] GSE212865 外周血芯片桥接")
    # 探针→symbol：GPL23159 family soft
    m = {}
    with gzip.open(os.path.join(DSET, "GPL_Annotations", "GPL23159_family.soft.gz"),
                   "rt", encoding="utf-8", errors="replace") as fh:
        txt = fh.read()
    sec = txt.split("!platform_table_begin", 1)[1].split("!platform_table_end", 1)[0]
    lines = sec.strip().split("\n")
    hdr = lines[0].split("\t")
    i_id = hdr.index("ID")
    # 注释长列（含 " // " 的基因指派文本）为无表头列，按内容识别
    sym_re = re.compile(r"\(([^()]+)\),\s*(?:mRNA|transcript)")
    for l in lines[1:]:
        c = l.split("\t")
        if len(c) <= i_id:
            continue
        ga = next((x for x in c if " // " in x), "")
        if not ga:
            continue
        first = ga.split(" /// ")[0]
        mm = sym_re.search(first)
        if mm:
            m[c[i_id]] = mm.group(1).strip()
    note(f"  GPL23159 探针→symbol {len(m)}")
    f = os.path.join(DSET, "GSE212865_data", "GSE212865_series_matrix.txt.gz")
    with gzip.open(f, "rt", encoding="utf-8", errors="replace") as fh:
        txt = fh.read()
    samples = re.search(r'^!Sample_geo_accession\t(.+)$', txt, re.M).group(1).replace('"', "").split("\t")
    body = txt.split("!series_matrix_table_begin\n", 1)[1].split("!series_matrix_table_end", 1)[0]
    lines = body.strip().split("\n")
    ids = [l.split("\t")[0].strip().strip('"') for l in lines[1:]]
    vals = np.array([[float(x) if x not in ("", "null", "NA") else np.nan
                      for x in l.split("\t")[1:]] for l in lines[1:]])
    mat = pd.DataFrame(vals, index=ids, columns=samples)
    mx = np.nanmax(vals)
    if mx > 100:
        note(f"  线性尺度(max={mx:.0f}) → log2(x+1)"); mat = np.log2(mat + 1.0)
    else:
        note(f"  已为 log 尺度(max={mx:.1f})")
    sym = pd.Series([m.get(p, "") for p in mat.index], index=mat.index)
    mat = mat.loc[sym != ""]; sym = sym[sym != ""]
    mat["sym"] = sym
    gm = mat.groupby("sym").mean()
    Z = gm.sub(gm.mean(axis=1), axis=0).div(gm.std(axis=1, ddof=1), axis=0)
    mech, nh = mech_score_from_z(Z, {g.upper() for g in mech_genes})
    note(f"  mech 基因命中 {nh}/23；样本 {mat.shape[1]}")
    df = s56[s56["cohort"] == "GSE212865"].copy()
    df["mech"] = mech.reindex(df["sample"]).values
    rows = []
    for target in ["UCS", "EIS", "MDI"]:
        r, p, n = spearman(df["mech"].values, df[target].values)
        rows.append({"layer": "b_pbmc_microarray", "cohort": "GSE212865", "target": target,
                     "spearman_rho": r, "p_two": p, "n": n,
                     "partial_rho_adj_composition": np.nan, "partial_p_two": np.nan, "n_partial": 0,
                     "composition_covariate": ""})
        note(f"  mech×{target}: ρ={r:+.3f} (p={p:.3g}, n={n})")
    return rows

# ---------------- (c) GSE158055 供体级 ----------------
def bridge_gse158055(mech_genes):
    note("[c] GSE158055 scRNA 供体级桥接（单核亚集 h5ad 内自洽重算臂分数+mech）")
    import anndata as ad
    a = ad.read_h5ad(os.path.join(RAW, "GSE158055_scRNA", "GSE158055_Mono_Processed.h5ad"))
    genes_upper = a.var_names.str.upper()
    man = pd.read_csv(os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "Mitoxyperilysis_Gene_Manifest_v1.0.csv"),
                      encoding="utf-8-sig")
    up = set(man.loc[man["arm"] == "upstream_collapse", "hgnc_symbol"].str.upper())
    ex = set(man.loc[man["arm"] == "execution_induction", "hgnc_symbol"].str.upper())
    mech_set = {g.upper() for g in mech_genes}
    note(f"  h5ad {a.shape}；臂命中 up {genes_upper.isin(up).sum()}/30 ex {genes_upper.isin(ex).sum()}/33 "
         f"mech {genes_upper.isin(mech_set).sum()}/23")
    keep = genes_upper.isin(up | ex | mech_set)
    X = a[:, keep].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    gnames = list(genes_upper[keep])
    obs = a.obs[["sampleID", "celltype"]].reset_index(drop=True)
    df = pd.DataFrame(X, columns=gnames)
    df["sampleID"] = obs["sampleID"].values
    df["donor"] = obs["sampleID"].str.replace(r"-\d+$", "", regex=True).values
    df["celltype"] = obs["celltype"].values
    # 供体×细胞类型均值（先按供体塌陷重复采样，再按供体×类型聚合）
    agg = df.groupby(["donor", "celltype"])[gnames].mean()
    ncell = df.groupby(["donor", "celltype"]).size().rename("n_cells")
    agg = agg.join(ncell)
    agg = agg[agg["n_cells"] >= 50]
    # 数据集内基因 z（跨供体×类型层面 z 化；log 空间假定与 h5ad 预处理一致，披露）
    Zc = agg[gnames].sub(agg[gnames].mean()).div(agg[gnames].std(ddof=1))
    def arm_mean(s):
        cols = [g for g in gnames if g in s]
        return Zc[cols].mean(axis=1)
    agg["UCS"] = arm_mean(up); agg["EIS"] = arm_mean(ex)
    agg["MDI"] = agg["EIS"] - agg["UCS"]
    agg["mech"] = arm_mean(mech_set)
    note(f"  供体×类型单元 {len(agg)}（供体 {agg.reset_index()['donor'].nunique()}）")
    rows = []
    for layer, sub in (("c_scRNA_donor_monocytes_all", agg),):
        for target in ["UCS", "EIS", "MDI"]:
            r, p, n = spearman(sub["mech"].values, sub[target].values)
            rows.append({"layer": layer, "cohort": "GSE158055", "target": target,
                         "spearman_rho": r, "p_two": p, "n": n,
                         "partial_rho_adj_composition": np.nan, "partial_p_two": np.nan, "n_partial": 0,
                         "composition_covariate": "供体×细胞类型单元内自洽重算；非 S67 合并"})
            note(f"  供体级 mech×{target}: ρ={r:+.3f} (p={p:.3g}, n={n})")
    return rows

# ---------------- (d) Visium 空间共定位（探索性，家族外） ----------------
def bridge_visium(mech_genes, seed=20260831):
    note("[d] Visium 空间共定位（探索性）")
    import anndata as ad
    a = ad.read_h5ad(VISIUM)
    genes = a.var_names.str.upper()
    keep = genes.isin({g.upper() for g in mech_genes})
    note(f"  Visium {a.shape}；mech 命中 {keep.sum()}/23")
    X = a[:, keep].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    obs = a.obs[["section", "condition", "up_score_z", "ex_score_z", "mdi_z"]].copy()
    coords = a.obsm["spatial"]
    del a
    df = pd.DataFrame(X, columns=genes[keep])
    df["section"] = obs["section"].values
    rows = []
    rng_master = np.random.default_rng(seed)
    pooled = {}
    for sec, idx in df.groupby("section").groups.items():
        idx = np.asarray(idx)
        sub = df.loc[idx]
        Zs = sub.drop(columns=["section"])
        Zs = Zs.sub(Zs.mean()).div(Zs.std(ddof=1).replace(0, np.nan))
        mech = Zs.mean(axis=1).values
        up = obs["up_score_z"].values[idx]; ex = obs["ex_score_z"].values[idx]; mdi = obs["mdi_z"].values[idx]
        n = len(idx)
        r_up, p_up, _ = spearman(mech, up)
        r_ex, p_ex, _ = spearman(mech, ex)
        r_mdi, p_mdi, _ = spearman(mech, mdi)
        # 双变量 Moran's I（kNN k=8 行标准化，999 置换）
        xy = coords[idx]
        d2 = ((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1)
        np.fill_diagonal(d2, np.inf)
        knn = np.argsort(d2, axis=1)[:, :8]
        W = np.zeros((n, n)); W[np.arange(n)[:, None], knn] = 1.0
        W = W / W.sum(1, keepdims=True)
        def bmi(x, y):
            xz = (x - x.mean()) / x.std(); yz = (y - y.mean()) / y.std()
            return float(xz @ W @ yz / n * n / n)  # = xz W yz / n（行标准化 W）
        def bmi_stat(x, y):
            xz = (x - x.mean()) / x.std(); yz = (y - y.mean()) / y.std()
            return float((xz * (W @ yz)).mean())
        i_ex = bmi_stat(mech, ex); i_up = bmi_stat(mech, up)
        B = 999
        ge_ex = ge_up = 1
        for b in range(B):
            perm = rng_master.permutation(n)
            if abs(bmi_stat(mech[perm], ex)) >= abs(i_ex): ge_ex += 1
            if abs(bmi_stat(mech[perm], up)) >= abs(i_up): ge_up += 1
        rows.append({"section": sec, "condition": obs["condition"].values[idx][0], "n_spots": n,
                     "rho_mech_UCS": r_up, "p_mech_UCS": p_up,
                     "rho_mech_EIS": r_ex, "p_mech_EIS": p_ex,
                     "rho_mech_MDI": r_mdi, "p_mech_MDI": p_mdi,
                     "bivarI_mech_EIS": i_ex, "perm_p_mech_EIS": ge_ex / (B + 1),
                     "bivarI_mech_UCS": i_up, "perm_p_mech_UCS": ge_up / (B + 1)})
        note(f"  [{sec}] n={n} ρ(mech,EIS)={r_ex:+.3f} I(mech,EIS)={i_ex:+.4f} (p={ge_ex/(B+1):.4f})")
    res = pd.DataFrame(rows)
    # Stouffer 合并（双侧）
    def stouffer(vals, ps):
        zs = stats.norm.isf(np.clip(ps, 1e-300, 1 - 1e-16) / 2) * np.sign(vals)
        z = zs.sum() / np.sqrt(len(zs))
        return z, 2 * stats.norm.sf(abs(z))
    for pair in [("rho_mech_UCS", "p_mech_UCS"), ("rho_mech_EIS", "p_mech_EIS"),
                 ("rho_mech_MDI", "p_mech_MDI"), ("bivarI_mech_EIS", "perm_p_mech_EIS"),
                 ("bivarI_mech_UCS", "perm_p_mech_UCS")]:
        z, p = stouffer(res[pair[0]].values, res[pair[1]].values)
        pooled[pair[0]] = {"stouffer_z": z, "p_two": p,
                           "n_positive": int((np.sign(res[pair[0]]) > 0).sum()), "k": len(res)}
        note(f"  [pooled] {pair[0]}: Stouffer z={z:+.2f} (p={p:.3g})，正值 {pooled[pair[0]]['n_positive']}/{len(res)}")
    fv = os.path.join(INTER, "M16_visium_mech_spatial.csv")
    res.to_csv(fv, index=False, encoding="utf-8-sig")
    note(f"[write] {fv}")
    return res, pooled

def main():
    note("== M16 step4 人体桥接层 2026-08-31（预注册判据 D） ==")
    os.makedirs(TAB, exist_ok=True); os.makedirs(INTER, exist_ok=True)
    mech_genes = load_mech()
    s56 = pd.read_csv(os.path.join(TAB, "Table_S56_M2_PerSample_Scores.csv"))
    rows = []
    rows += bridge_gse185263(mech_genes, s56)
    rows += bridge_gse212865(mech_genes, s56)
    rows += bridge_gse158055(mech_genes)
    df = pd.DataFrame(rows)
    # BH 家族 = (a)(b)(c) 全部主相关
    df = df.sort_values("p_two")
    m = len(df)
    df["bh_rank"] = np.arange(1, m + 1)
    df["p_two_bh"] = np.minimum.accumulate((df["p_two"] * m / df["bh_rank"])[::-1])[::-1].clip(upper=1)
    try:
        vres, pooled = bridge_visium(mech_genes)
        for k, v in pooled.items():
            row = {"layer": "d_visium_spatial_pooled", "cohort": "GSE271370", "target": k,
                   "spearman_rho": v["stouffer_z"], "p_two": v["p_two"], "n": v["k"],
                   "partial_rho_adj_composition": np.nan, "partial_p_two": np.nan, "n_partial": 0,
                   "composition_covariate": "Stouffer z（非 ρ）；空间层探索性，BH 家族外",
                   "p_two_bh": np.nan, "bh_rank": np.nan}
            df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
    except Exception as e:
        note(f"[d] Visium 层失败：{e}（如实披露，不重试伪造）")
    f91 = os.path.join(TAB, "Table_S91_M16_Bridge_Correlations.csv")
    df.to_csv(f91, index=False, encoding="utf-8-sig")
    note(f"[write] {f91}")
    with open(LOG, "w", encoding="utf-8") as f:
        f.write("\n".join(log) + "\n")
    print("[done]")

if __name__ == "__main__":
    main()
