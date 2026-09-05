# -*- coding: utf-8 -*-
"""
M16_step2_vili_scoring.py — M16C：4 套小鼠 VILI 数据集 mdi_v1.0 同标准重算 + mech_v1.0 评分
=====================================================================
预注册：M16_pre_registration_20260831.md（判据 A/C）
- 平台探针→小鼠 symbol：GEO 官方 GPL annot（本地 00_RAW_DATA/GEO_downloads，GPL339/5145/8321/1261）
- 小鼠 symbol → 人源臂/mech 基因：M16_human2mouse_orthologs.csv（HomoloGene 双向一对一，step1）
- 口径：数据集内 log2 空间逐基因 z（全样本、不分组）→ 臂均值；UCS=上游臂均值、EIS=执行臂均值、
  MDI=EIS−UCS；nomt 敏感性（剔除 MT-*）；mech=23 基因模块均值 z。
- 值域：max>100 判为线性尺度 → log2(x+1)（GSE2411/GSE9208）；GSE7742/GSE9368 已为 log2 尺度（披露）。
- 主对比：GSE2411 MV vs Control；GSE7742 WT vent vs WT ctrl；GSE9208 WT MV vs WT SpV；
  GSE9368 VALI vs Control。修饰层（家族外敏感性）：jnk1-/-、Nrf2-/-、rhPBEF 各层。
- 统计：Welch t + Mann-Whitney + Hedges' g（正值=通气组高）；单侧方向预测 UCS<0/EIS>0/MDI>0/mech>0；
  BH 家族=4 数据集×4 分数=16；合并=DL 随机效应 + Hartung–Knapp + 留一数据集。
输出：Table_S89_M16_VILI_Contrasts.csv、Table_S92_M16_Ortholog_Coverage_Audit.csv、
  _intermediate/M16_vili_per_sample_scores.csv、03_LOGS/M16_step2_log.txt
"""
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.abspath(__file__))
RAW = r"D:\!!!Research\!!!课题组\ARDS核心基因筛选策略方案\00_RAW_DATA\GEO_downloads"
TAB = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES", "SUPPLEMENTARY_Tables_CSV")
INTER = os.path.join(ROOT, "_intermediate")
LOG = os.path.join(ROOT, "03_LOGS", "M16_step2_log.txt")
MANIFEST = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "Mitoxyperilysis_Gene_Manifest_v1.0.csv")
ORTH = os.path.join(INTER, "M16_human2mouse_orthologs.csv")

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
log = []
def note(m=""):
    log.append(m); print(m, flush=True)

MT = {"MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CYB", "MT-ND1", "MT-ND2"}

# ---------------- GPL annot 解析（探针→小鼠 symbol） ----------------
def load_gpl_map(gpl):
    f = os.path.join(RAW, f"{gpl}.annot.gz")
    with gzip.open(f, "rt", encoding="utf-8", errors="replace") as fh:
        lines = [l for l in fh if not l.startswith(("^", "!", "#"))]
    hdr = lines[0].rstrip("\n").split("\t")
    i_id = hdr.index("ID"); i_sym = hdr.index("Gene symbol")
    m = {}
    for l in lines[1:]:
        c = l.rstrip("\n").split("\t")
        if len(c) <= max(i_id, i_sym):
            continue
        sym = c[i_sym].split("///")[0].strip()
        if sym:
            m[c[i_id].strip()] = sym
    note(f"[{gpl}] 探针→symbol {len(m)}")
    return m

# ---------------- series matrix 解析 ----------------
def load_series(gse):
    f = os.path.join(RAW, f"{gse}_series_matrix.txt.gz")
    with gzip.open(f, "rt", encoding="utf-8", errors="replace") as fh:
        txt = fh.read()
    samples = re.search(r'^!Sample_geo_accession\t(.+)$', txt, re.M).group(1).replace('"', "").split("\t")
    titles = re.search(r'^!Sample_title\t(.+)$', txt, re.M).group(1).replace('"', "").split("\t")
    body = txt.split("!series_matrix_table_begin\n", 1)[1].split("!series_matrix_table_end", 1)[0]
    lines = body.strip().split("\n")
    ids = [l.split("\t")[0].strip().strip('"') for l in lines[1:]]
    vals = np.array([[float(x) if x not in ("", "null", "NA") else np.nan
                      for x in l.split("\t")[1:]] for l in lines[1:]])
    mat = pd.DataFrame(vals, index=ids, columns=samples)
    mx = np.nanmax(vals)
    if mx > 100:
        note(f"[{gse}] 线性尺度(max={mx:.0f}) → log2(x+1)")
        mat = np.log2(mat + 1.0)
    else:
        note(f"[{gse}] 已为 log2 尺度(max={mx:.1f})，原样使用（披露）")
    return mat, dict(zip(samples, titles))

# ---------------- 分组（按官方样本标题，先于任何数值） ----------------
def groups(gse, titles):
    g = {}
    if gse == "GSE2411":
        for s, t in titles.items():
            if t.startswith("Control"): g[s] = "Control"
            elif t.startswith("MV+LPS"): g[s] = "MV+LPS"
            elif t.startswith("MV"): g[s] = "MV"
            elif t.startswith("LPS"): g[s] = "LPS"
    elif gse == "GSE7742":
        for s, t in titles.items():
            ko = "jnk1ko" in t.lower()
            vent = "ventilation" in t.lower()
            g[s] = ("KO_" if ko else "WT_") + ("vent" if vent else "ctrl")
    elif gse == "GSE9208":
        for s, t in titles.items():
            ko = t.startswith("KO")
            mv = "mechanical" in t.lower()
            g[s] = ("KO_" if ko else "WT_") + ("MV" if mv else "SpV")
    elif gse == "GSE9368":
        for s, t in titles.items():
            tl = t.lower()
            if tl.startswith("control"): g[s] = "Control"
            elif "rhpbef" in tl and "vali" in tl: g[s] = "rhPBEF+VALI"
            elif tl.startswith("vali"): g[s] = "VALI"
            elif tl.startswith("rhpbef"): g[s] = "rhPBEF"
    return g

# ---------------- 统计 ----------------
def hedges_g(ctrl, case):
    ctrl = np.asarray(ctrl, float); case = np.asarray(case, float)
    n1, n2 = len(ctrl), len(case)
    if n1 < 2 or n2 < 2:
        return np.nan
    d = (case.mean() - ctrl.mean()) / np.sqrt(((n1 - 1) * ctrl.var(ddof=1) + (n2 - 1) * case.var(ddof=1)) / (n1 + n2 - 2))
    return d * (1 - 3 / (4 * (n1 + n2) - 9))

def hedges_se(g, n1, n2):
    return np.sqrt((n1 + n2) / (n1 * n2) + g ** 2 / (2 * (n1 + n2)))

def one_sided(p_two, stat, direction):
    """方向性单侧 p：direction=-1 → H1: case<ctrl（stat 为 case-ctrl 读数，负值为支持方向）。"""
    if direction == -1:
        return p_two / 2 if stat < 0 else 1 - p_two / 2
    return p_two / 2 if stat > 0 else 1 - p_two / 2

DIRECTION = {"UCS": -1, "UCS_nomt": -1, "EIS": 1, "MDI": 1, "MDI_nomt": 1, "mech": 1}

def dl_hk_meta(gs, ses):
    """DerSimonian-Laird 随机效应 + Hartung–Knapp，单侧（方向=各层方向一致的符号）。"""
    gs = np.asarray(gs, float); ses = np.asarray(ses, float)
    k = len(gs); w = 1 / ses ** 2
    mu_fe = np.sum(w * gs) / np.sum(w)
    q = np.sum(w * (gs - mu_fe) ** 2)
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (q - (k - 1)) / c)
    wr = 1 / (ses ** 2 + tau2)
    mu = np.sum(wr * gs) / np.sum(wr)
    se_hk = np.sqrt(np.sum(wr * (gs - mu) ** 2) / ((k - 1) * np.sum(wr)))
    t = mu / se_hk
    p_two = 2 * stats.t.sf(abs(t), k - 1)
    sign = 1 if np.mean(np.sign(gs)) > 0 else -1
    p_one = one_sided(p_two, t * sign, 1)
    return dict(k=k, pooled=mu, se=se_hk, t=t, p_one=p_one, tau2=tau2,
                i2=max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0,
                dir_consistent=int(np.sum(np.sign(gs) == sign)))

def main():
    note("== M16 step2 VILI mdi_v1.0 同标准重算 2026-08-31（预注册 M16_pre_registration_20260831.md）==")
    os.makedirs(TAB, exist_ok=True); os.makedirs(INTER, exist_ok=True)

    man = pd.read_csv(MANIFEST, encoding="utf-8-sig")
    up = sorted(man.loc[man["arm"] == "upstream_collapse", "hgnc_symbol"].str.strip())
    ex = sorted(man.loc[man["arm"] == "execution_induction", "hgnc_symbol"].str.strip())
    mech = pd.read_csv(os.path.join(TAB, "Table_S88_M16_Mechanosensing_Module_v10.csv"), encoding="utf-8-sig")
    mech_genes = sorted(mech["hgnc_symbol"])
    orth = pd.read_csv(ORTH, encoding="utf-8-sig").fillna("")
    h2m = {r["human_symbol"].upper(): r["mouse_symbol"] for _, r in orth.iterrows() if r["mouse_symbol"]}

    SETS = [
        ("GSE2411", "GPL339", "MV", "Control", [("LPS", "Control"), ("MV+LPS", "Control"), ("MV+LPS", "LPS")]),
        ("GSE7742", "GPL5145", "WT_vent", "WT_ctrl", [("KO_vent", "KO_ctrl")]),
        ("GSE9208", "GPL8321", "WT_MV", "WT_SpV", [("KO_MV", "KO_SpV")]),
        ("GSE9368", "GPL1261", "VALI", "Control", [("rhPBEF+VALI", "VALI"), ("rhPBEF", "Control")]),
    ]

    all_rows, cover_rows, per_sample_all = [], [], []
    for gse, gpl, prim_case, prim_ctrl, mods in SETS:
        mat, titles = load_series(gse)
        grp = groups(gse, titles)
        pmap = load_gpl_map(gpl)
        # 探针→小鼠 symbol → 塌陷（max-mean，≤3 探针）
        sym = pd.Series([pmap.get(p, "") for p in mat.index], index=mat.index)
        mat = mat.loc[sym != ""]
        sym = sym[sym != ""]
        df = mat.copy(); df["sym"] = sym
        gm = df.groupby("sym").mean()  # 小鼠 symbol 级
        mu = gm.mean(axis=1)
        keep = mu.groupby(level=0)  # 单 symbol 已聚合
        # 人源映射层：仅保留有一对一同源的基因
        m2h = {v: k for k, v in h2m.items()}
        gm_h = gm.loc[[s for s in gm.index if s in m2h]]
        gm_h.index = [m2h[s] for s in gm_h.index]
        note(f"[{gse}] 人源可映射基因 {len(gm_h)}")
        # 覆盖率审计
        for arm_name, genes in (("UCS", up), ("EIS", ex), ("mech", mech_genes)):
            hit = [g for g in genes if g.upper() in set(gm_h.index.str.upper())]
            cover_rows.append({"gse": gse, "set": arm_name, "n_total": len(genes), "n_present": len(hit),
                               "missing": ";".join(sorted(set(genes) - set(hit)))})
        idx_up = [g for g in gm_h.index if g.upper() in {x.upper() for x in up}]
        idx_ex = [g for g in gm_h.index if g.upper() in {x.upper() for x in ex}]
        idx_up_nomt = [g for g in idx_up if g.upper() not in MT]
        idx_mech = [g for g in gm_h.index if g.upper() in {x.upper() for x in mech_genes}]
        Z = (gm_h - gm_h.mean(axis=1).values[:, None]) / gm_h.std(axis=1, ddof=1).values[:, None]
        sc = pd.DataFrame(index=gm_h.columns)
        sc["UCS"] = Z.loc[idx_up].mean()
        sc["EIS"] = Z.loc[idx_ex].mean()
        sc["MDI"] = sc["EIS"] - sc["UCS"]
        sc["UCS_nomt"] = Z.loc[idx_up_nomt].mean()
        sc["MDI_nomt"] = sc["EIS"] - sc["UCS_nomt"]
        sc["mech"] = Z.loc[idx_mech].mean()
        sc["group"] = pd.Series(grp)
        sc["gse"] = gse
        per_sample_all.append(sc.reset_index().rename(columns={"index": "sample"}))
        # 对比
        contrasts = [("primary", prim_case, prim_ctrl)] + [("sensitivity", a, b) for a, b in mods if a]
        for layer, case_g, ctrl_g in contrasts:
            ca = sc.loc[sc["group"] == case_g]; ct = sc.loc[sc["group"] == ctrl_g]
            if len(ca) < 2 or len(ct) < 2:
                note(f"[{gse}] {case_g} vs {ctrl_g} 样本不足，跳过"); continue
            for score in ["UCS", "EIS", "MDI", "mech", "UCS_nomt", "MDI_nomt"]:
                a, b = ca[score].values, ct[score].values
                t, p_t = stats.ttest_ind(a, b, equal_var=False)
                u, p_u = stats.mannwhitneyu(a, b, alternative="two-sided")
                g = hedges_g(b, a)
                d = DIRECTION.get(score, 1)
                all_rows.append({
                    "gse": gse, "layer": layer, "contrast": f"{case_g}_vs_{ctrl_g}",
                    "score": score, "n_case": len(a), "n_ctrl": len(b),
                    "mean_case": a.mean(), "mean_ctrl": b.mean(), "diff": a.mean() - b.mean(),
                    "hedges_g": g, "welch_t": t, "welch_p_two": p_t, "mw_p_two": p_u,
                    "welch_p_one": one_sided(p_t, t, d), "mw_p_one": one_sided(p_u, a.mean() - b.mean(), d),
                    "direction_hypothesis": { -1: "decrease", 1: "increase"}[d],
                    "n_genes_UCS": len(idx_up), "n_genes_EIS": len(idx_ex), "n_genes_mech": len(idx_mech),
                })
        note(f"[{gse}] 主对比 {prim_case} vs {prim_ctrl} 完成；臂基因覆盖 UCS {len(idx_up)}/30 EIS {len(idx_ex)}/33 mech {len(idx_mech)}/23")

    res = pd.DataFrame(all_rows)
    # BH 家族=4 数据集×4 分数（主对比层）
    prim = res[(res["layer"] == "primary") & res["score"].isin(["UCS", "EIS", "MDI", "mech"])].copy()
    prim = prim.sort_values("welch_p_one")
    m = len(prim)
    prim["bh_rank"] = np.arange(1, m + 1)
    prim["welch_p_one_bh"] = np.minimum.accumulate(
        (prim["welch_p_one"] * m / prim["bh_rank"])[::-1])[::-1].clip(upper=1)
    res = res.merge(prim[["gse", "contrast", "score", "welch_p_one_bh"]],
                    on=["gse", "contrast", "score"], how="left")

    # 合并（主对比层，REML+HK）
    meta_rows = []
    for score in ["UCS", "EIS", "MDI", "mech"]:
        sub = prim[prim["score"] == score]
        gs = sub["hedges_g"].values
        ses = [hedges_se(g, n1, n2) for g, n1, n2 in zip(sub["hedges_g"], sub["n_case"], sub["n_ctrl"])]
        mm = dl_hk_meta(gs, ses)
        loso_dirs = []
        for i in range(len(gs)):
            mm2 = dl_hk_meta(np.delete(gs, i), np.delete(ses, i))
            loso_dirs.append(int(np.sign(mm2["pooled"]) == np.sign(mm["pooled"])))
        meta_rows.append({"score": score, **mm, "loso_same_direction": f"{sum(loso_dirs)}/{len(loso_dirs)}"})
        note(f"[meta] {score}: pooled g={mm['pooled']:+.3f} (k={mm['k']}, I2={mm['i2']:.2f}, "
             f"单侧 p={mm['p_one']:.4f}, 方向一致 {mm['dir_consistent']}/{mm['k']}, LOSO {sum(loso_dirs)}/{len(loso_dirs)})")
    meta = pd.DataFrame(meta_rows)

    f89 = os.path.join(TAB, "Table_S89_M16_VILI_Contrasts.csv")
    res.to_csv(f89, index=False, encoding="utf-8-sig")
    f92 = os.path.join(TAB, "Table_S92_M16_Ortholog_Coverage_Audit.csv")
    pd.DataFrame(cover_rows).to_csv(f92, index=False, encoding="utf-8-sig")
    fmeta = os.path.join(INTER, "M16_vili_meta.csv")
    meta.to_csv(fmeta, index=False, encoding="utf-8-sig")
    fps = os.path.join(INTER, "M16_vili_per_sample_scores.csv")
    pd.concat(per_sample_all).to_csv(fps, index=False, encoding="utf-8-sig")
    note(f"[write] {f89}\n[write] {f92}\n[write] {fmeta}\n[write] {fps}")
    with open(LOG, "w", encoding="utf-8") as f:
        f.write("\n".join(log) + "\n")
    print("[done]")

if __name__ == "__main__":
    main()
