# -*- coding: utf-8 -*-
"""
M16_step3_gse2411_interaction.py — M16C 判据 B：GSE2411 力学×LPS 2×2 互作正式检验
=====================================================================
预注册：M16_pre_registration_20260831.md（判据 B）
- 设计：Control / MV / LPS / MV+LPS（各 n=6）；互作 β₃ = μ_MV+LPS − μ_MV − μ_LPS + μ_Control
- 分数：UCS / EIS / MDI / mech（mdi_v1.0 同口径，step2 逐样本分数）
- 推断：Welch-Satterthwaite 对比 t（主）+ 标签置换经验 p（B=10,000，固定种子，四组内整体洗牌）
- 方向预测：MDI β₃>0（力学放大解离）；EIS β₃>0；UCS β₃<0；mech β₃ 双向报告
- BH 家族 = 4 分数；基因级互作（臂基因+mech 基因）为描述性层
输出：Table_S90_M16_GSE2411_Interaction.csv、_intermediate/M16_gse2411_genelevel_interaction.csv、03_LOGS/M16_step3_log.txt
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES", "SUPPLEMENTARY_Tables_CSV")
INTER = os.path.join(ROOT, "_intermediate")
LOG = os.path.join(ROOT, "03_LOGS", "M16_step3_log.txt")

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
log = []
def note(m=""):
    log.append(m); print(m, flush=True)

DIRECTION = {"UCS": -1, "EIS": 1, "MDI": 1, "mech": 0}  # 0=双向

def welch_contrast(groups, order=("Control", "MV", "LPS", "MV+LPS")):
    """β₃ = μ4 − μ2 − μ3 + μ1，Welch-Satterthwaite df。"""
    vals = [np.asarray(groups[g], float) for g in order]
    ms = [v.mean() for v in vals]
    vs = [v.var(ddof=1) for v in vals]
    ns = [len(v) for v in vals]
    beta = ms[3] - ms[1] - ms[2] + ms[0]
    se = np.sqrt(vs[3] / ns[3] + vs[1] / ns[1] + vs[2] / ns[2] + vs[0] / ns[0])
    df_num = (vs[3] / ns[3] + vs[1] / ns[1] + vs[2] / ns[2] + vs[0] / ns[0]) ** 2
    df_den = sum((v / n) ** 2 / (n - 1) for v, n in zip(vs, ns))
    dfree = df_num / df_den if df_den > 0 else np.nan
    t = beta / se if se > 0 else np.nan
    p_two = 2 * stats.t.sf(abs(t), dfree) if np.isfinite(t) else np.nan
    return beta, se, t, dfree, p_two

def perm_p(values, labels, B=10000, seed=20260831):
    """四组标签整体洗牌的经验 p（双侧 + 方向性）。"""
    rng = np.random.default_rng(seed)
    order = ["Control", "MV", "LPS", "MV+LPS"]
    obs, *_ = welch_contrast({g: values[labels == g] for g in order})
    vals = np.asarray(values, float); lab = np.asarray(labels)
    cnt_two = cnt_pos = 0
    for _ in range(B):
        p = rng.permutation(lab)
        b, *_ = welch_contrast({g: vals[p == g] for g in order})
        if abs(b) >= abs(obs) - 1e-15: cnt_two += 1
        if b >= obs - 1e-15: cnt_pos += 1
    return obs, (1 + cnt_two) / (1 + B), (1 + cnt_pos) / (1 + B)

def one_sided(p_two, stat, d):
    if d == 0:
        return p_two
    if d == -1:
        return p_two / 2 if stat < 0 else 1 - p_two / 2
    return p_two / 2 if stat > 0 else 1 - p_two / 2

def main():
    note("== M16 step3 GSE2411 力学×LPS 互作正式检验（判据 B）2026-08-31 ==")
    sc = pd.read_csv(os.path.join(INTER, "M16_vili_per_sample_scores.csv"))
    sc = sc[sc["gse"] == "GSE2411"].copy()
    order = ["Control", "MV", "LPS", "MV+LPS"]
    note(f"[GSE2411] n={len(sc)}，各组 {sc['group'].value_counts().to_dict()}")

    rows = []
    for score in ["UCS", "EIS", "MDI", "mech"]:
        groups = {g: sc.loc[sc["group"] == g, score].values for g in order}
        beta, se, t, dfree, p_two = welch_contrast(groups)
        obs, p_perm_two, p_perm_pos = perm_p(sc[score].values, sc["group"].values)
        d = DIRECTION[score]
        p_one = one_sided(p_two, beta, d)
        # 置换层方向性：MDI/EIS/UCS 用方向化经验 p；mech 双向
        if d == 1:
            p_perm_dir = p_perm_pos
        elif d == -1:
            p_perm_dir = 1 - p_perm_pos + 1 / (10000 + 1)
        else:
            p_perm_dir = p_perm_two
        rows.append({"gse": "GSE2411", "score": score, "beta3_interaction": beta, "se": se,
                     "welch_t": t, "df_satterthwaite": dfree, "welch_p_two": p_two,
                     "welch_p_one": p_one, "perm_p_two": p_perm_two, "perm_p_directional": p_perm_dir,
                     "direction_hypothesis": {-1: "decrease", 1: "increase", 0: "two-sided"}[d],
                     "mean_Control": groups["Control"].mean(), "mean_MV": groups["MV"].mean(),
                     "mean_LPS": groups["LPS"].mean(), "mean_MVLPS": groups["MV+LPS"].mean(),
                     "n_per_group": 6, "perm_B": 10000, "seed": 20260831})
        note(f"[{score}] β₃={beta:+.4f} (SE {se:.4f}, t={t:+.3f}, df={dfree:.1f}) "
             f"welch_p_one={p_one:.4f} perm_dir={p_perm_dir:.4f}")
    res = pd.DataFrame(rows)
    res = res.sort_values("welch_p_one")
    m = len(res)
    res["bh_rank"] = np.arange(1, m + 1)
    res["welch_p_one_bh"] = np.minimum.accumulate((res["welch_p_one"] * m / res["bh_rank"])[::-1])[::-1].clip(upper=1)
    res["perm_p_directional_bh"] = np.minimum.accumulate(
        (res.sort_values("perm_p_directional")["perm_p_directional"] * m /
         np.arange(1, m + 1))[::-1])[::-1].clip(upper=1).values

    f90 = os.path.join(TAB, "Table_S90_M16_GSE2411_Interaction.csv")
    res.to_csv(f90, index=False, encoding="utf-8-sig")
    note(f"[write] {f90}")

    # 基因级互作（描述性）：从 step2 的中间产物重建基因级 z 矩阵
    sys.path.insert(0, ROOT)
    import importlib
    s2 = importlib.import_module("M16_step2_vili_scoring")
    mat, titles = s2.load_series("GSE2411")
    grp = s2.groups("GSE2411", titles)
    pmap = s2.load_gpl_map("GPL339")
    sym = pd.Series([pmap.get(p, "") for p in mat.index], index=mat.index)
    mat = mat.loc[sym != ""]; sym = sym[sym != ""]
    df = mat.copy(); df["sym"] = sym
    gm = df.groupby("sym").mean()
    man = pd.read_csv(os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "Mitoxyperilysis_Gene_Manifest_v1.0.csv"), encoding="utf-8-sig")
    orth = pd.read_csv(os.path.join(INTER, "M16_human2mouse_orthologs.csv"), encoding="utf-8-sig").fillna("")
    h2m = {r["human_symbol"].upper(): r["mouse_symbol"] for _, r in orth.iterrows() if r["mouse_symbol"]}
    m2h = {v: k for k, v in h2m.items()}
    gm_h = gm.loc[[s for s in gm.index if s in m2h]]
    gm_h.index = [m2h[s] for s in gm_h.index]
    geneset = sorted(set(gm_h.index.str.upper()))
    Z = (gm_h - gm_h.mean(axis=1).values[:, None]) / gm_h.std(axis=1, ddof=1).values[:, None]
    grows = []
    lab = pd.Series(grp)
    for gene in Z.index:
        groups = {g: Z.loc[gene, lab[lab == g].index].values for g in order}
        if min(len(v) for v in groups.values()) < 2:
            continue
        beta, se, t, dfree, p_two = welch_contrast(groups)
        grows.append({"human_ortholog": gene, "beta3": beta, "welch_t": t, "welch_p_two": p_two})
    gd = pd.DataFrame(grows).sort_values("welch_p_two")
    fg = os.path.join(INTER, "M16_gse2411_genelevel_interaction.csv")
    gd.to_csv(fg, index=False, encoding="utf-8-sig")
    note(f"[write] {fg}（{len(gd)} 基因，描述性层，未校正）")
    top = gd.head(10)
    for _, r in top.iterrows():
        note(f"  top互作: {r['human_ortholog']} β₃={r['beta3']:+.3f} p={r['welch_p_two']:.4f}")

    with open(LOG, "w", encoding="utf-8") as f:
        f.write("\n".join(log) + "\n")
    print("[done]")

if __name__ == "__main__":
    main()
