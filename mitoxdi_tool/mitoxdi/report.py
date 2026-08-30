# -*- coding: utf-8 -*-
"""
report.py — 对比统计、R² 对比、比例汇总与 ISED/混淆度解读报告（M15 工具核心之三）
================================================================================
- 对比口径：Mann-Whitney U 双侧 p、Cliff's delta（b 相对 a；2,000 bootstrap 95%CI）、
  Hedges' g（b 相对 a）± SE —— 与 mdi_lib 冻结实现一致。
- R²：队列内 obs vs comp 的 Pearson r²（UCS/EIS/MDI 三族）。
- ISED 解读（混淆度解读）输出：组成可解释方差占比、疾病效应衰减比例
  （g_obs -> g_resid）、组间中性粒/单核比例漂移。
"""
import numpy as np
import pandas as pd

from .mdi_core import mwu_cliff, hedges_g, clf_delta_ci


def contrast_stats(ctrl, case):
    """返回 dict：n_ctrl/n_case、MW 双侧 p、Cliff delta(±95%CI boot)、Hedges g±SE。"""
    ctrl = np.asarray(ctrl, float)
    case = np.asarray(case, float)
    p, delta = mwu_cliff(ctrl, case)
    delta_lo = delta_hi = np.nan
    try:
        _, delta_lo, delta_hi = clf_delta_ci(ctrl, case, n_boot=2000, seed=0)
    except Exception:
        pass
    g, gse = hedges_g(ctrl, case)
    return dict(n_ctrl=len(ctrl), n_case=len(case), MW_p=p, Cliff_delta=delta,
                Cliff_lo=delta_lo, Cliff_hi=delta_hi, Hedges_g=g, g_SE=gse)


def family_contrasts(df, group_col, ctrl_label, case_label, families):
    """对多族指标做同口径对比。返回 DataFrame（行=family）。"""
    rows = []
    g = df[group_col].astype(str)
    ctrl = df[g == ctrl_label]
    case = df[g == case_label]
    for fam in families:
        if fam not in df.columns:
            continue
        r = contrast_stats(ctrl[fam].dropna().values, case[fam].dropna().values)
        r["family"] = fam
        rows.append(r)
    return pd.DataFrame(rows)[["family", "n_ctrl", "n_case", "MW_p", "Cliff_delta",
                               "Cliff_lo", "Cliff_hi", "Hedges_g", "g_SE"]]


def r2_obs_comp(obs, comp):
    """obs/comp 对齐向量 -> (r, R2, n)。"""
    r = np.corrcoef(obs, comp)[0, 1]
    return float(r), float(r ** 2), int(len(obs))


def family_r2(df, families):
    """df 含 family 与 family_comp 列。返回 DataFrame（行=family）。"""
    rows = []
    for fam in families:
        if fam not in df.columns or fam + "_comp" not in df.columns:
            continue
        sub = df[[fam, fam + "_comp"]].dropna()
        r, r2, n = r2_obs_comp(sub[fam].values, sub[fam + "_comp"].values)
        rows.append(dict(family=fam, n=n, r=r, R2=r2))
    return pd.DataFrame(rows)


def attenuation(df, group_col, ctrl_label, case_label):
    """疾病效应衰减：obs vs resid 的 Hedges g 与衰减比例 = 1 - g_resid/g_obs。
    观测与残差方向相反（符号翻转）时衰减比例不适用，返回 nan。"""
    g = df[group_col].astype(str)
    ctrl = df[g == ctrl_label]
    case = df[g == case_label]
    g_obs, _ = hedges_g(ctrl["MDI"].values, case["MDI"].values)
    g_res, _ = hedges_g(ctrl["MDI_resid"].values, case["MDI_resid"].values)
    if abs(g_obs) > 1e-9 and np.sign(g_obs) == np.sign(g_res):
        att = 1.0 - g_res / g_obs
    else:
        att = np.nan  # 方向翻转：衰减比例不适用
    return dict(g_obs=g_obs, g_resid=g_res, attenuation=att)


def proportion_summary(prop, groups, ctrl_label, case_label):
    """组间主要细胞型比例均值（含 Neut/Mono/comp 合计）。"""
    out = pd.DataFrame(index=prop.index)
    neut_cols = [c for c in prop.columns if "neutrophil" in c.lower()]
    mono_cols = [c for c in prop.columns if c.lower() in ("monocytes", "monocytes c", "monocytes nc+i")]
    out["Neut"] = prop[neut_cols].sum(axis=1) if neut_cols else np.nan
    out["Mono"] = prop[mono_cols].sum(axis=1) if mono_cols else np.nan
    out["comp"] = out["Neut"] + out["Mono"]
    out["group"] = np.asarray(groups)
    return out.groupby("group").mean(numeric_only=True)


def write_demo_report(path, title, dataset, platform, n_samples, group_counts,
                      coverage, contrasts, r2_table, att, prop_sum, notes=""):
    """写出采纳演示的 Markdown 报告（ISED/混淆度解读）。"""
    lines = []
    lines.append(f"# {title}\n")
    lines.append(f"- 数据集：{dataset}；平台：{platform}；样本：{n_samples}")
    lines.append(f"- 组别：{group_counts}")
    lines.append(f"- 臂覆盖（mdi_v1.0）：上游 {coverage[0]}/30，执行 {coverage[1]}/33\n")
    lines.append("## 1. 逐族疾病效应（观测 / 组成预测 / 组成残差）\n\n")
    lines.append("| 指标 | n(对照) | n(病例) | MW p | Cliff δ | Hedges g (SE) |\n|---|---|---|---|---|---|\n")
    for _, r in contrasts.iterrows():
        lines.append(f"| {r['family']} | {r['n_ctrl']} | {r['n_case']} | {r['MW_p']:.4g} "
                     f"| {r['Cliff_delta']:+.3f} | {r['Hedges_g']:+.3f} ({r['g_SE']:.3f}) |\n")
    lines.append("\n## 2. 组成预测 R²（观测 vs 合成谱，Monaco|NNLS 主口径）\n\n")
    lines.append("| 指标 | n | r | R² |\n|---|---|---|---|\n")
    for _, r in r2_table.iterrows():
        lines.append(f"| {r['family']} | {r['n']} | {r['r']:+.3f} | {r['R2']:.3f} |\n")
    lines.append("\n## 3. 疾病效应衰减（g_obs -> g_resid）\n\n")
    att_txt = f"{att['attenuation']:.1%}" if np.isfinite(att["attenuation"]) else "不适用（方向翻转）"
    lines.append(f"- MDI 观测效应 g = {att['g_obs']:+.3f}；组成残差效应 g = {att['g_resid']:+.3f}；"
                 f"衰减比例 = {att_txt}")
    lines.append("\n## 4. 组间细胞比例（主细胞型均值）\n\n")
    lines.append(prop_sum.to_string())
    lines.append("\n## 5. ISED/混淆度解读（依据 GATE.md 判定门语言）\n")
    lines.append(notes)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
