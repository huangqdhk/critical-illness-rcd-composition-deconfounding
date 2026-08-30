# -*- coding: utf-8 -*-
"""
demo_GSE66099.py — M15 采纳演示 2：GSE66099 一键复算
================================================================================
数据集（已按红线 #1 经 GEO 官方记录核验，2026-08-27）：
  - GSE66099 "Unique Patients from the Genomics of Pediatric SIRS and Septic Shock
    Investigators (GPSSSI)"；提交 2015-02-19，更新 2020-11-30
  - 平台 GPL570 [HG-U133_Plus_2] Affymetrix；人全血微阵列（ICU Day 1）
  - 官方总体设计：276 个唯一患者（Day 1）= 6 个早期系列（GSE4607/GSE8121/GSE9692/
    GSE13904/GSE26378/GSE26440）合并；.CEL 文件以 gcRMA 统一重归一化
  - 本地文件：GSE66099_series_matrix.txt.gz（GEO 系列矩阵）
本脚本 = 工具"零下载"一键复算演示：系列矩阵 → GPL570 探针全塌陷（max-mean ≤3）→
mdi_v1.0 评分 → Monaco 29 型 NNLS 反卷积（主口径；OLS-CLS 交叉；ABIS Micro 11 型
平台匹配交叉）→ MDI_comp → MDI_resid → 逐族疾病效应（SepticShock/Sepsis/SIRS vs
Control）→ R² → 效应衰减 → 细胞比例漂移 → ISED/混淆度解读报告。
输出（_intermediate/）：
  M15_demo_GSE66099_per_sample.csv / _contrasts.csv / _r2.csv / _proportions.csv
报告：03_LOGS/M15_demo_GSE66099_report.md
用法（发布版）：python demo_GSE66099.py <series_matrix.txt.gz 路径> [--outdir 目录]
"""
import argparse
import os
import sys
import warnings

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.dirname(HERE)
sys.path.insert(0, TOOL)
sys.path.insert(0, os.path.join(TOOL, ".."))

import numpy as np
import pandas as pd

from mitoxdi import loaders as LD
from mitoxdi.mdi_core import load_manifest, arm_scores, parse_series_matrix, sample_characteristics
from mitoxdi.composition import (load_monaco_ref, load_abis_refs, deconvolve_nnls,
                                 deconvolve_olcls, synthetic_scores, residual_mdi,
                                 neutrophil_monocyte)
from mitoxdi.report import family_contrasts, family_r2, attenuation, write_demo_report
from mitoxdi import ROOT

OUTDIR = os.path.join(ROOT, "_intermediate")
LOGDIR = os.path.join(ROOT, "03_LOGS")

CONTRASTS = [("Control", "SepticShock"), ("Control", "Sepsis"), ("Control", "SIRS")]


def load_matrix(path):
    meta, samples, mat = parse_series_matrix(path)
    p2s = LD.read_probe2sym_csv(LD.GPL570_CSV)
    g = LD.collapse_all(mat, p2s)
    disease = sample_characteristics(meta, label="disease")
    groups = pd.Series(disease, index=samples)
    return g, groups


def main(matrix_path=None, outdir=None):
    outdir = outdir or OUTDIR
    os.makedirs(outdir, exist_ok=True)
    matrix_path = matrix_path or os.path.join(
        ROOT, "00_RAW_DATA", "GSE66099_Pediatric_Sepsis_WholeBlood",
        "GSE66099_series_matrix.txt.gz")

    _, arms, _ = load_manifest()
    UP, EX = arms["upstream_collapse"], arms["execution_induction"]
    mat, groups = load_matrix(matrix_path)
    known = groups.isin(["Control", "SIRS", "Sepsis", "SepticShock"])
    if not known.all():
        raise RuntimeError(f"未知组别：{sorted(set(groups[~known]))}")

    monaco = load_monaco_ref()
    _, abis_micro = load_abis_refs()

    # --- mdi_v1.0 评分（GPL570 直接符号塌陷矩阵，别名口径同 M2 冻结路径） ---
    sc = arm_scores(mat, arms)
    n_up, n_ex = int(sc["n_up"].iloc[0]), int(sc["n_ex"].iloc[0])

    # --- 反卷积与组成预测 ---
    prop_m, common_m = deconvolve_nnls(mat, monaco)
    prop_o, _ = deconvolve_olcls(mat, monaco)
    prop_a, _ = deconvolve_nnls(mat, abis_micro)
    sy_m = synthetic_scores(prop_m, monaco, UP, EX)
    sy_o = synthetic_scores(prop_o, monaco, UP, EX)
    resid, b0, b1 = residual_mdi(sc["MDI"].values, sy_m["MDI_comp"].values)
    out = pd.concat([sc, sy_m], axis=1)
    out["MDI_comp_olcls"] = sy_o["MDI_comp"].values
    out["MDI_resid"] = resid
    out["b0"], out["b1"] = b0, b1
    nm = neutrophil_monocyte(prop_m)
    out["Neut"], out["Mono"], out["comp"] = nm["Neut"], nm["Mono"], nm["comp"]
    out["group"] = groups.values
    out["sample"] = [str(i) for i in mat.index]
    out.to_csv(os.path.join(outdir, "M15_demo_GSE66099_per_sample.csv"), index=False)

    # --- 逐族疾病效应（三组 vs Control） ---
    fams = ["UCS", "EIS", "MDI", "MDI_nomt", "UCS_comp", "EIS_comp", "MDI_comp", "MDI_resid"]
    ctr_parts = []
    for ctrl, case in CONTRASTS:
        c = family_contrasts(out, "group", ctrl, case, fams)
        c["contrast"] = f"{case}_vs_{ctrl}"
        ctr_parts.append(c)
    ctr = pd.concat(ctr_parts, ignore_index=True)
    ctr.to_csv(os.path.join(outdir, "M15_demo_GSE66099_contrasts.csv"), index=False)

    # --- R²（观测 vs 合成；主对比口径 = 全样本） ---
    r2_main = family_r2(out, ["UCS", "EIS", "MDI"])
    o_out = out.copy()
    o_out["MDI_comp"] = out["MDI_comp_olcls"]
    r2_o = family_r2(o_out, ["UCS", "EIS", "MDI"])
    r2_o["method"] = "OLSCLS"
    sy_a = synthetic_scores(prop_a, abis_micro, UP, EX)
    a_out = out.copy()
    a_out["UCS_comp"], a_out["EIS_comp"], a_out["MDI_comp"] = sy_a["UCS_comp"], sy_a["EIS_comp"], sy_a["MDI_comp"]
    r2_a = family_r2(a_out, ["UCS", "EIS", "MDI"])
    r2_a["method"] = "ABIS_Micro_NNLS"
    r2_main["method"] = "Monaco_NNLS"
    r2 = pd.concat([r2_main, r2_o, r2_a], ignore_index=True)
    r2.to_csv(os.path.join(outdir, "M15_demo_GSE66099_r2.csv"), index=False)

    # --- 效应衰减（主对比 SepticShock vs Control） ---
    att = attenuation(out, "group", "Control", "SepticShock")

    # --- 比例漂移 ---
    prop_m["group"] = groups.values
    prop_sum = prop_m.groupby("group").mean(numeric_only=True)
    prop_sum.to_csv(os.path.join(outdir, "M15_demo_GSE66099_proportions.csv"))

    # --- 报告 ---
    print("== GSE66099 一键复算")
    print(f"   样本 {len(out)}；组别 {groups.value_counts().to_dict()}；"
          f"臂覆盖 up {n_up}/30 ex {n_ex}/33；Monaco 共同基因 {len(common_m)}")
    print(f"   主对比 SepticShock vs Control：MDI g={att['g_obs']:+.3f} -> resid g={att['g_resid']:+.3f} "
          f"（衰减 {att['attenuation']:.1%}）；MDI R²={r2_main[r2_main['family'] == 'MDI']['R2'].iloc[0]:.3f}")
    for _, r in ctr.iterrows():
        print(f"   {r['contrast']:<24} {r['family']:<12} n={r['n_ctrl']}v{r['n_case']} "
              f"MW p={r['MW_p']:.4g} Cliff={r['Cliff_delta']:+.3f} g={r['Hedges_g']:+.3f}")
    md_r2 = r2_main[r2_main["family"] == "MDI"]["R2"].iloc[0]
    ucs_r2 = r2_main[r2_main["family"] == "UCS"]["R2"].iloc[0]
    eis_r2 = r2_main[r2_main["family"] == "EIS"]["R2"].iloc[0]
    # Monaco 感知的比例汇总（同 demo_GSE157103 注释）
    pm_neut_cols = [c for c in prop_m.columns if "NEUTROPHIL" in c.upper()]
    pm_mono_cols = [c for c in prop_m.columns if c.upper().endswith("_MONO")]
    pm_summary = pd.DataFrame({
        "Neut": prop_m[pm_neut_cols].sum(axis=1),
        "Mono": prop_m[pm_mono_cols].sum(axis=1),
        "group": groups.values})
    pm_sum = pm_summary.groupby("group").mean(numeric_only=True)
    neut_d = pm_sum.loc["SepticShock", "Neut"] - pm_sum.loc["Control", "Neut"]
    mono_d = pm_sum.loc["SepticShock", "Mono"] - pm_sum.loc["Control", "Mono"]
    g_comp = ctr[(ctr["contrast"] == "SepticShock_vs_Control") & (ctr["family"] == "MDI_comp")]["Hedges_g"].iloc[0]
    g_res = ctr[(ctr["contrast"] == "SepticShock_vs_Control") & (ctr["family"] == "MDI_resid")]["Hedges_g"].iloc[0]
    g_obs = ctr[(ctr["contrast"] == "SepticShock_vs_Control") & (ctr["family"] == "MDI")]["Hedges_g"].iloc[0]
    notes = (
        f"- 组成预测值（Monaco|NNLS）解释 MDI 队列内方差的 R² = {md_r2:.3f}"
        f"（UCS {ucs_r2:.3f} / EIS {eis_r2:.3f}）——儿童脓毒症全血中组成漂移携带的 MDI 信号量级。\n"
        f"- 疾病效应三层读数（SepticShock vs Control）：观测 MDI g = {g_obs:+.3f}；"
        f"组成预测 MDI_comp g = {g_comp:+.3f}；组成残差 MDI_resid g = {g_res:+.3f}"
        f"（衰减 {att['attenuation']:.1%}）。\n"
        f"- 比例漂移（SepticShock − Control）：中性粒 {neut_d:+.3f}、单核 {mono_d:+.3f}。\n"
        f"- 解读（依据 GATE.md）：按 M10 检验 A 判定门语言解读三套读数"
        f"（未校正/组成校正/组成残差），如实报告 R² 与衰减，不评价格式化结论。\n"
        f"- 数据红线：本数据集经 GEO 官方记录核验（标题/设计/平台/文件，2026-08-27，"
        f"03_LOGS/M15_demo_GEO_verification_20260827.md），非本项目既有队列；"
        f"本复算仅陈述客观量（臂评分/组成预测/残差），不评价原论文结论。")
    write_demo_report(
        os.path.join(LOGDIR, "M15_demo_GSE66099_report.md"),
        "M15 采纳演示 2：GSE66099（儿童 SIRS/脓毒症/脓毒性休克全血）一键复算",
        "GSE66099 (GPSSSI, Wong et al.)", "microarray (GPL570, gcRMA)",
        f"{len(out)}（SepticShock 181 / Sepsis 18 / SIRS 30 / Control 47）",
        groups.value_counts().to_dict(), (n_up, n_ex), ctr, r2, att, prop_sum, notes)
    print(f"   DONE -> {os.path.join(outdir, 'M15_demo_GSE66099_per_sample.csv')}")
    return out, ctr, r2, att


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix", nargs="?", default=None, help="series_matrix.txt.gz 路径")
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args()
    main(a.matrix, a.outdir)
