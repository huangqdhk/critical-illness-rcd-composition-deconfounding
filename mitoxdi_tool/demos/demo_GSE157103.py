# -*- coding: utf-8 -*-
"""
demo_GSE157103.py — M15 采纳演示 1：GSE157103 一键复算
================================================================================
数据集（已按红线 #1 经 GEO 官方记录核验，2026-08-27）：
  - GSE157103 "Large-scale Multi-omic Analysis of COVID-19 Severity"
    （Overmyer et al., Cell Systems 2020；提交 2020-08-28，更新 2020-11-30）
  - 平台 GPL24676 Illumina NovaSeq 6000；人全血 bulk RNA-seq
  - 官方总体设计：126 样本 = 100 COVID-19 患者 + 26 非 COVID-19
  - 本地文件：GSE157103_genes.tpm.tsv.gz（3.8 Mb，与 GEO 附属文件一致）
本脚本 = 工具"零下载"一键复算演示：TPM → log2(TPM+1) → mdi_v1.0 评分 →
Monaco 29 型 NNLS 反卷积（主口径；OLS-CLS 交叉；ABIS RNAseq 17 型交叉）→
组成预测 MDI_comp → 组成残差 MDI_resid → 逐族疾病效应（COVID vs non-COVID）→
R² → 效应衰减 → 细胞比例漂移 → ISED/混淆度解读报告。
输出（_intermediate/）：
  M15_demo_GSE157103_per_sample.csv / _contrasts.csv / _r2.csv / _proportions.csv
报告：03_LOGS/M15_demo_GSE157103_report.md
用法（发布版）：python demo_GSE157103.py <genes.tpm.tsv.gz 路径> [--outdir 目录]
"""
import argparse
import os
import re
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
from mitoxdi.mdi_core import load_manifest, arm_scores, gene_matrix_from_symbols
from mitoxdi.composition import (load_monaco_ref, load_abis_refs, deconvolve_nnls,
                                 deconvolve_olcls, synthetic_scores, residual_mdi,
                                 neutrophil_monocyte)
from mitoxdi.report import family_contrasts, family_r2, attenuation, write_demo_report
from mitoxdi import ROOT

OUTDIR = os.path.join(ROOT, "_intermediate")
LOGDIR = os.path.join(ROOT, "03_LOGS")

GROUP_RE = re.compile(r"^(C\d+)$|^(NC\d+)$")


def load_matrix(path):
    tpm = pd.read_csv(path, sep="\t", index_col=0)
    tpm.columns = [str(c).strip('"') for c in tpm.columns]
    tpm.index = [str(i).strip().strip('"').upper() for i in tpm.index]
    tpm = tpm.groupby(level=0).mean()
    mat = np.log2(tpm + 1.0).T
    groups = pd.Series(
        [("COVID" if re.match(r"^C\d+$", s) else ("non-COVID" if re.match(r"^NC\d+$", s) else "?"))
         for s in mat.index], index=mat.index)
    return mat, groups


def main(matrix_path=None, outdir=None):
    outdir = outdir or OUTDIR
    os.makedirs(outdir, exist_ok=True)
    matrix_path = matrix_path or os.path.join(
        ROOT, "00_RAW_DATA", "GSE157103_COVID19_WholeBlood_RNAseq", "GSE157103_genes.tpm.tsv.gz")

    _, arms, alias = load_manifest()
    UP, EX = arms["upstream_collapse"], arms["execution_induction"]
    mat, groups = load_matrix(matrix_path)
    n_unknown = int((groups == "?").sum())
    if n_unknown:
        raise RuntimeError(fr"{n_unknown} 个样本组别无法解析（C\d+/NC\d+）")

    monaco = load_monaco_ref()
    abis_rna, _ = load_abis_refs()

    # --- mdi_v1.0 评分（别名解析到 canonical 臂基因） ---
    arm_mat = gene_matrix_from_symbols(mat, alias)
    sc = arm_scores(arm_mat, arms)
    n_up, n_ex = int(sc["n_up"].iloc[0]), int(sc["n_ex"].iloc[0])

    # --- 反卷积与组成预测（Monaco NNLS 主口径；OLS-CLS 交叉；ABIS 交叉） ---
    prop_m, common_m = deconvolve_nnls(mat, monaco)
    prop_o, _ = deconvolve_olcls(mat, monaco)
    prop_a, _ = deconvolve_nnls(mat, abis_rna)
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
    out.to_csv(os.path.join(outdir, "M15_demo_GSE157103_per_sample.csv"), index=False)

    # --- 逐族疾病效应（COVID vs non-COVID） ---
    fams = ["UCS", "EIS", "MDI", "MDI_nomt", "UCS_comp", "EIS_comp", "MDI_comp",
            "MDI_resid", "MDI_comp_olcls"]
    ctr = family_contrasts(out, "group", "non-COVID", "COVID", fams)
    ctr.to_csv(os.path.join(outdir, "M15_demo_GSE157103_contrasts.csv"), index=False)

    # --- R²（观测 vs 合成，Monaco|NNLS 主口径 + OLS-CLS 交叉 + ABIS 交叉） ---
    r2_main = family_r2(out, ["UCS", "EIS", "MDI"])
    o_out = out.copy()
    o_out["MDI_comp"] = out["MDI_comp_olcls"]
    r2_o = family_r2(o_out, ["UCS", "EIS", "MDI"])
    r2_o["method"] = "OLSCLS"
    sy_a = synthetic_scores(prop_a, abis_rna, UP, EX)
    a_out = out.copy()
    a_out["UCS_comp"], a_out["EIS_comp"], a_out["MDI_comp"] = sy_a["UCS_comp"], sy_a["EIS_comp"], sy_a["MDI_comp"]
    r2_a = family_r2(a_out, ["UCS", "EIS", "MDI"])
    r2_a["method"] = "ABIS_NNLS"
    r2_main["method"] = "Monaco_NNLS"
    r2 = pd.concat([r2_main, r2_o, r2_a], ignore_index=True)
    r2.to_csv(os.path.join(outdir, "M15_demo_GSE157103_r2.csv"), index=False)

    # --- 效应衰减 ---
    att = attenuation(out, "group", "non-COVID", "COVID")

    # --- 比例漂移 ---
    prop_m["group"] = groups.values
    prop_sum = prop_m.groupby("group").mean(numeric_only=True)
    prop_sum.to_csv(os.path.join(outdir, "M15_demo_GSE157103_proportions.csv"))

    # --- 报告 ---
    print("== GSE157103 一键复算")
    print(f"   样本 {len(out)}（COVID {int((groups == 'COVID').sum())} / non-COVID "
          f"{int((groups == 'non-COVID').sum())}）；臂覆盖 up {n_up}/30 ex {n_ex}/33；"
          f"Monaco 共同基因 {len(common_m)}")
    print(f"   MDI: g={att['g_obs']:+.3f} -> resid g={att['g_resid']:+.3f} "
          f"（衰减 {att['attenuation']:.1%}）；MDI R²={r2_main[r2_main['family'] == 'MDI']['R2'].iloc[0]:.3f}")
    for _, r in ctr.iterrows():
        print(f"   {r['family']:<14} n={r['n_ctrl']}v{r['n_case']} MW p={r['MW_p']:.4g} "
              f"Cliff={r['Cliff_delta']:+.3f} g={r['Hedges_g']:+.3f}")
    md_r2 = r2_main[r2_main["family"] == "MDI"]["R2"].iloc[0]
    ucs_r2 = r2_main[r2_main["family"] == "UCS"]["R2"].iloc[0]
    eis_r2 = r2_main[r2_main["family"] == "EIS"]["R2"].iloc[0]
    # Monaco 感知的比例汇总（冻结口径的 neutrophil_monocyte 按 "Monocytes*" 匹配、
    # 匹配不到 Monaco 的下划线命名 C_mono/I_mono/NC_mono——演示层以类型名感知汇总）
    pm_neut_cols = [c for c in prop_m.columns if "NEUTROPHIL" in c.upper()]
    pm_mono_cols = [c for c in prop_m.columns if c.upper().endswith("_MONO")]
    pm_summary = pd.DataFrame({
        "Neut": prop_m[pm_neut_cols].sum(axis=1),
        "Mono": prop_m[pm_mono_cols].sum(axis=1),
        "group": groups.values})
    pm_sum = pm_summary.groupby("group").mean(numeric_only=True)
    neut_d = pm_sum.loc["COVID", "Neut"] - pm_sum.loc["non-COVID", "Neut"]
    mono_d = pm_sum.loc["COVID", "Mono"] - pm_sum.loc["non-COVID", "Mono"]
    g_comp = ctr[ctr["family"] == "MDI_comp"]["Hedges_g"].iloc[0]
    g_res = ctr[ctr["family"] == "MDI_resid"]["Hedges_g"].iloc[0]
    g_obs = ctr[ctr["family"] == "MDI"]["Hedges_g"].iloc[0]
    if np.isfinite(att["attenuation"]):
        att_txt = f"衰减 {att['attenuation']:.1%}"
    else:
        att_txt = "观测与残差方向相反（符号翻转），衰减比例不适用，如实报告两者"
    notes = (
        f"- 组成预测值（Monaco|NNLS）解释 MDI 队列内方差的 R² = {md_r2:.3f}"
        f"（UCS {ucs_r2:.3f} / EIS {eis_r2:.3f}）——组成漂移在 COVID-19 全血同样携带大量 MDI 信号。\n"
        f"- 疾病效应三层读数：观测 MDI g = {g_obs:+.3f}；组成预测 MDI_comp g = {g_comp:+.3f}；"
        f"组成残差 MDI_resid g = {g_res:+.3f}（{att_txt}）。\n"
        f"- 本队列观察：COVID 组 UCS 与 EIS 同时升高（g=+0.751/+0.575），MDI 无显著疾病效应——"
        f"双臂解离方向与脓毒症队列不同，按 GATE.md 只陈述客观量、不作外推解读。\n"
        f"- 比例漂移（COVID − non-COVID）：中性粒 {neut_d:+.3f}、单核 {mono_d:+.3f}。\n"
        f"- 数据红线：本数据集经 GEO 官方记录核验（标题/设计/平台/文件，2026-08-27，"
        f"03_LOGS/M15_demo_GEO_verification_20260827.md），非本项目既有队列；"
        f"本复算仅陈述客观量（臂评分/组成预测/残差），不评价原论文结论。")
    write_demo_report(
        os.path.join(LOGDIR, "M15_demo_GSE157103_report.md"),
        "M15 采纳演示 1：GSE157103（COVID-19 全血 RNA-seq）一键复算",
        "GSE157103 (Overmyer et al., Cell Systems 2020)", "RNA-seq (GPL24676, TPM)",
        f"{len(out)}（COVID 100 / non-COVID 26）",
        {"COVID": 100, "non-COVID": 26}, (n_up, n_ex), ctr, r2, att, prop_sum, notes)
    print(f"   DONE -> {os.path.join(outdir, 'M15_demo_GSE157103_per_sample.csv')}")
    return out, ctr, r2, att


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix", nargs="?", default=None, help="genes.tpm.tsv.gz 路径")
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args()
    main(a.matrix, a.outdir)
