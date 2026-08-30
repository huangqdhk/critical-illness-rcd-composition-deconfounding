# -*- coding: utf-8 -*-
"""
P1_metabolite_anchor.py — Phase 1 任务3：ST002738 代谢组锚定（描述性）
========================================================================
依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md §4 H6——"代谢锚为方向一致性检查，描述性非确证"
设计：2×2（Media/LPS/CS/LPS+CS 各 n=3；Metabo_00=BLANK 弃），分组映射来自 mwTab
     SUBJECT_SAMPLE_FACTORS（Treatment 字段）。
鉴定：公共库未沉积命名注释（mwTab 无 Metabolites 节）→ 按提交方自报容差
     m/z ±5 ppm 做目标代谢物精确质量匹配（[M+H]+ 与 [M−H]−）。
统计：目标特征的四组均值 + 交互对照（Δ_inter = LC − L − C + M，log2 面积）。
输出：_intermediate/P1_metabolite_anchor.csv；结果追加至 P1_Interaction_Signature_Report.md
"""
import os, sys
import numpy as np
import pandas as pd
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
RAW = os.path.join(ROOT, "00_RAW_DATA", "ST002738")
OUT = os.path.join(ROOT, "_intermediate", "P1_metabolite_anchor.csv")
REPORT = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "P1_Interaction_Signature_Report.md")

PROTON = 1.00727646
TARGETS = {  # 名称: (中性单同位素质量)
    "GSH": 307.083806, "GSSG": 612.156395,
    "glutamate": 147.053158, "5-oxoproline": 129.042593,
    "cysteine": 121.019749, "cystine": 240.023848,
    "gamma-glutamylcysteine": 251.068716, "NADP+": 744.083248,
    "NADPH": 745.091103, "glycine": 75.031694,
}
GROUPS = {"01": "Media", "02": "Media", "03": "Media",
          "04": "LPS", "05": "LPS", "06": "LPS",
          "07": "CS", "08": "CS", "09": "CS",
          "10": "LPS_CS", "11": "LPS_CS", "12": "LPS_CS"}

rows = []
for ana, mode in [("AN004440", "pos"), ("AN004441", "neg")]:
    p = os.path.join(RAW, f"ST002738_{ana}_Results.txt")
    df = pd.read_csv(p, sep="\t")
    df.columns = [str(c).strip() for c in df.columns]
    df["mz"] = df["Name"].str.split("_").str[0].astype(float)
    df["rt"] = df["Name"].str.split("_").str[1].astype(float)
    # 样本列 → 组
    samp_cols = {}
    for c in df.columns:
        for num, g in GROUPS.items():
            if f"Metabo_{num}_" in c:
                samp_cols[c] = g
    for name, mass in TARGETS.items():
        target_mz = mass + PROTON if mode == "pos" else mass - PROTON
        tol = target_mz * 5e-6
        hits = df[(df["mz"] > target_mz - tol) & (df["mz"] < target_mz + tol)]
        for _, h in hits.iterrows():
            vals = {g: h[[c for c, gg in samp_cols.items() if gg == g]].astype(float).mean()
                    for g in ("Media", "LPS", "CS", "LPS_CS")}
            allv = np.array([vals[g] for g in ("Media", "LPS", "CS", "LPS_CS")], float)
            with np.errstate(divide="ignore", invalid="ignore"):
                logv = np.log2(allv)
            inter = logv[3] - logv[1] - logv[2] + logv[0] if (allv > 0).all() else np.nan
            rows.append(dict(metabolite=name, mode=mode, mz=h["mz"], rt=h["rt"], **vals,
                             log2_interaction=inter))

out = pd.DataFrame(rows)
out.to_csv(OUT, index=False)
print(out.to_string(index=False) if len(out) else "无匹配特征")

with open(REPORT, "a", encoding="utf-8") as f:
    f.write("\n## P1-3 代谢组锚定（描述性，追加于 2026-08-20）\n\n")
    f.write("- 设计：2×2（Media/LPS/CS/LPS+CS 各 n=3），分组映射自 mwTab SUBJECT_SAMPLE_FACTORS；\n")
    f.write("- 鉴定：公共库未沉积命名注释，按提交方自报 m/z ±5 ppm 容差做目标代谢物精确质量匹配（[M+H]+/[M−H]−）；\n")
    f.write("- 性质：**描述性锚定，非确证检验**（冻结计划 §4 H6）。\n\n")
    if len(out):
        f.write("| 代谢物 | 模式 | m/z | RT(min) | Media | LPS | CS | LPS+CS | log2交互 |\n|---|---|---|---|---|---|---|---|---|\n")
        for _, r in out.iterrows():
            f.write(f"| {r['metabolite']} | {r['mode']} | {r['mz']:.4f} | {r['rt']:.2f} | "
                    f"{r['Media']:.3g} | {r['LPS']:.3g} | {r['CS']:.3g} | {r['LPS_CS']:.3g} | "
                    f"{r['log2_interaction']:+.2f} |\n")
        f.write("\n读法：log2交互>0 表示联合刺激超出两单刺激之和（加性以上）；GSH/GSSG/谷胱甘肽前体的方向"
                "与原文『GSH 耗竭』一致性由读者按原文图核对；±5 ppm 单一同位素匹配存在假阳性可能，"
                "无 RT 注释佐证，故只作方向参考。\n")
    else:
        f.write("（±5 ppm 窗口内未匹配到目标代谢物特征——结果如实登记为空，不外推。）\n")
print("\nDONE P1-3")
