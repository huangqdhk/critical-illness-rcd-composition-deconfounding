# -*- coding: utf-8 -*-
"""mmc3 原始矩阵：补齐 mmc4 缺失的关键铁轴/氧化应激蛋白的肺组织 C/N 方向。"""
import pandas as pd
import io, sys, os, re
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
import numpy as np
from scipy.stats import mannwhitneyu

DL = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\00_RAW_DATA\Nie2021_Cell_MultiOrgan_Proteomics"
d = pd.read_excel(os.path.join(DL, "mmc3.xlsx"), sheet_name="A_Protein matrix")
d.columns = [str(c) for c in d.columns]
lung_cols = [c for c in d.columns if c.startswith("lung_")]
print("lung 列数:", len(lung_cols))
print("示例:", lung_cols[:8])

# 分组：127N..133N 系列为非COVID，127C..133C 为COVID（TMT 通道命名）
N_cols = [c for c in lung_cols if re.match(r"lung_b\d+_\d+N$", c)]
C_cols = [c for c in lung_cols if re.match(r"lung_b\d+_\d+C$", c)]
print(f"N列 {len(N_cols)}, C列 {len(C_cols)}")

# 关键基因
KEY = ["FTH1", "FTL", "GPX4", "SLC40A1", "ITPR1", "HMOX1", "GSR", "GCLM",
       "NQO1", "NFE2L2", "ATF4", "DDIT3", "HIF1A", "SOD2", "TFRC", "TOMM20",
       "VDAC1", "HSPA9", "CANX", "GCLC", "PRDX3", "TXN", "TXNRD2", "CAT",
       "SLC11A2", "SLC25A37", "SLC25A28", "IREB2", "ACO1", "NCOA4"]
import re
out = []
for g_ in KEY:
    rows = d[d["Gene name"].astype(str).str.upper() == g_]
    if len(rows) == 0:
        # 也接受含逗号/分号的形式
        rows = d[d["Gene name"].astype(str).str.upper().str.contains(
            r"^" + g_ + r"[;,]", regex=True)]
    if len(rows) == 0:
        print(f"{g_}: 未在矩阵中")
        continue
    vals = []
    for _, r in rows.iterrows():
        nv = pd.to_numeric(r[N_cols], errors="coerce").dropna()
        cv = pd.to_numeric(r[C_cols], errors="coerce").dropna()
        if len(nv) < 2 or len(cv) < 2:
            continue
        vals.append((len(nv), len(cv), nv.mean(), cv.mean(), cv.mean() / nv.mean()))
    if not vals:
        print(f"{g_}: 肺组织定量值不足")
        continue
    # 取覆盖最全的一行
    best = max(vals, key=lambda x: x[0] + x[1])
    out.append({"gene": g_, "n_N": best[0], "n_C": best[1],
                "mean_N": best[2], "mean_C": best[3], "C_N_ratio": best[4]})
    print(f"{g_}: N均值={best[2]:.3f} C均值={best[3]:.3f} C/N={best[4]:.3f} (nN={best[0]}, nC={best[1]})")

res = pd.DataFrame(out)
res.to_csv(r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\01_FIGURE_DATA_CSV\M4_Nie_lung_mmc3_key_genes.csv", index=False)
print("\n已保存 M4_Nie_lung_mmc3_key_genes.csv")
