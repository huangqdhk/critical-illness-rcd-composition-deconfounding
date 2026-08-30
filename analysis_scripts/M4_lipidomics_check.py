# -*- coding: utf-8 -*-
"""M4 脂质层：MTBLS6844 COVID 血清脂质组（Mild/Moderate/Severe）类别偏移 + oxylipin 代理。"""
import pandas as pd
import io, sys, os, re
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
import numpy as np
from scipy.stats import mannwhitneyu

BASE = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
outdir = BASE + r"\01_FIGURE_DATA_CSV"
os.makedirs(outdir, exist_ok=True)

# ---- 1. MTBLS6844 ----
p = BASE + r"\00_RAW_DATA\MTBLS6844_COVID19_Serum_Lipidomics\m_MTBLS6844_FIA-MS_positive__metabolite_profiling_v2_maf.tsv"
m = pd.read_csv(p, sep="\t")
groups = {}
for c in m.columns:
    mm = re.match(r"(MILD|MODERATE|SEVERE)\d+$", str(c))
    if mm:
        groups.setdefault(mm.group(1), []).append(c)
print("样本数:", {k: len(v) for k, v in groups.items()})

# 每类脂质总和（相对丰度）
class_rows = {}
for cls, sub in m.groupby("Class"):
    class_rows[cls] = sub[groups["MILD"] + groups["MODERATE"] + groups["SEVERE"]].sum(axis=0)
C = pd.DataFrame(class_rows).T
C = C.div(C.sum(axis=0), axis=1)  # 类别占比
print("脂质类别:", len(C))

rows = []
for cls in C.index:
    for g1, g2 in [("SEVERE", "MILD"), ("SEVERE", "MODERATE"), ("MODERATE", "MILD")]:
        a = C.loc[cls, groups[g1]].astype(float)
        b = C.loc[cls, groups[g2]].astype(float)
        if len(a) < 3 or len(b) < 3:
            continue
        w, pval = mannwhitneyu(a, b, alternative="two-sided")
        rows.append({"class": cls, "comparison": f"{g1}_vs_{g2}",
                     f"{g1}_median": a.median(), f"{g2}_median": b.median(),
                     "delta": a.median() - b.median(), "p": pval})
res = pd.DataFrame(rows)
res["BH"] = res.groupby("comparison")["p"].transform(
    lambda x: x * len(x) / x.rank(method="first")).clip(upper=1.0)
res.to_csv(outdir + r"\M4_MTBLS6844_lipid_class_shift.csv", index=False)
print("\n===== MTBLS6844: SEVERE vs MILD 显著 (BH<0.05) =====")
sig = res[(res["comparison"] == "SEVERE_vs_MILD") & (res["BH"] < 0.05)].sort_values("p")
print(sig.to_string(index=False, max_colwidth=18) if len(sig) else "(无显著类别)")
print("\n全部 SEVERE_vs_MILD（按 p 排序前 8）:")
print(res[res["comparison"] == "SEVERE_vs_MILD"].sort_values("p").head(8)
      .to_string(index=False, max_colwidth=18))

# ---- 2. Shen2020 oxylipin 代理：mmc6 论文官方差异表 ----
p2 = BASE + r"\00_RAW_DATA\Shen2020_Cell_Proteome_Metabolome\mmc6.xlsx"
print(f"\n===== Shen2020 oxylipin 代理（mmc6 论文官方差异表）=====")
oxy_rows = []
for sn in ["Metab Severe vs Healthy", "Metab non-Severe vs Healthy",
           "Metab Severe vs non-Severe"]:
    t = pd.read_excel(p2, sheet_name=sn)
    t.columns = [str(c).strip() for c in t.columns]
    feat_col = t.columns[0]
    sub = t[t[feat_col].astype(str).str.contains(
        "HODE|HETE|HETrE|isoprostane|MDA|HNE|oxo-ODE", case=False, na=False)]
    for _, r in sub.iterrows():
        fd_col = [c for c in t.columns if str(c).strip().lower() in ("fd", "log2 (fold change)")]
        oxy_rows.append({"sheet": sn, "metabolite": r[feat_col],
                         "fd": r[fd_col[0]] if fd_col else None,
                         "p": r.get("P_value"), "padj": r.get("P_value_adjust")})
oxy = pd.DataFrame(oxy_rows)
if len(oxy):
    print(oxy.to_string(index=False, max_colwidth=30))
oxy.to_csv(outdir + r"\M4_Shen2020_oxylipin_proxy.csv", index=False)
print("\n已保存 M4_Shen2020_oxylipin_proxy.csv")
