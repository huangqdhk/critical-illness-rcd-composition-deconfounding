# -*- coding: utf-8 -*-
"""
P0_manuscript_source_check_ext.py — §四-2 扩展版 manuscript-to-source 终检（2026-08-29）
==========================================================================================
在 P0_manuscript_source_check.py（17 项）基础上扩展 M1-M4 层关键句断言 + References 完整性 + 文稿侧包含性检查：
  M1 SMR 102 检验 BH≥0.65 / M1 pQTL BH≥0.752 / M2 S57 GSE32707 MDI q=0.0225 / D1b TOMM40 锚点 /
  S55 系列空间锚点（Stouffer z=6.33、残差化 4.04、SVG 9.36、21/23、CosMx 4/4）/
  S59a 双口径 33%/49% / Table 1 八行分值 / References 58 条完整性 / 查新重跑 08-29 切换
输出：04_AUDIT_GOVERNANCE/P0_Manuscript_Source_Check_Report_v2.md（原 08-21 报告保留不动）
"""
import os, re, sys, json
import pandas as pd
import numpy as np
from scipy.stats import norm
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
MS = os.path.join(ROOT, "111文稿_v5.md")
MSE = os.path.join(ROOT, "111文稿_v5_英文版.md")
INTER = os.path.join(ROOT, "_intermediate")
ST = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES", "SUPPLEMENTARY_Tables_CSV")
OUT = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "P0_Manuscript_Source_Check_Report_v2.md")

ms = open(MS, encoding="utf-8").read()
mse = open(MSE, encoding="utf-8").read()
checks = []
def add(tag, claim, verdict, detail):
    checks.append(dict(tag=tag, claim=claim, verdict=verdict, detail=detail))

# ================= M1 因果层 =================
smr = pd.read_csv(os.path.join(ST, "Table_S15f_SMR_HEIDI_Results.csv"))
n_smr, minbh = len(smr), smr["bh_p_pooled"].min()
ok = n_smr == 102 and abs(minbh - 0.6508) < 0.001 and minbh >= 0.65
add("M1-SMR102", "文稿：SMR/HEIDI 102 检验全部 BH≥0.65（最小 BH=0.6508）",
    "PASS" if ok else "MISMATCH", f"源表 n={n_smr}, min bh_p_pooled={minbh:.4f}")

pq = pd.read_csv(os.path.join(ST, "Table_S15e_pQTL_MR_Results.csv"))
minpq = pq["bh_p_pooled"].min()
itpr1 = pq[(pq["exposure"] == "UKB-PPP:ITPR1") & (pq["outcome"] == "IEU_sepsis_ieu-b-69")]
ok = minpq >= 0.752 and len(itpr1) and abs(itpr1.iloc[0]["OR"] - 1.142) < 0.001 and abs(itpr1.iloc[0]["pval"] - 0.406) < 0.001
add("M1-pQTL", "文稿：pQTL-MR 全部 BH≥0.752；ITPR1 Wald OR=1.142, p=0.406",
    "PASS" if ok else "MISMATCH", f"源表 n={len(pq)}, min bh={minpq:.4f}, ITPR1 OR={itpr1.iloc[0]['OR']:.3f}, p={itpr1.iloc[0]['pval']:.3f}" if len(itpr1) else "无ITPR1行")

# ================= M2 临床化 =================
s57 = pd.read_csv(os.path.join(ST, "Table_S57_M2_CrossDisease_Contrasts_Meta.csv"))
g = s57[(s57["cohort"] == "GSE32707") & (s57["score"] == "MDI") & (s57["contrast"] == "ARDSd0_vs_Control")]
ok = len(g) == 1 and abs(g.iloc[0]["BH_q"] - 0.0225) < 0.0005 and abs(g.iloc[0]["hedges_g"] - 0.772) < 0.005
add("M2-S57-GSE32707-MDI", "文稿：GSE32707 ARDS d0 MDI g=+0.77（p=0.015），BH 校正 q=0.0225",
    "PASS" if ok else "MISMATCH", f"源表 g={g.iloc[0]['hedges_g']:.3f}, MW_p={g.iloc[0]['MW_p']:.4f}, BH_q={g.iloc[0]['BH_q']:.5f}" if len(g) else "无记录")

# ================= D1b 三锚点（TOMM40） =================
s53 = pd.read_csv(os.path.join(ST, "Table_S53_D1b_Three_Anchor_Consistency.csv"))
t40 = s53[s53["row_id"] == "TOMM40"]
ok = len(t40) == 1 and t40.iloc[0]["GSE32707_ARDS_d0_vs_Control"] < 0 and t40.iloc[0]["GSE32707_ARDS_d0_vs_Control_BH_q"] < 0.05
add("D1b-TOMM40", "文稿：GSE32707 逐基因显著下调含 TOMM40（14 个 BH q<0.05 之一）",
    "PASS" if ok else "MISMATCH",
    f"源表 TOMM40 log2FC={t40.iloc[0]['GSE32707_ARDS_d0_vs_Control']:.3f}, BH_q={t40.iloc[0]['GSE32707_ARDS_d0_vs_Control_BH_q']:.4f}" if len(t40) else "无记录")

# ================= M1 空间（S55 系列） =================
a = pd.read_csv(os.path.join(ST, "Table_S55a_M1_Visium_Section_Spatial_Stats.csv"))
z_bv = a["z_bv"].sum() / np.sqrt(len(a))
neg_sig = int(((a["z_bv"] < 0) & (a["p_bv"] < 0.05)).sum())
ex_sig = int((a["p_ex"] < 0.05).sum())
ok = len(a) == 23 and abs(z_bv - 6.33) < 0.01 and neg_sig == 3 and ex_sig == 21
add("M1-S55a-bivariate", "文稿：双变量 Moran's I 23 切片 Stouffer 合并 z=+6.33；仅 3/23 显著负；执行臂 21/23 显著聚集",
    "PASS" if ok else "MISMATCH", f"源表复算 z={z_bv:.3f}, 负显著 {neg_sig}/23, 执行臂显著 {ex_sig}/23")

v = pd.read_csv(os.path.join(ST, "Table_S55v_M1_Visium_Bivariate_Residualized.csv"))
p_res = v["p_resid"].clip(lower=1e-3)
z_res = (np.sign(v["I_bv_resid"]) * norm.ppf(1 - p_res / 2)).sum() / np.sqrt(len(v))
ok = len(v) == 23 and abs(z_res - 4.04) < 0.01
add("M1-S55v-residualized", "文稿：计数残差化后合并 z=+4.04（p=5.4×10⁻⁵）",
    "PASS" if ok else "MISMATCH", f"源表复算（置换 p 下限 1e-3 截断）z={z_res:.3f}, n={len(v)}")

f = pd.read_csv(os.path.join(ST, "Table_S55f_M1_Visium_SVG_Enrichment.csv"))
exr = f[f["arm"] == "execution_induction"].iloc[0]
ok = abs(exr["stouffer_z"] - 9.36) < 0.01 and exr["n_sig_sections"] == 14 and exr["n_sections"] == 23
add("M1-S55f-SVG", "文稿：执行臂 SVG 富集 Stouffer z=+9.36（14/23 切片显著）",
    "PASS" if ok else "MISMATCH", f"源表 z={exr['stouffer_z']:.3f}, {exr['n_sig_sections']}/{exr['n_sections']}")

j = pd.read_csv(os.path.join(ST, "Table_S55j_M1_CosMx_Moran.csv"))
infj = j[j["feature"] == "inflammasome_module"]
ok = len(infj) == 4 and (infj["p_perm"] == 0.001).all() and infj["I"].between(0.022, 0.042).all()
add("M1-S55j-CosMx", "文稿：CosMx 炎症小体模块 4/4 TMA 显著聚集（置换 p=0.001）",
    "PASS" if ok else "MISMATCH", f"源表 4 TMA p_perm 均为 {sorted(infj['p_perm'].unique())}, I 范围 {infj['I'].min():.3f}~{infj['I'].max():.3f}")

# ================= M4 生化层（S59a 双口径） =================
sa = pd.read_csv(os.path.join(ST, "Table_S59a_M4_Protein_Transcript_Consistency.csv"))
A = sa[sa["dir_blood"].notna() & sa["dir_protein"].notna() & (sa["blood_tx_padj"] < 0.05)]
B = sa[sa["dir_balf"].notna() & sa["dir_protein"].notna() & (sa["balf_padj"] < 0.05)]
rA, rB = (A["dir_blood"] == A["dir_protein"]).mean(), (B["dir_balf"] == B["dir_protein"]).mean()
ok = len(A) == 30 and (A["dir_blood"] == A["dir_protein"]).sum() == 10 and len(B) == 37 and (B["dir_balf"] == B["dir_protein"]).sum() == 18
add("M4-S59a-dual", "文稿：口径 A 可判定 30 基因一致 10（33%）；口径 B 可判定 37 一致 18（49%）",
    "PASS" if ok else "MISMATCH", f"源表复算 A={len(A)}判/{int((A['dir_blood']==A['dir_protein']).sum())}致（{rA:.0%}）, B={len(B)}判/{int((B['dir_balf']==B['dir_protein']).sum())}致（{rB:.0%}）")

# ================= Table 1 八行分值 =================
s1 = pd.read_csv(os.path.join(ST, "Table_S1_Mitoxyperilysis_Score_by_CellType_Summary.csv"))
h = s1[s1["condition"] == "Healthy"]
pooled = (h.assign(w=h["mean_score"] * h["n_cells"]).groupby("cell_type")[["w", "n_cells"]].sum()
          .assign(pooled=lambda d: d["w"] / d["n_cells"])["pooled"].round(4))
expect = {"Macrophage": 0.4568, "Mono_c14": 0.4204, "DC": 0.3728, "Epithelial": 0.3717,
          "B_cell": 0.3397, "T_cell": 0.3161, "Club": 0.3154, "NK": 0.3050}
bad = {k: (pooled.get(k), v) for k, v in expect.items() if abs((pooled.get(k) or 9) - v) > 0.0005}
ok = not bad
add("Table1-eight", "文稿 Table 1：健康对照八细胞类型分值（分数据集评分后合并细胞池均值）0.4568/0.4204/0.3728/0.3717/0.3397/0.3161/0.3154/0.3050",
    "PASS" if ok else "MISMATCH", "八类型 n_cells 加权合并全部吻合" if ok else f"不符：{bad}")

# ================= References 完整性（2026-08-29 插入后合并复核） =================
def ref_profile(text):
    body = text.split("# References")[0]
    marks = [m for m in re.findall(r"\[[\d,\-\s]+\]", body) if m.strip() != "[0,1]"]
    expand = set()
    for m in marks:
        for part in m[1:-1].split(","):
            part = part.strip()
            if "-" in part:
                lo, hi = part.split("-"); expand |= set(range(int(lo), int(hi) + 1))
            elif part:
                expand.add(int(part))
    n_entries = len(re.findall(r"^\d+\. ", text.split("# References")[1], re.M))
    return marks, expand, n_entries

mz, ez, nz = ref_profile(ms)
me, ee, ne = ref_profile(mse)
ok = nz == 58 and ne == 58 and ez == set(range(1, 59)) and ee == set(range(1, 59)) and mz == me
add("REF-integrity", "References 插入后完整性：双版各 58 条、正文标记展开覆盖 1..58、双版标记序列一致",
    "PASS" if ok else "MISMATCH", f"中文 {nz} 条/英文 {ne} 条；覆盖一致={ez == ee == set(range(1,59))}；标记序列一致={mz == me}")

# ================= 文稿侧包含性（关键数字串必须在正文出现） =================
needles = ["0.6508", "0.752", "0.0225", "6.33", "9.36", "4.04", "0.4568", "0.3050", "33%", "49%",
           "21/23", "4/4 TMA", "p=0.001", "TOMM40", "1.142", "0.406"]
missing = [x for x in needles if x not in ms]
add("MS-presence", "文稿侧包含性：16 个关键数字锚点均在 v5 正文出现",
    "PASS" if not missing else "MISMATCH", "全部命中" if not missing else f"缺失：{missing}")

# ================= 查新（切换 08-29） =================
nov = json.load(open(os.path.join(INTER, "pubmed_novelty_rerun_20260829.json"), encoding="utf-8"))
ok = nov["Q1_term"]["count"] == "11" and nov["Q2_ARDS_lung"]["count"] == "0" and nov["Q3_sepsis"]["count"] == "1"
add("P0-novelty-0829", "查新重跑（2026-08-29）：Q1=11、Q2(ARDS/肺)=0、Q3(脓毒症)=1——空白维持",
    "PASS" if ok else "MISMATCH", f"Q1={nov['Q1_term']['count']}, Q2={nov['Q2_ARDS_lung']['count']}, Q3={nov['Q3_sepsis']['count']}")

# ================= 汇总 =================
n_pass = sum(1 for c in checks if c["verdict"] == "PASS")
n_mm = sum(1 for c in checks if c["verdict"] == "MISMATCH")
lines = ["# P0 Manuscript-to-Source 自动核对报告 v2（扩展版，2026-08-29）", "",
         "> 检查对象：`111文稿_v5.md`（References 已按 Vancouver 插入后）+ 英文版一致性。",
         "> 本版为 §四-2 扩展终检：在原 17 项（P0/P1/P2，见 08-21 版报告，保留不动）之上新增 M1-M4 层关键句、References 完整性与文稿侧包含性检查。",
         f"> 结果：**{n_pass} PASS / {n_mm} MISMATCH / {len(checks)} 项（本版新增项）**。",
         "", "| # | 断言 | 判定 | 详情 |", "|---|---|---|---|"]
for i, c in enumerate(checks, 1):
    lines.append(f"| {i} | {c['claim']} | {c['verdict']} | {c['detail']} |")
lines += ["", "## 结论", "",
          "- 本版覆盖 §四-2 指定全部锚点：SMR 102/BH≥0.65、pQTL BH≥0.752、S57 GSE32707 MDI q=0.0225、TOMM40、S55 系列空间锚点（Stouffer 6.33/残差化 4.04/SVG 9.36/21-23/4-4 TMA）、S59a 双口径 33%-49%、Table 1 八行分值。",
          "- 新增 References 插入后完整性断言（双版 58 条、覆盖 1..58、标记序列一致）与文稿侧包含性检查（16 锚点数字串）。",
          "- 查新断言已切换至 2026-08-29 重跑结果（空白维持：Q1=11、×ARDS/肺=0、×脓毒症=1）。",
          "- 连同原 17 项（08-21 版全 PASS），manuscript-to-source 核对合计覆盖 P0/P1/P2+M1-M4+References。",
          "- 若有 MISMATCH 项，须在正式稿中按源表修正后再投稿。", ""]
open(OUT, "w", encoding="utf-8").write("\n".join(lines))
print(f"{n_pass} PASS / {n_mm} MISMATCH / {len(checks)}")
for c in checks:
    print(f"[{c['verdict']}] {c['tag']}: {c['detail']}")
