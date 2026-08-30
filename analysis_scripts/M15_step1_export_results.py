# -*- coding: utf-8 -*-
"""
M15_step1_export_results.py — M10/M13/M14 结果导出至 01_FIGURE_DATA_CSV/Main 与
02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV（命名合规：Figure_Nx.csv / Table_S6x_*.csv）
"""
import os
import sys
import shutil

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
INTER = ROOT + r"\_intermediate"
MAIN = ROOT + r"\01_FIGURE_DATA_CSV\Main"
SUPP = ROOT + r"\02_SUPPLEMENTARY_TABLES\SUPPLEMENTARY_Tables_CSV"
GOV = ROOT + r"\04_AUDIT_GOVERNANCE"

log = []
def note(m=""):
    log.append(m); print(m, flush=True)

def read(p, **kw):
    return pd.read_csv(os.path.join(INTER, p), **kw)

# ======================================================================
# Figure 7 主图面板数据
# ======================================================================
note("== Figure 7 面板数据")

# --- 7A：M10 检验 A 组成预测 R²（Monaco|NNLS 主口径 + 4 套全口径） ---
r2 = read("M10A_r2_comparison.csv")
merged = read("M10A_r2_merged.csv")
cross = read("M10A_proportion_cross_reference.csv")
fig7a = r2.copy()
fig7a["source"] = "M10A_r2_comparison"
fig7a.to_csv(os.path.join(MAIN, "Figure_7A.csv"), index=False)
note(f"   Figure_7A.csv：{len(fig7a)} 行（4 套口径逐队列 R²）")

# --- 7B：M10 检验 A 残差疾病效应（k=16 层 + 合并 + LOSO） ---
layers = read("M10A_residual_layer_effects.csv")
loso = read("M10A_residual_meta_loso.csv")
meta = read("M10A_residual_meta.csv") if os.path.exists(os.path.join(INTER, "M10A_residual_meta.csv")) else None
fig7b = layers.copy()
fig7b.to_csv(os.path.join(MAIN, "Figure_7B.csv"), index=False)
note(f"   Figure_7B.csv：{len(fig7b)} 层残差效应")

# --- 7C：M10 检验 B 髓系 UCS 层 + 合并 ---
b_layers = read("M10B_myeloid_UCS_layers.csv")
b_meta = read("M10B_myeloid_UCS_meta.csv") if os.path.exists(os.path.join(INTER, "M10B_myeloid_UCS_meta.csv")) else None
fig7c = b_layers.copy()
fig7c.to_csv(os.path.join(MAIN, "Figure_7C.csv"), index=False)
note(f"   Figure_7C.csv：{len(fig7c)} 层")

# --- 7D：M13 目录回归散点数据 ---
cat = read("M13_catalog.csv")
fig7d = cat[["signature", "n_genes", "n_measured", "execution_share_ext", "myeloid_share",
             "g_185", "g_327", "r_UCS", "r_EIS", "dominant_arm"]].copy()
fig7d.to_csv(os.path.join(MAIN, "Figure_7D.csv"), index=False)
note(f"   Figure_7D.csv：{len(fig7d)} 签名目录")

# --- 7E：M14 腿2 core Gate 2 三层 ---
g32707 = read("M14_gate2_core_scores_GSE32707.csv")
g185 = read("M14_gate2_core_scores_GSE185263.csv")
scrna = read("M14_gate2_core_scRNA_localization.csv")
fig7e = pd.concat([
    g32707.assign(layer="GSE32707"),
    g185.assign(layer="GSE185263"),
    scrna.assign(layer="scRNA_donor"),
], ignore_index=True, sort=False)
fig7e.to_csv(os.path.join(MAIN, "Figure_7E.csv"), index=False)
note(f"   Figure_7E.csv：{len(fig7e)} 行（三层头对头）")

# --- 7F：M14 腿1 Torin 置换 ---
torin = read("M14_Torin_reversal_permutation.csv")
fig7f = torin.copy()
fig7f.to_csv(os.path.join(MAIN, "Figure_7F.csv"), index=False)
note(f"   Figure_7F.csv：{len(fig7f)} 行")

# --- 7G：M10 检验 C 区室对照 ---
c1 = read("M10C_sputum_blood.csv")
c2 = read("M10C_pbmc_ta.csv") if os.path.exists(os.path.join(INTER, "M10C_pbmc_ta.csv")) else pd.DataFrame()
c3 = read("M10C_balf_pbmc.csv")
fig7g = pd.concat([
    c1.assign(contrast="GSE148871_sputum_vs_blood"),
    c2.assign(contrast="GSE180578_TA_vs_PBMC"),
    c3.assign(contrast="BALF_vs_PBMC_cross_cohort"),
], ignore_index=True, sort=False)
fig7g.to_csv(os.path.join(MAIN, "Figure_7G.csv"), index=False)
note(f"   Figure_7G.csv：{len(fig7g)} 行")

# ======================================================================
# Supplementary tables（编号 S62 起）
# ======================================================================
note("\n== Supplementary tables（S62 起）")

# S62：IIAMD core 签名（小鼠+人源映射）
core = read("M14_IIAMD_core_v1.0.csv")
core.to_csv(os.path.join(SUPP, "Table_S62_M14_IIAMD_Core_v1.0.csv"), index=False)
note(f"   Table_S62（core {len(core)} 基因）")

# S63：M14 Gate 2 core 三层重跑（全表 = Figure 7E 同源 + n 列齐全）
fig7e.to_csv(os.path.join(SUPP, "Table_S63_M14_Gate2_Core_Rerun_ThreeLayers.csv"), index=False)

# S64：M10A 反卷积比例（4 套）
props = read("M10A_proportions.csv")
props.to_csv(os.path.join(SUPP, "Table_S64_M10A_Deconvolution_Proportions.csv"), index=False)
note(f"   Table_S64（{len(props)} 行）")

# S65：M10A R² 全口径 + 合并 + 交叉验证
r2.to_csv(os.path.join(SUPP, "Table_S65_M10A_R2_AllMethods.csv"), index=False)
merged.to_csv(os.path.join(SUPP, "Table_S65b_M10A_R2_Merged.csv"), index=False)
cross.to_csv(os.path.join(SUPP, "Table_S65c_M10A_Proportion_CrossReference.csv"), index=False)

# S66：M10A 残差
layers.to_csv(os.path.join(SUPP, "Table_S66_M10A_Residual_LayerEffects.csv"), index=False)
loso.to_csv(os.path.join(SUPP, "Table_S66b_M10A_Residual_Meta_LOSO.csv"), index=False)

# S67：M10B 三队列供者级臂评分与检验
b_scores = pd.concat([
    read("M10B_gse216009_arm_scores.csv").assign(cohort="GSE216009"),
    read("M10B_gse180578_arm_scores.csv").assign(cohort="GSE180578"),
    pd.read_csv(os.path.join(INTER, "P0_pseudobulk_arm_scores.csv")).assign(cohort=None),
], ignore_index=True, sort=False)
b_tests = pd.concat([
    read("M10B_gse216009_arm_tests.csv"),
    read("M10B_gse180578_arm_tests.csv"),
], ignore_index=True, sort=False)
b_scores.to_csv(os.path.join(SUPP, "Table_S67_M10B_Donor_Arm_Scores.csv"), index=False)
b_tests.to_csv(os.path.join(SUPP, "Table_S67b_M10B_Donor_Arm_Tests.csv"), index=False)
b_layers.to_csv(os.path.join(SUPP, "Table_S67c_M10B_Myeloid_UCS_Layers.csv"), index=False)
if b_meta is not None:
    b_meta.to_csv(os.path.join(SUPP, "Table_S67d_M10B_Myeloid_UCS_Meta.csv"), index=False)

# S68：M13 签名目录
cat.to_csv(os.path.join(SUPP, "Table_S68_M13_Signature_Catalog.csv"), index=False)
note(f"   Table_S68（{len(cat)} 签名）")

# S69：M13 零模型 + 反向审计
nulldf = read("M13_null_model.csv")
# 零模型压缩为摘要（10,000 行全存过大但可接受？——导出为压缩版 gzip 与摘要）
nulldf.to_csv(os.path.join(SUPP, "Table_S69a_M13_NullModel_10000.csv.gz"), index=False, compression="gzip")
nulldf.describe().T.reset_index().to_csv(os.path.join(SUPP, "Table_S69b_M13_NullModel_Summary.csv"), index=False)
audit = read("M13_reverse_audit.csv")
audit.to_csv(os.path.join(SUPP, "Table_S69c_M13_ReverseAudit.csv"), index=False)

# S70：M10C 区室对照
fig7g.to_csv(os.path.join(SUPP, "Table_S70_M10C_Compartment_Contrasts.csv"), index=False)

# S71：M14 全宇宙效应表（15,304 基因；b3/b_lps/b_cs/b_torin）
univ = read("M14_GSE235046_universe_effects.csv")
univ.to_csv(os.path.join(SUPP, "Table_S71_M14_GSE235046_Universe_Effects.csv.gz"),
            index=False, compression="gzip")
note(f"   Table_S71（{len(univ)} 基因，gzip）")

with open(os.path.join(ROOT, "03_LOGS", "M15_export_log.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(log))
note("\nDONE export")
