# -*- coding: utf-8 -*-
"""
M4_finalize_outputs.py — M4 产出固化：
1. 从分析输出生成正式结果文件 Figure_10A–10D（Main/）与 Table_S59a–d、S60；
2. 把分析中间文件移入 _intermediate；
3. 向 FIGURE_DATA_MANIFEST.csv 与 RESULTS_MANIFEST_v2.0.csv 登记（SHA256）。
（Figure_10E / Table_S59e 由 M4_HPA_localization_check.py 直接生成，本脚本一并登记。）
"""
import io, os, sys, shutil, hashlib
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
import pandas as pd
import numpy as np

BASE = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
MAIN = BASE + r"\01_FIGURE_DATA_CSV\Main"
TAB = BASE + r"\02_SUPPLEMENTARY_TABLES\SUPPLEMENTARY_Tables_CSV"
FIGD = BASE + r"\01_FIGURE_DATA_CSV"
INT = BASE + r"\_intermediate"

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(65536), b""):
            h.update(b)
    return h.hexdigest()

# ---------- 1. 读分析输出 ----------
two = pd.read_csv(FIGD + r"\M4_Nie_protein_direction_two_refs.csv")
key = pd.read_csv(FIGD + r"\M4_Nie_lung_mmc3_key_genes.csv")
lip = pd.read_csv(FIGD + r"\M4_MTBLS6844_lipid_class_shift.csv")
oxy = pd.read_csv(FIGD + r"\M4_Shen2020_oxylipin_proxy.csv")
gem = pd.read_csv(INT + r"\M4_GEM_v1.1\M4_GEM_group_comparison.csv")
iron = pd.read_csv(INT + r"\M4_GEM_v1.1\M4_GEM_iron_balance.csv")

# 蛋白显著性（mmc4 口径）并到铁轴表
prev = pd.read_csv(FIGD + r"\M4_Nie_lung_protein_vs_transcript.csv")
prev = prev[["gene", "protein_log2FC(Nie_lung)", "protein_padj"]]
prev.columns = ["gene", "mmc4_log2FC", "mmc4_padj"]

# ---------- 2. 生成正式文件 ----------
# Figure_10A / Table_S59a：蛋白-转录双口径方向一致
a = two.copy()
a["dir_blood"] = np.sign(pd.to_numeric(a["blood_tx_log2FC"], errors="coerce"))
a["dir_balf"] = np.sign(pd.to_numeric(a["balf_d"], errors="coerce"))
a["dir_protein"] = np.sign(pd.to_numeric(a["px_log2FC"], errors="coerce"))
a["gene_set_version"] = "Mitoxy-80_v1.0"
a["score_version"] = "m4px_v1.0"
a.to_csv(MAIN + r"\Figure_10A.csv", index=False)
a.to_csv(TAB + r"\Table_S59a_M4_Protein_Transcript_Consistency.csv", index=False)

# Figure_10B / Table_S59b：肺铁轴蛋白
b = key.merge(prev, on="gene", how="left")
b["gene_set_version"] = "Mitoxy-80_v1.0"
b["score_version"] = "m4px_v1.0"
b.to_csv(MAIN + r"\Figure_10B.csv", index=False)
b.to_csv(TAB + r"\Table_S59b_M4_Lung_Iron_Axis_Proteins.csv", index=False)

# Figure_10C / Table_S59c：脂质类别偏移（SEVERE_vs_MILD 全类别）
c = lip[lip["comparison"] == "SEVERE_vs_MILD"].copy()
c = c.rename(columns={"SEVERE_median": "severe_median", "MILD_median": "mild_median"})
c["gene_set_version"] = "Mitoxy-80_v1.0"
c["score_version"] = "m4lipid_v1.0"
c.to_csv(MAIN + r"\Figure_10C.csv", index=False)
lip["gene_set_version"] = "Mitoxy-80_v1.0"
lip["score_version"] = "m4lipid_v1.0"
lip.to_csv(TAB + r"\Table_S59c_M4_Lipid_Class_Shift.csv", index=False)

# Table_S59d：oxylipin 代理
oxy["gene_set_version"] = "Mitoxy-80_v1.0"
oxy["score_version"] = "m4lipid_v1.0"
oxy.to_csv(TAB + r"\Table_S59d_M4_Oxylipin_Proxy.csv", index=False)

# Figure_10D / Table_S60：GEM 通量组间比较（标注容量饱和反应）
cap = pd.read_csv(INT + r"\M4_GEM_v1.1\M4_GEM_capacity_matrix.csv", index_col=0)
saturated = set(cap.index[cap.std(axis=1) == 0])
g = gem.copy()
g["capacity_saturated"] = g["reaction"].isin(saturated)
g["score_version"] = "m4gem_v1.1"
g.to_csv(MAIN + r"\Figure_10D.csv", index=False)
g.to_csv(TAB + r"\Table_S60_M4_GEM_Flux_Comparison.csv", index=False)

print("正式文件已生成。")

# ---------- 3. 中间文件移入 _intermediate ----------
for f in ["M4_Nie_lung_protein_vs_transcript.csv",
          "M4_Nie_lung_mmc3_key_genes.csv",
          "M4_Nie_protein_direction_two_refs.csv",
          "M4_MTBLS6844_lipid_class_shift.csv",
          "M4_Shen2020_oxylipin_proxy.csv"]:
    src = FIGD + "\\" + f
    if os.path.exists(src) and not os.path.exists(INT + "\\" + f):
        shutil.move(src, INT + "\\" + f)
        print("moved:", f)
if os.path.exists(FIGD + r"\M4_GEM") and not os.path.exists(INT + r"\M4_GEM_v1.0_degenerate"):
    shutil.move(FIGD + r"\M4_GEM", INT + r"\M4_GEM_v1.0_degenerate")
    print("moved: M4_GEM (v1.0 退化输出) -> _intermediate")

# ---------- 4. 登记 manifest ----------
fig_rows = [
    ("Figure_10A", "Main/Figure_10A.csv", "Table_S59a_M4_Protein_Transcript_Consistency.csv",
     "生化层补证(M4)", "蛋白-转录双口径方向一致(全血33%/BALF 49%;上游臂蛋白反向升高)",
     "双向散点/方向一致热图", "口径A=GSE185263全血log2FC;口径B=BALF d_mean;判定<60%阈值"),
    ("Figure_10B", "Main/Figure_10B.csv", "Table_S59b_M4_Lung_Iron_Axis_Proteins.csv",
     "生化层补证(M4)", "尸检肺铁轴蛋白定量(TFRC/HMOX1/FTH1/FTL/SLC40A1/SLC25A28↑,NQO1↓)",
     "C/N比值条形图(带显著性)", "mmc3原始矩阵19N/19C;mmc4官方统计并注"),
    ("Figure_10C", "Main/Figure_10C.csv", "Table_S59c_M4_Lipid_Class_Shift.csv",
     "生化层补证(M4)", "血清脂质类别占比 SEVERE vs MILD(无显著类别;CE名义p=0.043)",
     "类别偏移条形图", "MTBLS6844 FIA-MS 53样本"),
    ("Figure_10D", "Main/Figure_10D.csv", "Table_S60_M4_GEM_Flux_Comparison.csv",
     "生化层补证(M4)", "GEM容量组间比较(38反应无一BH<0.05;铁失衡p=0.34)",
     "火山图/森林图", "m4gem_v1.1;饱和反应标注;v1.0退化见03_LOGS审计记录"),
]
fig_df = pd.read_csv(FIGD + r"\FIGURE_DATA_MANIFEST.csv", encoding="utf-8-sig", on_bad_lines="skip")
for panel, pf, src, title, content, plot, note in fig_rows:
    fig_df.loc[len(fig_df)] = [panel, pf, src, title, content, plot, note]
fig_df.to_csv(FIGD + r"\FIGURE_DATA_MANIFEST.csv", index=False, encoding="utf-8-sig")
print("FIGURE_DATA_MANIFEST 登记完成")

res_df = pd.read_csv(BASE + r"\04_AUDIT_GOVERNANCE\RESULTS_MANIFEST_v2.0.csv",
                     encoding="utf-8-sig", on_bad_lines="skip")
res_specs = [
    ("Figure_10A.csv", "01_FIGURE_DATA_CSV/Main", "figure_data", "mitoxy_80", "Mitoxy-80_v1.0", "m4px_v1.0",
     "Figure 10A: 蛋白-转录双口径方向一致"),
    ("Figure_10B.csv", "01_FIGURE_DATA_CSV/Main", "figure_data", "mitoxy_80", "Mitoxy-80_v1.0", "m4px_v1.0",
     "Figure 10B: 尸检肺铁轴蛋白"),
    ("Figure_10C.csv", "01_FIGURE_DATA_CSV/Main", "figure_data", "mitoxy_80", "Mitoxy-80_v1.0", "m4lipid_v1.0",
     "Figure 10C: 脂质类别偏移 SEVERE vs MILD"),
    ("Figure_10D.csv", "01_FIGURE_DATA_CSV/Main", "figure_data", "mitoxy_80", "Mitoxy-80_v1.0", "m4gem_v1.1",
     "Figure 10D: GEM 容量组间比较"),
    ("Figure_10E.csv", "01_FIGURE_DATA_CSV/Main", "figure_data", "mitoxy_80", "Mitoxy-80_v1.0", "m4hpa_v1.0",
     "Figure 10E: 80基因亚细胞定位核验"),
    ("Table_S59a_M4_Protein_Transcript_Consistency.csv", "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
     "supplementary_table", "mitoxy_80", "Mitoxy-80_v1.0", "m4px_v1.0", "S59a: 蛋白-转录方向一致双口径"),
    ("Table_S59b_M4_Lung_Iron_Axis_Proteins.csv", "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
     "supplementary_table", "mitoxy_80", "Mitoxy-80_v1.0", "m4px_v1.0", "S59b: 肺铁轴蛋白定量"),
    ("Table_S59c_M4_Lipid_Class_Shift.csv", "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
     "supplementary_table", "mitoxy_80", "Mitoxy-80_v1.0", "m4lipid_v1.0", "S59c: 脂质类别偏移(3对比)"),
    ("Table_S59d_M4_Oxylipin_Proxy.csv", "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
     "supplementary_table", "mitoxy_80", "Mitoxy-80_v1.0", "m4lipid_v1.0", "S59d: oxylipin 代理(mmc6)"),
    ("Table_S59e_M4_Subcellular_Annotation.csv", "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
     "supplementary_table", "mitoxy_80", "Mitoxy-80_v1.0", "m4hpa_v1.0", "S59e: 亚细胞定位核验"),
    ("Table_S60_M4_GEM_Flux_Comparison.csv", "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
     "supplementary_table", "mitoxy_80", "Mitoxy-80_v1.0", "m4gem_v1.1", "S60: GEM 通量组间比较"),
]
for fn, loc, cls, uni, gsv, scv, note in res_specs:
    path = BASE + "\\" + loc.replace("/", "\\") + "\\" + fn
    if not os.path.exists(path):
        print("跳过（未生成）:", fn)
        continue
    nrow = len(pd.read_csv(path)) - 1
    size = os.path.getsize(path)
    res_df.loc[len(res_df)] = [np.nan, fn, loc, cls, "canonical", uni, gsv, scv,
                               nrow, size, sha256(path), note, "M4 (2026-08-18)"]
res_df.to_csv(BASE + r"\04_AUDIT_GOVERNANCE\RESULTS_MANIFEST_v2.0.csv", index=False,
              encoding="utf-8-sig")
print("RESULTS_MANIFEST_v2.0 登记完成")
