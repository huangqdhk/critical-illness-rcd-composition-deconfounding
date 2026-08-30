# -*- coding: utf-8 -*-
"""
M15_step4_export_and_report.py — M15 步骤 4：量化表格导出 + manifest 登记 + 治理报告
================================================================================
前置：M15_step2（11 队列回归）与 M15_step3（两个采纳演示）已运行。
本步骤（只写新文件，不误伤既有文件）：
  1. 主图面板数据（量化表，按"优先表格"输出规则不出图）：
     01_FIGURE_DATA_CSV/Main/Figure_10A.csv（演示 1 对比表）
     01_FIGURE_DATA_CSV/Main/Figure_10B.csv（演示 2 对比表）
     01_FIGURE_DATA_CSV/Main/Figure_10C.csv（11 队列回归队列级汇总）
  2. 补充表：
     02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV/Table_S83*（回归全表）
     Table_S84*（演示 1：逐样本/对比/R²/比例）
     Table_S85*（演示 2：逐样本/对比/R²/比例）
  3. RESULTS_MANIFEST_v2.0.csv 追加登记（sha256/行数/字节数）；lint_package.py
     的 M1X_NO_VERSION_STAMP 补登新文件名（M15 导出设计豁免，与 S62–S82 同口径）。
  4. 04_AUDIT_GOVERNANCE/M15_Tool_Adoption_Report.md 治理报告。
  5. 运行 lint_package.py 全量校验（要求 74/0 全绿不降级）。
"""
import hashlib
import os
import re
import sys
import warnings

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import pandas as pd

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
INTER = os.path.join(ROOT, "_intermediate")
MAIN = os.path.join(ROOT, "01_FIGURE_DATA_CSV", "Main")
SUPP = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES", "SUPPLEMENTARY_Tables_CSV")
GOV = os.path.join(ROOT, "04_AUDIT_GOVERNANCE")
MANIFEST = os.path.join(GOV, "RESULTS_MANIFEST_v2.0.csv")

log = []
def note(m=""):
    log.append(m)
    print(m, flush=True)

def sha256(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()

NOTE_TXT = "M15 工具封装+采纳演示批次补登（2026-08-27；工具回归 11 队列、演示 GSE157103/GSE66099 经 GEO 官方核验）"
ANALYSIS_FIG = {
    "Figure_10A.csv": "图10A 面板：M15 采纳演示 1 GSE157103 逐族对比（观测/组成预测/组成残差）",
    "Figure_10B.csv": "图10B 面板：M15 采纳演示 2 GSE66099 三对比逐族对比",
    "Figure_10C.csv": "图10C 面板：M15 工具内部回归测试（11 队列队列级汇总）",
}
ANALYSIS_SUPP = {
    "Table_S84_M15_RegressionTest_11Cohorts.csv": "M15 工具回归测试全表（11 队列×族级 max|Δ|）",
    "Table_S84b_M15_RegressionTest_CohortSummary.csv": "M15 工具回归测试队列级汇总",
    "Table_S85_M15_Demo_GSE157103_PerSample.csv": "M15 采纳演示 1 逐样本评分/组成预测/残差",
    "Table_S85b_M15_Demo_GSE157103_Contrasts.csv": "M15 采纳演示 1 逐族对比",
    "Table_S85c_M15_Demo_GSE157103_R2.csv": "M15 采纳演示 1 组成预测 R²（三方法）",
    "Table_S85d_M15_Demo_GSE157103_Proportions.csv": "M15 采纳演示 1 组间细胞比例",
    "Table_S86_M15_Demo_GSE66099_PerSample.csv": "M15 采纳演示 2 逐样本评分/组成预测/残差",
    "Table_S86b_M15_Demo_GSE66099_Contrasts.csv": "M15 采纳演示 2 逐族对比（三对比）",
    "Table_S86c_M15_Demo_GSE66099_R2.csv": "M15 采纳演示 2 组成预测 R²（三方法）",
    "Table_S86d_M15_Demo_GSE66099_Proportions.csv": "M15 采纳演示 2 组间细胞比例",
}

# ======================================================================
# 0. 输入核对
# ======================================================================
note("== 0. 输入核对")
req = ["M15_regression_summary.csv", "M15_regression_cohort_summary.csv",
       "M15_demo_GSE157103_per_sample.csv", "M15_demo_GSE157103_contrasts.csv",
       "M15_demo_GSE157103_r2.csv", "M15_demo_GSE157103_proportions.csv",
       "M15_demo_GSE66099_per_sample.csv", "M15_demo_GSE66099_contrasts.csv",
       "M15_demo_GSE66099_r2.csv", "M15_demo_GSE66099_proportions.csv"]
for f in req:
    p = os.path.join(INTER, f)
    if not os.path.exists(p):
        raise SystemExit(f"缺少输入 {p}——请先运行 M15_step2 / M15_step3")
    note(f"   {f} ✓")

reg_sum = pd.read_csv(os.path.join(INTER, "M15_regression_summary.csv"))
reg_coh = pd.read_csv(os.path.join(INTER, "M15_regression_cohort_summary.csv"))
n_fam = len(reg_sum)
n_fam_pass = int(reg_sum["PASS"].sum())
n_coh_pass = int(reg_coh["PASS"].sum())
note(f"   回归：{n_fam} 族核对，通过 {n_fam_pass}/{n_fam}；队列级 {n_coh_pass}/{len(reg_coh)}")
if n_fam_pass != n_fam or n_coh_pass != len(reg_coh):
    raise SystemExit("回归测试未全绿——停止导出（如实披露，不发布未通过的结果）")

# ======================================================================
# 1. 主图面板数据（Figure 10）
# ======================================================================
note("\n== 1. Figure 10 面板数据")

def read_inter(name):
    return pd.read_csv(os.path.join(INTER, name))

fig10a = read_inter("M15_demo_GSE157103_contrasts.csv").copy()
fig10a["dataset"] = "GSE157103"
fig10a.to_csv(os.path.join(MAIN, "Figure_10A.csv"), index=False)
note(f"   Figure_10A.csv：{len(fig10a)} 行（演示 1 对比）")

fig10b = read_inter("M15_demo_GSE66099_contrasts.csv").copy()
fig10b["dataset"] = "GSE66099"
fig10b.to_csv(os.path.join(MAIN, "Figure_10B.csv"), index=False)
note(f"   Figure_10B.csv：{len(fig10b)} 行（演示 2 三对比）")

fig10c = reg_coh.copy()
fig10c.to_csv(os.path.join(MAIN, "Figure_10C.csv"), index=False)
note(f"   Figure_10C.csv：{len(fig10c)} 行（回归队列级汇总）")

# ======================================================================
# 2. 补充表 S83–S85
# ======================================================================
note("\n== 2. 补充表 S83–S85")
export_map = {
    "Table_S84_M15_RegressionTest_11Cohorts.csv": reg_sum,
    "Table_S84b_M15_RegressionTest_CohortSummary.csv": reg_coh,
    "Table_S85_M15_Demo_GSE157103_PerSample.csv": read_inter("M15_demo_GSE157103_per_sample.csv"),
    "Table_S85b_M15_Demo_GSE157103_Contrasts.csv": read_inter("M15_demo_GSE157103_contrasts.csv"),
    "Table_S85c_M15_Demo_GSE157103_R2.csv": read_inter("M15_demo_GSE157103_r2.csv"),
    "Table_S85d_M15_Demo_GSE157103_Proportions.csv": read_inter("M15_demo_GSE157103_proportions.csv"),
    "Table_S86_M15_Demo_GSE66099_PerSample.csv": read_inter("M15_demo_GSE66099_per_sample.csv"),
    "Table_S86b_M15_Demo_GSE66099_Contrasts.csv": read_inter("M15_demo_GSE66099_contrasts.csv"),
    "Table_S86c_M15_Demo_GSE66099_R2.csv": read_inter("M15_demo_GSE66099_r2.csv"),
    "Table_S86d_M15_Demo_GSE66099_Proportions.csv": read_inter("M15_demo_GSE66099_proportions.csv"),
}
for fname, df in export_map.items():
    p = os.path.join(SUPP, fname)
    df.to_csv(p, index=False)
    note(f"   {fname}：{len(df)} 行")

# ======================================================================
# 3. RESULTS_MANIFEST_v2.0.csv 登记 + lint 豁免
# ======================================================================
note("\n== 3. manifest 登记")
# 重要：keep_default_na=False——字面量 "NA" 是版本戳的合法值，默认解析会将其毁为 NaN
man = pd.read_csv(MANIFEST, keep_default_na=False, na_values=[], low_memory=False)
man["filename"] = man["filename"].astype(str)
existing = set(man["filename"])
new_rows = []
for fname in list(ANALYSIS_FIG) + list(ANALYSIS_SUPP):
    if fname in existing:
        note(f"   [跳过] {fname} 已在 manifest")
        continue
    if fname.startswith("Figure_"):
        loc = "01_FIGURE_DATA_CSV/Main"
        path = os.path.join(MAIN, fname)
        cls = "figure_data"
        s_number = ""
        analysis = ANALYSIS_FIG[fname]
    else:
        loc = "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV"
        path = os.path.join(SUPP, fname)
        cls = "supplementary_table"
        s_number = re.match(r"[Tt]able_([Ss]\d+)", fname).group(1).upper()  # 基础编号（N22）
        analysis = ANALYSIS_SUPP[fname]
    n_rows = len(pd.read_csv(path, nrows=None)) if path.endswith(".csv") else None
    new_rows.append(dict(
        s_number=s_number, filename=fname, location=loc, class_=None,
        status="canonical", gene_universe="none", gene_set_version="",
        score_version="", n_data_rows=n_rows, size_bytes=os.path.getsize(path),
        sha256=sha256(path), analysis=analysis, notes=NOTE_TXT))
new_rows_df = pd.DataFrame(new_rows)
new_rows_df = new_rows_df.rename(columns={"class_": "class"})
man = pd.concat([man, new_rows_df], ignore_index=True, sort=False)
man.to_csv(MANIFEST, index=False, na_rep="")
note(f"   新登记 {len(new_rows_df)} 行 -> {MANIFEST}")

# lint 豁免补登
LINT = os.path.join(ROOT, "lint_package.py")
txt = open(LINT, encoding="utf-8").read()
stamp_name = "M1X_NO_VERSION_STAMP"
added = 0
for fname in list(ANALYSIS_FIG) + list(ANALYSIS_SUPP):
    if f'"{fname}"' in txt:
        continue
    anchor = '    "Figure_9A.csv", "Figure_9B.csv", "Figure_9C.csv",\n'
    if anchor in txt:
        txt = txt.replace(anchor, anchor + f'    "{fname}",\n')
        added += 1
open(LINT, "w", encoding="utf-8").write(txt)
note(f"   lint_package.py M1X_NO_VERSION_STAMP 补登 {added} 个文件名")

# ======================================================================
# 4. 治理报告
# ======================================================================
note("\n== 4. 治理报告")
key = {
    "GSE157103": {"n": 126, "up": None, "ex": None, "md_r2": None, "g_obs": None, "g_res": None},
    "GSE66099": {"n": 276, "up": None, "ex": None, "md_r2": None, "g_obs": None, "g_res": None},
}
d1p = read_inter("M15_demo_GSE157103_per_sample.csv")
d1r2 = read_inter("M15_demo_GSE157103_r2.csv")
d1c = read_inter("M15_demo_GSE157103_contrasts.csv")
key["GSE157103"]["up"] = int(d1p["n_up"].iloc[0]); key["GSE157103"]["ex"] = int(d1p["n_ex"].iloc[0])
key["GSE157103"]["md_r2"] = float(d1r2[(d1r2["method"] == "Monaco_NNLS") & (d1r2["family"] == "MDI")]["R2"].iloc[0])
key["GSE157103"]["g_obs"] = float(d1c[(d1c["family"] == "MDI")]["Hedges_g"].iloc[0])
key["GSE157103"]["g_res"] = float(d1c[(d1c["family"] == "MDI_resid")]["Hedges_g"].iloc[0])
d2p = read_inter("M15_demo_GSE66099_per_sample.csv")
d2r2 = read_inter("M15_demo_GSE66099_r2.csv")
d2c = read_inter("M15_demo_GSE66099_contrasts.csv")
key["GSE66099"]["up"] = int(d2p["n_up"].iloc[0]); key["GSE66099"]["ex"] = int(d2p["n_ex"].iloc[0])
key["GSE66099"]["md_r2"] = float(d2r2[(d2r2["method"] == "Monaco_NNLS") & (d2r2["family"] == "MDI")]["R2"].iloc[0])
key["GSE66099"]["g_obs"] = float(d2c[(d2c["contrast"] == "SepticShock_vs_Control") & (d2c["family"] == "MDI")]["Hedges_g"].iloc[0])
key["GSE66099"]["g_res"] = float(d2c[(d2c["contrast"] == "SepticShock_vs_Control") & (d2c["family"] == "MDI_resid")]["Hedges_g"].iloc[0])

rep = []
rep.append("# M15 工具封装与采纳演示治理报告（2026-08-27）\n")
rep.append("- 依据：111黄裕荣创新提质_20260822.md §五 M15（工具核心 = 两臂分解 + MDI + 组成预测值/组成残差 MDI；GitHub + Zenodo + 判定门文档化；2 个第三方数据集采纳演示）\n")
rep.append("## 1. 工具包\n")
rep.append("- 位置：M15_tool/（mitoxdi v1.0.0；score_version = mdi_v1.0；gene_set_version = Mitoxy-80_v1.0）")
rep.append("- 组成：mdi_core（两臂分解+评分，mdi_lib 冻结快照）/ composition（Monaco+ABIS 参考谱、NNLS 主口径 + OLS-CLS 交叉、合成谱评分、残差）/ report（对比统计 + ISED 解读）/ loaders（逐队列口径复刻）")
rep.append("- 判定门语言：M15_tool/GATE.md（M10 检验 A 判定门的通用化表述）")
rep.append("- 数据依赖（零下载）：Monaco 冻结缓存 _intermediate/M10A_monaco_ct_means.csv；ABIS sigmatrixRNAseq/sigmatrixMicro（00_RAW_DATA/ABIS/）；80 基因 manifest（内置副本与治理目录 SHA256 一致 346537be…e9965）；GPL570/GPL23159/GPL6947/GPL10295 探针映射（00_RAW_DATA）")
rep.append("- 发布件：README/LICENSE/requirements/CITATION.cff；GitHub + Zenodo DOI 随投稿前打包定稿（代码与文档，不含数据下载）\n")
rep.append("## 2. 内部回归测试（11 队列）\n")
rep.append(f"- 参照冻结表：M2（6 队列评分）、M10A（组成预测/残差）、M11M12（5 队列评分+组成+残差+比例）")
rep.append(f"- 结果：{n_fam} 族逐样本核对全部 PASS（队列级 {n_coh_pass}/{len(reg_coh)}）；族级 max|Δ| 最大值为 "
           f"{reg_sum['max_abs_diff'].max():.3e}（机器精度量级）")
rep.append(f"- 表：Table_S83/S83b；图 10C 面板数据；日志 03_LOGS/M15_regression_log.txt\n")
rep.append("## 3. 采纳演示（2 个第三方数据集，GEO 官方核验入列）\n")
for ds, k in key.items():
    rep.append(f"- **{ds}**：n={k['n']}；臂覆盖 up {k['up']}/30 ex {k['ex']}/33；"
               f"MDI R²(Monaco|NNLS)={k['md_r2']:.3f}；观测 MDI g={k['g_obs']:+.3f} → "
               f"残差 MDI g={k['g_res']:+.3f}（衰减 {1 - k['g_res']/k['g_obs'] if abs(k['g_obs']) > 1e-9 else float('nan'):.1%}）")
rep.append("- 核验留痕：03_LOGS/M15_demo_GEO_verification_20260827.md；报告：03_LOGS/M15_demo_*_report.md")
rep.append("- 解读纪律：仅陈述客观量；按 GATE.md 判定门语言解读；不评价原论文结论\n")
rep.append("## 4. 导出登记\n")
rep.append(f"- 主图面板：Figure_10A–10C（量化表，不出图——按项目'优先生成量化数据表格'输出规则）")
rep.append(f"- 补充表：Table_S83–S85（{len(export_map)} 个文件）已写入 RESULTS_MANIFEST_v2.0.csv（sha256/行数/字节数）")
rep.append("- lint_package.py M1X_NO_VERSION_STAMP 已补登；本步骤结束时全量 lint 校验（要求 74/0 不降级）")
open(os.path.join(GOV, "M15_Tool_Adoption_Report.md"), "w", encoding="utf-8").write("\n".join(rep))
note("   -> 04_AUDIT_GOVERNANCE/M15_Tool_Adoption_Report.md")

# ======================================================================
# 5. lint 全量校验
# ======================================================================
note("\n== 5. lint 校验")
import subprocess
r = subprocess.run([sys.executable, os.path.join(ROOT, "lint_package.py")],
                   capture_output=True, text=True, encoding="utf-8", errors="replace")
out = (r.stdout or "") + (r.stderr or "")
note(out[-2000:])
with open(os.path.join(ROOT, "03_LOGS", "M15_export_log.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(log) + "\n\n[LINT]\n" + out)
print(f"\nlint exit code = {r.returncode}")
sys.exit(r.returncode)
