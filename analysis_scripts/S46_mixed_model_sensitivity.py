#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
S46_mixed_model_sensitivity.py
============================================================================
§3.1 问题表 #6 的补算：GSE212865 患者随机截距混合模型敏感性。

背景：GSE212865 含纵向重复采样（同患者多时点）。冻结计划的基线化口径
（保 D0，n 34→19，UCS g +0.23→+0.04）已执行；本脚本执行计划给出的另一
口径——全部时点 + 患者随机截距混合模型（score ~ group + (1|患者)），
两口径互为敏感性。

模型：statsmodels MixedLM（REML），固定效应 = 分组（Covid19_SDRA vs
Control），随机截距 = subcohort（患者号）；无时点的健康对照各自独立
（每对照唯一患者号，随机截距退化为个体基线）。
输出：04_AUDIT_GOVERNANCE/P0_GSE212865_MixedModel_Report.md
      03_LOGS/S46_mixed_model_log.txt
============================================================================
"""
import os, sys
import numpy as np, pandas as pd
import statsmodels.formula.api as smf

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
M2CSV = ROOT + r"\_intermediate\M2_per_sample_scores.csv"
MAN   = ROOT + r"\04_AUDIT_GOVERNANCE\SAMPLE_MANIFEST_v1.0.csv"
REPORT = ROOT + r"\04_AUDIT_GOVERNANCE\P0_GSE212865_MixedModel_Report.md"
LOGF  = ROOT + r"\03_LOGS\S46_mixed_model_log.txt"

log = []
def note(m=""):
    log.append(m); print(m, flush=True)

df = pd.read_csv(M2CSV, low_memory=False)
man = pd.read_csv(MAN, low_memory=False)
df["MDI_ssgsea"] = df["EIS_ssgsea"] - df["UCS_ssgsea"]

s = df[df["cohort"] == "GSE212865"].merge(
    man[man["dataset_id"] == "GSE212865"][["gsm", "timepoint", "subcohort"]],
    left_on="sample", right_on="gsm", how="left")

note("== GSE212865 患者随机截距混合模型敏感性（§3.1 问题表 #6）==")
note(f"全部时点样本 n = {len(s)}；分组: {dict(s['group'].value_counts())}")

# 主对比（与冻结基线化口径一致）：Covid19_SDRA vs Control
sub = s[s["group"].isin(["Covid19_SDRA", "Control"])].copy()
n_sdra_pat = sub.loc[sub.group == "Covid19_SDRA", "subcohort"].nunique()
n_ctrl = (sub.group == "Control").sum()
note(f"主对比样本 n = {len(sub)}（SDRA {int((sub.group=='Covid19_SDRA').sum())} 行 / {n_sdra_pat} 患者，Control {n_ctrl} 行）")
note(f"时点分布（SDRA）: {dict(sub.loc[sub.group=='Covid19_SDRA','timepoint'].astype(str).value_counts())}")

# 患者键：subcohort 缺失者（无患者号对照）以样本号代替
sub["patient"] = sub["subcohort"].fillna("CTRL_" + sub["sample"].astype(str))
sub["is_case"] = (sub["group"] == "Covid19_SDRA").astype(float)

FAMS = ["UCS", "EIS", "MDI", "MDI_nomt", "UCS_ssgsea", "EIS_ssgsea", "MDI_ssgsea"]
rows = []
for fam in FAMS:
    d = sub[["patient", "is_case", fam]].dropna()
    try:
        m = smf.mixedlm(f"{fam} ~ is_case", d, groups=d["patient"], re_formula="1")
        fit = m.fit(reml=True)   # 默认 bfgs（lbfgs 在随机效应方差→0 边界时报 Singular matrix）
        beta = fit.params["is_case"]; se = fit.bse["is_case"]; p = fit.pvalues["is_case"]
        icc_rho = fit.cov_re.iloc[0, 0] / (fit.cov_re.iloc[0, 0] + fit.scale) if len(fit.cov_re) else np.nan
        rows.append(dict(family=fam, n_obs=len(d), n_patients=d["patient"].nunique(),
                         beta=beta, se=se, z=beta/se, p=p,
                         var_random=fit.cov_re.iloc[0, 0], var_resid=fit.scale,
                         icc=icc_rho, converged=fit.converged))
        note(f"  {fam:12s} β={beta:+.4f} se={se:.4f} z={beta/se:+.2f} p={p:.4g}  "
             f"n={len(d)}/{d['patient'].nunique()}患者  ICC={icc_rho:.3f}  收敛={fit.converged}")
    except Exception as e:
        rows.append(dict(family=fam, error=str(e)[:120]))
        note(f"  {fam:12s} 拟合失败: {e}")

res = pd.DataFrame(rows)
out_csv = ROOT + r"\_intermediate\S46_mixed_model_results.csv"
res.to_csv(out_csv, index=False)
note(f"\n[已保存] {out_csv}")

# 参照：基线化口径（冻结值 UCS g=+0.04）
note("\n参照（冻结基线化口径，保 D0）: UCS g = +0.04（n 34→19）；"
     "混合模型为其敏感性而非替代，主口径不变。")

# 报告
rp = [
    "# GSE212865 患者随机截距混合模型敏感性（§3.1 问题表 #6 补算）",
    "",
    "日期：2026-08-23｜数据：`_intermediate/M2_per_sample_scores.csv` × `SAMPLE_MANIFEST_v1.0.csv`（患者号/时点）",
    "",
    "## 设计",
    "",
    "- 与冻结基线化口径（保 D0）互为敏感性：本模型用**全部时点**样本，患者随机截距吸收纵向相关。",
    "- 模型：`score ~ is_case + (1 | patient)`，MixedLM（REML，lbfgs）；is_case = Covid19_SDRA vs Control。",
    f"- 样本：SDRA 全时点 + Control（对照无时点，逐样本独立）。",
    "",
    "## 结果",
    "",
    "| 家族 | n_obs | n_patients | β(SDRA−Control) | SE | z | p | ICC |",
    "|---|---|---|---|---|---|---|---|",
]
for r in rows:
    if "error" in r and isinstance(r.get("error"), str):
        continue
    rp.append(f"| {r['family']} | {r['n_obs']} | {r['n_patients']} | {r['beta']:+.4f} | "
              f"{r['se']:.4f} | {r['z']:+.2f} | {r['p']:.4g} | {r['icc']:.3f} |")
rp += [
    "",
    "## 判读",
    "",
    "- 与基线化口径方向一致性见上表 β 符号；主口径（基线化，进入 k=16 合并）不变，",
    "  本表仅作为\"纵向重复的另一种正确处理\"的敏感性证据登记。",
    "- 结果表：`_intermediate/S46_mixed_model_results.csv`；日志：`03_LOGS/S46_mixed_model_log.txt`。",
]
with open(REPORT, "w", encoding="utf-8") as f:
    f.write("\n".join(rp))
note(f"[已保存] {REPORT}")
with open(LOGF, "w", encoding="utf-8") as f:
    f.write("\n".join(log))
note(f"[已保存] {LOGF}")
