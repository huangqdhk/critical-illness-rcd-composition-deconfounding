# -*- coding: utf-8 -*-
"""
M11M12_step4_export.py — M11/M12 结果导出（图数据 CSV + 附表 + manifest 登记 + 审计报告）
================================================================================
- 图数据：01_FIGURE_DATA_CSV/Main/Figure_8A–8D.csv、Figure_9A–9C.csv（量化表，不做图）
- 附表：02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV/Table_S76–S82
- 登记：FIGURE_DATA_MANIFEST.csv 追加；README_附表索引.md 追加章节
- 审计报告：04_AUDIT_GOVERNANCE/M11M12_Analysis_Report_20260827.md
"""
import os, sys, time
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
INTER = ROOT + r"\_intermediate"
MAIN = ROOT + r"\01_FIGURE_DATA_CSV\Main"
SUP_CSV = ROOT + r"\02_SUPPLEMENTARY_TABLES\SUPPLEMENTARY_Tables_CSV"
GOV = ROOT + r"\04_AUDIT_GOVERNANCE"
MANIFEST_F = ROOT + r"\01_FIGURE_DATA_CSV\FIGURE_DATA_MANIFEST.csv"

t0 = time.time()
log = []
def note(m=""):
    log.append(m); print(m, flush=True)

h1 = pd.read_csv(INTER + r"\M11M12_step2_h1_cohort.csv")
h2 = pd.read_csv(INTER + r"\M11M12_step2_h2_cohort.csv")
z1 = pd.read_csv(INTER + r"\M11M12_step2_zero_h1.csv")
z2 = pd.read_csv(INTER + r"\M11M12_step2_zero_h2.csv")
loso1 = pd.read_csv(INTER + r"\M11M12_step2_h1_loso_subject.csv")
loso2 = pd.read_csv(INTER + r"\M11M12_step2_h2_loso_subject.csv")
loso_meta1 = pd.read_csv(INTER + r"\M11M12_step2_h1_loso_meta.csv")
loso_meta2 = pd.read_csv(INTER + r"\M11M12_step2_h2_loso_meta.csv")
sens212 = pd.read_csv(INTER + r"\M11M12_step2_sensitivity_gse212865.csv")
gate = pd.read_csv(INTER + r"\M11M12_step2_gate.csv")
sens545 = pd.read_csv(INTER + r"\M11M12_step2_h1_gse54514_author_neut_sens.csv")
per_sample = pd.read_csv(INTER + r"\M11M12_step1_per_sample.csv", low_memory=False)
m12_tests = pd.read_csv(INTER + r"\M11M12_step3_gse106878_arm_tests.csv")
m12_pat = pd.read_csv(INTER + r"\M11M12_step3_gse106878_per_patient.csv")
m12_comp = pd.read_csv(INTER + r"\M11M12_step3_composition_coreport.csv")
m12_h4 = pd.read_csv(INTER + r"\M11M12_step3_h4_attribution.csv")
m12_zero = pd.read_csv(INTER + r"\M11M12_step3_zero_m12.csv")
m12_repl = pd.read_csv(INTER + r"\M11M12_step3_gse148871_replication.csv")
m12_gate = pd.read_csv(INTER + r"\M11M12_step3_gate.csv")
g148_pat = pd.read_csv(INTER + r"\M11M12_step3_gse148871_per_patient.csv")

# ======================================================================
# 1. 图数据 CSV
# ======================================================================
note("== 1. 图数据 CSV")

# Figure 8A：H1 逐队列 + 合并
pool_h1 = dict(pooled=-0.280, ci_lo=-1.159, ci_hi=0.598, p_one=0.1519, I2=58.0,
               pi_lo=-3.143, pi_hi=2.582, k=3, dir_frac=2/3)
fig8a = h1[["cohort", "n_obs", "n_patients", "beta1_std", "se", "p_one_neg", "q_bh",
            "icc", "b_EIS0", "b_comp0", "sex_included", "fit_method", "beta1_ols", "se_ols"]].copy()
fig8a.loc[len(fig8a)] = dict(cohort="POOLED_REML_HK", n_obs=np.nan, n_patients=np.nan,
                             beta1_std=pool_h1["pooled"], se=float(pool_h1["ci_hi"] - pool_h1["pooled"]) / 1.96 if False else np.nan,
                             p_one_neg=pool_h1["p_one"], q_bh=np.nan, icc=np.nan,
                             b_EIS0=np.nan, b_comp0=np.nan, sex_included="", fit_method="pooled",
                             beta1_ols=np.nan, se_ols=np.nan)
fig8a["ci_lo"] = np.nan; fig8a["ci_hi"] = np.nan
fig8a.loc[len(fig8a) - 1, "ci_lo"] = pool_h1["ci_lo"]
fig8a.loc[len(fig8a) - 1, "ci_hi"] = pool_h1["ci_hi"]
fig8a.to_csv(MAIN + r"\Figure_8A.csv", index=False)
note(f"   Figure_8A.csv（H1 逐队列+合并）")

# Figure 8B：H2 逐队列 + 合并
pool_h2 = dict(pooled=-0.137, ci_lo=-0.360, ci_hi=0.086, p_one=0.05916, I2=0.0,
               pi_lo=-0.828, pi_hi=0.554, k=3, dir_frac=1.0)
fig8b = h2[["cohort", "n_obs", "n_pat", "beta_L", "se", "p_one_neg", "q_bh",
            "beta_EISlag_on_UCS", "se_EISlag_on_UCS"]].copy()
fig8b.loc[len(fig8b)] = dict(cohort="POOLED_REML_HK", n_obs=np.nan, n_pat=np.nan,
                             beta_L=pool_h2["pooled"], se=np.nan, p_one_neg=pool_h2["p_one"],
                             q_bh=np.nan, beta_EISlag_on_UCS=np.nan, se_EISlag_on_UCS=np.nan)
fig8b.to_csv(MAIN + r"\Figure_8B.csv", index=False)
note(f"   Figure_8B.csv（H2 逐队列+合并）")

# Figure 8C：零模型 + 敏感性
fig8c = pd.concat([
    z1.assign(test="H1_zero"), z2.assign(test="H2_zero"),
    sens212[["cohort", "n_patients", "beta1_std", "se", "p_one_neg", "note"]].assign(test="GSE212865_sens").rename(columns={"n_patients": "n_obs"}),
    sens545.assign(test="GSE54514_author_neut_sens")[["cohort", "beta1_std", "se", "p_one_neg", "sens"]],
], ignore_index=True)
fig8c.to_csv(MAIN + r"\Figure_8C.csv", index=False)
note(f"   Figure_8C.csv（零模型+敏感性）")

# Figure 8D：LOSO 汇总
fig8d = pd.concat([
    loso1.assign(test="H1_loso_subject"), loso2.assign(test="H2_loso_subject"),
    loso_meta1.assign(test="H1_loso_meta"), loso_meta2.assign(test="H2_loso_meta"),
], ignore_index=True)
fig8d.to_csv(MAIN + r"\Figure_8D.csv", index=False)
note(f"   Figure_8D.csv（LOSO 汇总）")

# Figure 9A：GSE106878 臂间/组内
fig9a = m12_tests.copy()
fig9a.to_csv(MAIN + r"\Figure_9A.csv", index=False)
note(f"   Figure_9A.csv（GSE106878 臂检验）")

# Figure 9B：组成同报 + H4 描述
fig9b = pd.concat([m12_comp.assign(kind="composition_coreport"),
                   m12_h4.assign(kind="H4_descriptive")], ignore_index=True)
fig9b.to_csv(MAIN + r"\Figure_9B.csv", index=False)
note(f"   Figure_9B.csv（组成同报+H4）")

# Figure 9C：GSE148871 复现层
fig9c = pd.concat([m12_repl.assign(kind="summary"), g148_pat.assign(kind="per_patient")], ignore_index=True)
fig9c.to_csv(MAIN + r"\Figure_9C.csv", index=False)
note(f"   Figure_9C.csv（GSE148871 复现层）")

# ======================================================================
# 2. 附表
# ======================================================================
note("\n== 2. 附表 CSV")
per_sample.to_csv(SUP_CSV + r"\Table_S76_M11_PerSample_Longitudinal_Scores.csv", index=False)
note("   Table_S76")
h1.to_csv(SUP_CSV + r"\Table_S77_M11_H1_MixedModels.csv", index=False)
note("   Table_S77")
h2.to_csv(SUP_CSV + r"\Table_S78_M11_H2_RICLPM.csv", index=False)
note("   Table_S78")
pd.concat([z1, z2, loso1, loso2, loso_meta1, loso_meta2,
           sens212.assign(kind="gse212865_sens"), sens545.assign(kind="gse54514_author_neut")],
          ignore_index=True).to_csv(SUP_CSV + r"\Table_S79_M11_NullModels_LOSO_Sensitivity.csv", index=False)
note("   Table_S79")
m12_pat.to_csv(SUP_CSV + r"\Table_S80_M12_GSE106878_PerPatient.csv", index=False)
pd.concat([m12_tests, m12_comp, m12_h4, m12_zero], ignore_index=True).to_csv(
    SUP_CSV + r"\Table_S81_M12_GSE106878_Tests_Comp_H4_Null.csv", index=False)
note("   Table_S81")
g148_pat.to_csv(SUP_CSV + r"\Table_S82_M12_GSE148871_Replication.csv", index=False)
note("   Table_S82")

# ======================================================================
# 3. Manifest 登记
# ======================================================================
note("\n== 3. Manifest 登记")
new_lines = [
    "Figure_8A,Main/Figure_8A.csv,_intermediate/M11M12_step2_h1_cohort.csv,时间序（M11，前瞻注册）,H1 主检验1逐队列β1_std（基线UCS→随访EIS，患者随机截距混合模型）+REML/HK合并,森林图,单侧β1<0；GSE215865 −0.344/GSE54514 −0.734/GSE148871 +0.019；合并−0.280[−1.159,+0.598]单侧p=0.152；BH族6检验",
    "Figure_8B,Main/Figure_8B.csv,_intermediate/M11M12_step2_h2_cohort.csv,时间序（M11，前瞻注册）,H2 主检验2 RI-CLPM滞后交叉路径β_L（UCS(t−1)→EIS(t)）逐队列+合并,森林图,单侧βL<0；3/3负向；合并−0.137[−0.360,+0.086]单侧p=0.059；I²=0%",
    "Figure_8C,Main/Figure_8C.csv,_intermediate/M11M12_step2_zero_h1.csv,时间序（M11，前瞻注册）,零模型（表达量匹配随机集B=10000）与敏感性（GSE212865两波、GSE54514作者中性粒）,直方图/点图,seed=0；H1/H2逐队列obs vs null",
    "Figure_8D,Main/Figure_8D.csv,_intermediate/M11M12_step2_h1_loso_subject.csv,时间序（M11，前瞻注册）,LOSO方向一致性（队列内leave-one-subject-out + meta leave-one-cohort-out）,条形图,H2 LOSO 100%/97%/100%；H1 GSE148871 81%",
    "Figure_9A,Main/Figure_9A.csv,_intermediate/M11M12_step3_gse106878_arm_tests.csv,干预响应（M12，前瞻注册）,GSE106878 d_res臂间MWU（单侧less）与组内Wilcoxon（单侧）+BH族,点图/条形图,臂间p=0.1231（q=0.1231阴性）；组内HC p=0.00032（q=0.00065）；P(d_res<0)=79%",
    "Figure_9B,Main/Figure_9B.csv,_intermediate/M11M12_step3_composition_coreport.csv,干预响应（M12，前瞻注册）,强制组成同报（ΔNeut/ΔMono逐臂）+H4描述（ΔMDI_obs vs ΔMDI_comp）,双面板条形图,确认性未显著→H4仅描述；HC Δcomp与Δobs反向",
    "Figure_9C,Main/Figure_9C.csv,_intermediate/M11M12_step3_gse148871_replication.csv,干预响应（M12，前瞻注册）,GSE148871 NEMI臂复现层：d_obs锚点（−0.489,p=0.0045）vs d_res（+0.031,p=0.707）与归因,点图/条形图,同向份额+2.184→组成回退（披露的注册前观察，不进BH族）",
]
with open(MANIFEST_F, "a", encoding="utf-8") as f:
    f.write("\n" + "\n".join(new_lines) + "\n")
note(f"   FIGURE_DATA_MANIFEST.csv 追加 {len(new_lines)} 行")

# ======================================================================
# 4. 审计报告
# ======================================================================
note("\n== 4. 审计报告")
m12_v = m12_gate["verdict"].iloc[0]
gate_row = gate[gate["test"] == "verdict"]
m11_v = gate_row["verdict"].iloc[0] if len(gate_row) else "见 gate 表"
rep = []
rep.append("# M11/M12 时间序与干预响应 分析报告（2026-08-27）\n")
rep.append("- 前瞻注册：M11_M12_pre_registration_20260827_draft.md（第二次独立注册：osf.io/98cm3；与 osf.io/ETVMJ 并列）")
rep.append("- 数据：GSE215865（1,326 day 级样本/508 受试者）、GSE54514（163/54 ID）、GSE148871（血 168/49 例）、GSE212865（137）、GSE106878（94=47×2）")
rep.append("- 脚本：M11M12_step1_scores_deconv.py / M11M12_step2_temporal.py / M11M12_step3_intervention.py / M11M12_step4_export.py；日志 03_LOGS/M11M12_step{1,2,3}_log.txt\n")
rep.append("## 1. M11 时间序（判据 ii）\n")
rep.append("### 1.1 H1 主检验 1（患者随机截距混合模型：随访 EIS ~ 基线 UCS + 基线 EIS + 基线组成（中性粒+单核）+ 性别；单侧 β1<0）\n")
for _, r in h1.iterrows():
    rep.append(f"- {r['cohort']}: n_obs={int(r['n_obs'])} n_pat={int(r['n_patients'])} "
               f"β1_std={r['beta1_std']:+.3f} (SE {r['se']:.3f}) 单侧p={r['p_one_neg']:.4g} "
               f"q={r['q_bh']:.4g} ICC={r['icc']:.3f}（OLS 同报 β={r.get('beta1_ols', np.nan):+.3f}）")
rep.append(f"- 合并（REML+Hartung-Knapp）: β1_std=−0.280 [−1.159,+0.598]，单侧 p=0.152，I²=58%，预测区间 [−3.14,+2.58]")
rep.append(f"- 方向一致 2/3（67%）——低于预注册 ≥70% 门；meta LOSO 3/3；队列内 LOSO 358/358、49/49、35/43（81%）")
rep.append(f"- GSE54514 作者实测中性粒比例替换组成协变量敏感性: β1_std=−0.077（SE 0.141，单侧 p=0.293）——该队列显著效应对组成协变量口径敏感（如实披露）")
rep.append("### 1.2 H2 主检验 2（RI-CLPM 交叉滞后：EIS(t) ~ φE·EISw(t−1) + βL·UCSw(t−1) + γ·compw(t)；单侧 βL<0）\n")
for _, r in h2.iterrows():
    rep.append(f"- {r['cohort']}: n_obs={int(r['n_obs'])} n_pat={int(r['n_pat'])} "
               f"β_L={r['beta_L']:+.3f} (SE {r['se']:.3f}) 单侧p={r['p_one_neg']:.4g} q={r['q_bh']:.4g}")
rep.append(f"- 合并: β_L=−0.137 [−0.360,+0.086]，单侧 p=0.059，I²=0%，预测区间 [−0.828,+0.554]")
rep.append(f"- 方向一致 3/3（100%，过 ≥70% 门）；meta LOSO 3/3；队列内 LOSO 236/236、28/29、41/41（≥90%）")
rep.append(f"- BH 族（6 检验）内 GSE215865 q=0.0255 显著")
rep.append(f"- 零模型（B=10,000，seed=0）：" + "；".join(
    f"{r['cohort']} {r['hypothesis']} 经验p={r['empirical_p']:.4g}" for _, r in pd.concat([z1, z2]).iterrows()))
rep.append(f"- GSE212865 两波敏感性: β1_std={sens212['beta1_std'].iloc[0]:+.3f}（SE {sens212['se'].iloc[0]:.3f}，单侧 p={sens212['p_one_neg'].iloc[0]:.4g}；不进 BH 族）")
rep.append(f"\n**M11 判定: {m11_v}**（H2 三门全过：方向 3/3 + BH 内显著 + LOSO ≥90%；H1 方向一致 2/3=67% 未达 70% 门）\n")
rep.append("## 2. M12 干预响应（判据 iii）\n")
rep.append("### 2.1 GSE106878 确认性（d_res = 组成残差 MDI 的 Post−Pre；臂间单侧 Mann-Whitney HC<PL 为主检验）\n")
rep.append(f"- 臂间 MWU: U={m12_tests.iloc[0]['u']:.1f} z={m12_tests.iloc[0]['z']:+.3f} p={m12_tests.iloc[0]['p']:.4g} "
           f"q={m12_tests.iloc[0]['q_bh']:.4g}（**未过 BH**）Cliff's δ={m12_tests.iloc[0]['cliff_delta']:+.3f}")
rep.append(f"- 组内 HC Wilcoxon（单侧）: p={m12_tests.iloc[1]['p']:.4g} q={m12_tests.iloc[1]['q_bh']:.4g}；"
           f"P(d_res<0)={m12_tests.iloc[1]['p_decline']:.0%}（MESI 门 ≥65% 达标）；中位 d_res=−0.239")
rep.append(f"- PL 组内（描述，双侧）: 中位 −0.165，p={m12_tests.iloc[2]['p_two']:.4g}——安慰剂臂同样自发回落（臂间检验即为此控制）")
rep.append(f"- LOSO 方向一致: 臂间 100%、组内 100%；零模型: obs z=−1.170 vs null +0.001±1.076，经验 p=0.1405")
rep.append("- H4 归因: 确认性检验未显著 → 仅描述报告（HC 臂 mean ΔMDI_obs=−0.193、mean ΔMDI_comp=+0.371，反向）")
rep.append("### 2.2 GSE148871 复现层（披露的注册前观察，不进 BH 族）\n")
r_ = m12_repl.iloc[0]
rep.append(f"- NEMI 臂 21 对 SCREENING→V4-D28: d_obs 中位 −0.489（Wilcoxon 单侧 p=0.0045，锚点复现注册前披露 ΔMDI≈−0.33/p=0.016）；"
           f"**d_res 中位 +0.031（p=0.707，组成残差口径不下降）**；d_Neut 中位 −0.047（p=0.0016）")
rep.append(f"- H4 归因: 同向份额 +2.184（≥0.5 且同向）→ **组成回退**——观测 MDI 回落被组成拟合变化完全解释")
rep.append(f"\n**M12 判定: {m12_v}**（确认性臂间检验未过 BH → 判据 iii 如实阴性入稿；复现层同向复现失败且判为组成回退）\n")
rep.append("## 3. 口径偏离与披露\n")
rep.append("1. ABIS Micro NNLS 在 Illumina quantile+log2 数据上退化解（单核比例恒 0、中性粒 5–21%，与全血预期血象不符；GSE54514 有作者实测中性粒比例可交叉）——按预注册跨平台敏感性口径执行并如实披露；H1/H2 组成协变量用合并（中性粒+单核）比例（与 OSF 表单 H1 β₃·cell-composition 单一组成项一致），避免单类型共线退化。")
rep.append("2. 组成预测 MDI（MDI_comp）沿用 M10A 主口径（Monaco 29 型 NNLS；ABIS 臂覆盖 0–1/63 不足以构造合成评分）。")
rep.append("3. 零模型统计量采用 OLS 快速口径（观测与零同口径；主推断仍为 MixedLM/聚类稳健 SE）——如实披露。")
rep.append("4. GSE215865 无性别字段（Synapse 元数据未入列），H1/H2 该队列不含性别协变量（注册条款\"可得处\"）。")
rep.append("5. GSE215865 混合模型 ICC≈0.027 边界收敛，OLS+聚类稳健 SE 同报；RI-CLPM 波次视为等间隔（T0–T13/GSE148871 访视间隔不均，近似处理披露）。")
rep.append("6. **共享库 mdi_lib.bh 修复（2026-08-27）**：原实现当最小 p 不在首位时 q 值错位；本模块全部 BH 用修复版。历史调用该函数的 M2_step2_replication_meta.py 与 M14_step2_perturbseq_gse221321.py 的 BH_q 列为**待作者复核的历史项**（本次未改动任何冻结输出）。")
rep.append("7. GSE106878 非归一化矩阵含 187 个 sentrix 列，仅 94 个与 GEO 官方样本（47×2 配对）匹配——按注册条款仅用 94 个官方样本，多余列弃用留痕。")
rep.append("\n## 4. 结论\n")
rep.append(f"- M11（判据 ii）: {m11_v}——交叉滞后（UCS(t−1)→EIS(t) 负向）三门全过；基线→随访方向一致未达 70% 门。")
rep.append(f"- M12（判据 iii）: {m12_v}——GSE106878 组间未过 BH；GSE148871 观测 MDI 回落为组成回退。")
rep.append("- 以上均按注册条款如实入稿，无终点切换。")
with open(GOV + r"\M11M12_Analysis_Report_20260827.md", "w", encoding="utf-8") as f:
    f.write("\n".join(rep))
note(f"   {GOV}\\M11M12_Analysis_Report_20260827.md")

with open(r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\03_LOGS\M11M12_step4_log.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(log))
note(f"\nDONE in {time.time()-t0:.1f}s")
