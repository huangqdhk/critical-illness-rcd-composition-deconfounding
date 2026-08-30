# -*- coding: utf-8 -*-
"""
M2_step2_replication_meta.py — M2 第 2 步：包内回算组间对比 + 跨队列复制 + meta 合成
==============================================================================
预注册：M2_pre_registration_20260817.md
- 每队列主对比（疾病 vs 对照）：Mann-Whitney U + Cliff's delta + Hedges' g + BH
- 无对照队列（疾病-only）：单样本检验（MDI 均值 vs 0）+ EIS−UCS 符号
- meta：DerSimonian-Laird 随机效应（有对照队列的 Hedges' g）+ Stouffer z
- 自洽性锚点：GSE185263 MDI 方向必须与 Table S2（上游 −0.83/执行 +0.18）一致
输出：_intermediate/M2_contrasts.csv, M2_meta.csv
"""
import os, io, sys
import numpy as np
import pandas as pd
sys.path.insert(0, r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS")
import mdi_lib as L

ROOT = L.ROOT
INTER = ROOT + r"\_intermediate"
df = pd.read_csv(INTER + r"\M2_per_sample_scores.csv")
log = []
def note(m):
    log.append(m); print(m, flush=True)

# ---------------- 主对比定义（预注册：每队列 1 个主疾病对比 + 关键次对比）----------------
CONTRASTS = [
    # (cohort, case_groups, control_groups, contrast_id)
    ("GSE185263", ["Sepsis_COVID"], ["Control"], "SepsisCOVID_vs_Control"),
    ("GSE185263", ["Sepsis"], ["Control"], "Sepsis_vs_Control"),
    ("GSE32707",  ["ARDS_d0"],    ["Control"], "ARDSd0_vs_Control"),
    ("GSE32707",  ["Sepsis_d0"],  ["Control"], "Sepsisd0_vs_Control"),
    ("GSE212865", ["Covid19_SDRA"], ["Control"], "SDRA_vs_Control"),
    ("GSE212865", ["Covid19"],      ["Control"], "COVID_vs_Control"),
    ("GSE310929", ["Sepsis", "Sepsis - Shock", "Sepsis - ARDs"], ["Control"], "Sepsis_vs_Control"),
]

rows = []
for cohort, cases, ctrls, cid in CONTRASTS:
    sub = df[df["cohort"] == cohort]
    case = sub[sub["group"].isin(cases)]
    ctrl = sub[sub["group"].isin(ctrls)]
    if len(case) == 0 or len(ctrl) == 0:
        note(f"  {cid}: 空组 (case={len(case)}, ctrl={len(ctrl)})，跳过")
        continue
    for score in ("UCS", "EIS", "MDI"):
        p, delta = L.mwu_cliff(ctrl[score].values, case[score].values)
        g, gse = L.hedges_g(ctrl[score].values, case[score].values)
        rows.append(dict(cohort=cohort, contrast=cid, score=score,
                         n_case=len(case), n_ctrl=len(ctrl),
                         case_mean=case[score].mean(), ctrl_mean=ctrl[score].mean(),
                         diff=case[score].mean() - ctrl[score].mean(),
                         cliff_delta=delta, hedges_g=g, hedges_g_se=gse,
                         MW_p=p))
        note(f"  {cid} {score}: diff={case[score].mean()-ctrl[score].mean():+.3f} g={g:+.3f} p={p:.3g}")

ct = pd.DataFrame(rows)
# BH within each contrast-score family (per cohort, per contrast → 3 tests)
ct["BH_q"] = np.nan
for (cohort, cid), grp in ct.groupby(["cohort", "contrast"]):
    ct.loc[grp.index, "BH_q"] = L.bh(grp["MW_p"].values)
ct["gene_set_version"] = L.GENE_SET_VERSION
ct["score_version"] = L.SCORE_VERSION
ct.to_csv(INTER + r"\M2_contrasts.csv", index=False)

# ---------------- 自洽性锚点检查（预注册 §1 末条）----------------
note("\n== 自洽性锚点检查")
g185 = df[df["cohort"] == "GSE185263"]
sc = g185[g185["group"] == "Sepsis_COVID"]; co = g185[g185["group"] == "Control"]
note(f"  GSE185263 Sepsis_COVID vs Control: UCS Δ={sc['UCS'].mean()-co['UCS'].mean():+.3f} (期望<0), "
     f"EIS Δ={sc['EIS'].mean()-co['EIS'].mean():+.3f} (期望>0), MDI Δ={sc['MDI'].mean()-co['MDI'].mean():+.3f} (期望>0)")
anchor_ok = (sc["UCS"].mean() < co["UCS"].mean()) and (sc["EIS"].mean() > co["EIS"].mean())
note(f"  锚点检查: {'PASS' if anchor_ok else 'FAIL'}（Table S2: 上游臂 log2FC=-0.83<0、执行臂 +0.18>0）")

# ---------------- meta 合成（有对照队列的 Hedges' g）----------------
note("\n== 随机效应 meta（DerSimonian-Laird）")
meta_rows = []
for score in ("UCS", "EIS", "MDI"):
    sub = ct[(ct["score"] == score)]
    # 每队列取主对比：GSE185263 取 SepsisCOVID、GSE32707 取 ARDSd0、GSE310929 取 Sepsis
    keep = sub[sub["contrast"].isin(["SepsisCOVID_vs_Control", "ARDSd0_vs_Control", "Sepsis_vs_Control"])]
    keep = keep.drop_duplicates(subset=["cohort"], keep="first")
    if len(keep) >= 2:
        g, se, p, tau2, I2, Q, dfq, Qp = L.der_simonian_laird(keep["hedges_g"].values, keep["hedges_g_se"].values)
        z, pz = L.stouffer_z(keep["MW_p"].values)
        meta_rows.append(dict(score=score, n_cohorts=len(keep),
                              pooled_g=g, pooled_se=se, meta_p=p, tau2=tau2, I2=I2,
                              Q=Q, Q_df=dfq, Q_p=Qp, stouffer_z=z, stouffer_p=pz,
                              cohorts=";".join(f"{r.cohort}({r.contrast})" for r in keep.itertuples())))
        note(f"  {score}: pooled g={g:+.3f} p={p:.3g} I2={I2:.2f} 队列={list(keep['cohort'])}")

mt = pd.DataFrame(meta_rows)
mt["gene_set_version"] = L.GENE_SET_VERSION
mt["score_version"] = L.SCORE_VERSION
mt.to_csv(INTER + r"\M2_meta.csv", index=False)

# ---------------- 疾病-only 队列的单样本口径（GSE310929 脓毒症亚组 MDI vs 0）----------------
note("\n== 脓毒症亚组 MDI 单样本口径（z 空间，期望>0）")
sep = df[(df["cohort"] == "GSE310929") & (df["Disease Simplified"].isin(["Sepsis", "Sepsis - Shock", "Sepsis - ARDs"]))]
from scipy import stats as st
note(f"  脓毒症 n={len(sep)}: MDI mean={sep['MDI'].mean():+.3f}, 95%CI=[{sep['MDI'].mean()-1.96*sep['MDI'].std()/np.sqrt(len(sep)):+.3f}, {sep['MDI'].mean()+1.96*sep['MDI'].std()/np.sqrt(len(sep)):+.3f}], t={st.ttest_1samp(sep['MDI'], 0).statistic:+.2f}, p={st.ttest_1samp(sep['MDI'], 0).pvalue:.3g}")
note(f"  对照 n={len(df[(df['cohort']=='GSE310929')&(df['Disease Simplified']=='Control')])}")

# ---------------- GSE188309 CAP（疾病-only 单样本口径）----------------
note("\n== GSE188309 CAP MDI 单样本口径（z 空间，期望>0）")
cap = df[df["cohort"] == "GSE188309"]
if len(cap):
    note(f"  CAP n={len(cap)}: MDI mean={cap['MDI'].mean():+.3f}, UCS={cap['UCS'].mean():+.3f}, EIS={cap['EIS'].mean():+.3f}, "
         f"t={st.ttest_1samp(cap['MDI'], 0).statistic:+.2f}, p={st.ttest_1samp(cap['MDI'], 0).pvalue:.3g}")

# ---------------- GSE148871 COPD（血 vs 痰；治疗前后）----------------
note("\n== GSE148871 COPD 组织与治疗对比")
copd = df[df["cohort"] == "GSE148871"]
if len(copd):
    blood = copd[copd["tissue"] == "whole blood"]
    sput = copd[copd["tissue"] == "sputum"]
    for score in ("UCS", "EIS", "MDI"):
        p, delta = L.mwu_cliff(sput[score].values, blood[score].values)
        note(f"  Blood(n={len(blood)}) vs Sputum(n={len(sput)}) {score}: blood mean={blood[score].mean():+.3f}, sputum mean={sput[score].mean():+.3f}, delta={delta:+.3f}, p={p:.3g}")
    screen = copd[copd["visit"] == "SCREENING"]
    v4 = copd[copd["visit"] == "V4 DAY 28"]
    for score in ("UCS", "EIS", "MDI"):
        if len(screen) and len(v4):
            p, delta = L.mwu_cliff(screen[score].values, v4[score].values)
            note(f"  SCREENING(n={len(screen)}) vs V4D28(n={len(v4)}) {score}: Δ={v4[score].mean()-screen[score].mean():+.3f}, p={p:.3g}")

with open(INTER + r"\M2_step2_log.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(log))
print("\nDONE step2")
