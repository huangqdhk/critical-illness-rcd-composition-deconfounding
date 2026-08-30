# -*- coding: utf-8 -*-
"""
M10D_step4_pandisease_k20.py — M10 注册次终点 §1.5 剩余项：k≥20 泛疾病扩展
============================================================================
前瞻注册：M10_M13_M14_pre_registration_20260824.md（osf.io/ETVMJ）§1.5——
"≥20 个来源互斥、带对照的急性炎症队列（GSE310929 来源拆分 + 已持有队列 + C/D 类候选），
MDI 同口径回算 + 随机效应合并。"

【口径冻结——先于任何数值，写入日志】
- 层定义（来源互斥、带对照、急性炎症、全血 bulk）：
  * 冻结 16 层（M10A 检验 A 层口径，obs MDI 效应取自冻结表 M10A_residual_layer_effects.csv）：
    13 个 GSE310929 来源拆分子层（case∈{Sepsis, Sepsis-Shock, Sepsis-ARDs} vs Control）
    + GSE185263（Sepsis_COVID vs Control）+ GSE32707（ARDS_d0 vs Control）
    + GSE212865_base（Covid19_SDRA vs Control，D0）。
  * 新增 4 层（全部已入盘、已经 GEO 官方记录核验红线 #1、mdi_v1.0 同口径逐样本评分
    为既有冻结表，本脚本不重评分，只做层效应 + 合并）：
    - GSE66099（M15 采纳演示 2；儿童全血 GPL570）：case = Sepsis+SepticShock vs Control；
      来源：Table_S86_M15_Demo_GSE66099_PerSample.csv（冻结评分）。
    - GSE54514（M11 主队列 2；成人脓毒症 PAXgene 全血 Illumina）：case = 脓毒症存活+死亡
      的 Day_1 基线 vs 健康对照 Day_1；来源：M11M12_step1_per_sample.csv（冻结评分；
      GPL6947 平台臂覆盖 12/30+14/33，如实披露）。
    - GSE157103（M15 采纳演示 1；COVID-19 全血 RNA-seq）：case = COVID-19 vs non-COVID-19；
      GEO 官方设计 "26 non-COVID-19" 为**住院非 COVID 患者**（非健康对照，如实披露）；
      来源：Table_S85_M15_Demo_GSE157103_PerSample.csv（冻结评分）。
    - GSE232404（GSE310929 来源拆分子层，此前未被 M10A 使用）：case = "Sepsis - AKI" vs
      Control（n=5v5）；case 标签超出 M10A 冻结三标签集，作为层集扩展如实披露，且以
      不含本层的敏感性（k=19/k=18）交叉。
- 层效应：Hedges' g（mdi_lib.hedges_g，正值=case 高）+ Mann-Whitney p；合并 = REML τ² +
  Hartung-Knapp + 95% CI + 预测区间 + I²/Q + LOSO（与 M10A pool_hk 同款）。
- 主分析 k=20；敏感性 k=19（去 GSE232404）、k=18（去 GSE232404+GSE157103，仅健康对照
  新增层）、k=16（冻结锚点重算核对）。
- 本层只报告观测 MDI 同口径信号（注册原文"MDI 同口径回算"）；组成校正层不属于本终点。
输出：
  _intermediate/M10D_pandisease_k20_layers.csv
  _intermediate/M10D_pandisease_k20_meta.csv
  _intermediate/M10D_pandisease_k20_loso.csv
  02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV/Table_S87_M10D_Pandisease_k20_LayerEffects.csv
  02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV/Table_S87b_M10D_Pandisease_k20_Meta.csv
  Table_S87c_M10D_Pandisease_k20_LOSO.csv
日志：03_LOGS/M10D_pandisease_k20_log.txt
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import mdi_lib as L

ROOT = os.path.dirname(os.path.abspath(__file__))
INTER = os.path.join(ROOT, "_intermediate")
TAB = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES", "SUPPLEMENTARY_Tables_CSV")
LOGP = os.path.join(ROOT, "03_LOGS", "M10D_pandisease_k20_log.txt")

log = []
def note(m=""):
    log.append(m); print(m, flush=True)

def flush_log():
    with open(LOGP, "w", encoding="utf-8") as f:
        f.write("\n".join(log) + "\n")

# ---------------- meta 数学（P0/M10A 同款 REML + Hartung-Knapp） ----------------
def reml_tau2(g, se):
    g = np.asarray(g, float); v = np.asarray(se, float) ** 2
    def F(t):
        w = 1.0 / (v + t)
        mu = np.sum(w * g) / np.sum(w)
        return np.sum(w ** 2 * ((g - mu) ** 2 - v - t))
    if F(0.0) <= 0:
        return 0.0
    hi = 1.0
    while F(hi) > 0 and hi < 1e9:
        hi *= 4.0
    try:
        return float(brentq(F, 0.0, hi, xtol=1e-12))
    except Exception:
        return 0.0

def pool_hk(g, se):
    g = np.asarray(g, float); se = np.asarray(se, float); k = len(g)
    t2 = reml_tau2(g, se)
    w = 1.0 / (se ** 2 + t2)
    mu = float(np.sum(w * g) / np.sum(w))
    se_mu = float(np.sqrt(1.0 / np.sum(w)))
    se_hk = float(np.sqrt(np.sum(w * (g - mu) ** 2) / ((k - 1) * np.sum(w))))
    tstat = mu / se_hk
    wf = 1.0 / se ** 2
    mu_fe = float(np.sum(wf * g) / np.sum(wf))
    Q = float(np.sum(wf * (g - mu_fe) ** 2)); dfq = k - 1
    I2 = max(0.0, (Q - dfq) / Q) * 100 if Q > 0 else 0.0
    if k > 2:
        tp = float(stats.t.ppf(0.975, k - 2))
        half = tp * np.sqrt(t2 + se_mu ** 2)
        pi_lo, pi_hi = mu - half, mu + half
    else:
        pi_lo = pi_hi = np.nan
    return dict(k=k, pooled_g=mu, se=se_mu, se_hk=se_hk,
                ci_lo=mu - float(stats.t.ppf(0.975, k - 1)) * se_hk,
                ci_hi=mu + float(stats.t.ppf(0.975, k - 1)) * se_hk,
                p_pos=float(1 - stats.t.cdf(tstat, k - 1)),
                p_neg=float(stats.t.cdf(tstat, k - 1)),
                tau2=t2, I2=I2, Q=Q, Q_p=float(1 - stats.chi2.cdf(Q, dfq)),
                pi_lo=float(pi_lo), pi_hi=float(pi_hi),
                dir_pos=int((g > 0).sum()), dir_neg=int((g < 0).sum()))

def layer_effect(case_mdi, ctrl_mdi):
    g, gse = L.hedges_g(np.asarray(ctrl_mdi, float), np.asarray(case_mdi, float))
    p_mw = stats.mannwhitneyu(case_mdi, ctrl_mdi, alternative="two-sided").pvalue
    return g, gse, p_mw

def main():
    note("== M10D 泛疾病扩展（注册次终点 ETVMJ §1.5 剩余项）2026-08-27")
    note("口径冻结（先于任何数值，见脚本 docstring）：层集 = 冻结 16 层 obs MDI + 4 个新增"
         "入盘核验层（GSE66099 / GSE54514 / GSE157103 / GSE232404），层效应 Hedges' g，"
         "REML+Hartung-Knapp+PI+LOSO；主分析 k=20，敏感性 k=19/k=18，k=16 锚点重算核对。")

    # ---------------- 1. 冻结 16 层 ----------------
    frozen = pd.read_csv(os.path.join(INTER, "M10A_residual_layer_effects.csv"))
    frozen = frozen[["stratum", "n_case", "n_ctrl", "obs_g", "obs_g_se"]].copy()
    frozen["source"] = "frozen16"
    frozen["n_up"], frozen["n_ex"], frozen["MW_p"] = np.nan, np.nan, np.nan
    note(f"   冻结 16 层 obs MDI 载入（M10A 冻结表）；锚点合并 g={frozen['obs_g'].mean():+.3f}（算术，非合并口径）")

    # ---------------- 2. 新增层（复用既有冻结评分表，不重评分） ----------------
    rows = []

    # 2a GSE66099（M15 演示冻结表；case = Sepsis+SepticShock vs Control）
    t86 = pd.read_csv(os.path.join(TAB, "Table_S86_M15_Demo_GSE66099_PerSample.csv"))
    case = t86[t86["group"].isin(["Sepsis", "SepticShock"])]
    ctrl = t86[t86["group"] == "Control"]
    g, gse, pmw = layer_effect(case["MDI"], ctrl["MDI"])
    rows.append(dict(stratum="GSE66099", source="M15_demo_verified", n_case=len(case), n_ctrl=len(ctrl),
                     obs_g=g, obs_g_se=gse, MW_p=pmw,
                     n_up=int(case["n_up"].iloc[0]), n_ex=int(case["n_ex"].iloc[0]),
                     UCS_g=L.hedges_g(np.asarray(ctrl["UCS"], float), np.asarray(case["UCS"], float))[0],
                     EIS_g=L.hedges_g(np.asarray(ctrl["EIS"], float), np.asarray(case["EIS"], float))[0]))
    note(f"   GSE66099: Sepsis+SepticShock {len(case)} vs Control {len(ctrl)} | "
         f"MDI g={g:+.3f} (se={gse:.3f}) MW p={pmw:.3g} | n_up={rows[-1]['n_up']} n_ex={rows[-1]['n_ex']}")

    # 2b GSE54514（M11 冻结评分；case = 脓毒症（存活+死亡）Day_1 vs 健康对照 Day_1）
    m11 = pd.read_csv(os.path.join(INTER, "M11M12_step1_per_sample.csv"))
    g54 = m11[(m11["cohort"] == "GSE54514") & (m11["day_title"] == "Day_1")].copy()
    case = g54[g54["grp_title"].isin(["sepsis_survivor", "sepsis_nonsurvivor"])]
    ctrl = g54[g54["grp_title"] == "Control"]
    g, gse, pmw = layer_effect(case["MDI"], ctrl["MDI"])
    rows.append(dict(stratum="GSE54514", source="M11_verified", n_case=len(case), n_ctrl=len(ctrl),
                     obs_g=g, obs_g_se=gse, MW_p=pmw,
                     n_up=int(case["n_up"].iloc[0]), n_ex=int(case["n_ex"].iloc[0]),
                     UCS_g=L.hedges_g(np.asarray(ctrl["UCS"], float), np.asarray(case["UCS"], float))[0],
                     EIS_g=L.hedges_g(np.asarray(ctrl["EIS"], float), np.asarray(case["EIS"], float))[0]))
    note(f"   GSE54514: 脓毒症D1 {len(case)} vs 对照D1 {len(ctrl)} | "
         f"MDI g={g:+.3f} (se={gse:.3f}) MW p={pmw:.3g} | n_up={rows[-1]['n_up']} n_ex={rows[-1]['n_ex']}（GPL6947 覆盖受限，披露）")

    # 2c GSE157103（M15 演示冻结表；COVID vs non-COVID，对照=住院非 COVID 患者，披露）
    t85 = pd.read_csv(os.path.join(TAB, "Table_S85_M15_Demo_GSE157103_PerSample.csv"))
    case = t85[t85["group"] == "COVID"]
    ctrl = t85[t85["group"] == "non-COVID"]
    g, gse, pmw = layer_effect(case["MDI"], ctrl["MDI"])
    rows.append(dict(stratum="GSE157103", source="M15_demo_verified", n_case=len(case), n_ctrl=len(ctrl),
                     obs_g=g, obs_g_se=gse, MW_p=pmw,
                     n_up=int(case["n_up"].iloc[0]), n_ex=int(case["n_ex"].iloc[0]),
                     UCS_g=L.hedges_g(np.asarray(ctrl["UCS"], float), np.asarray(case["UCS"], float))[0],
                     EIS_g=L.hedges_g(np.asarray(ctrl["EIS"], float), np.asarray(case["EIS"], float))[0]))
    note(f"   GSE157103: COVID {len(case)} vs non-COVID {len(ctrl)} | "
         f"MDI g={g:+.3f} (se={gse:.3f}) MW p={pmw:.3g}（对照=住院非 COVID 患者，披露）")

    # 2d GSE232404（GSE310929 拆分子层；case="Sepsis - AKI" 5v5，标签扩展披露）
    m2 = pd.read_csv(os.path.join(INTER, "M2_per_sample_scores.csv"))
    g23 = m2[(m2["cohort"] == "GSE310929") & (m2["Dataset"] == "GSE232404")]
    case = g23[g23["group"] == "Sepsis - AKI"]
    ctrl = g23[g23["group"] == "Control"]
    g, gse, pmw = layer_effect(case["MDI"], ctrl["MDI"])
    rows.append(dict(stratum="310:GSE232404", source="GSE310929_subsource_extended", n_case=len(case), n_ctrl=len(ctrl),
                     obs_g=g, obs_g_se=gse, MW_p=pmw,
                     n_up=int(case["n_up"].iloc[0]), n_ex=int(case["n_ex"].iloc[0]),
                     UCS_g=L.hedges_g(np.asarray(ctrl["UCS"], float), np.asarray(case["UCS"], float))[0],
                     EIS_g=L.hedges_g(np.asarray(ctrl["EIS"], float), np.asarray(case["EIS"], float))[0]))
    note(f"   310:GSE232404: Sepsis-AKI {len(case)} vs Control {len(ctrl)} | "
         f"MDI g={g:+.3f} (se={gse:.3f}) MW p={pmw:.3g}（case 标签超出冻结三标签集 + n=5v5，敏感性交叉）")

    new = pd.DataFrame(rows)
    layers = pd.concat([frozen, new], ignore_index=True)
    note(f"\n   层集合计 {len(layers)}（冻结 16 + 新增 4）")

    # ---------------- 3. 合并（主 k=20 + 敏感性） ----------------
    def run_meta(sub, tag):
        m = pool_hk(sub["obs_g"].values, sub["obs_g_se"].values)
        note(f"   [{tag}] k={m['k']} pooled MDI g={m['pooled_g']:+.3f} "
             f"95%CI=[{m['ci_lo']:+.3f},{m['ci_hi']:+.3f}] 单侧p(pos)={m['p_pos']:.4g} "
             f"I2={m['I2']:.0f}% PI=[{m['pi_lo']:+.3f},{m['pi_hi']:+.3f}] "
             f"方向 {m['dir_pos']}/{m['k']} 为正")
        m["tag"] = tag
        return m

    meta_rows = []
    meta_rows.append(run_meta(layers, "k20_main"))
    meta_rows.append(run_meta(layers[layers["stratum"] != "310:GSE232404"], "k19_sens_noAKI"))
    meta_rows.append(run_meta(layers[~layers["stratum"].isin(["310:GSE232404", "GSE157103"])], "k18_sens_healthyctrl"))
    meta_rows.append(run_meta(layers[layers["source"] == "frozen16"], "k16_anchor_recompute"))
    meta = pd.DataFrame(meta_rows)

    # LOSO（主 k=20）
    loso_rows = []
    full_sign = np.sign(pool_hk(layers["obs_g"].values, layers["obs_g_se"].values)["pooled_g"])
    for i in range(len(layers)):
        keep = layers.drop(layers.index[i])
        m = pool_hk(keep["obs_g"].values, keep["obs_g_se"].values)
        loso_rows.append(dict(left_out=layers.iloc[i]["stratum"], k=m["k"], pooled_g=m["pooled_g"],
                              same_sign=bool(np.sign(m["pooled_g"]) == full_sign)))
    loso = pd.DataFrame(loso_rows)
    note(f"   LOSO（k=20）同向 {int(loso['same_sign'].sum())}/{len(loso)}")

    # ---------------- 4. 落盘 ----------------
    layers_out = layers[["stratum", "source", "n_case", "n_ctrl", "n_up", "n_ex",
                         "obs_g", "obs_g_se", "MW_p", "UCS_g", "EIS_g"]].rename(
        columns={"obs_g": "MDI_g", "obs_g_se": "MDI_g_se"})
    layers_out.to_csv(os.path.join(INTER, "M10D_pandisease_k20_layers.csv"), index=False)
    layers_out.to_csv(os.path.join(TAB, "Table_S87_M10D_Pandisease_k20_LayerEffects.csv"), index=False)
    meta_out = meta[["tag", "k", "pooled_g", "se", "se_hk", "ci_lo", "ci_hi", "p_pos", "p_neg",
                     "tau2", "I2", "Q", "Q_p", "pi_lo", "pi_hi", "dir_pos", "dir_neg"]]
    meta_out.to_csv(os.path.join(INTER, "M10D_pandisease_k20_meta.csv"), index=False)
    meta_out.to_csv(os.path.join(TAB, "Table_S87b_M10D_Pandisease_k20_Meta.csv"), index=False)
    loso.to_csv(os.path.join(INTER, "M10D_pandisease_k20_loso.csv"), index=False)
    loso.to_csv(os.path.join(TAB, "Table_S87c_M10D_Pandisease_k20_LOSO.csv"), index=False)
    flush_log()
    print("\nDONE")

if __name__ == "__main__":
    main()
