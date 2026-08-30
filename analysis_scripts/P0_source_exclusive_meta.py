# -*- coding: utf-8 -*-
"""
P0_source_exclusive_meta.py — Phase 0 任务4：来源互斥 meta 重算（H1–H3）
=========================================================================
依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md（OSF https://osf.io/C7RYD/，2026-08-20 冻结）
     §2 终点 / §3 排除与纳入 / §4 统计方法；判定门槛 §5。
输入：_intermediate/M2_per_sample_scores.csv（冻结 score_version=mdi_v1.0）
     04_AUDIT_GOVERNANCE/SAMPLE_MANIFEST_v1.0.csv（GSE212865 患者号与时点映射）
规则：
     1. GSE310929 按 Dataset 列拆回原始来源；GSE185263/GSE32707 来源行弃用
        （同一样本采用其直接处理版本，每生物学样本只计一次）；
     2. 内部来源 case={Sepsis, Sepsis-Shock, Sepsis-ARDs} vs ctrl=Control；
     3. 直并行用各自预注册主对比（GSE185263: Sepsis_COVID vs Control；
        GSE32707: ARDS_d0 vs Control；GSE212865: Covid19_SDRA vs Control，基线时点）；
     4. GSE188309（disease-only）单样本 t 检验，不入合并；GSE148871 无 case/ctrl，不入合并；
     5. 合并：REML + Hartung-Knapp，单侧 p，预测区间，LOSO 稳定性。
输出：_intermediate/P0_dataset_effects.csv / P0_meta_pooled.csv / P0_meta_loso.csv
     04_AUDIT_GOVERNANCE/P0_Source_Exclusive_Meta_Report.md
     03_LOGS/P0_source_exclusive_meta_log.txt
"""
import os, sys
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import brentq

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
M2CSV = ROOT + r"\_intermediate\M2_per_sample_scores.csv"
MAN = ROOT + r"\04_AUDIT_GOVERNANCE\SAMPLE_MANIFEST_v1.0.csv"
OLDMETA = ROOT + r"\_intermediate\M2_meta.csv"
INTER = ROOT + r"\_intermediate"
REPORT = ROOT + r"\04_AUDIT_GOVERNANCE\P0_Source_Exclusive_Meta_Report.md"
LOGF = ROOT + r"\03_LOGS\P0_source_exclusive_meta_log.txt"

sys.path.insert(0, ROOT)
import mdi_lib as L

log = []
def note(m=""):
    log.append(m)
    try:
        print(m, flush=True)
    except UnicodeEncodeError:
        print(m.encode("ascii", "replace").decode(), flush=True)

# ---------------- meta 数学：REML tau2 + Hartung-Knapp ----------------
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
    """REML 随机效应 + Hartung-Knapp 区间/单侧 p + 预测区间。g>0 = case>ctrl。"""
    g = np.asarray(g, float); se = np.asarray(se, float); k = len(g)
    t2 = reml_tau2(g, se)
    w = 1.0 / (se ** 2 + t2)
    mu = float(np.sum(w * g) / np.sum(w))
    se_mu = float(np.sqrt(1.0 / np.sum(w)))
    se_hk = float(np.sqrt(np.sum(w * (g - mu) ** 2) / ((k - 1) * np.sum(w))))
    tcrit = float(stats.t.ppf(0.975, k - 1))
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
                ci_lo=mu - tcrit * se_hk, ci_hi=mu + tcrit * se_hk,
                p_pos=float(1 - stats.t.cdf(tstat, k - 1)),
                p_neg=float(stats.t.cdf(tstat, k - 1)),
                tau2=t2, I2=I2, Q=Q, Q_p=float(1 - stats.chi2.cdf(Q, dfq)),
                pi_lo=float(pi_lo), pi_hi=float(pi_hi))

# ---------------- 数据与层（stratum）装配 ----------------
df = pd.read_csv(M2CSV, low_memory=False)
man = pd.read_csv(MAN, low_memory=False)
df["MDI_ssgsea"] = df["EIS_ssgsea"] - df["UCS_ssgsea"]

CASES_310 = ["Sepsis", "Sepsis - Shock", "Sepsis - ARDs"]
MIN_N = 3
strata = []          # (stratum_id, case_df, ctrl_df)
dropped = []

note("== 1. GSE310929 拆源（弃用其中 GSE185263/GSE32707 来源行——采用直并行版本）")
g3 = df[df["cohort"] == "GSE310929"]
for ds, s in g3.groupby("Dataset"):
    if ds in ("GSE185263", "GSE32707"):
        dropped.append((ds, "duplicate-of-direct", len(s)))
        continue
    case = s[s["group"].isin(CASES_310)]
    ctrl = s[s["group"] == "Control"]
    if len(case) >= MIN_N and len(ctrl) >= MIN_N:
        strata.append((f"310:{ds}", case, ctrl))
    else:
        dropped.append((f"310:{ds}", f"insufficient (case={len(case)}, ctrl={len(ctrl)})", len(s)))
note(f"   可合并来源层 {sum(1 for x in strata if x[0].startswith('310:'))} 个；"
     f"弃用 {len(dropped)} 个来源（case-only / control-only / 双计）")

note("\n== 2. 直并行（各自预注册主对比）")
for cohort, cgrp, label in [("GSE185263", ["Sepsis_COVID"], "SepsisCOVID_vs_Control"),
                            ("GSE32707", ["ARDS_d0"], "ARDSd0_vs_Control")]:
    s = df[df["cohort"] == cohort]
    strata.append((cohort, s[s["group"].isin(cgrp)], s[s["group"] == "Control"]))
    note(f"   {cohort}: {label}")

note("\n== 3. GSE212865 基线时点规则（患者号=subcohort，保 D0；健康对照无时点全保留）")
s = df[df["cohort"] == "GSE212865"].merge(
    man[man["dataset_id"] == "GSE212865"][["gsm", "timepoint", "subcohort"]],
    left_on="sample", right_on="gsm", how="left")
n0 = len(s)
s = s[s["timepoint"].isna() | (s["timepoint"].astype(str) == "D0")]
dup = s.duplicated(subset=["subcohort", "group"], keep=False) & s["subcohort"].notna()
if dup.any():
    note(f"   基线内仍有 {int(dup.sum())} 行同患者同分组重复，保留首行")
    s = s[~dup | ~s.duplicated(subset=["subcohort", "group"], keep="first")]
note(f"   {n0} → {len(s)}（基线化后）；SDRA vs Control")
strata.append(("GSE212865_base", s[s["group"] == "Covid19_SDRA"], s[s["group"] == "Control"]))

note(f"\n== 4. 不入合并的队列（按冻结计划）")
d188 = df[df["cohort"] == "GSE188309"]
note(f"   GSE188309 disease-only n={len(d188)}：单样本 t 检验，另报")
note(f"   GSE148871 无 case/ctrl 分组（血/痰×治疗访问），不入合并（记录）")

# ---------------- 数据集级效应 ----------------
FAMS_CONF = ["UCS", "EIS", "MDI"]
FAMS_SENS = ["MDI_nomt", "UCS_ssgsea", "EIS_ssgsea", "MDI_ssgsea"]
rows = []
for sid, case, ctrl in strata:
    for fam in FAMS_CONF + FAMS_SENS:
        g, gse = L.hedges_g(ctrl[fam].values, case[fam].values)
        rows.append(dict(stratum=sid, family=fam, n_case=len(case), n_ctrl=len(ctrl),
                         case_mean=case[fam].mean(), ctrl_mean=ctrl[fam].mean(),
                         hedges_g=g, hedges_g_se=gse))
eff = pd.DataFrame(rows)
eff["gene_set_version"] = L.GENE_SET_VERSION
eff["score_version"] = L.SCORE_VERSION
eff.to_csv(INTER + r"\P0_dataset_effects.csv", index=False)
note(f"\n== 5. 数据集级效应已写入 P0_dataset_effects.csv（{len(eff)} 行）")
for fam in FAMS_CONF:
    sub = eff[eff["family"] == fam].sort_values("hedges_g")
    signs = (sub["hedges_g"] > 0).sum()
    note(f"   {fam}: {len(sub)} 层，正方向 {signs} 层；"
         f"范围 [{sub['hedges_g'].min():+.2f}, {sub['hedges_g'].max():+.2f}]")

# ---------------- 合并 + LOSO ----------------
pooled_rows, loso_rows = [], []
for fam in FAMS_CONF + FAMS_SENS:
    sub = eff[eff["family"] == fam].dropna()
    res = pool_hk(sub["hedges_g"].values, sub["hedges_g_se"].values)
    res.update(family=fam, role=("confirmatory" if fam in FAMS_CONF else "sensitivity"),
               cohorts=";".join(sub["stratum"]))
    pooled_rows.append(res)
    if fam in FAMS_CONF + ["MDI_nomt"]:
        full_sign = np.sign(res["pooled_g"])
        stab = []
        for i in range(len(sub)):
            keep = sub.drop(sub.index[i])
            r2 = pool_hk(keep["hedges_g"].values, keep["hedges_g_se"].values)
            stab.append(np.sign(r2["pooled_g"]) == full_sign)
            loso_rows.append(dict(family=fam, left_out=sub.iloc[i]["stratum"],
                                  pooled_g=r2["pooled_g"], same_sign=bool(np.sign(r2["pooled_g"]) == full_sign)))
        res["loso_stable"] = all(stab)
        note(f"   POOLED {fam}: g={res['pooled_g']:+.3f} [{res['ci_lo']:+.3f},{res['ci_hi']:+.3f}] "
             f"HK 单侧p(pos)={res['p_pos']:.4g} 单侧p(neg)={res['p_neg']:.4g} "
             f"I2={res['I2']:.0f}% PI=[{res['pi_lo']:+.2f},{res['pi_hi']:+.2f}] "
             f"LOSO {'STABLE' if res['loso_stable'] else 'UNSTABLE'}")

pool_df = pd.DataFrame(pooled_rows)
loso_df = pd.DataFrame(loso_rows)
pool_df.to_csv(INTER + r"\P0_meta_pooled.csv", index=False)
loso_df.to_csv(INTER + r"\P0_meta_loso.csv", index=False)

# ---------------- GSE188309 单样本 + 锚点 ----------------
note("\n== 6. GSE188309 单样本口径（z 空间 vs 0）")
one188 = {}
for fam in FAMS_CONF + ["MDI_nomt"]:
    t = stats.ttest_1samp(d188[fam].dropna(), 0)
    one188[fam] = (float(d188[fam].mean()), float(t.statistic), float(t.pvalue))
    note(f"   {fam}: mean={d188[fam].mean():+.3f} t={t.statistic:+.2f} p={t.pvalue:.3g}")

note("\n== 7. 自洽性锚点（GSE185263 直并行，方向应与 Table S2 一致：UCS<0, EIS>0）")
a = eff[(eff["stratum"] == "GSE185263") & (eff["family"].isin(FAMS_CONF))]
anchor_ok = bool((a[a.family == "UCS"]["hedges_g"] < 0).all() and (a[a.family == "EIS"]["hedges_g"] > 0).all())
note(f"   锚点检查: {'PASS' if anchor_ok else 'FAIL'}")

# ---------------- 旧结果对照（原 M2_meta，已判定无效，仅并列展示） ----------------
def df_md(d, cols):
    if d is None or len(d) == 0:
        return "(empty)"
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, rr in d.iterrows():
        vals = []
        for c in cols:
            v = rr[c]
            vals.append(f"{v:.4g}" if isinstance(v, (int, float, np.floating)) else str(v))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)

old_txt = ""
if os.path.exists(OLDMETA):
    om = pd.read_csv(OLDMETA)
    old_txt = df_md(om, ["score", "n_cohorts", "pooled_g", "meta_p", "I2"])

# ---------------- 报告 ----------------
def fam_verdict(fam, hyp):
    r = pool_df[(pool_df.family == fam) & (pool_df.role == "confirmatory")].iloc[0]
    p = r["p_neg"] if hyp == "neg" else r["p_pos"]
    ok = (p < 0.05) and bool(r.get("loso_stable", False))
    return ("支持" if ok else "不支持"), r, p

v1, r1, p1 = fam_verdict("UCS", "neg")     # H1: g<0
v2, r2, p2 = fam_verdict("EIS", "pos")     # H2: g>0
v3, r3, p3 = fam_verdict("MDI", "pos")     # H3: g>0

os.makedirs(os.path.dirname(REPORT), exist_ok=True)
with open(REPORT, "w", encoding="utf-8") as f:
    f.write("# P0-4 来源互斥 Meta 重算报告（H1–H3）\n\n")
    f.write("- 日期：2026-08-20；冻结依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md（OSF https://osf.io/C7RYD/）\n")
    f.write(f"- 脚本：P0_source_exclusive_meta.py；score_version={L.SCORE_VERSION}；gene_set_version={L.GENE_SET_VERSION}\n")
    f.write(f"- 合并层数 k={int(r1['k'])}（GSE310929 拆源 {sum(1 for x in strata if x[0].startswith('310:'))} + 直并行 3）；"
            f"弃用来源 {len(dropped)} 个\n\n")
    f.write("## 判定（冻结计划 §5 口径）\n\n")
    f.write(f"| 假设 | 判定 | pooled g | 95%CI(HK) | 单侧p | LOSO |\n|---|---|---|---|---|---|\n")
    for h, v, r, p in [("H1 UCS g<0", v1, r1, p1), ("H2 EIS g>0", v2, r2, p2), ("H3 MDI g>0", v3, r3, p3)]:
        f.write(f"| {h} | **{v}** | {r['pooled_g']:+.3f} | [{r['ci_lo']:+.3f}, {r['ci_hi']:+.3f}] | {p:.4g} | "
                f"{'稳定' if r.get('loso_stable') else '不稳定'} |\n")
    f.write("\nPhase 0 总判定需结合 P0-5（供者级 pseudobulk）后给出；本表仅 H1–H3 层面。\n")
    f.write("\n## 数据集级效应（P0_dataset_effects.csv 摘要）\n\n")
    for fam in FAMS_CONF:
        sub = eff[eff["family"] == fam].sort_values("hedges_g")
        f.write(f"\n### {fam}\n\n| 层 | n_case | n_ctrl | g |\n|---|---|---|---|\n")
        for _, rr in sub.iterrows():
            f.write(f"| {rr['stratum']} | {rr['n_case']} | {rr['n_ctrl']} | {rr['hedges_g']:+.3f} |\n")
    f.write("\n## 敏感性族（MDI_nomt / ssGSEA）合并结果\n\n")
    f.write(df_md(pool_df[pool_df.role == "sensitivity"], ["family", "k", "pooled_g", "ci_lo", "ci_hi", "p_pos", "I2"]))
    f.write("\n\n## GSE188309 单样本检验（不入合并）\n\n| 族 | mean | t | p |\n|---|---|---|---|\n")
    for fam, (m, t, p) in one188.items():
        f.write(f"| {fam} | {m:+.3f} | {t:+.2f} | {p:.3g} |\n")
    f.write("\n注：GSE188309 的臂评分在其自身队列内整体 z 标准化，全队列均值恒为 0，"
            "故『单样本 vs 0』检验因评分构造而无信息量（均值/统计量全为 0 属必然，不构成『无差异』的证据）。"
            "该队列按冻结计划不入合并；此结构性限制在此如实披露。\n")
    f.write("\n## 与旧（无效）M2 meta 的并列对照\n\n")
    f.write(old_txt if isinstance(old_txt, str) and old_txt else "（旧 M2_meta.csv 不存在）")
    f.write("\n\n旧 k=3 meta 已被 P0-3 审计判定无效（GSE310929 全量包含另两队列），上表仅作对照展示，不具证据效力。\n")
    f.write("\n## 弃用来源清单\n\n| 来源 | 原因 | n |\n|---|---|---|\n")
    for ds, why, n in dropped:
        f.write(f"| {ds} | {why} | {n} |\n")
    f.write("\n## 注记\n\n- GSE310929 拆源层的臂评分继承自 GSE310929 整体运行的 z 标准化（冻结 score_version），\n"
            "  层内病例-对照效应（Hedges' g）对整体仿射变换不变，此近似在报告中如实披露。\n"
            "- GSE212865 已按患者号基线化（保 D0；健康对照无时点）。\n")

with open(LOGF, "w", encoding="utf-8") as f:
    f.write("\n".join(log))
print("\nDONE P0-4 ->", REPORT)
