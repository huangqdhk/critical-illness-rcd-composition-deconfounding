# -*- coding: utf-8 -*-
"""
M2_step3_outcomes_endotypes.py — M2 第 3 步：临床结局关联 + 增量价值 + 四象限内型
================================================================================
预注册：M2_pre_registration_20260817.md
- GSE310929（脓毒症图谱）：28 天死亡 Cox（未调整/年龄性别调整）+ Harrell C + log-rank 四分位
  + logistic OR + 校准（斜率/截距/Brier）+ DCA + 增量价值（基线=年龄+性别+MolecularSubtype；
  ΔC-index / IDI / 分类 NRI）+ 敏感性（剔除 GSE185263/GSE32707 来源样本）
- GSE188309（CAP）：logistic 住院死亡 + ROC-AUC + Hosmer-Lemeshow + DCA + 三分位趋势（数据就位后）
- 四象限内型：队列内 UCS/EIS 中位数切分（命名预注册固定）→ 结局梯度 + IGP（k=10, 500 置换）
输出：_intermediate/M2_outcomes.csv, M2_endotypes.csv
"""
import os, io, sys
import numpy as np
import pandas as pd
from scipy import stats as st
sys.path.insert(0, r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS")
import mdi_lib as L
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index as cindex
from sklearn.metrics import roc_auc_score

ROOT = L.ROOT
INTER = ROOT + r"\_intermediate"
df = pd.read_csv(INTER + r"\M2_per_sample_scores.csv")
log = []
def note(m):
    log.append(m); print(m, flush=True)

outcome_rows = []

def cox_model(d, covariates, dur="tte", ev="event"):
    """Cox PH；返回 (summary_df, cindex, n, n_events)。"""
    dd = d[[dur, ev] + list(covariates)].dropna()
    cph = CoxPHFitter()
    cph.fit(dd, duration_col=dur, event_col=ev)
    ci = cindex(dd[dur], -cph.predict_partial_hazard(dd), dd[ev])
    return cph, ci, len(dd), int(dd[ev].sum())

def logistic_summary(X, y, names):
    """sklearn LogisticRegression；返回每变量 OR/CI/p。"""
    from sklearn.linear_model import LogisticRegression
    from scipy.special import expit
    mod = LogisticRegression(penalty=None, max_iter=1000).fit(X, y)
    phat = mod.predict_proba(X)[:, 1]
    n = len(y)
    # 标准误（信息矩阵）
    Xd = np.column_stack([np.ones(n), X])
    W = phat * (1 - phat)
    H = (Xd * W[:, None]).T @ Xd
    cov = np.linalg.inv(H + 1e-10 * np.eye(H.shape[0]))
    se = np.sqrt(np.diag(cov))
    coefs = np.concatenate([[mod.intercept_[0]], mod.coef_[0]])
    out = []
    for nm, b, s in zip(["intercept"] + list(names), coefs, se):
        z = b / s
        p = 2 * (1 - st.norm.cdf(abs(z)))
        out.append(dict(term=nm, coef=b, se=s, OR=np.exp(b),
                        CI_lo=np.exp(b - 1.96 * s), CI_hi=np.exp(b + 1.96 * s), p=p))
    return out, phat

# ================================================================ GSE310929
note("== GSE310929 28 天死亡结局")
sep_all = df[df["cohort"] == "GSE310929"].copy()
# 结局构造（预注册：DaySurvivalEdited + TimetoEventEdited）
sep_all["event"] = (sep_all["DaySurvivalEdited"] == "Died").astype(float)
sep_all["tte"] = pd.to_numeric(sep_all["TimetoEventEdited"], errors="coerce")
sep_log = sep_all[sep_all["DaySurvivalEdited"].notna()].copy()          # 主分析：有 28d 状态即可
sep_out = sep_all[sep_all["DaySurvivalEdited"].notna() & sep_all["tte"].notna()].copy()  # Cox 亚组
sep_out.loc[sep_out["tte"] <= 0, "tte"] = 0.5
sep_out["tte"] = sep_out["tte"].clip(upper=28.0)
for d_ in (sep_log, sep_out):
    d_["MDI_sd"] = d_["MDI"] / d_["MDI"].std(ddof=1)
    d_["age"] = pd.to_numeric(d_["Age"], errors="coerce")
    d_["sex_male"] = (d_["Gender"].astype(str).str.lower() == "male").astype(float)
note(f"  主分析（logistic）n={len(sep_log)}（死亡 {int(sep_log['event'].sum())}）；"
     f"Cox 时间亚组 n={len(sep_out)}（事件 {int(sep_out['event'].sum())}）")

# --- 主分析：logistic 28d（未调整）
lr0_rows, phat0 = logistic_summary(sep_log[["MDI_sd"]].values, sep_log["event"].values, ["MDI_per_SD"])
m0 = [r for r in lr0_rows if r["term"] == "MDI_per_SD"][0]
auc0 = roc_auc_score(sep_log["event"], phat0)
note(f"  Logistic 未调整（n={len(sep_log)}）: OR={m0['OR']:.3f} (95%CI {m0['CI_lo']:.3f}-{m0['CI_hi']:.3f}), p={m0['p']:.3g}, AUC={auc0:.3f}")
outcome_rows.append(dict(cohort="GSE310929", model="Logistic_unadj_28d", term="MDI_per_SD",
                         HR_or_OR=m0["OR"], CI_lo=m0["CI_lo"], CI_hi=m0["CI_hi"],
                         p=m0["p"], cindex=auc0, n=len(sep_log), n_events=int(sep_log["event"].sum())))

# --- logistic 年龄性别调整
sep_adj = sep_log[sep_log["age"].notna()].copy()
lr1_rows, phat1 = logistic_summary(sep_adj[["MDI_sd", "age", "sex_male"]].values,
                                   sep_adj["event"].values, ["MDI_per_SD", "age", "sex_male"])
m1 = [r for r in lr1_rows if r["term"] == "MDI_per_SD"][0]
auc1 = roc_auc_score(sep_adj["event"], phat1)
note(f"  Logistic 调整（n={len(sep_adj)}）: OR={m1['OR']:.3f} (95%CI {m1['CI_lo']:.3f}-{m1['CI_hi']:.3f}), p={m1['p']:.3g}, AUC={auc1:.3f}")
outcome_rows.append(dict(cohort="GSE310929", model="Logistic_age_sex", term="MDI_per_SD",
                         HR_or_OR=m1["OR"], CI_lo=m1["CI_lo"], CI_hi=m1["CI_hi"],
                         p=m1["p"], cindex=auc1, n=len(sep_adj), n_events=int(sep_adj["event"].sum())))

# --- Cox（593 时间亚组，未调整 + 调整）
cph0, ci0, n0, ev0 = cox_model(sep_out, ["MDI_sd"])
hr0 = cph0.summary.loc["MDI_sd"]
note(f"  Cox 未调整（n={n0}，事件 {ev0}）: HR={np.exp(hr0['coef']):.3f} (95%CI {np.exp(hr0['coef']-1.96*hr0['se(coef)']):.3f}-{np.exp(hr0['coef']+1.96*hr0['se(coef)']):.3f}), p={hr0['p']:.3g}, C={ci0:.3f}")
outcome_rows.append(dict(cohort="GSE310929", model="Cox_unadj", term="MDI_per_SD",
                         HR_or_OR=np.exp(hr0["coef"]), CI_lo=np.exp(hr0["coef"]-1.96*hr0["se(coef)"]),
                         CI_hi=np.exp(hr0["coef"]+1.96*hr0["se(coef)"]), p=hr0["p"],
                         cindex=ci0, n=n0, n_events=ev0))
sub_adj2 = sep_out[sep_out["age"].notna()]
cph1, ci1, n1, ev1 = cox_model(sub_adj2, ["MDI_sd", "age", "sex_male"])
hr1 = cph1.summary.loc["MDI_sd"]
note(f"  Cox 年龄性别调整（n={n1}）: HR={np.exp(hr1['coef']):.3f}, p={hr1['p']:.3g}, C={ci1:.3f}")
outcome_rows.append(dict(cohort="GSE310929", model="Cox_age_sex", term="MDI_per_SD",
                         HR_or_OR=np.exp(hr1["coef"]), CI_lo=np.exp(hr1["coef"]-1.96*hr1["se(coef)"]),
                         CI_hi=np.exp(hr1["coef"]+1.96*hr1["se(coef)"]), p=hr1["p"],
                         cindex=ci1, n=n1, n_events=ev1))

# --- 四分位分层（logistic 主口径 + Cox 时间亚组；Q1 为参照）
sep_log["MDI_q"] = pd.qcut(sep_log["MDI"], 4, labels=["Q1_low", "Q2", "Q3", "Q4_high"])
note("  四分位 28d 死亡率（logistic 口径 n=%d）:" % len(sep_log))
quart_rows = []
for q in ["Q1_low", "Q2", "Q3", "Q4_high"]:
    qd = sep_log[sep_log["MDI_q"] == q]
    if q == "Q1_low":
        note(f"  {q}: n={len(qd)}, 死亡率={qd['event'].mean():.3f}, OR vs Q1=1.000 (参照)")
        quart_rows.append(dict(cohort="GSE310929", quartile=q, n=len(qd), mortality=qd['event'].mean(),
                               OR_vs_Q1=1.0, CI_lo=np.nan, CI_hi=np.nan, p=np.nan))
        continue
    # logistic vs Q1
    sub = sep_log[sep_log["MDI_q"].isin(["Q1_low", q])]
    dummy = (sub["MDI_q"] == q).astype(float).values
    lrq, _ = logistic_summary(dummy[:, None], sub["event"].values, [q])
    rq = [r for r in lrq if r["term"] == q][0]
    note(f"  {q}: n={len(qd)}, 死亡率={qd['event'].mean():.3f}, OR vs Q1={rq['OR']:.3f} (95%CI {rq['CI_lo']:.3f}-{rq['CI_hi']:.3f}), p={rq['p']:.3g}")
    quart_rows.append(dict(cohort="GSE310929", quartile=q, n=len(qd), mortality=qd['event'].mean(),
                           OR_vs_Q1=rq["OR"], CI_lo=rq["CI_lo"], CI_hi=rq["CI_hi"], p=rq["p"]))
# 趋势（有序四分位 logistic）
lr_tr, _ = logistic_summary(sep_log["MDI_q"].astype("category").cat.codes.values[:, None],
                            sep_log["event"].values, ["quartile_ord"])
rt = [r for r in lr_tr if r["term"] == "quartile_ord"][0]
note(f"  四分位有序 logistic 趋势: OR/级={rt['OR']:.3f}, p={rt['p']:.3g}")
quart_rows.append(dict(cohort="GSE310929", quartile="ordinal_trend", n=len(sep_log),
                       mortality=np.nan, OR_vs_Q1=rt["OR"], CI_lo=rt["CI_lo"], CI_hi=rt["CI_hi"], p=rt["p"]))
pd.DataFrame(quart_rows).to_csv(INTER + r"\M2_quartiles.csv", index=False)

# --- 校准 + DCA（logistic 主模型，n=2436）
intc, slope = L.calibration_curve_stats(sep_log["event"].values, phat0)
br = L.brier_score(sep_log["event"].values, phat0)
dca = L.dca_net_benefit(sep_log["event"].values, phat0)
note(f"  校准: 截距={intc:.3f}, 斜率={slope:.3f}, Brier={br:.4f}")
dca_rows = []
for t, (nb_m, nb_a, _) in dca.items():
    note(f"  DCA t={t}: NB_model={nb_m:+.4f} NB_all={nb_a:+.4f}")
    dca_rows.append(dict(cohort="GSE310929", threshold=t, NB_model=nb_m, NB_treat_all=nb_a, NB_treat_none=0.0))
cal_rows = [dict(cohort="GSE310929", calib_intercept=intc, calib_slope=slope, brier=br)]
pd.DataFrame(dca_rows).to_csv(INTER + r"\M2_dca.csv", index=False)
pd.DataFrame(cal_rows).to_csv(INTER + r"\M2_calibration.csv", index=False)

# --- 增量价值（基线=年龄+性别+MolecularSubtype）
sub_inc = sep_log[sep_log["age"].notna() & sep_log["MolecularSubtype"].notna()].copy()
note(f"\n== 增量价值（基线=age+sex+分子亚型 C1-C4；n={len(sub_inc)}）")
baseX = pd.get_dummies(sub_inc[["age", "sex_male", "MolecularSubtype"]],
                       columns=["MolecularSubtype"], drop_first=False).astype(float)
plusX = baseX.copy()
plusX["MDI_sd"] = sub_inc["MDI_sd"].values
from sklearn.linear_model import LogisticRegression
mb = LogisticRegression(penalty=None, max_iter=2000).fit(baseX.values, sub_inc["event"].values)
mp = LogisticRegression(penalty=None, max_iter=2000).fit(plusX.values, sub_inc["event"].values)
pb = mb.predict_proba(baseX.values)[:, 1]
pp = mp.predict_proba(plusX.values)[:, 1]
from sklearn.metrics import roc_auc_score
auc_b = roc_auc_score(sub_inc["event"], pb)
auc_p = roc_auc_score(sub_inc["event"], pp)
dAUC = auc_p - auc_b
ev = sub_inc["event"].values == 1
nev = ~ev
idi = (pp[ev].mean() - pb[ev].mean()) - (pp[nev].mean() - pb[nev].mean())
# 分类 NRI（风险三分位）
terc = np.quantile(pb, [1/3, 2/3])
def risk_cat(p):
    return np.where(p < terc[0], 0, np.where(p < terc[1], 1, 2))
cb_, cp_ = risk_cat(pb), risk_cat(pp)
nri_ev = (np.mean(cp_[ev] > cb_[ev]) - np.mean(cp_[ev] < cb_[ev]))
nri_nev = (np.mean(cp_[nev] < cb_[nev]) - np.mean(cp_[nev] > cb_[nev]))
nri = nri_ev + nri_nev
note(f"  AUC 基线={auc_b:.4f} → +MDI={auc_p:.4f}；ΔAUC={dAUC:+.4f}；IDI={idi:+.4f}；分类NRI={nri:+.4f}")
outcome_rows.append(dict(cohort="GSE310929", model="Increment", term="MDI_vs_base(age+sex+subtype)",
                         HR_or_OR=auc_p, CI_lo=auc_b, CI_hi=dAUC, p=np.nan, cindex=np.nan,
                         idi=idi, nri=nri,
                         n=len(sub_inc), n_events=int(sub_inc["event"].sum())))

# --- 敏感性：剔除 GSE185263/GSE32707 来源样本（logistic 主口径）
sens = sep_log[~sep_log["Dataset"].isin(["GSE185263", "GSE32707"])]
sens["MDI_sd"] = sens["MDI"] / sens["MDI"].std(ddof=1)
lr_s, phat_s = logistic_summary(sens[["MDI_sd"]].values, sens["event"].values, ["MDI_per_SD"])
ms = [r for r in lr_s if r["term"] == "MDI_per_SD"][0]
auc_s = roc_auc_score(sens["event"], phat_s)
note(f"\n== 敏感性（剔除 GSE185263/GSE32707 来源）: n={len(sens)}, 事件 {int(sens['event'].sum())}, OR={ms['OR']:.3f}, p={ms['p']:.3g}, AUC={auc_s:.3f}")
outcome_rows.append(dict(cohort="GSE310929_excl_internal", model="Logistic_unadj_28d", term="MDI_per_SD",
                         HR_or_OR=ms["OR"], CI_lo=ms["CI_lo"], CI_hi=ms["CI_hi"], p=ms["p"],
                         cindex=auc_s, n=len(sens), n_events=int(sens["event"].sum())))

# --- MDI_nomt 敏感性（logistic 主口径）
note("\n== MDI_nomt 敏感性（GSE310929）")
if "MDI_nomt" in sep_log.columns:
    sep_log["MDI_nomt_sd"] = sep_log["MDI_nomt"] / sep_log["MDI_nomt"].std(ddof=1)
    lr_n, phat_n = logistic_summary(sep_log[["MDI_nomt_sd"]].values, sep_log["event"].values, ["MDI_nomt_per_SD"])
    mn = [r for r in lr_n if r["term"] == "MDI_nomt_per_SD"][0]
    auc_n = roc_auc_score(sep_log["event"], phat_n)
    note(f"  MDI_nomt: OR={mn['OR']:.3f} (95%CI {mn['CI_lo']:.3f}-{mn['CI_hi']:.3f}), p={mn['p']:.3g}, AUC={auc_n:.3f}")
    outcome_rows.append(dict(cohort="GSE310929", model="Logistic_unadj_28d_nomt", term="MDI_nomt_per_SD",
                             HR_or_OR=mn["OR"], CI_lo=mn["CI_lo"], CI_hi=mn["CI_hi"], p=mn["p"],
                             cindex=auc_n, n=len(sep_log), n_events=int(sep_log["event"].sum())))

# ================================================================ GSE188309（若已加载）
cap = df[df["cohort"] == "GSE188309"]
if len(cap) and cap["mortality"].notna().any():
    note("\n== GSE188309 住院死亡（logistic）")
    cap = cap[cap["mortality"].notna()].copy()
    cap["MDI_sd"] = cap["MDI"] / cap["MDI"].std(ddof=1)
    lr_cap, phat_cap = logistic_summary(cap[["MDI_sd"]].values, cap["mortality"].values, ["MDI_per_SD"])
    mc = [r for r in lr_cap if r["term"] == "MDI_per_SD"][0]
    auc_cap = roc_auc_score(cap["mortality"], phat_cap)
    hl_c, hl_p = L.hosmer_lemeshow(cap["mortality"].values, phat_cap)
    note(f"  CAP: OR={mc['OR']:.3f} (95%CI {mc['CI_lo']:.3f}-{mc['CI_hi']:.3f}), p={mc['p']:.3g}, AUC={auc_cap:.3f}, HL p={hl_p:.3f}")
    dca_cap = L.dca_net_benefit(cap["mortality"].values, phat_cap)
    for t, (nb_m, nb_a, _) in dca_cap.items():
        note(f"  DCA t={t}: NB_model={nb_m:+.4f} NB_all={nb_a:+.4f}")
    # 年龄性别调整
    cap_a = cap[cap["age"].notna()].copy()
    cap_a["Sex"] = pd.to_numeric(cap_a["Sex"], errors="coerce")
    cap_a = cap_a[cap_a["Sex"].notna()]
    if len(cap_a) > 40:
        lr_cap2, _ = logistic_summary(cap_a[["MDI_sd", "age", "Sex"]].values, cap_a["mortality"].values,
                                      ["MDI_per_SD", "age", "Sex"])
        mc2 = [r for r in lr_cap2 if r["term"] == "MDI_per_SD"][0]
        note(f"  CAP 调整: OR={mc2['OR']:.3f}, p={mc2['p']:.3g}")
    # 三分位趋势
    cap["MDI_t"] = pd.qcut(cap["MDI"], 3, labels=["T1", "T2", "T3"])
    tab = pd.crosstab(cap["MDI_t"], cap["mortality"]).values.T
    zca, pca = L.cochran_armitage(tab)
    note(f"  三分位 Cochran-Armitage: z={zca:+.2f}, p={pca:.3g}；死亡率: {tab[1]/tab.sum(axis=0)}")
    outcome_rows.append(dict(cohort="GSE188309", model="Logistic_unadj", term="MDI_per_SD",
                             HR_or_OR=mc["OR"], CI_lo=mc["CI_lo"], CI_hi=mc["CI_hi"],
                             p=mc["p"], cindex=auc_cap, n=len(cap), n_events=int(cap["mortality"].sum())))
else:
    note("\n== GSE188309 未加载（注释待装），结局分析跳过")

# ================================================================ 四象限内型
note("\n== 四象限内型（队列内 UCS/EIS 中位数切分，命名预注册固定）")
endo_rows = []
QUAD_DEF = {  # (EIS高, UCS低)
    (True,  False): "Q1_execution_dominant",
    (True,  True):  "Q2_dual_high",
    (False, False): "Q3_dual_low",
    (False, True):  "Q4_collapse_dominant",
}
for cohort, sub0 in df.groupby("cohort"):
    sub = sub0.copy()
    eis_hi = sub["EIS"] >= sub["EIS"].median()
    ucs_lo = sub["UCS"] <= sub["UCS"].median()
    sub["quadrant"] = [QUAD_DEF[(bool(a), bool(b))] for a, b in zip(eis_hi, ucs_lo)]
    dist = sub["quadrant"].value_counts(normalize=True).to_dict()
    # IGP
    X = sub[["UCS", "EIS"]].values
    igp, igp_p = L.igp_permutation_pvalue(X, sub["quadrant"].values, k=10, n_perm=500)
    note(f"  {cohort}: 象限分布={ {k: round(v,3) for k,v in dist.items()} }")
    note(f"  {cohort}: IGP={igp:.3f}, 置换 p={igp_p:.3f}")
    # 结局关联（有结局字段的队列）
    if cohort == "GSE310929" and "DaySurvivalEdited" in sub.columns:
        q_out = []
        for q in sorted(sub["quadrant"].unique()):
            qd = sub[(sub["quadrant"] == q) & sub["DaySurvivalEdited"].notna()]
            mort = (qd["DaySurvivalEdited"] == "Died").mean()
            q_out.append(dict(cohort=cohort, quadrant=q, n=len(qd), mortality_rate=mort))
        ref = [q for q in q_out if q["quadrant"] == "Q1_execution_dominant"][0]
        for q in q_out:
            q["RR_vs_Q1"] = q["mortality_rate"] / ref["mortality_rate"]
            note(f"  {cohort} {q['quadrant']}: n={q['n']}, 28d死亡率={q['mortality_rate']:.3f}, RR vs Q1={q['RR_vs_Q1']:.2f}")
        endo_rows += q_out
    if cohort == "GSE188309" and "mortality" in sub.columns:
        q_out = []
        for q in sorted(sub["quadrant"].unique()):
            qd = sub[(sub["quadrant"] == q) & sub["mortality"].notna()]
            mort = qd["mortality"].mean()
            q_out.append(dict(cohort=cohort, quadrant=q, n=len(qd), mortality_rate=mort))
        for q in q_out:
            note(f"  {cohort} {q['quadrant']}: n={q['n']}, 住院死亡率={q['mortality_rate']:.3f}")
        endo_rows += q_out

od = pd.DataFrame(outcome_rows)
od["gene_set_version"] = L.GENE_SET_VERSION
od["score_version"] = L.SCORE_VERSION
od.to_csv(INTER + r"\M2_outcomes.csv", index=False)
ed = pd.DataFrame(endo_rows)
if len(ed):
    ed["gene_set_version"] = L.GENE_SET_VERSION
    ed["score_version"] = L.SCORE_VERSION
ed.to_csv(INTER + r"\M2_endotypes.csv", index=False)
with open(INTER + r"\M2_step3_log.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(log))
print("\nDONE step3")
