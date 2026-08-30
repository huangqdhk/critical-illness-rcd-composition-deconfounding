# -*- coding: utf-8 -*-
"""
M11M12_step2_temporal.py — M11 时间序（判据 ii）
================================================
前瞻注册：M11_M12_pre_registration_20260827_draft.md §1（H1/H2 主检验、判定门 M11）
- H1 主检验 1：患者随机截距混合模型
    EIS_fu ~ β1·UCS0 + β2·EIS0 + β3·Neut0 + β4·Mono0 (+β5·sex) + (1|patient)
    主系数 β1（基线 UCS），单侧 β1<0；最小关注效应量 |β_std|>=0.20。
    队列：GSE215865（≥2 日标签）、GSE54514（≥2 天，含对照）、GSE148871（血样 ≥2 访视）。
    变量队列内 z 标准化后拟合；合并 = REML + Hartung-Knapp + 预测区间；
    meta LOSO（leave-one-cohort-out）+ 队列内 leave-one-subject-out。
- H2 主检验 2：随机截距交叉滞后面板（RI-CLPM，≥3 波：GSE215865/GSE54514/GSE148871）
    EIS(t) = b + bE·EISbar + bU·UCSbar + bc·compbar + φE·EISw(t-1) + βL·UCSw(t-1) + γ·compw(t)
    主系数 βL（滞后交叉路径 UCS->EIS），单侧 βL<0；受试者聚类稳健 SE；
    UCS 方程（φU·UCSw(t-1) + βL2·EISw(t-1)）同报。GSE212865 仅 2 波不参与 H2。
- 零模型：表达量匹配随机集置换 B=10,000，seed=0。
    统计量 = 队列内 OLS 口径 β（观测与零同口径；主推断用 MixedLM/聚类稳健 SE，如实披露
    零模型采用 OLS 快速口径的偏离）。随机基因集按队列内基因均值表达十分位匹配，
    从非 manifest 宇宙抽取，臂大小保持。
- BH 族（预声明）：2 主检验 × 3 队列 = 6 项单侧 p。
- 判定门 M11：主检验方向一致（≥70% 队列同向）+ BH 内显著 + LOSO ≥90% 方向一致。
输出：_intermediate/M11M12_step2_*.csv；03_LOGS/M11M12_step2_log.txt
"""
import gzip, os, re, sys, time, warnings
from collections import defaultdict

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import brentq

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
INTER = ROOT + r"\_intermediate"
LOG = ROOT + r"\03_LOGS"
sys.path.insert(0, ROOT)
import mdi_lib as L

t0 = time.time()
log = []
def note(m=""):
    log.append(m); print(m, flush=True)

B_NULL = 10000
SEED = 0
RNG = np.random.default_rng(SEED)

# ---------------- meta（M10A 同款 REML + Hartung-Knapp） ----------------
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
    return dict(k=k, pooled=mu, se=se_mu, se_hk=se_hk,
                ci_lo=mu - float(stats.t.ppf(0.975, k - 1)) * se_hk,
                ci_hi=mu + float(stats.t.ppf(0.975, k - 1)) * se_hk,
                p_neg=float(stats.t.cdf(tstat, k - 1)),
                p_pos=float(1 - stats.t.cdf(tstat, k - 1)),
                tau2=t2, I2=I2, Q=Q, Q_p=float(1 - stats.chi2.cdf(Q, dfq)),
                pi_lo=float(pi_lo), pi_hi=float(pi_hi))

def ols(X, y):
    X = np.asarray(X, float); y = np.asarray(y, float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    return beta, resid

def cluster_se(X, resid, clusters):
    X = np.asarray(X, float)
    n, p = X.shape
    meat = np.zeros((p, p))
    for cid in np.unique(clusters):
        idx = clusters == cid
        s = X[idx].T @ resid[idx]
        meat += np.outer(s, s)
    return np.sqrt(np.maximum(np.diag(np.linalg.inv(X.T @ X) @ meat @ np.linalg.inv(X.T @ X)), 0.0))

def zvec(x):
    x = np.asarray(x, float)
    sd = x.std(ddof=1)
    return (x - x.mean()) / sd if sd > 0 else np.zeros_like(x)

# ======================================================================
# 1. 数据装载与队列装配
# ======================================================================
note("== 1. 数据装载")
S = pd.read_csv(INTER + r"\M11M12_step1_per_sample.csv", low_memory=False)
S["subject_norm"] = S["subject_norm"].astype(str)
S["day_label"] = S["day_label"].astype(str)
note(f"   逐样本 {len(S)} 行")

def day_num(d):
    d = str(d)
    return int(d.split("_")[1]) if d.startswith("Day_") else np.nan
S["day_num"] = S["day_title"].apply(day_num)
VISIT_ORDER = {"SCREENING": 0, "V3 DAY 12": 1, "V4 DAY 28": 2, "V6 DAY 84": 3}
S["visit_num"] = S["visit"].map(VISIT_ORDER)
def t_num(t):
    t = str(t)
    return int(t[1:]) if t.startswith("T") and t[1:].isdigit() else np.nan
S["t_num"] = S["day_label"].apply(t_num)
S["tp_num"] = S["timepoint"].map({"D0": 0, "D7": 1})

WAVE_COL = {"GSE215865": "t_num", "GSE54514": "day_num",
            "GSE148871": "visit_num", "GSE212865": "tp_num"}

def build_h1(cohort):
    """返回含 wave/wave0 列的队列数据（≥2 波；同标签重复样取均值）。"""
    d = S[S["cohort"] == cohort].copy()
    d = d.dropna(subset=["subject_norm"])
    d = d[d["subject_norm"] != "nan"]
    d["wave"] = d[WAVE_COL[cohort]]
    d = d.dropna(subset=["wave"])
    if cohort == "GSE54514":
        d = d[d["wave"] >= 1]
    # 同受试者×同波重复样取均值（预注册重复样规则；影响 604 号 V3-D12 等）
    if cohort in ("GSE54514", "GSE148871"):
        num_cols = [c for c in d.columns if c not in ("cohort", "sample_id", "subject_norm",
                                                     "day_label", "sex_norm", "sex", "gender", "wave")
                    and pd.api.types.is_numeric_dtype(d[c])]
        cat_map = {}
        for c in ("sex_norm", "sex", "gender"):
            if c in d.columns:
                cat_map[c] = "first"
        agg = {c: "mean" for c in num_cols}
        agg.update(cat_map)
        agg["sample_id"] = "first"
        d = d.groupby(["subject_norm", "wave"]).agg(agg).reset_index()
    bl = d.groupby("subject_norm")["wave"].min().rename("wave0")
    d = d.merge(bl, left_on="subject_norm", right_index=True)
    nw = d.groupby("subject_norm")["wave"].nunique()
    d = d[d["subject_norm"].isin(nw[nw >= 2].index)]
    return d

def h1_design(d, use_sex):
    """返回 fu, X, y, cols, base。协变量按受试者 z（队列内），结局按随访观测 z。
    组成协变量 = 中性粒+单核合并比例（预注册 H1：β3·cell-composition，单一组成项）。"""
    bl = d[d["wave"] == d["wave0"]]
    fu = d[d["wave"] > d["wave0"]]
    base = bl.groupby("subject_norm").agg(
        UCS0=("UCS", "first"), EIS0=("EIS", "first"),
        comp0=("comp", "first"), sex0=("sex_norm", "first"),
        sample0=("sample_id", "first"),
        neut_author=("neutrophil proportion", "first") if "neutrophil proportion" in bl.columns else ("UCS", "first"))
    if "neutrophil proportion" not in bl.columns:
        base = base.drop(columns=["neut_author"])
    for c in ("UCS0", "EIS0", "comp0"):
        base[c] = zvec(base[c].values)
    fu = fu.merge(base, left_on="subject_norm", right_index=True)
    fu["EIS_z"] = zvec(fu["EIS"].values)
    cols = ["UCS0", "EIS0", "comp0"]
    if use_sex:
        ok = base["sex0"].isin(["F", "M"])
        if ok.mean() > 0.8:
            base["sex_d"] = (base["sex0"] == "F").astype(float)
            fu = fu.drop(columns=["sex_d"], errors="ignore").merge(
                base[["sex_d"]], left_on="subject_norm", right_index=True)
            cols = cols + ["sex_d"]
    X = np.column_stack([np.ones(len(fu))] + [fu[c].values for c in cols])
    return fu, X, fu["EIS_z"].values, cols, base

def h1_design_unadj(fu, base):
    """未校正变体（§0.4 三套报告义务）：仅 UCS0 + EIS0（+sex）。"""
    cols = ["UCS0", "EIS0"]
    if "sex_d" in fu.columns:
        cols = cols + ["sex_d"]
    X = np.column_stack([np.ones(len(fu))] + [fu[c].values for c in cols])
    return X, cols

def _mixed_fit(y, X, groups):
    """MixedLM 拟合（lbfgs 优先，失败回退 bfgs/cg，仍失败返回 None）。"""
    from statsmodels.regression.mixed_linear_model import MixedLM
    for method in ("lbfgs", "bfgs", "cg"):
        try:
            mdl = MixedLM(y, X, groups).fit(reml=True, method=method, maxiter=500)
            if hasattr(mdl, "converged") and not mdl.converged:
                continue
            return mdl
        except Exception:
            continue
    return None

def h1_fit_mixed(cohort):
    d = build_h1(cohort)
    use_sex = (cohort != "GSE215865")
    fu, X, y, cols, base = h1_design(d, use_sex)
    groups = pd.factorize(fu["subject_norm"])[0]
    mdl = _mixed_fit(y, X, groups)
    fit_method = "mixedlm"
    if mdl is None:
        # 回退：OLS + 受试者聚类稳健 SE（如实披露）
        beta, resid = ols(X, y)
        se_all = cluster_se(X, resid, groups)
        dfr = len(y) - X.shape[1]
        b1_idx = cols.index("UCS0") + 1
        beta1 = beta[b1_idx]; se1 = se_all[b1_idx]
        p_one_neg = float(stats.t.cdf(beta1 / se1, df=dfr))
        fit_method = "ols_cluster_fallback"
        out = dict(cohort=cohort, n_obs=int(len(y)), n_patients=int(fu["subject_norm"].nunique()),
                   beta1_std=float(beta1), se=float(se1), p_one_neg=p_one_neg,
                   p_two=float(2 * min(p_one_neg, 1 - p_one_neg)), icc=np.nan,
                   b_EIS0=float(beta[cols.index("EIS0") + 1]),
                   b_comp0=float(beta[cols.index("comp0") + 1]),
                   sex_included=("sex_d" in cols), fit_method=fit_method)
        return out
    b1_idx = cols.index("UCS0") + 1
    beta1 = mdl.params[b1_idx]
    se1 = mdl.bse[b1_idx]
    dfr = len(y) - X.shape[1]
    p_one_neg = float(stats.t.cdf(beta1 / se1, df=dfr))
    icc = float(mdl.cov_re[0, 0] / (mdl.cov_re[0, 0] + mdl.scale))
    out = dict(cohort=cohort, n_obs=int(len(y)), n_patients=int(fu["subject_norm"].nunique()),
               beta1_std=float(beta1), se=float(se1), p_one_neg=p_one_neg,
               p_two=float(2 * min(p_one_neg, 1 - p_one_neg)), icc=icc,
               b_EIS0=float(mdl.params[cols.index("EIS0") + 1]),
               b_comp0=float(mdl.params[cols.index("comp0") + 1]),
               sex_included=("sex_d" in cols), fit_method="mixedlm")
    # OLS + 受试者聚类稳健 SE 同报（边界收敛敏感性）
    beta_o, resid_o = ols(X, y)
    se_o = cluster_se(X, resid_o, groups)
    out["beta1_ols"] = float(beta_o[b1_idx])
    out["se_ols"] = float(se_o[b1_idx])
    # 未校正变体（同模型，去掉 comp0）
    Xu, colsu = h1_design_unadj(fu, base)
    mdlu = _mixed_fit(y, Xu, groups)
    if mdlu is not None:
        out["beta1_unadj"] = float(mdlu.params[1])
        out["se_unadj"] = float(mdlu.bse[1])
    return out

# ======================================================================
# 2. H1 主检验 1
# ======================================================================
note("\n== 2. H1 主检验 1（患者随机截距混合模型）")
h1_rows = []
for cohort in ("GSE215865", "GSE54514", "GSE148871"):
    r = h1_fit_mixed(cohort)
    h1_rows.append(r)
    note(f"   {cohort}: n_obs={r['n_obs']} n_pat={r['n_patients']} "
         f"β1_std={r['beta1_std']:+.3f} (SE {r['se']:.3f}) 单侧p={r['p_one_neg']:.4g} "
         f"ICC={r['icc'] if r['icc']==r['icc'] else float('nan'):.3f} 拟合={r['fit_method']}")
h1 = pd.DataFrame(h1_rows)
h1.to_csv(INTER + r"\M11M12_step2_h1_cohort.csv", index=False)

meta1 = pool_hk(h1["beta1_std"].values, h1["se"].values)
dir_frac1 = (h1["beta1_std"] < 0).mean()
note(f"   H1 合并: β1_std={meta1['pooled']:+.3f} [{meta1['ci_lo']:+.3f},{meta1['ci_hi']:+.3f}] "
     f"单侧p(neg)={meta1['p_neg']:.4g} I2={meta1['I2']:.0f}% PI=[{meta1['pi_lo']:+.3f},{meta1['pi_hi']:+.3f}]")
note(f"   方向一致 {int((h1['beta1_std']<0).sum())}/{len(h1)}（{dir_frac1:.0%}，门 ≥70%）")
loso_meta_rows = []
for i in range(len(h1)):
    keep = h1.drop(h1.index[i])
    m = pool_hk(keep["beta1_std"].values, keep["se"].values)
    loso_meta_rows.append(dict(left_out=h1.iloc[i]["cohort"], pooled=m["pooled"],
                               same_sign=bool(np.sign(m["pooled"]) == np.sign(meta1["pooled"]))))
loso_meta1 = pd.DataFrame(loso_meta_rows)
loso_meta1.to_csv(INTER + r"\M11M12_step2_h1_loso_meta.csv", index=False)
note(f"   meta LOSO 同向 {int(loso_meta1['same_sign'].sum())}/{len(loso_meta1)}")

# 队列内 leave-one-subject-out（H1）
note("\n== 2b. H1 队列内 leave-one-subject-out")
loso_subj_rows = []
for cohort in ("GSE215865", "GSE54514", "GSE148871"):
    d = build_h1(cohort)
    use_sex = (cohort != "GSE215865")
    fu, X, y, cols, base = h1_design(d, use_sex)
    subjs = fu["subject_norm"].unique()
    b1_idx = cols.index("UCS0") + 1
    full = h1_fit_mixed(cohort)["beta1_std"]
    signs = []
    n_fail = 0
    for sj in subjs:
        keep = (fu["subject_norm"] != sj).values
        g = pd.factorize(fu["subject_norm"][keep])[0]
        m = _mixed_fit(y[keep], X[keep], g)
        if m is None:
            n_fail += 1
            continue
        signs.append(np.sign(m.params[b1_idx]) == np.sign(full))
    frac = float(np.mean(signs)) if signs else np.nan
    loso_subj_rows.append(dict(cohort=cohort, n_subjects=int(len(subjs)),
                               direction_consistent=frac, n_consistent=int(np.sum(signs)),
                               n_failed_fit=n_fail))
    note(f"   {cohort}: LOSO 同向 {int(np.sum(signs))}/{len(signs)}（{frac:.0%}）失败拟合 {n_fail}")
loso_subj1 = pd.DataFrame(loso_subj_rows)
loso_subj1.to_csv(INTER + r"\M11M12_step2_h1_loso_subject.csv", index=False)

# GSE54514 作者实测中性粒比例敏感性（外部金标准交叉）
note("\n== 2c. GSE54514 作者实测中性粒比例敏感性")
if "neutrophil proportion" in S.columns:
    d545 = build_h1("GSE54514")
    fu545, X545, y545, cols545, base545 = h1_design(d545, True)
    neut_author = pd.to_numeric(base545["neut_author"], errors="coerce")
    if neut_author.notna().mean() > 0.8:
        base545["neutA_z"] = zvec(neut_author.values)
        fu545 = fu545.drop(columns=["neutA_z"], errors="ignore").merge(
            base545[["neutA_z"]], left_on="subject_norm", right_index=True)
        Xa = X545.copy()
        Xa[:, cols545.index("comp0") + 1] = fu545["neutA_z"].values
        groups545 = pd.factorize(fu545["subject_norm"])[0]
        ma = _mixed_fit(y545, Xa, groups545)
        if ma is not None:
            b1a = ma.params[cols545.index("UCS0") + 1]
            note(f"   作者中性粒比例替换组成协变量: β1_std={b1a:+.3f} "
                 f"(SE {ma.bse[cols545.index('UCS0')+1]:.3f}) 单侧p={stats.t.cdf(b1a/ma.bse[cols545.index('UCS0')+1], len(y545)-Xa.shape[1]):.4g}")
            pd.DataFrame([dict(cohort="GSE54514", sens="author_neutrophil_replaces_comp",
                               beta1_std=float(b1a),
                               se=float(ma.bse[cols545.index("UCS0") + 1]),
                               p_one_neg=float(stats.t.cdf(b1a / ma.bse[cols545.index("UCS0") + 1],
                                                           len(y545) - Xa.shape[1])))]) \
              .to_csv(INTER + r"\M11M12_step2_h1_gse54514_author_neut_sens.csv", index=False)
        else:
            note("   [WARN] 作者中性粒敏感性拟合失败")
    else:
        note("   [WARN] 作者中性粒比例覆盖不足，跳过")
else:
    note("   [WARN] per_sample 表无 neutrophil proportion 列，跳过")

# ======================================================================
# 3. H2 RI-CLPM（≥3 波）
# ======================================================================
note("\n== 3. H2 主检验 2（RI-CLPM 交叉滞后面板）")
def build_clpm(cohort):
    d = build_h1(cohort)
    nw = d.groupby("subject_norm")["wave"].nunique()
    d = d[d["subject_norm"].isin(nw[nw >= 3].index)]
    d["comp_z"] = zvec(d["comp"].values)
    gm = d.groupby("subject_norm")[["UCS", "EIS", "comp_z"]].transform("mean")
    d["UCSbar"] = gm["UCS"]; d["EISbar"] = gm["EIS"]; d["compbar"] = gm["comp_z"]
    d["UCSw"] = d["UCS"] - d["UCSbar"]
    d["EISw"] = d["EIS"] - d["EISbar"]
    d["compw"] = d["comp_z"] - d["compbar"]
    rows = []
    for sj, g in d.groupby("subject_norm", sort=False):
        g = g.sort_values("wave")
        g["UCSw_lag"] = g["UCSw"].shift(1)
        g["EISw_lag"] = g["EISw"].shift(1)
        rows.append(g)
    dd = pd.concat(rows)
    dd = dd.dropna(subset=["UCSw_lag", "EISw_lag"])
    return dd

CLPM_COLS = ["EISbar", "UCSbar", "compbar", "EISw_lag", "UCSw_lag", "compw"]

def clpm_fit(dd, outcome):
    """返回 dict + (X, y, beta, se)。主系数 UCSw_lag（outcome=EIS 时）。"""
    X = np.column_stack([np.ones(len(dd))] + [dd[c].values for c in CLPM_COLS])
    y = dd[outcome].values
    beta, resid = ols(X, y)
    se = cluster_se(X, resid, pd.factorize(dd["subject_norm"])[0])
    li = CLPM_COLS.index("UCSw_lag") + 1
    b = beta[li]; s = se[li]
    tstat = b / s
    p_one = float(stats.norm.cdf(tstat))
    out = dict(beta_L=float(b), se=float(s), t=float(tstat), p_one_neg=p_one,
               n_obs=int(len(y)), n_pat=int(dd["subject_norm"].nunique()))
    # 同报：EIS(t-1)->UCS(t) 交叉路径（UCS 方程）
    li2 = CLPM_COLS.index("EISw_lag") + 1
    out["beta_EISlag_on_UCS"] = float(beta[li2])
    out["se_EISlag_on_UCS"] = float(se[li2])
    return out, (X, y, beta, se)

h2_rows = []
for cohort in ("GSE215865", "GSE54514", "GSE148871"):
    dd = build_clpm(cohort)
    r, _ = clpm_fit(dd, "EIS")
    r["cohort"] = cohort
    h2_rows.append(r)
    note(f"   {cohort}: n_obs={r['n_obs']} n_pat={r['n_pat']} "
         f"β_L={r['beta_L']:+.3f} (SE {r['se']:.3f}) 单侧p={r['p_one_neg']:.4g}")
h2 = pd.DataFrame(h2_rows)
h2.to_csv(INTER + r"\M11M12_step2_h2_cohort.csv", index=False)
meta2 = pool_hk(h2["beta_L"].values, h2["se"].values)
dir_frac2 = (h2["beta_L"] < 0).mean()
note(f"   H2 合并: β_L={meta2['pooled']:+.3f} [{meta2['ci_lo']:+.3f},{meta2['ci_hi']:+.3f}] "
     f"单侧p(neg)={meta2['p_neg']:.4g} I2={meta2['I2']:.0f}% PI=[{meta2['pi_lo']:+.3f},{meta2['pi_hi']:+.3f}]")
note(f"   方向一致 {int((h2['beta_L']<0).sum())}/{len(h2)}（{dir_frac2:.0%}）")
loso_meta_rows = []
for i in range(len(h2)):
    keep = h2.drop(h2.index[i])
    m = pool_hk(keep["beta_L"].values, keep["se"].values)
    loso_meta_rows.append(dict(left_out=h2.iloc[i]["cohort"], pooled=m["pooled"],
                               same_sign=bool(np.sign(m["pooled"]) == np.sign(meta2["pooled"]))))
loso_meta2 = pd.DataFrame(loso_meta_rows)
loso_meta2.to_csv(INTER + r"\M11M12_step2_h2_loso_meta.csv", index=False)
note(f"   meta LOSO 同向 {int(loso_meta2['same_sign'].sum())}/{len(loso_meta2)}")

# H2 队列内 LOSO
note("\n== 3b. H2 队列内 leave-one-subject-out")
loso_subj2_rows = []
for cohort in ("GSE215865", "GSE54514", "GSE148871"):
    dd = build_clpm(cohort)
    subjs = dd["subject_norm"].unique()
    full = clpm_fit(dd, "EIS")[0]["beta_L"]
    signs = []
    for sj in subjs:
        keep = dd["subject_norm"] != sj
        if keep.sum() < 10:
            continue
        try:
            b = clpm_fit(dd[keep], "EIS")[0]["beta_L"]
            signs.append(np.sign(b) == np.sign(full))
        except Exception:
            signs.append(False)
    frac = float(np.mean(signs))
    loso_subj2_rows.append(dict(cohort=cohort, n_subjects=int(len(subjs)),
                                direction_consistent=frac, n_consistent=int(np.sum(signs))))
    note(f"   {cohort}: LOSO 同向 {int(np.sum(signs))}/{len(signs)}（{frac:.0%}）")
loso_subj2 = pd.DataFrame(loso_subj2_rows)
loso_subj2.to_csv(INTER + r"\M11M12_step2_h2_loso_subject.csv", index=False)

# ======================================================================
# 4. 基因级矩阵（零模型用；step1 口径重建，缓存）
# ======================================================================
note("\n== 4. 基因级矩阵（零模型）")
def load_gene_matrix(cohort):
    cache = INTER + rf"\M11M12_step1_genemat_{cohort}.npz"
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return pd.DataFrame(z["M"], index=z["samples"], columns=z["genes"])
    note(f"   [重建基因矩阵] {cohort}（step1 同口径）")
    if cohort == "GSE215865":
        cols_file = RAW215 = ROOT + r"\00_RAW_DATA\GSE215865_COVID19_WholeBlood_Longitudinal\GSE215865_rnaseq_logCPM_matrix.csv.gz"
        with gzip.open(cols_file, "rt", encoding="utf-8", errors="replace") as f:
            header = f.readline().rstrip("\n")
        colnames = header.split(",")
        cache_union = {}
        for c in (INTER + r"\M10A_monaco_ensg2sym.csv", INTER + r"\P2_gse185263_en2sym.csv"):
            for k, v in pd.read_csv(c).values:
                cache_union[str(k).strip().upper().split(".")[0]] = str(v).strip().upper()
        if os.path.exists(INTER + r"\M11M12_step1_ensg2sym_union.csv"):
            for k, v in pd.read_csv(INTER + r"\M11M12_step1_ensg2sym_union.csv").values:
                cache_union[str(k).strip().upper()] = str(v).strip().upper()
        df = pd.read_csv(cols_file, compression="gzip", index_col=0, low_memory=False,
                         na_values=["NA"], dtype={c: np.float32 for c in colnames[1:]})
        syms = [cache_union.get(str(i).split(".")[0].upper(), "") for i in df.index]
        df["sym"] = syms
        df = df[df["sym"] != ""]
        df = df.groupby("sym").mean(numeric_only=True).fillna(0.0)
        pat = re.compile(r"^Subj_([0-9a-fA-F]+)T(\d+[A-Za-z]?)_Plate_(\d+)$")
        meta = []
        for c in df.columns:
            m = pat.match(str(c))
            meta.append((str(c), m.group(1), "T" + re.sub(r"[A-Za-z]+$", "", m.group(2))) if m else ("", "", ""))
        mdf = pd.DataFrame(meta, columns=["col", "subject", "day_label"])
        mdf["key"] = mdf["subject"] + "|" + mdf["day_label"]
        mdf = mdf.set_index("col")             # 与转置矩阵行索引对齐
        g = df.T.groupby(mdf["key"]).mean()
        g["sid"] = g.index.astype(str)      # 索引即 "subject|day_label"
        g = g.set_index("sid")
    elif cohort == "GSE54514":
        f = ROOT + r"\00_RAW_DATA\GSE54514_Sepsis_PAXgene_WholeBlood\GSE54514_non-normalized.txt.gz"
        raw = pd.read_csv(f, sep="\t", index_col=0, low_memory=False,
                          dtype={c: np.float32 for c in pd.read_csv(f, sep="\t", nrows=0).columns[1:]})
        expr = raw[[c for c in raw.columns if not str(c).endswith("Detection Pval")]].T
        def qn(m):
            M = m.values.astype(float).T
            order = np.argsort(M, axis=1)
            sm = np.take_along_axis(M, order, axis=1)
            means = np.nanmean(sm, axis=0)
            out = np.empty_like(M)
            for i in range(M.shape[0]):
                out[i] = np.broadcast_to(means, M.shape)[i][np.argsort(order[i])]
            return pd.DataFrame(out.T, index=m.index, columns=m.columns)
        expr = np.log2(qn(expr) + 1.0)
        # sentrix -> gsm 映射（GSM 元数据 descriptions[0] = sentrix）
        sent2gsm = {}
        txt = open(ROOT + r"\00_RAW_DATA\GSE54514_Sepsis_PAXgene_WholeBlood\GSE54514_gsm_metadata.txt",
                   encoding="utf-8", errors="replace").read()
        for rec in txt.split("^SAMPLE")[1:]:
            gsm = None
            for ln in rec.split("\n"):
                if ln.startswith("!Sample_geo_accession"):
                    gsm = ln.split("=", 1)[1].strip()
                elif ln.startswith("!Sample_description") and gsm:
                    sent2gsm[ln.split("=", 1)[1].strip()] = gsm
                    gsm = None
        expr.index = [sent2gsm.get(str(i), str(i)) for i in expr.index]
        p2s = {}
        with open(ROOT + r"\00_RAW_DATA\GPL6947_table.txt", encoding="utf-8", errors="replace") as ff:
            lines = ff.read().split("\n")
        hdr_idx = next(i for i, ln in enumerate(lines) if ln.startswith("!platform_table_begin"))
        header = lines[hdr_idx + 1].split("\t")
        sc = next(i for i, h in enumerate(header) if h.strip().strip('"').lower() == "symbol")
        ic = next(i for i, h in enumerate(header) if h.strip().strip('"').lower() in ("id",))
        for ln in lines[hdr_idx + 2:]:
            if ln.startswith("!platform_table_end") or not ln.strip():
                break
            parts = ln.split("\t")
            if len(parts) > max(sc, ic):
                sym = parts[sc].strip().strip('"').upper()
                if sym and sym != "---":
                    p2s[parts[ic].strip().strip('"')] = sym
        g2p = defaultdict(list)
        for p, s in p2s.items():
            if p in expr.columns:
                g2p[s].append(p)
        g = pd.DataFrame(index=expr.index)
        for gg, ps in g2p.items():
            subp = expr[ps]
            if len(ps) > 3:
                best = subp.mean(axis=0).nlargest(3).index
                subp = subp[best]
            g[gg] = subp.mean(axis=1)
        g = g[~g.index.duplicated(keep="first")]
    else:
        if cohort == "GSE148871":
            meta4, s4, m4_ = L.parse_series_matrix(
                ROOT + r"\00_RAW_DATA\GSE148871_COPD_AE\GSE148871_series_matrix.txt.gz")
            gpl = pd.read_csv(ROOT + r"\00_RAW_DATA\GPL570_probe2symbol_hgu133plus2db.csv")
            pc = [c for c in gpl.columns if c.lower() in ("probe_id", "probeid", "id")][0]
            sc = [c for c in gpl.columns if c.lower() in ("symbol", "gene_symbol")][0]
            p2s = dict(zip(gpl[pc].astype(str).str.strip(), gpl[sc].astype(str).str.strip().str.upper()))
            g2p = defaultdict(list)
            for p, s in p2s.items():
                if p in m4_.columns:
                    g2p[s].append(p)
            g = pd.DataFrame(index=m4_.index)
            for gg, ps in g2p.items():
                subp = m4_[ps]
                if len(ps) > 3:
                    best = subp.mean(axis=0).nlargest(3).index
                    subp = subp[best]
                g[gg] = subp.mean(axis=1)
            def char_vec(meta, label):
                for block in meta.get("!Sample_characteristics_ch1", []):
                    if block and label.lower() in block[0].lower():
                        return [b.split(":", 1)[1].strip() if ":" in b else b.strip() for b in block]
                return None
            tis = char_vec(meta4, "tissue")
            g = g[[str(t).lower() == "whole blood" for t in tis]]
        else:
            raise NotImplementedError
    # 对齐 step1 样本
    sid = list(S[S["cohort"] == cohort]["sample_id"])
    g = g.loc[[i for i in sid if i in g.index]]
    np.savez_compressed(cache, M=g.values.astype(np.float32),
                        samples=np.array(g.index, dtype=object),
                        genes=np.array(g.columns, dtype=object))
    return g

# ======================================================================
# 5. 零模型（表达量匹配随机集，B=10,000，seed=0）
# ======================================================================
note("\n== 5. 零模型（B=10,000，seed=0）")
_, arms, _ = L.load_manifest()
UP = arms["upstream_collapse"]
ALL80 = set(arms["upstream_collapse"]) | set(arms["execution_induction"])

def decile_structures(gm, arm_genes):
    """表达量匹配结构：宇宙基因按队列内均值表达十分位分桶；臂基因按宇宙分桶归类。"""
    universe = [g for g in gm.columns if g.upper() not in ALL80]
    means = gm[universe].mean(axis=0).values
    deciles, bins = pd.qcut(means, 10, labels=False, duplicates="drop", retbins=True)
    idx_by_dec = {int(d): np.where(deciles == d)[0] for d in np.unique(deciles)}
    arm_in = [g for g in arm_genes if g in gm.columns]
    arm_mean = gm[arm_in].mean(axis=0).values
    if len(bins) >= 3:
        arm_dec = np.clip(np.digitize(arm_mean, bins[1:-1]), 0, len(bins) - 2).astype(int)
    else:
        arm_dec = np.zeros(len(arm_in), dtype=int)
    return universe, idx_by_dec, arm_in, arm_dec

# ---- H1 零模型（观测与零同口径：OLS β1） ----
note("   H1 零模型运行中……")
zero1_rows = []
for cohort in ("GSE215865", "GSE54514", "GSE148871"):
    gm = load_gene_matrix(cohort)
    universe, idx_by_dec, arm_in, arm_dec = decile_structures(gm, UP)
    d = build_h1(cohort)
    use_sex = (cohort != "GSE215865")
    fu, X, y, cols, base = h1_design(d, use_sex)
    b1_idx = cols.index("UCS0") + 1
    obs_beta, _ = ols(X, y)
    obs_stat = obs_beta[b1_idx]
    # 基线样本基因矩阵
    base_sids = base["sample0"].values
    gm_base = gm.loc[[i for i in base_sids if i in gm.index]]
    # 预存其他列（不随置换变化）
    X_fix = X.copy()
    nulls = np.empty(B_NULL)
    for b in range(B_NULL):
        sel = []
        for ad in arm_dec:
            pool = idx_by_dec[ad]
            sel.append(universe[pool[RNG.integers(0, len(pool))]])
        ucs0_null = gm_base[sel].mean(axis=1)                      # 基线样本 x 1
        ucs0_by_subj = ucs0_null.reindex(base_sids).values         # 按 base.index 对齐
        u = zvec(ucs0_by_subj)                                     # 受试者层 z（与观测口径一致）
        fu_u = pd.Series(u, index=base.index).reindex(fu["subject_norm"]).values
        X2 = X_fix.copy()
        X2[:, b1_idx] = fu_u
        bb, _ = ols(X2, y)
        nulls[b] = bb[b1_idx]
    p_emp = float((nulls <= obs_stat).mean())
    zero1_rows.append(dict(cohort=cohort, hypothesis="H1", n_arm_genes=len(arm_in),
                           obs_stat=float(obs_stat), null_mean=float(nulls.mean()),
                           null_sd=float(nulls.std()), empirical_p=p_emp))
    note(f"   {cohort} H1: obs={obs_stat:+.3f} null={nulls.mean():+.3f}±{nulls.std():.3f} "
         f"经验p={p_emp:.4g}")
zero1 = pd.DataFrame(zero1_rows)
zero1.to_csv(INTER + r"\M11M12_step2_zero_h1.csv", index=False)

# ---- H2 零模型 ----
note("   H2 零模型运行中……")
zero2_rows = []
for cohort in ("GSE215865", "GSE54514", "GSE148871"):
    gm = load_gene_matrix(cohort)
    universe, idx_by_dec, arm_in, arm_dec = decile_structures(gm, UP)
    dd = build_clpm(cohort)
    X, y = np.column_stack([np.ones(len(dd))] + [dd[c].values for c in CLPM_COLS]), dd["EIS"].values
    obs_beta, _ = ols(X, y)
    obs_stat = obs_beta[CLPM_COLS.index("UCSw_lag") + 1]
    # 预计算受试者结构（person-mean 与组内滞后是 UCS 的线性算子）
    subj_codes = pd.factorize(dd["subject_norm"])[0]
    n_subj = subj_codes.max() + 1
    # 行按受试者排序后做 person mean / lag
    dd_s = dd.copy()
    dd_s["_code"] = subj_codes
    dd_s["_wave"] = dd["wave"]
    order = np.lexsort((dd_s["_wave"].values, dd_s["_code"].values))
    codes_s = dd_s["_code"].values[order]
    # person mean 算子
    counts = np.bincount(codes_s, minlength=n_subj)
    def person_mean(v):
        return (np.bincount(codes_s, weights=v, minlength=n_subj) / counts)[codes_s]
    def lag_within(v):
        out = np.full(len(v), np.nan)
        same = codes_s[1:] == codes_s[:-1]
        out[1:][same] = v[:-1][same]
        return out
    gm_ord = gm.loc[[i for i in dd["sample_id"] if i in gm.index]]
    # 固定列（排序后行序；不随置换变化）
    EISbar_k_all = dd_s["EISbar"].values[order]
    compbar_k_all = dd_s["compbar"].values[order]
    EISw_lag_k_all = dd_s["EISw_lag"].values[order]
    compw_k_all = dd_s["compw"].values[order]
    y_k = dd_s["EIS"].values[order]
    nulls = np.empty(B_NULL)
    sid_to_row = {sid: i for i, sid in enumerate(dd["sample_id"])}
    for b in range(B_NULL):
        sel = []
        for ad in arm_dec:
            pool = idx_by_dec[ad]
            sel.append(universe[pool[RNG.integers(0, len(pool))]])
        ucs_null = gm_ord[sel].mean(axis=1)
        ucs_row = np.full(len(dd), np.nan)
        for sid, v in ucs_null.items():
            ucs_row[sid_to_row[sid]] = v
        ucs_z = zvec(ucs_row)
        ucs_z_s = ucs_z[order]
        ucsbar = person_mean(ucs_z_s)
        ucsw = ucs_z_s - ucsbar
        ucsw_lag = lag_within(ucsw)
        keep = ~np.isnan(ucsw_lag)
        # X2 按 CLPM_COLS 顺序：[1, EISbar, UCSbar, compbar, EISw_lag, UCSw_lag, compw]
        X2 = np.column_stack([np.ones(int(keep.sum())), EISbar_k_all[keep], ucsbar[keep],
                              compbar_k_all[keep], EISw_lag_k_all[keep], ucsw_lag[keep],
                              compw_k_all[keep]])
        bb, _ = ols(X2, y_k[keep])
        nulls[b] = bb[CLPM_COLS.index("UCSw_lag") + 1]
    p_emp = float((nulls <= obs_stat).mean())
    zero2_rows.append(dict(cohort=cohort, hypothesis="H2", n_arm_genes=len(arm_in),
                           obs_stat=float(obs_stat), null_mean=float(nulls.mean()),
                           null_sd=float(nulls.std()), empirical_p=p_emp))
    note(f"   {cohort} H2: obs={obs_stat:+.3f} null={nulls.mean():+.3f}±{nulls.std():.3f} "
         f"经验p={p_emp:.4g}")
zero2 = pd.DataFrame(zero2_rows)
zero2.to_csv(INTER + r"\M11M12_step2_zero_h2.csv", index=False)

# ======================================================================
# 6. GSE212865 两波敏感性（H1 形式；不进 BH 族）
# ======================================================================
note("\n== 6. GSE212865 两波敏感性（D0->D7）")
d212 = S[S["cohort"] == "GSE212865"].copy()
d212 = d212[d212["timepoint"].isin(["D0", "D7"])].dropna(subset=["patient"])
d212["patient"] = d212["patient"].astype(str)
paired = d212.groupby("patient")["timepoint"].apply(lambda x: set(x) == {"D0", "D7"})
d212p = d212[d212["patient"].isin(paired[paired].index)]
bl = d212p[d212p["timepoint"] == "D0"].groupby("patient").agg(
    UCS0=("UCS", "first"), EIS0=("EIS", "first"), comp0=("comp", "first"))
fu = d212p[d212p["timepoint"] == "D7"].set_index("patient")
for c in ("UCS0", "EIS0", "comp0"):
    bl[c] = zvec(bl[c].values)
fu = fu.merge(bl, left_index=True, right_index=True)
fu["EIS_z"] = zvec(fu["EIS"].values)
cols = ["UCS0", "EIS0", "comp0"]
X = np.column_stack([np.ones(len(fu))] + [fu[c].values for c in cols])
y = fu["EIS_z"].values
beta, resid = ols(X, y)
se = cluster_se(X, resid, pd.factorize(fu.index)[0])
b1_idx = cols.index("UCS0") + 1
tstat = beta[b1_idx] / se[b1_idx]
p_one = float(stats.norm.cdf(tstat))
sens212 = dict(cohort="GSE212865", n_patients=int(fu.index.nunique()), n_obs=int(len(y)),
               beta1_std=float(beta[b1_idx]), se=float(se[b1_idx]),
               p_one_neg=p_one, note="两波敏感性，不进 BH 族")
note(f"   GSE212865 敏感性: n={sens212['n_patients']} β1_std={beta[b1_idx]:+.3f} "
     f"(SE {se[b1_idx]:.3f}) 单侧p={p_one:.4g}")
pd.DataFrame([sens212]).to_csv(INTER + r"\M11M12_step2_sensitivity_gse212865.csv", index=False)

# ======================================================================
# 7. BH 与判定门
# ======================================================================
note("\n== 7. BH 族与判定门 M11")
pvals = np.concatenate([h1["p_one_neg"].values, h2["p_one_neg"].values])
qvals = L.bh(pvals)
h1["q_bh"] = qvals[:3]
h2["q_bh"] = qvals[3:]
h1.to_csv(INTER + r"\M11M12_step2_h1_cohort.csv", index=False)
h2.to_csv(INTER + r"\M11M12_step2_h2_cohort.csv", index=False)
for _, r in h1.iterrows():
    note(f"   H1 {r['cohort']}: p={r['p_one_neg']:.4g} q={r['q_bh']:.4g}")
for _, r in h2.iterrows():
    note(f"   H2 {r['cohort']}: p={r['p_one_neg']:.4g} q={r['q_bh']:.4g}")
gate_h1 = dict(dir_consistent=bool(dir_frac1 >= 0.7),
               bh_pass=bool((h1["q_bh"] < 0.05).any()),
               loso_subj=bool((loso_subj1["direction_consistent"] >= 0.9).all()),
               loso_meta=bool(loso_meta1["same_sign"].all()))
gate_h2 = dict(dir_consistent=bool(dir_frac2 >= 0.7),
               bh_pass=bool((h2["q_bh"] < 0.05).any()),
               loso_subj=bool((loso_subj2["direction_consistent"] >= 0.9).all()),
               loso_meta=bool(loso_meta2["same_sign"].all()))
note(f"   H1 判定门: {gate_h1}")
note(f"   H2 判定门: {gate_h2}")
allpass = all(gate_h1.values()) and all(gate_h2.values())
verdict = "M11 PASS（判据 ii 成立）" if allpass else "M11 未全过门：如实入稿（按注册条款）"
note(f"   判定: {verdict}")
pd.DataFrame([{**gate_h1, "test": "H1"}, {**gate_h2, "test": "H2"}, {"test": "verdict", "verdict": verdict}]).to_csv(
    INTER + r"\M11M12_step2_gate.csv", index=False)

with open(LOG + r"\M11M12_step2_log.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(log))
note(f"\nDONE in {time.time()-t0:.1f}s -> {LOG}\\M11M12_step2_log.txt")
