# -*- coding: utf-8 -*-
"""
M10B_step3_cross_cohort_meta.py — M10 检验 B 跨队列合并判定
===========================================================
前瞻注册：M10_M13_M14_pre_registration_20260824.md §1.3
- 队列：GSE158055（PBMC）、GSE145926（BALF，把握度不足作描述）、GSE216009（全血）、
  GSE180578（PBMC，RAW 重建）
- 主检验：髓系（单核/巨噬）细胞类型内部 UCS 的跨队列方向一致性与随机效应合并
  （REML + Hartung-Knapp，P0 同款）
- 判定门：>=3 个独立 scRNA 队列中髓系内部 UCS 一致下降（>=2/3 队列方向一致且合并显著）
  -> 支持细胞内在成分（H_int）
输出：_intermediate/M10B_myeloid_UCS_meta.csv / M10B_myeloid_UCS_layers.csv
"""
import os, sys, time, warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
INTER = ROOT + r"\_intermediate"
GOV = ROOT + r"\04_AUDIT_GOVERNANCE"

sys.path.insert(0, ROOT)
import mdi_lib as L
from scipy import stats
from scipy.optimize import brentq

t0 = time.time()
log = []
def note(m=""):
    log.append(m); print(m, flush=True)

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
    return dict(k=k, pooled_g=mu, se_hk=se_hk,
                ci_lo=mu - float(stats.t.ppf(0.975, k - 1)) * se_hk,
                ci_hi=mu + float(stats.t.ppf(0.975, k - 1)) * se_hk,
                p_neg=float(stats.t.cdf(tstat, k - 1)),
                p_pos=float(1 - stats.t.cdf(tstat, k - 1)),
                I2=I2, Q=Q, Q_p=float(1 - stats.chi2.cdf(Q, dfq)),
                pi_lo=float(pi_lo), pi_hi=float(pi_hi))

# ---------------- 各队列髓系 UCS 效应 ----------------
note("== 1. 髓系 UCS 层效应")
MYELOID_CT = {
    "GSE158055": ["Mono_c14"],
    "GSE145926": ["Mono_c14", "Macrophage"],
    "GSE216009": ["Classical_monocytes", "Non-classical_monocytes",
                  "Mature_neutrophils", "S100A8-9_hi_neutrophils"],
    "GSE180578": ["Classical monocytes", "Intermediate monocytes", "Non-classical monocytes"],
}

def load_scores(path, cohort, cond_map=None):
    d = pd.read_csv(path)
    d["cohort"] = cohort
    if "condition" in d.columns and cond_map:
        d["condition"] = d["condition"].map(cond_map)
    return d

p0 = pd.read_csv(INTER + r"\P0_pseudobulk_arm_scores.csv")
p0["cohort"] = p0["dataset"]
g216 = pd.read_csv(INTER + r"\M10B_gse216009_arm_scores.csv")
g180 = pd.read_csv(INTER + r"\M10B_gse180578_arm_scores.csv")
g180["cohort"] = "GSE180578"
# 供者独立性：covid 仅保留 timepoint==1；healthy 全保留
g180["timepoint"] = g180["timepoint"].astype(str)
tp_num = pd.to_numeric(g180["timepoint"], errors="coerce")
g180 = g180[(g180["condition"] == "healthy") | (tp_num == 1.0)]
# GSE180578 表只存 raw——按 P0 口径在细胞型内重建 z（UCS/EIS/MDI）
def zc180(s):
    s = np.asarray(s, float); sd = s.std(ddof=1)
    return (s - s.mean()) / sd if sd > 0 else s * 0
g180["UCS"] = np.nan; g180["EIS"] = np.nan; g180["MDI"] = np.nan
for ct in g180["cell_type"].unique():
    m = g180["cell_type"] == ct
    g180.loc[m, "UCS"] = zc180(g180.loc[m, "UCS_raw"])
    g180.loc[m, "EIS"] = zc180(g180.loc[m, "EIS_raw"])
    g180.loc[m, "MDI"] = g180.loc[m, "EIS"] - g180.loc[m, "UCS"]
allsc = pd.concat([p0[["cohort", "sampleID", "cell_type", "condition", "UCS", "EIS", "MDI"]],
                   g216[["cohort", "sampleID", "cell_type", "condition", "UCS", "EIS", "MDI"]],
                   g180[["cohort", "sampleID", "cell_type", "condition", "UCS", "EIS", "MDI"]]],
                  ignore_index=True)
allsc["sampleID"] = allsc["sampleID"].astype(str)
# 条件标签统一
cond_map = {"COVID_severe": "case", "Healthy": "ctrl", "Sepsis": "case", "Control": "ctrl",
            "covid": "case", "healthy": "ctrl"}
allsc["grp"] = allsc["condition"].map(lambda x: cond_map.get(str(x), str(x)))
allsc = allsc[allsc["grp"].isin(["case", "ctrl"])]

layers = []
for cohort, cts in MYELOID_CT.items():
    sub = allsc[allsc["cohort"] == cohort]
    for ct in cts:
        s = sub[sub["cell_type"].astype(str) == ct]
        a = s[s["grp"] == "case"]["UCS"].dropna().values
        b = s[s["grp"] == "ctrl"]["UCS"].dropna().values
        if len(a) < 3 or len(b) < 3:
            note(f"   {cohort} {ct}: 不足（{len(a)}v{len(b)}）——跳过")
            continue
        g, gse = L.hedges_g(b, a)   # g>0 = case 更高；UCS 下降 = g<0
        layers.append(dict(cohort=cohort, cell_type=ct, n_case=len(a), n_ctrl=len(b),
                           hedges_g=g, se=gse))
        note(f"   {cohort} {ct:<28} n={len(a)}v{len(b)} UCS g={g:+.3f} se={gse:.3f}")
lay = pd.DataFrame(layers)
lay.to_csv(INTER + r"\M10B_myeloid_UCS_layers.csv", index=False)

# ---------------- 合并（主口径：单核/巨噬，每队列一个代表层） ----------------
note("\n== 2. 合并（每队列单核/巨噬 UCS）")
# 代表层：GSE158055 Mono_c14；GSE145926 Mono_c14（若有）；GSE216009 Classical_monocytes；
# GSE180578 Classical monocytes
rep_sel = {("GSE158055", "Mono_c14"), ("GSE145926", "Mono_c14"),
           ("GSE216009", "Classical_monocytes"), ("GSE180578", "Classical monocytes")}
rep = lay[lay.apply(lambda r: (r["cohort"], r["cell_type"]) in rep_sel, axis=1)]
note(f"   代表层 {len(rep)} 个：")
for _, r in rep.iterrows():
    note(f"     {r['cohort']} {r['cell_type']}: g={r['hedges_g']:+.3f}")
m = pool_hk(rep["hedges_g"].values, rep["se"].values)
dir_frac = (rep["hedges_g"] < 0).mean()
note(f"   合并：g={m['pooled_g']:+.3f} 95%CI(HK)=[{m['ci_lo']:+.3f},{m['ci_hi']:+.3f}] "
     f"单侧p(neg)={m['p_neg']:.4g} I2={m['I2']:.0f}% PI=[{m['pi_lo']:+.3f},{m['pi_hi']:+.3f}]")
note(f"   方向一致（g<0）: {int((rep['hedges_g']<0).sum())}/{len(rep)}（{dir_frac:.0%}）")
# LOSO
loso = []
full_sign = np.sign(m["pooled_g"])
for i in range(len(rep)):
    keep = rep.drop(rep.index[i])
    r2 = pool_hk(keep["hedges_g"].values, keep["se"].values)
    loso.append(dict(left_out=rep.iloc[i]["cohort"] + "|" + rep.iloc[i]["cell_type"],
                     pooled_g=r2["pooled_g"], same_sign=bool(np.sign(r2["pooled_g"]) == full_sign)))
loso_df = pd.DataFrame(loso)
loso_df.to_csv(INTER + r"\M10B_myeloid_UCS_loso.csv", index=False)
note(f"   LOSO 同向 {int(loso_df['same_sign'].sum())}/{len(loso_df)}")

meta_row = dict(family="UCS_myeloid", k=len(rep), pooled_g=m["pooled_g"],
                ci_lo=m["ci_lo"], ci_hi=m["ci_hi"], p_neg=m["p_neg"], p_pos=m["p_pos"],
                I2=m["I2"], pi_lo=m["pi_lo"], pi_hi=m["pi_hi"],
                direction_consistent_frac=float(dir_frac), loso_stable=bool(loso_df["same_sign"].all()))
pd.DataFrame([meta_row]).to_csv(INTER + r"\M10B_myeloid_UCS_meta.csv", index=False)

# 敏感性：加中性粒细胞层（GSE216009 两个中性粒层 + GSE180578 无中性粒）
note("\n== 3. 敏感性：单核+中性粒全髓系层")
m2l = lay.copy()
m3 = pool_hk(m2l["hedges_g"].values, m2l["se"].values)
note(f"   全髓系层（k={len(m2l)}）：g={m3['pooled_g']:+.3f} [{m3['ci_lo']:+.3f},{m3['ci_hi']:+.3f}] "
     f"单侧p(neg)={m3['p_neg']:.4g} 方向一致 {(m2l['hedges_g']<0).mean():.0%}")

with open(ROOT + r"\03_LOGS\M10B_cross_cohort_log.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(log))
note(f"\nDONE in {time.time()-t0:.1f}s")
