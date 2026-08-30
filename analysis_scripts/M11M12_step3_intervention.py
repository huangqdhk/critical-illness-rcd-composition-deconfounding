# -*- coding: utf-8 -*-
"""
M11M12_step3_intervention.py — M12 干预响应（判据 iii）
=======================================================
前瞻注册：M11_M12_pre_registration_20260827_draft.md §2（H3/H4 主检验、判定门 M12 方案 A）
- H3（确认性，GSE106878，注册后首次触碰）：
    每患者 d_res = MDI_resid(post) − MDI_resid(pre)。
    主检验 = 臂间单侧 Mann-Whitney（氢化可的松 vs 安慰剂；H1: HC d_res < PL d_res）；
    组内配对单侧 Wilcoxon（HC 臂，H1: median(d_res) < 0）同报（同一 BH 族）。
    强制组成同报：同臂同一时点的细胞组成变化（中性粒/单核比例）同款检验。
    MESI：组内 P(d_res < 0) ≥ 0.65。
- H4（仅确认性下降检验显著时判定）：
    ΔMDI_obs ~ ΔMDI_comp 回归；组成拟合变化与观测下降同向且 ΔMDI_comp/ΔMDI_obs ≥ 0.5
    → "组成回退"；比值 < 0.5（含 ≤0）→ "状态可逆"。
- 复现层（披露的注册前观察，不进 BH 族）：GSE148871 NEMI 臂 SCREENING→V4-D28
    同款组内 Wilcoxon（d_res 与 d_obs 同报；d_obs 应复现 ΔMDI≈−0.33, p=0.016 锚点）。
- 零模型：表达量匹配随机集置换 B=10,000，seed=0（双臂同替换；统计量 = 臂间 MWU 单侧 z）。
- LOSO：leave-one-patient-out 方向一致性（≥90%）。
- 判定门 M12：GSE106878 组间 MWU 过 BH + GSE148871 同向复现 → 判据 iii 成立；
  GSE106878 阴性 → 如实阴性入稿。
输出：_intermediate/M11M12_step3_*.csv；03_LOGS/M11M12_step3_log.txt
"""
import gzip, os, re, sys, time, warnings
from collections import defaultdict

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import pandas as pd
from scipy import stats

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

def zvec(x):
    x = np.asarray(x, float)
    sd = x.std(ddof=1)
    return (x - x.mean()) / sd if sd > 0 else np.zeros_like(x)

_, arms, _ = L.load_manifest()
UP, EX = arms["upstream_collapse"], arms["execution_induction"]
ALL80 = set(UP) | set(EX)

S = pd.read_csv(INTER + r"\M11M12_step1_per_sample.csv", low_memory=False)
S["subject_norm"] = S["subject_norm"].astype(str)

# ======================================================================
# 1. GSE106878 患者级 d_res
# ======================================================================
note("== 1. GSE106878 患者级变化分")
d = S[S["cohort"] == "GSE106878"].copy()
d["patient"] = d["individual"].astype(str)
d["tp"] = d["timepoint"].astype(str)
note(f"   样本 {len(d)}；臂 {dict(d['treatment'].value_counts())}")
piv = d.pivot_table(index=["patient", "treatment"],
                    columns="tp",
                    values=["MDI", "MDI_resid", "MDI_comp", "UCS", "EIS", "Neut", "Mono", "comp"],
                    aggfunc="first").reset_index()
piv.columns = ["_".join(c).strip("_") for c in piv.columns]
pairs = piv.dropna(subset=["MDI_resid_Post(24h)", "MDI_resid_Pre"]).copy()
note(f"   完全配对患者 {len(pairs)}（HC {int((pairs['treatment']=='hydrocortisone').sum())} / "
     f"PL {int((pairs['treatment']=='placebo').sum())}）")
for c in ("MDI", "MDI_resid", "MDI_comp", "UCS", "EIS", "Neut", "Mono", "comp"):
    pairs["d_" + c] = pairs[f"{c}_Post(24h)"] - pairs[f"{c}_Pre"]
hc = pairs[pairs["treatment"] == "hydrocortisone"]
pl = pairs[pairs["treatment"] == "placebo"]

# 主检验：臂间单侧 MWU（HC < PL）
u_stat, p_mwu = stats.mannwhitneyu(hc["d_MDI_resid"], pl["d_MDI_resid"], alternative="less")
n1, n2 = len(hc), len(pl)
z_mwu = (u_stat - n1 * n2 / 2) / np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
cliff = stats.mannwhitneyu(hc["d_MDI_resid"], pl["d_MDI_resid"],
                           alternative="two-sided") if False else None
# Cliff's delta：P(HC < PL) - P(HC > PL)
dd_ = 0.0
for a in hc["d_MDI_resid"]:
    dd_ += np.sum(pl["d_MDI_resid"] > a) - np.sum(pl["d_MDI_resid"] < a)
cliff_delta = dd_ / (n1 * n2)
# 组内 HC 单侧 Wilcoxon
w_stat, p_wil = stats.wilcoxon(hc["d_MDI_resid"], alternative="less")
p_decline = float((hc["d_MDI_resid"] < 0).mean())
note(f"   臂间 MWU（单侧 less）: U={u_stat:.1f} z={z_mwu:+.3f} p={p_mwu:.4g} "
     f"Cliff's δ={cliff_delta:+.3f}")
note(f"   组内 HC Wilcoxon（单侧）: p={p_wil:.4g} | P(d_res<0)={p_decline:.0%} "
     f"（MESI 门 ≥65%）| 中位 d_res={hc['d_MDI_resid'].median():+.3f}")
note(f"   PL 组内（描述，双侧）: 中位 d_res={pl['d_MDI_resid'].median():+.3f} "
     f"p={stats.wilcoxon(pl['d_MDI_resid']).pvalue:.4g}")

# BH 族（2 检验）
q = L.bh([p_mwu, p_wil])
note(f"   BH 族: 组间 q={q[0]:.4g}，组内 q={q[1]:.4g}")

# 组成同报（双侧）
comp_rows = []
for arm_name, arm_df in (("hydrocortisone", hc), ("placebo", pl)):
    for cvar in ("Neut", "Mono", "comp", "MDI_comp"):
        c = arm_df["d_" + cvar]
        w2, p2 = stats.wilcoxon(c)
        comp_rows.append(dict(arm=arm_name, variable="d_" + cvar, n=len(c),
                              median=float(c.median()), wilcoxon_p_two=float(p2)))
comp = pd.DataFrame(comp_rows)
note(f"   组成同报：HC d_Neut 中位 {comp[(comp['arm']=='hydrocortisone')&(comp['variable']=='d_Neut')]['median'].iloc[0]:+.3f}"
     f"（p={comp[(comp['arm']=='hydrocortisone')&(comp['variable']=='d_Neut')]['wilcoxon_p_two'].iloc[0]:.4g}）；"
     f"PL d_Neut 中位 {comp[(comp['arm']=='placebo')&(comp['variable']=='d_Neut')]['median'].iloc[0]:+.3f}"
     f"（p={comp[(comp['arm']=='placebo')&(comp['variable']=='d_Neut')]['wilcoxon_p_two'].iloc[0]:.4g}）")

# LOSO（leave-one-patient-out；方向 = 检验统计量符号与全数据一致）
mwu_z_loso = []
wil_dir_loso = []
w_full, _ = stats.wilcoxon(hc["d_MDI_resid"], alternative="less")
dir_full_wil = w_full - len(hc) * (len(hc) + 1) / 4
for i in range(len(pairs)):
    keep = pairs.drop(pairs.index[i])
    hc_k = keep[keep["treatment"] == "hydrocortisone"]
    pl_k = keep[keep["treatment"] == "placebo"]
    if len(hc_k) < 3 or len(pl_k) < 3:
        continue
    u_k, _ = stats.mannwhitneyu(hc_k["d_MDI_resid"], pl_k["d_MDI_resid"], alternative="less")
    z_k = (u_k - len(hc_k) * len(pl_k) / 2) / np.sqrt(len(hc_k) * len(pl_k) * (len(hc_k) + len(pl_k) + 1) / 12)
    mwu_z_loso.append(np.sign(z_k) == np.sign(z_mwu))
    w_k, _ = stats.wilcoxon(hc_k["d_MDI_resid"], alternative="less")
    dir_k = w_k - len(hc_k) * (len(hc_k) + 1) / 4
    wil_dir_loso.append(np.sign(dir_k) == np.sign(dir_full_wil))
loso_frac_mwu = float(np.mean(mwu_z_loso))
loso_frac_wil = float(np.mean(wil_dir_loso))
note(f"   LOSO 方向一致: 臂间 MWU {np.mean(mwu_z_loso):.0%}、组内 Wilcoxon {np.mean(wil_dir_loso):.0%}")

# ======================================================================
# 2. H4 归因（仅当确认性下降检验显著）
# ======================================================================
note("\n== 2. H4 归因")
confirmatory_pass = (q[0] < 0.05)
h4_rows = []
if confirmatory_pass:
    # 回归：ΔMDI_obs ~ ΔMDI_comp（全患者 + 分臂）
    for name, sub in (("all", pairs), ("hydrocortisone", hc), ("placebo", pl)):
        X = np.column_stack([np.ones(len(sub)), sub["d_MDI_comp"].values])
        y = sub["d_MDI"].values
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = np.corrcoef(sub["d_MDI_comp"], sub["d_MDI"])[0, 1]
        h4_rows.append(dict(group=name, n=len(sub),
                            slope_obs_on_comp=float(beta[1]), intercept=float(beta[0]), r=float(r),
                            mean_dMDI_comp=float(sub["d_MDI_comp"].mean()),
                            mean_dMDI_obs=float(sub["d_MDI"].mean())))
    share_hc = h4_rows[1]["mean_dMDI_comp"] / h4_rows[1]["mean_dMDI_obs"] if h4_rows[1]["mean_dMDI_obs"] != 0 else np.nan
    same_dir = (h4_rows[1]["mean_dMDI_comp"] < 0) and (h4_rows[1]["mean_dMDI_obs"] < 0)
    adjudication = ("组成回退" if (same_dir and share_hc >= 0.5) else "状态可逆")
    note(f"   HC 臂: mean ΔMDI_obs={h4_rows[1]['mean_dMDI_obs']:+.3f}, "
         f"mean ΔMDI_comp={h4_rows[1]['mean_dMDI_comp']:+.3f}, "
         f"同向份额={share_hc:+.3f} → 判定: {adjudication}")
    for h in h4_rows:
        h["adjudication"] = adjudication if h["group"] == "hydrocortisone" else ""
else:
    note("   确认性检验未显著 → H4 仅描述报告（不执行判定）")
    adjudication = "未触发（确认性检验未显著）"
    for name, sub in (("all", pairs), ("hydrocortisone", hc), ("placebo", pl)):
        X = np.column_stack([np.ones(len(sub)), sub["d_MDI_comp"].values])
        y = sub["d_MDI"].values
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = np.corrcoef(sub["d_MDI_comp"], sub["d_MDI"])[0, 1]
        h4_rows.append(dict(group=name, n=len(sub),
                            slope_obs_on_comp=float(beta[1]), intercept=float(beta[0]), r=float(r),
                            mean_dMDI_comp=float(sub["d_MDI_comp"].mean()),
                            mean_dMDI_obs=float(sub["d_MDI"].mean()), adjudication=adjudication))
h4 = pd.DataFrame(h4_rows)

# ======================================================================
# 3. GSE148871 复现层（披露的注册前观察，不进 BH 族）
# ======================================================================
note("\n== 3. GSE148871 复现层（NEMI 臂 SCREENING->V4-D28）")
g = S[S["cohort"] == "GSE148871"].copy()
g = g[g["treatment"] == "All NEMI"]
g = g[g["visit"].isin(["SCREENING", "V4 DAY 28"])]
g["patient"] = g["subject_norm"]
gpiv = g.pivot_table(index="patient", columns="visit",
                     values=["MDI", "MDI_resid", "MDI_comp", "Neut", "Mono", "comp"],
                     aggfunc="first").reset_index()
gpiv.columns = ["_".join(c).strip("_") for c in gpiv.columns]
gpr = gpiv.dropna(subset=["MDI_resid_SCREENING", "MDI_resid_V4 DAY 28"]).copy()
for c in ("MDI", "MDI_resid", "MDI_comp", "Neut", "Mono", "comp"):
    gpr["d_" + c] = gpr[f"{c}_V4 DAY 28"] - gpr[f"{c}_SCREENING"]
note(f"   NEMI 配对 {len(gpr)} 例")
w_res, p_res = stats.wilcoxon(gpr["d_MDI_resid"], alternative="less")
w_obs, p_obs = stats.wilcoxon(gpr["d_MDI"], alternative="less")
note(f"   组内 Wilcoxon（单侧）: d_res 中位={gpr['d_MDI_resid'].median():+.3f} p={p_res:.4g} | "
     f"d_obs 中位={gpr['d_MDI'].median():+.3f} p={p_obs:.4g}（锚点：注册前披露 ΔMDI≈−0.33, p=0.016）")
note(f"   组成同报: d_Neut 中位={gpr['d_Neut'].median():+.3f} p={stats.wilcoxon(gpr['d_Neut']).pvalue:.4g}")
# H4 归因（复现层）
X = np.column_stack([np.ones(len(gpr)), gpr["d_MDI_comp"].values])
y = gpr["d_MDI"].values
beta, *_ = np.linalg.lstsq(X, y, rcond=None)
share_g = gpr["d_MDI_comp"].mean() / gpr["d_MDI"].mean() if gpr["d_MDI"].mean() != 0 else np.nan
adj_g = ("组成回退" if (gpr["d_MDI_comp"].mean() < 0 and gpr["d_MDI"].mean() < 0 and share_g >= 0.5)
         else "状态可逆")
note(f"   GSE148871 归因: 同向份额={share_g:+.3f} → {adj_g}")
repl = dict(cohort="GSE148871", arm="NEMI", n=int(len(gpr)),
            median_d_res=float(gpr["d_MDI_resid"].median()),
            wilcoxon_p_one=float(p_res), median_d_obs=float(gpr["d_MDI"].median()),
            wilcoxon_p_obs=float(p_obs), p_decline_res=float((gpr["d_MDI_resid"] < 0).mean()),
            share_comp=float(share_g) if np.isfinite(share_g) else np.nan,
            adjudication=adj_g, note="披露的注册前观察，不进 BH 族")

# ======================================================================
# 4. 零模型（双臂同替换，表达量匹配，B=10,000，seed=0）
# ======================================================================
note("\n== 4. 零模型（B=10,000，seed=0）")
# GSE106878 基因级矩阵（step1 同口径，缓存）
gm_cache = INTER + r"\M11M12_step1_genemat_GSE106878.npz"
if os.path.exists(gm_cache):
    zz = np.load(gm_cache, allow_pickle=True)
    gm = pd.DataFrame(zz["M"], index=zz["samples"], columns=zz["genes"])
else:
    f = ROOT + r"\00_RAW_DATA\GSE106878_Sepsis_Hydrocortisone_CORTICUS\GSE106878_non-normalized_data.txt.gz"
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
    # sentrix -> gsm 映射（GSE106878 的 Sample_description[1] = sentrix）
    sent2gsm = {}
    txt = open(ROOT + r"\00_RAW_DATA\GSE106878_Sepsis_Hydrocortisone_CORTICUS\GSE106878_gsm_metadata.txt",
               encoding="utf-8", errors="replace").read()
    for rec in txt.split("^SAMPLE")[1:]:
        gsm = None; descs = []
        for ln in rec.split("\n"):
            if ln.startswith("!Sample_geo_accession"):
                gsm = ln.split("=", 1)[1].strip()
            elif ln.startswith("!Sample_description"):
                descs.append(ln.split("=", 1)[1].strip())
        if gsm and len(descs) > 1:
            sent2gsm[descs[1]] = gsm
        elif gsm and descs:
            sent2gsm[descs[0]] = gsm
    expr.index = [sent2gsm.get(str(i), str(i)) for i in expr.index]
    p2s = {}
    with open(ROOT + r"\00_RAW_DATA\GPL10295_table.txt", encoding="utf-8", errors="replace") as ff:
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
    gmx = pd.DataFrame(index=expr.index)
    for gg, ps in g2p.items():
        subp = expr[ps]
        if len(ps) > 3:
            best = subp.mean(axis=0).nlargest(3).index
            subp = subp[best]
        gmx[gg] = subp.mean(axis=1)
    sid = list(S[S["cohort"] == "GSE106878"]["sample_id"])
    gm = gmx.loc[[i for i in sid if i in gmx.index]]
    np.savez_compressed(gm_cache, M=gm.values.astype(np.float32),
                        samples=np.array(gm.index, dtype=object),
                        genes=np.array(gm.columns, dtype=object))
note(f"   基因矩阵 {gm.shape[0]} 样本 x {gm.shape[1]} 基因")
# MDI_comp 固定（step1 口径）
mdicomp = S[S["cohort"] == "GSE106878"].set_index("sample_id")["MDI_comp"]
# 观测 MWU z 已算（z_mwu）
# 零模型：双臂同替换
universe = [g for g in gm.columns if g.upper() not in ALL80]
means = gm[universe].mean(axis=0).values
deciles, bins = pd.qcut(means, 10, labels=False, duplicates="drop", retbins=True)
idx_by_dec = {int(d): np.where(deciles == d)[0] for d in np.unique(deciles)}
up_in = [g for g in UP if g in gm.columns]
ex_in = [g for g in EX if g in gm.columns]
up_mean = gm[up_in].mean(axis=0).values
ex_mean = gm[ex_in].mean(axis=0).values
if len(bins) >= 3:
    up_dec = np.clip(np.digitize(up_mean, bins[1:-1]), 0, len(bins) - 2).astype(int)
    ex_dec = np.clip(np.digitize(ex_mean, bins[1:-1]), 0, len(bins) - 2).astype(int)
else:
    up_dec = np.zeros(len(up_in), dtype=int)
    ex_dec = np.zeros(len(ex_in), dtype=int)
# 患者臂映射
pat_arm = dict(zip(pairs["patient"], pairs["treatment"]))
sid2pat = dict(zip(d["sample_id"], d["patient"]))
sid2tp = dict(zip(d["sample_id"], d["tp"]))
# 观测统计量（z_mwu）
null_z = np.empty(B_NULL)
for b in range(B_NULL):
    sel_up = [universe[idx_by_dec[ad][RNG.integers(0, len(idx_by_dec[ad]))]] for ad in up_dec]
    sel_ex = [universe[idx_by_dec[ad][RNG.integers(0, len(idx_by_dec[ad]))]] for ad in ex_dec]
    ucs_n = zvec(gm[sel_up].mean(axis=1).values)
    eis_n = zvec(gm[sel_ex].mean(axis=1).values)
    mdi_n = eis_n - ucs_n
    # 残差化（refit OLS MDI_null ~ MDI_comp）
    comp_vals = mdicomp.reindex(gm.index).values
    A = np.column_stack([np.ones(len(mdi_n)), comp_vals])
    beta_n, *_ = np.linalg.lstsq(A, mdi_n, rcond=None)
    resid_n = mdi_n - A @ beta_n
    # d_res per patient
    res_s = pd.Series(resid_n, index=gm.index)
    pat_rows = {}
    for sid_, v in res_s.items():
        p_ = sid2pat.get(sid_); tp_ = sid2tp.get(sid_)
        if p_ and tp_:
            pat_rows.setdefault(p_, {})[tp_] = v
    hc_v = []; pl_v = []
    for p_, vals in pat_rows.items():
        if "Pre" in vals and "Post(24h)" in vals:
            (hc_v if pat_arm[p_] == "hydrocortisone" else pl_v).append(vals["Post(24h)"] - vals["Pre"])
    if len(hc_v) < 2 or len(pl_v) < 2:
        null_z[b] = 0.0
        continue
    u_b, _ = stats.mannwhitneyu(hc_v, pl_v, alternative="less")
    n1b, n2b = len(hc_v), len(pl_v)
    null_z[b] = (u_b - n1b * n2b / 2) / np.sqrt(n1b * n2b * (n1b + n2b + 1) / 12)
p_emp = float((null_z <= z_mwu).mean())
note(f"   零模型: obs z={z_mwu:+.3f} null z={null_z.mean():+.3f}±{null_z.std():.3f} 经验p={p_emp:.4g}")
zero = pd.DataFrame(dict(n_arm_genes=[len(up_in) + len(ex_in)],
                         obs_z=[float(z_mwu)], null_mean=[float(null_z.mean())],
                         null_sd=[float(null_z.std())], empirical_p=[p_emp]))

# ======================================================================
# 5. 判定门 M12
# ======================================================================
note("\n== 5. 判定门 M12")
repl_same_direction = (repl["median_d_res"] < 0)
gate = dict(confirmatory_between_arm_BH=bool(q[0] < 0.05),
            within_arm_co_report_BH=bool(q[1] < 0.05),
            loso_mwu=bool(loso_frac_mwu >= 0.9),
            loso_wilcoxon=bool(loso_frac_wil >= 0.9),
            replication_gse148871_same_direction=repl_same_direction,
            mesi_p_decline=float(p_decline))
if gate["confirmatory_between_arm_BH"] and gate["replication_gse148871_same_direction"]:
    verdict = "M12 PASS（判据 iii 成立）"
elif not gate["confirmatory_between_arm_BH"]:
    verdict = "M12 阴性：GSE106878 确认性检验未过 BH → 按注册条款如实阴性入稿"
else:
    verdict = "M12 未全过门（复现层不同向）→ 如实入稿"
note(f"   判定: {verdict}")

# ======================================================================
# 6. 导出
# ======================================================================
note("\n== 6. 导出")
arm_tests = pd.DataFrame([
    dict(test="between_arm_MWU_one_sided_less", n_hc=n1, n_pl=n2, u=float(u_stat),
         z=float(z_mwu), p=float(p_mwu), q_bh=float(q[0]), cliff_delta=float(cliff_delta)),
    dict(test="within_arm_Wilcoxon_HC", n=len(hc), stat=float(w_stat), p=float(p_wil),
         q_bh=float(q[1]), p_decline=p_decline),
    dict(test="within_arm_Wilcoxon_PL_descriptive", n=len(pl),
         p_two=float(stats.wilcoxon(pl["d_MDI_resid"]).pvalue),
         median=float(pl["d_MDI_resid"].median())),
])
arm_tests.to_csv(INTER + r"\M11M12_step3_gse106878_arm_tests.csv", index=False)
pairs[["patient", "treatment"] + [f"d_{c}" for c in ("MDI", "MDI_resid", "MDI_comp", "UCS", "EIS", "Neut", "Mono", "comp")]].to_csv(
    INTER + r"\M11M12_step3_gse106878_per_patient.csv", index=False)
comp.to_csv(INTER + r"\M11M12_step3_composition_coreport.csv", index=False)
h4.to_csv(INTER + r"\M11M12_step3_h4_attribution.csv", index=False)
pd.DataFrame([repl]).to_csv(INTER + r"\M11M12_step3_gse148871_replication.csv", index=False)
zero.to_csv(INTER + r"\M11M12_step3_zero_m12.csv", index=False)
pd.DataFrame([{**gate, "verdict": verdict}]).to_csv(INTER + r"\M11M12_step3_gate.csv", index=False)
gpr[["patient"] + [f"d_{c}" for c in ("MDI", "MDI_resid", "MDI_comp", "Neut", "Mono", "comp")]].to_csv(
    INTER + r"\M11M12_step3_gse148871_per_patient.csv", index=False)
note("   已导出 step3 全部结果表")

with open(LOG + r"\M11M12_step3_log.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(log))
note(f"\nDONE in {time.time()-t0:.1f}s -> {LOG}\\M11M12_step3_log.txt")
