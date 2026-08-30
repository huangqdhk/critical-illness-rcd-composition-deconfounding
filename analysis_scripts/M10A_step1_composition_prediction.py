# -*- coding: utf-8 -*-
"""
M10A_step1_composition_prediction.py — M10 检验 A：组成预测值 vs 观测值（主终点）
================================================================================
前瞻注册：M10_M13_M14_pre_registration_20260824.md（osf.io/ETVMJ）§1.2
- 参考谱：口径 A = Monaco GSE107011 29 种免疫细胞 TPM（含中性粒细胞，log2(TPM+1) 细胞型均值）；
  口径 B = ABIS sigmatrixRNAseq 17 细胞型（log2(+1)）。
- 反卷积：NNLS（scipy 主口径）；方法交叉 = OLS 后置零截断（离线环境无法运行 CIBERSORTx，
  以 OLS-CLS 替代，如实披露偏离注册）。两口径 x 两方法 = 4 套比例估计，全部报告。
- 合成谱：X_comp(s) = sum_c pi_c(s)*mu_c（log2 空间）；对合成谱按 mdi_v1.0 冻结定义
  （臂均值 -> 队列内全样本 z，不分组）计算组成预测 UCS/EIS/MDI。
- 主检验 1：真实 MDI(s) 对 MDI_comp(s) 线性回归 R^2（队列内逐队列 + 跨队列合并），
  预注册阈值 R^2 >= 0.5 判组成解释成立。
- 主检验 2：回归残差 MDI（观测 MDI 减去 b0+b1*MDI_comp 拟合值）的疾病效应
  （Hedges' g），在来源互斥层（k=16，P0 同款层定义）REML+Hartung-Knapp 随机效应合并
  + 预测区间 + LOSO；报告层向一致性占比。
- 判定门（预注册）：A-FAIL（Track A）= R^2>=0.5 且校正后疾病效应不显著；
  A-PASS（Track B）= R^2<0.5 且残差效应 >=70% 层方向一致且合并显著。中间情形如实报告。
输出：
  _intermediate/M10A_monaco_ct_means.csv（参考谱细胞型均值，log2 空间）
  _intermediate/M10A_proportions.csv（4 套比例估计长表）
  _intermediate/M10A_synthetic_scores.csv（逐样本组成预测 UCS/EIS/MDI，4 套）
  _intermediate/M10A_r2_comparison.csv（逐队列 R^2 + 合并）
  _intermediate/M10A_residual_layer_effects.csv（残差层效应）
  _intermediate/M10A_residual_meta.csv（残差合并 + LOSO）
  04_AUDIT_GOVERNANCE/M10A_Composition_Test_Report.md
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
RAW = ROOT + r"\00_RAW_DATA"
INTER = ROOT + r"\_intermediate"
GOV = ROOT + r"\04_AUDIT_GOVERNANCE"
M2CSV = INTER + r"\M2_per_sample_scores.csv"
MAN = GOV + r"\SAMPLE_MANIFEST_v1.0.csv"
MONACO_F = RAW + r"\GSE107011_Monaco_ImmuneRef\GSE107011_Processed_data_TPM.txt.gz"
ABIS_F = RAW + r"\ABIS\sigmatrixRNAseq.txt"

sys.path.insert(0, ROOT)
import mdi_lib as L

t0 = time.time()
log = []
def note(m=""):
    log.append(m); print(m, flush=True)

man_df, arms, alias = L.load_manifest()
UP, EX = arms["upstream_collapse"], arms["execution_induction"]

# ---------------- meta 数学（P0 同款 REML + Hartung-Knapp） ----------------
from scipy import stats
from scipy.optimize import brentq
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
                pi_lo=float(pi_lo), pi_hi=float(pi_hi))

# ======================================================================
# 1. Monaco 参考谱（29 细胞型，log2(TPM+1)）
# ======================================================================
note("== 1. Monaco 参考谱")
mon = pd.read_csv(MONACO_F, sep="\t", index_col=0)
note(f"   TPM 矩阵 {mon.shape[0]} 基因 x {mon.shape[1]} 样本")
cols = [str(c) for c in mon.columns]
def parse_ct(c):
    return c.split("_", 1)[1] if "_" in c else c
ctypes = sorted(set(parse_ct(c) for c in cols))
ctypes = [c for c in ctypes if c != "PBMC"]   # 排除未分选 PBMC 混合样本列（README 29 型口径）
note(f"   细胞类型 {len(ctypes)}（已排除 PBMC 混合列）: {ctypes}")
mon.columns = [parse_ct(c) for c in cols]
mon = mon[[c for c in mon.columns if c in ctypes]]
# 同类型多供者取均值
mu_tpm = mon.groupby(level=0, axis=1).mean()
# ENSG（去版本）-> symbol
cache = INTER + r"\M10A_monaco_ensg2sym.csv"
ensg_raw = [str(i).split(".")[0] for i in mu_tpm.index]
if os.path.exists(cache):
    e2s = dict(pd.read_csv(cache).values)
    e2s = {str(k).upper(): str(v).upper() for k, v in e2s.items()}
else:
    import mygene
    mg = mygene.MyGeneInfo()
    res = mg.querymany(ensg_raw, scopes="ensembl.gene", fields="symbol", species="human",
                       returnall=True, verbose=False)
    e2s = {}
    for r in res.get("out", []):
        if r.get("symbol"):
            e2s[str(r["query"]).upper()] = str(r["symbol"]).upper()
    pd.DataFrame({"ensg": list(e2s.keys()), "symbol": list(e2s.values())}).to_csv(cache, index=False)
note(f"   ENSG->symbol 映射 {len(e2s)}/{len(ensg_raw)}")
mu = mu_tpm.copy()
mu["sym"] = [e2s.get(e, "") for e in ensg_raw]
mu = mu[mu["sym"] != ""]
mu = mu.groupby("sym").mean()          # 同符号取均值
mu = mu.T                                # 细胞型 x 基因
note(f"   参考谱（符号化）: {mu.shape[0]} 细胞型 x {mu.shape[1]} 基因")
mu_log = np.log2(mu + 1.0)
mu_log.to_csv(INTER + r"\M10A_monaco_ct_means.csv")
UP_MON = [g for g in UP if g in mu.columns]
EX_MON = [g for g in EX if g in mu.columns]
note(f"   manifest 臂覆盖：上游 {len(UP_MON)}/{len(UP)}，执行 {len(EX_MON)}/{len(EX)}")

# ======================================================================
# 2. ABIS 参考谱（17 细胞型）
# ======================================================================
note("\n== 2. ABIS 参考谱")
abis = pd.read_csv(ABIS_F, sep="\t", index_col=0)
abis.index = [str(i).strip().strip('"\'').upper() for i in abis.index]
abis.columns = [str(c).strip().strip('"\'') for c in abis.columns]
abis = abis.groupby(level=0).mean()
abis_log = np.log2(abis + 1.0).T        # 细胞型 x 基因
UP_ABIS = [g for g in UP if g in abis_log.columns]
EX_ABIS = [g for g in EX if g in abis_log.columns]
note(f"   ABIS {abis_log.shape[0]} 细胞型 x {abis_log.shape[1]} 基因；"
     f"臂覆盖：上游 {len(UP_ABIS)}/{len(UP)}，执行 {len(EX_ABIS)}/{len(EX)}")

# ======================================================================
# 3. 队列矩阵加载（样本 x 基因 log2 空间，符号列）
# ======================================================================
note("\n== 3. 队列矩阵加载")
def fast_collapse(probe_mat, probe2sym_csv, wanted_symbols, max_probes=3):
    """probe_mat: samples x probes(列=probe id)；返回 samples x genes（max-mean <=3 探针）。"""
    pmap = pd.read_csv(probe2sym_csv)
    pmap.columns = [str(c).upper().strip() for c in pmap.columns]
    pc = [c for c in pmap.columns if c in ("PROBEID", "ID", "PROBE_ID")][0]
    sc = [c for c in pmap.columns if c in ("SYMBOL", "GENE_SYMBOL", "GENE.SYMBOL")][0]
    pmap = pmap[[pc, sc]]
    pmap[pc] = pmap[pc].astype(str).str.strip()
    pmap[sc] = pmap[sc].astype(str).str.strip().str.upper()
    g2p = {}
    for pid, sym in pmap.itertuples(index=False, name=None):
        if sym in wanted_symbols:
            g2p.setdefault(sym, []).append(pid)
    out = pd.DataFrame(index=probe_mat.index)
    for g, probes in g2p.items():
        ps = [p for p in probes if p in probe_mat.columns]
        if not ps:
            continue
        sub = probe_mat[ps]
        if len(ps) > max_probes:
            best = sub.mean(axis=0).nlargest(max_probes).index
            sub = sub[best]
        out[g] = sub.mean(axis=1)
    return out

# 3.1 GSE185263
cnt = pd.read_csv(RAW + r"\GSE185263_Lung_ARDS\GSE185263_raw_counts.csv", index_col=0)
cnt.index = [str(i).upper() for i in cnt.index]
cnt.columns = [str(c) for c in cnt.columns]
cache2 = INTER + r"\P2_gse185263_en2sym.csv"
e2s2 = dict(pd.read_csv(cache2).values)
e2s2 = {str(k).upper(): str(v).upper() for k, v in e2s2.items()}
sub = cnt.loc[[e for e in e2s2 if e in cnt.index]].copy()
sub.index = [e2s2[e] for e in sub.index]
sub = sub.groupby(level=0).mean()
cpm = sub / sub.sum(axis=0) * 1e6
g185 = np.log2(cpm + 1).T
g185_groups = pd.read_csv(GOV + r"\GSE185263_groups.csv")
g185_groups.columns = ["sample_id", "group"]
g185_groups["sample_id"] = g185_groups["sample_id"].astype(str)
note(f"   GSE185263: {g185.shape[0]} 样本 x {g185.shape[1]} 基因")

# 3.2 GSE32707
m7 = pd.read_csv(RAW + r"\GSE32707_data\GSE32707_gene_symbol_log2_matrix.csv.gz", index_col=0)
m7.index = [str(i).upper() for i in m7.index]
m7.columns = [str(c) for c in m7.columns]
g327 = m7.groupby(level=0).mean().T
note(f"   GSE32707: {g327.shape[0]} 样本 x {g327.shape[1]} 基因")

# 3.3 GSE212865
meta2, s2, m2_ = L.parse_series_matrix(RAW + r"\GSE212865_data\GSE212865_series_matrix.txt.gz")
g212 = fast_collapse(m2_, RAW + r"\GPL23159_probe2symbol_clariomsdb.csv", set(mu.columns))
note(f"   GSE212865: {g212.shape[0]} 样本 x {g212.shape[1]} 基因（max-mean 塌陷）")

# 3.4 GSE188309
meta8, s8, m8_ = L.parse_series_matrix(RAW + r"\GSE188309_CAP_WholeBlood\GSE188309_series_matrix.txt.gz")
g188 = fast_collapse(m8_, RAW + r"\GPL23159_probe2symbol_clariomsdb.csv", set(mu.columns))
note(f"   GSE188309: {g188.shape[0]} 样本 x {g188.shape[1]} 基因")

# 3.5 GSE310929
tsv = RAW + r"\GSE310929_Sepsis\GSE310929_AllSampleExpressionSubmitted.tsv\GSE310929_AllSampleExpressionSubmitted.tsv"
m3 = pd.read_csv(tsv, sep="\t", index_col=0)
m3.index = [str(i).strip('"').upper() for i in m3.index]
m3.columns = [str(c).strip('"') for c in m3.columns]
m3 = m3.groupby(level=0).mean().T
want = [g for g in mu.columns if g in m3.columns]
g310 = m3[want]
note(f"   GSE310929: {g310.shape[0]} 样本 x {g310.shape[1]} 基因（提交方分位数标准化+ComBat 矩阵，NNLS 近似性披露）")

# 3.6 GSE148871（血样本，R^2 敏感性用）
meta4, s4, m4_ = L.parse_series_matrix(RAW + r"\GSE148871_COPD_AE\GSE148871_series_matrix.txt.gz")
g148 = fast_collapse(m4_, RAW + r"\GPL570_probe2symbol_hgu133plus2db.csv", set(mu.columns))
tissue4 = L.sample_characteristics(meta4, label="tissue")
s4t = pd.Series(tissue4, index=s4)
blood4 = [s for s in s4 if "blood" in str(s4t.get(s, "")).lower()]
g148b = g148.loc[[s for s in blood4 if s in g148.index]]
note(f"   GSE148871 blood: {g148b.shape[0]} 样本 x {g148b.shape[1]} 基因")

# ======================================================================
# 4. 反卷积（4 套：Monaco/ABIS x NNLS/OLS-CLS）
# ======================================================================
note("\n== 4. 反卷积与组成预测评分")
def zc(s):
    s = np.asarray(s, float); sd = s.std(ddof=1)
    return (s - s.mean()) / sd if sd > 0 else np.zeros_like(s)

PROPS_F = INTER + r"\M10A_proportions.csv"
SYNTH_F = INTER + r"\M10A_synthetic_scores.csv"
cohorts = [("GSE185263", g185), ("GSE32707", g327), ("GSE212865", g212),
           ("GSE188309", g188), ("GSE310929", g310), ("GSE148871_blood", g148b)]
refs = [("Monaco", mu_log, UP_MON, EX_MON), ("ABIS", abis_log, UP_ABIS, EX_ABIS)]
methods = ["NNLS", "OLSCLS"]
if os.path.exists(PROPS_F) and os.path.exists(SYNTH_F):
    props = pd.read_csv(PROPS_F)
    synth = pd.read_csv(SYNTH_F)
    props["sample"] = props["sample"].astype(str)
    synth["sample"] = synth["sample"].astype(str)
    note("   [缓存复用] 读取既有 M10A_proportions.csv / M10A_synthetic_scores.csv")
else:
    def deconvolve(X, ref_log, method):
        """X: samples x genes(log2)；ref_log: celltypes x genes(log2)。返回 proportions (samples x celltypes)。"""
        common = [g for g in X.columns if g in ref_log.columns]
        A = ref_log[common].values.T        # genes x celltypes
        Xc = X[common].values               # samples x genes
        AtA = A.T @ A
        n_ct = A.shape[1]
        rows = []
        if method == "NNLS":
            from scipy.linalg import cholesky, solve_triangular
            Lc = cholesky(AtA + np.eye(n_ct) * 1e-10, lower=True)
            from scipy.optimize import nnls
            for i in range(Xc.shape[0]):
                Atx = A.T @ Xc[i]
                rhs = solve_triangular(Lc, Atx, lower=True)
                pi, _ = nnls(Lc.T, rhs)
                s = pi.sum()
                rows.append((pi / s) if s > 0 else np.full(n_ct, 1.0 / n_ct))
        else:  # OLS-CLS
            for i in range(Xc.shape[0]):
                pi = np.linalg.lstsq(AtA + np.eye(n_ct) * 1e-6, A.T @ Xc[i], rcond=None)[0]
                pi = np.clip(pi, 0, None)
                s = pi.sum()
                rows.append((pi / s) if s > 0 else np.full(n_ct, 1.0 / n_ct))
        return pd.DataFrame(rows, index=X.index, columns=ref_log.index), common

    def synthetic_arm_scores(prop, ref_log, up_genes, ex_genes):
        """X_comp = pi @ mu；臂均值（log2 空间）-> 队列内 z -> MDI_comp。"""
        up_in = [g for g in up_genes if g in ref_log.columns]
        ex_in = [g for g in ex_genes if g in ref_log.columns]
        mu_up = ref_log[up_in].values          # celltypes x genes
        mu_ex = ref_log[ex_in].values
        P = prop.values                        # samples x celltypes
        uc = P @ mu_up.mean(axis=1)
        ec = P @ mu_ex.mean(axis=1)
        out = pd.DataFrame(index=prop.index)
        out["UCS_comp_raw"] = uc; out["EIS_comp_raw"] = ec
        out["UCS_comp"] = zc(uc); out["EIS_comp"] = zc(ec)
        out["MDI_comp"] = out["EIS_comp"] - out["UCS_comp"]
        out["n_up_ref"] = len(up_in); out["n_ex_ref"] = len(ex_in)
        return out

    prop_rows, synth_rows = [], []
    for cname, X in cohorts:
        for rname, ref_log, up_in, ex_in in refs:
            for method in methods:
                key = f"{cname}|{rname}|{method}"
                prop, common = deconvolve(X, ref_log, method)
                prop_rows.append(prop.assign(cohort=cname, reference=rname, method=method))
                sc = synthetic_arm_scores(prop, ref_log, up_in, ex_in)
                sc = sc.assign(cohort=cname, reference=rname, method=method)
                synth_rows.append(sc)
                note(f"   {key:<34} n_common={len(common):<5} prop_sum_range="
                     f"[{prop.sum(axis=1).min():.3f},{prop.sum(axis=1).max():.3f}] "
                     f"arm_cov=({sc['n_up_ref'].iloc[0]},{sc['n_ex_ref'].iloc[0]})")

    props = pd.concat(prop_rows).reset_index(names="sample")
    synth = pd.concat(synth_rows).reset_index(names="sample")
    props.to_csv(PROPS_F, index=False)
    synth.to_csv(SYNTH_F, index=False)

# ======================================================================
# 5. R^2：观测（M2 冻结 mdi_v1.0）vs 组成预测（逐队列 + 跨队列合并）
# ======================================================================
note("\n== 5. R^2 对比")
obs = pd.read_csv(M2CSV, low_memory=False)
obs["sample"] = obs["sample"].astype(str)

r2_rows = []
for cname, X in cohorts:
    o = obs[obs["cohort"] == ("GSE148871" if cname == "GSE148871_blood" else cname)]
    for rname, _, _, _ in refs:
        for method in methods:
            sc = synth[(synth["cohort"] == cname) & (synth["reference"] == rname) & (synth["method"] == method)]
            sc = sc.set_index("sample")
            # 对齐样本
            common_samples = sorted(set(o["sample"]) & set(sc.index))
            if len(common_samples) < 10:
                continue
            oo = o.set_index("sample").loc[common_samples]
            ss = sc.loc[common_samples]
            for fam in ("UCS", "EIS", "MDI"):
                r = np.corrcoef(oo[fam].values, ss[f"{fam}_comp"].values)[0, 1]
                r2_rows.append(dict(cohort=cname, reference=rname, method=method,
                                    family=fam, n=len(common_samples), r=r, R2=r ** 2))
r2 = pd.DataFrame(r2_rows)
r2.to_csv(INTER + r"\M10A_r2_comparison.csv", index=False)
note("   逐队列 R^2（Monaco|NNLS 主口径，MDI 族）：")
main_r2 = r2[(r2["reference"] == "Monaco") & (r2["method"] == "NNLS") & (r2["family"] == "MDI")]
for _, r in main_r2.iterrows():
    note(f"     {r['cohort']:<18} n={r['n']:<5} r={r['r']:+.3f} R^2={r['R2']:.3f}")

# 跨队列合并（Fisher z meta + 合并样本相关系数）
note("\n   跨队列合并（Monaco|NNLS 主口径）:")
merged_rows = []
for fam in ("UCS", "EIS", "MDI"):
    sub = r2[(r2["reference"] == "Monaco") & (r2["method"] == "NNLS") & (r2["family"] == fam)]
    zs = np.arctanh(np.clip(sub["r"].values, -0.999, 0.999))
    z_pool = zs.mean()
    r_pool = np.tanh(z_pool)
    # 合并样本口径
    xs, ys = [], []
    for cname, X in cohorts:
        o = obs[obs["cohort"] == ("GSE148871" if cname == "GSE148871_blood" else cname)].set_index("sample")
        sc = synth[(synth["cohort"] == cname) & (synth["reference"] == "Monaco") & (synth["method"] == "NNLS")].set_index("sample")
        cmn = sorted(set(o.index) & set(sc.index))
        xs.append(o.loc[cmn, fam].values); ys.append(sc.loc[cmn, f"{fam}_comp"].values)
    r_all = np.corrcoef(np.concatenate(xs), np.concatenate(ys))[0, 1]
    merged_rows.append(dict(family=fam, fisher_z_meta_r=r_pool, n_cohorts=len(sub),
                            pooled_sample_r=r_all, pooled_sample_R2=r_all ** 2))
    note(f"     {fam}: Fisher-z meta r={r_pool:+.3f}（k={len(sub)}）; 合并样本 r={r_all:+.3f}, R^2={r_all**2:.3f}")
pd.DataFrame(merged_rows).to_csv(INTER + r"\M10A_r2_merged.csv", index=False)

# 两参考口径比例交叉验证（Monaco-NNLS vs ABIS-NNLS，匹配细胞型 Spearman）
note("\n   两参考口径比例交叉验证（Monaco-NNLS vs ABIS-NNLS）:")
TYPE_MAP = {"Monocytes C": "C_mono", "NK": "NK", "T CD8 Memory": "CD8_CM",
            "T CD8 Naive": "CD8_naive", "T CD4 Naive": "CD4_naive", "T CD4 Memory": "CD4_TE",
            "B Memory": "B_SM", "B Naive": "B_naive", "Neutrophils LD": "Neutrophils",
            "Basophils LD": "Basophils", "pDCs": "pDC", "mDCs": "mDC",
            "Plasmablasts": "Plasmablasts", "MAIT": "MAIT", "T gd Vd2": "VD2+",
            "T gd non-Vd2": "VD2-"}
cross_rows = []
for cname, X in cohorts:
    pm = props[(props["cohort"] == cname) & (props["reference"] == "Monaco") & (props["method"] == "NNLS")]
    pa = props[(props["cohort"] == cname) & (props["reference"] == "ABIS") & (props["method"] == "NNLS")]
    if len(pm) == 0 or len(pa) == 0:
        continue
    pm = pm.set_index("sample"); pa = pa.set_index("sample")
    cmn = sorted(set(pm.index) & set(pa.index))
    rs_list = []
    for a_t, m_t in TYPE_MAP.items():
        if a_t in pa.columns and m_t in pm.columns:
            r = np.corrcoef(pa.loc[cmn, a_t], pm.loc[cmn, m_t])[0, 1]
            rs_list.append((m_t, r))
    cross_rows.append(dict(cohort=cname, n_samples=len(cmn), n_types=len(rs_list),
                           median_r=float(np.median([r for _, r in rs_list])),
                           min_r=float(np.min([r for _, r in rs_list]))))
    note(f"     {cname:<18} n={len(cmn):<5} 匹配类型 {len(rs_list)} 中位 r={np.median([r for _, r in rs_list]):+.3f}")
pd.DataFrame(cross_rows).to_csv(INTER + r"\M10A_proportion_cross_reference.csv", index=False)

# ======================================================================
# 6. 主检验 2：残差 MDI 疾病效应（来源互斥层 k=16，P0 同款层定义）
# ======================================================================
note("\n== 6. 残差 MDI 疾病效应（来源互斥层）")
sc_main = synth[(synth["reference"] == "Monaco") & (synth["method"] == "NNLS")].set_index("sample")
# 按队列做 OLS MDI ~ MDI_comp，得残差
resid_parts = []
for cname, X in cohorts:
    cohort_name = "GSE148871" if cname == "GSE148871_blood" else cname
    o = obs[obs["cohort"] == cohort_name].set_index("sample")
    sc = sc_main[sc_main["cohort"] == cname]
    cmn = sorted(set(o.index) & set(sc.index))
    oo = o.loc[cmn]; ss = sc.loc[cmn]
    if len(cmn) < 10:
        continue
    A = np.column_stack([np.ones(len(cmn)), ss["MDI_comp"].values])
    beta, *_ = np.linalg.lstsq(A, oo["MDI"].values, rcond=None)
    resid = oo["MDI"].values - A @ beta
    rdf = pd.DataFrame({"sample": cmn, "MDI_obs": oo["MDI"].values,
                        "MDI_comp": ss["MDI_comp"].values, "MDI_resid": resid,
                        "b0": beta[0], "b1": beta[1]})
    resid_parts.append(rdf)
resid_all = pd.concat(resid_parts)
resid_all.to_csv(INTER + r"\M10A_residual_per_sample.csv", index=False)

obs_full = obs.copy()
obs_full["sample"] = obs_full["sample"].astype(str)
res_map = resid_all.set_index("sample")["MDI_resid"].to_dict()
obs_full["MDI_resid"] = obs_full["sample"].map(res_map)

CASES_310 = ["Sepsis", "Sepsis - Shock", "Sepsis - ARDs"]
MIN_N = 3
strata = []
note("   层装配：")
g3 = obs_full[obs_full["cohort"] == "GSE310929"]
for ds, s in g3.groupby("Dataset"):
    if ds in ("GSE185263", "GSE32707"):
        continue
    case = s[s["group"].isin(CASES_310)]
    ctrl = s[s["group"] == "Control"]
    if len(case) >= MIN_N and len(ctrl) >= MIN_N:
        strata.append((f"310:{ds}", case, ctrl))
for cohort, cgrp in [("GSE185263", ["Sepsis_COVID"]), ("GSE32707", ["ARDS_d0"])]:
    s = obs_full[obs_full["cohort"] == cohort]
    strata.append((cohort, s[s["group"].isin(cgrp)], s[s["group"] == "Control"]))
s2x = obs_full[obs_full["cohort"] == "GSE212865"].merge(
    pd.read_csv(MAN)[lambda d: d["dataset_id"] == "GSE212865"][["gsm", "timepoint", "subcohort"]],
    left_on="sample", right_on="gsm", how="left")
s2x = s2x[s2x["timepoint"].isna() | (s2x["timepoint"].astype(str) == "D0")]
dup = s2x.duplicated(subset=["subcohort", "group"], keep=False) & s2x["subcohort"].notna()
if dup.any():
    s2x = s2x[~dup | ~s2x.duplicated(subset=["subcohort", "group"], keep="first")]
strata.append(("GSE212865_base", s2x[s2x["group"] == "Covid19_SDRA"], s2x[s2x["group"] == "Control"]))

note(f"   可合并层 {len(strata)} 个（含残差可得层）")
layer_rows = []
for sid, case, ctrl in strata:
    case_r = case["MDI_resid"].dropna().values
    ctrl_r = ctrl["MDI_resid"].dropna().values
    if len(case_r) < MIN_N or len(ctrl_r) < MIN_N:
        continue
    g, gse = L.hedges_g(ctrl_r, case_r)
    # 观测值同层效应（对照锚点，应与 P0 一致）
    g_obs, gse_obs = L.hedges_g(ctrl["MDI"].values, case["MDI"].values)
    p_mw = stats.mannwhitneyu(case_r, ctrl_r, alternative="two-sided").pvalue
    layer_rows.append(dict(stratum=sid, n_case=len(case_r), n_ctrl=len(ctrl_r),
                           resid_g=g, resid_g_se=gse, resid_MW_p=p_mw,
                           obs_g=g_obs, obs_g_se=gse_obs))
layers = pd.DataFrame(layer_rows)
layers.to_csv(INTER + r"\M10A_residual_layer_effects.csv", index=False)
note(f"   残差可得层 {len(layers)}；方向（resid_g>0）占比 {int((layers['resid_g']>0).sum())}/{len(layers)} "
     f"({(layers['resid_g']>0).mean():.0%})")

meta_res = pool_hk(layers["resid_g"].values, layers["resid_g_se"].values)
meta_obs = pool_hk(layers["obs_g"].values, layers["obs_g_se"].values)
note(f"   [观测 MDI 同层合并（锚点）] g={meta_obs['pooled_g']:+.3f} "
     f"95%CI=[{meta_obs['ci_lo']:+.3f},{meta_obs['ci_hi']:+.3f}] 单侧p={meta_obs['p_pos']:.4g} I2={meta_obs['I2']:.0f}%")
note(f"   [残差 MDI 合并]         g={meta_res['pooled_g']:+.3f} "
     f"95%CI=[{meta_res['ci_lo']:+.3f},{meta_res['ci_hi']:+.3f}] 单侧p(pos)={meta_res['p_pos']:.4g} "
     f"单侧p(neg)={meta_res['p_neg']:.4g} I2={meta_res['I2']:.0f}% PI=[{meta_res['pi_lo']:+.3f},{meta_res['pi_hi']:+.3f}]")
loso_res = []
full_sign = np.sign(meta_res["pooled_g"])
for i in range(len(layers)):
    keep = layers.drop(layers.index[i])
    r2m = pool_hk(keep["resid_g"].values, keep["resid_g_se"].values)
    loso_res.append(dict(family="MDI_resid", left_out=layers.iloc[i]["stratum"],
                         pooled_g=r2m["pooled_g"], same_sign=bool(np.sign(r2m["pooled_g"]) == full_sign)))
loso_df = pd.DataFrame(loso_res)
loso_df.to_csv(INTER + r"\M10A_residual_meta_loso.csv", index=False)
note(f"   LOSO 同向 {int(loso_df['same_sign'].sum())}/{len(loso_df)}")

# ======================================================================
# 7. 判定与报告
# ======================================================================
main_mdi = r2[(r2["reference"] == "Monaco") & (r2["method"] == "NNLS") & (r2["family"] == "MDI")]
n_R2_ge_half = int((main_mdi["R2"] >= 0.5).sum())
pool_R2 = merged_rows[-1]["pooled_sample_R2"]
res_sig_pos = meta_res["p_pos"] < 0.05
res_sig_neg = meta_res["p_neg"] < 0.05
dir_frac = (layers["resid_g"] > 0).mean()
verdict = ""
if pool_R2 >= 0.5 and not (res_sig_pos and dir_frac >= 0.7):
    verdict = "A-FAIL（R^2>=0.5 且校正后疾病效应不显著/不一致）-> Track A 方向"
elif pool_R2 < 0.5 and res_sig_pos and dir_frac >= 0.7:
    verdict = "A-PASS（R^2<0.5 且残差效应显著、方向一致）-> Track B 方向"
else:
    verdict = "中间情形：按支持方向取权重报告，不强制二分"
note(f"\n== 7. 判定门 M10 检验 A：{verdict}")
note(f"   合并样本 R^2={pool_R2:.3f}（阈值 0.5）；逐队列 R^2>=0.5 的队列 {n_R2_ge_half}/{len(main_mdi)}")
note(f"   残差合并 g={meta_res['pooled_g']:+.3f}，方向一致占比 {dir_frac:.0%}，单侧 p={meta_res['p_pos']:.4g}")

report = []
report.append("# M10 检验 A 报告：组成预测值 vs 观测值\n")
report.append("- 日期：2026-08-26；依据：M10_M13_M14_pre_registration_20260824.md §1.2（osf.io/ETVMJ）")
report.append("- 参考谱：口径 A Monaco（29 细胞型，含中性粒细胞，log2(TPM+1)）；口径 B ABIS（17 细胞型，log2(+1)）")
report.append("- 方法：NNLS 主口径；OLS-CLS 交叉（离线环境无法运行 CIBERSORTx，如实披露偏离注册）；4 套比例估计全部报告")
report.append("- 合成谱按 mdi_v1.0 冻结定义评分（臂均值->队列内全样本 z）\n")
report.append("## 1. 逐队列 R^2（真实 MDI vs 组成预测 MDI；Monaco|NNLS 主口径）\n\n")
report.append("| 队列 | n | r | R^2 |\n|---|---|---|---|\n")
for _, r in main_mdi.iterrows():
    report.append(f"| {r['cohort']} | {r['n']} | {r['r']:+.3f} | {r['R2']:.3f} |\n")
report.append(f"\n合并样本 R^2 = {pool_R2:.3f}（Fisher-z meta 见 M10A_r2_merged.csv）\n")
report.append("\n## 2. 主检验 2：残差 MDI 疾病效应（来源互斥层）\n\n")
report.append(f"- 观测 MDI 同层合并（锚点，应与 P0 g=+1.50 一致）：g={meta_obs['pooled_g']:+.3f} "
              f"[{meta_obs['ci_lo']:+.3f},{meta_obs['ci_hi']:+.3f}] 单侧 p={meta_obs['p_pos']:.4g}\n")
report.append(f"- **残差 MDI 合并：g={meta_res['pooled_g']:+.3f} [{meta_res['ci_lo']:+.3f},{meta_res['ci_hi']:+.3f}]，"
              f"单侧 p(pos)={meta_res['p_pos']:.4g}，p(neg)={meta_res['p_neg']:.4g}，I2={meta_res['I2']:.0f}%，"
              f"预测区间 [{meta_res['pi_lo']:+.3f},{meta_res['pi_hi']:+.3f}]**\n")
report.append(f"- 层向一致性：{int((layers['resid_g']>0).sum())}/{len(layers)} 层为正（{dir_frac:.0%}）\n")
report.append(f"- LOSO 同向：{int(loso_df['same_sign'].sum())}/{len(loso_df)}\n")
report.append("\n## 3. 判定\n\n")
report.append(f"**{verdict}**\n")
report.append("\n（4 套比例估计与全部口径的 R^2 见 M10A_proportions.csv / M10A_r2_comparison.csv）")
with open(GOV + r"\M10A_Composition_Test_Report.md", "w", encoding="utf-8") as f:
    f.write("\n".join(report))
with open(ROOT + r"\03_LOGS\M10A_composition_test_log.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(log))
note(f"\nDONE in {time.time()-t0:.1f}s -> {GOV}\\M10A_Composition_Test_Report.md")
