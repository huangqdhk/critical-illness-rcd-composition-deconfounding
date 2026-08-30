# -*- coding: utf-8 -*-
"""
M13_step1_confounding_audit.py — M13 RCD 去混淆审计（三子任务）
==============================================================
前瞻注册：M10_M13_M14_pre_registration_20260824.md §2（osf.io/ETVMJ）
- 子任务1 目录：MSigDB v7.5.1（HALLMARK/C2.CP-KEGG-REACTOME/C5 GO）PCD/线粒体相关签名 +
  已发表 PCD 签名；逐签名：执行臂占比（外部本体）/上游臂占比（外部本体）/髓系基因占比/
  臂内方向一致性/疾病效应（Hedges' g，GSE185263 Sepsis_COVID vs Control；GSE32707 敏感性）/
  净评分主导臂（与 UCS/EIS 的相关）。
- 子任务2 预注册回归（主检验）：签名疾病效应 g ~ 执行臂占比 + 髓系基因占比（OLS）；
  执行臂系数单侧检验（H1: beta>0）；LOSO 稳定性（>=90% 方向一致）；
  表达量匹配随机等大小基因集零模型 B=10,000（seed=0）给出经验 p；
  去循环化：臂定义主口径 = 外部本体（GO/REACTOME OXPHOS/线粒体被膜 vs 炎症小体/gasdermin），
  敏感性 = 项目 80 基因 manifest 臂。
- 子任务3 反向审计：已发表 mitoxyperilysis 策展基因集 vs IIAMD v1.0（GSE235046 实验锚定）
  Fisher 检验 + baseMean 匹配随机集零模型 B=10,000。
  措辞纪律：只陈述基因集一致性，不评价他人论文结论对错。
输出：_intermediate/M13_catalog.csv / M13_regression.csv / M13_null_model.csv /
      M13_reverse_audit.csv；04_AUDIT_GOVERNANCE/M13_Confounding_Audit_Report.md
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
MSIG = RAW + r"\MSigDB_v7.5.1_Hs"
SEED = 0
B_NULL = 10000

sys.path.insert(0, ROOT)
import mdi_lib as L

t0 = time.time()
log = []
def note(m=""):
    log.append(m); print(m, flush=True)

def zc(s):
    s = np.asarray(s, float); sd = s.std(ddof=1)
    return (s - s.mean()) / sd if sd > 0 else s * 0

# ======================================================================
# 1. 数据：GSE185263 / GSE32707 基因矩阵 + Monaco 髓系 + manifest
# ======================================================================
note("== 1. 数据加载")
cnt = pd.read_csv(RAW + r"\GSE185263_Lung_ARDS\GSE185263_raw_counts.csv", index_col=0)
cnt.index = [str(i).upper() for i in cnt.index]
cnt.columns = [str(c) for c in cnt.columns]
e2s = dict(pd.read_csv(INTER + r"\P2_gse185263_en2sym.csv").values)
e2s = {str(k).upper(): str(v).upper() for k, v in e2s.items()}
sub = cnt.loc[[e for e in e2s if e in cnt.index]].copy()
sub.index = [e2s[e] for e in sub.index]
sub = sub.groupby(level=0).mean()
cpm = sub / sub.sum(axis=0) * 1e6
g185 = np.log2(cpm + 1)                    # genes x samples
grp = pd.read_csv(GOV + r"\GSE185263_groups.csv")
grp.columns = ["sample_id", "group"]
grp["sample_id"] = grp["sample_id"].astype(str)
case185 = grp[grp["group"] == "Sepsis_COVID"]["sample_id"].tolist()
ctrl185 = grp[grp["group"] == "Control"]["sample_id"].tolist()
case185 = [s for s in case185 if s in g185.columns]
ctrl185 = [s for s in ctrl185 if s in g185.columns]
X185 = g185.T                                 # samples x genes
# 每基因 z（跨样本）+ 疾病 log2FC
g_mean_case = X185.loc[case185].mean(axis=0)
g_mean_ctrl = X185.loc[ctrl185].mean(axis=0)
lfc185 = (g_mean_case - g_mean_ctrl)
Z185 = X185.apply(zc, axis=0)                # samples x genes z
note(f"   GSE185263: {X185.shape[0]} 样本 x {X185.shape[1]} 基因（case={len(case185)}, ctrl={len(ctrl185)}）")

m7 = pd.read_csv(RAW + r"\GSE32707_data\GSE32707_gene_symbol_log2_matrix.csv.gz", index_col=0)
m7.index = [str(i).upper() for i in m7.index]
m7.columns = [str(c) for c in m7.columns]
g327 = m7.groupby(level=0).mean().T
meta7, _, _ = L.parse_series_matrix(RAW + r"\GSE32707_data\GSE32707_series_matrix.txt.gz")
gmap = {"untreated": "Control", "SIRS Day 0": "SIRS_d0", "Sepsis Day 0": "Sepsis_d0",
        "Sepsis Day 7": "Sepsis_d7", "se/ARDS Day 0": "ARDS_d0", "se/ARDS Day 7": "ARDS_d7"}
m7df = pd.DataFrame({"gsm": meta7["!Sample_geo_accession"][0],
                     "group": [gmap.get(s, None) for s in meta7["!Sample_source_name_ch1"][0]]}).set_index("gsm")
case327 = m7df[m7df["group"] == "ARDS_d0"].index.tolist()
ctrl327 = m7df[m7df["group"] == "Control"].index.tolist()
Z327 = g327.apply(zc, axis=0)
note(f"   GSE32707: {g327.shape[0]} 样本 x {g327.shape[1]} 基因（ARDS_d0={len(case327)}, Control={len(ctrl327)}）")

monaco = pd.read_csv(INTER + r"\M10A_monaco_ct_means.csv", index_col=0)  # celltypes x genes
MYELOID_TYPES = {"C_mono", "I_mono", "NC_mono", "Neutrophils", "mDC", "pDC", "Basophils"}
myeloid_ct = [c for c in MYELOID_TYPES if c in monaco.index]
argmax_ct = monaco.idxmax(axis=0)
myeloid_genes = set(argmax_ct[argmax_ct.isin(myeloid_ct)].index)
note(f"   髓系基因（Monaco argmax）: {len(myeloid_genes)}")

man, arms, alias = L.load_manifest()
UP_MAN = set(arms["upstream_collapse"])
EX_MAN = set(arms["execution_induction"])

# ======================================================================
# 2. MSigDB 签名选择 + 外部本体臂定义
# ======================================================================
note("\n== 2. MSigDB 加载与选择")
def load_gmt(path):
    sets = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            name, desc, genes = parts[0], parts[1], [g.upper() for g in parts[2:] if g]
            sets[name] = set(genes)
    return sets

gmt_all = {}
for fn in ("HALLMARK.gmt", "C2_CP_KEGG_REACTOME.gmt", "C5_GO_BP_CC.gmt"):
    s = load_gmt(os.path.join(MSIG, fn))
    gmt_all.update(s)
    note(f"   {fn}: {len(s)} 集")
note(f"   合计 {len(gmt_all)} 集")

KW = ["APOPT", "FERROPT", "PYROPT", "NECROPT", "PANOPT", "CUPROPT", "GASDERMIN",
      "INFLAMMASOM", "CASPASE", "CELL_DEATH", "APOPTOTIC",
      "OXIDATIVE_PHOSPHORYLATION", "MITOCHONDRI", "RESPIRATORY_CHAIN",
      "ELECTRON_TRANSPORT", "OXPHOS"]
ARM_UP_NAMES = {"GOCC_MITOCHONDRIAL_ENVELOPE", "GOCC_MITOCHONDRIAL_INNER_MEMBRANE",
                "GOCC_MITOCHONDRIAL_MATRIX", "GOCC_RESPIRASOME",
                "GOBP_OXIDATIVE_PHOSPHORYLATION"}
ARM_EX_NAMES = {"GOBP_PYROPTOSIS", "GOBP_INFLAMMASOME_COMPLEX_ASSEMBLY",
                "GOCC_INFLAMMASOME_COMPLEX", "REACTOME_PYROPTOSIS"}
ARM_EX_MANUAL = {"GSDMA", "GSDMB", "GSDMC", "GSDMD", "GSDME", "PYCARD", "CASP1",
                 "CASP4", "CASP5", "IL1B", "IL18", "NLRP1", "NLRP3", "NLRC4", "AIM2"}

arm_up_ext = set()
arm_ex_ext = set(ARM_EX_MANUAL)
for n, gs in gmt_all.items():
    if n in ARM_UP_NAMES:
        arm_up_ext |= gs
    if n in ARM_EX_NAMES:
        arm_ex_ext |= gs
note(f"   外部本体臂：上游（OXPHOS/线粒体被膜等）{len(arm_up_ext)} 基因；"
     f"执行（炎症小体/gasdermin/焦亡）{len(arm_ex_ext)} 基因")

selected = {}
seen = set()
for n, gs in gmt_all.items():
    if n in ARM_UP_NAMES or n in ARM_EX_NAMES:
        continue                                # 臂定义集本身不入签名池（防循环）
    if not any(kw in n.upper() for kw in KW):
        continue
    if not (10 <= len(gs) <= 600):
        continue
    fz = frozenset(gs)
    if fz in seen:
        continue
    seen.add(fz)
    selected[n] = gs
note(f"   入选签名 {len(selected)} 个（关键词 + 10<=n<=600 + 去重）")

# 项目 6 PCD 评分来源（MSigDB 同源：KEGG/Reactome 焦亡/铁死亡等已在 C2 内；此处补登记已发表签名口径）
# Apoptosis = HALLMARK_APOPTOSIS / KEGG_APOPTOSIS 等已含；不重复添加。

# ======================================================================
# 3. 目录指标
# ======================================================================
note("\n== 3. 目录指标")
def hedges_g_2(a, b):
    n1, n2 = len(a), len(b)
    s1, s2 = np.var(a, ddof=1), np.var(b, ddof=1)
    sp = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    g = (np.mean(b) - np.mean(a)) / sp
    g = g * (1 - 3 / (4 * (n1 + n2) - 9))
    se = np.sqrt(1 / n1 + 1 / n2 + g * g / (2 * (n1 + n2)))
    return g, se

# UCS/EIS 项目臂评分（GSE185263，供"主导臂"相关）
ucs_score185 = zc(X185[[g for g in UP_MAN if g in X185.columns]].mean(axis=1))
eis_score185 = zc(X185[[g for g in EX_MAN if g in X185.columns]].mean(axis=1))
case_pos = [X185.index.get_loc(s) for s in case185]
ctrl_pos = [X185.index.get_loc(s) for s in ctrl185]
case_pos327 = [g327.index.get_loc(s) for s in case327]
ctrl_pos327 = [g327.index.get_loc(s) for s in ctrl327]

rows = []
for name, gs in selected.items():
    meas = [g for g in gs if g in X185.columns]
    if len(meas) < 10:
        continue
    # 评分
    score = zc(X185[meas].mean(axis=1))
    g, gse = hedges_g_2(score[ctrl_pos], score[case_pos])
    # GSE32707 敏感性
    meas327 = [g for g in gs if g in g327.columns]
    g327eff = np.nan
    if len(meas327) >= 10:
        s327 = zc(g327[meas327].mean(axis=1))
        g327eff, _ = hedges_g_2(s327[ctrl_pos327], s327[case_pos327])
    # 臂占比（外部本体）
    up_ov = len(set(meas) & arm_up_ext)
    ex_ov = len(set(meas) & arm_ex_ext)
    my_ov = len(set(meas) & myeloid_genes)
    # 臂内方向一致性（GSE185263 疾病 log2FC）
    lfc = lfc185[meas]
    up_lfc = lfc[[g for g in meas if g in arm_up_ext]]
    ex_lfc = lfc[[g for g in meas if g in arm_ex_ext]]
    up_cons = float((up_lfc < 0).mean()) if len(up_lfc) else np.nan
    ex_cons = float((ex_lfc > 0).mean()) if len(ex_lfc) else np.nan
    # 主导臂（与项目 UCS/EIS 评分相关）
    r_ucs = np.corrcoef(score, ucs_score185)[0, 1]
    r_eis = np.corrcoef(score, eis_score185)[0, 1]
    rows.append(dict(signature=name, n_genes=len(gs), n_measured=len(meas),
                     upstream_share_ext=up_ov / len(meas),
                     execution_share_ext=ex_ov / len(meas),
                     myeloid_share=my_ov / len(meas),
                     up_direction_consistency=up_cons,
                     ex_direction_consistency=ex_cons,
                     g_185=g, g_185_se=gse, g_327=g327eff,
                     r_UCS=r_ucs, r_EIS=r_eis,
                     dominant_arm="UCS" if abs(r_ucs) > abs(r_eis) else "EIS"))
cat = pd.DataFrame(rows)
cat.to_csv(INTER + r"\M13_catalog.csv", index=False)
note(f"   目录 {len(cat)} 个签名；g_185 分布："
     f"min={cat['g_185'].min():+.2f} med={cat['g_185'].median():+.2f} max={cat['g_185'].max():+.2f}")
note(f"   execution_share_ext 分布：med={cat['execution_share_ext'].median():.3f} "
     f"max={cat['execution_share_ext'].max():.3f}（>0 的签名 {(cat['execution_share_ext']>0).sum()}）")
note(f"   myeloid_share 分布：med={cat['myeloid_share'].median():.3f}")

# ======================================================================
# 4. 预注册回归（主检验）
# ======================================================================
note("\n== 4. 预注册回归")
import statsmodels.api as sm

def fit_reg(df, xcols, y="g_185"):
    X = sm.add_constant(df[xcols].astype(float))
    fit = sm.OLS(df[y].astype(float), X).fit()
    return fit

fit_main = fit_reg(cat, ["execution_share_ext", "myeloid_share"])
b_exec = fit_main.params["execution_share_ext"]
p_exec_one = fit_main.pvalues["execution_share_ext"] / 2 if b_exec > 0 else 1 - fit_main.pvalues["execution_share_ext"] / 2
note(f"   主回归（g ~ exec_share + myeloid_share）：beta_exec={b_exec:+.3f} "
     f"se={fit_main.bse['execution_share_ext']:.3f} 单侧p={p_exec_one:.4g} "
     f"beta_myeloid={fit_main.params['myeloid_share']:+.3f} p={fit_main.pvalues['myeloid_share']:.3g} "
     f"R2={fit_main.rsquared:.3f} n={len(cat)}")

# LOSO
loso_signs = []
for i in range(len(cat)):
    sub = cat.drop(cat.index[i])
    f = fit_reg(sub, ["execution_share_ext", "myeloid_share"])
    loso_signs.append(np.sign(f.params["execution_share_ext"]))
loso_frac = float((np.sign(np.array(loso_signs)) == np.sign(b_exec)).mean())
note(f"   LOSO：beta_exec 方向一致 {int(np.sum(np.sign(loso_signs)==np.sign(b_exec)))}/{len(cat)} "
     f"（{loso_frac:.0%}；阈值 >=90%）")

# 敏感性：manifest 臂定义
cat["execution_share_man"] = [len(set(gs) & EX_MAN) / max(1, len([g for g in gs if g in X185.columns]))
                              for _, gs in [(r.signature, selected[r.signature]) for r in cat.itertuples()]]
cat["upstream_share_man"] = [len(set(gs) & UP_MAN) / max(1, len([g for g in gs if g in X185.columns]))
                             for _, gs in [(r.signature, selected[r.signature]) for r in cat.itertuples()]]
fit_man = fit_reg(cat, ["execution_share_man", "myeloid_share"])
note(f"   敏感性（manifest 臂定义）：beta_exec={fit_man.params['execution_share_man']:+.3f} "
     f"单侧p={fit_man.pvalues['execution_share_man']/2 if fit_man.params['execution_share_man']>0 else 1-fit_man.pvalues['execution_share_man']/2:.4g}")

# ======================================================================
# 5. 零模型（表达量匹配随机等大小基因集，B=10,000）
# ======================================================================
note("\n== 5. 零模型（B=10,000, seed=0）")
genes185 = list(X185.columns)
mean_expr = X185.mean(axis=0)
deciles = pd.qcut(mean_expr.rank(method="first"), 10, labels=False)
dec_map = dict(zip(genes185, deciles))
lfc_map = lfc185.to_dict()
z_rows = {g: Z185[g].values for g in genes185}
case_idx = [i for i, s in enumerate(X185.index) if s in case185]
ctrl_idx = [i for i, s in enumerate(X185.index) if s in ctrl185]

rng = np.random.default_rng(SEED)
null_rows = []
pool_by_dec = {d: [g for g, dd in dec_map.items() if dd == d] for d in range(10)}
sig_meas_list = [[g for g in selected[r.signature] if g in X185.columns] for r in cat.itertuples()]
gene_pos = {g: i for i, g in enumerate(genes185)}
Zmat = Z185.values                            # samples x genes（列序=genes185）
Zmat_t = Zmat.T                               # genes x samples
for it in range(B_NULL):
    sig_genes_i = sig_meas_list[it % len(cat)]   # 大小分布 = 实际签名
    ms = [pool_by_dec[dec_map[g]][int(rng.integers(len(pool_by_dec[dec_map[g]])))] for g in sig_genes_i]
    idx = [gene_pos[g] for g in ms]
    sc = Zmat_t[idx].mean(axis=0)
    g_null, _ = hedges_g_2(sc[ctrl_idx], sc[case_idx])
    ex_share = len(set(ms) & arm_ex_ext) / len(ms)
    my_share = len(set(ms) & myeloid_genes) / len(ms)
    null_rows.append(dict(g_185=g_null, execution_share_ext=ex_share, myeloid_share=my_share))
    if (it + 1) % 1000 == 0:
        note(f"   零模型 {it+1}/{B_NULL}（{time.time()-t0:.0f}s）")
null_df = pd.DataFrame(null_rows)
null_df.to_csv(INTER + r"\M13_null_model.csv", index=False)
fit_null = fit_reg(null_df, ["execution_share_ext", "myeloid_share"])
null_beta = fit_null.params["execution_share_ext"]
# 自助法经验 p（1000 次重拟合 null 数据）
boot_betas = []
rng2 = np.random.default_rng(SEED)
for b in range(1000):
    idx = rng2.integers(len(null_df), size=len(null_df))
    fb = fit_reg(null_df.iloc[idx], ["execution_share_ext", "myeloid_share"])
    boot_betas.append(fb.params["execution_share_ext"])
boot_betas = np.array(boot_betas)
emp_p_boot = float((1 + (boot_betas >= b_exec).sum()) / (len(boot_betas) + 1))
note(f"   bootstrap 经验 p（1000 次重拟合）= {emp_p_boot:.4g}；null beta 分布 "
     f"mean={boot_betas.mean():+.3f} sd={boot_betas.std():.3f}")

# ======================================================================
# 6. 反向审计（子任务3）
# ======================================================================
note("\n== 6. 反向审计（已发表策展集 vs IIAMD v1.0）")
from scipy.stats import fisher_exact
sig_full = pd.read_csv(GOV + r"\P1_IIAMD_signature_v1.0.csv")
sig_full["geneSymbol"] = sig_full["geneSymbol"].astype(str).str.upper()
sig_genes = set(sig_full["geneSymbol"])
sig_up = set(sig_full.loc[sig_full["direction"] == "up", "geneSymbol"])
sig_dn = set(sig_full.loc[sig_full["direction"] == "down", "geneSymbol"])
univ = pd.read_csv(INTER + r"\M14_GSE235046_universe_effects.csv")
univ["geneSymbol"] = univ["geneSymbol"].astype(str).str.upper()
univ = univ.set_index("geneSymbol")
univ["bm_decile"] = pd.qcut(univ["baseMean"].rank(method="first"), 10, labels=False)
universe = set(univ.index)
note(f"   IIAMD v1.0：up={len(sig_up)} down={len(sig_dn)}（小鼠符号）；宇宙 {len(universe)} 基因")

curated_sets = {
    "PMC13189447_melanoma_BAX_BAK1_BID": ["BAX", "BAK1", "BID"],
    "PMC13189447_melanoma_+mTORC2": ["BAX", "BAK1", "BID", "MTOR", "RICTOR"],
    "PMC13189447_melanoma_+mTORC2_full": ["BAX", "BAK1", "BID", "MTOR", "RICTOR", "MLST8", "MAPKAP1"],
}
audit_rows = []
for name, genes in curated_sets.items():
    genes_u = [g.upper() for g in genes]
    in_univ = [g for g in genes_u if g in universe]
    in_sig = [g for g in in_univ if g in sig_genes]
    n_overlap = len(in_sig)
    # Fisher
    n_sig_u = len(sig_genes & universe)
    n_non = len(universe) - n_sig_u
    table = [[n_overlap, len(in_univ) - n_overlap], [n_sig_u - n_overlap, n_non - (len(in_univ) - n_overlap)]]
    oddsr, p_fisher = fisher_exact(table, alternative="greater")
    # baseMean 匹配随机集零模型 B=10,000
    nulls = []
    rng3 = np.random.default_rng(SEED)
    pool = univ[~univ.index.isin(set(in_univ))]
    for it in range(B_NULL):
        rs = []
        for g in in_univ:
            dd = univ.loc[g, "bm_decile"]
            layer = pool[pool["bm_decile"] == dd]
            if len(layer):
                rs.append(layer.index[int(rng3.integers(len(layer)))])
        nulls.append(len(set(rs) & sig_genes))
    nulls = np.array(nulls)
    p_emp = float((1 + (nulls >= n_overlap).sum()) / (B_NULL + 1))
    audit_rows.append(dict(paper=name, genes=",".join(genes_u), n_genes=len(genes_u),
                           n_in_universe=len(in_univ), n_overlap_IIAMD=n_overlap,
                           expected_overlap=float(len(in_univ) * n_sig_u / len(universe)),
                           fisher_OR=oddsr, fisher_p=p_fisher,
                           null_mean=float(nulls.mean()), null_max=int(nulls.max()),
                           empirical_p=p_emp, n_perm=B_NULL))
    note(f"   {name}: overlap={n_overlap}/{len(in_univ)} (期望 {len(in_univ)*n_sig_u/len(universe):.2f}) "
         f"Fisher OR={oddsr:.3g} p={p_fisher:.3g} | 匹配零模型 mean={nulls.mean():.2f} "
         f"empirical p={p_emp:.4g}")
    # 方向登记
    for g in in_sig:
        direction = "up" if g in sig_up else ("down" if g in sig_dn else "NA")
        audit_rows.append(dict(paper=name, genes=g, n_genes=np.nan, n_in_universe=np.nan,
                               n_overlap_IIAMD=np.nan, expected_overlap=np.nan,
                               fisher_OR=np.nan, fisher_p=np.nan, null_mean=np.nan,
                               null_max=np.nan, empirical_p=np.nan, n_perm=np.nan,
                               direction=direction))
audit = pd.DataFrame(audit_rows)
audit.to_csv(INTER + r"\M13_reverse_audit.csv", index=False)

# 结直肠癌论文（PMC13376262）：35 基因进入 ML、StepCox 定终签名；补充材料仅含图（S1-S3），
# 基因表未见文本化清单 -> 留痕待补（作者侧可从论文图表/通讯作者处获得后补算）。
note("   [留痕] PMC13376262（结直肠癌）：补充材料 DataSheet1.docx 仅含图 S1-S3，无基因表；")
note("   35 基因/最终 StepCox 签名文本化清单暂缺——该文审计列为待补项，不编造数据。")

# ======================================================================
# 7. 报告
# ======================================================================
report = []
report.append("# M13 RCD 去混淆审计报告\n")
report.append("- 日期：2026-08-26；依据：M10_M13_M14_pre_registration_20260824.md §2（osf.io/ETVMJ）")
report.append(f"- 签名目录：MSigDB v7.5.1 HALLMARK/C2.CP(KEGG+REACTOME)/C5.GO 关键词筛选 "
              f"+ 10<=n<=600 + 去重 -> {len(cat)} 个签名")
report.append(f"- 臂定义（去循环化主口径）：外部本体 上游=GOCC 线粒体被膜/内膜/基质/RESPIRASOME+GOBP OXPHOS"
              f"（{len(arm_up_ext)} 基因）；执行=GOBP 焦亡/炎症小体装配+GOCC 炎症小体复合物+REACTOME 焦亡"
              f"+gasdermin/CASP1 家族（{len(arm_ex_ext)} 基因）")
report.append(f"- 髓系基因：Monaco 29 型 argmax 于髓系类型（{len(myeloid_genes)} 基因）\n")
report.append("## 1. 预注册回归（主检验）\n\n")
report.append(f"- 模型：签名疾病效应 g（GSE185263 Sepsis_COVID vs Control）~ 执行臂占比 + 髓系基因占比")
report.append(f"- **beta_exec={b_exec:+.3f}（se={fit_main.bse['execution_share_ext']:.3f}），单侧 p={p_exec_one:.4g}**")
report.append(f"- beta_myeloid={fit_main.params['myeloid_share']:+.3f}（p={fit_main.pvalues['myeloid_share']:.3g}）；"
              f"R^2={fit_main.rsquared:.3f}，n={len(cat)}")
report.append(f"- LOSO：beta_exec 方向一致 {int(np.sum(np.sign(loso_signs)==np.sign(b_exec)))}/{len(cat)}（{loso_frac:.0%}，阈值>=90%）")
report.append(f"- 零模型（B=10,000 表达量匹配随机等大小基因集）：null beta_exec={null_beta:+.3f}；"
              f"bootstrap 经验 p={emp_p_boot:.4g}")
report.append(f"- 敏感性（manifest 臂定义）：beta_exec={fit_man.params['execution_share_man']:+.3f}\n")
report.append("## 2. 反向审计\n\n")
report.append("| 论文 | 策展集 | overlap/可测 | 期望 | Fisher p | 匹配零模型 p |")
report.append("|---|---|---|---|---|---|")
for _, r in audit[audit["n_overlap_IIAMD"].notna()].iterrows():
    report.append(f"| {r['paper']} | {r['genes']} | {int(r['n_overlap_IIAMD'])}/{int(r['n_in_universe'])} "
                  f"| {r['expected_overlap']:.2f} | {r['fisher_p']:.3g} | {r['empirical_p']:.4g} |")
report.append("\n（措辞纪律：只陈述基因集一致性，不评价他人论文结论对错。）\n")
report.append("## 3. 判定门 M13\n\n")
report.append("- 判定门：预注册回归显著 + LOSO 稳定（>=90% 方向一致）+ 外部本体定义下同向。判定见正文。")
with open(GOV + r"\M13_Confounding_Audit_Report.md", "w", encoding="utf-8") as f:
    f.write("\n".join(report))
with open(ROOT + r"\03_LOGS\M13_confounding_audit_log.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(log))
note(f"\nDONE in {time.time()-t0:.1f}s -> {GOV}\\M13_Confounding_Audit_Report.md")
