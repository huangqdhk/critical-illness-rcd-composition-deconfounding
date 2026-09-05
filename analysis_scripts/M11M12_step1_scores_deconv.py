# -*- coding: utf-8 -*-
"""
M11M12_step1_scores_deconv.py — M11/M12 步骤 1：五队列评分 + 组成反卷积 + 组成残差 MDI
================================================================================
前瞻注册：M11_M12_pre_registration_20260827.md（§1.3/§2.4 平台口径锁定）
口径（冻结）：
- 评分：mdi_v1.0（臂均值 -> 队列内全样本 z；mdi_lib.arm_scores）
- GSE215865：官方 logCPM（log2 空间）；ENSG 去版本号 -> HUGO（本地缓存 + mygene 补缺）；
  反卷积主口径 ABIS RNAseq 17 型 NNLS（Monaco 29 型 NNLS 交叉）；重复样（受试者×日标签）取均值后评分。
- GSE54514：non-normalized -> quantile + log2(x+1)（lumi lumiN 等价，预注册口径）；
  ILMN->HUGO 用 GEO GPL6947 官方平台表；max-mean <=3 探针塌陷；ABIS Micro 11 型 NNLS（跨平台敏感性口径）。
- GSE148871（血）：GPL570 系列矩阵 -> 探针塌陷；ABIS Micro NNLS；与 M2 冻结 mdi_v1.0 评分交叉核对。
- GSE212865（2 波敏感性）：GPL23159 系列矩阵；ABIS Micro NNLS。
- GSE106878：non-normalized nuID -> quantile + log2(x+1)；nuID->HUGO 用 GEO GPL10295 官方平台表；
  ABIS Micro NNLS（跨平台敏感性口径）。
- 组成残差：逐队列 OLS MDI ~ MDI_comp -> resid（沿用 M10A 口径）。
输出：
  _intermediate/M11M12_step1_per_sample.csv（逐样本长表）
  _intermediate/M11M12_step1_proportions.csv（ABIS 主口径细胞比例）
  _intermediate/M11M12_step1_reference_means.csv（ABIS 参考谱 log2 均值，缓存）
  03_LOGS/M11M12_step1_log.txt
"""
import gzip, os, re, sys, time, warnings
from collections import Counter, defaultdict

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.linalg import cholesky, solve_triangular

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
RAW = ROOT + r"\00_RAW_DATA"
INTER = ROOT + r"\_intermediate"
GOV = ROOT + r"\04_AUDIT_GOVERNANCE"
LOG = ROOT + r"\03_LOGS"
sys.path.insert(0, ROOT)
import mdi_lib as L

t0 = time.time()
log = []
def note(m=""):
    log.append(m); print(m, flush=True)

man_df, arms, alias = L.load_manifest()
UP, EX = arms["upstream_collapse"], arms["execution_induction"]
ALL80 = set(UP) | set(EX)

# ======================================================================
# 0. ABIS 参考谱（RNAseq 17 型 / Micro 11 型），log2(+1)
# ======================================================================
note("== 0. ABIS 参考谱")
REF_CACHE = INTER + r"\M11M12_step1_reference_means.csv"
if os.path.exists(REF_CACHE):
    ref_cache = pd.read_csv(REF_CACHE)
    note("   [缓存复用] M11M12_step1_reference_means.csv")
else:
    rows = []
    for fname, tag in [(r"ABIS\sigmatrixRNAseq.txt", "ABIS_rna"),
                       (r"ABIS\sigmatrixMicro.txt", "ABIS_micro")]:
        df = pd.read_csv(RAW + "\\" + fname, sep="\t", index_col=0)
        df.index = [str(i).strip().strip('"\'').upper() for i in df.index]
        df.columns = [str(c).strip().strip('"\'') for c in df.columns]
        df = df.groupby(level=0).mean()
        lg = np.log2(df + 1.0)
        rows.append(pd.DataFrame(lg.stack()).reset_index().rename(
            columns={"level_0": "gene", "level_1": "celltype", 0: "log2val"}).assign(reference=tag))
    ref_cache = pd.concat(rows)
    ref_cache.to_csv(REF_CACHE, index=False)
    note(f"   ABIS 参考谱缓存已建：{len(ref_cache)} 行")

def ref_matrix(tag):
    """返回 (celltypes x genes) DataFrame。"""
    sub = ref_cache[ref_cache["reference"] == tag]
    return sub.pivot(index="celltype", columns="gene", values="log2val")

ABIS_R = ref_matrix("ABIS_rna")     # 17 型
ABIS_M = ref_matrix("ABIS_micro")   # 11 型
# Monaco 29 型（M10A 缓存，log2(TPM+1)，细胞型 x 基因）——组成预测 MDI 主口径
MON = pd.read_csv(INTER + r"\M10A_monaco_ct_means.csv", index_col=0)
MON.index = [str(i).strip().upper() for i in MON.index]
MON.columns = [str(c).strip().upper() for c in MON.columns]
note(f"   ABIS RNAseq: {ABIS_R.shape[0]} 型 x {ABIS_R.shape[1]} 基因")
note(f"   ABIS Micro:  {ABIS_M.shape[0]} 型 x {ABIS_M.shape[1]} 基因 | 类型: {list(ABIS_M.index)}")
note(f"   Monaco 29 型: {MON.shape[0]} x {MON.shape[1]} 基因")
note(f"   臂覆盖：ABIS_R up/ex = {len([g for g in UP if g in ABIS_R.columns])}/{len([g for g in EX if g in ABIS_R.columns])}"
     f" | ABIS_M = {len([g for g in UP if g in ABIS_M.columns])}/{len([g for g in EX if g in ABIS_M.columns])}"
     f" | Monaco = {len([g for g in UP if g in MON.columns])}/{len([g for g in EX if g in MON.columns])}")

# ======================================================================
# 1. 反卷积（NNLS 主口径）与合成评分（mdi_v1.0 同口径）
# ======================================================================
def zc(s):
    s = np.asarray(s, float); sd = s.std(ddof=1)
    return (s - s.mean()) / sd if sd > 0 else np.zeros_like(s)

def deconvolve_nnls(X, ref_log):
    """X: samples x genes(log2); ref_log: celltypes x genes(log2)。返回 proportions。"""
    common = [g for g in X.columns if g in ref_log.columns]
    A = ref_log[common].values.T
    Xc = X[common].values
    AtA = A.T @ A
    n_ct = A.shape[1]
    Lc = cholesky(AtA + np.eye(n_ct) * 1e-10, lower=True)
    rows = []
    for i in range(Xc.shape[0]):
        rhs = solve_triangular(Lc, A.T @ Xc[i], lower=True)
        pi, _ = nnls(Lc.T, rhs)
        s = pi.sum()
        rows.append((pi / s) if s > 0 else np.full(n_ct, 1.0 / n_ct))
    return pd.DataFrame(rows, index=X.index, columns=ref_log.index), common

def synthetic_scores(prop, ref_log, up_genes, ex_genes):
    up_in = [g for g in up_genes if g in ref_log.columns]
    ex_in = [g for g in ex_genes if g in ref_log.columns]
    mu_up = ref_log[up_in].values
    mu_ex = ref_log[ex_in].values
    P = prop.values
    uc = P @ mu_up.mean(axis=1)
    ec = P @ mu_ex.mean(axis=1)
    out = pd.DataFrame(index=prop.index)
    out["UCS_comp"] = zc(uc); out["EIS_comp"] = zc(ec)
    out["MDI_comp"] = out["EIS_comp"] - out["UCS_comp"]
    out["n_up_ref"] = len(up_in); out["n_ex_ref"] = len(ex_in)
    return out

def neutrophil_monocyte(prop):
    """返回 DataFrame: Neut / Mono / comp（Neut+Mono 合计比例）。"""
    out = pd.DataFrame(index=prop.index)
    neut_cols = [c for c in prop.columns if "neutrophil" in c.lower()]
    mono_cols = [c for c in prop.columns if c.lower() in ("monocytes", "monocytes c", "monocytes nc+i")]
    out["Neut"] = prop[neut_cols].sum(axis=1) if neut_cols else np.nan
    out["Mono"] = prop[mono_cols].sum(axis=1) if mono_cols else np.nan
    out["comp"] = out["Neut"] + out["Mono"]
    return out

def score_and_resid(mat_gene, ref_main, tag):
    """mat_gene: samples x genes(log2, 大写符号)。
    - 比例：ref_main（ABIS，预注册锁定口径）NNLS + Monaco NNLS 交叉；
    - 组成预测 MDI：Monaco（全基因组臂覆盖）主口径；
    - 残差：逐队列 OLS MDI ~ MDI_comp。"""
    sc = L.arm_scores(mat_gene, arms)
    prop_main, common = deconvolve_nnls(mat_gene, ref_main)
    prop_mon, common_mon = deconvolve_nnls(mat_gene, MON)
    nm = neutrophil_monocyte(prop_main)
    nm_mon = neutrophil_monocyte(prop_mon).rename(columns=lambda c: c + "_monaco")
    sy = synthetic_scores(prop_mon, MON, UP, EX)
    sy_abis = synthetic_scores(prop_main, ref_main, UP, EX)
    sy_abis = sy_abis.rename(columns=lambda c: c + "_abis")
    out = pd.concat([sc, nm, nm_mon, sy, sy_abis], axis=1)
    A = np.column_stack([np.ones(len(out)), out["MDI_comp"].values])
    beta, *_ = np.linalg.lstsq(A, out["MDI"].values, rcond=None)
    out["MDI_resid"] = out["MDI"].values - A @ beta
    out["b0"] = beta[0]; out["b1"] = beta[1]
    out["comp_ref"] = tag
    out["n_common_deconv"] = len(common)
    out["n_common_monaco"] = len(common_mon)
    return out, prop_main, prop_mon

def collapse_maxmean(probe_mat, probe2sym, max_probes=3):
    """probe_mat: samples x probes；probe2sym: dict probe->symbol(大写)。全探针塌陷。"""
    g2p = defaultdict(list)
    for p, s in probe2sym.items():
        if p in probe_mat.columns:
            g2p[s].append(p)
    out = pd.DataFrame(index=probe_mat.index)
    for g, ps in g2p.items():
        sub = probe_mat[ps]
        if len(ps) > max_probes:
            best = sub.mean(axis=0).nlargest(max_probes).index
            sub = sub[best]
        out[g] = sub.mean(axis=1)
    return out

def quantile_normalize(mat):
    """列分位数归一化（等价 lumiN）。mat: samples x features。"""
    M = mat.values.astype(float).T   # features x samples
    order = np.argsort(M, axis=1)
    sortedM = np.take_along_axis(M, order, axis=1)
    means = np.nanmean(sortedM, axis=0)
    Q = np.broadcast_to(means, M.shape)
    out = np.empty_like(M)
    for i in range(M.shape[0]):
        out[i] = Q[i][np.argsort(order[i])]
    return pd.DataFrame(out.T, index=mat.index, columns=mat.columns)

# ======================================================================
# 2. GSM 元数据解析
# ======================================================================
def parse_gsm_metadata(path):
    txt = open(path, encoding="utf-8", errors="replace").read()
    recs = txt.split("^SAMPLE")
    out = []
    for r in recs[1:]:
        d = {}
        for ln in r.split("\n"):
            ln = ln.strip()
            if ln.startswith("!Sample_geo_accession"):
                d["gsm"] = ln.split("=", 1)[1].strip()
            elif ln.startswith("!Sample_title"):
                d["title"] = ln.split("=", 1)[1].strip()
            elif ln.startswith("!Sample_description"):
                d.setdefault("descriptions", []).append(ln.split("=", 1)[1].strip())
            elif ln.startswith("!Sample_characteristics_ch1"):
                kv = ln.split("=", 1)[1].strip()
                if ":" in kv:
                    k, v = kv.split(":", 1)
                    d[k.strip().lower()] = v.strip()
        out.append(d)
    return pd.DataFrame(out)

# ======================================================================
# 3. GSE215865（M11 主队列 1）
# ======================================================================
note("\n== 3. GSE215865")
cols_file = RAW + r"\GSE215865_COVID19_WholeBlood_Longitudinal\GSE215865_rnaseq_logCPM_matrix.csv.gz"
with gzip.open(cols_file, "rt", encoding="utf-8", errors="replace") as f:
    header = f.readline().rstrip("\n")
colnames = header.split(",")
note(f"   logCPM 列数 {len(colnames)}（含 ENSG 首列）")
# ENSG 映射缓存
cache_union = {}
for c in (INTER + r"\M10A_monaco_ensg2sym.csv", INTER + r"\P2_gse185263_en2sym.csv"):
    for k, v in pd.read_csv(c).values:
        cache_union[str(k).strip().upper().split(".")[0]] = str(v).strip().upper()
note(f"   本地 ENSG->sym 缓存 {len(cache_union)} 条")
need_syms = set(UP) | set(EX) | set(ABIS_R.columns) | set(ABIS_M.columns)
def ensg2sym(ensg):
    e = str(ensg).split(".")[0].upper()
    return cache_union.get(e, "")
# 逐行流式读矩阵（只读一次）：
df215 = pd.read_csv(cols_file, compression="gzip", index_col=0, low_memory=False,
                    na_values=["NA"], dtype={c: np.float32 for c in colnames[1:]})
df215 = df215.astype(np.float32)
note(f"   矩阵 {df215.shape[0]} 基因 x {df215.shape[1]} 样本")
ensg_list = [str(i) for i in df215.index]
sym_list = [ensg2sym(e) for e in ensg_list]
missing = sorted({ensg_list[i] for i in range(len(ensg_list)) if (not sym_list[i]) and ensg_list[i] in need_syms})
if missing:
    note(f"   缓存缺失需 mygene 查询的基因 {len(missing)}: {missing}")
    import mygene
    mg = mygene.MyGeneInfo()
    res = mg.querymany(missing, scopes="ensembl.gene", fields="symbol", species="human",
                       returnall=True, verbose=False)
    for r in res.get("out", []):
        if r.get("symbol"):
            cache_union[str(r["query"]).split(".")[0].upper()] = str(r["symbol"]).upper()
    sym_list = [ensg2sym(e) for e in ensg_list]
    pd.DataFrame({"ensg": list(cache_union.keys()), "symbol": list(cache_union.values())}).to_csv(
        INTER + r"\M11M12_step1_ensg2sym_union.csv", index=False)
n_sym = sum(1 for s in sym_list if s)
note(f"   符号化覆盖 {n_sym}/{len(ensg_list)}")
df215["sym"] = sym_list
df215 = df215[df215["sym"] != ""]
df215 = df215.groupby("sym").mean(numeric_only=True)
df215 = df215.fillna(0.0)   # NA logCPM = 零计数基因 -> log2(0+1)=0
note(f"   符号化后 {df215.shape[0]} 基因")
# 列 -> subject/label/plate；重复样（受试者×日标签）取均值
pat = re.compile(r"^Subj_([0-9a-fA-F]+)T(\d+[A-Za-z]?)_Plate_(\d+)$")
meta215 = []
for c in df215.columns:
    m = pat.match(str(c))
    if m:
        meta215.append((str(c), m.group(1), "T" + re.sub(r"[A-Za-z]+$", "", m.group(2)), m.group(2)))
    else:
        meta215.append((str(c), "", "", ""))
mdf = pd.DataFrame(meta215, columns=["col", "subject", "day_label", "raw_label"])
note(f"   可解析列 {int((mdf['subject'] != '').sum())}/{len(mdf)}")
df215_t = df215.T
df215_t.index.name = "col"
m215 = mdf.set_index("col")
# 受试者×日标签重复 -> 均值
m215["key"] = m215["subject"] + "|" + m215["day_label"]
dup_cols = m215[m215.duplicated(subset=["key"], keep=False) & (m215["subject"] != "")]
n_dup = dup_cols["key"].nunique()
note(f"   受试者×日标签重复样本 {n_dup} 组 -> 取均值")
merged = df215_t.groupby(m215["key"]).mean()
merged.index = merged.index.str.split("|", expand=True)
merged.index.names = ["subject", "day_label"]
day_order = sorted(set(merged.index.get_level_values(1)), key=lambda x: int(re.sub(r"\D", "", x)))
note(f"   day 级标签 {len(day_order)}: {day_order}")
per_subj = merged.groupby(level=0).size()
note(f"   受试者 {len(per_subj)}；>=2 天 {int((per_subj >= 2).sum())}、>=3 天 {int((per_subj >= 3).sum())}")
g215 = merged.loc[per_subj[per_subj >= 1].index]  # 全受试者
# 评分（全 day 级样本）+ ABIS RNAseq NNLS
s215, p215, pm215 = score_and_resid(g215, ABIS_R, "ABIS_rna_NNLS")
s215["cohort"] = "GSE215865"
s215.index = g215.index
note(f"   GSE215865 评分: n={len(s215)} | n_up={s215['n_up'].iloc[0]} n_ex={s215['n_ex'].iloc[0]} | "
     f"b1(MDI~MDI_comp)={s215['b1'].iloc[0]:+.3f} | n_common={s215['n_common_deconv'].iloc[0]}")

# ======================================================================
# 4. GSE54514（M11 主队列 2）
# ======================================================================
note("\n== 4. GSE54514")
raw54 = pd.read_csv(RAW + r"\GSE54514_Sepsis_PAXgene_WholeBlood\GSE54514_non-normalized.txt.gz",
                    sep="\t", index_col=0, low_memory=False,
                    dtype={c: np.float32 for c in pd.read_csv(
                        RAW + r"\GSE54514_Sepsis_PAXgene_WholeBlood\GSE54514_non-normalized.txt.gz",
                        sep="\t", nrows=0).columns[1:]})
expr_cols = [c for c in raw54.columns if not str(c).endswith("Detection Pval")]
g545_raw = raw54[expr_cols].T
note(f"   表达矩阵 {g545_raw.shape[0]} 样本 x {g545_raw.shape[1]} 探针")
# 样本注释
gsm54 = parse_gsm_metadata(RAW + r"\GSE54514_Sepsis_PAXgene_WholeBlood\GSE54514_gsm_metadata.txt")
note(f"   GSM 记录 {len(gsm54)}")
# title: "PAXgene whole blood, {grp}, Day_{n}, ID={id}"
def title_parse(t):
    m = re.search(r",\s*(.+?),\s*Day_(\d+),\s*ID=(\S+)", str(t))
    if m:
        return m.group(1).strip(), "Day_" + m.group(2), m.group(3)
    return None, None, None
g54p = gsm54["title"].apply(title_parse)
gsm54["grp_title"] = g54p.str[0]
gsm54["day_title"] = g54p.str[1]
gsm54["id_title"] = g54p.str[2]
gsm54["sentrix"] = gsm54["descriptions"].apply(lambda ds: ds[0] if ds else "")
note(f"   title 解析: ID {gsm54['id_title'].nunique()} | group {Counter(gsm54['grp_title'])}")
# characteristics
note(f"   characteristics 字段: {sorted(set(k for d in gsm54.to_dict('records') for k in d if k not in ('gsm','title','descriptions')))}")
# 列 sentrix -> GSM
sentrix2gsm = dict(zip(gsm54["sentrix"], gsm54["gsm"]))
g545_raw["gsm"] = g545_raw.index.map(lambda c: sentrix2gsm.get(str(c), ""))
n_matched = int((g545_raw["gsm"] != "").sum())
note(f"   矩阵列->GSM 匹配 {n_matched}/{g545_raw.shape[0]}")
g545_raw = g545_raw[g545_raw["gsm"] != ""].set_index("gsm", drop=True)
# quantile + log2
q545 = quantile_normalize(g545_raw)
g545_log = np.log2(q545 + 1.0)
note("   quantile + log2 完成")
# GPL6947 表 -> probe2sym
gpl6947 = pd.read_csv(RAW + r"\GPL6947_table.txt", sep="\t", skiprows=lambda i: i < 0, dtype=str, on_bad_lines="skip", comment="!", low_memory=False)
# 重新读：平台表在 !platform_table_begin 之后
p2s = {}
with open(RAW + r"\GPL6947_table.txt", encoding="utf-8", errors="replace") as f:
    lines = f.read().split("\n")
hdr_idx = next(i for i, ln in enumerate(lines) if ln.startswith("!platform_table_begin"))
header = lines[hdr_idx + 1].split("\t")
sym_col = next(i for i, h in enumerate(header) if h.strip().strip('"').lower() == "symbol")
id_col = next(i for i, h in enumerate(header) if h.strip().strip('"').lower() in ("id",))
n_mapped = 0
for ln in lines[hdr_idx + 2:]:
    if ln.startswith("!platform_table_end") or not ln.strip():
        break
    parts = ln.split("\t")
    if len(parts) > max(sym_col, id_col):
        s = parts[sym_col].strip().strip('"').upper()
        pid = parts[id_col].strip().strip('"')
        if s and s != "---":
            p2s[pid] = s
            n_mapped += 1
note(f"   GPL6947 探针映射 {n_mapped} 条")
g545_gene = collapse_maxmean(g545_log, p2s)
note(f"   基因级矩阵 {g545_gene.shape[0]} 样本 x {g545_gene.shape[1]} 基因")
s545, p545, pm545 = score_and_resid(g545_gene, ABIS_M, "ABIS_micro_NNLS")
s545["cohort"] = "GSE54514"
s545.index = g545_gene.index
# 注释并入
anno545 = gsm54.set_index("gsm")
for k in ("grp_title", "day_title", "id_title", "gender", "age (years)", "neutrophil proportion",
          "disease status", "group_day", "group_id"):
    if k in anno545.columns:
        s545[k] = anno545.reindex(s545.index)[k].values
note(f"   GSE54514 评分: n={len(s545)} | n_up={s545['n_up'].iloc[0]} n_ex={s545['n_ex'].iloc[0]} | "
     f"b1={s545['b1'].iloc[0]:+.3f}")

# ======================================================================
# 5. GSE148871（M11 队列 3 + M12 复现层；仅血样）
# ======================================================================
note("\n== 5. GSE148871")
meta4, s4, m4_ = L.parse_series_matrix(RAW + r"\GSE148871_COPD_AE\GSE148871_series_matrix.txt.gz")
note(f"   系列矩阵 {m4_.shape[0]} 样本 x {m4_.shape[1]} 探针")
# characteristics 块结构
note("   characteristics 块:")
for b in meta4.get("!Sample_characteristics_ch1", [])[:12]:
    if b:
        note(f"     {b[0]} | 示例: {b[1] if len(b)>1 else ''}")
def char_vec(meta, label):
    for block in meta.get("!Sample_characteristics_ch1", []):
        if block and label.lower() in block[0].lower():
            return [b.split(":", 1)[1].strip() if ":" in b else b.strip() for b in block]
    return None
tissue4 = char_vec(meta4, "tissue")
arm4 = char_vec(meta4, "treatment") or char_vec(meta4, "arm")
visit4 = char_vec(meta4, "visit")
subj4 = char_vec(meta4, "subject") or char_vec(meta4, "patient") or char_vec(meta4, "donor")
sex4 = char_vec(meta4, "gender") or char_vec(meta4, "sex")
note(f"   tissue: {Counter(tissue4) if tissue4 else 'N/A'}")
note(f"   treatment: {Counter(arm4) if arm4 else 'N/A'}")
note(f"   visit: {Counter(visit4) if visit4 else 'N/A'}")
note(f"   subject n: {len(set(subj4)) if subj4 else 'N/A'} | sex: {Counter(sex4) if sex4 else 'N/A'}")
# 探针塌陷（GPL570 全探针）
gpl570 = pd.read_csv(RAW + r"\GPL570_probe2symbol_hgu133plus2db.csv")
pc = [c for c in gpl570.columns if c.lower() in ("probe_id", "probeid", "id")][0]
sc = [c for c in gpl570.columns if c.lower() in ("symbol", "gene_symbol")][0]
p2s570 = dict(zip(gpl570[pc].astype(str).str.strip(), gpl570[sc].astype(str).str.strip().str.upper()))
g148_gene = collapse_maxmean(m4_, p2s570)
note(f"   基因级 {g148_gene.shape[0]} 样本 x {g148_gene.shape[1]} 基因")
s148, p148, pm148 = score_and_resid(g148_gene, ABIS_M, "ABIS_micro_NNLS")
s148["cohort"] = "GSE148871"
s148.index = g148_gene.index
s148["tissue"] = tissue4; s148["treatment"] = arm4; s148["visit"] = visit4
s148["subject"] = subj4; s148["sex"] = sex4
s148 = s148[s148["tissue"].astype(str).str.lower().str.contains("blood")]
note(f"   GSE148871 血样 n={len(s148)} | n_up={s148['n_up'].iloc[0]} n_ex={s148['n_ex'].iloc[0]}")
# M2 锚点核对
m2 = pd.read_csv(INTER + r"\M2_per_sample_scores.csv", low_memory=False)
m2["sample"] = m2["sample"].astype(str)
m2b = m2[m2["cohort"] == "GSE148871"].set_index("sample")
common148 = sorted(set(s148.index) & set(m2b.index))
r_mdi = np.corrcoef(s148.loc[common148, "MDI"].values, m2b.loc[common148, "MDI"].values)[0, 1]
r_ucs = np.corrcoef(s148.loc[common148, "UCS"].values, m2b.loc[common148, "UCS"].values)[0, 1]
note(f"   M2 锚点核对（n={len(common148)}）: MDI r={r_mdi:+.4f} UCS r={r_ucs:+.4f}")

# ======================================================================
# 6. GSE212865（2 波敏感性）
# ======================================================================
note("\n== 6. GSE212865")
meta2, s2, m2_ = L.parse_series_matrix(RAW + r"\GSE212865_data\GSE212865_series_matrix.txt.gz")
note(f"   系列矩阵 {m2_.shape[0]} 样本 x {m2_.shape[1]} 探针")
gpl23159 = pd.read_csv(RAW + r"\GPL23159_probe2symbol_clariomsdb.csv")
pc = [c for c in gpl23159.columns if c.lower() in ("probe_id", "probeid", "id")][0]
sc = [c for c in gpl23159.columns if c.lower() in ("symbol", "gene_symbol")][0]
p2s23159 = dict(zip(gpl23159[pc].astype(str).str.strip(), gpl23159[sc].astype(str).str.strip().str.upper()))
g212_gene = collapse_maxmean(m2_, p2s23159)
note(f"   基因级 {g212_gene.shape[0]} 样本 x {g212_gene.shape[1]} 基因")
s212, p212, pm212 = score_and_resid(g212_gene, ABIS_M, "ABIS_micro_NNLS")
s212["cohort"] = "GSE212865"
s212.index = g212_gene.index
# manifest 注释（患者/时点/亚组）
man = pd.read_csv(GOV + r"\SAMPLE_MANIFEST_v1.0.csv", low_memory=False)
m212 = man[man["dataset_id"] == "GSE212865"].set_index("gsm")
s212["patient"] = m212.reindex(s212.index)["subcohort"].values
s212["timepoint"] = m212.reindex(s212.index)["timepoint"].values
# 组别：从 M2 取 group（Covid19_SDRA / Control）
m2s212 = m2[m2["cohort"] == "GSE212865"].set_index("sample")
s212["group"] = m2s212.reindex(s212.index)["group"].values
note(f"   GSE212865 n={len(s212)} | 组: {Counter(s212['group'].dropna())} | 时点: {Counter(s212['timepoint'].dropna())}")

# ======================================================================
# 7. GSE106878（M12 确认队列）
# ======================================================================
note("\n== 7. GSE106878")
raw878 = pd.read_csv(RAW + r"\GSE106878_Sepsis_Hydrocortisone_CORTICUS\GSE106878_non-normalized_data.txt.gz",
                     sep="\t", index_col=0, low_memory=False,
                     dtype={c: np.float32 for c in pd.read_csv(
                         RAW + r"\GSE106878_Sepsis_Hydrocortisone_CORTICUS\GSE106878_non-normalized_data.txt.gz",
                         sep="\t", nrows=0).columns[1:]})
expr_cols = [c for c in raw878.columns if not str(c).endswith("Detection Pval")]
g878_raw = raw878[expr_cols].T
note(f"   表达矩阵 {g878_raw.shape[0]} 样本 x {g878_raw.shape[1]} nuID")
gsm878 = parse_gsm_metadata(RAW + r"\GSE106878_Sepsis_Hydrocortisone_CORTICUS\GSE106878_gsm_metadata.txt")
note(f"   GSM 记录 {len(gsm878)}")
gsm878["sentrix"] = gsm878["descriptions"].apply(lambda ds: ds[1] if len(ds) > 1 else (ds[0] if ds else ""))
sentrix2gsm878 = dict(zip(gsm878["sentrix"], gsm878["gsm"]))
g878_raw["gsm"] = g878_raw.index.map(lambda c: sentrix2gsm878.get(str(c), ""))
note(f"   矩阵列->GSM 匹配 {int((g878_raw['gsm'] != '').sum())}/{g878_raw.shape[0]}")
g878_raw = g878_raw[g878_raw["gsm"] != ""].set_index("gsm", drop=True)
q878 = quantile_normalize(g878_raw)
g878_log = np.log2(q878 + 1.0)
note("   quantile + log2 完成")
# GPL10295 表 -> nuID->symbol
p2s878 = {}
with open(RAW + r"\GPL10295_table.txt", encoding="utf-8", errors="replace") as f:
    lines = f.read().split("\n")
hdr_idx = next(i for i, ln in enumerate(lines) if ln.startswith("!platform_table_begin"))
header = lines[hdr_idx + 1].split("\t")
sym_col = next(i for i, h in enumerate(header) if h.strip().strip('"').lower() == "symbol")
id_col = next(i for i, h in enumerate(header) if h.strip().strip('"').lower() in ("id",))
for ln in lines[hdr_idx + 2:]:
    if ln.startswith("!platform_table_end") or not ln.strip():
        break
    parts = ln.split("\t")
    if len(parts) > max(sym_col, id_col):
        s = parts[sym_col].strip().strip('"').upper()
        pid = parts[id_col].strip().strip('"')
        if s and s != "---":
            p2s878[pid] = s
note(f"   GPL10295 nuID 映射 {len(p2s878)} 条")
g878_gene = collapse_maxmean(g878_log, p2s878)
note(f"   基因级 {g878_gene.shape[0]} 样本 x {g878_gene.shape[1]} 基因")
s878, p878, pm878 = score_and_resid(g878_gene, ABIS_M, "ABIS_micro_NNLS")
s878["cohort"] = "GSE106878"
s878.index = g878_gene.index
anno878 = gsm878.set_index("gsm")
for k in ("individual", "timepoint", "treatment", "gender", "age", "acth", "survival (28 days)"):
    if k in anno878.columns:
        s878[k] = anno878.reindex(s878.index)[k].values
note(f"   GSE106878 n={len(s878)} | 臂: {Counter(s878['treatment'])} | 时点: {Counter(s878['timepoint'])}")
note(f"   n_up={s878['n_up'].iloc[0]} n_ex={s878['n_ex'].iloc[0]} | b1={s878['b1'].iloc[0]:+.3f}")

# ======================================================================
# 8. 汇总导出
# ======================================================================
note("\n== 8. 汇总导出")
pieces = []
for s in (s215, s545, s148, s212, s878):
    d = s.copy()
    if isinstance(s.index, pd.MultiIndex):
        d["sample_id"] = [f"{i[0]}|{i[1]}" for i in s.index]
        d["subject"] = [i[0] for i in s.index]
        d["day_label"] = [i[1] for i in s.index]
    else:
        d["sample_id"] = [str(i) for i in s.index]
    d = d.reset_index(drop=True)
    pieces.append(d)
all_s = pd.concat(pieces, ignore_index=True)
# 统一标识列：subject / sex 归一化
all_s["subject_norm"] = all_s["subject"].fillna(
    all_s.get("id_title")).fillna(all_s.get("patient")).fillna(all_s.get("individual"))
all_s["sex_norm"] = all_s.get("sex").fillna(all_s.get("gender"))
# 列顺序整理：标识列前置
ident_cols = ["cohort", "sample_id", "subject_norm", "day_label"]
other = [c for c in all_s.columns if c not in ident_cols]
all_s = all_s[ident_cols + other]
all_s.to_csv(INTER + r"\M11M12_step1_per_sample.csv", index=False)
note(f"   逐样本长表 {len(all_s)} 行 -> M11M12_step1_per_sample.csv")
prop_rows = []
for p, cohort, refname in [(p215, "GSE215865", "ABIS_rna"), (pm215, "GSE215865", "Monaco"),
                            (p545, "GSE54514", "ABIS_micro"), (pm545, "GSE54514", "Monaco"),
                            (p148, "GSE148871", "ABIS_micro"), (pm148, "GSE148871", "Monaco"),
                            (p212, "GSE212865", "ABIS_micro"), (pm212, "GSE212865", "Monaco"),
                            (p878, "GSE106878", "ABIS_micro"), (pm878, "GSE106878", "Monaco")]:
    pr = p.copy()
    pr = pr.reset_index().rename(columns={"index": "sample"})
    pr["cohort"] = cohort
    pr["reference"] = refname
    prop_rows.append(pr)
props = pd.concat(prop_rows, ignore_index=True)
props.to_csv(INTER + r"\M11M12_step1_proportions.csv", index=False)
note(f"   比例长表 {len(props)} 行 -> M11M12_step1_proportions.csv")

with open(LOG + r"\M11M12_step1_log.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(log))
note(f"\nDONE in {time.time()-t0:.1f}s -> {LOG}\\M11M12_step1_log.txt")
