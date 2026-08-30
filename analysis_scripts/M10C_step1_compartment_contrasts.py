# -*- coding: utf-8 -*-
"""
M10C_step1_compartment_contrasts.py — M10 检验 C：区室对照（描述性，不设硬判定）
===============================================================================
前瞻注册：M10_M13_M14_pre_registration_20260824.md §1.4
- C1 痰 vs 血（GSE148871，COPD AE，患者配对）：UCS/EIS/MDI 配对差异 + Wilcoxon/配对t
- C2 血 vs 气道（GSE180578，COVID PBMC vs 气管吸出物 TA，RAW 全基因组矩阵）：
  同患者同时点配对臂评分（全部细胞合并口径，与痰/血口径一致）
- C3 BALF vs PBMC 跨队列（GSE145926 vs GSE158055，severe 供者，P0 pseudobulk arm 评分）
- C4 肺 Visium（GSE271370）：引用已有 Table S55b（MDI 组间 p=0.715），无新增计算
输出：_intermediate/M10C_sputum_blood.csv / M10C_pbmc_ta.csv / M10C_balf_pbmc.csv
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
MANIFEST = ROOT + r"\04_AUDIT_GOVERNANCE\Mitoxyperilysis_Gene_Manifest_v1.0.csv"

sys.path.insert(0, ROOT)
import mdi_lib as L

t0 = time.time()
log = []
def note(m=""):
    log.append(m); print(m, flush=True)

def zc(s):
    s = np.asarray(s, float); sd = s.std(ddof=1)
    return (s - s.mean()) / sd if sd > 0 else s * 0

man = pd.read_csv(MANIFEST)
UP = set(man.loc[man["arm"] == "upstream_collapse", "hgnc_symbol"].astype(str).str.upper())
EX = set(man.loc[man["arm"] == "execution_induction", "hgnc_symbol"].astype(str).str.upper())

# ======================================================================
# C1 痰 vs 血（GSE148871）
# ======================================================================
note("== C1 痰 vs 血（GSE148871，患者配对）")
m2 = pd.read_csv(INTER + r"\M2_per_sample_scores.csv", low_memory=False)
g = m2[m2["cohort"] == "GSE148871"].copy()
g["sample"] = g["sample"].astype(str)
g["subject"] = g["subject"].astype(str)
pairs = []
for subj, s in g.groupby("subject"):
    b = s[s["tissue"].astype(str).str.lower().str.contains("blood")]
    sp = s[s["tissue"].astype(str).str.lower().str.contains("sputum")]
    if len(b) == 1 and len(sp) == 1:
        pairs.append((subj, b.iloc[0], sp.iloc[0]))
note(f"   配对患者 {len(pairs)} 例")
from scipy.stats import wilcoxon, ttest_rel, mannwhitneyu
# 非配对全样本口径（与 v4 正文一致的主口径）
blood_all = g[g["tissue"].astype(str).str.lower().str.contains("blood")]
sp_all = g[g["tissue"].astype(str).str.lower().str.contains("sputum")]
note(f"   全样本：血 {len(blood_all)} / 痰 {len(sp_all)}；配对 {len(pairs)} 例")
rows = []
for fam in ("UCS", "EIS", "MDI"):
    d_blood = [b[fam] for _, b, _ in pairs]
    d_sp = [s[fam] for _, _, s in pairs]
    diffs = [s[fam] - b[fam] for _, b, s in pairs]
    w, pw = wilcoxon(diffs) if len(diffs) >= 5 and any(x != 0 for x in diffs) else (np.nan, np.nan)
    t, pt = ttest_rel(d_sp, d_blood) if len(diffs) >= 5 else (np.nan, np.nan)
    u, pu = mannwhitneyu(sp_all[fam].dropna(), blood_all[fam].dropna(), alternative="two-sided")
    rows.append(dict(family=fam, n_pairs=len(pairs),
                     n_blood_all=len(blood_all), n_sputum_all=len(sp_all),
                     blood_mean=float(np.mean(d_blood)), sputum_mean=float(np.mean(d_sp)),
                     diff_mean=float(np.mean(diffs)) if diffs else np.nan,
                     blood_all_mean=float(blood_all[fam].mean()),
                     sputum_all_mean=float(sp_all[fam].mean()),
                     all_diff=float(sp_all[fam].mean() - blood_all[fam].mean()),
                     MW_all_p=float(pu),
                     wilcoxon_p=float(pw), paired_t_p=float(pt)))
c1 = pd.DataFrame(rows)
c1.to_csv(INTER + r"\M10C_sputum_blood.csv", index=False)
for _, r in c1.iterrows():
    note(f"   {r['family']}: 全样本 血 {r['blood_all_mean']:+.2f} 痰 {r['sputum_all_mean']:+.2f} "
         f"Δ={r['all_diff']:+.3f} MW p={r['MW_all_p']:.3g} | 配对(n={r['n_pairs']}) "
         f"Δ={r['diff_mean']:+.3f} Wilcoxon p={r['wilcoxon_p']:.3g}")

# ======================================================================
# C2 血 vs 气道（GSE180578，PBMC vs TA，RAW 矩阵，同患者配对）
# ======================================================================
note("\n== C2 血 vs 气道（GSE180578 PBMC vs TA）")
import scanpy as sc
from collections import defaultdict
MT_GENES = {"MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB",
            "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6"}
a = sc.read_h5ad(RAW + r"\GSE180578_COVID19_PBMC_Trachea\GSE180578_cillo_covid19_study_aggregrated_annotated_data.h5ad")
obs = a.obs
hmap = defaultdict(set)
for bc, row in obs[["patient_id", "disease", "sample_type", "timepoint"]].iterrows():
    key = str(bc).split("_")[0]
    hmap[key].add((str(row["patient_id"]), str(row["disease"]),
                   str(row["sample_type"]), row["timepoint"]))
# TA 样本 = GSM 标题含 ETA 的样本目录
OUT_DIR = INTER + r"\gse180578_gex"
import gzip as _gzip
import scipy.io as sio
import io as _io

def read_sample(sp_dir):
    fp_m = os.path.join(sp_dir, "matrix.mtx")
    fp_b = os.path.join(sp_dir, "barcodes.tsv")
    fp_f = os.path.join(sp_dir, "features.tsv")
    if not (os.path.exists(fp_m) and os.path.exists(fp_b) and os.path.exists(fp_f)):
        return None
    bcs = [x.strip() for x in open(fp_b, encoding="utf-8", errors="replace")]
    syms = [x.split("\t")[1].upper() if "\t" in x else x.upper()
            for x in open(fp_f, encoding="utf-8", errors="replace")]
    raw_mtx = open(fp_m, encoding="utf-8", errors="replace").read()
    filtered = "\n".join(l for l in raw_mtx.splitlines() if not l.startswith("%metadata_json"))
    X = sio.mmread(_io.StringIO(filtered)).tocsc()
    if X.shape[0] != len(syms):
        X = X.T.tocsc()
    return bcs, syms, X

# 收集 TA（含 ETA 标题的 GSM）与同患者同时点 PBMC
want_genes = sorted((UP | EX | MT_GENES))
ta_recs = {}     # (patient, tp) -> (n, sum_vector)
pb_recs = {}
for gsm_dir in sorted(os.listdir(OUT_DIR)):
    p = os.path.join(OUT_DIR, gsm_dir)
    if not os.path.isdir(p):
        continue
    is_eta = "ETA" in gsm_dir.upper()
    for samp in sorted(os.listdir(p)):
        sp_dir = os.path.join(p, samp)
        if not os.path.isdir(sp_dir):
            continue
        got = read_sample(sp_dir)
        if got is None:
            continue
        bcs, syms, X = got
        sym_upper = [s.upper() for s in syms]
        gene_pos = {g: [i for i, s in enumerate(sym_upper) if s == g] for g in want_genes}
        lib = np.asarray(X.sum(axis=0)).ravel()
        # 每细胞注释
        cell_ann = {}
        for j, bc in enumerate(bcs):
            hits = hmap.get(bc.split("-")[0])
            if not hits:
                continue
            patients = {x[0] for x in hits}
            if len(patients) != 1:
                continue
            pid, disease, stype, tp = next(iter(hits))
            tps = str(tp).strip().upper()
            tkey = "NA" if tps in ("NA", "NAN", "NONE", "") else str(int(float(tp)))
            cell_ann[j] = (pid, tkey)
        if not cell_ann:
            continue
        idx = np.array(list(cell_ann.keys()), dtype=int)
        # 向量化：逐基因 log1p(CPM1e4)
        vmat = np.zeros((len(want_genes), len(idx)))
        for gi, g in enumerate(want_genes):
            pos = gene_pos.get(g, [])
            if pos:
                vmat[gi] = np.asarray(X[pos, :][:, idx].sum(axis=0)).ravel() / lib[idx] * 1e4
        vmat = np.log1p(vmat)
        ag = defaultdict(lambda: [0, np.zeros(len(want_genes))])
        for k, j in enumerate(idx):
            key = cell_ann[j]
            ag[key][0] += 1
            ag[key][1] += vmat[:, k]
        for (pid, tp), (n, s) in ag.items():
            if n < 30:
                continue
            mu = s / n
            up_m = float(np.mean([mu[gi] for gi, g in enumerate(want_genes) if g in UP]))
            ex_m = float(np.mean([mu[gi] for gi, g in enumerate(want_genes) if g in EX]))
            rec = dict(n_cells=n, UCS_raw=up_m, EIS_raw=ex_m)
            if is_eta:
                ta_recs[(pid, tp)] = rec
            else:
                pb_recs[(pid, tp)] = rec

note(f"   TA 记录 {len(ta_recs)}；PBMC 记录 {len(pb_recs)}")
paired = [(k, ta_recs[k], pb_recs[k]) for k in ta_recs if k in pb_recs]
note(f"   同患者同时点配对 {len(paired)}")
rows2 = []
for fam in ("UCS_raw", "EIS_raw"):
    dta = [t[fam] for _, t, _ in paired]
    dpb = [b[fam] for _, _, b in paired]
    diffs = [t[fam] - b[fam] for _, t, b in paired]
    if len(diffs) >= 4 and any(x != 0 for x in diffs):
        w, pw = wilcoxon(diffs)
        t, pt = ttest_rel(dta, dpb)
        rows2.append(dict(family=fam, n_pairs=len(paired),
                          TA_mean=float(np.mean(dta)), PBMC_mean=float(np.mean(dpb)),
                          diff_mean=float(np.mean(diffs)),
                          wilcoxon_p=float(pw), paired_t_p=float(pt)))
c2 = pd.DataFrame(rows2)
c2.to_csv(INTER + r"\M10C_pbmc_ta.csv", index=False)
for _, r in c2.iterrows():
    note(f"   {r['family']}: TA {r['TA_mean']:+.3f} PBMC {r['PBMC_mean']:+.3f} "
         f"Δ={r['diff_mean']:+.3f} Wilcoxon p={r['wilcoxon_p']:.3g}")
# 全部 TA（未配对）描述性
if ta_recs:
    ta_all = pd.DataFrame(ta_recs.values())
    note(f"   TA 全部（未配对，n={len(ta_all)}）：UCS_raw mean {ta_all['UCS_raw'].mean():.3f}，"
         f"EIS_raw mean {ta_all['EIS_raw'].mean():.3f}")

# ======================================================================
# C3 BALF vs PBMC 跨队列（描述性）
# ======================================================================
note("\n== C3 BALF vs PBMC（GSE145926 vs GSE158055，severe 供者）")
p0 = pd.read_csv(INTER + r"\P0_pseudobulk_arm_scores.csv")
p0["dataset"] = p0["dataset"].astype(str)
sev = p0[p0["condition"] == "COVID_severe"]
rows3 = []
for fam in ("UCS_raw", "EIS_raw"):
    balf = sev[sev["dataset"] == "GSE145926"][fam].dropna()
    pbmc = sev[sev["dataset"] == "GSE158055"][fam].dropna()
    rows3.append(dict(family=fam, BALF_n=len(balf), BALF_mean=float(balf.mean()),
                      PBMC_n=len(pbmc), PBMC_mean=float(pbmc.mean()),
                      diff=float(balf.mean() - pbmc.mean())))
c3 = pd.DataFrame(rows3)
c3.to_csv(INTER + r"\M10C_balf_pbmc.csv", index=False)
for _, r in c3.iterrows():
    note(f"   {r['family']}: BALF {r['BALF_mean']:+.3f}（n={r['BALF_n']}）vs "
         f"PBMC {r['PBMC_mean']:+.3f}（n={r['PBMC_n']}）Δ={r['diff']:+.3f}（跨队列描述性）")
note("\n== C4 肺 Visium：引用 Table S55b（臂评分/双变量 I 组间 Kruskal-Wallis 均不显著；"
     "MDI p=0.715）——无新增计算。")
with open(ROOT + r"\03_LOGS\M10C_compartment_log.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(log))
note(f"\nDONE in {time.time()-t0:.1f}s")
