# -*- coding: utf-8 -*-
"""
M10B_step1_gse180578_extract.py — GSE180578 RAW GEX 解包 + 条形码映射 + 供者×细胞型 arm 评分
========================================================================================
M10 检验 B 队列之一。RAW 矩阵 = CellRanger 3.1.0 全基因组（33,538 基因）；
h5ad = 逐细胞注释（cell_types 11 类、patient_id、disease、sample_type、timepoint、outcome）。
口径（沿用 P0 供者级 pseudobulk）：
  - 每细胞 log1p(CPM*1e4)（与 combined_processed.h5ad 的 log1p 标准化同尺）
  - arm 评分 = 供者×细胞型内 arm 基因均值 -> 细胞型内 z -> Welch t（BH 族 = celltype x score）
  - donor x celltype >= 50 细胞、双侧各 >= 3 供者；主对比 covid(timepoint=1) vs healthy；PBMC 区室
输出：
  _intermediate/M10B_gse180578_arm_scores.csv（供者级）
  _intermediate/M10B_gse180578_arm_tests.csv（Welch t + BH）
"""
import os, sys, glob, tarfile, warnings, time
import numpy as np
import pandas as pd
import scipy.sparse as sp

warnings.filterwarnings("ignore")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
RAW = ROOT + r"\00_RAW_DATA"
INTER = ROOT + r"\_intermediate"
GSE_RAW = RAW + r"\GSE180578_COVID19_PBMC_Trachea\GSE180578_RAW"
OUT_DIR = INTER + r"\gse180578_gex"
H5AD_F = RAW + r"\GSE180578_COVID19_PBMC_Trachea\GSE180578_cillo_covid19_study_aggregrated_annotated_data.h5ad"
MANIFEST = ROOT + r"\04_AUDIT_GOVERNANCE\Mitoxyperilysis_Gene_Manifest_v1.0.csv"
MT_GENES = {"MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB",
            "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6"}

sys.path.insert(0, ROOT)
import mdi_lib as L

t0 = time.time()
log = []
def note(m=""):
    log.append(m); print(m, flush=True)

def zc(s):
    s = np.asarray(s, float); sd = s.std(ddof=1)
    return (s - s.mean()) / sd if sd > 0 else s * 0

def bh(pvals):
    p = np.asarray(pvals, float); out = np.full_like(p, np.nan)
    m = ~np.isnan(p); ps = p[m]
    if len(ps) == 0: return out
    order = np.argsort(ps); n = len(ps); adj = np.empty(n)
    for k in range(n - 1, -1, -1):
        adj[order[k]] = min(ps[order[k]] * n / (k + 1), adj[order[k + 1]] if k + 1 < n else 1.0)
    out[m] = adj; return out

# ---------------- 1. 解包 48 个 GEX tar.gz ----------------
note("== 1. 解包 GEX tars")
os.makedirs(OUT_DIR, exist_ok=True)
tars = sorted(glob.glob(GSE_RAW + r"\*GEX*.tar.gz"))
note(f"   GEX tars: {len(tars)}")
for i, f in enumerate(tars):
    d = os.path.join(OUT_DIR, os.path.basename(f).replace(".tar.gz", ""))
    if not os.path.isdir(d):
        with tarfile.open(f, "r:gz") as t:
            t.extractall(d)
    note(f"   [{i+1}/{len(tars)}] {os.path.basename(f)}")

# ---------------- 2. h5ad 注释映射（barcode -> patient/celltype/condition） ----------------
note("\n== 2. h5ad 注释映射")
import scanpy as sc
a = sc.read_h5ad(H5AD_F)
obs = a.obs
from collections import defaultdict
hmap = defaultdict(set)
for bc, row in obs[["patient_id", "cell_types", "disease", "sample_type", "timepoint"]].iterrows():
    key = str(bc).split("_")[0]           # 去掉 "_1_1" 后缀
    hmap[key].add((str(row["patient_id"]), str(row["cell_types"]),
                   str(row["disease"]), str(row["sample_type"]),
                   row["timepoint"]))
n_amb = sum(1 for k, v in hmap.items() if len({x[0] for x in v}) > 1)
note(f"   h5ad 去后缀 barcode {len(hmap)}；跨患者碰撞 {n_amb}（{n_amb/len(hmap):.1%}）——歧义细胞剔除")

# ---------------- 3. manifest 基因（arm + MT） ----------------
man = pd.read_csv(MANIFEST)
up_genes = set(man.loc[man["arm"] == "upstream_collapse", "hgnc_symbol"].astype(str).str.upper())
ex_genes = set(man.loc[man["arm"] == "execution_induction", "hgnc_symbol"].astype(str).str.upper())
note(f"   arm 基因：上游 {len(up_genes)} / 执行 {len(ex_genes)}")

# ---------------- 4. 逐样本 RAW 矩阵 -> 供者×细胞型 arm 评分 ----------------
note("\n== 4. RAW 矩阵聚合")
MIN_CELLS, MIN_DONORS = 50, 3
recs = []
stats = []
for gsm_dir in sorted(os.listdir(OUT_DIR)):
    p = os.path.join(OUT_DIR, gsm_dir)
    if not os.path.isdir(p):
        continue
    for samp in sorted(os.listdir(p)):
        sp_dir = os.path.join(p, samp)
        if not os.path.isdir(sp_dir):
            continue
        def _find(name):
            for cand in (name, name + ".gz"):
                fp = os.path.join(sp_dir, cand)
                if os.path.exists(fp):
                    return fp
            return None
        mtx_f = _find("matrix.mtx")
        bc_f = _find("barcodes.tsv")
        ft_f = _find("features.tsv")
        if not (mtx_f and bc_f and ft_f):
            continue
        import gzip as _gzip
        def _open_lines(fp):
            if fp.endswith(".gz"):
                return [x.decode().strip() for x in _gzip.open(fp, "rb")]
            return [x.strip() for x in open(fp, encoding="utf-8", errors="replace")]
        bcs = _open_lines(bc_f)
        fts = [x.split("\t") for x in _open_lines(ft_f)]
        syms = [f[1].upper() if len(f) > 1 else f[0].upper() for f in fts]
        import scipy.io as sio
        import io as _io
        if mtx_f.endswith(".gz"):
            raw_mtx = _gzip.open(mtx_f, "rt", encoding="utf-8", errors="replace").read()
        else:
            raw_mtx = open(mtx_f, encoding="utf-8", errors="replace").read()
        # CellRanger 3.1 的 %metadata_json 注释行与 scipy fast-mmread 不兼容，过滤后读取
        filtered = "\n".join(l for l in raw_mtx.splitlines() if not l.startswith("%metadata_json"))
        X = sio.mmread(_io.StringIO(filtered)).tocsc()
        if X.shape[0] != len(syms):
            X = X.T.tocsc()
        lib = np.asarray(X.sum(axis=0)).ravel()          # 每细胞文库
        want = [i for i, s in enumerate(syms) if s in up_genes or s in ex_genes or s in MT_GENES]
        Xw = X[want, :].tocsc()
        up_idx = [pos for pos, i in enumerate(want) if syms[i] in up_genes]
        ex_idx = [pos for pos, i in enumerate(want) if syms[i] in ex_genes]
        up_nomt = [pos for pos, i in enumerate(want) if syms[i] in up_genes and syms[i] not in MT_GENES]
        # 细胞 -> 注释（歧义 barcode：跨患者碰撞则剔除该细胞）
        cell_map = {}
        n_hit, n_miss, n_amb_drop = 0, 0, 0
        for j, bc in enumerate(bcs):
            key = bc.split("-")[0]
            hits = hmap.get(key)
            if not hits:
                n_miss += 1
                continue
            patients = {x[0] for x in hits}
            if len(patients) != 1:
                n_amb_drop += 1
                continue
            cell_map[j] = next(iter(hits)); n_hit += 1
        stats.append(dict(gsm=gsm_dir, sample=samp, n_cells=X.shape[1],
                          n_barcode_hit=n_hit, n_miss=n_miss, n_amb_drop=n_amb_drop))
        # 聚合：仅 PBMC 区室
        ag = {}
        for j, ann in cell_map.items():
            pid, ct, disease, stype, tp = ann
            if stype != "PBMC":
                continue
            l = lib[j]
            if l <= 0:
                continue
            v = np.asarray(Xw[:, j].todense()).ravel() / l * 1e4
            v = np.log1p(v)
            k = (pid, ct, disease, tp)
            if k not in ag:
                ag[k] = [0, np.zeros(len(want))]
            ag[k][0] += 1
            ag[k][1] += v
        for (pid, ct, disease, tp), (n, acc) in ag.items():
            if n < MIN_CELLS:
                continue
            mu = acc / n
            recs.append(dict(cohort="GSE180578", sampleID=pid, cell_type=ct,
                             condition=disease, timepoint=("NA" if pd.isna(tp) else str(tp)),
                             n_cells=n,
                             UCS_raw=float(mu[up_idx].mean()) if up_idx else np.nan,
                             EIS_raw=float(mu[ex_idx].mean()) if ex_idx else np.nan,
                             UCS_nomt_raw=float(mu[up_nomt].mean()) if up_nomt else np.nan))
note(f"   样本级命中率：中位 {np.median([s['n_barcode_hit']/max(1,s['n_cells']) for s in stats]):.0%}")
pd.DataFrame(stats).to_csv(INTER + r"\M10B_gse180578_barcode_stats.csv", index=False)

if not recs:
    note("[ERROR] 无任何 barcode 命中——终止（检查 h5ad/RAW 条形码格式）")
    raise SystemExit(1)
df = pd.DataFrame(recs)
note(f"   供者×细胞型记录 {len(df)}；供者 {df['sampleID'].nunique()}")
# 去重（同供者同细胞型可能跨 GSM 重复出现——同一时点样本只计一次；若同供者同细胞型出现多行则合并求和）
grp = df.groupby(["sampleID", "cell_type", "condition", "timepoint"], as_index=False).agg(
    n_cells=("n_cells", "sum"),
    UCS_raw=("UCS_raw", lambda s: np.average(s, weights=df.loc[s.index, "n_cells"])),
    EIS_raw=("EIS_raw", lambda s: np.average(s, weights=df.loc[s.index, "n_cells"])),
    UCS_nomt_raw=("UCS_nomt_raw", lambda s: np.average(s, weights=df.loc[s.index, "n_cells"])))
df = grp
note(f"   合并后记录 {len(df)}")

# ---------------- 5. 细胞型内 z + Welch t（主对比 covid vs healthy） ----------------
note("\n== 5. 供者级检验")
CASE, CTRL = "covid", "healthy"
# 主对比：covid 仅保留 timepoint==1（供者独立性）；healthy 全部
tp = df["timepoint"].astype(str)
mask_case = (df["condition"] == CASE) & (tp == "1")
mask_ctrl = df["condition"] == CTRL
d = df[mask_case | mask_ctrl].copy()
d["arm_z"] = np.nan
for ct in d["cell_type"].unique():
    m = d["cell_type"] == ct
    d.loc[m, "UCS"] = zc(d.loc[m, "UCS_raw"])
    d.loc[m, "EIS"] = zc(d.loc[m, "EIS_raw"])
    d.loc[m, "MDI"] = d.loc[m, "EIS"] - d.loc[m, "UCS"]
    d.loc[m, "UCS_nomt"] = zc(d.loc[m, "UCS_nomt_raw"])
    d.loc[m, "MDI_nomt"] = d.loc[m, "EIS"] - d.loc[m, "UCS_nomt"]
tests = []
from scipy.stats import ttest_ind
for ct in d["cell_type"].unique():
    sub = d[d["cell_type"] == ct]
    for score in ("UCS", "EIS", "MDI", "MDI_nomt"):
        a_ = sub.loc[sub["condition"] == CASE, score].values
        b_ = sub.loc[sub["condition"] == CTRL, score].values
        if len(a_) >= MIN_DONORS and len(b_) >= MIN_DONORS:
            t, p = ttest_ind(a_, b_, equal_var=False)
            tests.append(dict(cohort="GSE180578", cell_type=ct, score=score,
                              n_case=len(a_), n_ctrl=len(b_),
                              diff=float(a_.mean() - b_.mean()), p=float(p)))
tests = pd.DataFrame(tests)
tests["BH_q"] = bh(tests["p"].values) if len(tests) else np.nan
note("   检验结果：")
for _, r in tests.iterrows():
    note(f"     {r['cell_type']:<16} {r['score']:<9} diff={r['diff']:+.3f} p={r['p']:.3g} q={r['BH_q']:.3g}")

df.to_csv(INTER + r"\M10B_gse180578_arm_scores.csv", index=False)
tests.to_csv(INTER + r"\M10B_gse180578_arm_tests.csv", index=False)
with open(ROOT + r"\03_LOGS\M10B_gse180578_log.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(log))
note(f"\nDONE in {time.time()-t0:.1f}s")
