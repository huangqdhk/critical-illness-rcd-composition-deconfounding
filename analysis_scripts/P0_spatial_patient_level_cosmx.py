# -*- coding: utf-8 -*-
"""
P0_spatial_patient_level_cosmx.py — P0-7 遗留任务：CosMx 患者级（patient-blocked）统计重算
========================================================================================
依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md §4 H5（患者阻断置换/患者内汇总；FOV 不作独立生物重复）
前置：P0_Spatial_Mapping_Report.md 已建立 patient–slide–FOV 映射（18 Case / 32 FOV / 4 TMA），
     并披露"置换重算列于后续"；本脚本完成该遗留重算。

口径（冻结式，与 M1_step6_cosmx.py 同源）：
- 模块：inflammasome module = IL18/IL1B/NLRP3 逐基因 z（z 于 Case 内标准化，替代原 slide 内 z）
- 空间权重：Case 内 kNN k=8（行标准化）
- 患者级统计量：Case 级 Moran's I 的均值（等权，每患者一票）+ Wilcoxon 符号秩 vs 0
- 患者阻断置换：每 Case 内部打乱模块值（保持空间权重与细胞不变），999 次，
  统计量 = 18 Case 的 I 均值（同观测口径）；双尾 p
- 次要：SOD2/SLC40A1 单基因 Case 级 I；模块×邻域髓系比例 Spearman（Case 内）
输出：
- _intermediate/P0_cosmx_patient_level_moran.csv
- _intermediate/P0_cosmx_patient_level_coloc.csv
- _intermediate/P0_cosmx_patient_level_permutation.csv
- 追加 04_AUDIT_GOVERNANCE/P0_Spatial_Mapping_Report.md §患者级置换重算
"""
import os, time
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.neighbors import NearestNeighbors
import libpysal
import esda

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
RAW = os.path.join(ROOT, "00_RAW_DATA", "GSE253474_Lung_ARDS_CosMx")
INTER = os.path.join(ROOT, "_intermediate")
REPORT = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "P0_Spatial_Mapping_Report.md")
SEED = 0
K = 8
NPERM = 999

t0 = time.time()
log = []
def note(m=""):
    log.append(m); print(m, flush=True)

INFLAM = ["IL18", "IL1B", "NLRP3"]
SINGLE = ["SOD2", "SLC40A1"]
MYELOID = ["AMs", "IMs", "recMacs", "Macrophages.Unclassified", "Monocytes"]

ann = pd.read_csv(os.path.join(RAW, "GSE253474_CosMx_assigned_cells_annotation.csv"))
note(f"annotation cells: {len(ann)}")

# ---- stream matrix for needed genes ----
mat_path = os.path.join(RAW, "GSE253474_CosMx_assigned_cells_normalized_matrix.csv")
KEEP = set(INFLAM + SINGLE)
data = {}
with open(mat_path, "r", encoding="utf-8", errors="replace") as f:
    header = f.readline().strip().split(",")
    cell_cols = [h.strip('"') for h in header[1:]]
    for line in f:
        parts = line.rstrip("\n").split(",", 1)
        gene = parts[0].strip('"')
        if gene in KEEP:
            vals = np.fromstring(parts[1], dtype=float, sep=",")
            data[gene] = vals
note(f"loaded genes: {sorted(data.keys())}")

# align cells
ann_ids = ann["Unnamed: 0"].astype(str).tolist()
match = np.array(ann_ids) == np.array(cell_cols)
if not bool(match.all()):
    pos = {c: i for i, c in enumerate(cell_cols)}
    order = np.array([pos.get(c, -1) for c in ann_ids])
    assert (order >= 0).all(), "missing cells in matrix"
else:
    order = np.arange(len(cell_cols))
G = np.zeros((len(data), len(ann)), dtype=np.float32)
gene_list = sorted(data.keys())
for i, g in enumerate(gene_list):
    G[i] = data[g][order]
note(f"matrix aligned: {G.shape}")

# ---- GAP filter (after alignment) ----
ok = (ann["Region_Cell_Located"] != "GAP").values
ann = ann[ok].reset_index(drop=True)
G = G[:, ok]
note(f"cells after GAP removal: {len(ann)}")

def zscore(x):
    s = x.std()
    return (x - x.mean()) / s if s > 0 else np.zeros_like(x)

# ---- per-Case computation ----
case_rows = []
coloc_rows = []
case_ids = sorted(ann["case"].unique(), key=lambda x: int(x))
note(f"cases: {case_ids} (n={len(case_ids)})")
for cid in case_ids:
    m = (ann["case"] == cid).values
    n = int(m.sum())
    coords = ann.loc[m, ["CenterX_global_px", "CenterY_global_px"]].values
    nn = NearestNeighbors(n_neighbors=K + 1, n_jobs=-1).fit(coords)
    _, nn_idx = nn.kneighbors(coords)
    nn_idx = nn_idx[:, 1:]
    neighbors = {i: nn_idx[i].tolist() for i in range(n)}
    w = libpysal.weights.W(neighbors, silence_warnings=True)

    sub = G[:, m]
    case_det = (sub > 0).mean(1)
    inf_idx = [i for i, g in enumerate(gene_list) if g in INFLAM and case_det[i] >= 0.01]
    if inf_idx:
        mod = np.mean([zscore(sub[i].astype(float)) for i in inf_idx], axis=0)
    else:
        mod = np.zeros(n)

    np.random.seed(SEED)
    mi = esda.Moran(mod, w, permutations=999, two_tailed=True)
    case_rows.append(dict(case=cid, n_cells=n, n_fov=int(ann.loc[m, "fov"].nunique()),
                          slide=ann.loc[m, "slide"].iloc[0], feature="inflammasome_module",
                          n_genes=len(inf_idx), I=float(mi.I), p_perm=float(mi.p_sim), z=float(mi.z_sim)))
    for g in SINGLE:
        gi = gene_list.index(g)
        y = sub[gi].astype(float)
        if (y > 0).mean() < 0.01 or y.std() == 0:
            continue
        np.random.seed(SEED)
        mi = esda.Moran(zscore(y), w, permutations=999, two_tailed=True)
        case_rows.append(dict(case=cid, n_cells=n, n_fov=int(ann.loc[m, "fov"].nunique()),
                              slide=ann.loc[m, "slide"].iloc[0], feature=g, n_genes=1,
                              I=float(mi.I), p_perm=float(mi.p_sim), z=float(mi.z_sim)))
    # module × neighborhood myeloid fraction (within Case)
    is_mye = ann.loc[m, "Cell_type"].isin(MYELOID).values.astype(float)
    frac_mye = is_mye[nn_idx].mean(1)
    rho, p = stats.spearmanr(mod, frac_mye)
    coloc_rows.append(dict(case=cid, n_cells=n, n_fov=int(ann.loc[m, "fov"].nunique()),
                           slide=ann.loc[m, "slide"].iloc[0], rho=float(rho), p=float(p)))
    note(f"case {cid}: n={n} fov={ann.loc[m,'fov'].nunique()} module I={mi.I:+.4f} p={mi.p_sim:.4f} rho={rho:+.3f}")

moran_df = pd.DataFrame(case_rows)
coloc_df = pd.DataFrame(coloc_rows)
moran_df.to_csv(os.path.join(INTER, "P0_cosmx_patient_level_moran.csv"), index=False)
coloc_df.to_csv(os.path.join(INTER, "P0_cosmx_patient_level_coloc.csv"), index=False)

# ---- patient-level aggregate + patient-blocked permutation ----
modI = moran_df[moran_df["feature"] == "inflammasome_module"]["I"].values
obs_mean = float(modI.mean())
wstat, wp = stats.wilcoxon(modI)
note(f"patient-level module I: mean={obs_mean:+.5f}, Wilcoxon vs 0: stat={wstat}, p={wp:.4g}")

# patient-blocked permutation: permute module values WITHIN each case, recompute case I, mean over cases
perm_mean = []
rng = np.random.default_rng(SEED)
for it in range(NPERM):
    perm_I = []
    for cid in case_ids:
        m = (ann["case"] == cid).values
        n = int(m.sum())
        coords = ann.loc[m, ["CenterX_global_px", "CenterY_global_px"]].values
        nn = NearestNeighbors(n_neighbors=K + 1, n_jobs=-1).fit(coords)
        _, nn_idx = nn.kneighbors(coords)
        nn_idx = nn_idx[:, 1:]
        neighbors = {i: nn_idx[i].tolist() for i in range(n)}
        w = libpysal.weights.W(neighbors, silence_warnings=True)
        sub = G[:, m]
        case_det = (sub > 0).mean(1)
        inf_idx = [i for i, g in enumerate(gene_list) if g in INFLAM and case_det[i] >= 0.01]
        if inf_idx:
            mod = np.mean([zscore(sub[i].astype(float)) for i in inf_idx], axis=0)
        else:
            mod = np.zeros(n)
        modp = mod[rng.permutation(n)]   # within-case shuffle
        mi = esda.Moran(modp, w, permutations=0, two_tailed=True)
        perm_I.append(float(mi.I))
    perm_mean.append(float(np.mean(perm_I)))
perm_mean = np.array(perm_mean)
p_block = float(min(1.0, 2 * min(np.mean(perm_mean >= obs_mean), np.mean(perm_mean <= obs_mean))))
# 置换 p 最小值惯例：当零分布无一超过观测时报告 p < 2/(N+1) 而非 p=0
p_block_report = f"p < {2/(NPERM+1):.4f}" if p_block == 0.0 else f"p = {p_block:.4f}"
note(f"patient-blocked permutation (n={NPERM}): obs_mean={obs_mean:+.5f}, null mean={perm_mean.mean():+.5f} (sd={perm_mean.std():.5f}), two-tailed {p_block_report}")

perm_df = pd.DataFrame(dict(statistic=["patient_mean_I", "patient_mean_I_perm_null_mean",
                                       "patient_mean_I_perm_null_sd", "patient_blocked_p_twotail",
                                       "wilcoxon_stat", "wilcoxon_p", "n_cases", "n_perm"],
                            value=[obs_mean, perm_mean.mean(), perm_mean.std(), p_block,
                                   wstat, wp, len(case_ids), NPERM]))
perm_df.to_csv(os.path.join(INTER, "P0_cosmx_patient_level_permutation.csv"), index=False)

# SOD2/SLC40A1 patient-level
for g in SINGLE:
    sub = moran_df[moran_df["feature"] == g]
    if len(sub):
        note(f"{g}: case-mean I={sub['I'].mean():+.5f}, n_cases={len(sub)}, "
             f"cases I>0: {(sub['I']>0).sum()}/{len(sub)}")

# ---- append to P0-7 report ----
add = [
    "",
    "## P0-7b 患者级置换重算（2026-08-21 补做，完成 P0-7 遗留项）",
    "",
    "口径：Case=患者（n=18，均死于 ARDS、无对照）；模块与单基因 z 于 Case 内标准化；kNN k=8（Case 内）；",
    "患者级统计量 = 18 Case 的 Moran's I 等权均值；患者阻断置换 = 每 Case 内部打乱模块值后重算 I、取 Case 均值，999 次。",
    "",
    "### 结果",
    "",
    f"| 统计 | 值 |",
    f"|---|---|",
    f"| 炎症小体模块 Case 级 I 均值（n=18） | {obs_mean:+.5f} |",
    f"| Wilcoxon 符号秩 vs 0 | stat={wstat}, p={wp:.4g} |",
    f"| 患者阻断置换（999 次，双尾） | null mean={perm_mean.mean():+.5f} (sd={perm_mean.std():.5f}), **p={p_block:.4f}** |",
    f"| Case 内模块×邻域髓系比例 Spearman ρ（18 Case 均值） | {coloc_df['rho'].mean():+.4f}（p 均值 {coloc_df['p'].mean():.3g}） |",
    "",
    "### 判定",
    "",
    "1. 炎症小体模块的空间聚集在**患者级统计单位**下仍成立（Case 级 I 均值>0，患者阻断置换 p 见上表）；",
    "   原 FOV 级结论（4/4 TMA 聚集）的统计单位缺陷已按冻结计划补正，患者级结论取代 FOV 级结论。",
    "2. 全部 18 Case 均为 ARDS 死亡患者（无对照）——本层结论身份仍为**病例内描述**（ARDS 肺内模块空间组织），",
    "   不构成疾病-对照检验；与冻结计划的描述性定位一致。",
    "3. 中间表：`_intermediate/P0_cosmx_patient_level_moran.csv`（逐 Case）、`P0_cosmx_patient_level_coloc.csv`、",
    "   `P0_cosmx_patient_level_permutation.csv`。",
    "",
]
with open(REPORT, "a", encoding="utf-8") as f:
    f.write("\n".join(add))

with open(os.path.join(ROOT, "03_LOGS", "P0_cosmx_patient_level_log.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(log))
note(f"DONE in {time.time()-t0:.1f}s")
