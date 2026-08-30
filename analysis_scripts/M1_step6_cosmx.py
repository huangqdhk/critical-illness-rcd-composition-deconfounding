# -*- coding: utf-8 -*-
"""
M1 Step 6: CosMx (GSE253474) supporting analysis at single-cell resolution.

Panel reality (checked): 960 real genes; 80-gene overlap = 5 execution
(IL18, IL1B, NLRP3, SLC40A1, SOD2) + 9 non-arm; upstream = 0/30. TFRC absent,
TNFSF13B present.  All 18 cases died of ARDS (no non-ARDS control).

Analyses (pre-registered C1-C4):
  C1: execution-module & single-gene spatial autocorrelation (Moran's I, per slide)
  C2: module score vs neighborhood myeloid fraction (Spearman, per slide)
  C3: module score across viral-region niches (POS/ADJ/NEG) per cell type
  C4: TNFSF13B cellular source + neighborhood enrichment
"""
import os, csv
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.neighbors import NearestNeighbors
import libpysal
import esda

RAW = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\00_RAW_DATA\GSE253474_Lung_ARDS_CosMx"
OUT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\02_SUPPLEMENTARY_TABLES\SUPPLEMENTARY_Tables_CSV"
SEED = 0
K = 8

# keep-set: 14 overlapping + TNFSF13B + context genes
OVERLAP = ["ATG5", "BECN1", "DDIT3", "HIF1A", "IL18", "IL1B", "MAP1LC3B",
           "NFKB1", "NLRP3", "PPARG", "RELA", "SLC40A1", "SOD2", "SQSTM1"]
CONTEXT = ["TNFSF13B", "CD68", "CD163", "LYZ", "S100A8", "S100A9", "CD14",
           "FCGR3A", "ITGAX", "EPCAM", "KRT5", "KRT8", "PDPN", "COL1A1",
           "VWF", "PECAM1", "MRC1", "MARCO", "FABP4", "SPP1", "C1QA",
           "C1QB", "APOE", "CCL2", "CXCL10", "ISG15", "IFI27", "MX1", "OAS1"]
KEEP = sorted(set(OVERLAP + CONTEXT))
print("keep set size:", len(KEEP))

# ---- load annotation ----
ann = pd.read_csv(os.path.join(RAW, "GSE253474_CosMx_assigned_cells_annotation.csv"))
print("cells:", len(ann))
MYELOID = ["AMs", "IMs", "recMacs", "Macrophages.Unclassified", "Monocytes"]

# ---- stream matrix rows for KEEP genes ----
mat_path = os.path.join(RAW, "GSE253474_CosMx_assigned_cells_normalized_matrix.csv")
data = {}
with open(mat_path, "r", encoding="utf-8", errors="replace") as f:
    header = f.readline().strip().split(",")
    col_names = [h.strip('"') for h in header]
    # first col is gene id col
    cell_cols = col_names[1:]
    for line in f:
        parts = line.rstrip("\n").split(",", 1)
        gene = parts[0].strip('"')
        if gene in KEEP:
            vals = np.fromstring(parts[1], dtype=float, sep=",")
            if len(vals) != len(cell_cols):
                print("row length mismatch", gene, len(vals), len(cell_cols))
                continue
            data[gene] = vals
            print("loaded", gene)
print("loaded genes:", sorted(data.keys()))
missing = [g for g in KEEP if g not in data]
print("not found in matrix (will be dropped):", missing)

# align cells: matrix columns should match annotation cell_ID order?  Check via index col
# annotation 'Unnamed: 0' == cell id string like TMA_A_8_1_2
ann_ids = ann["Unnamed: 0"].astype(str).tolist()
match = np.array(ann_ids) == np.array(cell_cols)
print("annotation row order matches matrix col order:", bool(match.all()))
if not bool(match.all()):
    # map by id
    pos = {c: i for i, c in enumerate(cell_cols)}
    order = np.array([pos.get(c, -1) for c in ann_ids])
    assert (order >= 0).all(), "missing cells in matrix"
else:
    order = np.arange(len(cell_cols))

G = np.zeros((len(data), len(ann)), dtype=np.float32)
gene_list = sorted(data.keys())
for i, g in enumerate(gene_list):
    G[i] = data[g][order]

# ---- QC: drop GAP-region cells for main analyses (keep all for composition) ----
ann["region3"] = ann["Region_Cell_Located"]
ok = (ann["region3"] != "GAP").values
print("cells after GAP removal:", ok.sum())

# detection rates
det = (G > 0).mean(1)
det_df = pd.DataFrame({"gene": gene_list, "detection_rate_all": det,
                       "detection_rate_nonGAP": (G[:, ok] > 0).mean(1)})
det_df.to_csv(os.path.join(OUT, "Table_S55i_M1_CosMx_Gene_Detection.csv"), index=False)
print(det_df.to_string())

# ---- per-slide analysis ----
def zscore(x):
    s = x.std()
    return (x - x.mean()) / s if s > 0 else np.zeros_like(x)

slide_rows = []
coloc_rows = []
module_vals = np.full(len(ann), np.nan)
for slide in sorted(ann["slide"].unique()):
    m = (ann["slide"] == slide).values & ok
    n = m.sum()
    coords = ann.loc[m, ["CenterX_global_px", "CenterY_global_px"]].values
    nn = NearestNeighbors(n_neighbors=K + 1, n_jobs=-1).fit(coords)
    _, nn_idx = nn.kneighbors(coords)
    nn_idx = nn_idx[:, 1:]
    neighbors = {i: nn_idx[i].tolist() for i in range(n)}
    w = libpysal.weights.W(neighbors, silence_warnings=True)

    # execution module: mean z over inflammasome genes with det>=1% in slide
    sub = G[:, m]
    slide_det = (sub > 0).mean(1)
    inf_genes = [i for i, g in enumerate(gene_list) if g in ("IL18", "IL1B", "NLRP3")
                 and slide_det[i] >= 0.01]
    if inf_genes:
        mod = np.mean([zscore(sub[i].astype(float)) for i in inf_genes], axis=0)
    else:
        mod = np.zeros(n)
    module_vals[m] = mod

    np.random.seed(SEED)
    mi = esda.Moran(mod, w, permutations=999, two_tailed=True)
    slide_rows.append(dict(slide=slide, n_cells=n, feature="inflammasome_module",
                           n_genes=len(inf_genes), I=mi.I, p_perm=mi.p_sim, z=mi.z_sim))
    for g in ["SOD2", "SLC40A1", "TNFSF13B", "NLRP3", "IL1B", "IL18"]:
        gi = gene_list.index(g)
        y = sub[gi].astype(float)
        if (y > 0).mean() < 0.01 or y.std() == 0:
            continue
        np.random.seed(SEED)
        mi = esda.Moran(zscore(y), w, permutations=999, two_tailed=True)
        slide_rows.append(dict(slide=slide, n_cells=n, feature=g, n_genes=1,
                               I=mi.I, p_perm=mi.p_sim, z=mi.z_sim))

    # C2: neighborhood myeloid fraction vs module score
    is_mye = ann.loc[m, "Cell_type"].isin(MYELOID).values.astype(float)
    frac_mye = is_mye[nn_idx].mean(1)
    rho, p = stats.spearmanr(mod, frac_mye)
    # bivariate moran (module x myeloid fraction), manual permutation
    Wy = w.sparse.dot(frac_mye)
    I_bv = float((mod @ Wy) / (mod @ mod))
    rng = np.random.default_rng(SEED)
    sims = np.empty(999)
    for i in range(999):
        yp = frac_mye[rng.permutation(n)]
        sims[i] = float((mod @ w.sparse.dot(yp)) / (mod @ mod))
    p_bv = float(min(1.0, 2 * min(np.mean(sims >= I_bv), np.mean(sims <= I_bv))))
    coloc_rows.append(dict(slide=slide, n_cells=n, rho=rho, p=rho_p if False else p,
                           I_bv=I_bv, p_bv=p_bv))
    print(f"{slide}: module I={mi.I:+.4f} p={mi.p_sim:.4f} | rho(mye)={rho:+.3f} p={p:.3g}")

pd.DataFrame(slide_rows).to_csv(os.path.join(OUT, "Table_S55j_M1_CosMx_Moran.csv"), index=False)
pd.DataFrame(coloc_rows).to_csv(os.path.join(OUT, "Table_S55k_M1_CosMx_Myeloid_Coloc.csv"), index=False)

# ---- C2b: cell-type composition of module-high cells ----
ann["module_score"] = module_vals
high = {}
for slide in ann["slide"].unique():
    m = ann["slide"] == slide
    thr = ann.loc[m, "module_score"].quantile(0.75)
    ann.loc[m, "module_high"] = ann.loc[m, "module_score"] >= thr
comp_rows = []
for ct, g in ann[ok].groupby("Cell_type"):
    comp_rows.append(dict(cell_type=ct, n_cells=len(g),
                          module_mean=g["module_score"].mean(),
                          frac_module_high=g["module_high"].mean()))
comp = pd.DataFrame(comp_rows).sort_values("module_mean", ascending=False)
comp.to_csv(os.path.join(OUT, "Table_S55l_M1_CosMx_Module_by_CellType.csv"), index=False)
print(comp.head(12).to_string())

# ---- C3: viral-region niche comparison (per cell type) ----
rows3 = []
ct_list = sorted(ct for ct in ann["Cell_type"].unique() if isinstance(ct, str))
for ct in ct_list:
    g = ann[(ann["Cell_type"] == ct) & ok]
    if len(g) < 50:
        continue
    for reg in ["POS", "ADJ", "NEG"]:
        a = g.loc[g["region3"] == reg, "module_score"].dropna()
        if len(a) >= 10:
            rows3.append(dict(cell_type=ct, region=reg, n=len(a), mean=a.mean(), median=a.median()))
        else:
            rows3.append(dict(cell_type=ct, region=reg, n=len(a), mean=np.nan, median=np.nan))
    # POS vs NEG test
    a = g.loc[g["region3"] == "POS", "module_score"].dropna()
    b = g.loc[g["region3"] == "NEG", "module_score"].dropna()
    if len(a) >= 10 and len(b) >= 10:
        u, p = stats.mannwhitneyu(a, b)
        rows3.append(dict(cell_type=ct, region="POS_vs_NEG", n=len(a) + len(b),
                          mean=u, median=p, test="MWU_u_and_p"))
niche = pd.DataFrame(rows3)
niche.to_csv(os.path.join(OUT, "Table_S55m_M1_CosMx_ViralRegion_Niche.csv"), index=False)

# ---- C4: TNFSF13B source + neighborhood ----
gi = gene_list.index("TNFSF13B")
ann["tnfsf13b"] = G[gi]
tn_pos = (ann["tnfsf13b"] > 0) & ok
src_rows = []
for ct, g in ann[ok].groupby("Cell_type"):
    src_rows.append(dict(cell_type=ct, n_cells=len(g), frac_TNFSF13B_pos=(g["tnfsf13b"] > 0).mean(),
                         mean_logCPM=g["tnfsf13b"].mean()))
src = pd.DataFrame(src_rows).sort_values("frac_TNFSF13B_pos", ascending=False)
src.to_csv(os.path.join(OUT, "Table_S55n_M1_CosMx_TNFSF13B_Source.csv"), index=False)
print(src.head(10).to_string())

# neighborhood enrichment of TNFSF13B+ cells for macrophage subtypes (per slide)
enr_rows = []
for slide in ann["slide"].unique():
    m = (ann["slide"] == slide).values & ok
    coords = ann.loc[m, ["CenterX_global_px", "CenterY_global_px"]].values
    n = m.sum()
    nn = NearestNeighbors(n_neighbors=K + 1, n_jobs=-1).fit(coords)
    _, nn_idx = nn.kneighbors(coords)
    nn_idx = nn_idx[:, 1:]
    pos = (ann.loc[m, "tnfsf13b"] > 0).values
    if pos.sum() < 20:
        continue
    for ct in MYELOID + ["Fibroblasts", "AT2", "AT1"]:
        lab = ann.loc[m, "Cell_type"].isin([ct]).values.astype(float)
        obs = lab[nn_idx[pos]].mean()
        # permutation
        rng = np.random.default_rng(SEED)
        perms = np.array([lab[rng.permutation(n)[nn_idx[pos]]].mean()
                          for _ in range(999)])
        z = (obs - perms.mean()) / (perms.std() + 1e-12)
        p = (1 + (perms >= obs).sum()) / 1000
        enr_rows.append(dict(slide=slide, n_TNFSF13B_pos=int(pos.sum()), neighbor_type=ct,
                             obs_frac=obs, perm_mean=perms.mean(), z=z, p=p))
enr = pd.DataFrame(enr_rows)
enr.to_csv(os.path.join(OUT, "Table_S55o_M1_CosMx_TNFSF13B_Neighborhood.csv"), index=False)
print(enr.to_string())
print("step 6 done")
