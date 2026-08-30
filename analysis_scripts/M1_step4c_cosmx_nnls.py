# -*- coding: utf-8 -*-
"""
M1 Step 4c: Validation deconvolution with lung-specific CosMx signatures
(27 cell types, mean logCPM per type) on shared genes (~940).

Cross-platform caveat documented; used only as validation reference for the
scRNA-signature NNLS (step 4b) and for lung-structural cell-type mapping.
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy.optimize import nnls
from scipy.sparse import issparse

RAW = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\00_RAW_DATA\GSE253474_Lung_ARDS_CosMx"
INT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate"
OUT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\02_SUPPLEMENTARY_TABLES\SUPPLEMENTARY_Tables_CSV"

# ---- load full CosMx matrix (960 rows x 98850 cells) ----
mat_path = os.path.join(RAW, "GSE253474_CosMx_assigned_cells_normalized_matrix.csv")
genes = []
rows_list = []
with open(mat_path, "r", encoding="utf-8", errors="replace") as f:
    header = f.readline().strip().split(",")
    cell_cols = [h.strip('"') for h in header[1:]]
    for line in f:
        parts = line.rstrip("\n").split(",", 1)
        genes.append(parts[0].strip('"'))
        rows_list.append(np.fromstring(parts[1], dtype=np.float32, sep=","))
G = np.vstack(rows_list)
print("CosMx full matrix:", G.shape, "n genes:", len(genes))

ann = pd.read_csv(os.path.join(RAW, "GSE253474_CosMx_assigned_cells_annotation.csv"))
ann_ids = ann["Unnamed: 0"].astype(str).tolist()
match = np.array(ann_ids) == np.array(cell_cols)
print("order match:", bool(match.all()))
if not bool(match.all()):
    pos = {c: i for i, c in enumerate(cell_cols)}
    order = np.array([pos[c] for c in ann_ids])
    G = G[:, order]

# mean per cell type
types = sorted(t for t in ann["Cell_type"].unique() if isinstance(t, str))
means = np.zeros((len(genes), len(types)))
for j, t in enumerate(types):
    m = (ann["Cell_type"] == t).values
    means[:, j] = G[:, m].mean(1)
means_df = pd.DataFrame(means, index=genes, columns=types)
means_df.to_csv(os.path.join(INT, "M1_cosmx_type_means.csv"))
print("saved M1_cosmx_type_means.csv")

# ---- NNLS deconvolution of Visium on shared genes ----
adata = ad.read_h5ad(os.path.join(INT, "M1_visium_merged.h5ad"))
adata.var_names_make_unique()
vgenes = list(adata.var_names)
shared = sorted(set(genes) & set(vgenes))
print("shared genes:", len(shared))
S = means_df.loc[shared].values.astype(np.float64)
S = S / (S.sum(0, keepdims=True) + 1e-12)  # normalize signature per type to composition-like scale
idx = [vgenes.index(g) for g in shared]

props_all = []
for sec in sorted(adata.obs["section"].unique()):
    m = (adata.obs["section"] == sec).values
    X = adata.X[m][:, idx]
    if issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)
    X = X / (X.sum(1, keepdims=True) + 1e-12)
    P = np.zeros((X.shape[0], len(types)))
    for i in range(X.shape[0]):
        w, _ = nnls(S, X[i])
        s = w.sum()
        P[i] = w / s if s > 0 else 0
    props_all.append(P)
    print(sec, "done")

P_all = np.vstack(props_all)
props = pd.DataFrame(P_all, columns=types, index=adata.obs.index)
props.to_csv(os.path.join(INT, "M1_visium_nnls_props_cosmx.csv"))
summary = props.copy()
summary["condition"] = adata.obs["condition"].values
summ = summary.groupby("condition")[types].mean().T
summ.to_csv(os.path.join(OUT, "Table_S55u_M1_Visium_NNLS_CosMx_by_Condition.csv"))
print(summ.round(4).to_string())
print("step 4c done")
