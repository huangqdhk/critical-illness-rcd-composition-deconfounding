# -*- coding: utf-8 -*-
"""
M1 Step 4a: Extract 8 cell-type signatures from the package scRNA object
(combined_processed.h5ad, backed mode) for Visium deconvolution.

Reference: GSE145926+GSE158055 merged, cell_type = {Mono_c14, T_cell, B_cell,
NK, Epithelial, DC, Club, Macrophage} (manuscript's canonical 8 types).

Outputs:
  - M1_scRNA_type_means.csv  (mean log1pCPM per type per gene)
  - M1_scRNA_type_markers.csv (top-50 markers per type, ranked by z-mean difference)
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

INT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate"
os.makedirs(INT, exist_ok=True)

a = ad.read_h5ad(r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\00_RAW_DATA\_scrna_work\combined_processed.h5ad", backed="r")
cell_types = sorted(a.obs["cell_type"].unique())
print("cell types:", cell_types)
groups = {ct: np.where(a.obs["cell_type"].values == ct)[0] for ct in cell_types}
genes = list(a.var_names)

# accumulate sums per cell type (chunked over cells to bound memory)
n_genes = len(genes)
sums = {ct: np.zeros(n_genes, dtype=np.float64) for ct in cell_types}
ns = {ct: len(idx) for ct, idx in groups.items()}
CHUNK = 20000
n_cells = a.shape[0]
for start in range(0, n_cells, CHUNK):
    end = min(start + CHUNK, n_cells)
    X = a.X[start:end]
    if sparse.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)
    for ct in cell_types:
        idx = groups[ct]
        sel = idx[(idx >= start) & (idx < end)] - start
        if len(sel):
            sums[ct] += X[sel].sum(0)
    print(f"chunk {start}-{end} done")

a.file.close()

means = pd.DataFrame({ct: sums[ct] / max(ns[ct], 1) for ct in cell_types}, index=genes)
means.to_csv(os.path.join(INT, "M1_scRNA_type_means.csv"))

# z-score genes across types, rank top-50 markers per type
z = (means - means.mean(1).values[:, None]) / (means.std(1).values[:, None] + 1e-12)
markers = {}
for ct in cell_types:
    top = z[ct].sort_values(ascending=False).head(50).index.tolist()
    markers[ct] = top
pd.DataFrame(markers).to_csv(os.path.join(INT, "M1_scRNA_type_markers.csv"), index=False)
print("saved M1_scRNA_type_means.csv / M1_scRNA_type_markers.csv")
for ct in cell_types:
    print(ct, "top5:", markers[ct][:5])
