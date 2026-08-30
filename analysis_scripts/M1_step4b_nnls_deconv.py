# -*- coding: utf-8 -*-
"""
M1 Step 4b: NNLS deconvolution of Visium spots using the package scRNA
8 cell-type signatures (marker-based gene subset).

Signature source: M1_scRNA_type_means.csv / M1_scRNA_type_markers.csv
(step 4a, combined_processed.h5ad, cell_type column).
Gene subset for NNLS: union of top-50 markers per type + 80-gene manifest genes.
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy.optimize import nnls
from scipy import stats

INT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate"
OUT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\02_SUPPLEMENTARY_TABLES\SUPPLEMENTARY_Tables_CSV"

means = pd.read_csv(os.path.join(INT, "M1_scRNA_type_means.csv"), index_col=0)
markers = pd.read_csv(os.path.join(INT, "M1_scRNA_type_markers.csv"))
types = list(markers.columns)
man = pd.read_csv(r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\04_AUDIT_GOVERNANCE\Mitoxyperilysis_Gene_Manifest_v1.0.csv", encoding="utf-8-sig")
man["sym"] = man["hgnc_symbol"].fillna(man["gene_symbol"])

adata = ad.read_h5ad(os.path.join(INT, "M1_visium_merged.h5ad"))
adata.var_names_make_unique()
all_genes = list(adata.var_names)

genes_use = set()
for t in types:
    genes_use.update([g for g in markers[t].dropna() if g in all_genes][:50])
genes_use.update([g for g in man["sym"] if g in all_genes])
genes_use = sorted(genes_use & set(means.index))
print("NNLS gene subset:", len(genes_use))

S = means.loc[genes_use, types].values.astype(np.float64)
# scale signatures to Visium log1pCPM scale: per-gene quantile-matched
idx = [all_genes.index(g) for g in genes_use]

props_all = []
for sec in sorted(adata.obs["section"].unique()):
    m = (adata.obs["section"] == sec).values
    X = adata.X[m][:, idx]
    from scipy.sparse import issparse
    if issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)
    P = np.zeros((X.shape[0], len(types)))
    for i in range(X.shape[0]):
        w, _ = nnls(S, X[i])
        s = w.sum()
        P[i] = w / s if s > 0 else 0
    props_all.append(P)
    print(sec, "done", P.shape)

P_all = np.vstack(props_all)
props = pd.DataFrame(P_all, columns=types, index=adata.obs.index)
props.to_csv(os.path.join(INT, "M1_visium_nnls_props_scRNA.csv"))

# summary per condition (per-spot means)
summary = props.copy()
summary["condition"] = adata.obs["condition"].values
summ = summary.groupby("condition")[types].mean().T
summ.to_csv(os.path.join(OUT, "Table_S55p_M1_Visium_NNLS_CellType_by_Condition.csv"))
print(summ.round(4).to_string())

# marker-score validation: correlation between NNLS prop and marker mean-z per type
scores = {}
for t in types:
    mk = [g for g in markers[t].dropna() if g in all_genes][:50]
    mi = [all_genes.index(g) for g in mk]
    X = adata.X[:, mi]
    if issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)
    Z = (X - X.mean(0)) / (X.std(0) + 1e-12)
    scores[t] = Z.mean(1)
scores = pd.DataFrame(scores, index=adata.obs.index)
corrs = {}
for t in types:
    rho, p = stats.spearmanr(props[t], scores[t])
    corrs[t] = dict(rho=rho, p=p)
corr_df = pd.DataFrame(corrs).T.reset_index().rename(columns={"index": "cell_type"})
corr_df.to_csv(os.path.join(OUT, "Table_S55q_M1_Visium_NNLS_MarkerScore_Correlation.csv"), index=False)
print(corr_df.round(4).to_string())
print("step 4b done")
