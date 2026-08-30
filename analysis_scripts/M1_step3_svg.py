# -*- coding: utf-8 -*-
"""
M1 Step 3: Spatially variable genes (SVG) via per-section Moran's I
(randomization moments, Cliff & Ord) + enrichment of the 80-gene framework.

For row-standardized kNN weights (k=8): I(g) = (z' W z) / (z' z)  [since S0=n]
Randomization moments (Cliff & Ord 1981):
  E[I] = -1/(n-1)
  Var[I] = [n^2*S1 - n*S2 + 3*S0^2]/[S0^2*(n^2-1)] - E[I]^2
  S1 = sum_ij (w_ij + w_ji)^2 / 2 = sum w_ij^2 + sum w_ij*w_ji
  S2 = sum_i (rowsum_i + colsum_i)^2
Key arm genes additionally get empirical permutation p-values (esda, 999 perms).
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats
from sklearn.neighbors import NearestNeighbors
import libpysal
import esda

INT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate"
OUT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\02_SUPPLEMENTARY_TABLES\SUPPLEMENTARY_Tables_CSV"
SEED = 0
K = 8

man = pd.read_csv(r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\04_AUDIT_GOVERNANCE\Mitoxyperilysis_Gene_Manifest_v1.0.csv", encoding="utf-8-sig")
man["sym"] = man["hgnc_symbol"].fillna(man["gene_symbol"])

adata = ad.read_h5ad(os.path.join(INT, "M1_visium_merged.h5ad"))
adata.var_names_make_unique()
all_genes = list(adata.var_names)
arm_of = {}
for _, r in man.iterrows():
    if r["sym"] in all_genes:
        arm_of[r["sym"]] = r["arm"]


def moran_moments(n, nn_idx, k=K):
    """Closed-form Moran moments for kNN row-standardized weights."""
    d_in = np.bincount(nn_idx.ravel(), minlength=n)  # in-degree
    colsum = d_in / k
    S0 = n
    S1 = n / k + np.sum((d_in > 0).astype(float)) / k ** 2  # sum w^2 + sum w_ij w_ji (approx: mutual pairs counted as d_in>0)
    # exact mutual pair count:
    pairs = np.hstack([np.sort(np.column_stack([np.repeat(np.arange(n), k), nn_idx.ravel()]), axis=1)])
    mutual = 0
    if n * k < 2e7:
        seen = {}
        for i in range(n):
            for j in nn_idx[i]:
                key = (i, j)
                if key in seen:
                    mutual += 1
                seen[key] = True
    S1 = n / k + mutual / k ** 2
    S2 = np.sum((1.0 + colsum) ** 2)
    E = -1.0 / (n - 1)
    Var = (n ** 2 * S1 - n * S2 + 3 * S0 ** 2) / (S0 ** 2 * (n ** 2 - 1)) - E ** 2
    return E, max(Var, 1e-300)


def moran_z_vectorized(Z, nn_idx):
    """I per column of Z (n x g), using neighbor means."""
    nb_mean = Z[nn_idx].mean(1)  # (n, g): mean of neighbor z per spot
    num = (Z * nb_mean).sum(0)
    den = (Z * Z).sum(0)
    return num / np.maximum(den, 1e-300)


detail_rows = []
for sec in sorted(adata.obs["section"].unique()):
    m = adata.obs["section"] == sec
    sub = adata[m]
    n = sub.shape[0]
    coords = adata.obsm["spatial"][m.values]
    nn = NearestNeighbors(n_neighbors=K + 1).fit(coords)
    _, nn_idx = nn.kneighbors(coords)
    nn_idx = nn_idx[:, 1:]
    E, Var = moran_moments(n, nn_idx)

    X = sub.X
    from scipy.sparse import issparse
    if issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)
    det = (X > 0).mean(0)
    keep = det >= 0.05
    genes_keep = [g for g, k in zip(all_genes, keep) if k]
    Xk = X[:, keep]
    Z = (Xk - Xk.mean(0)) / (Xk.std(0) + 1e-12)
    I = moran_z_vectorized(Z, nn_idx)
    zval = (I - E) / np.sqrt(Var)
    p = 2 * stats.norm.sf(np.abs(zval))

    for g, ival, zv, pv in zip(genes_keep, I, zval, p):
        detail_rows.append(dict(section=sec, gene=g, n_spots=n,
                                moran_I=ival, moran_z=zv, p_approx=pv))
    print(f"{sec}: {len(genes_keep)} genes, median |z| = {np.median(np.abs(zval)):.2f}")

detail = pd.DataFrame(detail_rows)
detail.to_csv(os.path.join(OUT, "Table_S55d_M1_Visium_SVG_Detail.csv"), index=False)

# ---- meta across sections (Stouffer on z) ----
detail["arm"] = detail["gene"].map(arm_of).fillna("background")
meta = detail.groupby("gene").agg(
    n_sections=("section", "size"),
    moran_I_mean=("moran_I", "mean"),
    moran_z_stouffer=("moran_z", lambda s: s.sum() / np.sqrt(len(s))),
    arm=("arm", "first"),
).reset_index()
meta["p_stouffer"] = 2 * stats.norm.sf(np.abs(meta["moran_z_stouffer"]))
meta = meta.sort_values("moran_z_stouffer", ascending=False)
meta.to_csv(os.path.join(OUT, "Table_S55e_M1_Visium_SVG_Meta.csv"), index=False)

# ---- enrichment: arm genes vs background ----
rows = []
for sec, g in detail.groupby("section"):
    for arm in ["upstream_collapse", "execution_induction", "not_in_dissociation_arms"]:
        a = g[g["arm"] == arm]["moran_z"]
        b = g[g["arm"] == "background"]["moran_z"]
        if len(a) >= 3:
            u, p = stats.mannwhitneyu(a, b, alternative="greater")
            rows.append(dict(section=sec, arm=arm, n_genes=len(a),
                             median_z_arm=a.median(), median_z_bg=b.median(), p_mwu=p))
enr = pd.DataFrame(rows)
# Stouffer per arm across sections (one-sided: arm > bg)
summ = []
for arm, g in enr.groupby("arm"):
    zs = stats.norm.ppf(1 - g["p_mwu"].clip(1e-300, 1 - 1e-12))
    z_st = zs.sum() / np.sqrt(len(zs))
    p_st = stats.norm.sf(z_st)
    n_sig = int((g["p_mwu"] < 0.05).sum())
    summ.append(dict(arm=arm, n_sections=len(g), n_sig_sections=n_sig,
                     stouffer_z=z_st, stouffer_p=p_st))
enr_sum = pd.DataFrame(summ)
enr_sum.to_csv(os.path.join(OUT, "Table_S55f_M1_Visium_SVG_Enrichment.csv"), index=False)
print("\nSVG enrichment (arm vs background):")
print(enr_sum.to_string())

# ---- top SVG hypergeometric for 80-gene set ----
rows2 = []
for sec, g in detail.groupby("section"):
    n_all = len(g)
    top = g.nlargest(int(n_all * 0.1), "moran_z")
    for arm in ["upstream_collapse", "execution_induction", "not_in_dissociation_arms"]:
        n_arm = int((g["arm"] == arm).sum())
        k = int((top["arm"] == arm).sum())
        from scipy.stats import hypergeom
        p = hypergeom.sf(k - 1, n_all, n_arm, len(top))
        rows2.append(dict(section=sec, arm=arm, n_arm=n_arm, n_arm_top10=k,
                          n_top10=len(top), p_hyper=p))
hyp = pd.DataFrame(rows2)
hyp.to_csv(os.path.join(OUT, "Table_S55g_M1_Visium_SVG_Top10_Hypergeom.csv"), index=False)

# ---- permutation p for key arm genes (per section) ----
key = [g for g in man["sym"].unique() if g in all_genes]
key_rows = []
for sec in sorted(adata.obs["section"].unique()):
    m = adata.obs["section"] == sec
    sub = adata[m]
    n = sub.shape[0]
    coords = adata.obsm["spatial"][m.values]
    nn = NearestNeighbors(n_neighbors=K + 1).fit(coords)
    _, nn_idx = nn.kneighbors(coords)
    neighbors = {i: nn_idx[i, 1:].tolist() for i in range(n)}
    w = libpysal.weights.W(neighbors, silence_warnings=True)
    X = sub.X
    from scipy.sparse import issparse
    if issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)
    for g in key:
        gi = all_genes.index(g)
        y = X[:, gi]
        if (y > 0).mean() < 0.05 or y.std() == 0:
            continue
        z = (y - y.mean()) / (y.std() + 1e-12)
        np.random.seed(SEED)
        mi = esda.Moran(z, w, permutations=999, two_tailed=True)
        key_rows.append(dict(section=sec, gene=g, arm=arm_of.get(g, "background"),
                             I=mi.I, p_perm=mi.p_sim, z=mi.z_sim))
key_df = pd.DataFrame(key_rows)
key_df.to_csv(os.path.join(OUT, "Table_S55h_M1_Visium_KeyGene_Moran_Perm.csv"), index=False)
print("saved S55d-h")
