# -*- coding: utf-8 -*-
"""
M1 Step 5: TNFSF13B-TFRC axis, spatial neighborhood enrichment (Visium).

For each section: TNFSF13B+ / TFRC+ = top-quartile log1p(CPM) spots.
Observed: mean fraction of TFRC+ among kNN(k=8) neighbors of TNFSF13B+ spots.
Null: 999 label permutations (fixed seed).  Meta: Stouffer across sections.
Also: bivariate Moran (z-scores) + cell-type attribution via NNLS proportions
(needs M1_visium_nnls_props_scRNA.csv from step 4b).
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import issparse

INT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate"
OUT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\02_SUPPLEMENTARY_TABLES\SUPPLEMENTARY_Tables_CSV"
SEED = 0
K = 8
NPERM = 999

adata = ad.read_h5ad(os.path.join(INT, "M1_visium_merged.h5ad"))
adata.var_names_make_unique()
genes = list(adata.var_names)
i_tn = genes.index("TNFSF13B")
i_tf = genes.index("TFRC")

rows = []
for sec in sorted(adata.obs["section"].unique()):
    m = (adata.obs["section"] == sec).values
    X = adata.X[m][:, [i_tn, i_tf]]
    if issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)
    n = X.shape[0]
    coords = adata.obsm["spatial"][m]
    nn = NearestNeighbors(n_neighbors=K + 1).fit(coords)
    _, nn_idx = nn.kneighbors(coords)
    nn_idx = nn_idx[:, 1:]

    thr_tn = np.quantile(X[:, 0], 0.75)
    thr_tf = np.quantile(X[:, 1], 0.75)
    pos_tn = X[:, 0] >= thr_tn
    pos_tf = X[:, 1] >= thr_tf
    n_tn, n_tf = int(pos_tn.sum()), int(pos_tf.sum())
    if n_tn < 20 or n_tf < 20:
        rows.append(dict(section=sec, n_spots=n, n_TNFSF13B=n_tn, n_TFRC=n_tf,
                         obs=0.25, z=np.nan, p=np.nan, direction="TNFSF13B->TFRC"))
        continue

    # TNFSF13B+ -> TFRC+ neighbors
    obs = pos_tf[nn_idx[pos_tn]].mean()
    rng = np.random.default_rng(SEED)
    perms = np.empty(NPERM)
    for i in range(NPERM):
        perms[i] = pos_tf[rng.permutation(n)[nn_idx[pos_tn]]].mean()
    z = (obs - perms.mean()) / (perms.std() + 1e-12)
    p = float(min(1.0, 2 * min(np.mean(perms >= obs), np.mean(perms <= obs))))
    rows.append(dict(section=sec, n_spots=n, n_TNFSF13B=n_tn, n_TFRC=n_tf,
                     obs=obs, perm_mean=perms.mean(), z=z, p=p,
                     direction="TNFSF13B->TFRC"))

    # reverse
    obs2 = pos_tn[nn_idx[pos_tf]].mean()
    perms2 = np.empty(NPERM)
    for i in range(NPERM):
        perms2[i] = pos_tn[rng.permutation(n)[nn_idx[pos_tf]]].mean()
    z2 = (obs2 - perms2.mean()) / (perms2.std() + 1e-12)
    p2 = float(min(1.0, 2 * min(np.mean(perms2 >= obs2), np.mean(perms2 <= obs2))))
    rows.append(dict(section=sec, n_spots=n, n_TNFSF13B=n_tn, n_TFRC=n_tf,
                     obs=obs2, perm_mean=perms2.mean(), z=z2, p=p2,
                     direction="TFRC->TNFSF13B"))

    # bivariate Moran on z-scores
    zt = (X[:, 0] - X[:, 0].mean()) / (X[:, 0].std() + 1e-12)
    zf = (X[:, 1] - X[:, 1].mean()) / (X[:, 1].std() + 1e-12)
    nb_mean = zf[nn_idx].mean(1)
    I_bv = float((zt * nb_mean).sum() / (zt * zt).sum())
    rows.append(dict(section=sec, n_spots=n, n_TNFSF13B=n_tn, n_TFRC=n_tf,
                     obs=I_bv, perm_mean=np.nan, z=np.nan, p=np.nan,
                     direction="bivariate_Moran_z"))
    print(f"{sec}: TN->TF obs={obs:.3f} z={z:+.2f} p={p:.3f} | I_bv={I_bv:+.3f}")

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "Table_S55s_M1_Visium_TNFSF13B_TFRC_Neighborhood.csv"), index=False)

# meta per direction
for d in ["TNFSF13B->TFRC", "TFRC->TNFSF13B", "bivariate_Moran_z"]:
    g = df[(df["direction"] == d) & df["z"].notna()]
    if len(g) == 0:
        g = df[(df["direction"] == d)]
        zst = np.nan; pst = np.nan
    else:
        zst = g["z"].sum() / np.sqrt(len(g))
        pst = 2 * stats.norm.sf(abs(zst))
    n_enr = int((g["z"] > 0).sum()) if "z" in g else np.nan
    print(f"{d}: Stouffer z={zst:+.2f} p={pst:.3g} n_pos_z={n_enr}/{len(g)}")

# ---- cell-type attribution via NNLS props (if available) ----
props_path = os.path.join(INT, "M1_visium_nnls_props_scRNA.csv")
if os.path.exists(props_path):
    props = pd.read_csv(props_path, index_col=0).loc[adata.obs.index]
    X = adata.X[:, [i_tn, i_tf]]
    if issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)
    att = {}
    for gname, col in [("TNFSF13B", X[:, 0]), ("TFRC", X[:, 1])]:
        att[gname] = {}
        for t in props.columns:
            rho, p = stats.spearmanr(col, props[t].values)
            att[gname][t] = rho
    att_df = pd.DataFrame(att)
    att_df.to_csv(os.path.join(OUT, "Table_S55t_M1_Visium_TNFSF13B_TFRC_CellType_Attribution.csv"))
    print(att_df.round(3).to_string())
print("step 5 done")
