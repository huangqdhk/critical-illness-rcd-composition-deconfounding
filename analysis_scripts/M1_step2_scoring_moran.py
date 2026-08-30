# -*- coding: utf-8 -*-
"""
M1 Step 2: Two-arm spatial scoring + spatial autocorrelation + spatial domains
for GSE271370 Visium (23 sections).

Pre-registered (M1_pre_registration_20260817.md):
  score_version = "score_genes_v1_spatial" (AddModuleScore-like, 25 bins, seed=0)
  S1: per-arm global Moran's I (kNN k=8, row-standardized, 999 perms)
  S2: bivariate Moran's I (upstream x execution)  <- main compartmentalization test
  S3: spatial domains (HVG2000 -> PCA30 + neighbor-avg PCA30 -> KMeans k=5),
      Kruskal-Wallis across domains + domain-level Spearman rho
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats
from scipy.sparse import issparse
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import libpysal
import esda

INT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate"
OUT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\02_SUPPLEMENTARY_TABLES\SUPPLEMENTARY_Tables_CSV"
SEED = 0
K_NEIGH = 8
N_PERMS = 999
NPCS = 30
NDOMAINS = 5
N_BINS = 25

rng = np.random.default_rng(SEED)

man = pd.read_csv(r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\04_AUDIT_GOVERNANCE\Mitoxyperilysis_Gene_Manifest_v1.0.csv", encoding="utf-8-sig")
man["sym"] = man["hgnc_symbol"].fillna(man["gene_symbol"])
up_genes = man.loc[man.arm == "upstream_collapse", "sym"].tolist()
ex_genes = man.loc[man.arm == "execution_induction", "sym"].tolist()

adata = ad.read_h5ad(os.path.join(INT, "M1_visium_merged.h5ad"))
adata.var_names_make_unique()
all_genes = list(adata.var_names)
up_det = [g for g in up_genes if g in all_genes]
ex_det = [g for g in ex_genes if g in all_genes]
print(f"upstream detectable: {len(up_det)}/30 | execution detectable: {len(ex_det)}/33")


def dense_for(genes, sub):
    """Return (n_spots, n_genes) float64 dense slice for given genes."""
    idx = [all_genes.index(g) for g in genes]
    X = sub[:, idx].X
    if issparse(X):
        X = X.toarray()
    return np.asarray(X, dtype=np.float64)


def add_module_score(X_pool, set_idx, pool_mean, pool_std, n_bins=N_BINS, seed=SEED):
    """AddModuleScore-like scoring.
    X_pool: (n_spots, n_pool_genes) log1p expression; set_idx: indices into pool columns.
    Returns score vector (n_spots,)."""
    z = (X_pool - pool_mean) / pool_std
    bins = np.quantile(pool_mean, np.linspace(0, 1, n_bins + 1))
    bin_ids = np.digitize(pool_mean, bins) - 1
    bin_ids = np.clip(bin_ids, 0, n_bins - 1)
    set_idx = np.asarray(set_idx)
    ctrl_idx = []
    for i in set_idx:
        pool_candidates = np.where(bin_ids == bin_ids[i])[0]
        pool_candidates = pool_candidates[pool_candidates != i]
        if len(pool_candidates) == 0:
            pool_candidates = np.arange(X_pool.shape[1])
        ctrl_idx.append(int(rng.choice(pool_candidates)))
    ctrl_idx = np.array(ctrl_idx)
    return z[:, set_idx].mean(1) - z[:, ctrl_idx].mean(1)


def make_weights(coords, k=K_NEIGH):
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    dist, ind = nn.kneighbors(coords)
    neighbors = {i: ind[i, 1:].tolist() for i in range(len(coords))}
    w = libpysal.weights.W(neighbors, silence_warnings=True)
    w.transform = "r"
    return w, ind[:, 1:]


def moran_row(y, w, perms=N_PERMS, seed=SEED):
    np.random.seed(seed)
    mi = esda.Moran(y, w, permutations=perms, two_tailed=True)
    return mi.I, mi.p_sim, mi.z_sim


def moran_bv_row(x, y, w, perms=N_PERMS, seed=SEED):
    """Bivariate Moran's I = (x'W y)/(x'x) with row-standardized W; manual
    permutation test (two-tailed) for full control of the null."""
    Wy = w.sparse.dot(y)
    I_obs = float((x @ Wy) / (x @ x))
    rng = np.random.default_rng(seed)
    sims = np.empty(perms)
    for i in range(perms):
        yp = y[rng.permutation(len(y))]
        sims[i] = float((x @ w.sparse.dot(yp)) / (x @ x))
    p_two = float(min(1.0, 2 * min(np.mean(sims >= I_obs), np.mean(sims <= I_obs))))
    z = (I_obs - sims.mean()) / (sims.std() + 1e-12)
    return I_obs, p_two, z


def hvg_dispersion(sub, top=2000):
    """Variance-based HVG on log1p data (Seurat v1 style, deterministic)."""
    X = sub.X
    if issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)
    mu = X.mean(0)
    var = X.var(0)
    # standardize by mean: dispersion = var / (mu + 1e-6)
    disp = var / (mu + 1e-6)
    order = np.argsort(-disp)
    return order[:top]


rows = []
dom_rows = []
obs_score = {}
adata.obs["up_score_z"] = np.nan
adata.obs["ex_score_z"] = np.nan
adata.obs["mdi_z"] = np.nan
adata.obs["domain"] = -1
adata.obs["up_score_raw"] = np.nan
adata.obs["ex_score_raw"] = np.nan

for sec in sorted(adata.obs["section"].unique()):
    m = adata.obs["section"] == sec
    sub = adata[m].copy()
    n = sub.shape[0]
    cond = sub.obs["condition"].iloc[0]

    # detection filter for background pool (>=5% spots)
    Xs = sub.X
    if issparse(Xs):
        det = np.asarray((Xs > 0).sum(0)).ravel() / n
    else:
        det = (Xs > 0).mean(0)
    pool = np.where(det >= 0.05)[0]
    X_pool = dense_for([all_genes[i] for i in pool], sub)
    pool_mean = X_pool.mean(0)
    pool_std = X_pool.std(0) + 1e-12

    up_idx = [i for i, g in enumerate(pool) if all_genes[g] in up_det]
    ex_idx = [i for i, g in enumerate(pool) if all_genes[g] in ex_det]
    if len(up_idx) == 0 or len(ex_idx) == 0:
        print("skip", sec, "missing arm genes")
        continue
    up_raw = add_module_score(X_pool, up_idx, pool_mean, pool_std)
    ex_raw = add_module_score(X_pool, ex_idx, pool_mean, pool_std)
    up_z = (up_raw - up_raw.mean()) / (up_raw.std() + 1e-12)
    ex_z = (ex_raw - ex_raw.mean()) / (ex_raw.std() + 1e-12)
    mdi = ex_z - up_z

    coords = adata.obsm["spatial"][m.values]
    w, nn_idx = make_weights(coords)

    I_up, p_up, z_up_mi = moran_row(up_z, w)
    I_ex, p_ex, z_ex_mi = moran_row(ex_z, w)
    I_bv, p_bv, z_bv = moran_bv_row(up_z, ex_z, w)

    # neighborhood cross-correlation (lag-1): corr(own arm, neighbor-mean of other arm)
    ex_nb = ex_z[nn_idx].mean(1)
    up_nb = up_z[nn_idx].mean(1)
    ncc_up_to_ex = np.corrcoef(up_z, ex_nb)[0, 1]
    ncc_ex_to_up = np.corrcoef(ex_z, up_nb)[0, 1]

    # per-spot Spearman
    rho_spot, p_spot = stats.spearmanr(up_z, ex_z)
    # residualized (log10 total_counts + n_genes)
    tc = np.log10(adata.obs["n_counts"].values[m.values].astype(float))
    ng = adata.obs["n_genes"].values[m.values].astype(float)
    A = np.column_stack([np.ones(n), tc, ng])
    beta_up = np.linalg.lstsq(A, up_z, rcond=None)[0]
    beta_ex = np.linalg.lstsq(A, ex_z, rcond=None)[0]
    up_res = up_z - A @ beta_up
    ex_res = ex_z - A @ beta_ex
    rho_res, p_res = stats.spearmanr(up_res, ex_res)

    # spatial domains (non-circular w.r.t. arm scores)
    hv = hvg_dispersion(sub, top=2000)
    X_hv = dense_for([all_genes[i] for i in hv], sub)
    X_hv = (X_hv - X_hv.mean(0)) / (X_hv.std(0) + 1e-12)
    pca = PCA(n_components=min(NPCS, n, X_hv.shape[1]), random_state=SEED)
    pc = pca.fit_transform(X_hv)
    pc_nb = pc[nn_idx].mean(1)
    feat = np.hstack([pc, pc_nb])
    km = KMeans(n_clusters=NDOMAINS, n_init=10, random_state=SEED).fit(feat)
    dom = km.labels_

    # domain-level stats
    dom_df = pd.DataFrame({"domain": dom, "up_z": up_z, "ex_z": ex_z})
    g = dom_df.groupby("domain").agg(up_mean=("up_z", "mean"), ex_mean=("ex_z", "mean"),
                                     n_spots=("up_z", "size")).reset_index()
    if len(g) > 1:
        h_up, p_kw_up = stats.kruskal(*[dom_df.up_z[dom_df.domain == d].values for d in sorted(set(dom))])
        h_ex, p_kw_ex = stats.kruskal(*[dom_df.ex_z[dom_df.domain == d].values for d in sorted(set(dom))])
        rho_dom, p_dom = stats.spearmanr(g["up_mean"], g["ex_mean"])
    else:
        p_kw_up = p_kw_ex = rho_dom = p_dom = np.nan

    rows.append(dict(
        section=sec, condition=cond, n_spots=n,
        n_up_genes=len(up_idx), n_ex_genes=len(ex_idx),
        I_up=I_up, p_up=p_up, z_up=z_up_mi,
        I_ex=I_ex, p_ex=p_ex, z_ex=z_ex_mi,
        I_bv=I_bv, p_bv=p_bv, z_bv=z_bv,
        ncc_up_to_ex=ncc_up_to_ex, ncc_ex_to_up=ncc_ex_to_up,
        rho_spot=rho_spot, p_spot=p_spot, rho_res=rho_res, p_res=p_res,
        p_kw_up=p_kw_up, p_kw_ex=p_kw_ex, rho_dom=rho_dom, p_dom=p_dom,
        gene_set_version="Mitoxy-80_v1.0", score_version="score_genes_v1_spatial",
    ))
    for d, r in g.iterrows():
        dom_rows.append(dict(section=sec, condition=cond, domain=int(r.domain),
                             n_spots=int(r.n_spots), up_mean=r.up_mean, ex_mean=r.ex_mean,
                             gene_set_version="Mitoxy-80_v1.0", score_version="score_genes_v1_spatial"))

    adata.obs.loc[m, "up_score_z"] = up_z
    adata.obs.loc[m, "ex_score_z"] = ex_z
    adata.obs.loc[m, "mdi_z"] = mdi
    adata.obs.loc[m, "domain"] = dom
    adata.obs.loc[m, "up_score_raw"] = up_raw
    adata.obs.loc[m, "ex_score_raw"] = ex_raw
    print(f"{sec:12s} {cond:18s} n={n:5d} I_up={I_up:+.3f}(p={p_up:.4f}) "
          f"I_ex={I_ex:+.3f}(p={p_ex:.4f}) I_bv={I_bv:+.4f}(p={p_bv:.4f}) "
          f"rho_spot={rho_spot:+.3f} rho_dom={rho_dom:+.3f}")

sec_df = pd.DataFrame(rows)
sec_df.to_csv(os.path.join(OUT, "Table_S55a_M1_Visium_Section_Spatial_Stats.csv"), index=False)
pd.DataFrame(dom_rows).to_csv(os.path.join(OUT, "Table_S55c_M1_Visium_Domain_Stats.csv"), index=False)
adata.write_h5ad(os.path.join(INT, "M1_visium_scored.h5ad"))
print("saved Table_S55a / Table_S55c / M1_visium_scored.h5ad")

# ---- cross-section summaries ----
sdf = sec_df.copy()
# Stouffer on z_bv (two-sided)
z_st = sdf["z_bv"].sum() / np.sqrt(len(sdf))
p_st = 2 * stats.norm.sf(abs(z_st))
n_neg = int((sdf["I_bv"] < 0).sum())
n_neg_sig = int(((sdf["I_bv"] < 0) & (sdf["p_bv"] < 0.05)).sum())
n_pos_sig = int(((sdf["I_bv"] > 0) & (sdf["p_bv"] < 0.05)).sum())
n_up_sig = int(((sdf["I_up"] > 0) & (sdf["p_up"] < 0.05)).sum())
n_ex_sig = int(((sdf["I_ex"] > 0) & (sdf["p_ex"] < 0.05)).sum())
print("\n===== SUMMARY =====")
print(f"Stouffer z_bv = {z_st:+.3f}, p = {p_st:.4g}")
print(f"sections with I_bv<0: {n_neg}/23; I_bv<0 & p<0.05: {n_neg_sig}/23; "
      f"I_bv>0 & p<0.05: {n_pos_sig}/23")
print(f"sections with I_up>0 & p<0.05: {n_up_sig}/23; I_ex>0 & p<0.05: {n_ex_sig}/23")

# condition-level: per-section means
grp = adata.obs.groupby(["section", "condition"]).agg(
    up_mean=("up_score_z", "mean"), ex_mean=("ex_score_z", "mean"),
    mdi_mean=("mdi_z", "mean")).reset_index()
sec_df2 = sec_df.merge(grp, on="section", how="left")
cond_sum = sec_df2.groupby("condition").agg(
    n_sections=("section", "size"),
    up_mean_of_section_means=("up_mean", "mean"),
    ex_mean_of_section_means=("ex_mean", "mean"),
    mdi_mean=("mdi_mean", "mean"),
    I_bv_mean=("I_bv", "mean"),
    I_up_mean=("I_up", "mean"),
    I_ex_mean=("I_ex", "mean"),
).reset_index()
if cond_sum["n_sections"].nunique() > 1:
    h1, p1 = stats.kruskal(*[sec_df2.up_mean[sec_df2.condition == c].values for c in sorted(set(sec_df2.condition))])
    h2, p2 = stats.kruskal(*[sec_df2.ex_mean[sec_df2.condition == c].values for c in sorted(set(sec_df2.condition))])
    h3, p3 = stats.kruskal(*[sec_df2.mdi_mean[sec_df2.condition == c].values for c in sorted(set(sec_df2.condition))])
    cond_sum["p_kw_up"] = p1; cond_sum["p_kw_ex"] = p2; cond_sum["p_kw_mdi"] = p3
    print("\ncondition KW (per-section means): up p=", p1, "ex p=", p2, "mdi p=", p3)
cond_sum.to_csv(os.path.join(OUT, "Table_S55b_M1_Visium_Condition_Summary.csv"), index=False)
print("saved Table_S55b")
