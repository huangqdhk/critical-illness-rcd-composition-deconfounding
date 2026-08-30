# -*- coding: utf-8 -*-
"""
M1 Step 2b: cross-section summaries from Table_S55a (no recomputation).
Stouffer meta of bivariate Moran's I; sign counts; condition-level Kruskal-Wallis.
Saves Table_S55b.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import anndata as ad

INT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate"
OUT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\02_SUPPLEMENTARY_TABLES\SUPPLEMENTARY_Tables_CSV"

sdf = pd.read_csv(os.path.join(OUT, "Table_S55a_M1_Visium_Section_Spatial_Stats.csv"))
assert len(sdf) == 23

z_st = sdf["z_bv"].sum() / np.sqrt(len(sdf))
p_st = 2 * stats.norm.sf(abs(z_st))
n_neg = int((sdf["I_bv"] < 0).sum())
n_neg_sig = int(((sdf["I_bv"] < 0) & (sdf["p_bv"] < 0.05)).sum())
n_pos_sig = int(((sdf["I_bv"] > 0) & (sdf["p_bv"] < 0.05)).sum())
n_up_sig = int(((sdf["I_up"] > 0) & (sdf["p_up"] < 0.05)).sum())
n_ex_sig = int(((sdf["I_ex"] > 0) & (sdf["p_ex"] < 0.05)).sum())
n_dom_neg = int((sdf["rho_dom"] < 0).sum())
n_dom_neg_sig = int(((sdf["rho_dom"] < 0) & (sdf["p_dom"] < 0.05)).sum())

print("===== SUMMARY =====")
print(f"Stouffer z_bv = {z_st:+.3f}, p = {p_st:.4g}")
print(f"I_bv<0: {n_neg}/23 | I_bv<0 & p<0.05: {n_neg_sig}/23 | I_bv>0 & p<0.05: {n_pos_sig}/23")
print(f"I_up>0 & p<0.05: {n_up_sig}/23 | I_ex>0 & p<0.05: {n_ex_sig}/23")
print(f"domain rho<0: {n_dom_neg}/23 | rho<0 & p<0.05: {n_dom_neg_sig}/23")
print("median I_bv by condition:")
print(sdf.groupby("condition")["I_bv"].median().to_string())
print("median I_up / I_ex by condition:")
print(sdf.groupby("condition")[["I_up", "I_ex"]].median().to_string())

adata = ad.read_h5ad(os.path.join(INT, "M1_visium_scored.h5ad"))
grp = adata.obs.groupby(["section", "condition"]).agg(
    up_mean=("up_score_z", "mean"), ex_mean=("ex_score_z", "mean"),
    mdi_mean=("mdi_z", "mean")).reset_index()
grp = grp.drop(columns=["condition"])
sec_df2 = sdf.merge(grp, on="section", how="left")
cond_sum = sec_df2.groupby("condition").agg(
    n_sections=("section", "size"),
    up_mean_of_section_means=("up_mean", "mean"),
    ex_mean_of_section_means=("ex_mean", "mean"),
    mdi_mean=("mdi_mean", "mean"),
    I_bv_mean=("I_bv", "mean"),
    I_up_mean=("I_up", "mean"),
    I_ex_mean=("I_ex", "mean"),
).reset_index()
conds = sorted(set(sec_df2.condition))
h1, p1 = stats.kruskal(*[sec_df2.up_mean[sec_df2.condition == c].values for c in conds])
h2, p2 = stats.kruskal(*[sec_df2.ex_mean[sec_df2.condition == c].values for c in conds])
h3, p3 = stats.kruskal(*[sec_df2.mdi_mean[sec_df2.condition == c].values for c in conds])
h4, p4 = stats.kruskal(*[sec_df2.I_bv[sec_df2.condition == c].values for c in conds])
cond_sum["p_kw_up"] = p1; cond_sum["p_kw_ex"] = p2
cond_sum["p_kw_mdi"] = p3; cond_sum["p_kw_Ibv"] = p4
print("\ncondition KW (per-section): up p=%.3f ex p=%.3f mdi p=%.3f I_bv p=%.3f" % (p1, p2, p3, p4))
cond_sum.to_csv(os.path.join(OUT, "Table_S55b_M1_Visium_Condition_Summary.csv"), index=False)

# pairwise post-hoc for mdi/up/ex if KW p<0.05
post = []
for var in ["up_mean", "ex_mean", "mdi_mean", "I_bv"]:
    if {"up_mean": p1, "ex_mean": p2, "mdi_mean": p3, "I_bv": p4}[var] < 0.05:
        from itertools import combinations
        for c1, c2 in combinations(conds, 2):
            a = sec_df2[var][sec_df2.condition == c1].values
            b = sec_df2[var][sec_df2.condition == c2].values
            u, p = stats.mannwhitneyu(a, b)
            post.append(dict(variable=var, group1=c1, group2=c2, p_mwu=p,
                             median1=np.median(a), median2=np.median(b)))
post = pd.DataFrame(post)
if len(post):
    # BH within variable
    post["p_bh"] = post.groupby("variable")["p_mwu"].transform(
        lambda s: s * len(s) / s.rank())
    post.to_csv(os.path.join(OUT, "Table_S55r_M1_Visium_Condition_Posthoc.csv"), index=False)
print("step 2b done")
