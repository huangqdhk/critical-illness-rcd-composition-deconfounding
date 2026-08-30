#!/usr/bin/env python3
"""
M3 腿4 步骤8：SeismicGWAS 空间链升级（位点→空间域→细胞类型→基因）
====================================================================
输入：Table_S19f_GWAS_Loci.csv（提示性位点 p<1e-5，FinnGen ARDS）
      M1_visium_scored.h5ad（23 切片 93,869 spots：X 表达 + domain（Banksy-lite k=5））
      M1_cosmx_type_means.csv（CosMx 基因×27 类型均值）
      M1_scRNA_type_means.csv（scRNA 8 类型均值，CosMx 缺席基因回退）
方法（预注册 m3seis_v1.0）：
  Visium 域归因：逐切片逐域基因均值 → 域内 z（切片内域间标准化）→ 跨切片 Stouffer z；
    逐切片域间 Kruskal-Wallis p → Fisher 合并。
  CosMx 细胞类型归因：基因×27 类型均值最大者为主类型（面板缺席回退 scRNA 8 类型）。
输出：Table_S19i_SeismicGWAS_Spatial_CellType_Chain.csv
      （另含 80 基因框架全域归因表，供"位点→域"联动核验）
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = Path(r"E:\SCI\proj")
INT = BASE / "_intermediate"
OUT = BASE / "02_SUPPLEMENTARY_TABLES" / "SUPPLEMENTARY_Tables_CSV"
GENE_SET_VERSION = "Mitoxy-80_v1.0"
SCORE_VERSION = "m3seis_v1.0"


def load_loci():
    f = OUT / "Table_S19f_GWAS_Loci.csv"
    if not f.exists():
        f = OUT / "Table_S19_GWAS_Mitoxy_Integration.csv"
    loci = pd.read_csv(OUT / "Table_S19f_GWAS_Loci.csv")
    return loci


def stouffer(ps):
    ps = np.array(ps, dtype=float)
    ps = ps[(ps > 0) & (ps < 1)]
    if len(ps) == 0:
        return np.nan, np.nan
    zs = -stats.norm.ppf(ps)
    z = zs.sum() / np.sqrt(len(ps))
    return z, 2 * (1 - stats.norm.cdf(abs(z)))


def visium_attribution(a, genes):
    """返回 DataFrame：gene -> top_domain, stouffer_z, kw_fisher_p, n_sections"""
    rows = []
    X = a.X
    obs = a.obs.reset_index(drop=True)
    gene_idx = {g: i for i, g in enumerate(a.var_names)}
    for gene in genes:
        if gene not in gene_idx:
            continue
        col = X[:, gene_idx[gene]].toarray().ravel() if hasattr(X, "toarray") else X[:, gene_idx[gene]]
        sub = obs.copy()
        sub["expr"] = col
        # 逐切片：域均值 + 域间 KW
        zs_by_domain, kw_ps = {}, []
        for sec, sdf in sub.groupby("section"):
            if sdf["domain"].nunique() < 2:
                continue
            dm = sdf.groupby("domain")["expr"].mean()
            if dm.std() == 0:
                zdom = pd.Series(0.0, index=dm.index)
            else:
                zdom = (dm - dm.mean()) / dm.std()
            for d, zv in zdom.items():
                zs_by_domain.setdefault(d, []).append(zv)
            groups = [sdf[sdf["domain"] == d]["expr"].values for d in sorted(sdf["domain"].unique())]
            if all(len(g) > 1 for g in groups) and len(groups) > 1:
                kw_ps.append(stats.kruskal(*groups).pvalue)
        if not zs_by_domain:
            continue
        domain_z = {d: np.mean(v) * np.sqrt(len(v)) for d, v in zs_by_domain.items()}
        top_domain = max(domain_z, key=domain_z.get)
        fisher_chi = -2 * np.sum(np.log([p for p in kw_ps if p > 0])) if kw_ps else np.nan
        kw_fisher_p = float(stats.chi2.sf(fisher_chi, 2 * len(kw_ps))) if kw_ps else np.nan
        rows.append({
            "gene": gene, "visium_top_domain": int(top_domain),
            "visium_stouffer_z": round(domain_z[top_domain], 4),
            "visium_kw_fisher_p": round(kw_fisher_p, 6),
            "n_slices": len(zs_by_domain.get(top_domain, [])),
        })
    return pd.DataFrame(rows)


def celltype_attribution(genes):
    cos = pd.read_csv(INT / "M1_cosmx_type_means.csv", index_col=0)
    scr = pd.read_csv(INT / "M1_scRNA_type_means.csv", index_col=0)
    rows = []
    for gene in genes:
        rec = {"gene": gene}
        if gene in cos.index:
            r = cos.loc[gene]
            rec["cosmx_top_celltype"] = r.idxmax()
            rec["cosmx_top_mean"] = round(float(r.max()), 4)
            rec["cosmx_mean_min"] = round(float(r.min()), 4)
            rec["attribution_source"] = "CosMx_27types"
        elif gene in scr.index:
            r = scr.loc[gene]
            rec["cosmx_top_celltype"] = r.idxmax()
            rec["cosmx_top_mean"] = round(float(r.max()), 4)
            rec["cosmx_mean_min"] = round(float(r.min()), 4)
            rec["attribution_source"] = "scRNA_8types_fallback"
        else:
            rec["cosmx_top_celltype"] = np.nan
            rec["cosmx_top_mean"] = np.nan
            rec["cosmx_mean_min"] = np.nan
            rec["attribution_source"] = "not_measured"
        rows.append(rec)
    return pd.DataFrame(rows)


def main():
    print("=" * 70)
    print("M3 step8: SeismicGWAS spatial chain (m3seis_v1.0)")
    print("=" * 70)
    a = ad.read_h5ad(INT / "M1_visium_scored.h5ad")
    loci = load_loci()
    print(f"loci: {len(loci)}")

    manifest = pd.read_csv(BASE / "04_AUDIT_GOVERNANCE" / "Mitoxyperilysis_Gene_Manifest_v1.0.csv")
    fw80 = set(manifest["gene_symbol"])

    # 位点内基因集合
    locus_genes = {}
    for _, r in loci.iterrows():
        gs = [g.strip() for g in str(r["genes"]).split(",") if g.strip()]
        locus_genes[r["locus_id"]] = gs
    all_genes = sorted(set(g for gs in locus_genes.values() for g in gs))

    vis = visium_attribution(a, all_genes)
    ct = celltype_attribution(all_genes)
    print(f"Visium attribution: {len(vis)} genes; celltype attribution: {len(ct)} genes")

    # 80 基因框架全域归因（联动核验用）
    vis80 = visium_attribution(a, sorted(fw80))
    ct80 = celltype_attribution(sorted(fw80))
    fw = vis80.merge(ct80, on="gene", how="outer")
    fw["in_80_framework"] = "Yes"
    fw["locus_id"] = ""
    fw["lead_rsid"] = ""

    chain = []
    for _, r in loci.iterrows():
        for g in locus_genes[r["locus_id"]]:
            v = vis[vis["gene"] == g]
            c = ct[ct["gene"] == g]
            chain.append({
                "locus_id": r["locus_id"], "chrom": r["chrom"],
                "lead_rsid": r["lead_rsid"], "lead_pos": r["lead_pos"],
                "lead_pval": r["lead_pval"], "gene": g,
                "in_80_framework": "Yes" if g in fw80 else "No",
                "visium_top_domain": v["visium_top_domain"].iloc[0] if len(v) else np.nan,
                "visium_stouffer_z": v["visium_stouffer_z"].iloc[0] if len(v) else np.nan,
                "visium_kw_fisher_p": v["visium_kw_fisher_p"].iloc[0] if len(v) else np.nan,
                "cosmx_top_celltype": c["cosmx_top_celltype"].iloc[0] if len(c) else np.nan,
                "cosmx_top_mean": c["cosmx_top_mean"].iloc[0] if len(c) else np.nan,
                "attribution_source": c["attribution_source"].iloc[0] if len(c) else "not_measured",
            })
    chain_df = pd.DataFrame(chain)
    fw_out = pd.concat([fw, chain_df[chain_df["in_80_framework"] == "No"]], ignore_index=True)
    chain_df["gene_set_version"] = GENE_SET_VERSION
    chain_df["score_version"] = SCORE_VERSION
    fw_out["gene_set_version"] = GENE_SET_VERSION
    fw_out["score_version"] = SCORE_VERSION
    chain_df.to_csv(OUT / "Table_S19i_SeismicGWAS_Spatial_CellType_Chain.csv", index=False)
    fw_out.to_csv(INT / "M3_seismic_80gene_spatial_attribution.csv", index=False)
    print("\n===== Locus -> domain -> celltype -> gene chain =====")
    print(chain_df[["locus_id", "lead_rsid", "gene", "in_80_framework",
                    "visium_top_domain", "visium_stouffer_z", "cosmx_top_celltype"]].to_string(index=False))
    print(f"\nSaved {len(chain_df)} chain rows + {len(fw_out)} 80-gene attribution rows")


if __name__ == "__main__":
    main()
