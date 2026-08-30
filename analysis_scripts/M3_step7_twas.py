#!/usr/bin/env python3
"""
M3 腿3 步骤7：TWAS（FUSION 口径）GTEx v8 肺+全血权重 × FinnGen ARDS
====================================================================
权重：M3_step6 导出的 80 基因×2 组织（top1/enet/susie/lasso，cv rsq）
主口径：best-CV（cv rsq 最大且 >0 的模型，FUSION 默认）；敏感性：top1。
统计量：Z_twas = w'z / sqrt(w' Σ w)
  z 为 FinnGen z-score（对齐权重效应等位 a1）；
  Σ 为 1000G EUR 参考基因型（按权重 a1 编码）的 Pearson 相关矩阵。
输出：Table_S15g_TWAS_Results.csv；score_version = m3twas_v1.0
"""
import gzip, sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = Path(r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS")
RAW = BASE / "00_RAW_DATA"
INT = BASE / "_intermediate"
OUT = BASE / "02_SUPPLEMENTARY_TABLES" / "SUPPLEMENTARY_Tables_CSV"
FINNGEN = RAW / "GWAS" / "finngen_R10_J10_ARDS.gz"
EUR = RAW / "1000Genomes_Reference" / "EUR"
COMP = str.maketrans("ACGT", "TGCA")
GENE_SET_VERSION = "Mitoxy-80_v1.0"
SCORE_VERSION = "m3twas_v1.0"


def load_weights(tissue):
    f = INT / f"M3_twas_weights_{tissue}.csv"
    if not f.exists():
        return None
    return pd.read_csv(f, dtype={"rsid": str})


def choose_best(wdf):
    """best-CV 主口径：在 enet/lasso/susie 中取 cv_rsq 最大且 >0 的模型
    （GTEx v8 权重 'top1' 列为全 SNP 边际效应向量，非单 SNP 模型，不参与
    best-CV 选择；top1 敏感性改为单 SNP 测试）。"""
    per_gene = {}
    for gene, sub in wdf.groupby("gene"):
        cand = sub[sub["model"].isin(["enet", "lasso", "susie"])]
        cand = cand.dropna(subset=["cv_rsq"])
        cand = cand[cand["cv_rsq"] > 0]
        if len(cand) == 0:
            continue
        best_model = cand.loc[cand["cv_rsq"].idxmax(), "model"]
        per_gene[gene] = sub[sub["model"] == best_model]
    return per_gene


def choose_top1(wdf):
    """top1 敏感性：每基因取 |top1 边际权重| 最大的单个 SNP。"""
    out = {}
    for g, sub in wdf.groupby("gene"):
        t = sub[sub["model"] == "top1"]
        if len(t) == 0:
            continue
        top = t.loc[t["weight"].abs().idxmax()]
        out[g] = pd.DataFrame([top])
    return out


def load_1000g(needed_rsids):
    """读 EUR bed 中目标变体块（int8 剂量，2=两份 a1 等位）。
    返回 G[局部行], row_map(bim行号->局部行), rs2row, bim。"""
    bim = []
    with open(str(EUR) + ".bim") as fh:
        for line in fh:
            p = line.split()
            bim.append((p[1], p[4].upper(), p[5].upper()))  # rsid, a1, a2
    rs2row = {b[0]: i for i, b in enumerate(bim)}
    wanted = sorted({rs2row[rs] for rs in needed_rsids if rs in rs2row})
    row_map = {r: i for i, r in enumerate(wanted)}
    n_samples = 503
    n_bytes = (n_samples + 3) // 4
    G = np.zeros((len(wanted), n_samples), dtype=np.int8)
    with open(str(EUR) + ".bed", "rb") as fh:
        for r in wanted:
            fh.seek(3 + r * n_bytes)
            block = fh.read(n_bytes)
            if len(block) < n_bytes:
                continue
            buf = np.frombuffer(block, dtype=np.uint8)
            li = row_map[r]
            for j in range(n_samples):
                code = (buf[j // 4] >> (6 - 2 * (j % 4))) & 0b11
                if code == 0b00:
                    G[li, j] = 0
                elif code == 0b01:
                    G[li, j] = 1
                elif code == 0b10:
                    G[li, j] = 2
                else:
                    G[li, j] = -1  # missing
    return G, row_map, rs2row, bim


def extract_finngen_keys(keys):
    """keys: set of 'chr:pos'；返回 dict key->row。"""
    hits = {}
    with gzip.open(FINNGEN, "rt") as fh:
        for line in fh:
            if line.startswith("#chrom"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            k = f"{p[0]}:{p[1]}"
            if k in keys:
                hits[k] = p
    return hits


def run_twas(tissue, gene_models, G, row_map, rs2row, bim, fg_hits):
    rows = []
    detail = []
    for gene, wdf in sorted(gene_models.items()):
        model = wdf["model"].iloc[0]
        # 全量对齐到权重效应等位 a1：w 用原始权重；z 对齐 a1；基因型按 a1 编码
        w, z, gvecs, snp_info = [], [], [], []
        n_matched = 0
        for _, r in wdf.iterrows():
            key = f"{r['chr']}:{r['pos_hg38']}"
            p = fg_hits.get(key)
            if p is None:
                continue
            n_matched += 1
            ref, alt = p[2].upper(), p[3].upper()
            a1 = str(r["a1"]).upper()
            z_alt = float(p[8]) / float(p[9])
            s = 0
            if a1 == alt:
                s = 1
            elif a1 == ref:
                s = -1
            elif a1 == alt.translate(COMP):
                s = 1
            elif a1 == ref.translate(COMP):
                s = -1
            else:
                continue
            row = rs2row.get(r["rsid"])
            if row is None or row not in row_map:
                continue
            b_a1 = bim[row][1]
            b_a2 = bim[row][2]
            if a1 == b_a1:
                flip = False
            elif a1 == b_a2:
                flip = True
            else:
                continue
            gv = G[row_map[row]].astype(float)
            if flip:
                gv = 2 - gv
            miss = gv < 0
            if miss.all():
                continue
            gv[miss] = gv[~miss].mean() if (~miss).any() else 0.0
            w.append(float(r["weight"]))
            z.append(s * z_alt)
            gvecs.append(gv)
            snp_info.append(r["rsid"])
        if len(w) == 0:
            continue
        w = np.array(w)
        z = np.array(z)
        Gm = np.stack(gvecs)
        Gm = Gm - Gm.mean(axis=1, keepdims=True)
        Sig = (Gm @ Gm.T) / (Gm.shape[1] - 1)
        std = np.sqrt(np.diag(Sig))
        std[std == 0] = 1.0
        Sig = Sig / np.outer(std, std)
        num = float(w @ z)
        den = float(np.sqrt(w @ Sig @ w))
        twas_z = num / den if den > 0 else np.nan
        twas_p = 2 * (1 - stats.norm.cdf(abs(twas_z))) if np.isfinite(twas_z) else np.nan
        rows.append({
            "tissue": tissue, "gene": gene, "model": model,
            "n_snps_weight": len(wdf), "n_snps_matched_gwas": n_matched,
            "n_snps_ld": len(w), "twas_z": twas_z, "twas_p": twas_p,
        })
        for i in range(len(w)):
            detail.append({
                "tissue": tissue, "gene": gene, "model": model,
                "rsid": snp_info[i], "weight": w[i],
                "gwas_z_a1": z[i],
            })
    return pd.DataFrame(rows), pd.DataFrame(detail)


def main():
    print("=" * 70)
    print("M3 step7: TWAS (m3twas_v1.0)")
    print("=" * 70)
    # 先收集全部所需 rsid
    needed = set()
    tissue_weights = {}
    for tissue in ("GTExv8.EUR.Lung", "GTExv8.EUR.Whole_Blood"):
        wdf = load_weights(tissue)
        tissue_weights[tissue] = wdf
        if wdf is not None:
            needed.update(wdf["rsid"].dropna().astype(str))
    print(f"needed rsids: {len(needed)}")
    G, row_map, rs2row, bim = load_1000g(needed)
    print(f"1000G EUR: {len(row_map)}/{len(needed)} variants loaded")

    all_res, all_det = [], []
    for tissue in ("GTExv8.EUR.Lung", "GTExv8.EUR.Whole_Blood"):
        wdf = tissue_weights[tissue]
        if wdf is None:
            print(f"{tissue}: no weights file, skip")
            continue
        print(f"[{tissue}] one FinnGen pass for all model SNPs ...")
        keys = set()
        for _, r in wdf.iterrows():
            keys.add(f"{r['chr']}:{r['pos_hg38']}")
        print(f"  {len(keys)} SNP keys")
        fg_hits = extract_finngen_keys(keys)
        for tag, chooser in (("best", choose_best), ("top1", choose_top1)):
            gene_models = chooser(wdf)
            print(f"  [{tissue}/{tag}] {len(gene_models)} genes")
            res, det = run_twas(tissue, gene_models, G, row_map, rs2row, bim, fg_hits)
            res["method"] = tag
            det["method"] = tag
            all_res.append(res)
            all_det.append(det)
            print(f"  [{tissue}/{tag}] genes tested: {len(res)}")

    res = pd.concat(all_res, ignore_index=True)
    det = pd.concat(all_det, ignore_index=True)
    res["bh_p_pooled"] = multipletests(res["twas_p"], method="fdr_bh")[1]
    for tissue, sub in res.groupby("tissue"):
        res.loc[sub.index, "bh_p_tissue"] = multipletests(sub["twas_p"], method="fdr_bh")[1]
    res["gene_set_version"] = GENE_SET_VERSION
    res["score_version"] = SCORE_VERSION
    det["gene_set_version"] = GENE_SET_VERSION
    det["score_version"] = SCORE_VERSION
    res.to_csv(OUT / "Table_S15g_TWAS_Results.csv", index=False)
    det.to_csv(INT / "M3_twas_snp_detail.csv", index=False)
    print("\n===== TWAS (best-CV) TOP =====")
    best = res[res["method"] == "best"].sort_values("twas_p")
    print(best[["tissue", "gene", "model", "n_snps_ld", "twas_z", "twas_p", "bh_p_pooled"]].head(15).to_string(index=False))
    print(f"\nTWAS best-CV: genes tested={len(best)}, min p={best['twas_p'].min():.3g}, "
          f"BH<0.05 pooled={(best['bh_p_pooled'] < 0.05).sum()}")


if __name__ == "__main__":
    main()
