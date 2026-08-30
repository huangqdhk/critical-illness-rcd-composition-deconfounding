#!/usr/bin/env python3
"""
M3 腿2 步骤5：SMR / HEIDI（eQTLGen 全血 + GTEx v8 肺 eQTL × FinnGen ARDS）
=============================================================================
层1 eQTLGen（n=31,470，全血）：Z→b 转换 b=Z/sqrt(2p(1-p)(N+Z²))，p 取 1000G EUR
    AssessedAllele 频率（近似，如实披露）；SNP 位置 hg19 → FinnGen 按 rsID 匹配。
层2 GTEx v8 肺 signif pairs（q<0.05）：slope/slope_se 原生；hg38 → chr:pos 匹配。
SMR：top cis-eQTL 为工具；b_SMR=b_GWAS/b_eQTL，delta 法 SE。
HEIDI：d_i=b_GWAS,i−b_SMR·b_eQTL,i；var(d_i)=se²_GWAS,i+b_SMR²·se²_eQTL,i；
      T=Σd²/var ~ χ²(m−1)；m≥3 才执行，否则 NA；p<0.05 → 异质性（连锁/多效）。
判定（预注册）：SMR BH<0.05 且 HEIDI 通过 → 共享因果支持；其余不成立。
score_version = m3smr_v1.0
输出：Table_S15f_SMR_HEIDI_Results.csv / Table_S15f_SMR_HEIDI_GeneDetail.csv
"""
import gzip, math, subprocess, sys
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
EQTLGEN = RAW / "GWAS" / "eQTLGen" / "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
GTEX_LUNG = RAW / "GWAS" / "gtex_v8_eQTL" / "Lung.v8.signif_variant_gene_pairs.txt.gz"
PLINK = r"E:\SCI\plink19\plink.exe"
EUR = RAW / "1000Genomes_Reference" / "EUR"
COMP = str.maketrans("ACGT", "TGCA")
N_EQTLGEN = 31470
CIS_MB = 1_000_000
GENE_SET_VERSION = "Mitoxy-80_v1.0"
SCORE_VERSION = "m3smr_v1.0"


def load_manifest():
    m = pd.read_csv(BASE / "04_AUDIT_GOVERNANCE" / "Mitoxyperilysis_Gene_Manifest_v1.0.csv")
    return m[["gene_symbol", "ensembl_gene_id", "module", "arm"]]


def parse_eqtlgen(ensg_set):
    """流式解析 eQTLGen：80 基因 cis ±1Mb 行。"""
    rows = []
    with gzip.open(EQTLGEN, "rt") as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 15:
                continue
            if p[7] not in ensg_set:
                continue
            try:
                snpchr, snppos = int(p[2]), int(p[3])
                genechr, genepos = int(p[9]), int(p[10])
            except ValueError:
                continue
            if snpchr != genechr or abs(snppos - genepos) > CIS_MB:
                continue
            rows.append({
                "gene_ensg": p[7],
                "rsid": p[1], "chr": snpchr, "pos_hg19": snppos,
                "assessed": p[4], "other": p[5],
                "zscore": float(p[6]), "p_eqtl": float(p[0]),
                "gene_pos": genepos,
            })
    df = pd.DataFrame(rows)
    if len(df):
        df["gene_symbol"] = df["gene_ensg"].map(
            {r.ensembl_gene_id: r.gene_symbol for r in load_manifest().itertuples()})
    return df


def parse_gtex_lung_api(manifest, gene_pos):
    """GTEx Portal API v2：按位点查询 GTEx v8 肺显著 eQTL（80 基因 TSS±1Mb）。
    NES 与 slope 成逐基因常数比例（NES=slope/SD(Y)），SMR 方向与 HEIDI 均对该
    缩放不变（推导见 M3_README）；se' = |NES|/qnorm(p/2, lower.tail=F)。"""
    import json as _json
    rows = []
    sym2ensg = dict(zip(manifest["gene_symbol"], manifest["ensembl_gene_id"]))
    for gene in sorted(gene_pos):
        chrom, tss = gene_pos[gene]
        lo, hi = tss - 1_200_000, tss + 1_200_000
        page = 1
        while True:
            url = ("https://gtexportal.org/api/v2/association/singleTissueEqtlByLocation"
                   f"?tissueSiteDetailId=Lung&chromosome=chr{chrom}&start={lo}&end={hi}"
                   f"&datasetId=gtex_v8&itemsPerPage=250&page={page}")
            cmd = ["curl.exe", "-s", "--max-time", "120", url]
            r = subprocess.run(cmd, capture_output=True, text=True,
                               encoding="utf-8", errors="replace")
            if r.returncode != 0:
                break
            try:
                js = _json.loads(r.stdout)
            except _json.JSONDecodeError:
                break
            data = js.get("singleTissueEqtl") or []
            if not data:
                break
            for d in data:
                gs = d.get("geneSymbol")
                if gs not in sym2ensg:
                    continue
                nes = d.get("nes")
                pval = d.get("pValue")
                if nes is None or pval is None or pval <= 0 or pval >= 1:
                    continue
                z = stats.norm.ppf(pval / 2)
                rows.append({
                    "gene_symbol": gs, "gene_ensg": sym2ensg[gs],
                    "rsid": d.get("snpId"), "chr": int(chrom),
                    "pos_hg38": int(d.get("pos")),
                    "variant_id": d.get("variantId"),
                    "p_eqtl": float(pval), "b_eqtl": float(nes),
                    "se_eqtl": abs(float(nes)) / abs(z),
                })
            pg = js.get("paging_info") or {}
            npages = pg.get("numberOfPages", 0)
            if page >= npages:
                break
            page += 1
    return pd.DataFrame(rows)


def get_eqtlgen_af(rsids):
    """plink --freq 计算 1000G EUR 频率（rsid → (A1, MAF)）。"""
    if not rsids:
        return {}
    lst = INT / "M3_smr_af_snps.txt"
    with open(lst, "w") as f:
        for r in sorted(rsids):
            f.write(r + "\n")
    out = INT / "M3_smr_af"
    cmd = [PLINK, "--bfile", str(EUR), "--extract", str(lst), "--freq", "--out", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-300:])
    af = {}
    with open(str(out) + ".frq") as f:
        next(f)
        for line in f:
            p = line.split()
            if len(p) >= 5:
                af[p[1]] = (p[2], float(p[4]))
    return af


def extract_finngen(lookup):
    """单次流式提取 FinnGen 行。lookup: {'rsid:xxx' 或 'pos:chr:pos' -> key}"""
    want_rsid = {k.split(":", 1)[1] for k in lookup if k.startswith("rsid:")}
    want_pos = {k.split(":", 1)[1] for k in lookup if k.startswith("pos:")}
    hits = {}
    with gzip.open(FINNGEN, "rt") as fh:
        for line in fh:
            if line.startswith("#chrom"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            for tok in p[4].split(","):
                if tok in want_rsid:
                    hits[f"rsid:{tok}"] = p
            if f"{p[0]}:{p[1]}" in want_pos:
                hits[f"pos:{p[0]}:{p[1]}"] = p
    return hits


def harmonize_gwas(p, ea_gwas_info):
    """返回 (beta_gwas, se_gwas) 相对 eQTL 效应等位。"""
    # FinnGen: ref/alt，beta 相对 alt
    ref, alt = p[2].upper(), p[3].upper()
    beta, se = float(p[8]), float(p[9])
    ea, oa = ea_gwas_info
    ea, oa = ea.upper(), oa.upper()
    if ea == alt:
        return beta, se
    if ea == ref:
        return -beta, se
    if ea == alt.translate(COMP):
        return beta, se
    if ea == ref.translate(COMP):
        return -beta, se
    return None


def run_layer(name, pairs, lookup_key, af_map, fg_hits):
    """pairs: DataFrame 带基因、eQTL 统计、效应等位；逐基因 SMR+HEIDI。"""
    out = []
    detail = []
    for (gene, chrom), gdf in pairs.groupby(["gene_symbol", "chr"], sort=False):
        recs = []
        for _, r in gdf.iterrows():
            key = lookup_key(r)
            p = fg_hits.get(key)
            if p is None:
                continue
            ea_info = r["ea_info"]
            h = harmonize_gwas(p, ea_info)
            if h is None:
                continue
            beta_g, se_g = h
            recs.append({
                "rsid": r["rsid"], "pos_hg38": r.get("pos_hg38"),
                "p_eqtl": r["p_eqtl"], "b_eqtl": r["b_eqtl"], "se_eqtl": r["se_eqtl"],
                "b_gwas": beta_g, "se_gwas": se_g,
            })
        if not recs:
            continue
        rec = pd.DataFrame(recs).sort_values("p_eqtl")
        top = rec.iloc[0]
        b_smr = top["b_gwas"] / top["b_eqtl"]
        se_smr = math.sqrt(
            (top["se_gwas"] / top["b_eqtl"]) ** 2 +
            (top["b_gwas"] * top["se_eqtl"] / top["b_eqtl"] ** 2) ** 2)
        p_smr = 2 * (1 - stats.norm.cdf(abs(b_smr / se_smr)))
        m = len(rec)
        heidi_p = heidi_T = np.nan
        if m >= 3:
            d = rec["b_gwas"] - b_smr * rec["b_eqtl"]
            vard = rec["se_gwas"] ** 2 + (b_smr * rec["se_eqtl"]) ** 2
            heidi_T = float((d ** 2 / vard).sum())
            heidi_p = float(stats.chi2.sf(heidi_T, m - 1))
        for _, r in rec.iterrows():
            detail.append({
                "layer": name, "gene_symbol": gene, "chr": chrom,
                "rsid": r["rsid"], "pos_hg38": r["pos_hg38"],
                "p_eqtl": r["p_eqtl"], "b_eqtl": r["b_eqtl"], "se_eqtl": r["se_eqtl"],
                "b_gwas": r["b_gwas"], "se_gwas": r["se_gwas"],
            })
        out.append({
            "layer": name, "gene_symbol": gene, "chr": chrom,
            "top_snp": top["rsid"], "n_snps_matched": m,
            "top_p_eqtl": top["p_eqtl"], "top_b_eqtl": top["b_eqtl"],
            "top_se_eqtl": top["se_eqtl"],
            "top_b_gwas": top["b_gwas"], "top_se_gwas": top["se_gwas"],
            "b_smr": b_smr, "se_smr": se_smr, "p_smr": p_smr,
            "heidi_T": heidi_T, "heidi_df": m - 1 if m >= 3 else np.nan,
            "heidi_p": heidi_p,
            "heidi_pass": bool(m >= 3 and heidi_p >= 0.05),
        })
    return pd.DataFrame(out), pd.DataFrame(detail)


def get_gene_positions(manifest):
    """80 基因 hg38 基因座中点（来自导出的 FUSION 权重 SNP 坐标范围；
    权重围绕基因 ±500kb，中点误差 <20kb，查询窗取 ±1.2Mb 覆盖 TSS±1Mb）。"""
    sym2ensg = dict(zip(manifest["gene_symbol"], manifest["ensembl_gene_id"]))
    out = {}
    for tissue in ("GTExv8.EUR.Lung", "GTExv8.EUR.Whole_Blood"):
        f = INT / f"M3_twas_weights_{tissue}.csv"
        if not f.exists():
            continue
        w = pd.read_csv(f, dtype={"rsid": str})
        w["ensg"] = w["gene"]
        for ensg, sub in w.groupby("ensg"):
            if ensg in out:
                continue
            chrom = sub["chr"].iloc[0]
            lo, hi = sub["pos_hg38"].min(), sub["pos_hg38"].max()
            out[ensg] = (int(chrom), (int(lo) + int(hi)) // 2)
    res = {}
    for sym, ensg in sym2ensg.items():
        if ensg in out:
            res[sym] = out[ensg]
    return res


def main():
    print("=" * 70)
    print("M3 step5: SMR/HEIDI (m3smr_v1.0)")
    print("=" * 70)
    manifest = load_manifest()
    ensg_set = set(manifest["ensembl_gene_id"])
    symbols = set(manifest["gene_symbol"])

    # ---- 层1 eQTLGen ----
    print("[1/4] parsing eQTLGen (whole blood) ...")
    eq = parse_eqtlgen(ensg_set)
    print(f"  rows: {len(eq)}, genes: {eq['gene_symbol'].nunique() if len(eq) else 0}")
    eq.to_csv(INT / "M3_smr_eqtlgen_rows.csv", index=False)

    # ---- 层2 GTEx v8 肺（API） ----
    print("[2/4] querying GTEx v8 Lung eQTL via Portal API ...")
    gene_pos = get_gene_positions(manifest)
    print(f"  genes with hg38 TSS: {len(gene_pos)}/80")
    gl = parse_gtex_lung_api(manifest, gene_pos)
    print(f"  rows: {len(gl)}, genes: {gl['gene_symbol'].nunique() if len(gl) else 0}")
    gl.to_csv(INT / "M3_smr_gtex_lung_rows.csv", index=False)

    # ---- eQTLGen AF + Z→b ----
    af_map = get_eqtlgen_af(sorted(set(eq["rsid"]))) if len(eq) else {}
    print(f"[3/4] 1000G EUR AF lookup: {len(af_map)} rsids")
    if len(eq):
        conv = []
        for _, r in eq.iterrows():
            af = af_map.get(r["rsid"])
            if af is None:
                continue
            a1, maf = af
            if a1 not in (r["assessed"], r["other"]):
                continue
            z = r["zscore"]
            p = maf
            b = z / math.sqrt(2 * p * (1 - p) * (N_EQTLGEN + z ** 2))
            conv.append({**r, "b_eqtl": b, "se_eqtl": 1 / math.sqrt(
                2 * p * (1 - p) * (N_EQTLGEN + z ** 2))})
        eq = pd.DataFrame(conv)
        eq["ea_info"] = list(zip(eq["assessed"], eq["other"]))
        eq["pos_hg38"] = np.nan  # eQTLGen 无 hg38 坐标，按 rsid 匹配
        print(f"  eQTLGen converted rows: {len(eq)}")

    # GTEx 肺层：variant_id chr_pos_ref_alt_b38 -> (ea=alt, oa=ref)
    if len(gl):
        sp = gl["variant_id"].str.split("_")
        gl["ea_info"] = list(zip(sp.str[3], sp.str[2]))
        gl = gl.dropna(subset=["pos_hg38"])

    # ---- FinnGen 提取 ----
    lookup = {}
    for _, r in eq.iterrows():
        lookup[f"rsid:{r['rsid']}"] = True
    for _, r in gl.iterrows():
        lookup[f"pos:{r['chr']}:{r['pos_hg38']}"] = True
    print(f"[4/4] streaming FinnGen for {len(lookup)} keys ...")
    fg_hits = extract_finngen(lookup)
    print(f"  FinnGen hits: {len(fg_hits)}")

    res1, det1 = run_layer("eQTLGen_whole_blood", eq,
                           lambda r: f"rsid:{r['rsid']}", af_map, fg_hits)
    res2, det2 = run_layer("GTEx_v8_Lung", gl,
                           lambda r: f"pos:{r['chr']}:{r['pos_hg38']}", af_map, fg_hits)
    res = pd.concat([res1, res2], ignore_index=True)
    det = pd.concat([det1, det2], ignore_index=True)

    # BH（合并 + 分层双口径）
    res["bh_p_pooled"] = multipletests(res["p_smr"], method="fdr_bh")[1]
    for layer, sub in res.groupby("layer"):
        res.loc[sub.index, "bh_p_layer"] = multipletests(sub["p_smr"], method="fdr_bh")[1]
    res["gene_set_version"] = GENE_SET_VERSION
    res["score_version"] = SCORE_VERSION
    det["gene_set_version"] = GENE_SET_VERSION
    det["score_version"] = SCORE_VERSION
    res.to_csv(OUT / "Table_S15f_SMR_HEIDI_Results.csv", index=False)
    det.to_csv(OUT / "Table_S15f_SMR_HEIDI_GeneDetail.csv", index=False)
    print("\n===== SMR/HEIDI RESULTS =====")
    cols = ["layer", "gene_symbol", "top_snp", "n_snps_matched", "top_p_eqtl",
            "b_smr", "p_smr", "bh_p_pooled", "heidi_p", "heidi_pass"]
    print(res[cols].sort_values("p_smr").to_string(index=False))
    print(f"\nSMR BH<0.05 (pooled): {(res['bh_p_pooled'] < 0.05).sum()}")
    print(f"SMR BH<0.05 & HEIDI pass: {((res['bh_p_pooled'] < 0.05) & (res['heidi_pass'] == True)).sum()}")


if __name__ == "__main__":
    main()
