#!/usr/bin/env python3
"""
M3 附加核查：SLC7A11 全血 eQTL（GTEx v8 API）× FinnGen ARDS 共定位
用于甄别 TWAS 弥散模型信号（TWAS z=+5.13 但 top1 p=0.64）。
eQTL 侧：NES 缩放（beta=nes, se=|nes|/z）；sdY=1；MAF 取 1000G EUR。
"""
import gzip, json, math, subprocess, sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = Path(r"E:\SCI\proj")
INT = BASE / "_intermediate"
FINNGEN = BASE / "00_RAW_DATA" / "GWAS" / "finngen_R10_J10_ARDS.gz"
PLINK = r"E:\SCI\plink19\plink.exe"
EUR = BASE / "00_RAW_DATA" / "1000Genomes_Reference" / "EUR"
LO, HI = 137600000, 138900000  # SLC7A11 chr4 hg38 ±~600kb


def gtex_api(lo, hi):
    rows = []
    page = 1
    while True:
        url = ("https://gtexportal.org/api/v2/association/singleTissueEqtlByLocation"
               f"?tissueSiteDetailId=Whole_Blood&chromosome=chr4&start={lo}&end={hi}"
               f"&datasetId=gtex_v8&itemsPerPage=250&page={page}")
        r = subprocess.run(["curl.exe", "-s", "--max-time", "120", url],
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        if r.returncode != 0:
            break
        js = json.loads(r.stdout)
        data = js.get("singleTissueEqtl") or []
        if not data:
            break
        rows.extend(data)
        npg = js.get("paging_info", {}).get("numberOfPages", 0)
        if page >= npg:
            break
        page += 1
    return pd.DataFrame([d for d in rows if d.get("geneSymbol") == "SLC7A11"])


def finngen_region():
    rows = []
    with gzip.open(FINNGEN, "rt") as fh:
        for line in fh:
            if line.startswith("#chrom"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 12 or p[0] != "4":
                continue
            pos = int(p[1])
            if LO <= pos <= HI:
                rows.append(p)
    return rows


def af_lookup(rsids):
    if not rsids:
        return {}
    lst = INT / "M3_slc7a11_af_snps.txt"
    with open(lst, "w") as f:
        for r in sorted(rsids):
            f.write(r + "\n")
    out = INT / "M3_slc7a11_af"
    subprocess.run([PLINK, "--bfile", str(EUR), "--extract", str(lst),
                    "--freq", "--out", str(out)],
                   capture_output=True, text=True, encoding="utf-8", errors="replace")
    af = {}
    with open(str(out) + ".frq") as f:
        next(f)
        for line in f:
            p = line.split()
            if len(p) >= 5:
                af[p[1]] = (p[2], float(p[4]))
    return af


def main():
    eq = gtex_api(LO, HI)
    print("SLC7A11 whole-blood eQTL rows:", len(eq))
    fg = finngen_region()
    fg_map = {}
    for p in fg:
        fg_map[int(p[1])] = p
    print("FinnGen region rows:", len(fg_map))
    af = af_lookup(set(eq["snpId"].dropna()))
    print("AF lookup:", len(af))

    snp, beta_e, se_e, maf_e = [], [], [], []
    beta_o, se_o, maf_o = [], [], []
    used = 0
    for _, r in eq.iterrows():
        pos = int(r["pos"])
        p = fg_map.get(pos)
        if p is None:
            continue
        rs = r["snpId"]
        a = af.get(rs)
        if a is None:
            continue
        nes, pval = r["nes"], r["pValue"]
        if nes is None or pval is None or pval <= 0 or pval >= 1 or nes == 0:
            continue
        z = stats.norm.ppf(pval / 2)
        se_e_val = abs(nes) / abs(z)
        # FinnGen 侧按 alt 对齐 NES 方向：variantId chr4_pos_ref_alt_b38
        vid = r["variantId"].split("_")
        ref, alt = vid[2].upper(), vid[3].upper()
        fg_ref, fg_alt = p[2].upper(), p[3].upper()
        if alt == fg_alt:
            s = 1
        elif alt == fg_ref:
            s = -1
        elif alt == fg_alt.translate(str.maketrans("ACGT", "TGCA")):
            s = 1
        elif alt == fg_ref.translate(str.maketrans("ACGT", "TGCA")):
            s = -1
        else:
            continue
        try:
            b_o, se_o_val = float(p[8]), float(p[9])
        except ValueError:
            continue
        maf_o_val = min(float(p[10]), 1 - float(p[10]))
        snp.append(f"chr4:{pos}")
        beta_e.append(float(nes))
        se_e.append(se_e_val)
        maf_e.append(min(a[1], 1 - a[1]))
        beta_o.append(s * b_o)
        se_o.append(se_o_val)
        maf_o.append(maf_o_val)
        used += 1
    print("shared SNPs for coloc:", used)
    d = pd.DataFrame({
        "snp": snp, "beta_e": beta_e, "se_e": se_e, "maf_e": maf_e,
        "beta_o": beta_o, "se_o": se_o, "maf_o": maf_o,
    })
    d.to_csv(INT / "M3_coloc_SLC7A11_wholeblood.csv", index=False)
    print("saved -> _intermediate/M3_coloc_SLC7A11_wholeblood.csv")
    print(d.head(3).to_string())


if __name__ == "__main__":
    main()
