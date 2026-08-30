#!/usr/bin/env python3
"""
M3 腿1 步骤1：构建 pQTL 工具变量（cis-only）
=============================================
数据源（预注册 M3_pre_registration_20260818.md）：
  - deCODE pQTL 全基因组汇总统计（ITPR1/SOD2/TFRC；空格分隔；
    ID=hg19 坐标、GENPOS=hg38 坐标；BETA/SE 相对 ALLELE1）
  - UKB-PPP (Sun2023) MOESM3 ST15：显著 pQTL（p<1.7e-5）
  - 1000G Phase3 EUR（PLINK bfile，hg19）用于 rsID 映射与 LD clumping

口径：cis = 变体位于基因 hg38 [start-1Mb, end+1Mb]；p<5e-8；
LD clump r2<0.001 kb=10000（PLINK 1.9）；F>10。
输出：_intermediate/M3_pqtl_instruments_deCODE.csv（候选+clump 归属）
      _intermediate/M3_pqtl_instruments_final.csv（clump 索引 + UKB-PPP 敏感性行）
      _intermediate/M3_ukbppp_cis_pqtl.csv（UKB-PPP 覆盖度核查表）
score_version = m3pqtl_v1.0
"""
import gzip, os, subprocess, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = Path(r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS")
RAW = BASE / "00_RAW_DATA"
INT = BASE / "_intermediate"
INT.mkdir(parents=True, exist_ok=True)
PLINK = r"E:\SCI\plink19\plink.exe"
EUR = RAW / "1000Genomes_Reference" / "EUR"
GENE_SET_VERSION = "Mitoxy-80_v1.0"
SCORE_VERSION = "m3pqtl_v1.0"

# 基因坐标（Ensembl GRCh38，基因体）
GENES = {
    "ITPR1": {"chr": 3, "start": 4493696, "end": 5026555},
    "SOD2":  {"chr": 6, "start": 159679069, "end": 159762529},
    "TFRC":  {"chr": 3, "start": 196049527, "end": 196082199},
}
DECODE_DIRS = {
    "ITPR1": RAW / "deCODE" / "ITPR1_Q14643_OID30216_v1_Cardiometabolic_II",
    "SOD2":  RAW / "deCODE" / "SOD2_P04179_OID21114_v1_Neurology",
    "TFRC":  RAW / "deCODE" / "TFRC_P02786_OID20369_v1_Cardiometabolic",
}

COMP = str.maketrans("ACGT", "TGCA")
PVAL_THRESH = 5e-8
WINDOW = 1_000_000


def load_bim(chroms):
    """Stream EUR.bim, keep chr in `chroms`: (chr,pos) -> dict {allele_pair: rsid}."""
    m = {}
    with open(str(EUR) + ".bim") as fh:
        for line in fh:
            p = line.split()
            if len(p) < 6:
                continue
            try:
                c = int(p[0])
            except ValueError:
                continue
            if c not in chroms:
                continue
            pos = int(p[3])
            a1, a2 = p[4].upper(), p[5].upper()
            key = (c, pos)
            pair = tuple(sorted((a1, a2)))
            d = m.setdefault(key, {})
            if pair in d:
                d[pair] = None  # 冲突：同 chr:pos 同等位对多 rsID -> 弃用
            else:
                d[pair] = p[1]
    return m


def lookup_rsid(bim_map, chrom, pos, a0, a1):
    """deCODE 变体（hg19 pos + A0/A1）-> 1000G rsID；允许链翻转。"""
    a0, a1 = a0.upper(), a1.upper()
    cand = {tuple(sorted((a0, a1))), tuple(sorted((a0.translate(COMP), a1.translate(COMP))))}
    d = bim_map.get((chrom, pos))
    if not d:
        return None
    hits = [rs for pair, rs in d.items() if pair in cand and rs is not None]
    return hits[0] if len(hits) == 1 else None


def parse_decode(gene):
    """解析 deCODE cis 区：返回候选变体表（p<5e-8）。"""
    info = GENES[gene]
    lo, hi = info["start"] - WINDOW, info["end"] + WINDOW
    rows = []
    for chrom_file in sorted(DECODE_DIRS[gene].glob("discovery_chr*.gz")):
        fchrom = chrom_file.name.split("_chr")[1].split("_")[0]
        if fchrom != str(info["chr"]):
            continue
        with gzip.open(chrom_file, "rt") as fh:
            for line in fh:
                if line.startswith("CHROM"):
                    continue
                p = line.split()
                if len(p) < 14:
                    continue
                try:
                    genpos = int(p[1])
                except ValueError:
                    continue
                if not (lo <= genpos <= hi):
                    continue
                try:
                    log10p = float(p[12])
                except ValueError:
                    continue
                if log10p < 7.30103:  # -log10(5e-8)
                    continue
                idf = p[2].split(":")
                if len(idf) < 3:
                    continue
                rows.append({
                    "gene": gene, "chr": int(p[0]), "pos_hg38": genpos,
                    "pos_hg19": int(idf[1]), "allele0": p[3], "allele1": p[4],
                    "a1freq": float(p[5]), "n": int(p[7]),
                    "beta": float(p[9]), "se": float(p[10]),
                    "pval": 10 ** (-log10p),
                })
    df = pd.DataFrame(rows)
    df["F"] = (df["beta"] / df["se"]) ** 2
    df = df[df["F"] > 10].reset_index(drop=True)
    return df


def clump_gene(gene, df, bim_map):
    """plink --clump，返回每个候选的 rsid 与 clump 归属。"""
    recs = []
    for _, r in df.iterrows():
        rs = lookup_rsid(bim_map, r["chr"], r["pos_hg19"], r["allele0"], r["allele1"])
        recs.append((rs, r))
    mapped = [(rs, r) for rs, r in recs if rs is not None]
    n_miss = len(recs) - len(mapped)
    print(f"  [{gene}] candidates={len(recs)} rsid-mapped={len(mapped)} unmatched={n_miss}")
    if not mapped:
        return df.assign(rsid="", clump_index=False, clump_rsid="")
    clump_in = INT / f"M3_clump_{gene}.txt"
    with open(clump_in, "w") as f:
        f.write("SNP P\n")
        for rs, r in mapped:
            f.write(f"{rs} {r['pval']:.6e}\n")
    out = INT / f"M3_clump_{gene}"
    cmd = [PLINK, "--bfile", str(EUR), "--clump", str(clump_in),
           "--clump-p1", "5e-8", "--clump-p2", "5e-8",
           "--clump-r2", "0.001", "--clump-kb", "10000",
           "--clump-snp-field", "SNP", "--clump-field", "P",
           "--out", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    if r.returncode != 0:
        print("  plink stderr:", r.stderr[-500:])
        raise RuntimeError(f"plink clump failed for {gene}")
    clumped_file = Path(str(out) + ".clumped")
    index_snps = set()
    with open(clumped_file) as f:
        for line in f:
            parts = line.split()
            if not parts or parts[0] == "CHR":
                continue
            index_snps.add(parts[2])
    df = df.copy()
    df["rsid"] = [rs for rs, _ in recs]
    df["clump_index"] = df["rsid"].isin(index_snps)
    return df


def ukbppp_cis():
    """UKB-PPP ST15 中 8 目标基因的 cis pQTL 覆盖度核查。"""
    import openpyxl
    targets = {"ITPR1", "VDAC1", "CANX", "SLC40A1", "GPX4", "SOD2", "TFRC", "FTH1"}
    wb = openpyxl.load_workbook(
        RAW / "Sun2023_Nature_pQTL" / "41586_2023_6592_MOESM3_ESM.xlsx", read_only=True)
    ws = wb["ST15"]
    hdr, rows = None, []
    for r in ws.iter_rows(values_only=True):
        if hdr is None:
            if r[0] is not None and "Variant ID" in str(r[0]):
                hdr = [str(x) if x is not None else "" for x in r]
                idx = {h: i for i, h in enumerate(hdr)}
            continue
        if len(r) <= idx["Assay Target"] or r[idx["Assay Target"]] not in targets:
            continue
        rows.append({
            "gene": r[idx["Assay Target"]],
            "variant_id": r[idx["Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)"]],
            "chrom": r[idx["CHROM"]], "pos_hg38": r[idx["GENPOS (hg38)"]],
            "rsid": r[idx["rsID"]],
            "a1freq": r[idx["A1FREQ (discovery)"]],
            "beta": r[idx["BETA (discovery, wrt. A1)"]],
            "se": r[idx["SE (discovery)"]],
            "log10p": r[idx["log10(p) (discovery)"]],
            "cis_trans": r[idx["cis/trans"]],
        })
    wb.close()
    df = pd.DataFrame(rows)
    df.to_csv(INT / "M3_ukbppp_cis_pqtl.csv", index=False)
    cis = df[df["cis_trans"] == "cis"].copy()
    cis["F"] = (cis["beta"] / cis["se"]) ** 2
    cis["pval"] = 10 ** (-cis["log10p"].astype(float))
    # variant_id 形如 "chr:pos:A0:A1:imp:v1" -> 拆出等位
    sp = cis["variant_id"].str.split(":", regex=False)
    cis["allele0"] = sp.str[2]
    cis["allele1"] = sp.str[3]
    print("  UKB-PPP 8-gene coverage: cis hits =",
          cis["gene"].tolist() if len(cis) else "NONE")
    print("  genes WITHOUT significant cis pQTL:",
          sorted(targets - set(cis["gene"])))
    return cis


def main():
    print("=" * 70)
    print("M3 step1: build cis-pQTL instruments (m3pqtl_v1.0)")
    print("=" * 70)
    bim_map = load_bim({3, 6})
    print(f"bim map: {sum(len(d) for d in bim_map.values())} variants (chr3/6)")
    parts = []
    for gene in GENES:
        df = parse_decode(gene)
        df = clump_gene(gene, df, bim_map)
        parts.append(df)
    decode_all = pd.concat(parts, ignore_index=True)
    decode_all["source"] = "deCODE"
    decode_all["gene_set_version"] = GENE_SET_VERSION
    decode_all["score_version"] = SCORE_VERSION
    decode_all.to_csv(INT / "M3_pqtl_instruments_deCODE.csv", index=False)

    ukb = ukbppp_cis()
    if len(ukb):
        ukb_out = ukb[["gene", "variant_id", "chrom", "pos_hg38", "rsid",
                       "allele0", "allele1", "a1freq", "beta", "se", "pval", "F"]].copy()
        ukb_out = ukb_out.rename(columns={"chrom": "chr"})
        ukb_out["source"] = "UKB-PPP"
        ukb_out["clump_index"] = True  # 单变体敏感性
        ukb_out["gene_set_version"] = GENE_SET_VERSION
        ukb_out["score_version"] = SCORE_VERSION
        final = pd.concat([
            decode_all[decode_all["clump_index"]],
            ukb_out,
        ], ignore_index=True, sort=False)
    else:
        final = decode_all[decode_all["clump_index"]].copy()
    final.to_csv(INT / "M3_pqtl_instruments_final.csv", index=False)
    print("\nFinal instruments:")
    print(final[["gene", "rsid", "source", "beta", "se", "F"]].to_string(index=False))
    print(f"\nSaved -> _intermediate/M3_pqtl_instruments_final.csv ({len(final)} rows)")


if __name__ == "__main__":
    main()
