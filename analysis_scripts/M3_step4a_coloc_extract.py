#!/usr/bin/env python3
"""
M3 腿1 步骤4a：共定位输入提取（lead±100kb）
===========================================
对每个（基因×结局）组合提取区域全量 SNP：
  deCODE（hg38 GENPOS + hg19 ID 坐标，A1FREQ/N/BETA/SE）
  FinnGen ARDS（本地流式，hg38，ref/alt/beta/sebeta/af_alt）
  IEU 脓毒症（API 区间查询，hg19）
等位 harmonize（非回文：方向对齐；回文 MAF<=0.42 保留、>0.42 剔除）
输出：_intermediate/M3_coloc_{gene}_{outcome}.csv（snp 键 + 双侧统计）
"""
import gzip, json, subprocess, sys
from pathlib import Path

import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = Path(r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS")
RAW = BASE / "00_RAW_DATA"
INT = BASE / "_intermediate"
FINNGEN = RAW / "GWAS" / "finngen_R10_J10_ARDS.gz"
TOKEN_FILE = BASE / ".opengwas_token"
COMP = str.maketrans("ACGT", "TGCA")
WINDOW = 100_000
NCASE_FG, NCONTROL_FG = 357, 406536
NCASE_IEU, NCONTROL_IEU = 10154, 454764

DECODE_DIRS = {
    "ITPR1": RAW / "deCODE" / "ITPR1_Q14643_OID30216_v1_Cardiometabolic_II",
    "SOD2":  RAW / "deCODE" / "SOD2_P04179_OID21114_v1_Neurology",
    "TFRC":  RAW / "deCODE" / "TFRC_P02786_OID20369_v1_Cardiometabolic",
}
GENES = {"ITPR1": 3, "SOD2": 6, "TFRC": 3}


def get_leads():
    """deCODE 每基因 top 变体（LOG10P 最大者）。"""
    cand = pd.read_csv(INT / "M3_pqtl_instruments_deCODE.csv")
    leads = {}
    for g, sub in cand.groupby("gene"):
        top = sub.loc[sub["pval"].idxmin()]
        leads[g] = {"pos_hg38": int(top["pos_hg38"]), "pos_hg19": int(top["pos_hg19"]),
                    "rsid": top["rsid"]}
    return leads


def extract_decode(gene, lead):
    """deCODE lead±100kb 全量行（hg38 窗口过滤）。"""
    chrom = GENES[gene]
    lo, hi = lead["pos_hg38"] - WINDOW, lead["pos_hg38"] + WINDOW
    rows = []
    for f in sorted(DECODE_DIRS[gene].glob("discovery_chr*.gz")):
        if f.name.split("_chr")[1].split("_")[0] != str(chrom):
            continue
        with gzip.open(f, "rt") as fh:
            for line in fh:
                if line.startswith("CHROM"):
                    continue
                p = line.split()
                if len(p) < 14:
                    continue
                try:
                    gp = int(p[1])
                except ValueError:
                    continue
                if not (lo <= gp <= hi):
                    continue
                idf = p[2].split(":")
                if len(idf) < 3:
                    continue
                rows.append({
                    "pos_hg38": gp, "pos_hg19": int(idf[1]),
                    "a0": p[3], "a1": p[4],
                    "maf": min(float(p[5]), 1 - float(p[5])),
                    "n": int(p[7]), "beta": float(p[9]), "se": float(p[10]),
                })
    return pd.DataFrame(rows)


def extract_finngen(chrom, lo, hi):
    """FinnGen 区域流式提取（hg38）。"""
    rows = []
    with gzip.open(FINNGEN, "rt") as fh:
        for line in fh:
            if line.startswith("#chrom"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            if p[0] != str(chrom):
                continue
            pos = int(p[1])
            if not (lo <= pos <= hi):
                continue
            rows.append({
                "pos_hg38": pos, "ref": p[2], "alt": p[3],
                "beta": float(p[8]), "se": float(p[9]), "af_alt": float(p[10]),
            })
    return pd.DataFrame(rows)


def extract_ieu(chrom, lo, hi):
    """IEU API 区间查询（hg19）；10kb 分块避免结果截断。"""
    token = TOKEN_FILE.read_text().strip()
    out = []
    chunk = 10_000
    start = lo
    while start <= hi:
        end = min(start + chunk - 1, hi)
        q = f"id=ieu-b-69&variant={chrom}:{start}-{end}"
        cmd = ["curl.exe", "-s", "-X", "POST",
               "-H", f"Authorization: Bearer {token}",
               "--max-time", "300", f"https://api.opengwas.io/api/associations?{q}"]
        r = subprocess.run(cmd, capture_output=True, text=True,
                           encoding="utf-8", errors="replace")
        if r.returncode != 0:
            raise RuntimeError(r.stderr[-300:])
        data = json.loads(r.stdout)
        out.extend(data)
        start = end + 1
    return pd.DataFrame([{
        "pos_hg19": int(d["position"]), "ea": d.get("ea"), "nea": d.get("nea"),
        "beta": d["beta"], "se": d["se"], "eaf": d.get("eaf"),
    } for d in out])


def is_pal(a, b):
    return {a.upper(), b.upper()} in ({"A", "T"}, {"C", "G"})


def main():
    leads = get_leads()
    print("leads:", {g: leads[g]["rsid"] for g in leads})
    for gene, lead in leads.items():
        chrom = GENES[gene]
        dec = extract_decode(gene, lead)
        lo, hi = lead["pos_hg38"] - WINDOW, lead["pos_hg38"] + WINDOW
        fg = extract_finngen(chrom, lo, hi)
        lo19, hi19 = lead["pos_hg19"] - WINDOW, lead["pos_hg19"] + WINDOW
        ieu = extract_ieu(chrom, lo19, hi19)
        print(f"[{gene}] deCODE={len(dec)} FinnGen={len(fg)} IEU={len(ieu)}")

        for outcome, od, key19 in (
            ("FinnGen_R10_ARDS", fg, False),
            ("IEU_sepsis_ieu-b-69", ieu, True),
        ):
            od = od.copy()
            if outcome == "FinnGen_R10_ARDS":
                od["pos_hg19"] = od["pos_hg38"]  # 占位，实际用 hg38 键
            merged = dec.merge(
                od, left_on="pos_hg19" if key19 else "pos_hg38",
                right_on="pos_hg19" if key19 else "pos_hg38", suffixes=("_e", "_o"))
            if len(merged) == 0:
                print(f"  {outcome}: 0 shared SNPs")
                continue
            keep, flipped = [], 0
            for _, r in merged.iterrows():
                a1_e = str(r["a1"]).upper()
                if outcome == "FinnGen_R10_ARDS":
                    oa, ob = str(r["ref"]).upper(), str(r["alt"]).upper()
                    beta_o, se_o = r["beta_o"], r["se_o"]
                    maf_o = min(r["af_alt"], 1 - r["af_alt"])
                else:
                    oa, ob = str(r["nea"]).upper(), str(r["ea"]).upper()
                    beta_o, se_o = r["beta_o"], r["se_o"]
                    maf_o = r["eaf"]
                    if maf_o is not None:
                        maf_o = min(maf_o, 1 - maf_o)
                flip = False
                if a1_e == ob:
                    pass
                elif a1_e == oa:
                    flip = True
                elif a1_e == ob.translate(COMP):
                    pass
                elif a1_e == oa.translate(COMP):
                    flip = True
                else:
                    continue
                if is_pal(oa, ob):
                    if maf_o is None or maf_o > 0.42 or r["maf"] > 0.42:
                        continue  # 回文且不可判向
                if flip:
                    beta_o = -beta_o
                    flipped += 1
                snp = f"{chrom}:{r['pos_hg38']}"
                keep.append({
                    "snp": snp, "pos_hg38": r["pos_hg38"],
                    "pqtl_beta": r["beta_e"], "pqtl_se": r["se_e"],
                    "pqtl_maf": r["maf"], "pqtl_n": r["n"],
                    "out_beta": beta_o, "out_se": se_o,
                    "out_maf": maf_o,
                })
            out_df = pd.DataFrame(keep)
            out_df.to_csv(INT / f"M3_coloc_{gene}_{outcome}.csv", index=False)
            print(f"  {outcome}: {len(out_df)} shared (flipped {flipped})")


if __name__ == "__main__":
    main()
