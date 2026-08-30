#!/usr/bin/env python3
"""
M3 腿1 步骤2：结局数据提取与 harmonization
===========================================
结局：IEU OpenGWAS ieu-b-69（脓毒症，hg19，经官方 API 逐 SNP 查询）
      FinnGen R10 J10_ARDS（本地 2.1GB gz 流式提取，hg38）
暴露：M3_step1 输出的 cis-pQTL 工具（deCODE 主 + UKB-PPP 敏感性）
harmonization 口径与既有 twmr_finngen.R / Table S15a/b 一致：
rsID + 等位对匹配（含链翻转）；回文用 EAF 差<=0.2 判定；含糊剔除。
输出：_intermediate/M3_pqtl_harmonized.csv（长表，逐 工具×结局 一行）
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
EAF_TOL = 0.2


def load_instruments():
    df = pd.read_csv(INT / "M3_pqtl_instruments_final.csv")
    return df


def extract_finngen(inst):
    """流式读取 FinnGen ARDS gz，按 (chrom,pos) 提取目标变体行（rsids 列可为空）。"""
    keys = {(str(int(r["chr"])), str(int(r["pos_hg38"]))) for _, r in inst.iterrows()
            if pd.notna(r.get("chr")) and pd.notna(r.get("pos_hg38"))}
    hits = {}
    with gzip.open(FINNGEN, "rt") as fh:
        for line in fh:
            if line.startswith("#chrom"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            key = (p[0], p[1])
            if key in keys:
                hits[key] = {
                    "chrom": p[0], "pos": p[1], "ref": p[2], "alt": p[3],
                    "rsid": p[4], "pval": float(p[6]), "beta": float(p[8]),
                    "sebeta": float(p[9]), "af_alt": float(p[10]),
                }
    return hits


def query_ieu(rsids):
    """OpenGWAS API 批量逐 SNP 查询（POST /api/associations，query 参数）。"""
    token = TOKEN_FILE.read_text().strip()
    params = [("id", "ieu-b-69")]
    for r in rsids:
        params.append(("variant", r))
    qs = "&".join(f"{k}={v}" for k, v in params)
    cmd = [
        "curl.exe", "-s", "-X", "POST",
        "-H", f"Authorization: Bearer {token}",
        "--max-time", "300",
        f"https://api.opengwas.io/api/associations?{qs}",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    if r.returncode != 0:
        raise RuntimeError(f"IEU API failed: {r.stderr[-300:]}")
    try:
        data = json.loads(r.stdout)
    except json.JSONDecodeError:
        raise RuntimeError(f"IEU API bad JSON: {r.stdout[:300]}")
    return {d["rsid"]: d for d in data}


def harmonize(expo, out):
    """单 SNP 等位 harmonize。expo: EA/otherA/beta/se/eaf; out: EA/otherA/beta/se/eaf。
    返回 (keep, flip) 或 None=剔除。"""
    e1, e2 = expo["EA"].upper(), expo["OA"].upper()
    o1, o2 = out["EA"].upper(), out["OA"].upper()
    flip = None
    if e1 == o1 and e2 == o2:
        flip = False
    elif e1 == o2 and e2 == o1:
        flip = True
    else:
        o1c, o2c = o1.translate(COMP), o2.translate(COMP)
        if e1 == o1c and e2 == o2c:
            flip = False
        elif e1 == o2c and e2 == o1c:
            flip = True
        else:
            return None
    pal = {e1, e2} in ({"A", "T"}, {"C", "G"})
    if pal and flip is not None:
        if pd.isna(expo["eaf"]) or pd.isna(out["eaf"]):
            return None
        if abs(expo["eaf"] - out["eaf"]) <= EAF_TOL:
            pass
        elif abs(expo["eaf"] - (1 - out["eaf"])) <= EAF_TOL:
            flip = not flip
        else:
            return None
    return flip


def main():
    print("=" * 70)
    print("M3 step2: outcome extraction + harmonization")
    print("=" * 70)
    inst = load_instruments()
    rsids = sorted(inst["rsid"].dropna().unique().tolist())
    print(f"instruments: {len(inst)} rows, {len(rsids)} unique rsids")

    fg = extract_finngen(inst)
    print(f"FinnGen matched: {len(fg)}/{len(inst)} instrument rows")
    ieu = query_ieu(rsids)
    print(f"IEU matched: {len(ieu)}/{len(rsids)}")

    rows = []
    for _, e in inst.iterrows():
        fg_key = (str(int(e["chr"])), str(int(e["pos_hg38"])))
        for outcome, od, src in (
            ("IEU_sepsis_ieu-b-69", ieu.get(e["rsid"]), "IEU"),
            ("FinnGen_R10_ARDS", fg.get(fg_key), "FinnGen"),
        ):
            if od is None:
                continue
            if src == "IEU":
                out_info = {"EA": od["ea"], "OA": od["nea"], "beta": od["beta"],
                            "se": od["se"], "eaf": od.get("eaf"), "n": od.get("n")}
            else:
                out_info = {"EA": od["alt"], "OA": od["ref"], "beta": od["beta"],
                            "se": od["sebeta"], "eaf": od["af_alt"], "n": 406893}
            flip = harmonize(
                {"EA": e["allele1"], "OA": e["allele0"], "beta": e["beta"],
                 "se": e["se"], "eaf": e["a1freq"]},
                out_info,
            )
            if flip is None:
                continue
            beta_o = -out_info["beta"] if flip else out_info["beta"]
            eaf_o = 1 - out_info["eaf"] if flip else out_info["eaf"]
            rows.append({
                "gene": e["gene"], "source": e["source"], "rsid": e["rsid"],
                "outcome": outcome,
                "EA_e": e["allele1"], "OA_e": e["allele0"],
                "beta_e": e["beta"], "se_e": e["se"], "eaf_e": e["a1freq"], "F_e": e["F"],
                "EA_o": out_info["EA"], "OA_o": out_info["OA"],
                "beta_o": beta_o, "se_o": out_info["se"], "eaf_o": eaf_o,
                "harmonize": "flip" if flip else "same", "n_outcome": out_info.get("n"),
            })
    h = pd.DataFrame(rows)
    h.to_csv(INT / "M3_pqtl_harmonized.csv", index=False)
    print("\nHarmonized table:")
    print(h[["gene", "source", "rsid", "outcome", "beta_e", "beta_o", "se_o", "harmonize"]].to_string(index=False))
    print(f"\nSaved {len(h)} rows -> _intermediate/M3_pqtl_harmonized.csv")
    # 汇总统计
    for oc in ["IEU_sepsis_ieu-b-69", "FinnGen_R10_ARDS"]:
        sub = h[h["outcome"] == oc]
        print(f"{oc}: {len(sub)} instrument-outcome rows, {sub['gene'].nunique()} genes")


if __name__ == "__main__":
    main()
