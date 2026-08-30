# -*- coding: utf-8 -*-
"""
M1 Step 11: 恢复 RESULTS_MANIFEST 中被 step9 pandas 往返抹掉的 'NA' 版本戳（2026-08-17）
====================================================================================
背景：M1_step9_register_manifest.py 用 `pd.read_csv`（默认 na_values 含字面量 'NA'）
读取并 `to_csv` 重写了整个 manifest，把原本合法的版本戳 'NA'（"不依赖基因集/评分
版本"，见 RESULTS_MANIFEST.md 约定）全部转成空串，导致 lint check_versions 对
170 行报"版本戳 != manifest"。本脚本按表格文件的真实版本格恢复 manifest 值：

- canonical CSV：若 manifest 的 gene_set_version/score_version 为空且表中对应格
  非空（'NA' 或真实戳），以表为准恢复；
- canonical TXT：无版本页脚（原以 manifest 'NA' 跳过页脚校验），gsv/sv 为空时恢复 'NA'；
- GSE67530_beta_for_GrimAge.csv.gz：rawdata_cleanup 移入子目录，文件字节未变
  （SHA256 已核对一致），仅修正 location。

只改 manifest 登记字段，不改任何表文件。
"""
import csv
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "04_AUDIT_GOVERNANCE" / "RESULTS_MANIFEST_v1.0.csv"


def read_rows(p):
    raw = p.read_bytes()
    enc = "utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf-8"
    with open(p, "r", encoding=enc, newline="") as f:
        return [r for r in csv.reader(f)]


def main():
    man = read_rows(MANIFEST)
    hdr, mrows = man[0], man[1:]
    idx = {c: hdr.index(c) for c in hdr}
    n_csv = n_txt = 0
    for r in mrows:
        if r[idx["status"]] != "canonical":
            continue
        fn = r[idx["filename"]]
        loc = r[idx["location"]]
        if not loc:
            continue
        p = ROOT / loc / fn
        if not p.exists():
            continue
        if fn.lower().endswith(".csv"):
            rows = read_rows(p)
            if not rows:
                continue
            hdr_i = 0
            while hdr_i < len(rows) and rows[hdr_i] and str(rows[hdr_i][0]).startswith("#"):
                hdr_i += 1  # '#' 注释行（如 Table_S27b）
            if hdr_i >= len(rows) or "gene_set_version" not in rows[hdr_i] \
                    or "score_version" not in rows[hdr_i]:
                continue
            gi, si = rows[hdr_i].index("gene_set_version"), rows[hdr_i].index("score_version")
            tv = None
            for row in rows[hdr_i + 1:]:
                if row and str(row[0]).startswith("#"):
                    continue
                tv = row
                break
            if tv is None:
                continue
            changed = False
            for col_i, mf in ((gi, "gene_set_version"), (si, "score_version")):
                cell = tv[col_i] if col_i < len(tv) else ""
                if r[idx[mf]] == "" and cell:
                    r[idx[mf]] = cell
                    changed = True
            if changed:
                n_csv += 1
                print(f"[CSV] {fn}: gsv/sv -> {r[idx['gene_set_version']]!r}/{r[idx['score_version']]!r}")
        elif fn.lower().endswith(".txt"):
            if r[idx["gene_set_version"]] == "" and r[idx["score_version"]] == "":
                r[idx["gene_set_version"]] = "NA"
                r[idx["score_version"]] = "NA"
                n_txt += 1
                print(f"[TXT] {fn}: gsv/sv -> NA/NA")

    # gz location 修正（文件字节未变，SHA256 已核对一致）
    for r in mrows:
        if r[idx["filename"]] == "GSE67530_beta_for_GrimAge.csv.gz":
            r[idx["location"]] = "00_RAW_DATA/GSE67530_beta_for_GrimAge.csv"
            print("[GZ] location -> 00_RAW_DATA/GSE67530_beta_for_GrimAge.csv")

    with open(MANIFEST, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(hdr)
        w.writerows(mrows)
    print(f"恢复完成：CSV {n_csv} 行、TXT {n_txt} 行、gz location 1 行；manifest 共 {len(mrows)} 行")


if __name__ == "__main__":
    main()
