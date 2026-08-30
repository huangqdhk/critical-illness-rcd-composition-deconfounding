# -*- coding: utf-8 -*-
"""修复剩余 7 项 lint 失败：4 个 Figure CSV 补版本列；GSE67530 gz 位置；重写 manifest 行。"""
import csv, hashlib, os, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
MAN = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "RESULTS_MANIFEST_v2.0.csv")

def sha256_of(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

def read_rows(p):
    raw = open(p, "rb").read()
    enc = "utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf-8"
    with open(p, encoding=enc, newline="") as f:
        return list(csv.reader(f)), enc

def write_rows(p, rows, enc):
    with open(p, "w", encoding=enc, newline="") as f:
        csv.writer(f, lineterminator="\r\n").writerows(rows)

STAMPS = {
    "Figure_1B.csv": ("01_FIGURE_DATA_CSV/Main", "NA", "NA"),
    "Figure_1C.csv": ("01_FIGURE_DATA_CSV/Main", "NA", "NA"),
    "Figure_2E.csv": ("01_FIGURE_DATA_CSV/Main", "Mitoxy-80_v1.0", "NA"),
    "Figure_S2F.csv": ("01_FIGURE_DATA_CSV/Supplementary", "NA", "NA"),
}
for fn, (loc, g, s) in STAMPS.items():
    p = os.path.join(ROOT, loc, fn)
    rows, enc = read_rows(p)
    if "gene_set_version" not in rows[0]:
        rows[0] = rows[0] + ["gene_set_version", "score_version"]
        for r in rows[1:]:
            r += [g, s]
        write_rows(p, rows, enc)
        print(f"[补列] {fn} -> ({g},{s})")
    else:
        print(f"[跳过] {fn} 已有版本列")

# manifest 行更新
with open(MAN, encoding="utf-8-sig", newline="") as f:
    man = list(csv.DictReader(f))
    FIELDS = man[0].keys()
by_fn = {r["filename"]: r for r in man}
for fn, (loc, g, s) in STAMPS.items():
    r = by_fn[fn]
    r["gene_set_version"], r["score_version"] = g, s
    p = os.path.join(ROOT, loc, fn)
    r["sha256"] = sha256_of(p)
    r["size_bytes"] = str(os.path.getsize(p))
    r["n_data_rows"] = str(sum(1 for _ in open(p, "rb")) - 1)
    print(f"[manifest] {fn} 版本/校验和更新")
# GSE67530 gz：实际位于 00_RAW_DATA/GSE67530_beta_for_GrimAge.csv 目录下
r = by_fn.get("GSE67530_beta_for_GrimAge.csv.gz")
if r:
    r["location"] = "00_RAW_DATA/GSE67530_beta_for_GrimAge.csv"
    print(f"[manifest] GSE67530 gz location -> {r['location']}")

with open(MAN, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\r\n")
    w.writeheader()
    w.writerows(man)
print("DONE fix-remaining")
