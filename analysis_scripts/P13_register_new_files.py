# -*- coding: utf-8 -*-
"""
P13_register_new_files.py — 将 111提质计划13 新增产出登记入 RESULTS_MANIFEST_v2.0.csv
（幂等：同名行先删后增；sha256/size 实算）
"""
import csv, hashlib, io, os, sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = os.path.dirname(os.path.abspath(__file__))
MAN = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "RESULTS_MANIFEST_v2.0.csv")

NEW = [
    ("P13_Table4_PreReg_Traceability.md", "04_AUDIT_GOVERNANCE", "governance",
     "P13 Table 4 预注册溯源核对报告（17 行逐条；行 1/2/9 重大修正）",
     "111提质计划13 §四 必办事项", "", ""),
    ("Table4_criteria_verdicts.csv", "01_FIGURE_DATA_CSV", "results_table",
     "Table 4（17 条预设判定标准×结果×判决 + 预注册出处/注册状态溯源列）",
     "111提质计划13 §四.3；已同步登记 FIGURE_DATA_MANIFEST.csv",
     "Mitoxy-80_v1.0", "mdi_v1.0"),
]

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

with io.open(MAN, encoding="utf-8-sig", newline="") as f:
    rows = list(csv.DictReader(f))
    fields = list(rows[0].keys())

names = {n for n, *_ in NEW}
body = [r for r in rows if r["filename"] not in names]
for name, loc, cls, desc, note, gsv, scv in NEW:
    p = os.path.join(ROOT, loc, name)
    nrow = 0
    if name.endswith(".csv"):
        with io.open(p, encoding="utf-8-sig") as f:
            nrow = sum(1 for _ in f) - 1
    body.append({
        "s_number": "", "filename": name, "location": loc, "class": cls,
        "status": "canonical", "gene_universe": "none",
        "gene_set_version": gsv, "score_version": scv,
        "n_data_rows": str(nrow) if nrow else "",
        "size_bytes": str(os.path.getsize(p)), "sha256": sha(p),
        "analysis": desc, "notes": note, "class_": ""})
with io.open(MAN, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, lineterminator="\r\n")
    w.writeheader()
    w.writerows(body)
print("registered:", [n for n, *_ in NEW])
