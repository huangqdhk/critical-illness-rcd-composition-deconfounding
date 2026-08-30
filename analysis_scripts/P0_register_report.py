# -*- coding: utf-8 -*-
"""登记 P0_Execution_Completion_Report_20260821.md 入 manifest。"""
import csv, hashlib, os, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
MAN = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "RESULTS_MANIFEST_v2.0.csv")
p = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "P0_Execution_Completion_Report_20260821.md")
h = hashlib.sha256()
with open(p, "rb") as f:
    for b in iter(lambda: f.read(1 << 20), b""):
        h.update(b)
with open(MAN, encoding="utf-8-sig", newline="") as f:
    rows = list(csv.DictReader(f))
    FIELDS = rows[0].keys()
exists = any(r["filename"] == "P0_Execution_Completion_Report_20260821.md" for r in rows)
if not exists:
    rows.append({
        "s_number": "", "filename": "P0_Execution_Completion_Report_20260821.md",
        "location": "04_AUDIT_GOVERNANCE", "class": "audit_report", "status": "canonical",
        "gene_universe": "none", "gene_set_version": "NA", "score_version": "NA",
        "n_data_rows": "", "size_bytes": str(os.path.getsize(p)),
        "sha256": h.hexdigest(),
        "analysis": "分析计划执行完成状态总报告（2026-08-21；Phase 0-2 完成 / Phase 3-4 湿实验边界声明）",
        "notes": "",
    })
    print("registered new row")
with open(MAN, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\r\n")
    w.writeheader()
    w.writerows(rows)
print("DONE")
