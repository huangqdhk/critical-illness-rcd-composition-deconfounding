# -*- coding: utf-8 -*-
"""同步 RESULTS_MANIFEST.md 的 sha/size 入 v2.0 manifest。"""
import csv, hashlib, os, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
MAN = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "RESULTS_MANIFEST_v2.0.csv")
p = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "RESULTS_MANIFEST.md")
h = hashlib.sha256()
with open(p, "rb") as f:
    for b in iter(lambda: f.read(1 << 20), b""):
        h.update(b)
with open(MAN, encoding="utf-8-sig", newline="") as f:
    rows = list(csv.DictReader(f))
    FIELDS = rows[0].keys()
for r in rows:
    if r["filename"] == "RESULTS_MANIFEST.md":
        r["sha256"] = h.hexdigest()
        r["size_bytes"] = str(os.path.getsize(p))
        print("updated RESULTS_MANIFEST.md:", r["size_bytes"], r["sha256"][:16])
with open(MAN, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\r\n")
    w.writeheader()
    w.writerows(rows)
print("DONE")
