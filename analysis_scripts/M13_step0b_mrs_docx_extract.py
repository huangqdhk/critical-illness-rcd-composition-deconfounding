# -*- coding: utf-8 -*-
"""Extract the 35-gene MRS list from PMC13376262 DataSheet1.docx."""
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import docx  # noqa: E402

doc = docx.Document(r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\00_RAW_DATA\PMC13376262_DataSheet1.docx")
all_text = []
for p in doc.paragraphs:
    t = p.text.strip()
    if t:
        all_text.append(t)
# tables
tables_text = []
for ti, tb in enumerate(doc.tables):
    rows = []
    for row in tb.rows:
        cells = [c.text.strip() for c in row.cells]
        rows.append(cells)
    tables_text.append((ti, len(rows), rows))
print("paragraphs:", len(all_text), "tables:", len(tables_text))
for ti, nr, rows in tables_text:
    print(f"--- table {ti} ({nr} rows) head: {rows[:3]}")

# look for a column of gene symbols (e.g. "gene", "Genes", single words of 2-10 letters)
gene_re = re.compile(r"^[A-Za-z][A-Za-z0-9-]{1,15}$")
candidates = {}
for ti, nr, rows in tables_text:
    if nr < 25 or nr > 500:
        continue
    for ci in range(min(12, len(rows[0]))):
        col = [r[ci] if ci < len(r) else "" for r in rows]
        hit = [c for c in col if gene_re.match(c)]
        frac = len(hit) / max(1, len(col))
        if frac > 0.6 and len(hit) > 10:
            candidates[(ti, ci)] = (frac, hit)
            print(f"candidate table {ti} col {ci}: {len(hit)} genes, frac {frac:.2f}")
            print("  ", hit[:45])
print("DONE")
