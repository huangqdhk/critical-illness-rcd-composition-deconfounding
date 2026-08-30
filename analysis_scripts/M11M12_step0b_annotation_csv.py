# -*- coding: utf-8 -*-
"""M11M12_step0b_annotation_csv.py — 由 GSM 官方元数据生成逐样本注释 CSV（审计脚本）。"""
import sys, re
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import pandas as pd

RAW = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\00_RAW_DATA"

def parse_gsm(path):
    txt = open(path, encoding="utf-8", errors="replace").read()
    recs = txt.split("^SAMPLE")
    out = []
    for r in recs[1:]:
        d = {}
        for ln in r.split("\n"):
            ln = ln.strip()
            if ln.startswith("!Sample_geo_accession"):
                d["gsm"] = ln.split("=", 1)[1].strip()
            elif ln.startswith("!Sample_title"):
                d["title"] = ln.split("=", 1)[1].strip()
            elif ln.startswith("!Sample_description"):
                d.setdefault("descriptions", []).append(ln.split("=", 1)[1].strip())
            elif ln.startswith("!Sample_characteristics_ch1"):
                kv = ln.split("=", 1)[1].strip()
                if ":" in kv:
                    k, v = kv.split(":", 1)
                    d[k.strip().lower()] = v.strip()
        out.append(d)
    return pd.DataFrame(out)

# GSE54514
g = parse_gsm(RAW + r"\GSE54514_Sepsis_PAXgene_WholeBlood\GSE54514_gsm_metadata.txt")
def tparse(t):
    m = re.search(r",\s*(.+?),\s*Day_(\d+),\s*ID=(\S+)", str(t))
    return (m.group(1), "Day_" + m.group(2), m.group(3)) if m else ("", "", "")
g[["group", "day", "patient_id"]] = g["title"].apply(lambda t: pd.Series(tparse(t)))
g["sentrix"] = g["descriptions"].apply(lambda ds: ds[0] if ds else "")
g[["gsm", "title", "patient_id", "group", "day", "gender", "age (years)",
   "neutrophil proportion", "sentrix"]].to_csv(
    RAW + r"\GSE54514_Sepsis_PAXgene_WholeBlood\GSE54514_sample_annotation.csv", index=False)
print("GSE54514:", g.shape[0], "行")

# GSE106878
g2 = parse_gsm(RAW + r"\GSE106878_Sepsis_Hydrocortisone_CORTICUS\GSE106878_gsm_metadata.txt")
g2["sentrix"] = g2["descriptions"].apply(lambda ds: ds[1] if len(ds) > 1 else (ds[0] if ds else ""))
cols = [c for c in ["gsm", "title", "individual", "timepoint", "treatment", "gender", "age",
                    "acth", "survival (28 days)", "pre-treatment ifng il10 ratio", "sentrix"]
        if c in g2.columns]
g2[cols].to_csv(RAW + r"\GSE106878_Sepsis_Hydrocortisone_CORTICUS\GSE106878_sample_annotation.csv", index=False)
print("GSE106878:", g2.shape[0], "行 | 臂:", g2["treatment"].value_counts().to_dict(),
      "| 时点:", g2["timepoint"].value_counts().to_dict())
