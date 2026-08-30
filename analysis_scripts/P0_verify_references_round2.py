# -*- coding: utf-8 -*-
"""
P0_verify_references_round2.py — 修复首轮误匹配的参考文献条目
========================================================================
仅重查首轮结果明显错误/可疑的编号；用更精确检索式。正确命中保持不动。
"""
import json, sys, time, urllib.request, urllib.parse
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def esearch(term, retmax=5):
    q = urllib.parse.urlencode({"db": "pubmed", "term": term, "retmax": retmax, "retmode": "json"})
    with urllib.request.urlopen(f"{EUTILS}/esearch.fcgi?{q}", timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))

def esummary(pmid):
    q = urllib.parse.urlencode({"db": "pubmed", "id": pmid, "retmode": "json"})
    with urllib.request.urlopen(f"{EUTILS}/esummary.fcgi?{q}", timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))["result"][pmid]

# 修复列表：编号 -> 精确检索式（含作者/标题锚定）
FIXES = {
    2: "Ranieri VM[Author] AND Berlin Definition acute respiratory distress syndrome JAMA 2012",
    4: "Thompson BT[Author] AND Acute Respiratory Distress Syndrome NEJM 2017",
    5: "Ware LB[Author] AND Matthay MA[Author] AND acute respiratory distress syndrome NEJM 2000",
    6: "mitochondrial ultrastructure alveolar epithelial cells type II acute respiratory distress syndrome electron microscopy prognosis",
    7: "sepsis mitochondrial dysfunction transcriptome peripheral blood mononuclear cells bioenergetics downregulation",
    8: "SARS-CoV-2 infection downregulates mitochondrial genes OXPHOS signature peripheral blood COVID-19",
    12: "apoptosis alveolar epithelial cell death acute lung injury review",
    13: "necroptosis RIPK3 MLKL acute lung injury pneumonia",
    14: "pyroptosis NLRP3 inflammasome gasdermin acute lung injury ARDS",
    15: "ferroptosis acute lung injury pulmonary mechanism review",
    16: "caspase inhibitor clinical trial sepsis human apoptosis inhibition failed",
    17: "necroptosis inhibition RIPK1 clinical sepsis ARDS therapeutic limitation",
    18: "Zhang K[Author] AND cell atlas chromatin accessibility 25 adult human tissues",
    24: "Ren X[Author] AND COVID-19 immune features large-scale single-cell transcriptome atlas Cell 2021",
    27: "microarray gene expression COVID-19 severe high-risk patients dynamic profiling",
    33: "FinnGen provides genetic insights from a well-phenotyped isolated population Nature",
    55: "sepsis whole blood transcriptome mitochondrial energy metabolism dysregulation",
    56: "critical illness mitochondrial dysfunction bioenergetic failure review",
    57: "PPARG LXR macrophage lipid handling hypoxia inducible factor sepsis immune tolerance",
    63: "Nemeth E[Author] AND hepcidin ferroportin iron efflux internalization Science 2004",
}

rows = []
for num, term in FIXES.items():
    try:
        res = esearch(term)
        ids = res["esearchresult"]["idlist"]
        if not ids:
            rows.append(dict(ref=num, pmid="", title="", journal="", year="", status="NO_HIT", term=term))
            print(f"[{num}] NO HIT", flush=True)
            continue
        for pmid in ids[:3]:
            sm = esummary(pmid)
            title = sm.get("title", "")[:150]
            rows.append(dict(ref=num, pmid=pmid, title=title,
                             journal=sm.get("fulljournalname", "") or sm.get("source", ""),
                             year=(sm.get("pubdate", "") or "")[:4], status="CAND", term=term))
            print(f"[{num}] PMID={pmid} ({rows[-1]['year']}) :: {title[:100]}", flush=True)
        time.sleep(0.4)
    except Exception as e:
        rows.append(dict(ref=num, pmid="", title="", journal="", year="", status=f"ERROR {e}", term=term))
        print(f"[{num}] ERROR {e}", flush=True)

import pandas as pd
df = pd.DataFrame(rows)
out = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate\P0_reference_verification_round2.csv"
df.to_csv(out, index=False, encoding="utf-8-sig")
print("saved", out)
