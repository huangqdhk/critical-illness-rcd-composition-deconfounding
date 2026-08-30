# -*- coding: utf-8 -*-
"""
P0_verify_references_round4.py — 第四轮：剩余条目的经典文献锚定检索
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

FIXES = {
    6: "mitochondrial ultrastructure alveolar epithelium ARDS lung injury electron microscopy 2020:2024[dp]",
    7: "Brealey D[Author] AND mitochondrial dysfunction severity outcome septic shock",
    8: "COVID-19 mitochondria oxidative phosphorylation downregulation host response 2020:2022[dp]",
    12: "Tang PS[Author] AND Fas-mediated apoptosis alveolar epithelium lung injury",
    13: "RIPK3 necroptosis ARDS acute lung injury review 2018:2024[dp]",
    55: "sepsis blood transcriptome mitochondrial respiratory chain downregulation 2015:2021[dp]",
    56: "Singer M[Author] AND mitochondrial dysfunction sepsis multi-organ failure",
    57: "PPARG macrophage lipid metabolism sepsis immunometabolism review 2018:2024[dp]",
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
        for pmid in ids[:4]:
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
out = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate\P0_reference_verification_round4.csv"
df.to_csv(out, index=False, encoding="utf-8-sig")
print("saved", out)
