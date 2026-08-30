# -*- coding: utf-8 -*-
"""
P0_verify_references_round5c.py — 第五轮c：[6] 高精度补查（relevance，2026-08-29）
"""
import json, sys, time, urllib.request, urllib.parse
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def esearch(term, retmax=15):
    q = urllib.parse.urlencode({"db": "pubmed", "term": term, "retmax": retmax,
                                "sort": "relevance", "retmode": "json"})
    with urllib.request.urlopen(f"{EUTILS}/esearch.fcgi?{q}", timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))

def esummary(pmid):
    q = urllib.parse.urlencode({"db": "pubmed", "id": pmid, "retmode": "json"})
    with urllib.request.urlopen(f"{EUTILS}/esummary.fcgi?{q}", timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))["result"][pmid]

TERMS = [
    '"alveolar type II"[tiab] AND mitochondr*[tiab] AND (ARDS[tiab] OR "acute respiratory distress"[tiab] OR "acute lung injury"[tiab])',
    '"type II alveolar epithelial"[tiab] AND mitochondr*[tiab] AND (death[tiab] OR mortality[tiab] OR prognosis[tiab] OR outcome*[tiab] OR survival[tiab])',
    'mitochondr*[tiab] AND (swelling[tiab] OR vacuolization[tiab] OR cristae[tiab] OR ultrastructur*[tiab]) AND ("respiratory distress"[tiab] OR "diffuse alveolar damage"[tiab] OR DAD[tiab])',
    '"lung injury"[tiab] AND "electron microscopy"[tiab] AND mitochondr*[tiab] AND (human[tiab] OR patients[tiab] OR autops*[tiab] OR biops*[tiab])',
    'alveolar[ti] AND mitochondr*[ti] AND (injur*[ti] OR distress[ti] OR damage[ti])',
]

rows = []
for ti, term in enumerate(TERMS):
    tag = f"6.{chr(97+ti)}"
    try:
        res = esearch(term)
        ids = res["esearchresult"]["idlist"]
        print(f"[{tag}] count={res['esearchresult']['count']}", flush=True)
        for pmid in ids:
            sm = esummary(pmid)
            title = sm.get("title", "")[:170]
            rows.append(dict(term_tag=tag, pmid=pmid, title=title,
                             journal=sm.get("fulljournalname", "") or sm.get("source", ""),
                             year=(sm.get("pubdate", "") or "")[:4], term=term))
            print(f"    PMID={pmid} ({rows[-1]['year']}) :: {title[:120]}", flush=True)
            time.sleep(0.35)
    except Exception as e:
        print(f"[{tag}] ERROR {e}", flush=True)
    time.sleep(0.4)

import pandas as pd
df = pd.DataFrame(rows)
out = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate\P0_reference_verification_round5c.csv"
df.to_csv(out, index=False, encoding="utf-8-sig")
print("saved", out, "rows:", len(df))
