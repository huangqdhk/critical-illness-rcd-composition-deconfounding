# -*- coding: utf-8 -*-
"""
P0_verify_references_round5d.py — 第五轮d：[8][12][13][55] 候选摘要核实 + [12][55] 补充检索（2026-08-29）
"""
import json, sys, time, re, urllib.request, urllib.parse
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def esearch(term, retmax=6):
    q = urllib.parse.urlencode({"db": "pubmed", "term": term, "retmax": retmax,
                                "sort": "relevance", "retmode": "json"})
    with urllib.request.urlopen(f"{EUTILS}/esearch.fcgi?{q}", timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))

def esummary(pmid):
    q = urllib.parse.urlencode({"db": "pubmed", "id": pmid, "retmode": "json"})
    with urllib.request.urlopen(f"{EUTILS}/esummary.fcgi?{q}", timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))["result"][pmid]

def efetch_abstract(pmid, maxlen=1100):
    q = urllib.parse.urlencode({"db": "pubmed", "id": pmid, "rettype": "abstract", "retmode": "text"})
    with urllib.request.urlopen(f"{EUTILS}/efetch.fcgi?{q}", timeout=60) as r:
        return re.sub(r"\s+", " ", r.read().decode("utf-8", errors="replace"))[:maxlen]

# ---- 1. 补充检索（[12] 凋亡×ARDS 经典；[55] 脓毒症血转录组×线粒体下调）----
SEARCHES = {
    12: [
        'apoptosis[ti] AND ("acute lung injury"[ti] OR "acute respiratory distress"[ti] OR ARDS[ti])',
        '"alveolar epithelial"[tiab] AND apoptosis[tiab] AND (ARDS[tiab] OR "lung injury"[tiab]) AND (human[tiab] OR patients[tiab])',
        'Matute-Bello[au] AND apoptosis[ti]',
    ],
    55: [
        'sepsis[ti] AND (transcriptome[ti] OR transcriptomic[ti] OR microarray[ti] OR "gene expression"[ti])',
        'sepsis[tiab] AND "peripheral blood"[tiab] AND (mitochondr*[tiab] OR "oxidative phosphorylation"[tiab] OR "respiratory chain"[tiab]) AND (downregul*[tiab] OR suppress*[tiab] OR reduced[tiab] OR decreased[tiab])',
        '"critical illness"[tiab] AND mitochondri*[ti] AND (blood[tiab] OR leukocyte*[tiab] OR muscle[tiab]) AND transcript*[tiab]',
    ],
}
supp_rows = []
for num, terms in SEARCHES.items():
    for ti, term in enumerate(terms):
        tag = f"{num}.s{ti+1}"
        try:
            res = esearch(term)
            ids = res["esearchresult"]["idlist"]
            print(f"[{tag}] count={res['esearchresult']['count']}", flush=True)
            for pmid in ids:
                sm = esummary(pmid)
                supp_rows.append(dict(ref=num, term_tag=tag, pmid=pmid,
                                      title=sm.get("title", "")[:160],
                                      journal=sm.get("fulljournalname", "") or sm.get("source", ""),
                                      year=(sm.get("pubdate", "") or "")[:4]))
                print(f"    PMID={pmid} ({supp_rows[-1]['year']}) :: {supp_rows[-1]['title'][:110]}", flush=True)
                time.sleep(0.35)
        except Exception as e:
            print(f"[{tag}] ERROR {e}", flush=True)
        time.sleep(0.4)

import pandas as pd
pd.DataFrame(supp_rows).to_csv(
    r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate\P0_reference_verification_round5d_search.csv",
    index=False, encoding="utf-8-sig")

# ---- 2. 全部候选 + 🟡 待确认条目的摘要 ----
ABSTRACTS = {
    8:  ["42487974", "34594344"],
    12: ["42495753", "41383584", "41007859"],
    13: ["42030446", "40188179", "40023302"],
    55: ["42011032", "41704334", "29078758"],
    14: ["40555044"], 15: ["41383584"], 16: ["19740426"], 17: ["22195746"], 56: ["24185508"],
}
print("\n" + "=" * 25, "CANDIDATE ABSTRACTS", "=" * 25, flush=True)
for num, pmids in ABSTRACTS.items():
    for pmid in pmids:
        try:
            print(f"\n===== [{num}] PMID {pmid} =====\n{efetch_abstract(pmid)}", flush=True)
        except Exception as e:
            print(f"\n===== [{num}] PMID {pmid} FETCH ERROR {e}", flush=True)
        time.sleep(0.4)
