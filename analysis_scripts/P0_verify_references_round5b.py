# -*- coding: utf-8 -*-
"""
P0_verify_references_round5b.py — 第五轮b：relevance 排序 + 经典作者锚点（2026-08-29）
变化：esearch 加 sort=relevance（默认按日期，导致前轮被 2026 文刷屏）；[57] 加核受体经典
综述锚点；对入围短名单 efetch 摘要核实内容匹配度。
"""
import json, sys, time, re, urllib.request, urllib.parse
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def esearch(term, retmax=8):
    q = urllib.parse.urlencode({"db": "pubmed", "term": term, "retmax": retmax,
                                "sort": "relevance", "retmode": "json"})
    with urllib.request.urlopen(f"{EUTILS}/esearch.fcgi?{q}", timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))

def esummary(pmid):
    q = urllib.parse.urlencode({"db": "pubmed", "id": pmid, "retmode": "json"})
    with urllib.request.urlopen(f"{EUTILS}/esummary.fcgi?{q}", timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))["result"][pmid]

def efetch_abstract(pmid):
    q = urllib.parse.urlencode({"db": "pubmed", "id": pmid, "rettype": "abstract", "retmode": "text"})
    with urllib.request.urlopen(f"{EUTILS}/efetch.fcgi?{q}", timeout=60) as r:
        return r.read().decode("utf-8", errors="replace")

TERMS = {
    6: [
        'ARDS AND "alveolar type II" AND mitochondri* AND (ultrastructur* OR "electron microscop*")',
        '"acute respiratory distress syndrome" AND "type II pneumocytes" AND mitochondri*',
        'ARDS AND alveolar AND mitochondri* AND (prognosis OR mortality OR outcome) AND (human OR patients)',
        '(diffuse alveolar damage) AND mitochondri* AND (ultrastructur* OR microscop*)',
        '"granular pneumocyte" AND mitochondri* AND (injury OR hyperoxia OR distress)',
        'Crouser ED[au] AND mitochondri*[ti]',
        'Weibel ER[au] AND (distress[ti] OR injury[ti])',
        'mitochondri*[ti] AND (ARDS[ti] OR "acute lung injury"[ti]) AND (prognos*[tiab] OR mortalit*[tiab] OR outcome*[tiab])',
    ],
    57: [
        'Bensinger[au] AND Tontonoz[au]',
        'Castrillo[au] AND Tontonoz[au]',
        '"liver X receptor"[ti] AND macrophage*[ti]',
        'NR1H3[ti] AND macrophage*[tiab]',
        '"HIF-2alpha"[ti] AND macrophage*[ti]',
        'EPAS1[ti] AND (macrophage*[tiab] OR myeloid[tiab])',
        'Imtiyaz[au] AND HIF',
        'macrophage*[ti] AND (immunometabolism[tiab] OR "metabolic reprogramming"[tiab]) AND (nuclear receptor*[tiab] OR hypoxia[tiab] OR LXR[tiab] OR PPAR[tiab])',
    ],
}

rows = []
for num, terms in TERMS.items():
    for ti, term in enumerate(terms):
        tag = f"{num}.{chr(65+ti)}"
        try:
            res = esearch(term)
            ids = res["esearchresult"]["idlist"]
            n_hits = res["esearchresult"]["count"]
            print(f"[{tag}] count={n_hits}", flush=True)
            if not ids:
                rows.append(dict(ref=num, term_tag=tag, pmid="", title="", journal="", year="", status="NO_HIT", term=term))
            for pmid in ids[:6]:
                sm = esummary(pmid)
                title = sm.get("title", "")[:170]
                rows.append(dict(ref=num, term_tag=tag, pmid=pmid, title=title,
                                 journal=sm.get("fulljournalname", "") or sm.get("source", ""),
                                 year=(sm.get("pubdate", "") or "")[:4], status="CAND", term=term))
                print(f"    PMID={pmid} ({rows[-1]['year']}) :: {title[:110]}", flush=True)
                time.sleep(0.35)
        except Exception as e:
            rows.append(dict(ref=num, term_tag=tag, pmid="", title="", journal="", year="", status=f"ERROR {e}", term=term))
            print(f"[{tag}] ERROR {e}", flush=True)
        time.sleep(0.4)

# 入围短名单：抓摘要核实
SHORTLIST = {
    6: ["42592019", "42432720", "42486791"],   # ARDS 线粒体综述类，核实是否论及 AT2 超微结构/预后
    57: ["42439344", "42291298", "42311122"],
}
# 加上作者锚点命中的第一波结果（脚本运行后由 CAND 行补充：Bensinger/Castrillo/Imtiyaz 的 PMID）
seen = set()
for r in rows:
    if r["pmid"] and r["status"] == "CAND":
        for au_anchor in ("6.F", "6.G", "57.A", "57.B", "57.G"):
            if r["term_tag"] == au_anchor and r["pmid"] not in seen:
                SHORTLIST.setdefault(r["ref"], []).append(r["pmid"]); seen.add(r["pmid"])

print("\n" + "="*30, "SHORTLIST ABSTRACTS", "="*30, flush=True)
import pandas as pd
df = pd.DataFrame(rows)
out = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate\P0_reference_verification_round5b.csv"
df.to_csv(out, index=False, encoding="utf-8-sig")
print("saved", out, "rows:", len(df), flush=True)

for num, pmids in SHORTLIST.items():
    for pmid in pmids:
        try:
            ab = efetch_abstract(pmid)
            ab_clean = re.sub(r"\s+", " ", ab)[:1400]
            print(f"\n----- [{num}] PMID {pmid} -----\n{ab_clean}", flush=True)
            time.sleep(0.4)
        except Exception as e:
            print(f"\n----- [{num}] PMID {pmid} FETCH ERROR {e}", flush=True)
