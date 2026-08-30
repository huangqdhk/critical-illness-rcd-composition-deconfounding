# -*- coding: utf-8 -*-
"""
P0_geo_citations.py — 从 GEO 获取各数据集 Series 记录的科学论文 PMID/标题（文本解析版）
GDS efetch 返回文本格式（忽略 retmode=xml），从中解析 "PubMed" 行。
"""
import json, os, re, time, urllib.request, urllib.parse, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ACC = ["GSE145926", "GSE158055", "GSE185263", "GSE212865", "GSE32707",
       "GSE165659", "GSE67530", "GSE271370", "GSE310929", "GSE188309",
       "GSE148871", "GSE253474", "GSE235046"]

def get(url, timeout=60):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read().decode("utf-8", errors="replace")

rows = []
for acc in ACC:
    try:
        q = urllib.parse.urlencode({"db": "gds", "term": f"{acc}[ACCN] AND gse[ETYP]", "retmode": "json"})
        res = json.loads(get(f"{EUTILS}/esearch.fcgi?{q}"))
        uids = res["esearchresult"]["idlist"]
        if not uids:
            rows.append(dict(accession=acc, uid="", title="", pmid="", note="no GSE record found"))
            continue
        uid = uids[0]
        q3 = urllib.parse.urlencode({"db": "gds", "id": uid, "retmode": "text"})
        txt = get(f"{EUTILS}/efetch.fcgi?{q3}")
        title = txt.split("\n")[0].strip().rstrip(".")[:250]
        m = re.search(r"PubMed\s+(\d+)", txt)
        pmid = m.group(1) if m else ""
        rows.append(dict(accession=acc, uid=uid, title=title, pmid=pmid, note="" if pmid else "no PMID in record"))
        print(f"{acc}: PMID={pmid} | {title[:100]}", flush=True)
        time.sleep(0.4)
    except Exception as e:
        rows.append(dict(accession=acc, uid="", title="", pmid="", note=f"ERROR {e}"))
        print(f"{acc}: ERROR {e}", flush=True)

import pandas as pd
df = pd.DataFrame(rows)
out = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate\P0_geo_series_pmids.csv"
df.to_csv(out, index=False, encoding="utf-8-sig")
print("saved", out)
