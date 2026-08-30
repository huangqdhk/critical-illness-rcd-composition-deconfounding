# -*- coding: utf-8 -*-
"""
P0_geo_elink.py — 用 NCBI ELink 将 GDS Series UID 链接到 PubMed PMID
"""
import json, sys, time, urllib.request, urllib.parse
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def get(url, timeout=60):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read().decode("utf-8")

ACC = ["GSE145926", "GSE158055", "GSE185263", "GSE212865", "GSE32707",
       "GSE165659", "GSE67530", "GSE271370", "GSE310929", "GSE188309",
       "GSE148871", "GSE253474", "GSE235046"]

rows = []
for acc in ACC:
    try:
        q = urllib.parse.urlencode({"db": "gds", "term": f"{acc}[ACCN] AND gse[ETYP]", "retmode": "json"})
        res = json.loads(get(f"{EUTILS}/esearch.fcgi?{q}"))
        uids = res["esearchresult"]["idlist"]
        if not uids:
            rows.append(dict(accession=acc, gds_uid="", pmid="", note="no GSE record"))
            continue
        uid = uids[0]
        q2 = urllib.parse.urlencode({"dbfrom": "gds", "db": "pubmed", "id": uid, "retmode": "json"})
        link = json.loads(get(f"{EUTILS}/elink.fcgi?{q2}"))
        pmids = []
        for ls in link.get("linksets", []):
            for lidb in ls.get("linksetdbs", []):
                if lidb.get("linkname") in ("gds_pubmed", "gds_pubmed_citedin"):
                    pmids += lidb.get("links", [])
        pmids = list(dict.fromkeys(pmids))
        rows.append(dict(accession=acc, gds_uid=uid, pmid=",".join(pmids) if pmids else "",
                         note="" if pmids else "no pubmed link"))
        print(f"{acc}: PMIDs={pmids}", flush=True)
        time.sleep(0.4)
    except Exception as e:
        rows.append(dict(accession=acc, gds_uid="", pmid="", note=f"ERROR {e}"))
        print(f"{acc}: ERROR {e}", flush=True)

import pandas as pd
df = pd.DataFrame(rows)
out = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate\P0_geo_elink_pmids.csv"
df.to_csv(out, index=False, encoding="utf-8-sig")
print("saved", out)
