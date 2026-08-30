# -*- coding: utf-8 -*-
"""
P0_refs_gbt_fetch.py — GB/T 7714-2025 落盘准备：批量抓取全部条目著录数据（2026-08-29）
输出 _intermediate/P0_refs_gbt_metadata.json：
  - 全部 PMID 条目的作者/题名/刊名(ISO缩写)/年/卷/期/页（esummary）
  - 5 个 GEO accession 的题名/发布日期（db=gds esummary）
  - 新增条目核实：Yao 2023 (36747806)、Overmyer 2020 (esearch)
"""
import json, sys, time, urllib.request, urllib.parse
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def get(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read().decode("utf-8", errors="replace")

def esummary_ids(db, ids):
    q = urllib.parse.urlencode({"db": db, "id": ",".join(ids), "retmode": "json"})
    return json.loads(get(f"{EUTILS}/esummary.fcgi?{q}"))["result"]

def esearch(term, db="pubmed", retmax=5):
    q = urllib.parse.urlencode({"db": db, "term": term, "retmax": retmax, "retmode": "json"})
    return json.loads(get(f"{EUTILS}/esearch.fcgi?{q}"))["esearchresult"]

# ---------- 1. 新增条目核实 ----------
print("== 新增条目核实 ==", flush=True)
yao = esummary_ids("pubmed", ["36747806"])["36747806"]
print("Yao 36747806:", yao.get("source"), yao.get("pubdate")[:4], "|", yao.get("title")[:90], flush=True)

over_res = esearch("Overmyer[au] AND (Cell Syst[jour] OR 'cell systems'[jour])")
print("Overmyer Cell Syst 命中:", over_res["idlist"], flush=True)
over_uid = over_res["idlist"]
over_pmid = None
if over_uid:
    sm = esummary_ids("pubmed", over_uid[:3])
    for uid in over_uid[:3]:
        d = sm[uid]
        if "COVID" in d.get("title", "") or "covid" in json.dumps(d).lower():
            over_pmid = uid
            print("  选定:", uid, "|", d.get("source"), d.get("pubdate")[:4], "|", d.get("title")[:90], flush=True)
time.sleep(0.4)

# ---------- 2. 全部 PMID 条目著录数据 ----------
PMID_MAP = {
    # 编号槽 -> PMID（'A:' 前缀=作者-年份条目）
    1:"26903337", 2:"22797452", 3:"30872586", 4:"28792873", 5:"10793167",
    6:"7075161", 7:"12133657", 8:"42487974", 9:"17178908", 10:"22850819",
    11:"22992592", 12:"12974968", 13:"42030446", 14:"40555044", 15:"41383584",
    16:"19740426", 17:"22195746", 19:"37307952", 21:"32461654",
    23:"32398875", 24:"33657410", 25:"35027333", 32:"29195078", 33:"36653562",
    34:"27543902", 35:"34475573", 36:"37962376",
    54:"32416070", 55:"29078758", 56:"24185508", 57:"18650918", 58:"20644254",
    60:"22632970", 61:"28985560", 62:"27535533", 63:"15514116", 64:"22461369",
    66:"38377798", 68:"36830965", 70:"34857953", 71:"37794186", 72:"32913098",
    73:"32492406", 75:"33503446", 76:"29457794",
    "A:Wang":"41317732", "A:NCCD":"29362479", "A:Casadio":"41540269", "A:Kapp":"16613834",
    "N:Yao":"36747806",
}
if over_pmid:
    PMID_MAP["N:Overmyer"] = over_pmid

slots = list(PMOD := PMID_MAP.keys())
all_pmids = list(PMID_MAP.values())
meta = {}
for i in range(0, len(all_pmids), 40):
    batch = all_pmids[i:i+40]
    res = esummary_ids("pubmed", batch)
    for slot in slots:
        pmid = PMID_MAP[slot]
        if pmid in res:
            d = res[pmid]
            meta[str(slot)] = dict(
                pmid=pmid,
                authors=[f"{a.get('name','')}" for a in d.get("authors", [])],
                title=d.get("title", "").rstrip("."),
                journal=d.get("source", ""),
                year=(d.get("pubdate", "") or "")[:4],
                volume=d.get("volume", ""), issue=d.get("issue", ""),
                pages=d.get("pages", ""),
                doi=(d.get("elocationid", "") or ""),
            )
    time.sleep(0.5)
print(f"esummary 完成：{len(meta)} 条", flush=True)

# ---------- 3. GEO accession 元数据（db=gds） ----------
GEO_ACCS = ["GSE165659", "GSE212865", "GSE271370", "GSE310929", "GSE148871"]
geo_meta = {}
for acc in GEO_ACCS:
    try:
        r = esearch(acc + "[ACCN]", db="gds", retmax=3)
        ids = r["idlist"]
        if ids:
            d = esummary_ids("gds", ids[:1])[ids[0]]
            geo_meta[acc] = dict(uid=ids[0], title=d.get("title", ""), summary=(d.get("summary", "") or "")[:200],
                                 pdat=d.get("pdat", ""), n_samples=d.get("n_samples", ""))
            print(f"{acc}: {d.get('pdat','')} | {d.get('title','')[:80]}", flush=True)
        else:
            print(f"{acc}: gds 无命中", flush=True)
    except Exception as e:
        print(f"{acc}: ERROR {e}", flush=True)
    time.sleep(0.5)

out = dict(pubmed=meta, geo=geo_meta, overmyer_pmid=over_pmid,
           yao=dict(journal=yao.get("source"), year=yao.get("pubdate")[:4], title=yao.get("title")))
with open(r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate\P0_refs_gbt_metadata.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("saved _intermediate/P0_refs_gbt_metadata.json", flush=True)
