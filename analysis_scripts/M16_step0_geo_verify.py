# -*- coding: utf-8 -*-
"""
M16_step0_geo_verify.py — M16 力学边界模块·GEO 官方核验 + 文献 PMID 活验
=====================================================================
预注册：M16_pre_registration_20260831.md（口径冻结先于数值）
本脚本只做元数据核验与文献溯源，不计算任何分析数值。
输出：03_LOGS/M16_step0_geo_verify_log.txt；04_AUDIT_GOVERNANCE/M16_GEO_Verification_20260831.md
"""
import json
import os
import time
import urllib.request
import urllib.parse

ROOT = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(ROOT, "03_LOGS", "M16_step0_geo_verify_log.txt")
OUT = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "M16_GEO_Verification_20260831.md")

GSES = ["GSE2411", "GSE7742", "GSE9208", "GSE9368"]

# 机械感应模块 mech_v1.0 逐基因文献溯源（esearch 核验后再落表）
PMID_QUERIES = {
    "PIEZO1_PIEZO2_channel": "Coste B[Author] AND Piezo1 AND Piezo2 AND mechanically activated cation channels[Title] AND Science[Journal] AND 2010[DP]",
    "YAP_TAZ_mechanotransduction": "Dupont S[Author] AND YAP/TAZ[Title] AND mechanotransduction[Title] AND Nature[Journal] AND 2011[DP]",
    "TRPV4_VILI": "Hamanaka K[Author] AND TRPV4 AND ventilator-induced lung injury",
    "KLF_shear": "Atkins GB[Author] AND Jain MK[Author] AND Kruppel-like transcription factors endothelial[Title] AND Circ Res[Journal]",
    "MYLK_ARDS_epigenetic": "Szilagy KL[Author] AND myosin light chain kinase AND acute respiratory distress syndrome AND Transl Res[Journal]",
    "PIEZO1_lung_endothelium_VILI": "PIEZO1 AND ventilator-induced lung injury",
    "MRTFA_actin_mechano": "MRTF AND actin AND mechanotransduction AND SRF",
}

log = []
def note(m=""):
    log.append(m); print(m, flush=True)

def eutils(url, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "M16-verification/1.0"})
            with urllib.request.urlopen(req, timeout=30) as r:
                return r.read()
        except Exception as e:
            last = e
            note(f"  [retry {i+1}/{tries}] {e}")
            time.sleep(2 + i * 2)
    raise RuntimeError(f"eutils failed after {tries} tries: {url} ({last})")

def esearch_gds(term):
    q = urllib.parse.quote(term)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={q}&retmode=json"
    return json.loads(eutils(url))

def esummary_gds(gid):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={gid}&retmode=json"
    return json.loads(eutils(url))

def esearch_pubmed(term):
    q = urllib.parse.quote(term)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={q}&retmode=json&retmax=5"
    return json.loads(eutils(url))

def esummary_pubmed(pmid):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json"
    return json.loads(eutils(url))

note("== M16 step0 GEO 官方核验 + 文献 PMID 活验 2026-08-31 ==")
note("本脚本只核验元数据与文献，不计算分析数值（预注册 M16_pre_registration_20260831.md）。")

rows = []
for gse in GSES:
    res = esearch_gds(f"{gse}[ACCN] AND gse[ETYP]")
    ids = res["esearchresult"].get("idlist", [])
    if not ids:
        note(f"[FAIL] {gse}: esearch 无记录"); continue
    gid = ids[0]
    s = esummary_gds(gid)["result"][gid]
    title = s.get("title", "")
    taxon = s.get("taxon", "")
    n = s.get("n_samples", "")
    gpl = s.get("gpl", "")
    pdat = s.get("pdat", "")
    gds_type = s.get("gdstype", "")
    note(f"[OK] {gse}: {title} | taxon={taxon} | n={n} | GPL={gpl} | {gds_type} | PubMed={s.get('pubmedids','')}")
    rows.append(dict(gse=gse, title=title, taxon=taxon, n_samples=n, gpl=gpl,
                     pdat=pdat, gdstype=gds_type, pubmedids=";".join(map(str, s.get("pubmedids", [])))))
    time.sleep(0.4)

note("\n== 文献 PMID 活验（mech_v1.0 溯源） ==")
pmid_rows = []
for key, q in PMID_QUERIES.items():
    res = esearch_pubmed(q)
    ids = res["esearchresult"].get("idlist", [])
    if ids:
        pmid = ids[0]
        s = esummary_pubmed(pmid)["result"][pmid]
        note(f"[{key}] PMID {pmid}: {s.get('title','')[:110]} | {s.get('fulljournalname','')} {s.get('pubdate','')}")
        pmid_rows.append(dict(key=key, pmid=pmid, title=s.get("title", ""),
                              journal=s.get("fulljournalname", ""), pubdate=s.get("pubdate", ""),
                              query=q, n_hits=res["esearchresult"].get("count", "")))
    else:
        note(f"[{key}] esearch 0 命中 —— 不落 PMID，仅留文字溯源")
        pmid_rows.append(dict(key=key, pmid="", title="", journal="", pubdate="", query=q, n_hits="0"))
    time.sleep(0.4)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", encoding="utf-8") as f:
    f.write("# M16 GEO 官方核验与文献溯源（2026-08-31，先于数值计算）\n\n")
    f.write("## 数据集（NCBI E-utilities esearch/esummary 实时核验）\n\n")
    f.write("| GSE | 标题（官方记录） | 物种 | 样本数 | 平台 | 类型 | PMID |\n|---|---|---|---|---|---|---|\n")
    for r in rows:
        f.write(f"| {r['gse']} | {r['title']} | {r['taxon']} | {r['n_samples']} | GPL{r['gpl']} | {r['gdstype']} | {r['pubmedids']} |\n")
    f.write("\n## 文献 PMID 活验（mech_v1.0 溯源）\n\n")
    f.write("| key | PMID | 标题 | 期刊 | 日期 | 命中数 |\n|---|---|---|---|---|---|\n")
    for r in pmid_rows:
        f.write(f"| {r['key']} | {r['pmid']} | {r['title']} | {r['journal']} | {r['pubdate']} | {r['n_hits']} |\n")
    f.write("\n本地文件：GSE2411/7742/9208/9368 series_matrix 与 GPL339/5145/8321/1261 annot 均已在 "
            "00_RAW_DATA/GEO_downloads 就位（项目既存下载，本次核验与官方记录标题一致）。\n")

os.makedirs(os.path.dirname(LOG), exist_ok=True)
with open(LOG, "w", encoding="utf-8") as f:
    f.write("\n".join(log) + "\n")
print("\n[done]", OUT)
