# -*- coding: utf-8 -*-
"""
P0_verify_references_round5.py — 第五轮：[6]/[57] 专项候选检索（2026-08-29）
策略变化：字段限定 [ti]/[tiab] + 旧式词汇（pneumocyte、MeSH）+ 经典作者锚点；[57] 面向
NR1H3(LXRα)/PPARG/EPAS1(HIF-2α) 三因子的髓系功能文献而非基因名检索。
"""
import json, sys, time, urllib.request, urllib.parse
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def esearch(term, retmax=6):
    q = urllib.parse.urlencode({"db": "pubmed", "term": term, "retmax": retmax, "retmode": "json"})
    with urllib.request.urlopen(f"{EUTILS}/esearch.fcgi?{q}", timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))

def esummary(pmid):
    q = urllib.parse.urlencode({"db": "pubmed", "id": pmid, "retmode": "json"})
    with urllib.request.urlopen(f"{EUTILS}/esummary.fcgi?{q}", timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))["result"][pmid]

TERMS = {
    6: [
        # A. 字段限定、放宽为"线粒体+超微/EM+AT2 细胞"（不强制 ARDS 同现）
        'mitochondr*[tiab] AND (ultrastructure[tiab] OR microscopy[tiab] OR ultrastructural[tiab]) '
        'AND ("type II pneumocyte*"[tiab] OR "type 2 pneumocyte*"[tiab] OR "alveolar epithelial type II"[tiab] '
        'OR "alveolar type II"[tiab] OR "type II alveolar"[tiab])',
        # B. 标题层：线粒体 × ARDS/ALI
        'mitochondr*[ti] AND (ARDS[ti] OR "acute respiratory distress"[ti] OR "acute lung injury"[ti])',
        # C. MeSH：ARDS 病理 × 线粒体病理
        '"Respiratory Distress Syndrome, Adult"[MeSH] AND Mitochondria[MeSH]',
        # D. 旧词汇：pneumocyte 标题
        'pneumocyte*[ti] AND mitochondr*[ti]',
        # E. 主张核心：AT2 线粒体退变 × 预后/结局
        'mitochondr*[tiab] AND (alveolar[tiab] OR pneumocyte*[tiab]) AND (ARDS[tiab] OR "respiratory distress"[tiab]) '
        'AND (prognos*[tiab] OR outcome*[tiab] OR survival[tiab] OR mortalit*[tiab])',
        # F. 经典作者锚点：ARDS 肺超微结构经典文献（Bachofen/Weibel 系）
        'Bachofen E[au] AND (ultrastructure[ti] OR structural[ti])',
        # G. COVID 尸检 EM（AT2 线粒体改变的时代文献）
        '(COVID[ti] OR SARS-CoV-2[ti]) AND ultrastructure[ti] AND mitochondr*[ti]',
        # H. 肺损伤模型 AT2 线粒体超微结构（放宽病种）
        '"lung injury"[ti] AND mitochondr*[ti] AND (ultrastructur*[ti] OR microscopy[ti])',
    ],
    57: [
        # A. LXR × 巨噬细胞脂质处理（综述优先）
        '(NR1H3[tiab] OR "liver X receptor"[tiab] OR LXR[tiab]) AND macrophage*[tiab] '
        'AND (lipid[tiab] OR cholesterol[tiab])',
        # B. PPARγ × 巨噬细胞抗炎许可/替代激活
        'PPARG[tiab] AND macrophage*[ti] AND (anti-inflammatory[tiab] OR "alternative activation"[tiab] OR polarization[tiab])',
        # C. HIF-2α/EPAS1 × 髓系
        '(EPAS1[tiab] OR "HIF-2alpha"[tiab] OR "HIF2alpha"[tiab] OR "hypoxia-inducible factor 2"[tiab]) '
        'AND (macrophage*[tiab] OR myeloid[tiab] OR monocyte*[tiab])',
        # D. 巨噬细胞免疫代谢综述（可同时覆盖三因子）
        'macrophage*[ti] AND immunometabolism[ti]',
        # E. 巨噬细胞极化转录控制综述
        'transcriptional[ti] AND macrophage*[ti] AND polarization[ti]',
        # F. 三因子同现
        '(LXR[tiab] AND PPAR[tiab] AND (HIF[tiab] OR hypoxia[tiab]) AND macrophage*[tiab])',
        # G. 脓毒症髓系代谢重编程
        'sepsis[tiab] AND (immunometabolism[tiab] OR "metabolic reprogramming"[tiab]) AND macrophage*[tiab]',
        # H. MeSH：PPARγ × 巨噬细胞 × 炎症
        '"PPAR gamma"[MeSH] AND Macrophages[MeSH] AND inflammation',
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
            if not ids:
                rows.append(dict(ref=num, term_tag=tag, pmid="", title="", journal="", year="", status="NO_HIT", term=term))
                print(f"[{tag}] NO HIT (count={n_hits})", flush=True)
            else:
                print(f"[{tag}] count={n_hits}", flush=True)
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

import pandas as pd
df = pd.DataFrame(rows)
out = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate\P0_reference_verification_round5.csv"
df.to_csv(out, index=False, encoding="utf-8-sig")
print("saved", out, "rows:", len(df))
