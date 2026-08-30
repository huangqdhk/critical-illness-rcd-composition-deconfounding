# -*- coding: utf-8 -*-
"""
P0_compile_reference_list.py — 合并四轮核实结果，生成正式参考文献清单与核实报告
============================================================================
- 已核实条目写入 v4 文稿 References；
- 无法唯一锚定的条目标注"待作者确认"，绝不虚构；
- 从未在 v4 正文引用的编号（19,20,22,26,28-31,37-53,58,59）登记为"未引用、无法重建"
  （v1 文稿不在项目内，其 1-64 清单在所有版本中均为占位）。
"""
import json, sys, time, urllib.request, urllib.parse
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def esearch(term, retmax=4):
    q = urllib.parse.urlencode({"db": "pubmed", "term": term, "retmax": retmax, "retmode": "json"})
    with urllib.request.urlopen(f"{EUTILS}/esearch.fcgi?{q}", timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))

def esummary(pmid):
    q = urllib.parse.urlencode({"db": "pubmed", "id": pmid, "retmode": "json"})
    with urllib.request.urlopen(f"{EUTILS}/esummary.fcgi?{q}", timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))["result"][pmid]

# 最终映射：编号 -> (PMID, 状态) ；状态 VERIFIED / CANDIDATE / AUTHOR_CONFIRM / GEO_ACCESSION / UNUSED
MAPPING = {
    1: ("26903337", "VERIFIED"),   # Bellani JAMA 2016 LUNG SAFE
    2: ("22797452", "VERIFIED"),   # Ranieri JAMA 2012 Berlin
    3: ("30872586", "VERIFIED"),   # Matthay Nat Rev Dis Primers 2019
    4: ("28792873", "VERIFIED"),   # Thompson NEJM 2017
    5: ("10793167", "VERIFIED"),   # Ware & Matthay NEJM 2000
    6: ("", "AUTHOR_CONFIRM"),     # ATII 线粒体超微结构退变与预后
    7: ("12133657", "VERIFIED"),   # Brealey Lancet 2002 脓毒症线粒体功能障碍
    8: ("", "AUTHOR_CONFIRM"),     # COVID 方向呼吸链亚基转录抑制
    9: ("17178908", "VERIFIED"),   # Szabadkai JCB 2006
    10: ("22850819", "VERIFIED"),  # Rizzuto NRMCB 2012
    11: ("22992592", "VERIFIED"),  # Rowland & Voeltz NRMCB 2012
    12: ("", "AUTHOR_CONFIRM"),    # 凋亡-ARDS
    13: ("", "AUTHOR_CONFIRM"),    # 坏死性凋亡-ARDS
    14: ("40555044", "CANDIDATE"), # 焦亡 ARDS 综述 2025
    15: ("41383584", "CANDIDATE"), # 铁死亡 ALI/ARDS 综述 2025
    16: ("19740426", "CANDIDATE"), # VX-166 caspase 抑制剂脓毒症 2009
    17: ("22195746", "CANDIDATE"), # Duprez Immunity 2011 RIP 坏死 SIRS
    18: ("", "GEO_ACCESSION"),     # GSE165659（预印本，无 PMID）
    19: ("", "UNUSED"),
    20: ("", "UNUSED"),
    21: ("32461654", "VERIFIED"),  # Karczewski Nature 2020 gnomAD
    22: ("", "UNUSED"),
    23: ("32398875", "VERIFIED"),  # Liao Nat Med 2020
    24: ("33657410", "VERIFIED"),  # Ren Cell 2021
    25: ("35027333", "VERIFIED"),  # Baghela EBioMedicine 2022
    26: ("", "UNUSED"),
    27: ("", "GEO_ACCESSION"),     # GSE212865（无 PMID）
    28: ("", "UNUSED"), 29: ("", "UNUSED"), 30: ("", "UNUSED"), 31: ("", "UNUSED"),
    32: ("29195078", "VERIFIED"),  # Subramanian Cell 2017
    33: ("36653562", "VERIFIED"),  # Kurki Nature 2023 FinnGen
    34: ("27543902", "VERIFIED"),  # Szilágyi AJRCMB 2016
    35: ("34475573", "VERIFIED"),  # Võsa Nat Genet 2021
    36: ("37962376", "VERIFIED"),  # JASPAR 2024
    37: ("", "UNUSED"), 38: ("", "UNUSED"), 39: ("", "UNUSED"), 40: ("", "UNUSED"),
    41: ("", "UNUSED"), 42: ("", "UNUSED"), 43: ("", "UNUSED"), 44: ("", "UNUSED"),
    45: ("", "UNUSED"), 46: ("", "UNUSED"), 47: ("", "UNUSED"), 48: ("", "UNUSED"),
    49: ("", "UNUSED"), 50: ("", "UNUSED"), 51: ("", "UNUSED"), 52: ("", "UNUSED"),
    53: ("", "UNUSED"),
    54: ("32416070", "VERIFIED"),  # Blanco-Melo Cell 2020
    55: ("", "AUTHOR_CONFIRM"),    # 脓毒症转录组线粒体方向
    56: ("24185508", "CANDIDATE"), # Singer 2014 线粒体功能障碍脓毒症多器官衰竭
    57: ("", "AUTHOR_CONFIRM"),    # NR1H3/PPARG/EPAS1
    58: ("", "UNUSED"), 59: ("", "UNUSED"),
    60: ("22632970", "VERIFIED"),  # Dixon Cell 2012
    61: ("28985560", "VERIFIED"),  # Stockwell Cell 2017
    62: ("27535533", "VERIFIED"),  # Lek Nature 2016
    63: ("15514116", "VERIFIED"),  # Nemeth Science 2004
    64: ("22461369", "VERIFIED"),  # Dolinay AJRCCM 2012
}

# 最后一搏：给 AUTHOR_CONFIRM 条目补充候选（不写入正式清单，仅登记为候选）
LAST_TRY = {
    6: "alveolar epithelial cell mitochondria ultrastructure ARDS prognosis",
    8: "COVID-19 peripheral blood mitochondrial gene downregulation transcriptome",
    12: "apoptosis lung injury ARDS review epithelial death",
    13: "necroptosis acute respiratory distress syndrome RIPK3",
    55: "sepsis whole blood gene expression mitochondrial downregulation",
    57: "macrophage LXR PPARG hypoxia lipid sepsis",
}
print("=== 最后候选检索 ===", flush=True)
last_cand = {}
for num, term in LAST_TRY.items():
    try:
        res = esearch(term)
        ids = res["esearchresult"]["idlist"]
        cands = []
        for pmid in ids[:3]:
            sm = esummary(pmid)
            cands.append(f"{pmid}: {sm.get('title','')[:90]} ({sm.get('pubdate','')[:4]})")
        last_cand[num] = cands
        for c in cands:
            print(f"[{num}] {c}", flush=True)
        time.sleep(0.4)
    except Exception as e:
        last_cand[num] = [f"ERROR {e}"]
        print(f"[{num}] ERROR {e}", flush=True)

# 取已核实条目完整引文信息
print("=== 拉取正式清单引文详情 ===", flush=True)
rows = []
for num in sorted(MAPPING):
    pmid, status = MAPPING[num]
    if pmid and status in ("VERIFIED", "CANDIDATE"):
        sm = esummary(pmid)
        auth = sm.get("authors", [])[:6]
        astr = ", ".join(a["name"] for a in auth)
        if len(sm.get("authors", [])) > 6:
            astr += ", et al"
        rows.append(dict(ref=num, pmid=pmid, status=status,
                         title=sm.get("title", ""),
                         journal=sm.get("fulljournalname", "") or sm.get("source", ""),
                         year=(sm.get("pubdate", "") or "")[:4],
                         volume=sm.get("volume", ""), issue=sm.get("issue", ""),
                         pages=sm.get("pages", ""),
                         authors_short=astr,
                         doi=next((a["value"] for a in sm.get("articleids", []) if a.get("idtype") == "doi"), "")))
        time.sleep(0.35)
    else:
        rows.append(dict(ref=num, pmid=pmid, status=status, title="", journal="",
                         year="", volume="", issue="", pages="", authors_short="", doi=""))

import pandas as pd
df = pd.DataFrame(rows)
df.to_csv(r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate\P0_reference_final_list.csv",
          index=False, encoding="utf-8-sig")

LIST_CTX = {
    6: "ARDS ATII 线粒体超微结构退变与预后独立相关",
    8: "MT 编码呼吸链亚基转录抑制与 mtROS 增加（COVID 方向）",
    12: "凋亡参与肺泡细胞死亡",
    13: "坏死性凋亡参与 ARDS",
    55: "脓毒症转录组呼吸链/线粒体下调",
    57: "NR1H3/PPARG/EPAS1 髓系驱动候选",
}
report = ["# P0 参考文献 1-64 重建与 PubMed 核实报告（2026-08-21）", "",
          "> 依据：《创新提质方案 v2》§8.3『补齐全部参考文献』。",
          "> 核实接口：NCBI E-utilities（PubMed esearch/esummary）；GEO 数据集条目经 ELink（GDS→PubMed）确认。",
          "> 背景披露：v1 文稿（含 1-64 条原始清单）不在项目内；所有现存文稿版本中 References 均为占位。",
          "> 因此本清单为**按 v4 正文引用上下文重建 + PubMed 逐条核实**的结果；无法唯一锚定的条目标注『待作者确认』，绝不虚构。",
          "", "## 状态分类", "",
          "- **VERIFIED**：检索命中且与引用上下文明确匹配（含作者/标题锚定）。",
          "- **CANDIDATE**：命中与上下文基本匹配，但存在同主题替代文献，作者应最终确认。",
          "- **AUTHOR_CONFIRM**：多轮检索无法唯一锚定到上下文所指文献（v1 原始选择不可恢复），须作者提供原始条目。",
          "- **GEO_ACCESSION**：数据集源文献无 PMID（预印本/未索引），沿用 v4 现行口径以 GEO accession 引用。",
          "- **UNUSED**：该编号在 v4 正文中从未被引用，且 v1 清单不可恢复——建议在正式稿中删除或由作者重新指派。",
          "", "## 核实结果", "",
          "| # | 状态 | PMID | 期刊/年 | 标题 |",
          "|---|---|---|---|---|"]
for _, r in df.iterrows():
    t = (r["title"] or "—")[:95].replace("|", "/")
    j = f"{r['journal']} {r['year']}" if r["journal"] else "—"
    pm = r["pmid"] or "—"
    report.append(f"| {r['ref']} | {r['status']} | {pm} | {j} | {t} |")
report += ["", "## 待作者确认条目的最后候选（未写入正式清单）", ""]
for num, cands in sorted(last_cand.items()):
    report.append(f"- **[{num}]**（{LIST_CTX.get(num, '')}）：")
    for c in cands:
        report.append(f"  - {c}")
report += ["", "## 说明", "",
           "- [14]/[15]/[16]/[17]/[56] 为 CANDIDATE：主题匹配但非唯一解，投稿前请作者按原始意图确认。",
           "- [18]/[27] 源数据集的科学论文无 PubMed 索引（GSE165659 sci-ATAC 为预印本图谱；GSE212865 未索引），按 GEO accession 引用。",
           "- UNUSED 编号（19,20,22,26,28-31,37-53,58,59 共 29 个）：v4 正文未引用；正式英文稿建议压缩编号或由作者补回原条目。",
           "- 本报告与 `_intermediate/P0_reference_final_list.csv` 同步生成；检索过程留痕于 P0_verify_references*.py 与 GEO ELink 输出。"]

out = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\04_AUDIT_GOVERNANCE\P0_References_1to64_Verification_Report.md"
with open(out, "w", encoding="utf-8") as f:
    f.write("\n".join(report))
print("saved", out)
