# -*- coding: utf-8 -*-
"""
P0_rerun_novelty_search.py — §8.3 硬阻断项：投稿前重跑查新检索（2026-08-21）
==================================================================================
依据：《111黄裕荣创新提质_v2.md》§8.3 与 P0_FROZEN_ANALYSIS_PLAN_v1.0.md
任务：按预注册检索式在 PubMed 重跑系统查新并存档检索日志。
检索式（与 P0_Novelty_Search_Table_v1.0.md 基线一致）：
  1) "mitoxyperilysis"[tiab] OR "mitoxyperiosis"[tiab]           —— 术语总体
  2) ("mitoxyperilysis"[tiab] OR "mitoxyperiosis"[tiab]) AND
     ("acute respiratory distress"[tiab] OR "ARDS"[tiab] OR lung[tiab]) —— 与 ARDS/肺交叉
  3) "mitoxyperilysis"[tiab] OR "mitoxyperiosis"[tiab] AND sepsis[tiab] —— 与脓毒症交叉
输出：_intermediate/pubmed_novelty_rerun_20260821.json（原始响应）
      04_AUDIT_GOVERNANCE/P0_Novelty_Search_Rerun_Log_20260821.md（检索日志）
不修改 P0_Novelty_Search_Table_v1.0.md（该表为冻结基线；重跑结论记录于新日志并追加）
"""
import json, os, time, urllib.request, urllib.parse
import xml.etree.ElementTree as ET

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
INTER = os.path.join(ROOT, "_intermediate")
OUT_LOG = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "P0_Novelty_Search_Rerun_Log_20260821.md")
JSON_OUT = os.path.join(INTER, "pubmed_novelty_rerun_20260821.json")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def esearch(term, retmax=50):
    q = urllib.parse.urlencode({"db": "pubmed", "term": term, "retmax": retmax, "retmode": "json"})
    with urllib.request.urlopen(f"{EUTILS}/esearch.fcgi?{q}", timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))

def esummary(pmids):
    q = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(pmids), "retmode": "json"})
    with urllib.request.urlopen(f"{EUTILS}/esummary.fcgi?{q}", timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))

def efetch_abstract(pmids):
    """Batch fetch abstracts (XML) for a list of PMIDs."""
    q = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"})
    with urllib.request.urlopen(f"{EUTILS}/efetch.fcgi?{q}", timeout=120) as r:
        return r.read().decode("utf-8")

QUERIES = {
    "Q1_term": '"mitoxyperilysis"[tiab] OR "mitoxyperiosis"[tiab]',
    "Q2_ARDS_lung": '("mitoxyperilysis"[tiab] OR "mitoxyperiosis"[tiab]) AND ("acute respiratory distress"[tiab] OR ARDS[tiab] OR lung[tiab])',
    "Q3_sepsis": '("mitoxyperilysis"[tiab] OR "mitoxyperiosis"[tiab]) AND sepsis[tiab]',
}

log = [f"# 查新检索重跑日志（2026-08-21）",
       "",
       f"> 检索时间：{time.strftime('%Y-%m-%d %H:%M:%S')}（本地时区）",
       "> 检索人：分析作者 + AI 辅助检索；检索接口：NCBI E-utilities（PubMed）",
       "> 依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md 与《创新提质方案 v2》§8.3『重跑查新检索并更新查新表』",
       "> 基线：`P0_Novelty_Search_Table_v1.0.md`（2026-08-20 版，Mitoxyperilysis×ARDS=0 条）",
       ""]

all_results = {}
for qk, qv in QUERIES.items():
    try:
        res = esearch(qv)
        ids = res["esearchresult"]["idlist"]
        count = res["esearchresult"]["count"]
        all_results[qk] = {"term": qv, "count": count, "ids": ids}
        log.append(f"## {qk}")
        log.append("")
        log.append(f"- 检索式：`{qv}`")
        log.append(f"- 命中：**{count}** 条")
        log.append(f"- PMID 列表：{', '.join(ids) if ids else '（无）'}")
        log.append("")
        time.sleep(0.4)
    except Exception as e:
        all_results[qk] = {"term": qv, "error": str(e)}
        log.append(f"## {qk}")
        log.append("")
        log.append(f"- 检索式：`{qv}`")
        log.append(f"- **检索失败**：{e}")
        log.append("")

# Q1 全量明细（标题/期刊/年份）
q1_ids = all_results.get("Q1_term", {}).get("ids", [])
if q1_ids:
    try:
        sm = esummary(q1_ids)
        log.append("## Q1 逐条明细（PubMed esummary）")
        log.append("")
        log.append("| PMID | 标题 | 期刊 | 年 | DOI |")
        log.append("|---|---|---|---|---|")
        for pid in q1_ids:
            it = sm["result"].get(pid, {})
            title = (it.get("title") or "")[:140].replace("|", "/")
            jour = it.get("fulljournalname", "") or it.get("source", "")
            year = (it.get("pubdate") or "")[:4]
            doi = "; ".join(it.get("articleids", []) and [a.get("value") for a in it.get("articleids", []) if a.get("idtype") == "doi"] or [])
            log.append(f"| {pid} | {title} | {jour} | {year} | {doi} |")
        log.append("")
        time.sleep(0.4)
    except Exception as e:
        log.append(f"- esummary 失败：{e}")
        log.append("")

# Q2/Q3 交叉：逐条取标题确认是否真为 ARDS/肺/脓毒症语境
for qk in ("Q2_ARDS_lung", "Q3_sepsis"):
    ids = all_results.get(qk, {}).get("ids", [])
    if ids:
        try:
            sm = esummary(ids)
            log.append(f"## {qk} 逐条明细（PubMed esummary）")
            log.append("")
            log.append("| PMID | 标题 | 期刊 | 年 |")
            log.append("|---|---|---|---|")
            for pid in ids:
                it = sm["result"].get(pid, {})
                title = (it.get("title") or "")[:160].replace("|", "/")
                jour = it.get("fulljournalname", "") or it.get("source", "")
                year = (it.get("pubdate") or "")[:4]
                log.append(f"| {pid} | {title} | {jour} | {year} |")
            log.append("")
            time.sleep(0.4)
        except Exception as e:
            log.append(f"- esummary 失败：{e}")
            log.append("")

log.append("## 结论判定（对照基线表）")
log.append("")
q2_count = all_results.get("Q2_ARDS_lung", {}).get("count", "?")
q3_count = all_results.get("Q3_sepsis", {}).get("count", "?")
log.append(f"- 术语总命中（Q1）：{all_results.get('Q1_term', {}).get('count', '?')} 条（基线 2026-08-20 为 10 条 Title/Abstract 记录口径）")
log.append(f"- Mitoxyperilysis×ARDS/肺（Q2）：{q2_count} 条")
log.append(f"- Mitoxyperilysis×脓毒症（Q3）：{q3_count} 条")
log.append("")
log.append("> 判读原则：命中条数 >0 不等于『直接实验证据』；须逐条核对文章类型（原始研究/评论/综述）与疾病语境后再更新查新表。")
log.append("> 本日志为检索过程留痕；查新表主文件更新另见 `P0_Novelty_Search_Table_v1.1.md`（如新增记录）。")

with open(JSON_OUT, "w", encoding="utf-8") as f:
    json.dump(all_results, f, ensure_ascii=False, indent=2)
with open(OUT_LOG, "w", encoding="utf-8") as f:
    f.write("\n".join(log))
print(f"查新检索日志已写入：{OUT_LOG}")
print(f"Q1 count={all_results.get('Q1_term', {}).get('count')}, Q2 count={q2_count}, Q3 count={q3_count}")
print("DONE P0 novelty re-run")
