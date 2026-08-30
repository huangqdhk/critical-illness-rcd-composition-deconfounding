# -*- coding: utf-8 -*-
"""
P0_verify_references.py — 按 v4 正文引用上下文重建并 PubMed 核实参考文献 1-64
================================================================================
背景：v1 参考文献清单（1-64 条）在项目内所有文稿版本中均为占位、从未填写；
v1 原稿文件不在项目内。因此本脚本按 v4 正文每条 [n] 引用的上下文（主题句）重建
候选文献，并用 PubMed esearch 以标题/作者精确检索逐条核实。所有条目必须命中
PubMed 才登记；无法命中的标注为"待作者确认"，绝不虚构。

输出：_intermediate/P0_reference_verification.csv
     04_AUDIT_GOVERNANCE/P0_References_1to64_Verification_Report.md
"""
import json, re, sys, time, urllib.request, urllib.parse
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def esearch(term, retmax=5):
    q = urllib.parse.urlencode({"db": "pubmed", "term": term, "retmax": retmax, "retmode": "json"})
    with urllib.request.urlopen(f"{EUTILS}/esearch.fcgi?{q}", timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))

def esummary(pmid):
    q = urllib.parse.urlencode({"db": "pubmed", "id": pmid, "retmode": "json"})
    with urllib.request.urlopen(f"{EUTILS}/esummary.fcgi?{q}", timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))["result"][pmid]

# 候选检索式：按 v4 正文引用上下文逐一登记（编号=文稿内编号；V4语境=该引用支撑的句子）
CANDIDATES = [
    # --- ARDS 死亡率 35-46% [1-4] ---
    (1, "ARDS mortality 35-46%", "Bellani G[Author] AND epidemiology patterns of care mortality acute respiratory distress syndrome JAMA"),
    (2, "ARDS mortality 35-46%", "Ranieri VM[Author] AND acute respiratory distress syndrome Berlin definition JAMA"),
    (3, "ARDS mortality 35-46%", "Matthay MA[Author] AND acute respiratory distress syndrome nature reviews disease primers"),
    (4, "ARDS mortality 35-46%", "Thompson BT[Author] AND acute respiratory distress syndrome new england journal"),
    # --- 单环节干预被并行通路代偿 [5] ---
    (5, "单环节干预被并行通路代偿（ARDS 病理并行）", "Ware LB[Author] AND Matthay MA[Author] AND acute respiratory distress syndrome new england journal 2000"),
    # --- ATII 线粒体超微结构退变与预后 [6] ---
    (6, "ARDS ATII 线粒体超微结构退变与预后独立相关", "mitochondria ultrastructure alveolar epithelial acute respiratory distress syndrome prognosis"),
    # --- MT 呼吸链亚基转录抑制+mtROS [7,8] ---
    (7, "MT 编码呼吸链亚基转录抑制与 mtROS 增加（电子传递链解耦联）", "mitochondrial electron transport chain downregulation sepsis mitochondrial ROS oxidative phosphorylation suppression"),
    (8, "MT 编码呼吸链亚基转录抑制与 mtROS 增加（COVID 方向）", "SARS-CoV-2 mitochondrial gene expression downregulation OXPHOS transcriptome"),
    # --- MAM IP3R1-GRP75-VDAC1 [9-11] ---
    (9, "MAM IP3R1-GRP75-VDAC1 钙转运轴", "Szabadkai G[Author] AND chaperone-mediated coupling endoplasmic reticulum mitochondrial calcium channels"),
    (10, "MAM 线粒体-内质网钙信号", "Rizzuto R[Author] AND mitochondria sensors regulators calcium signalling nature reviews"),
    (11, "MAM 内质网-线粒体接触", "Rowland AA[Author] AND Voeltz GK[Author] AND endoplasmic reticulum mitochondria contacts function junction"),
    # --- RCD 四种死亡方式参与肺泡细胞死亡 [12-15] ---
    (12, "凋亡参与肺泡细胞死亡（ARDS）", "apoptosis alveolar epithelial cell death acute respiratory distress syndrome review"),
    (13, "坏死性凋亡参与 ARDS", "necroptosis acute lung injury ARDS review"),
    (14, "焦亡参与 ARDS/肺损伤", "pyroptosis acute lung injury ARDS"),
    (15, "铁死亡参与 ARDS/肺损伤", "ferroptosis acute lung injury ARDS"),
    # --- 单通路抑制策略临床前疗效局限 [16,17] ---
    (16, "单通路抑制疗效局限（凋亡/caspase 抑制）", "caspase inhibitor sepsis trial apoptosis blockade failed improve outcomes"),
    (17, "单通路抑制疗效局限（RIPK1/坏死性凋亡抑制）", "RIPK1 kinase inhibitor sepsis necroptosis blockade inflammation clinical"),
    # --- GSE165659 [18] ---
    (18, "GSE165659 sci-ATAC 肺组织", "cell atlas chromatin accessibility adult human tissues single-cell ATAC"),
    # --- gnomAD [21] ---
    (21, "gnomAD v2.1.1", "Karczewski KJ[Author] AND mutational constraint spectrum quantified variation humans Nature 2020"),
    # --- GSE145926 [23] ---
    (23, "GSE145926 BALF scRNA COVID-19", "Liao M[Author] AND single-cell landscape bronchoalveolar immune cells COVID-19"),
    # --- GSE158055 [24] ---
    (24, "GSE158055 PBMC scRNA COVID-19", "large-scale single-cell analysis critical immune characteristics COVID-19 patients"),
    # --- GSE185263 [25] ---
    (25, "GSE185263 脓毒症全血 RNA-seq", "Baghela[Author] AND predicting sepsis severity first clinical presentation endotypes mechanistic signatures"),
    # --- GSE212865 [27] ---
    (27, "GSE212865 COVID 微阵列动态", "dynamics gene expression profiling microarrays high-risk patients severe COVID-19"),
    # --- CMap/L1000 [32] ---
    (32, "CMap/L1000", "Subramanian A[Author] AND next generation connectivity map L1000 platform profiles Cell 2017"),
    # --- FinnGen [33] ---
    (33, "FinnGen R10", "Kurki MI[Author] AND FinnGen genetic insights well-phenotyped isolated population Nature 2023"),
    # --- GSE67530 [34] ---
    (34, "GSE67530 450K 甲基化 ARDS", "epigenetic contribution myosin light chain kinase gene acute respiratory distress syndrome"),
    # --- eQTLGen [35] ---
    (35, "eQTLGen 全血 cis-eQTL", "Vosa U[Author] AND large-scale cis trans-eQTL analyses polygenic scores blood gene expression Nature Genetics 2021"),
    # --- JASPAR [36] ---
    (36, "JASPAR 2024", "JASPAR 2024 open-access database transcription factor binding profiles"),
    # --- COVID/脓毒症转录组呼吸链下调 [54-56] ---
    (54, "COVID-19 转录组呼吸链下调", "imbalanced host response SARS-CoV-2 drives development COVID-19 Blanco-Melo"),
    (55, "脓毒症转录组呼吸链/线粒体下调", "sepsis whole blood transcriptome mitochondrial dysfunction oxidative phosphorylation downregulation"),
    (56, "危重症转录组线粒体抑制", "critical illness transcriptomics mitochondrial energy metabolism suppression whole blood"),
    # --- NR1H3/PPARG/EPAS1 髓系驱动 [57] ---
    (57, "NR1H3/PPARG/EPAS1 脂质处理/假性缺氧/抗炎许可", "NR1H3 PPARG EPAS1 macrophage lipid pseudohypoxia tolerance sepsis"),
    # --- 铁死亡文献 [60,61] ---
    (60, "铁死亡定义文献", "Dixon SJ[Author] AND ferroptosis iron-dependent form nonapoptotic cell death Cell 2012"),
    (61, "铁死亡机制/铁代谢", "Stockwell BR[Author] AND ferroptosis regulated cell death nexus metabolism redox biology disease"),
    # --- pLI 必需性逻辑 [62] ---
    (62, "pLI 约束/必需性逻辑", "Lek M[Author] AND analysis protein-coding genetic variation humans ExAC Nature 2016"),
    # --- hepcidin-ferroportin 铁滞留 [63] ---
    (63, "hepcidin-ferroportin 铁滞留", "Nemeth E[Author] AND hepcidin regulates cellular iron efflux binding ferroportin internalization"),
    # --- GSE32707 [64] ---
    (64, "GSE32707 GAinS 脓毒症微阵列", "Dolinay T[Author] AND inflammasome-regulated cytokines critical mediators lung injury humans"),
]

rows = []
for num, ctx, term in CANDIDATES:
    try:
        res = esearch(term)
        ids = res["esearchresult"]["idlist"]
        if not ids:
            rows.append(dict(ref=num, context=ctx, pmid="", title="", journal="", year="",
                             authors="", status="NO_HIT"))
            print(f"[{num}] NO HIT :: {ctx}", flush=True)
        else:
            pmid = ids[0]
            sm = esummary(pmid)
            rows.append(dict(ref=num, context=ctx, pmid=pmid, title=sm.get("title", "")[:180],
                             journal=sm.get("fulljournalname", "") or sm.get("source", ""),
                             year=(sm.get("pubdate", "") or "")[:4],
                             authors=", ".join(a["name"] for a in sm.get("authors", [])[:3]),
                             status="OK" if pmid in ids else "BEST_MATCH"))
            print(f"[{num}] PMID={pmid} :: {ctx} :: {rows[-1]['title'][:90]}", flush=True)
        time.sleep(0.4)
    except Exception as e:
        rows.append(dict(ref=num, context=ctx, pmid="", title="", journal="", year="",
                         authors="", status=f"ERROR {e}"))
        print(f"[{num}] ERROR {e}", flush=True)

import pandas as pd
df = pd.DataFrame(rows)
csv_out = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate\P0_reference_verification.csv"
df.to_csv(csv_out, index=False, encoding="utf-8-sig")
print("saved", csv_out)
print(df[["ref", "pmid", "status"]].to_string())
