# -*- coding: utf-8 -*-
"""
P2_ortholog_map.py — Phase 2 前置：IIAMD v1.0 签名 小鼠→人 ortholog 映射
==========================================================================
依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md §4 H6——"仅保留高可信一对一 mouse-human ortholog；
     无法映射基因报告，不以近似同源替代"
方法：mygene.info 在线查询（scopes=symbol, species=mouse, fields=homologene）
     一对一判定（双向）：小鼠符号→恰好 1 个人类同源 且 该人类符号不被多个小鼠符号命中
输出：_intermediate/P2_IIAMD_human_v1.0.csv（mouseSymbol,humanSymbol,direction,log2FC,padj）
     03_LOGS/P2_ortholog_log.txt
"""
import os, sys
import pandas as pd
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import mygene

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
SIG = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "P1_IIAMD_signature_v1.0.csv")
OUT = os.path.join(ROOT, "_intermediate", "P2_IIAMD_human_v1.0.csv")
LOGF = os.path.join(ROOT, "03_LOGS", "P2_ortholog_log.txt")

log = []
def note(m=""):
    log.append(m); print(m, flush=True)

sig = pd.read_csv(SIG)
genes = sorted(set(sig["geneSymbol"]))
note(f"签名基因（小鼠符号）: {len(genes)}")

mg = mygene.MyGeneInfo()
res = mg.querymany(genes, scopes="symbol", species="mouse", fields="homologene",
                   returnall=True, verbose=False)
hit = {r["query"]: r for r in res.get("out", []) if "homologene" in r}
note(f"mygene 返回含 homologene 的查询: {len(hit)}/{len(genes)}")

# 提取人类同源（taxid 9606）
mouse2human = {}
ambiguous = []
for q, r in hit.items():
    homologs = r["homologene"].get("genes", [])
    hs = [g for g in homologs if g[0] == 9606]
    if len(hs) == 1:
        mouse2human[q] = hs[0][1]
    elif len(hs) > 1:
        ambiguous.append(q)
note(f"单人类同源(GeneID): {len(mouse2human)}；多义: {len(ambiguous)}；无: {len(genes)-len(mouse2human)-len(ambiguous)}")

# homologene 给的是 Entrez GeneID —— 批量转 HGNC 符号
ids = sorted(set(mouse2human.values()))
ginfo = mg.getgenes(ids, fields="symbol")
id2sym = {str(g["_id"]): str(g["symbol"]).upper() for g in ginfo if g.get("symbol")}
mouse2human = {m: id2sym[str(i)] for m, i in mouse2human.items() if str(i) in id2sym}
note(f"GeneID→HGNC 符号转换成功: {len(mouse2human)}（无法转符号 {len(ids)-len(id2sym)} 个 ID）")

# 双向一对一：同一人类符号被多个小鼠符号命中时全部剔除
from collections import Counter
cnt = Counter(mouse2human.values())
multi_human = [h for h, c in cnt.items() if c > 1]
one2one = {m: h for m, h in mouse2human.items() if cnt[h] == 1}
note(f"剔除 {sum(cnt[h] for h in multi_human)} 个命中同一人类符号的多小鼠基因 → 一对一 {len(one2one)}")

mapped = sig[sig["geneSymbol"].isin(one2one)].copy()
mapped["humanSymbol"] = mapped["geneSymbol"].map(one2one)
mapped = mapped.rename(columns={"geneSymbol": "mouseSymbol"})
mapped = mapped[["mouseSymbol", "humanSymbol", "direction", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"]]
mapped.to_csv(OUT, index=False)
up = (mapped["direction"] == "up").sum(); dn = (mapped["direction"] == "down").sum()
note(f"输出: {OUT}（up {up} / down {dn}，映射率 {len(mapped)}/{len(sig)} = {len(mapped)/len(sig):.1%}）")

# 锚定节点映射（供 P2 锚定面板用）
ANCHOR = ["Bax","Bak1","Bid","Rictor","Rhoa","Ninj1","Rptor","Mtor","Akt1","Tlr4","Myd88","Tnf","Casp1","Gsdmd","Mlkl","Gpx4","Slc7a11","Acsl4","Nlrp3"]
am = mg.querymany(ANCHOR, scopes="symbol", species="mouse", fields="homologene", returnall=True, verbose=False)
ah = {}
aid2sym = {}
aids = []
for r in am.get("out", []):
    if "homologene" in r:
        hs = [g for g in r["homologene"]["genes"] if g[0] == 9606]
        if len(hs) == 1:
            ah[r["query"]] = str(hs[0][1]); aids.append(str(hs[0][1]))
if aids:
    ag = mg.getgenes(sorted(set(aids)), fields="symbol")
    aid2sym = {g["_id"]: str(g["symbol"]).upper() for g in ag if g.get("symbol")}
ah = {m: aid2sym[i] for m, i in ah.items() if i in aid2sym}
note("锚定节点 human: " + str(ah))
pd.DataFrame([{"mouseSymbol": k, "humanSymbol": v} for k, v in ah.items()]).to_csv(
    os.path.join(ROOT, "_intermediate", "P2_anchor_human.csv"), index=False)

with open(LOGF, "w", encoding="utf-8") as f:
    f.write("\n".join(log))
print("DONE P2 ortholog map")
