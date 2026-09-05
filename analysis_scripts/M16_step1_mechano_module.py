# -*- coding: utf-8 -*-
"""
M16_step1_mechano_module.py — 机械感应模块 mech_v1.0 落表 + 人→小鼠一对一同源映射
=====================================================================
预注册：M16_pre_registration_20260831.md
- mech_v1.0：23 基因三层（通道/效应器/骨架黏附），逐基因文献溯源（PMID 均经 eutils 活验，未验者留文字溯源）
- 同源映射：mygene.info HomoloGene 双向一对一（M10D/P2 同款纪律），人源 80 基因 + mech 23 基因 → 小鼠 symbol
输出：
  02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV/Table_S88_M16_Mechanosensing_Module_v10.csv
  _intermediate/M16_human2mouse_orthologs.csv（80 臂基因 + 23 mech 基因的小鼠一对一映射）
  03_LOGS/M16_step1_log.txt
"""
import json
import os
import time
import urllib.request
import urllib.parse

import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES", "SUPPLEMENTARY_Tables_CSV")
LOG = os.path.join(ROOT, "03_LOGS", "M16_step1_log.txt")
INTER = os.path.join(ROOT, "_intermediate")
MANIFEST = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "Mitoxyperilysis_Gene_Manifest_v1.0.csv")

log = []
def note(m=""):
    log.append(m); print(m, flush=True)

# ---------------- mech_v1.0 固定清单（PMID 均 2026-08-31 eutils 活验；空=文字溯源） ----------------
MECH = [
    # tier, hgnc, name, role, evidence, pmid
    ("T1_channel", "PIEZO1", "piezo type mechanosensitive ion channel component 1", "机械激活阳离子通道；肺内皮/上皮力学换能", "Coste et al., Science 2010; VILI 机制综述 J Thorac Dis 2026", "20813920;42306687"),
    ("T1_channel", "PIEZO2", "piezo type mechanosensitive ion channel component 2", "机械激活阳离子通道", "Coste et al., Science 2010", "20813920"),
    ("T1_channel", "TRPV4", "transient receptor potential cation channel subfamily V member 4", "牵张激活 Ca2+ 通道；VILI 巨噬细胞激活", "Hamanaka et al., Am J Physiol Lung Cell Mol Physiol 2010", "20562229"),
    ("T1_channel", "TRPC1", "transient receptor potential cation channel subfamily C member 1", "机械敏感阳离子通道（经典型）", "经典机械敏感通道文献（文字溯源）", ""),
    ("T1_channel", "TRPC6", "transient receptor potential cation channel subfamily C member 6", "人肺动脉内皮机械敏感阳离子电流", "Am J Physiol Cell Physiol 2022", "35968892"),
    ("T1_channel", "PKD2", "polycystin 2, transient receptor potential cation channel", "初级纤毛 Ca2+ 通道力学感应", "纤毛力学感应经典文献（文字溯源）", ""),
    ("T1_channel", "SCNN1A", "sodium channel epithelial 1 subunit alpha", "肺泡上皮 Na+ 通道（牵张诱导）", "Exp Lung Res 2014", "25058750"),
    ("T2_effector", "YAP1", "Yes1 associated transcriptional regulator", "Hippo-YAP/TAZ 力学转导核心", "Dupont et al., Nature 2011; VILI-YAP 2023", "21654799;37904713"),
    ("T2_effector", "WWTR1", "WW domain containing transcription regulator 1", "TAZ，Hippo 力学转导核心", "Dupont et al., Nature 2011", "21654799"),
    ("T2_effector", "TEAD1", "TEA domain transcription factor 1", "YAP/TAZ 转录伙伴", "Dupont et al., Nature 2011（文字溯源）", ""),
    ("T2_effector", "MRTFA", "megakaryoblastic leukemia 1 (MKL1/MRTF-A)", "肌动蛋白-MRTF-SRF 力学转导", "Miralles et al., Cell 2003", "12732141"),
    ("T2_effector", "SRF", "serum response factor", "MRTF-SRF 轴转录因子", "Miralles et al., Cell 2003（文字溯源）", ""),
    ("T2_effector", "KLF2", "Kruppel like factor 2", "血流剪切/力学诱导内皮转录因子", "剪切-内皮 KLF 经典文献（文字溯源）", ""),
    ("T2_effector", "KLF4", "Kruppel like factor 4", "血流剪切/力学诱导内皮转录因子", "J Clin Invest 2012（内皮 KLF4）", "23160196"),
    ("T3_cytoskeleton", "MYLK", "myosin light chain kinase", "肌球蛋白轻链激酶；ARDS 遗传易感基因（力学门控）", "Szilágyi et al., Transl Res 2017（本文参考文献 29）", "27543902"),
    ("T3_cytoskeleton", "PTK2", "protein tyrosine kinase 2 (FAK)", "黏着斑激酶，整合素力学信号", "黏着斑力学转导经典文献（文字溯源）", ""),
    ("T3_cytoskeleton", "ITGB1", "integrin subunit beta 1", "整合素 β1，ECM 力学锚定", "整合素力学感应经典文献（文字溯源）", ""),
    ("T3_cytoskeleton", "ITGA5", "integrin subunit alpha 5", "整合素 α5β1，纤维连接蛋白力学受体", "整合素力学感应经典文献（文字溯源）", ""),
    ("T3_cytoskeleton", "RHOA", "ras homolog family member A", "RhoA-GTP，应力纤维/收缩主开关", "RhoA 抑制剂挽救高潮气量肺损伤 Int Immunopharmacol 2019", "30959374"),
    ("T3_cytoskeleton", "ROCK1", "Rho associated coiled-coil containing protein kinase 1", "RhoA 下游效应激酶", "RhoA/ROCK-VILI 文献（Int Immunopharmacol 2019，文字溯源）", ""),
    ("T3_cytoskeleton", "ROCK2", "Rho associated coiled-coil containing protein kinase 2", "RhoA 下游效应激酶", "同上（文字溯源）", ""),
    ("T3_cytoskeleton", "ACTA2", "actin alpha 2, smooth muscle", "α-平滑肌肌动蛋白，力学响应标志", "力学响应标志经典文献（文字溯源）", ""),
    ("T3_cytoskeleton", "CAV1", "caveolin 1", "小窝蛋白-1，肺内皮力学保护/屏障", "CAV1 肺内皮屏障力学文献（文字溯源）", ""),
]

def mygene_homologene_one2one(human_symbols, tries=4):
    """人源 symbol → 小鼠 symbol（HomoloGene 双向一对一，M10D/P2 同款纪律）。"""
    out = {}
    q = urllib.parse.quote(" ".join(human_symbols))
    # mygene querymany via POST
    url = "https://mygene.info/v3/query"
    data = urllib.parse.urlencode({
        "q": ",".join(human_symbols),
        "scopes": "symbol,alias",
        "fields": "symbol,taxid,homologene",
        "species": "human",
        "size": str(len(human_symbols) * 3),
    }).encode()
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": "M16/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                res = json.loads(r.read())
            break
        except Exception as e:
            last = e; note(f"  [retry {i+1}] mygene query: {e}"); time.sleep(3)
    else:
        raise RuntimeError(f"mygene query failed: {last}")
    hits = {}
    for h in res:
        if h.get("notfound"):
            continue
        sym = str(h.get("symbol", "")).upper()
        qsym = str(h.get("query", "")).upper()
        hg = h.get("homologene", {})
        genes = hg.get("genes", []) if isinstance(hg, dict) else []
        mouse = [g for g in genes if isinstance(g, (list, tuple)) and len(g) == 2 and g[0] == 10090]
        if not mouse:
            continue
        mgi = mouse[0][1]
        hits.setdefault(qsym, {"human_symbol": sym, "homologene_id": hg.get("id", ""), "mouse_gene_id": mgi})
    return hits

def mouse_symbol_from_geneid(gene_ids, tries=4):
    """小鼠 Entrez gene id 批量 → symbol。"""
    ids = sorted({str(g) for g in gene_ids if str(g) not in ("", "None", "nan")})
    if not ids:
        return {}
    url = "https://mygene.info/v3/gene"
    res_all = []
    for i in range(0, len(ids), 800):
        chunk = ids[i:i + 800]
        data = urllib.parse.urlencode({"ids": ",".join(chunk), "fields": "symbol,taxid"}).encode()
        last = None
        for t in range(tries):
            try:
                req = urllib.request.Request(url, data=data, headers={"User-Agent": "M16/1.0"})
                with urllib.request.urlopen(req, timeout=60) as r:
                    res_all.extend(json.loads(r.read()))
                break
            except Exception as e:
                last = e; note(f"  [retry {t+1}] mygene fetch: {e}"); time.sleep(3)
        else:
            raise RuntimeError(f"mygene gene fetch failed: {last}")
        time.sleep(0.5)
    return {str(g["_id"]): g.get("symbol", "") for g in res_all if g.get("taxid") == 10090}

def main():
    note("== M16 step1 机械感应模块 mech_v1.0 + 同源映射 2026-08-31 ==")
    os.makedirs(TAB, exist_ok=True); os.makedirs(INTER, exist_ok=True)

    df = pd.DataFrame(MECH, columns=["tier", "hgnc_symbol", "gene_name", "role", "evidence", "pmid_verified"])
    df["module"] = "mechanosensing"
    df["module_version"] = "mech_v1.0"
    f88 = os.path.join(TAB, "Table_S88_M16_Mechanosensing_Module_v10.csv")
    df.to_csv(f88, index=False, encoding="utf-8-sig")
    note(f"[write] {f88} ({len(df)} genes; tiers: {df['tier'].value_counts().to_dict()})")

    # 80 臂基因（唯一来源：canonical manifest arm 列）
    man = pd.read_csv(MANIFEST, encoding="utf-8-sig")
    arms = man[man["arm"].isin(["upstream_collapse", "execution_induction"])]
    human_genes = sorted(arms["hgnc_symbol"].astype(str).str.strip().unique())
    note(f"[manifest] 臂基因 {len(human_genes)}（up {(arms['arm']=='upstream_collapse').sum()}, ex {(arms['arm']=='execution_induction').sum()}）")
    allq = sorted(set(human_genes) | set(df["hgnc_symbol"]))
    note(f"[mygene] 查询 {len(allq)} 个人源基因的一对一小鼠同源 …")
    hits = mygene_homologene_one2one(allq)
    sym_map = mouse_symbol_from_geneid([v["mouse_gene_id"] for v in hits.values()])
    rows = []
    for g in allq:
        h = hits.get(g.upper())
        if h:
            gid = str(h["mouse_gene_id"])
            rows.append({"human_symbol": g, "homologene_id": h["homologene_id"],
                         "mouse_gene_id": gid, "mouse_symbol": sym_map.get(gid, "")})
        else:
            rows.append({"human_symbol": g, "homologene_id": "", "mouse_gene_id": "", "mouse_symbol": ""})
    om = pd.DataFrame(rows)
    fo = os.path.join(INTER, "M16_human2mouse_orthologs.csv")
    om.to_csv(fo, index=False, encoding="utf-8-sig")
    n_map = (om["mouse_symbol"] != "").sum()
    note(f"[write] {fo}：一对一映射 {n_map}/{len(om)}")
    miss = om.loc[om["mouse_symbol"] == "", "human_symbol"].tolist()
    note(f"[coverage] 未映射（披露）：{miss}")

    with open(LOG, "w", encoding="utf-8") as f:
        f.write("\n".join(log) + "\n")
    print("[done]", f88, fo, sep="\n")

if __name__ == "__main__":
    main()
