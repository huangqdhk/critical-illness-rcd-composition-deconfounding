# -*- coding: utf-8 -*-
"""
M4_HPA_localization_check.py — 80 基因亚细胞定位核验（方案 M4 步骤 4，指南 §5.3 干侧）。
口径：UniProt reviewed（Swiss-Prot）亚细胞定位注释（HPA 数据汇入口径），
判定 MAM/线粒体/内质网等定位类别，作"框架合法性补充"（不参与 M4 主判定）。
输出：01_FIGURE_DATA_CSV/Main/Figure_10E.csv；02_SUPPLEMENTARY_TABLES/.../Table_S59e_M4_Subcellular_Annotation.csv
"""
import io, os, sys, urllib.request, urllib.parse, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
import pandas as pd

BASE = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
man = pd.read_csv(BASE + r"\04_AUDIT_GOVERNANCE\Mitoxyperilysis_Gene_Manifest_v1.0.csv")
genes = man["gene_symbol"].tolist()

def fetch_exact(cands):
    query = " OR ".join(f"gene_exact:{s}" for s in cands) + " AND reviewed:true"
    fields = "accession,gene_primary,cc_subcellular_location,protein_name"
    url = ("https://rest.uniprot.org/uniprotkb/stream?format=tsv&query="
           + urllib.parse.quote(query) + "&fields=" + fields)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read().decode("utf-8")

def candidates(g):
    sym = g["gene_symbol"]
    cands = [sym]
    al = str(g.get("aliases", "") or "")
    cands += [a.strip() for a in al.split(",") if a.strip() and a.strip() != "nan"]
    if sym.startswith("MT-"):
        cands.append(sym[3:])
    return list(dict.fromkeys(cands))

cache = BASE + r"\_intermediate\uniprot_80genes.tsv"
rows_all = []
for i, (_, g) in enumerate(man.iterrows()):
    sym = g["gene_symbol"]
    hit = None
    for attempt in range(3):
        try:
            txt = fetch_exact(candidates(g))
            hit = pd.read_csv(io.StringIO(txt), sep="\t") if txt.strip() else pd.DataFrame()
            break
        except Exception:
            if attempt == 2:
                raise
            time.sleep(5)
    # 取 gene_primary 在候选集内的第一行
    if len(hit):
        hit["ok"] = hit["Gene Names (primary)"].isin(candidates(g))
        hit = hit[hit["ok"]]
    if len(hit):
        rows_all.append({"gene": sym, "uniprot": hit["Entry"].iloc[0],
                         "loc": str(hit["Subcellular location [CC]"].iloc[0])})
    else:
        rows_all.append({"gene": sym, "uniprot": "", "loc": ""})
    if (i + 1) % 20 == 0:
        print(f"  {i+1}/80 genes done")
cache_df = pd.DataFrame(rows_all)
cache_df.to_csv(cache, index=False)
print("UniProt 逐基因查询完成:", (cache_df["uniprot"] != "").sum(), "/80 命中")

def classify(loc):
    if not isinstance(loc, str):
        return "未注释"
    l = loc.lower()
    tags = []
    if "mitochondrion" in l:
        tags.append("线粒体")
    if "endoplasmic reticulum" in l:
        tags.append("内质网")
    if "cell membrane" in l or "membrane" in l and "endoplasmic" not in l:
        tags.append("膜")
    if "cytoplasm" in l or "cytosol" in l:
        tags.append("胞质")
    if "nucleus" in l:
        tags.append("核")
    if "secreted" in l:
        tags.append("分泌")
    return "+".join(tags) if tags else "其他"

rows = []
for _, g in man.iterrows():
    sym = g["gene_symbol"]
    hit = cache_df[cache_df["gene"] == sym]
    if len(hit) == 0 or not hit["uniprot"].iloc[0]:
        rows.append({"gene": sym, "arm": g["arm"], "uniprot": "", "subcellular": "未找到", "is_mito_er": ""})
        continue
    loc = hit["loc"].iloc[0]
    cls = classify(loc)
    mito = "线粒体" in cls
    er = "内质网" in cls
    rows.append({"gene": sym, "arm": g["arm"],
                 "uniprot": hit["uniprot"].iloc[0],
                 "subcellular": cls,
                 "is_mito_er": ("是" if (mito or er) else "否"),
                 "subcellular_raw": loc[:150]})
out = pd.DataFrame(rows)
print("\n模块级统计（含线粒体/内质网注释的比例）:")
stat = out.assign(has_mito_er=out["is_mito_er"] == "是").groupby("arm")["has_mito_er"].agg(["sum", "count"])
stat["pct"] = (stat["sum"] / stat["count"] * 100).round(0)
print(stat.to_string())
print("\n上游臂未注释或非线粒体/内质网的基因:")
print(out[(out["arm"] == "upstream_collapse") & (out["is_mito_er"] != "是")][["gene", "subcellular"]].to_string(index=False))

# 输出（Figure_10E + Table_S59e）
out["gene_set_version"] = "Mitoxy-80_v1.0"
out["score_version"] = "m4hpa_v1.0"
out.to_csv(BASE + r"\01_FIGURE_DATA_CSV\Main\Figure_10E.csv", index=False)
out.to_csv(BASE + r"\02_SUPPLEMENTARY_TABLES\SUPPLEMENTARY_Tables_CSV\Table_S59e_M4_Subcellular_Annotation.csv", index=False)
print("\n已写 Figure_10E.csv / Table_S59e_M4_Subcellular_Annotation.csv")
