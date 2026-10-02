# -*- coding: utf-8 -*-
"""
M4_HPA_localization_check.py — 80 基因亚细胞定位核验（方案 M4 步骤 4，指南 §5.3 干侧）。

口径：UniProt **reviewed（Swiss-Prot）人源（organism_id:9606）** 亚细胞定位注释
（HPA 数据汇入口径），判定 MAM/线粒体/内质网等定位类别，作"框架合法性补充"
（不参与 M4 主判定）。

⚠ 2026-09-20 修复（M2，本脚本为唯一根因所在）
------------------------------------------------------------------
本脚本原查询只有 ``AND reviewed:true``，**缺 organism_id:9606**；而 ``:49-54``
又"取基因名命中候选集的第一行"，物种由 UniProt 返回顺序决定 → 80 行里
**32 行**拿到的是牛/大鼠/小鼠等非人源 accession（如 CALR=P15253 牛、
CANX=P24643 牛、HSPA9=O35501 大鼠、SIGMAR1=Q60492 小鼠、VDAC1=A0A6P7EFR0）。
预注册 M4_pre_registration:15/:36 写明"UniProt reviewed 注释（**HPA 汇入口径**）"
——HPA 是人源专属库，故这不是"没想到要限定物种"，而是**代码没实现计划里的口径**，
属计划外实现偏差（M2 原复核记录见 归档\07_审计报告备份\可视化文档_合并前_20260928/README_审稿意见.md §1 M2）。

修法两处：
  ① 查询加 ``AND organism_id:9606``（口径落地）；
  ② 命中多行时不再按返回顺序取第一行，而是**优先 gene_primary 与基因符号完全
     相等**的那行（返回顺序取决于 UniProt 内部排序，不可作为口径）。

输出：01_FIGURE_DATA_CSV/Main/**Figure_5J**.csv（4J 读的就是这个文件名）
      02_SUPPLEMENTARY_TABLES/.../Table_S59e_M4_Subcellular_Annotation.csv
"""
import io, os, sys, urllib.request, urllib.parse, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
import pandas as pd

BASE = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"

#: 人源口径（HPA 汇入口径）。写成常量而不是散在查询串里，便于 lint 断言。
ORGANISM_ID = 9606

man = pd.read_csv(BASE + r"\04_AUDIT_GOVERNANCE\Mitoxyperilysis_Gene_Manifest_v1.0.csv")
genes = man["gene_symbol"].tolist()

def fetch_exact(cands):
    query = ("(" + " OR ".join(f"gene_exact:{s}" for s in cands) + ")"
             + " AND reviewed:true AND organism_id:%d" % ORGANISM_ID)
    fields = "accession,gene_primary,organism_id,cc_subcellular_location,protein_name"
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

def pick_row(hit, cands, sym):
    """候选集过滤 + 确定性取行：gene_primary == 基因符号者优先。

    原实现是 ``hit.iloc[0]``，也就是"UniProt 谁排前面就取谁"。加了物种限定
    之后同一基因仍可能有多条人源条目，故显式给出优先级，不再依赖返回顺序。
    """
    hit = hit[hit["Gene Names (primary)"].isin(cands)].copy()
    if not len(hit):
        return None
    hit["_exact"] = (hit["Gene Names (primary)"] == sym).astype(int)
    hit = hit.sort_values("_exact", ascending=False, kind="stable")
    return hit.iloc[0]

#: 缓存文件名带 _human 后缀：旧缓存 ``uniprot_80genes.tsv`` 装的是未限物种的
#: 结果，沿用旧名会静默复用错误数据（正是本次 M2 漏网的机制之一）。
cache = BASE + r"\_intermediate\uniprot_80genes_human.tsv"
#: 复核期冻结的人源查询结果（_m2_verify/uniprot_80genes_human.tsv），
#: 可离线复现；schema 与本脚本不同，需改名后落成正式缓存。
FROZEN = BASE + r"\_m2_verify\uniprot_80genes_human.tsv"

rows_all = []
if os.path.exists(cache):
    print("[cache] 复用", cache)
    cache_df = pd.read_csv(cache)
elif os.path.exists(FROZEN):
    print("[cache] 由冻结结果落盘（复核期查询，离线可复现）", FROZEN)
    fz = pd.read_csv(FROZEN)
    cache_df = pd.DataFrame({"gene": fz["gene"],
                             "uniprot": fz["uniprot_human"].fillna(""),
                             "loc": fz["loc_human"].fillna(""),
                             "organism": fz.get("organism", "")})
    cache_df.to_csv(cache, index=False)
else:
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
        row = pick_row(hit, candidates(g), sym) if len(hit) else None
        if row is not None:
            rows_all.append({"gene": sym, "uniprot": row["Entry"],
                             "loc": str(row["Subcellular location [CC]"]),
                             "organism": str(row.get("Organism", ""))})
        else:
            rows_all.append({"gene": sym, "uniprot": "", "loc": "",
                             "organism": ""})
        if (i + 1) % 20 == 0:
            print(f"  {i+1}/80 genes done")
    cache_df = pd.DataFrame(rows_all)
    cache_df.to_csv(cache, index=False)
print("UniProt 逐基因查询完成（organism_id=%d）: %d/80 命中"
      % (ORGANISM_ID, int((cache_df["uniprot"].fillna("") != "").sum())))

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
    if len(hit) == 0 or not str(hit["uniprot"].fillna("").iloc[0]):
        rows.append({"gene": sym, "arm": g["arm"], "uniprot": "", "subcellular": "未找到", "is_mito_er": ""})
        continue
    loc = hit["loc"].fillna("").iloc[0]
    cls = classify(loc)
    mito = "线粒体" in cls
    er = "内质网" in cls
    rows.append({"gene": sym, "arm": g["arm"],
                 "uniprot": str(hit["uniprot"].iloc[0]),
                 "subcellular": cls,
                 "is_mito_er": ("是" if (mito or er) else "否"),
                 "subcellular_raw": str(loc)[:150]})
out = pd.DataFrame(rows)
print("\n模块级统计（含线粒体/内质网注释的比例）:")
stat = out.assign(has_mito_er=out["is_mito_er"] == "是").groupby("arm")["has_mito_er"].agg(["sum", "count"])
stat["pct"] = (stat["sum"] / stat["count"] * 100).round(0)
print(stat.to_string())
print("\n上游臂未注释或非线粒体/内质网的基因:")
print(out[(out["arm"] == "upstream_collapse") & (out["is_mito_er"] != "是")][["gene", "subcellular"]].to_string(index=False))

# 输出（Figure_5J + Table_S59e）
# ⚠ 原脚本写 Figure_10E.csv，而 4J 面板读的是 Figure_5J.csv —— 输出名与消费名
# 脱钩，导致这条腿永远不会被重跑刷新（改名是 09-08 手工做的，没有脚本）。
out["gene_set_version"] = "Mitoxy-80_v1.0"
#: 版本号升位：口径由"未限物种"变为"人源"，产物必须能与旧版区分。
out["score_version"] = "m4hpa_v1.1_human"
out.to_csv(BASE + r"\01_FIGURE_DATA_CSV\Main\Figure_5J.csv", index=False)
out.to_csv(BASE + r"\02_SUPPLEMENTARY_TABLES\SUPPLEMENTARY_Tables_CSV\Table_S59e_M4_Subcellular_Annotation.csv", index=False)
print("\n已写 Figure_5J.csv / Table_S59e_M4_Subcellular_Annotation.csv")
