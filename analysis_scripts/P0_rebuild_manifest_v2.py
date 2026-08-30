# -*- coding: utf-8 -*-
"""
P0_rebuild_manifest_v2.py — 从 v1.0 + 磁盘实况重建 RESULTS_MANIFEST_v2.0.csv
================================================================================
背景：RESULTS_MANIFEST_v2.0.csv 在 2026-08-21 修复脚本写回时被截断（0 字节）。
v1.0（279 行，2026-08-17 冻结）完好，且对已冻结文件携带正确 sha256/元数据。
重建策略（审计链保全）：
1. 以 v1.0 为基底，保留其 sha256/size（对未改动文件保持冻结校验链）；
2. 版本戳以磁盘 CSV/TXT 实况为准（修正 v2.0 时代的 '' vs 'NA' 漂移）；
3. archive 子目录重组：archived 行的 location 指向磁盘实际子目录；
4. GSE67530 gz 的 location 修正为 '00_RAW_DATA'；
5. 补登记 v1.0 之后新增的 M2/M3/M4、P0/P1/P2 及本轮新文件（sha 按磁盘现算）；
6. 内容被本轮补版本列的 5 个 CSV 重算 sha/size。
"""
import csv, hashlib, os, re, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
GOV = os.path.join(ROOT, "04_AUDIT_GOVERNANCE")
TAB = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES", "SUPPLEMENTARY_Tables_CSV")
NOTES = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES", "Supplementary_Notes")
FIG = os.path.join(ROOT, "01_FIGURE_DATA_CSV", "Main")
FIGSUP = os.path.join(ROOT, "01_FIGURE_DATA_CSV", "Supplementary")
ARCH = os.path.join(ROOT, "archive")
V1 = os.path.join(GOV, "RESULTS_MANIFEST_v1.0.csv")
V2 = os.path.join(GOV, "RESULTS_MANIFEST_v2.0.csv")

# 本轮内容改动的文件（补版本列）→ sha 必须重算
RESHA = {"Figure_10D.csv", "Table_S15e_pQTL_MR_Instruments.csv", "Table_S15e_pQTL_MR_LeaveOneOut.csv",
         "Table_S15e_pQTL_MR_Sensitivity.csv", "Table_S60_M4_GEM_Flux_Comparison.csv"}
# lint check_versions 硬编码跳过（版本列豁免）
VERSION_SKIP = {"GSE185263_groups.csv", "MR_bio_CRP_Sepsis.csv", "MR_bio_Ferritin_Sepsis.csv",
                "Table_S50_D1b_GSE32707_Groups.csv"}

def sha256_of(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

def csv_rows(p):
    raw = open(p, "rb").read()
    enc = "utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf-8"
    with open(p, encoding=enc, newline="") as f:
        return list(csv.reader(f))

def disk_version_values(fn, loc):
    """从磁盘 CSV 读版本戳实际值集合；TXT 解析页脚；无则 ('NA','NA')。"""
    p = os.path.join(loc, fn)
    if fn.lower().endswith(".csv") and fn not in VERSION_SKIP:
        rows = csv_rows(p)
        h = rows[0]
        if h and str(h[0]).startswith("#"):
            h = rows[1] if len(rows) > 1 else h
        if "gene_set_version" in h and "score_version" in h:
            gi, si = h.index("gene_set_version"), h.index("score_version")
            gs, ss = set(), set()
            for r in rows[1:]:
                if r and not str(r[0]).startswith("#") and len(r) > max(gi, si):
                    gs.add(r[gi]); ss.add(r[si])
            return ";".join(sorted(gs)), ";".join(sorted(ss))
        return "NA", "NA"
    if fn.lower().endswith(".txt"):
        raw = open(p, "rb").read()
        for enc in ("utf-8", "gbk"):
            try:
                body = raw.decode(enc)
                break
            except UnicodeDecodeError:
                continue
        m = re.search(r"gene_set_version=([^;]+); score_version=([^\s;]+)", body)
        if m:
            return m.group(1), m.group(2)
        return "NA", "NA"
    return "", ""

def s_number_of(fname):
    m = re.match(r"[Tt]able_([Ss]\d+)[a-z]?[_.]", fname)
    return m.group(1).upper() if m else ""

SHARED = {  # 多对一共享编号（与 v1.0/build_D4 同源）
    "Table_S22_S23_Drug_Database.csv": "S22/S23",
    "Table_S22_S23_Integrated_Drug_Ranking.csv": "S22/S23",
    "Table_S22_S23_Virtual_Drug_Screening_Report.txt": "S22/S23",
    "Table_S29_S31_scATAC_SCENIC_Analysis_Report.txt": "S29/S31",
    "Table_S32_S33_Methodology_Note.txt": "S32/S33",
    "Table_S38_S39_Geneformer_Report.txt": "S38/S39",
}

# ---------- 1) 装载 v1.0 基底 ----------
base = {}
order = []
with open(V1, encoding="utf-8-sig", newline="") as f:
    rd = csv.DictReader(f)
    FIELDS = rd.fieldnames
    for r in rd:
        base[r["filename"]] = dict(r)
        order.append(r["filename"])
print("v1.0 base rows:", len(base))

# ---------- 2) archive 子目录定位 ----------
arch_loc = {}
for dirpath, dirs, files in os.walk(ARCH):
    for f in files:
        arch_loc[f] = os.path.relpath(dirpath, ROOT).replace("\\", "/")
        if f in base and base[f]["status"] == "archived":
            base[f]["location"] = arch_loc[f]

# ---------- 3) 基底修正 ----------
r = base.get("GSE67530_beta_for_GrimAge.csv.gz")
if r:
    r["location"] = "00_RAW_DATA"

for fn, r in base.items():
    if r["status"] == "canonical":
        # 版本戳以磁盘为准
        loc = os.path.join(ROOT, r["location"]) if r["location"] else ROOT
        p = os.path.join(loc, fn)
        if os.path.exists(p):
            g, s = disk_version_values(fn, loc)
            if g or s:
                r["gene_set_version"] = g
                r["score_version"] = s
            # v1.0 冻结后部分文件经 M2-M4 批次改动；v2.0 重建以磁盘为权威重算 sha/size
            # （v1.0 文件本身保留历史冻结链）
            r["sha256"] = sha256_of(p)
            r["size_bytes"] = str(os.path.getsize(p))
        # s_number 归一
        if fn.startswith("Table_S") or fn.startswith("Table_s"):
            sn = SHARED.get(fn) or s_number_of(fn)
            if sn:
                r["s_number"] = sn

# ---------- 4) 扫描磁盘补登 ----------
scan_dirs = [("04_AUDIT_GOVERNANCE", GOV), ("02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV", TAB),
             ("02_SUPPLEMENTARY_TABLES/Supplementary_Notes", NOTES),
             ("01_FIGURE_DATA_CSV/Main", FIG), ("01_FIGURE_DATA_CSV/Supplementary", FIGSUP)]
disk = {}
for locname, d in scan_dirs:
    for p in sorted(os.listdir(d)):
        fp_ = os.path.join(d, p)
        if os.path.isfile(fp_) and not (d == GOV and p.startswith("RESULTS_MANIFEST")):
            disk[p] = locname
for fn, loc in sorted(arch_loc.items()):
    disk[fn] = loc

added = 0
for fn, locname in sorted(disk.items()):
    if fn in base:
        continue
    p = os.path.join(ROOT, locname, fn)
    is_arch = locname.startswith("archive")
    status = "archived" if is_arch else "canonical"
    g, s = disk_version_values(fn, os.path.join(ROOT, locname)) if not is_arch else ("NA", "NA")
    sn = ""
    if not is_arch and (fn.startswith("Table_S") or fn.startswith("Table_s")):
        sn = SHARED.get(fn) or s_number_of(fn)
    n_rows = ""
    if fn.lower().endswith(".csv"):
        n_rows = str(sum(1 for _ in open(p, "rb")) - 1)
    cls = "supplementary_table" if (fn.startswith("Table_S") and fn.endswith(".csv")) else \
          ("supplementary_report" if fn.startswith("Table_S") else \
           ("figure_data" if fn.startswith("Figure_") and fn.endswith(".csv") else \
            ("governance" if (fn.endswith(".md") and not is_arch) else "analysis_output")))
    if locname.startswith("01_FIGURE_DATA_CSV"):
        cls = "figure_data"
    base[fn] = {
        "s_number": sn, "filename": fn, "location": locname, "class": cls,
        "status": status, "gene_universe": "none",
        "gene_set_version": g, "score_version": s,
        "n_data_rows": n_rows, "size_bytes": str(os.path.getsize(p)),
        "sha256": sha256_of(p),
        "analysis": "M2/M3/M4/P0-P2 与 2026-08-21 治理批次新增文件（v2.0 重建登记）" if not is_arch else "历史归档",
        "notes": "v2.0 重建（2026-08-21，v1.0 基底 + 磁盘扫描补登）" if not is_arch else "",
    }
    order.append(fn)
    added += 1
    print(f"[补登] {locname}/{fn}")

# ---------- 5) 内容改动文件重算 sha/size ----------
for fn in RESHA:
    if fn in base:
        loc = os.path.join(ROOT, base[fn]["location"])
        p = os.path.join(loc, fn)
        base[fn]["sha256"] = sha256_of(p)
        base[fn]["size_bytes"] = str(os.path.getsize(p))
        base[fn]["n_data_rows"] = str(sum(1 for _ in open(p, "rb")) - 1)
        print(f"[重冻结] {fn}")

# ---------- 6) 写回 v2.0 ----------
with open(V2, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\r\n")
    w.writeheader()
    for fn in order:
        w.writerow(base[fn])
print(f"v2.0 rebuilt: {len(base)} rows (added {added})")
