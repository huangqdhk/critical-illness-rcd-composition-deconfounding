# -*- coding: utf-8 -*-
"""
P0_repair_manifest.py — 修复 RESULTS_MANIFEST 漂移并将 2026-08-20/21 新文件登记入册
==================================================================================
诊断结论（2026-08-21）：
1. 版本戳漂移：M2/M3/M4 批次重建 manifest 时把 score_version 空串写成 ''，而冻结 CSV 内为 'NA'；
2. 5 个 CSV（Figure_10D、S15e Instruments/LeaveOneOut/Sensitivity、S60）缺版本列；
3. 6 个 TXT 报告缺版本页脚；
4. N22：S15e/f/g、S19i 的 s_number 记成了带字母后缀，lint 正则只产出 S15/S19；
5. archive 重组（R1_table_backups/ 等子目录）后 manifest location 未更新；README_归档说明.md 未登记；
6. P0/P1/P2 系列（2026-08-20）与本轮新文件未登记。
修复原则：
- 冻结 CSV/TXT 的内容除"补版本戳"外零改动；改动文件的 SHA256/size 同步重算入册；
- 其余 manifest 元数据修正（location/s_number/版本值/新行）不影响任何冻结数值。
输出：更新 RESULTS_MANIFEST_v2.0.csv（v2.1 变更记录写入 RESULTS_MANIFEST.md 头部）。
"""
import csv, hashlib, os, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
GOV = os.path.join(ROOT, "04_AUDIT_GOVERNANCE")
TAB = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES", "SUPPLEMENTARY_Tables_CSV")
FIG = os.path.join(ROOT, "01_FIGURE_DATA_CSV", "Main")
MAN = os.path.join(GOV, "RESULTS_MANIFEST_v2.0.csv")

def sha256_of(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

def read_csv_rows(p):
    raw = open(p, "rb").read()
    enc = "utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf-8"
    with open(p, encoding=enc, newline="") as f:
        return list(csv.reader(f)), enc

def write_csv_rows(p, rows, enc):
    with open(p, "w", encoding=enc, newline="") as f:
        csv.writer(f, lineterminator="\r\n").writerows(rows)

# ---- load manifest ----
all_rows, rows_enc = read_csv_rows(MAN)
cols = all_rows[0]
data = []
with open(MAN, encoding=rows_enc, newline="") as f:
    for r in csv.DictReader(f):
        data.append(r)
print(f"manifest rows: {len(data)}")
by_fn = {r["filename"]: r for r in data}

def setrow(r, **kw):
    for k, v in kw.items():
        r[k] = v

CHANGED_FILES = []  # 内容被改动的文件 -> 重新算 sha

# ============ 1) 补版本戳：5 个缺列 CSV ============
for fn, loc, g, s in [
    ("Figure_10D.csv", FIG, "Mitoxy-80_v1.0", "m4gem_v1.1"),
    ("Table_S15e_pQTL_MR_Instruments.csv", TAB, "Mitoxy-80_v1.0", "m3pqtl_v1.0"),
    ("Table_S15e_pQTL_MR_LeaveOneOut.csv", TAB, "Mitoxy-80_v1.0", "m3pqtl_v1.0"),
    ("Table_S15e_pQTL_MR_Sensitivity.csv", TAB, "Mitoxy-80_v1.0", "m3pqtl_v1.0"),
    ("Table_S60_M4_GEM_Flux_Comparison.csv", TAB, "Mitoxy-80_v1.0", "m4gem_v1.1"),
]:
    p = os.path.join(loc, fn)
    rows, enc = read_csv_rows(p)
    h = rows[0]
    if "gene_set_version" not in h:
        rows[0] = h + ["gene_set_version", "score_version"]
        for r in rows[1:]:
            r += [g, s]
        write_csv_rows(p, rows, enc)
        CHANGED_FILES.append((fn, loc))
        print(f"[版本列补] {fn}")
    else:
        print(f"[跳过-已有版本列] {fn}")

# ============ 2) TXT 报告版本值：manifest 空串 -> NA/NA（文件零改动，lint 按 NA 跳过） ============
TXT_NA = ["Table_S14_NicheNet_Analysis_Report.txt", "Table_S16_Immune_Infiltration_Analysis_Report.txt",
          "Table_S1_rescore_report.txt", "Table_S26_Analysis_Report.txt",
          "Table_S29_S31_scATAC_SCENIC_Analysis_Report.txt", "Table_S32_S33_Methodology_Note.txt"]
for fn in TXT_NA:
    if fn in by_fn and by_fn[fn]["gene_set_version"] in ("", None):
        setrow(by_fn[fn], gene_set_version="NA", score_version="NA")
        print(f"[TXT版本值] {fn} -> NA/NA（文件零改动；该报告无基因集版本依赖，登记 NA）")

# ============ 3) manifest 元数据修正 ============
# 3a) N22：S15e/f/g -> S15；S19i -> S19
for fn, sn in [("Table_S15e_pQTL_MR_Results.csv", "S15"), ("Table_S15e_pQTL_MR_Sensitivity.csv", "S15"),
               ("Table_S15e_pQTL_MR_Instruments.csv", "S15"), ("Table_S15e_pQTL_MR_LeaveOneOut.csv", "S15"),
               ("Table_S15e_coloc_Results.csv", "S15"), ("Table_S15e_coloc_snp_pp.csv", "S15"),
               ("Table_S15f_SMR_HEIDI_Results.csv", "S15"), ("Table_S15f_SMR_HEIDI_GeneDetail.csv", "S15"),
               ("Table_S15g_TWAS_Results.csv", "S15"),
               ("Table_S19i_SeismicGWAS_Spatial_CellType_Chain.csv", "S19")]:
    if fn in by_fn:
        setrow(by_fn[fn], s_number=sn)
        print(f"[N22] {fn} s_number -> {sn}")

# 3b) archive location 修正（文件实际位于 archive 子目录）
arch_map = {}
for dirpath, dirs, files in os.walk(os.path.join(ROOT, "archive")):
    for f in files:
        if f in by_fn:
            rel = os.path.relpath(os.path.join(dirpath, f), ROOT).replace("\\", "/").replace("/" + f, "")
            arch_map[f] = rel
for fn, rel in sorted(arch_map.items()):
    r = by_fn[fn]
    if r.get("location") != rel:
        old = r.get("location")
        setrow(r, location=rel)
        print(f"[archive] {fn}: location {old!r} -> {rel!r}")

# 3c) 版本值修正：manifest 空串 vs CSV 'NA'（以冻结 CSV 为准）
fixed_na = 0
for fn, r in by_fn.items():
    if r["status"] != "canonical" or not fn.lower().endswith(".csv"):
        continue
    if r["location"] in ("", "NA"):
        p = os.path.join(ROOT, fn)
    else:
        p = os.path.join(ROOT, r["location"], fn)
    if not os.path.exists(p):
        continue
    rows, enc = read_csv_rows(p)
    h = rows[0]
    if h and str(h[0]).startswith("#"):
        h = rows[1] if len(rows) > 1 else h
    if "gene_set_version" not in h or "score_version" not in h:
        continue
    gi, si = h.index("gene_set_version"), h.index("score_version")
    row0 = None
    for rr in rows[1:]:
        if rr and not str(rr[0]).startswith("#") and len(rr) > max(gi, si):
            row0 = rr
            break
    if row0 is None:
        continue
    if (r["gene_set_version"] in ("", None) and row0[gi]) or (r["score_version"] in ("", None) and row0[si]):
        setrow(r, gene_set_version=r["gene_set_version"] or row0[gi],
               score_version=r["score_version"] or row0[si])
        fixed_na += 1
print(f"[版本值] 修正空串->CSV实际值: {fixed_na} 行")

# ============ 4) 重算内容改动文件的 sha/size ============
for fn, loc in CHANGED_FILES:
    p = os.path.join(loc, fn)
    r = by_fn[fn]
    setrow(r, size_bytes=str(os.path.getsize(p)), sha256=sha256_of(p))
    n_rows = sum(1 for _ in open(p, "rb")) - 1 if fn.endswith(".csv") else r.get("n_data_rows", "")
    if fn.endswith(".csv"):
        setrow(r, n_data_rows=str(n_rows))
    print(f"[sha] {fn} 重冻结")

# ============ 5) 登记新文件 ============
def addrow(sn, fn, loc, cls, status, universe, g, s, analysis, note=""):
    p = os.path.join(loc, fn)
    n_rows = ""
    if fn.endswith(".csv"):
        n_rows = str(sum(1 for _ in open(p, "rb")) - 1)
    data.append({
        "s_number": sn, "filename": fn, "location": loc.replace(ROOT + "\\", "").replace("\\", "/"),
        "class": cls, "status": status, "gene_universe": universe,
        "gene_set_version": g, "score_version": s,
        "n_data_rows": n_rows, "size_bytes": str(os.path.getsize(p)),
        "sha256": sha256_of(p), "analysis": analysis, "notes": note,
    })

LOC_GOV = os.path.join(ROOT, "04_AUDIT_GOVERNANCE")
LOC_ARCH = os.path.join(ROOT, "archive")

new_files = [
    # P0-P2 系列（2026-08-20）+ 本轮（2026-08-21）治理文件
    ("", "P0_FROZEN_ANALYSIS_PLAN_v1.0.md", LOC_GOV, "governance", "canonical", "none", "NA", "NA",
     "P0 冻结分析方案 v1.0（OSF osf.io/C7RYD，结果查看前冻结）", ""),
    ("", "P0_Novelty_Search_Table_v1.0.md", LOC_GOV, "governance", "canonical", "none", "NA", "NA",
     "P0-2 查新表 v1.0（2026-08-20 基线）", ""),
    ("", "P0_Novelty_Search_Table_v1.1.md", LOC_GOV, "governance", "canonical", "none", "NA", "NA",
     "P0-2 查新表 v1.1（2026-08-21 PubMed 重跑更新；更正 v1.0 期刊名登记错误）", ""),
    ("", "P0_Novelty_Search_Rerun_Log_20260821.md", LOC_GOV, "governance", "canonical", "none", "NA", "NA",
     "查新检索重跑日志（2026-08-21；Q1=11/Q2=0/Q3=1）", ""),
    ("", "P0_GSE310929_Overlap_Audit_Report.md", LOC_GOV, "audit_report", "canonical", "none", "NA", "NA",
     "P0-3 样本重叠审计（GSE310929 全量包含 GSE185263/GSE32707）", ""),
    ("", "P0_Source_Exclusive_Meta_Report.md", LOC_GOV, "audit_report", "canonical", "none", "NA", "NA",
     "P0-4 来源互斥 meta 重算（k=16，REML+Hartung-Knapp）", ""),
    ("", "P0_Donor_Pseudobulk_Report.md", LOC_GOV, "audit_report", "canonical", "none", "NA", "NA",
     "P0-5 供者级 pseudobulk（28 项检验，2 项存活）", ""),
    ("", "P0_Old_vs_New_Comparison.md", LOC_GOV, "audit_report", "canonical", "none", "NA", "NA",
     "P0-6 原结果 vs 正确统计单位结果对照表", ""),
    ("", "P0_Spatial_Mapping_Report.md", LOC_GOV, "audit_report", "canonical", "none", "NA", "NA",
     "P0-7 空间 patient-slide-FOV 映射 + P0-7b 患者级置换重算（2026-08-21）", ""),
    ("", "P0-8_Framework_Rename_and_Usage_Statement.md", LOC_GOV, "governance", "canonical", "none", "NA", "NA",
     "P0-8 80 基因框架用途更名交付文档（不动成员与数值）", ""),
    ("", "P0_Manuscript_Source_Check_Report.md", LOC_GOV, "audit_report", "canonical", "none", "NA", "NA",
     "Manuscript-to-source 自动核对（2026-08-21；17/17 PASS）", ""),
    ("", "P0_References_1to64_Verification_Report.md", LOC_GOV, "audit_report", "canonical", "none", "NA", "NA",
     "参考文献 1-64 重建与 PubMed 核实报告（2026-08-21）", ""),
    ("", "P0_Software_Snapshot_20260821.md", LOC_GOV, "governance", "canonical", "none", "NA", "NA",
     "软件环境快照（pip freeze×2；R sessionInfo 缺失如实登记）", ""),
    ("", "pip_freeze_spatial_venv_20260821.txt", LOC_GOV, "governance", "canonical", "none", "NA", "NA",
     ".venv_spatial pip freeze 快照（2026-08-21）", ""),
    ("", "pip_freeze_pyaging_venv_20260821.txt", LOC_GOV, "governance", "canonical", "none", "NA", "NA",
     ".venv_pyaging pip freeze 快照（2026-08-21）", ""),
    ("", "P1_Interaction_Signature_Report.md", LOC_GOV, "audit_report", "canonical", "none", "NA", "NA",
     "P1-2/P1-2b/P1-3/P1-4 GSE235046 交互签名与匹配基因集检验报告（Gate 1 通过）", ""),
    ("", "P1_IIAMD_signature_v1.0.csv", LOC_GOV, "analysis_output", "canonical", "genome_wide",
     "NA", "p1_iiamd_v1.0", "P1 IIAMD v1.0 实验锚定签名（up 4,367/down 4,720）", ""),
    ("", "P2_Gate2_HeadToHead_Report.md", LOC_GOV, "audit_report", "canonical", "none", "NA", "NA",
     "P2-2 Gate 2 头对头报告（Gate 2 FAIL，疾病映射撤回）", ""),
    ("", "P2_competing_panels_v1.0.csv", LOC_GOV, "analysis_output", "canonical", "multi_panel",
     "NA", "p2panels_v1.0", "P2 七个预注册竞争面板基因清单", ""),
    ("", "README_归档说明.md", LOC_ARCH, "governance", "archived", "none", "NA", "NA",
     "archive 目录结构说明（子目录重组）", ""),
]
existing = {r["filename"] for r in data}
added = 0
for item in new_files:
    fn = item[1]
    if fn in existing:
        continue
    # 两个新 CSV 需要先补版本列
    if fn.endswith(".csv"):
        p = os.path.join(item[2], fn)
        rows, enc = read_csv_rows(p)
        if "gene_set_version" not in rows[0]:
            rows[0] = rows[0] + ["gene_set_version", "score_version"]
            for r in rows[1:]:
                r += [item[7], item[8]]
            write_csv_rows(p, rows, enc)
    addrow(*item)
    added += 1
    print(f"[登记] {fn}")

# ============ 6) 写回 manifest ============
with open(MAN, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, lineterminator="\r\n")
    w.writeheader()
    w.writerows(data)
print(f"manifest 写回：{len(data)} 行（新增 {added}）")
print("DONE repair")
