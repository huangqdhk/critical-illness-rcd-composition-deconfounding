# -*- coding: utf-8 -*-
"""
M15_step4b_manifest_repair.py — 修复 M15_step4 重写造成的 manifest 损伤
================================================================================
事故：M15_step4 以 pd.read_csv 默认 na_values 读取 manifest，把字面量 "NA" 的
版本戳单元格（gene_set_version/score_version）解析为 NaN，回写后丢失。
修复口径（基于"今晨 lint 74/0 全绿"这一事实）：
  1. 以 keep_default_na=False, na_values=[] 读取（不再二次损伤）。
  2. 对 canonical CSV 行（跳过 MISSING_OK/硬编码跳过/豁免名单）按 lint 同口径从
     文件本身读取第 0 行版本戳，写回 manifest——lint 全绿时 manifest 值与文件戳必然一致，
     故该重建忠实于原值（含字面 "NA"）。
  3. 对 canonical .txt 行：空版本戳单元格恢复为字面 "NA"。
  4. 与 RESULTS_MANIFEST_v1.0.csv 交叉核对：凡 v1.0 为字面 "NA" 而 v2.0 现为空的
     单元格，恢复 "NA"（仅 gene_set_version/score_version 两列有此类，已含于 2/3）。
  5. 新登记行 s_number 由 "S84b" 类修正为基础编号 "S84"（N22 校验要求）。
  6. 登记 04_AUDIT_GOVERNANCE/M15_Tool_Adoption_Report.md（治理目录全部文件须登记）。
输出：修复后 manifest；03_LOGS/M15_manifest_repair_log.txt。
"""
import csv
import hashlib
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import pandas as pd

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
MANIFEST = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "RESULTS_MANIFEST_v2.0.csv")
V1 = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "RESULTS_MANIFEST_v1.0.csv")

MISSING_OK = {
    "Table_S11_Mitoxyperilysis_TF_PerCell.csv",
    "Table_S11_TF_Activity_Matrix_AllCells.csv",
    "Table_S1_Mitoxyperilysis_Score_by_CellType_prev38.csv",
    "Table_S1_Mitoxyperilysis_Score_by_CellType_Summary_prev38.csv",
    "GSE67530_beta_for_GrimAge.csv.gz",
}
HARD_SKIP = {"GSE185263_groups.csv", "MR_bio_CRP_Sepsis.csv", "MR_bio_Ferritin_Sepsis.csv",
             "Table_S50_D1b_GSE32707_Groups.csv"}
# 豁免名单从 lint_package.py 同步读取
lint_txt = open(os.path.join(ROOT, "lint_package.py"), encoding="utf-8").read()
m = re.search(r"M1X_NO_VERSION_STAMP = \{(.*?)\n\}", lint_txt, re.S)
EXEMPT = set(re.findall(r'"([^"]+\.csv)"', m.group(1)))

log = []
def note(s=""):
    log.append(s)
    print(s, flush=True)

def sha256_of(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()

# ======================================================================
note("== 1. 读取（keep_default_na=False）")
v2 = pd.read_csv(MANIFEST, keep_default_na=False, na_values=[], low_memory=False)
v2["filename"] = v2["filename"].astype(str)
note(f"   v2 行数 {len(v2)}")

# ======================================================================
note("\n== 2. canonical CSV 行版本戳重建（lint 同口径）")
n_rebuilt = 0
for idx, r in v2.iterrows():
    fn = r["filename"]
    if r["status"] != "canonical":
        continue
    if not fn.lower().endswith(".csv"):
        continue
    if fn in MISSING_OK or fn in HARD_SKIP or fn in EXEMPT:
        continue
    p = os.path.join(ROOT, r["location"].replace("/", "\\"), fn)
    if not os.path.exists(p):
        note(f"   [跳过] 文件不存在: {fn}")
        continue
    raw = open(p, "rb").read()
    enc = "utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf-8"
    gi = si = -1
    row0 = None
    try:
        with open(p, "r", encoding=enc, newline="") as f:
            rd = csv.reader(f)
            hdr = next(rd)
            if hdr and str(hdr[0]).startswith("#"):
                hdr = next(rd)
            if "gene_set_version" not in hdr or "score_version" not in hdr:
                note(f"   [缺列] {fn}（不做重建）")
                continue
            gi, si = hdr.index("gene_set_version"), hdr.index("score_version")
            for row in rd:
                if row and str(row[0]).startswith("#"):
                    continue
                row0 = row
                break
    except StopIteration:
        continue
    if row0 is None:
        continue
    gsv = row0[gi] if len(row0) > gi else ""
    sv = row0[si] if len(row0) > si else ""
    if str(r["gene_set_version"]) != gsv or str(r["score_version"]) != sv:
        v2.at[idx, "gene_set_version"] = gsv
        v2.at[idx, "score_version"] = sv
        n_rebuilt += 1
note(f"   重建 {n_rebuilt} 行")

# ======================================================================
note("\n== 3. canonical .txt 行：空版本戳 -> 字面 'NA'")
n_txt = 0
for idx, r in v2.iterrows():
    if r["status"] != "canonical" or not str(r["filename"]).lower().endswith(".txt"):
        continue
    if str(r["gene_set_version"]) == "":
        v2.at[idx, "gene_set_version"] = "NA"
        n_txt += 1
    if str(r["score_version"]) == "":
        v2.at[idx, "score_version"] = "NA"
        n_txt += 1
note(f"   .txt 行恢复 'NA' {n_txt} 处")

# ======================================================================
note("\n== 4. v1.0 交叉核对（字面 NA 恢复兜底）")
v1 = pd.read_csv(V1, keep_default_na=False, na_values=[], low_memory=False)
v1["filename"] = v1["filename"].astype(str)
v1_map = {r["filename"]: r for _, r in v1.iterrows()}
n_v1 = 0
for idx, r in v2.iterrows():
    old = v1_map.get(r["filename"])
    if old is None:
        continue
    for c in ("gene_set_version", "score_version"):
        if str(old[c]) == "NA" and str(r[c]) == "":
            v2.at[idx, c] = "NA"
            n_v1 += 1
note(f"   v1.0 兜底恢复 {n_v1} 处")

# ======================================================================
note("\n== 5. 新行 s_number 基础编号修正（N22）")
n_snum = 0
for idx, r in v2.iterrows():
    fn = r["filename"]
    mm = re.match(r"[Tt]able_([Ss]\d+)", fn)
    if mm and r["s_number"] not in ("", mm.group(1).upper()):
        v2.at[idx, "s_number"] = mm.group(1).upper()
        n_snum += 1
note(f"   修正 {n_snum} 行 s_number")

# ======================================================================
note("\n== 6. 登记 M15_Tool_Adoption_Report.md")
if "M15_Tool_Adoption_Report.md" not in set(v2["filename"]):
    p = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "M15_Tool_Adoption_Report.md")
    new_row = pd.DataFrame([dict(
        s_number="", filename="M15_Tool_Adoption_Report.md",
        location="04_AUDIT_GOVERNANCE", class_="governance", status="canonical",
        gene_universe="none", gene_set_version="NA", score_version="NA",
        n_data_rows="", size_bytes=os.path.getsize(p),
        sha256=sha256_of(p),
        analysis="M15 工具封装与采纳演示治理报告（2026-08-27）",
        notes="M15 批次补登")])
    v2 = pd.concat([v2, new_row], ignore_index=True, sort=False)
    note("   已登记 M15_Tool_Adoption_Report.md")
else:
    note("   已在 manifest")

# ======================================================================
note("\n== 7. 回写")
v2.to_csv(MANIFEST, index=False, na_rep="")
note(f"   已写 {MANIFEST}（{len(v2)} 行）")
open(os.path.join(ROOT, "03_LOGS", "M15_manifest_repair_log.txt"), "w", encoding="utf-8").write("\n".join(log))
print("\nDONE repair")
