# -*- coding: utf-8 -*-
"""
N7_renumber_supp_tables_20261002.py — 附表编号迁移（旧编号 → 最终工作簿编号）

背景（2026-10-02，承接 RENUMBER_STATUS_20260928.md 的执行清单）:
  · 编号唯一权威来源 = 最终投稿工作簿 `Supplementary_Tables_S1-S93.xlsx`
    （投稿包 05_SupplementaryTables_Notes/，2026-10-01 版）Index 页的源文件清单；
    该编号已与文稿 v5、Zenodo 大表记录（concept DOI 10.5281/zenodo.23018056）一致。
  · 本脚本把 02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV 的表文件按
    工作簿编号批量重命名；**内容零改动**（逐件 SHA256 前后断言一致）。

匹配规则:
  1) 后缀匹配：Table_S<旧组><字母>_<余部> → Table_S<新组><字母>_<余部>
     （余部含下划线后全串，如 `S15a_..._Instruments` 只唯一命中 `S10a_..._Instruments`）
  2) 工作簿 Location 列的 4 个 Zenodo 大表按 Zenodo 记录实际文件名改名
  3) 合并命名 Table_S22_S23_* → 工作簿 S16 的 Table_S16_Drug_Database /
     Table_S16_Integrated_Drug_Ranking
  4) 未出现于最终工作簿的文件不改名（保留旧编号；对照表 note 列逐条说明）

同步更新（与改名同一原子批次）:
  · 04_AUDIT_GOVERNANCE/RESULTS_MANIFEST_v2.0.csv 的 filename / s_number 两列
  · 01_FIGURE_DATA_CSV/FIGURE_DATA_MANIFEST.csv 的 `source_file(01_RESULTS_TABLES)` 列
  · 输出对照表 `02_SUPPLEMENTARY_TABLES/附表最终编号对照表_20261002.csv`

备份：_review_tmp/supp_pre_renumber_20261002/（重命名前整目录复制，仅首次）
台账：_review_tmp/supp_renumber_ledger_20261002.csv（逐件 旧名/新名/SHA256）

用法:  python N7_renumber_supp_tables_20261002.py  [--dry-run]
"""
from __future__ import annotations

import csv
import hashlib
import io
import os
import re
import shutil
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(HERE, "02_SUPPLEMENTARY_TABLES", "SUPPLEMENTARY_Tables_CSV")
SUPP_DIR = os.path.join(HERE, "02_SUPPLEMENTARY_TABLES")
MANIFEST = os.path.join(HERE, "04_AUDIT_GOVERNANCE", "RESULTS_MANIFEST_v2.0.csv")
FIGMAN = os.path.join(HERE, "01_FIGURE_DATA_CSV", "FIGURE_DATA_MANIFEST.csv")
REVIEW = os.path.join(HERE, "_review_tmp")
BACKUP = os.path.join(REVIEW, "supp_pre_renumber_20261002")
LEDGER = os.path.join(REVIEW, "supp_renumber_ledger_20261002.csv")
XLSX = (r"E:\SCI\CELL DEATH AND DIFFERENTIATION投稿上传_UPLOAD_READY"
        r"\CELL DEATH AND DIFFERENTIATION投稿上传_UPLOAD_READY"
        r"\05_SupplementaryTables_Notes\Supplementary_Tables_S1-S93.xlsx")

DRY = "--dry-run" in sys.argv

# 特例（工作簿 Location 列 / 合并命名）：旧名 → 新名
SPECIAL = {
    "Table_S15f_SMR_HEIDI_GeneDetail.csv": "Table_S10f_SMR_HEIDI_GeneDetail.csv",
    "Table_S55d_M1_Visium_SVG_Detail.csv": "Table_S28d_M1_Visium_SVG_Detail.csv",
    "Table_S7k_Velocity_Pseudotime.csv": "Table_S5k_Velocity_Pseudotime.csv",
    "Table_S26k_CellType_Scores.csv": "Table_S79k_CellType_Scores.csv",
    "Table_S22_S23_Drug_Database.csv": "Table_S16_Drug_Database.csv",
    "Table_S22_S23_Integrated_Drug_Ranking.csv": "Table_S16_Integrated_Drug_Ranking.csv",
}

# 最终工作簿列出、但本目录无独立源文件的表（说明用，不改名、不新建）
NO_LOCAL_FILE = [
    ("Table_S66_WetLab_WB_Densitometry.csv", "Table S66", "随最终工作簿交付（湿实验；源数据在 111湿实验原始数据/ 与图 8 源包）"),
    ("Table_S67_WetLab_Summary.csv", "Table S67", "同上"),
    ("Table_S68_WetLab_Stats.csv", "Table S68", "同上"),
]


def sha256_of(path: str, block: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(block), b""):
            h.update(b)
    return h.hexdigest()


def newline_of(path: str) -> str:
    with open(path, "rb") as f:
        return "\r\n" if b"\r\n" in f.read(1 << 16) else "\n"


def key_of(fn: str):
    m = re.match(r"^Table_S(\d+)([a-z]*)_(.+)$", fn)
    return (m.group(2), m.group(3)) if m else None


def num_of(fn: str):
    m = re.match(r"^Table_S(\d+)", fn)
    return "S" + m.group(1) if m else ""


# ---------- 1) 从最终工作簿 Index 推导映射 ----------
print("== 1) 读取最终工作簿 Index（唯一权威编号来源） ==")
print("   ", XLSX)
print("    workbook sha256 =", sha256_of(XLSX)[:16], "…")
import openpyxl
wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
new_by_suffix, dup = {}, []
for r in wb["Index"].iter_rows(min_row=2, values_only=True):
    if not r[0]:
        continue
    for fn in [x.strip() for x in str(r[3] or "").split(";") if x.strip() and x.strip() != "—"]:
        k = key_of(fn)
        if k is None:
            continue
        if k in new_by_suffix:
            dup.append((k, new_by_suffix[k], fn))
        new_by_suffix[k] = fn
assert not dup, "工作簿源文件后缀不唯一: %s" % dup[:3]
print("    Index 表数 93；可解析源文件键 %d（后缀唯一）" % len(new_by_suffix))

files = sorted(os.listdir(TAB))
mapping, kept = {}, []
for fn in files:
    if not fn.lower().endswith((".csv", ".gz")):
        continue
    if fn in SPECIAL:
        mapping[fn] = SPECIAL[fn]
        continue
    k = key_of(fn)
    if k and k in new_by_suffix:
        mapping[fn] = new_by_suffix[k]
    else:
        kept.append(fn)

# 断言：特例目标与后缀表无冲突
for a, b in SPECIAL.items():
    assert a not in mapping or mapping[a] == b
n_same = sum(1 for a, b in mapping.items() if a == b)
print("    本目录表文件 %d：改名 %d（其中 %d 个新旧同名）、保留原名 %d"
      % (len([f for f in files if f.lower().endswith((".csv", ".gz"))]),
         len(mapping), n_same, len(kept)))

targets = list(mapping.values())
assert len(targets) == len(set(targets)), "新名冲突"
for t in targets:
    if t not in files:
        pass  # 目标不存在 → 正常
    else:
        assert mapping.get(t) == t, "链式/占用冲突: " + t

if DRY:
    print("\n[dry-run] 全量改名清单（旧 → 新）：")
    for a in sorted(mapping):
        flag = "  " if mapping[a] != a else "= "
        print("   %s%s  ->  %s" % (flag, a, mapping[a]))
    print("\n[dry-run] 保留原名 %d 个：" % len(kept))
    for a in kept:
        print("      ", a)
    sys.exit(0)

# ---------- 2) 备份（仅首次） ----------
print()
print("== 2) 备份 ==")
os.makedirs(REVIEW, exist_ok=True)
if not os.path.isdir(BACKUP):
    shutil.copytree(TAB, BACKUP)
    print("   已备份 %d 个文件 → %s" % (len(os.listdir(BACKUP)), BACKUP))
else:
    print("   备份已存在（跳过）：%s（%d 文件）" % (BACKUP, len(os.listdir(BACKUP))))

# ---------- 3) 批量改名（逐件 SHA256 前后断言） ----------
print()
print("== 3) 批量改名（内容零改动） ==")
ledger_rows, n_ren, n_skip = [], 0, 0
for old in sorted(mapping):
    new = mapping[old]
    src, dst = os.path.join(TAB, old), os.path.join(TAB, new)
    if old == new:
        ledger_rows.append((old, new, sha256_of(src), "同名"))
        continue
    if not os.path.exists(src) and os.path.exists(dst):
        ledger_rows.append((old, new, sha256_of(dst), "已改（重跑跳过）"))
        n_skip += 1
        continue
    assert os.path.exists(src), "缺源文件: " + old
    assert not os.path.exists(dst), "目标已存在: " + new
    h1 = sha256_of(src)
    os.rename(src, dst)
    h2 = sha256_of(dst)
    assert h1 == h2, "SHA256 变化: " + old
    ledger_rows.append((old, new, h1, "改名"))
    n_ren += 1
print("   改名 %d 个（重跑跳过 %d 个）；SHA256 逐件一致 ✓" % (n_ren, n_skip))

with open(LEDGER, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f, lineterminator="\r\n")
    w.writerow(["旧文件名", "新文件名", "sha256(未变)", "处理"])
    w.writerows(ledger_rows)
print("   台账 →", LEDGER)

# ---------- 4) 输出对照表（含保留原名与无源文件项） ----------
print()
print("== 4) 对照表 ==")
# 旧组 → 新组（组级换算，供人工快速查阅）
group_map = {}
for old in sorted(mapping):
    if mapping[old] == old:
        continue
    og = re.match(r"^Table_S(\d+)", old).group(1)
    ng = re.match(r"^Table_S(\d+)", mapping[old]).group(1)
    group_map.setdefault("S%s" % og, set()).add("S%s" % ng)
out_csv = os.path.join(SUPP_DIR, "附表最终编号对照表_20261002.csv")
with open(out_csv, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f, lineterminator="\r\n")
    w.writerow(["旧文件名", "新文件名", "最终表号", "处理", "备注"])
    for old in sorted(mapping):
        new = mapping[old]
        note = ""
        if old in SPECIAL:
            note = "工作簿 Location/Zenodo 命名特例"
        if old == new:
            note = note or "编号未变"
        w.writerow([old, new, num_of(new), "改名" if old != new else "编号未变", note])
    for old in kept:
        w.writerow([old, old, num_of(old), "保留原名", "未纳入最终工作簿（2026-10-01 版）；保留旧编号供溯源"])
    for fn, tab, note in NO_LOCAL_FILE:
        w.writerow(["（无独立文件）", fn, tab, "工作簿内交付", note])
print("   →", out_csv)
print("   组级换算：", "；".join("%s→%s" % (k, "/".join(sorted(v)))
                                  for k, v in sorted(group_map.items(), key=lambda kv: int(kv[0][1:]))))

# ---------- 5) 更新 RESULTS_MANIFEST_v2.0.csv ----------
print()
print("== 5) 更新 manifest（filename / s_number） ==")
tabloc = "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV"
raw = open(MANIFEST, encoding="utf-8-sig", newline="").read()
rows = list(csv.DictReader(io.StringIO(raw)))
FIELDS = list(rows[0].keys())
n_up, other_loc = 0, []
for r in rows:
    fn = r["filename"]
    if fn in mapping and mapping[fn] != fn:
        if r["location"] != tabloc:
            other_loc.append((fn, r["location"], r["status"]))
            continue
        r["filename"] = mapping[fn]
        r["s_number"] = num_of(mapping[fn])
        n_up += 1
print("   更新 %d 行；同名单但 location 在别处（不动）：%d" % (n_up, len(other_loc)))
for fn, loc, st in other_loc:
    print("      [留] %s | location=%s | status=%s" % (fn, loc, st))
with open(MANIFEST, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator=newline_of(MANIFEST))
    w.writeheader()
    w.writerows(rows)
print("   写回", MANIFEST)

# ---------- 6) 更新图数据包溯源列 ----------
print()
print("== 6) 更新 FIGURE_DATA_MANIFEST.csv 溯源列 ==")
SRCKEY = "source_file(01_RESULTS_TABLES)"
fraw = open(FIGMAN, encoding="utf-8-sig", newline="").read()
frows = list(csv.DictReader(io.StringIO(fraw)))
FF = list(frows[0].keys())
assert SRCKEY in FF
n_fig = 0
for r in frows:
    v = (r.get(SRCKEY) or "").strip()
    if v in mapping and mapping[v] != v:
        r[SRCKEY] = mapping[v]
        n_fig += 1
print("   溯源列更新 %d 处（共 %d 行）" % (n_fig, len(frows)))
with open(FIGMAN, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=FF, lineterminator=newline_of(FIGMAN))
    w.writeheader()
    w.writerows(frows)
print("   写回", FIGMAN)

print()
print("== 完成 ==")
print("下一步：① python lint_package.py（须全绿）② python P0_build_release_dir.py（重建发布目录）")
