# -*- coding: utf-8 -*-
"""
N6_integrate_wetlab_figdata_20261002.py — 把正文湿实验三图的面板源数据并入图数据包

背景（2026-10-02，经作者确认执行）:
  · 正文 Figure 8（湿实验，8 面板）、Figure S15（5 面板）、Figure S16（6 面板）的源数据
    此前按 09-28 台账以"外部登记（ext:）"方式留在包外；本次按该台账**命名口径①**入包：
      图 8 面板 -> Main/Figure_8w{A..H}.csv   （w = wet-lab；与包内既有同名不同内容的
                                               Figure_8A–8G.csv（正文 Figure 6 世代）区隔）
      S15/S16  -> Supplementary/Figure_S15{A..G}.csv、Figure_S16{A..F}.csv（沿用原文件名）
  · 本包 FIGURE_CROSSWALK.csv 以 `可视化/00_对照表/FIGURE_CROSSWALK.csv`（09-30 版，
    出图管线唯一真源）为基准重建：19 行 ext: 路径改写为包内路径；同时吸收该版本对
    S10/S11 的既有修正。FIGURE_DATA_MANIFEST.csv 由 crosswalk 重建，描述列按
    package_file 从旧 manifest 继承（文件级匹配，避免世代标签错位）。
  · 数值零改动：拷入的 21 个 CSV 与源文件逐字节一致（脚本内 MD5 断言）。

用法：python N6_integrate_wetlab_figdata_20261002.py  [--dry-run]
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

HERE = os.path.dirname(os.path.abspath(__file__))          # 项目根
PKG = os.path.join(HERE, "01_FIGURE_DATA_CSV")
VISX = os.path.join(HERE, "可视化", "00_对照表", "FIGURE_CROSSWALK.csv")
BACKUP = os.path.join(HERE, "_review_tmp", "figdata_pre_integration_20261002")

WET8 = os.path.join(HERE, "Figure_8", "Figure_8")
WETS = os.path.join(HERE, "Figure_S15-16待可视化数据")

DRY = "--dry-run" in sys.argv


def md5(p: str) -> str:
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def newline_of(p: str) -> str:
    with open(p, "rb") as f:
        return "\r\n" if b"\r\n" in f.read(1 << 16) else "\n"


# ---------- 1) 拷贝清单 ----------
COPIES = []
for L in "ABCDEFGH":
    COPIES.append((os.path.join(WET8, f"Figure_8{L}.csv"),
                   os.path.join(PKG, "Main", f"Figure_8w{L}.csv")))
for L in "ABCDEFG":
    COPIES.append((os.path.join(WETS, f"Figure_S15{L}.csv"),
                   os.path.join(PKG, "Supplementary", f"Figure_S15{L}.csv")))
for L in "ABCDEF":
    COPIES.append((os.path.join(WETS, f"Figure_S16{L}.csv"),
                   os.path.join(PKG, "Supplementary", f"Figure_S16{L}.csv")))

print("== 1) 拷贝 %d 个面板 CSV（内容零改动） ==" % len(COPIES))
os.makedirs(BACKUP, exist_ok=True)
for src, dst in COPIES:
    assert os.path.isfile(src), "源缺失: " + src
    if os.path.exists(dst):
        same = md5(src) == md5(dst)
        print("   [skip] 已存在 %s （内容一致=%s）" % (os.path.relpath(dst, HERE), same))
        continue
    if DRY:
        print("   [dry] %s -> %s" % (os.path.relpath(src, HERE), os.path.relpath(dst, HERE)))
        continue
    shutil.copy2(src, dst)
    assert md5(src) == md5(dst)
    print("   [ok] %s -> %s  md5=%s" % (os.path.relpath(src, HERE), os.path.relpath(dst, HERE), md5(dst)))

# ---------- 2) crosswalk 重建（以 09-30 可视化版为基准） ----------
print()
print("== 2) 重建 FIGURE_CROSSWALK.csv ==")
for f in ("FIGURE_CROSSWALK.csv", "FIGURE_DATA_MANIFEST.csv", "README_图数据包说明.md"):
    p = os.path.join(PKG, f)
    b = os.path.join(BACKUP, f)
    if os.path.isfile(p) and not os.path.isfile(b) and not DRY:
        shutil.copy2(p, b)
print("   旧件备份目录 %s（已存在则不覆盖）" % os.path.relpath(BACKUP, HERE))

rows = list(csv.DictReader(open(VISX, encoding="utf-8-sig")))
fields = list(rows[0].keys())

EXT_MAP_FIG8 = re.compile(r"^ext:Figure_8/Figure_8([A-H])\.csv$")
EXT_MAP_SUP = re.compile(r"^ext:Figure_S15-16/(Figure_S1[56][A-G])\.csv$")

changed = 0
for r in rows:
    m = EXT_MAP_FIG8.match(r["csv_relpath"])
    if m:
        L = m.group(1)
        r["csv_relpath"] = "Main/Figure_8w%s.csv" % L
        r["pkg_panel"] = "Figure_8w%s" % L
        changed += 1
        continue
    m = EXT_MAP_SUP.match(r["csv_relpath"])
    if m:
        r["csv_relpath"] = "Supplementary/%s.csv" % m.group(1)
        changed += 1
print("   ext: 行改写 %d 行" % changed)

planned = {os.path.relpath(dst, PKG).replace(os.sep, "/") for _, dst in COPIES}
missing = [r["csv_relpath"] for r in rows
           if r["csv_relpath"] not in planned
           and not os.path.isfile(os.path.join(PKG, r["csv_relpath"].replace("/", os.sep)))]
assert not missing, "crosswalk 指向不存在的文件: %s" % missing[:5]
assert not any(r["csv_relpath"].startswith("ext:") for r in rows), "仍有 ext: 行"
print("   全部 %d 行路径校验通过；无 ext: 残留" % len(rows))

if not DRY:
    with open(os.path.join(PKG, "FIGURE_CROSSWALK.csv"), "w", encoding="utf-8-sig",
              newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator=newline_of(VISX))
        w.writeheader()
        w.writerows(rows)

# ---------- 3) manifest 重建 ----------
print()
print("== 3) 重建 FIGURE_DATA_MANIFEST.csv ==")
mp = os.path.join(PKG, "FIGURE_DATA_MANIFEST.csv")
src_manifest = os.path.join(BACKUP, "FIGURE_DATA_MANIFEST.csv")
if not os.path.isfile(src_manifest):
    src_manifest = mp
old = list(csv.DictReader(open(src_manifest, encoding="utf-8-sig")))
old_hdr = list(old[0].keys())
SRCKEY = "source_file(01_RESULTS_TABLES)"
by_file = {}
for r in old:
    if r["panel"] == "Table_4":
        continue
    prev = by_file.get(r["package_file"])
    # 同一文件若有多行，优先取源表列非空者
    if prev is None or (not (prev.get(SRCKEY) or "").strip() and (r.get(SRCKEY) or "").strip()):
        by_file[r["package_file"]] = r
print("   继承源：%s（行数 %d，表头 %s）" % (os.path.relpath(src_manifest, HERE), len(old), ",".join(old_hdr)))

def ms_label(ms_figure: str, ms_panel: str) -> str:
    m = re.match(r"^Figure_(S?\d+)$", ms_figure)
    tok = m.group(1)
    letters = ms_panel[len(tok):] if ms_panel.startswith(tok) else ms_panel
    return "Figure_%s%s" % (tok, letters.upper())

WET_SRC = re.compile(r"^(Figure_(?:8w?[A-H]|S15[A-G]|S16[A-F]))\.csv$")
new_rows = []
filled_by_file = 0
for r in rows:
    panel = ms_label(r["ms_figure"], r["ms_panel"])
    pf = r["csv_relpath"]
    o = by_file.get(pf)
    if o:
        filled_by_file += 1
    if o and (o.get("panel_content") or o.get(SRCKEY)):
        src = o.get(SRCKEY, "")
        title = o.get("figure_title", "") or r.get("ms_title", "")
        content = o.get("panel_content", "") or r.get("panel_content", "")
        plot = o.get("suggested_plot", "") or r.get("suggested_plot", "")
        notes = o.get("notes", "")
    else:
        src = ""
        if r["csv_relpath"].startswith("Main/Figure_8w") or "S15" in r["csv_relpath"] or "S16" in r["csv_relpath"]:
            src = (r.get("pkg_source_file", "") or os.path.basename(pf)) + "（湿实验源包）"
        title = r.get("ms_title", "")
        content = r.get("panel_content", "")
        plot = r.get("suggested_plot", "")
        notes = ""
    if r["ms_panel"] == "S15E":
        notes = (notes + "；" if notes else "") + "2026-09-30 二次改版后 S15 为 5 面板，本行为其第 5 面板（源文件为 Figure_S15F.csv）"
    new_rows.append({
        "panel": panel,
        "package_file": pf,
        "source_file(01_RESULTS_TABLES)": src,
        "figure_title": title,
        "panel_content": content,
        "suggested_plot": plot,
        "notes": notes,
    })
# Table 4 行原样保留在末尾
t4 = [r for r in old if r["panel"] == "Table_4"]
if t4:
    t4 = t4[0]
    new_rows.append({
        "panel": "Table_4",
        "package_file": t4["package_file"],
        "source_file(01_RESULTS_TABLES)": t4["source_file(01_RESULTS_TABLES)"],
        "figure_title": t4["figure_title"],
        "panel_content": t4["panel_content"],
        "suggested_plot": t4["suggested_plot"],
        "notes": t4["notes"],
    })
print("   → 新 manifest 行数 %d（其中按文件继承描述 %d 行；Table 4 行保留）"
      % (len(new_rows), filled_by_file))
wet = [x for x in new_rows if x["panel"].startswith(("Figure_8", "Figure_S15", "Figure_S16"))]
print("   新增湿实验行 %d：%s" % (len(wet), ", ".join(x["panel"] for x in wet)))

if not DRY:
    with open(mp, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=old_hdr, lineterminator=newline_of(mp))
        w.writeheader()
        w.writerows(new_rows)

print()
print("== 完成 ==")
print("（dry-run，未写盘）" if DRY else "文件已写入；请随后运行 P0_build_release_dir.py 重建发布目录。")
