# -*- coding: utf-8 -*-
"""
N5b_build_source_data_xlsx_crosswalk.py — 按**正文图号**生成 Source Data xlsx

为什么要重做（2026-09-20）
--------------------------
CDD 的 Guide to Authors 明写 "Include spread-sheet data in supplementary
materials"（GTA L1227），而 `01_FIGURE_DATA_CSV/README_图数据包说明.md:61`
承诺的 `投稿…/Source_Data/Main|Supplementary/` 此前**在任何层级都不存在**
（审稿复核标记为"投稿前最大的真实缺口"）。

旧脚本 `N5_build_source_data_xlsx_20260905.py` 按**文件名里的数字**分组，
而本库数据包有**三套图号世代**（数据包 9 张主图 / 正文 7 张主图 + 14 张附图 /
21 张交付图），包内 Figure_5A–5J 其实对应正文 Figure 4A–4J。直接沿用会产出
"Figure_5_SourceData.xlsx 里装的是正文 Figure 4 的数据"，把三套世代的混乱
原样带进投稿包。

本脚本改以 `可视化/00_对照表/FIGURE_CROSSWALK.csv` 为**唯一真源**（该表把
ms_figure × ms_panel 映射到 csv_relpath），因此：
  · 正文图号 = 输出文件名里的图号（Main/Figure_4_SourceData.xlsx 装 4A–4J）；
  · sheet 名 = 正文面板字母（Panel_J），与图面字母逐一对应；
  · 包内编号只作为溯源列写进 sheet，不参与命名。

数值零改动：逐字节读 CSV 后原样写出（保守布尔列文本，避免 Excel 往返
True→1.0 漂移；剥离前导 '#' 注释行并单列 notes sheet）。

用法::

    python N5b_build_source_data_xlsx_crosswalk.py
"""

from __future__ import annotations

import csv
import glob
import io
import os
import re
import sys

import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")

ROOT = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(ROOT, "01_FIGURE_DATA_CSV")
OUT = os.path.join(ROOT, "投稿CELL DEATH AND DIFFERENTIATION", "Source_Data")
XWALK = os.path.join(ROOT, "可视化", "00_对照表", "FIGURE_CROSSWALK.csv")

#: crosswalk 里 `csv_relpath` 的外部源前缀（正文 Figure 8 的湿实验数据包不在
#: 01_FIGURE_DATA_CSV 下）。遇到即跳过并登记，见 main() 的说明。
EXT_PREFIX = "ext:"

#: Index 页的口径说明（审稿复核 m5 之②）。这些说明**不能**写进 CSV 本体：面板
#: 脚本用纯 pandas 读这些 CSV，行首 '#' 注释会顶掉表头解析；因此登记在工作簿的
#: Index 页（与 CSV 同处一份交付物，编辑/审稿人一眼可见）。
#: 键 = (folder, ms_figure, ms_panel)。
INDEX_NOTES = {
    ("Main", "Figure_2", "2D"):
        "列名 Cuprotosis 为上游源表既有拼写（= Cuproptosis），已登记于 "
        "01_FIGURE_DATA_CSV/README_图数据包说明.md；图面与图注一律用 Cuproptosis。",
    ("Supplementary", "Figure_S4", "S4H"):
        "列名 Cuprotosis 为上游源表既有拼写（= Cuproptosis），已登记于 "
        "01_FIGURE_DATA_CSV/README_图数据包说明.md；图面与图注一律用 Cuproptosis。",
}


def read_panel(fp: str) -> tuple[pd.DataFrame, list[str]]:
    """读面板 CSV：剥离前导 '#' 注释行（注释进 notes sheet，不当表头）。"""
    with open(fp, encoding="utf-8-sig", newline="") as fh:
        lines = fh.read().splitlines()
    notes, body = [], []
    for ln in lines:
        if ln.startswith("#") and not body:
            notes.append(ln)
        else:
            body.append(ln)
    return pd.read_csv(io.StringIO("\n".join(body)), low_memory=False), notes


def preserve_bool_text(df: pd.DataFrame) -> pd.DataFrame:
    """布尔列按 CSV 原文写出（Excel 往返会把 True/False 漂成 1.0/0.0）。"""
    for c in df.columns:
        s = df[c]
        if s.dtype == bool:
            df[c] = s.astype(object).map({True: "True", False: "False"})
        elif s.dtype == object:
            vals = set(s.dropna().astype(str).unique())
            if vals and vals <= {"True", "False"}:
                df[c] = s.astype(object).map({True: "True", False: "False"}).fillna("")
    return df


def panel_sort_key(ms_panel: str) -> tuple:
    m = re.match(r"([A-Za-z]+)(\d*)([a-z']*)", str(ms_panel))
    return (m.group(1), m.group(2), m.group(3)) if m else (str(ms_panel), "", "")


def main() -> int:
    rows = list(csv.DictReader(open(XWALK, encoding="utf-8-sig")))
    print(f"[source_data] crosswalk 行数 = {len(rows)}")

    buckets: dict[tuple[str, str], list[dict]] = {}
    missing, no_panel, external = [], 0, []
    for r in rows:
        ms_fig = (r.get("ms_figure") or "").strip()
        ms_panel = (r.get("ms_panel") or "").strip()
        rel = (r.get("csv_relpath") or "").strip()
        if not ms_fig:
            continue
        if not ms_panel:
            no_panel += 1
            continue
        if ms_fig.startswith("Figure_S"):
            folder = "Supplementary"
        elif ms_fig.startswith("Figure_"):
            folder = "Main"
        else:
            continue
        # `ext:` = 源在 01_FIGURE_DATA_CSV 之外（正文 Figure 8 的湿实验数据包
        # E:\SCI\Figure_8）。这类行不进本工作簿：它们的 Source Data 需在
        # Figure 8 进正文时另行产出。**显式登记后跳过**，不能当作路径错误，
        # 也不能当成 01_FIGURE_DATA_CSV 下的相对路径去 join。
        if rel.startswith(EXT_PREFIX):
            external.append((ms_fig, ms_panel, rel))
            continue
        fp = os.path.join(SRC, rel.replace("/", os.sep))
        if not os.path.exists(fp):
            missing.append((ms_fig, ms_panel, rel))
            continue
        buckets.setdefault((folder, ms_fig), []).append(
            {"panel": ms_panel, "path": fp, "pkg": (r.get("pkg_panel") or "").strip(),
             "rel": rel})

    # 同一面板号重复指向同一个 CSV 时去重（例如某些面板共用一张源表）
    made = []
    for (folder, ms_fig), items in sorted(buckets.items()):
        items.sort(key=lambda d: panel_sort_key(d["panel"]))
        seen: set[str] = set()
        outdir = os.path.join(OUT, folder)
        os.makedirs(outdir, exist_ok=True)
        out = os.path.join(outdir, f"{ms_fig}_SourceData.xlsx")
        with pd.ExcelWriter(out, engine="openpyxl") as w:
            # 索引页：面板 → 数据包编号 / 相对路径（三套图号世代的对照就在这一页）
            # + 口径说明列（m5：拼写等已在册的说明，避免写进 CSV 破坏表头解析）
            pd.DataFrame([{"panel": it["panel"], "package_panel": it["pkg"],
                           "csv_relpath": it["rel"],
                           "note": INDEX_NOTES.get((folder, ms_fig, it["panel"]),
                                                   "")}
                          for it in items]).to_excel(
                w, sheet_name="Index", index=False)
            for it in items:
                key = it["panel"]
                if key in seen:
                    continue
                seen.add(key)
                df, notes = read_panel(it["path"])
                df = preserve_bool_text(df)
                df.to_excel(w, sheet_name=f"Panel_{key}", index=False)
                if notes:
                    pd.DataFrame({"note": notes}).to_excel(
                        w, sheet_name=f"Panel_{key}_notes", index=False)
                # 同源附属小表（如 Figure_S2E_fit_stats.csv）：并进同一工作簿，
                # 各自成页 —— 审稿复核 m5 指出面板表里混入过量纲不同的汇总行，
                # 拆出后既不丢数据，也不污染面板表的单一量纲。
                stem = os.path.splitext(os.path.basename(it["path"]))[0]
                for extra in sorted(glob.glob(os.path.join(
                        os.path.dirname(it["path"]), f"{stem}_*.csv"))):
                    suffix = os.path.basename(extra)[len(stem) + 1:-4]
                    edf, enotes = read_panel(extra)
                    pd.DataFrame(edf).to_excel(
                        w, sheet_name=f"Panel_{key}_{suffix}", index=False)
                    if enotes:
                        pd.DataFrame({"note": enotes}).to_excel(
                            w, sheet_name=f"Panel_{key}_{suffix}_notes",
                            index=False)
                    print(f"  [extra] {os.path.relpath(extra, ROOT)} -> "
                          f"{os.path.basename(out)} : Panel_{key}_{suffix}")
        made.append((out, len(seen)))

    for p, n in made:
        print(f"  {os.path.relpath(p, ROOT)}  <- {n} panels, "
              f"{os.path.getsize(p) // 1024} KB")
    if no_panel:
        print(f"[source_data] 跳过 {no_panel} 行（crosswalk 无 ms_panel）")
    if external:
        figs = sorted({f for f, _p, _r in external})
        print(f"[source_data] 跳过 {len(external)} 行外部源（{', '.join(figs)}；"
              f"csv_relpath 带 ext: 前缀）：其源数据不在 {os.path.relpath(SRC, ROOT)} "
              f"下，Source Data 需在 Figure 8 进正文时另行产出")
        for f, p, rel in external:
            print(f"    {f} {p} -> {rel}")
    if missing:
        print("[source_data] !! crosswalk 指向的 CSV 不存在：")
        for f, p, rel in missing:
            print(f"    {f} {p} -> {rel}")
        return 1
    print(f"[source_data] 共 {len(made)} 个 xlsx 已写入 {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
