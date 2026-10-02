# -*- coding: utf-8 -*-
"""2026-09-05 CDD 投稿适配：由 01_FIGURE_DATA_CSV 逐图生成组合 Excel 源数据文件
（Nature Portfolio / Cell Death & Differentiation Source Data 惯例：每图一个 xlsx、每面板一个 sheet）。
数值零改动；QC 修订版（当日下午严格质控后）：
  1) 布尔列保真：CSV 中 'True'/'False' 文本列在 Excel 往返中会漂移为 1.0/0.0 浮点
     （Figure_8A sex_included、Figure_8D same_sign、Figure_S8E heidi_pass），现固定写字符串；
  2) 注释首行处理：Figure_S12D.csv 首行为 '#' 注释行（非表头）且 CRLF——原版误把注释当表头
     导致 sheet 损坏（1×1），现剥离注释行入 Panel_X_notes sheet、其余行按真表头解析。
"""
import os, re, glob, io
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(ROOT, "01_FIGURE_DATA_CSV")
OUT = os.path.join(ROOT, "投稿CELL DEATH AND DIFFERENTIATION", "Source_Data")

def panel_key(fp):
    m = re.search(r"Figure_S?\d+([A-Z]+)\.csv$", os.path.basename(fp))
    return m.group(1)

def read_panel(fp):
    """读面板 CSV：剥离前导 '#' 注释行（返回 df, notes）。"""
    with open(fp, encoding="utf-8-sig", newline="") as f:
        raw = f.read()
    lines = raw.splitlines()
    notes, body = [], []
    for ln in lines:
        if ln.startswith("#") and not body:
            notes.append(ln)
        else:
            body.append(ln)
    df = pd.read_csv(io.StringIO("\n".join(body)), low_memory=False)
    return df, notes

def preserve_bool_text(df):
    """布尔列按 CSV 原文写出（避免 Excel 往返 True→1.0 漂移）。"""
    for c in df.columns:
        s = df[c]
        if s.dtype == bool:
            df[c] = s.astype(object).map({True: "True", False: "False"})
        elif s.dtype == object:
            vals = set(s.dropna().astype(str).unique())
            if vals and vals <= {"True", "False"}:
                df[c] = s.astype(object).map({True: "True", False: "False"}).fillna("")
    return df

def build(folder, prefix):
    indir = os.path.join(SRC, folder)
    outdir = os.path.join(OUT, folder)
    os.makedirs(outdir, exist_ok=True)
    groups = {}
    for fp in glob.glob(os.path.join(indir, f"{prefix}*.csv")):
        m = re.match(rf"{prefix}(\d+)([A-Z]+)\.csv$", os.path.basename(fp))
        if m:
            groups.setdefault(int(m.group(1)), []).append(fp)
    made = []
    for num in sorted(groups):
        files = sorted(groups[num], key=panel_key)
        tag = f"{prefix}{num}"
        out = os.path.join(outdir, f"{tag}_SourceData.xlsx")
        with pd.ExcelWriter(out, engine="openpyxl") as w:
            for fp in files:
                letter = panel_key(fp)
                df, notes = read_panel(fp)
                df = preserve_bool_text(df)
                df.to_excel(w, sheet_name=f"Panel_{letter}", index=False)
                if notes:
                    pd.DataFrame({"note": notes}).to_excel(
                        w, sheet_name=f"Panel_{letter}_notes", index=False)
        made.append((out, len(files)))
    return made

def main():
    all_made = build("Main", "Figure_") + build("Supplementary", "Figure_S")
    t4 = os.path.join(SRC, "Table4_criteria_verdicts.csv")
    if os.path.exists(t4):
        out = os.path.join(OUT, "Main", "Table_4_SourceData.xlsx")
        with pd.ExcelWriter(out, engine="openpyxl") as w:
            df, _ = read_panel(t4)
            preserve_bool_text(df).to_excel(w, sheet_name="Table_4", index=False)
        all_made.append((out, 1))
    for p, n in all_made:
        print(f"{os.path.relpath(p, ROOT)}  <- {n} panels, {os.path.getsize(p)//1024} KB")

if __name__ == "__main__":
    main()
