# -*- coding: utf-8 -*-
"""
M1 Step 10: lint 收尾登记修正（2026-08-17）
================================================================
M1 交付后 lint 全包校验暴露的登记层缺口，一次性修正（只改登记元数据与表尾
版本列，不改任何数值列）：

1. 21 行 S55 manifest 记录 location 误为 '02_SUPPLEMENTARY_TABLES'，
   实际文件在 '02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV' → 修正；
2. 19 张 S55 表缺 gene_set_version/score_version 列（M1_README §5 承诺
   "全部新表带版本列"，与 lint check_versions 口径冲突）→ 按 manifest
   常量值补列（S55a/S55c 已有，跳过）；
3. 未登记文件补登记：M1_spatial_sample_manifest.csv（GOV）、
   rawdata_cleanup_20260817.log/.py（GOV）、archive/ 下 13 个 M1 脚本
   （archived）；
4. 对全部改动文件重算 size_bytes + sha256 并重写 manifest（重冻结）。

执行后运行 lint_package.py 验证全绿。
"""
import csv
import hashlib
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GOV = ROOT / "04_AUDIT_GOVERNANCE"
TAB = ROOT / "02_SUPPLEMENTARY_TABLES" / "SUPPLEMENTARY_Tables_CSV"
ARCH = ROOT / "archive"
MANIFEST = GOV / "RESULTS_MANIFEST_v1.0.csv"

TAB_LOC = "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV"
GOV_LOC = "04_AUDIT_GOVERNANCE"


def sha256_of(p, block=1024 * 1024):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while True:
            b = f.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def read_rows_csv(p):
    raw = p.read_bytes()
    enc = "utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf-8"
    with open(p, "r", encoding=enc, newline="") as f:
        return [r for r in csv.reader(f)]


def write_rows_csv(p, rows):
    with open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerows(rows)


def append_version_cols(p, gsv, sv):
    """表尾追加 gene_set_version/score_version 常量列（不改动已有单元格）。"""
    rows = read_rows_csv(p)
    assert rows, f"empty file: {p}"
    hdr = rows[0]
    if "gene_set_version" in hdr and "score_version" in hdr:
        return False
    rows[0] = hdr + ["gene_set_version", "score_version"]
    for i in range(1, len(rows)):
        rows[i] = rows[i] + [gsv, sv]
    write_rows_csv(p, rows)
    return True


def main():
    man = read_rows_csv(MANIFEST)
    hdr, mrows = man[0], man[1:]
    idx = {c: hdr.index(c) for c in hdr}
    by_fn = {r[idx["filename"]]: r for r in mrows}

    # ---- 1+2: S55 location 修正 + 19 表补版本列 + 重算哈希 ----
    s55 = [r for r in mrows if r[idx["s_number"]] == "S55"]
    assert len(s55) == 21, f"S55 rows = {len(s55)}"
    changed_files = []
    for r in s55:
        fn = r[idx["filename"]]
        r[idx["location"]] = TAB_LOC
        p = TAB / fn
        assert p.exists(), f"missing {p}"
        if append_version_cols(p, r[idx["gene_set_version"]], r[idx["score_version"]]):
            changed_files.append((p, r, fn))
    print(f"[S55] location 修正 21 行；补版本列 {len(changed_files)} 张")

    # ---- 3a: M1_spatial_sample_manifest.csv 补空版本列 + 登记 ----
    msp = GOV / "M1_spatial_sample_manifest.csv"
    append_version_cols(msp, "", "")
    n_msp = sum(1 for _ in open(msp, encoding="utf-8")) - 1
    assert n_msp == 139, f"M1 sample manifest rows = {n_msp}"

    # ---- 3b: 未登记文件补登记 ----
    new_rows = []
    reg = set(by_fn)

    def add_row(fn, loc, cls, status, analysis, notes, universe="none", gsv="", sv="", s_number=""):
        p = ({"04_AUDIT_GOVERNANCE": GOV, "archive": ARCH,
              "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV": TAB}[loc]) / fn
        n_lines = sum(1 for _ in open(p, "rb"))
        n_data = n_lines - 1 if fn.lower().endswith((".csv",)) else n_lines
        return [s_number, fn, loc, cls, status, universe, gsv, sv, str(n_data),
                str(p.stat().st_size), sha256_of(p), analysis, notes]

    if "M1_spatial_sample_manifest.csv" not in reg:
        new_rows.append(add_row("M1_spatial_sample_manifest.csv", GOV_LOC, "analysis_output",
                                "canonical", "M1 空间数据样本清单（GSE271370 23 切片 + GSE253474 116 FOV）", ""))
    for fn, cls, analysis, notes in [
        ("rawdata_cleanup_20260817.py", "audit_output", "2026-08-17 原始数据清理脚本", ""),
        ("rawdata_cleanup_20260817.log", "audit_report", "2026-08-17 原始数据清理记录", ""),
    ]:
        if fn not in reg:
            new_rows.append(add_row(fn, GOV_LOC, cls, "canonical", analysis, notes))
    for fn in sorted(p.name for p in ARCH.iterdir() if p.is_file() and p.name.startswith("M1_")
                     and p.name not in reg):
        new_rows.append(add_row(fn, "archive", "analysis_output", "archived", "",
                                "M1 管线辅助/检查脚本，2026-08-17 归档（见 M1_README_分析流程.md）"))
    print(f"[登记] 新增 manifest 行 {len(new_rows)}")

    # ---- 4: 更新改动文件哈希 + 写回 manifest ----
    for p, r, fn in changed_files:
        r[idx["size_bytes"]] = str(p.stat().st_size)
        r[idx["sha256"]] = sha256_of(p)
        print(f"[SHA256] 重冻结 {fn}")
    man = [hdr] + mrows + new_rows
    write_rows_csv(MANIFEST, man)
    print(f"manifest 重写完成：{len(mrows)} + {len(new_rows)} = {len(man) - 1} 行")


if __name__ == "__main__":
    main()
