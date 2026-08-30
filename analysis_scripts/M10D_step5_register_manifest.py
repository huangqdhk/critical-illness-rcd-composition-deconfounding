# -*- coding: utf-8 -*-
"""
M10D_step5_register_manifest.py — M10D 收尾批次 manifest 登记
================================================================
登记 Table S83b（猪同源敏感性）+ Table S87/S87b/S87c（k≥20 泛疾病扩展）
到 RESULTS_MANIFEST_v2.0.csv（追加，不改动既有行），并输出 lint 白名单
需要的 M1X_NO_VERSION_STAMP 增补清单（由本脚本打印，人工并入 lint_package.py）。
"""
import csv
import hashlib
import os
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
MAN = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "RESULTS_MANIFEST_v2.0.csv")
TAB = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES", "SUPPLEMENTARY_Tables_CSV")

NEW = [
    ("S83b", "Table_S83b_M10D_Pig_Ensembl_Sensitivity.csv",
     "M10D 猪层同源敏感性复算（Ensembl ortholog_one2one 49 条 vs 降级注释逐字口径）",
     "2026-08-27 M10D 补做批次（Ensembl 服务恢复后 one2one 复算；方向与主口径一致）"),
    ("S87", "Table_S87_M10D_Pandisease_k20_LayerEffects.csv",
     "M10D 泛疾病扩展 k=20 层效应表（冻结 16 层 + 新增 4 层；MDI Hedges' g）",
     "2026-08-27 M10D 补做批次（注册 ETVMJ §1.5 剩余项；新增层来源与限制逐行披露）"),
    ("S87b", "Table_S87b_M10D_Pandisease_k20_Meta.csv",
     "M10D 泛疾病扩展合并表（k=20 主 + k=19/k=18 敏感性 + k=16 锚点重算）",
     "2026-08-27 M10D 补做批次（REML+Hartung-Knapp+PI；k=16 锚点 +1.504 与冻结 +1.50 互核）"),
    ("S87c", "Table_S87c_M10D_Pandisease_k20_LOSO.csv",
     "M10D 泛疾病扩展 LOSO 留一法（k=20 主分析）",
     "2026-08-27 M10D 补做批次（LOSO 20/20 同向）"),
]


def sha256_of(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def n_rows(p):
    with open(p, "r", encoding="utf-8", newline="") as f:
        return sum(1 for _ in csv.reader(f)) - 1  # 减表头


def main():
    rows = []
    with open(MAN, "r", encoding="utf-8", newline="") as f:
        rd = csv.DictReader(f)
        cols = rd.fieldnames
        rows = list(rd)
    existing = {r["filename"] for r in rows}
    for s_no, fn, desc, note in NEW:
        if fn in existing:
            print("skip (already):", fn)
            continue
        p = os.path.join(TAB, fn)
        rows.append({
            "s_number": s_no, "filename": fn,
            "location": "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
            "class": "supplementary_table", "status": "canonical",
            "gene_universe": "none", "gene_set_version": "", "score_version": "",
            "n_data_rows": n_rows(p), "size_bytes": os.path.getsize(p),
            "sha256": sha256_of(p), "analysis": desc, "notes": note, "class_": "",
        })
        print("registered:", fn, rows[-1]["n_data_rows"], "rows")
    with open(MAN, "w", encoding="utf-8", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=cols)
        wr.writeheader()
        wr.writerows(rows)
    print("manifest rows now:", len(rows))
    print()
    print("M1X_NO_VERSION_STAMP additions (paste into lint_package.py):")
    for s_no, fn, _d, _n in NEW:
        print('    "%s",' % fn)


if __name__ == "__main__":
    main()
