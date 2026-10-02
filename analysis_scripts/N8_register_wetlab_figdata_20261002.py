# -*- coding: utf-8 -*-
"""
N8_register_wetlab_figdata_20261002.py — 补齐 2026-10-02 两项登记缺口

背景：N6（湿实验三图入包）只更新了 01_FIGURE_DATA_CSV/FIGURE_DATA_MANIFEST.csv，
未把 21 个新图源文件登记进 04_AUDIT_GOVERNANCE/RESULTS_MANIFEST_v2.0.csv；
另 归档/ 顶层留有 2026-09-28 的 Table_S26k 留档副本未登记。本脚本补齐：

  1) 登记 21 个湿实验图源（canonical）：
     Main/Figure_8wA–8wH（8）+ Supplementary/Figure_S15A–G（7）+ Figure_S16A–F（6）
  2) 登记 归档/Table_S26k_CellType_Scores.csv（archived；canonical 已于 N7 重编号为 S79k）
  3) lint_package.py：把 21 个新文件加入 M1X_NO_VERSION_STAMP 豁免
     （湿实验源包 CSV 无 gene_set_version/score_version 列，manifest 记 NA——设计内豁免）

幂等：同名行先删后增；重复运行结果一致。
用法：python N8_register_wetlab_figdata_20261002.py
"""
import csv, hashlib, io, os, sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
MAN = os.path.join(HERE, "04_AUDIT_GOVERNANCE", "RESULTS_MANIFEST_v2.0.csv")
LINT = os.path.join(HERE, "lint_package.py")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def nrows(p):
    with io.open(p, encoding="utf-8-sig", newline="") as f:
        return sum(1 for _ in f) - 1


WET = []
for L in "ABCDEFGH":
    WET.append(("Main", "Figure_8w%s.csv" % L, "8%s" % L, "Figure_8",
                "Figure_8\\Figure_8\\Figure_8%s.csv" % L))
for L in "ABCDEFG":
    WET.append(("Supplementary", "Figure_S15%s.csv" % L, "S15%s" % L, "Figure_S15",
                "Figure_S15-16待可视化数据\\Figure_S15%s.csv" % L))
for L in "ABCDEF":
    WET.append(("Supplementary", "Figure_S16%s.csv" % L, "S16%s" % L, "Figure_S16",
                "Figure_S15-16待可视化数据\\Figure_S16%s.csv" % L))

print("== 1) 登记 %d 个湿实验图源（canonical） ==" % len(WET))
with io.open(MAN, encoding="utf-8-sig", newline="") as f:
    rows = list(csv.DictReader(f))
    fields = list(rows[0].keys())

new_names = {fn for _s, fn, _p, _fig, _src in WET} | {"Table_S26k_CellType_Scores.csv"}
body = [r for r in rows if r["filename"] not in new_names]

for sub, fn, panel, fig, src in WET:
    p = os.path.join(HERE, "01_FIGURE_DATA_CSV", sub, fn)
    body.append({
        "s_number": "", "filename": fn, "location": "01_FIGURE_DATA_CSV/" + sub,
        "class": "figure_data", "status": "canonical", "gene_universe": "none",
        "gene_set_version": "NA", "score_version": "NA",
        "n_data_rows": str(nrows(p)), "size_bytes": str(os.path.getsize(p)),
        "sha256": sha(p),
        "analysis": "ms_panel=%s（湿实验 %s）；2026-10-02 N6 湿实验三图入包" % (panel, fig),
        "notes": "源包 %s（逐字节复制、MD5 断言见 N6_integrate_wetlab_figdata_20261002.py）" % src,
        "class_": ""})

print("== 2) 登记 归档/Table_S26k_CellType_Scores.csv（archived） ==")
p = os.path.join(HERE, "归档", "Table_S26k_CellType_Scores.csv")
body.append({
    "s_number": "S26", "filename": "Table_S26k_CellType_Scores.csv", "location": "归档",
    "class": "supplementary_table", "status": "archived", "gene_universe": "none",
    "gene_set_version": "NA", "score_version": "NA",
    "n_data_rows": str(nrows(p)), "size_bytes": str(os.path.getsize(p)),
    "sha256": sha(p),
    "analysis": "Zenodo 上传件归档留档（2026-09-28）",
    "notes": "canonical 件 2026-10-02 重编号为 Table_S79k_CellType_Scores.csv；本副本保留原名供溯源",
    "class_": ""})

with io.open(MAN, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, lineterminator="\r\n")
    w.writeheader()
    w.writerows(body)
print("   manifest 行数 %d → %d" % (len(rows), len(body)))

print("== 3) lint 豁免：M1X_NO_VERSION_STAMP += 21 湿实验文件 ==")
src = open(LINT, encoding="utf-8").read()
if "2026-10-02 N6 湿实验三图入包" not in src:
    anchor = "    \"Figure_7D.csv\",\n"
    assert src.count(anchor) == 1, "锚点不唯一"
    block = ("    # 2026-10-02 N6 湿实验三图入包：源包 CSV 无 gene_set_version/score_version 列，\n"
             "    # manifest 记 NA（设计内豁免；同上 M10–M16 导出设计）\n")
    for _s, fn, _p, _fig, _src in WET:
        block += '    "%s",\n' % fn
    src = src.replace(anchor, anchor + block, 1)
    open(LINT, "w", encoding="utf-8", newline="").write(src)
    print("   已写入豁免块（21 项）")
else:
    print("   豁免块已存在（跳过）")

print("== 完成 ==")
