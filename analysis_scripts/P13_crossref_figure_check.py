# -*- coding: utf-8 -*-
"""
P13_crossref_figure_check.py — 111提质计划13 §三/§十一：交叉引用与图件一致性检查（只读）
======================================================================================
检查项：
1. 英文稿 `Section N(.M)` / `§N(.M)` 引用 ↔ Results ### N. 小节标题（悬空引用检测）
2. 中文稿 `第 N 节` / `§N.M` 引用 ↔ Results ### N. 小节标题
3. 图注复核清单（提质计划 §3.3）：Figure 6–11 图注中的 (Section N) 与新编号一致
4. 主图/补图 panel ↔ 01_FIGURE_DATA_CSV/ 下 CSV 文件与 FIGURE_DATA_MANIFEST.csv 登记一致性
5. 中英文稿各自 Results 小节序列完整性（1–17 无缺/无重）
输出：03_LOGS/P13_crossref_figure_check_20260901.md（仅报告，不改任何文件）
"""
import io, os, re, sys, csv
from collections import OrderedDict

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = os.path.dirname(os.path.abspath(__file__))
EN = io.open(os.path.join(ROOT, "111文稿_v5_英文版.md"), encoding="utf-8").read()
# 2026-09-05：中文稿已随 2026-09-04 归档批次移入 归档/（弃用只读），路径同步；输出加日期不覆盖 09-01 报告
# 2026-09-21：归档区分类整理，中文稿移入 归档/04_文稿与计划备份/，路径再同步。
CN = io.open(os.path.join(ROOT, "归档", "04_文稿与计划备份", "111文稿_v5.md"),
             encoding="utf-8").read()
OUT = os.path.join(ROOT, "03_LOGS", "P13_crossref_figure_check_20260916.md")
rep = []
def log(s=""): rep.append(s); print(s)

def results_headings(txt, marker):
    """提取 Results 与 Discussion 之间的 ### N. 小节号"""
    i = txt.index("# Results"); j = txt.index("# Discussion")
    body = txt[i:j]
    return sorted(set(int(m.group(1)) for m in re.finditer(r"^### (\d+)\. ", body, re.M)))

# ---------- 1/2/5 小节引用 vs 标题 ----------
for name, txt, pats in [
    ("EN", EN, [r"Section (\d+)(?:\.(\d+))?", r"§(\d+)(?:\.(\d+))?"]),
    ("CN", CN, [r"第\s*(\d+)\s*节", r"§(\d+)(?:\.(\d+))?"]),
]:
    heads = results_headings(txt, name)
    log(f"## {name} Results 小节序列: {heads}")
    missing = [n for n in range(1, 18) if n not in heads]
    log(f"- 序列完整性: {'OK (1–17 全)' if not missing else '缺 ' + str(missing)}")
    refs = set()
    for pat in pats:
        for m in re.finditer(pat, txt):
            top = int(m.group(1))
            sub = m.group(2) if m.lastindex >= 2 else None
            refs.add((top, sub))
    dangling = sorted({t for t, s in refs if t not in heads})
    log(f"- 顶层引用 {sorted({t for t,s in refs})}; 悬空顶层引用: {dangling if dangling else '无'}")
    # 小节级引用粗查（11.x 是否落在 §11 的已知子节 1–7）
    subrefs = sorted({(t, int(s)) for t, s in refs if s})
    log(f"- 小节级引用: {subrefs}")
log("")

# ---------- 3 图注 Section 复核 ----------
# 2026-09-16 A+B 压缩后：仅新 Figure 5（旧 6+7 合并）图注带 "(Sections 10–11)" 标签
# 修复：本文档 # References 在 # Figure Legends 之前，原切片 EN[i:j] 得空串——改为切到文末
expect = {"Figure 5": "Sections 10–11"}
i = EN.index("# Figure Legends")
leg_en = EN[i:]
log("## EN 图注 (Section N) 复核")
for fig, sec in expect.items():
    m2 = re.search(re.escape(fig) + r"[. ].{0,400}?\((Sections? [\d–—-]+)\)", leg_en, re.S)
    found = m2.group(1) if m2 else None
    ok = "OK" if found == sec else "MISMATCH"
    log(f"- {fig}: 图注 Section={found} 期望={sec} -> {ok}")
# CN 图注 第N节（CN 稿为归档只读旧图号快照——2026-09-16 A+B 压缩未同步中文归档版，仅复核其冻结状态的两条有效标签）
i = CN.index("# Figure Legends")
leg_cn = CN[i:]
expect_cn = {"Figure 6": "10", "Figure 7": "11"}  # 归档 CN 稿（2026-09-06 九图体系）仅存此两条节标签
log("## CN 图注（第 N 节）复核（归档旧编号，仅供参考）")
for fig, sec in expect_cn.items():
    m2 = re.search(re.escape(fig) + r"[. ].{0,400}?第\s*(\d+)\s*节", leg_cn, re.S)
    found = m2.group(1) if m2 else None
    ok = "OK" if found == sec else "MISMATCH"
    log(f"- {fig}: 图注 第{found}节 期望={sec} -> {ok}")
# S12/S13 由 §16 引用（2026-09-16 A+B 附图级联：旧 S10→S12、旧 S11→S13）
m = re.search(r"### 16\.[^\n]*S13", EN[EN.index("## Part VI"):EN.index("# Discussion")])
log(f"- EN §16 标题引用 Figure S13: {'OK' if m else 'CHECK'}")
log("")

# ---------- 4 panel ↔ CSV ↔ manifest ----------
figdir = os.path.join(ROOT, "01_FIGURE_DATA_CSV")
files = set()
for dp, dn, fn in os.walk(figdir):
    for f in fn:
        if f.endswith(".csv") and f.startswith("Figure_"):
            files.add(f[:-4])
man = {}
with io.open(os.path.join(figdir, "FIGURE_DATA_MANIFEST.csv"), encoding="utf-8-sig", newline="") as f:
    for row in csv.DictReader(f):
        man[row["panel"]] = row["package_file"]
log("## panel ↔ CSV ↔ manifest")
missing_file, missing_man = [], []
# 2026-09-16：CN 稿为归档旧图号快照（弃用只读），其面板引用不再对照新编号 CSV——仅核查 EN 终稿
for name, leg, pref in [("EN", leg_en, None)]:
    for m in re.finditer(r"\*\*(Figure |图 )(S?\d+)[. ]", leg):
        num = m.group(2)
        # 抓该图注行内的 panel 字母
        seg = leg[m.start(): leg.find("\n\n", m.start()) if leg.find("\n\n", m.start()) > 0 else len(leg)]
        letters = sorted(set(re.findall(r"\(([A-Z])(?=[,)–—])", seg)))
        for L in letters:
            key = f"Figure_{num}{L}"
            if key not in files and f"S{num}{L}" not in files and key not in man:
                missing_file.append(f"{name}:{key}")
log(f"- manifest 登记 panel 数: {len(man)}; 磁盘 Figure_* CSV 数: {len(files)}")
log(f"- 图注提及但无 CSV/manifest 登记的 panel: {missing_file if missing_file else '无（或均为示意面板）'}")
# 反向：manifest 有登记但磁盘缺文件
rev = [p for p, pkg in man.items()
       if p.startswith("Figure") and not os.path.exists(os.path.join(figdir, pkg))]
log(f"- manifest 登记但磁盘缺文件: {rev if rev else '无'}")

io.open(OUT, "w", encoding="utf-8", newline="\n").write("\n".join(rep))
print("\nreport ->", OUT)
