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
CN = io.open(os.path.join(ROOT, "111文稿_v5.md"), encoding="utf-8").read()
OUT = os.path.join(ROOT, "03_LOGS", "P13_crossref_figure_check_20260901.md")
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
expect = {"Figure 6": "Section 10", "Figure 7": "Section 11", "Figure 8": "Section 12",
          "Figure 9": "Section 13", "Figure 10": "Section 14", "Figure 11": "Section 15"}
i = EN.index("# Figure Legends"); j = EN.index("# References")
leg_en = EN[i:j]
log("## EN 图注 (Section N) 复核")
for fig, sec in expect.items():
    m = re.search(re.escape(fig) + r"[. ].{0,400}?\((?:见 )?" + sec.replace("Section", "Section ") + r"\)", leg_en, re.S)
    m2 = re.search(re.escape(fig) + r"[. ].{0,400}?\(Section (\d+)\)", leg_en, re.S)
    found = m2.group(1) if m2 else None
    ok = "OK" if found == sec.split()[-1] else "MISMATCH"
    log(f"- {fig}: 图注 Section={found} 期望={sec.split()[-1]} -> {ok}")
# CN 图注 第N节
i = CN.index("# Figure Legends"); j = CN.index("# References")
leg_cn = CN[i:j]
expect_cn = {"Figure 6": "10", "Figure 7": "11", "Figure 8": "12",
             "Figure 9": "13", "Figure 10": "14", "图 11": "15"}
log("## CN 图注（第 N 节）复核")
for fig, sec in expect_cn.items():
    m2 = re.search(re.escape(fig) + r"[. ].{0,400}?第\s*(\d+)\s*节", leg_cn, re.S)
    found = m2.group(1) if m2 else None
    ok = "OK" if found == sec else "MISMATCH"
    log(f"- {fig}: 图注 第{found}节 期望={sec} -> {ok}")
# S11/S12 由 §16 引用
m = re.search(r"### 16\.[^\n]*S11", EN[EN.index("## Part VI"):EN.index("# Discussion")])
log(f"- EN §16 标题引用 Figure S11: {'OK' if m else 'CHECK'}")
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
for name, leg, pref in [("EN", leg_en, None), ("CN", leg_cn, None)]:
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
