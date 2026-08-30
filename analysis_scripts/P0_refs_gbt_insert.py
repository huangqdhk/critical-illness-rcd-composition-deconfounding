# -*- coding: utf-8 -*-
"""
P0_refs_gbt_insert.py — GB/T 7714-2025 参考文献落盘（2026-08-29）
功能：
  1. 按首次出现序重编号（数字标记 + 4 条作者-年份 + 3 处 prose 引用 + 2 条新增条目）
  2. 双版标记替换（中文/英文，保护行尾与 p∈[0,1]）
  3. 生成 GB/T 7714-2025 参考文献表（[J] 期刊 / [DS/OL] 数据集 / [PP/OL] 预印本 / [DB/OL] 数据库）
  4. 替换双版 References 占位节
  5. 一致性终检
"""
import json, re, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
ZH = BASE + r"\111文稿_v5.md"
EN = BASE + r"\111文稿_v5_英文版.md"
META = json.load(open(BASE + r"\_intermediate\P0_refs_gbt_metadata.json", encoding="utf-8"))
CITE_DATE = "2026-08-29"

# ---------- 1. 首次出现序（按行号+列号扫描中文版；prose 事件按行号注入） ----------
prose_events = {138: "Kapp", 158: "Yao", 168: "Overmyer"}  # 行号 -> 槽名
zh_text = open(ZH, encoding="utf-8", newline="").read()
zh_lines = zh_text.split("\n")

events = []  # (line, col, slot)
for i, line in enumerate(zh_lines, 1):
    for m in re.finditer(r"\[[\d,\-\s]+\]", line):
        if m.group(0).strip() == "[0,1]":
            continue
        for part in m.group(0)[1:-1].split(","):
            part = part.strip()
            if "-" in part:
                a, b = part.split("-")
                for n in range(int(a), int(b) + 1):
                    events.append((i, m.start(), str(n)))
            elif part:
                events.append((i, m.start(), part))
    for m in re.finditer(r"\[Galluzzi et al\. 2018\]", line):
        events.append((i, m.start(), "Galluzzi"))
    for m in re.finditer(r"\[Wang 2025[；;][^\]]*\]|\[Wang 2025\]", line):
        slot = "Wang" if m.group(0) in ("[Wang 2025]",) else "Wang+Casadio"
        if slot == "Wang":
            events.append((i, m.start(), "Wang"))
        else:
            events.append((i, m.start(), "Wang"))
            events.append((i, m.start() + 2, "Casadio"))
    for m in re.finditer(r"\[PMID 41317732\]", line):
        events.append((i, m.start(), "Wang"))  # 已编号则不重复计
    if i in prose_events:
        events.append((i, 0, prose_events[i]))

order, mapping = [], {}
for _, _, slot in sorted(events, key=lambda x: (x[0], x[1])):
    if slot not in mapping:
        mapping[slot] = len(order) + 1
        order.append(slot)
print("首次出现序（%d 条）：" % len(order))
print("  " + ", ".join(f"{s}→{mapping[s]}" for s in order))

# ---------- 2. 标记替换 ----------
def compress(nums):
    nums = sorted(set(nums))
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(str(nums[i]) if i == j else f"{nums[i]}-{nums[j]}")
        i = j + 1
    return "[" + ",".join(out) + "]"

def rewrite_numeric(text):
    def cb(m):
        if m.group(0).strip() == "[0,1]":
            return m.group(0)
        nums = []
        for part in m.group(0)[1:-1].split(","):
            part = part.strip()
            if "-" in part:
                a, b = part.split("-")
                nums += [mapping[str(n)] for n in range(int(a), int(b) + 1)]
            elif part:
                nums.append(mapping[part])
        return compress(nums)
    return re.sub(r"\[[\d,\-\s]+\]", cb, text)

REPL_ZH = [
    ("[Galluzzi et al. 2018]", "[%d]" % mapping["Galluzzi"]),
    ("[Wang 2025；Casadio *Nat Cell Biol* 2026, PMID 41540269]",
     "[%d,%d]" % (mapping["Wang"], mapping["Casadio"])),
    ("[Wang 2025]", "[%d]" % mapping["Wang"]),
    ("[PMID 41317732]", "[%d]" % mapping["Wang"]),
    ("Kapp & Tibshirani 2007", "Kapp & Tibshirani 2007 [%d]" % mapping["Kapp"]),
    ("Yao et al. bioRxiv 2023, PMID 36747806", "Yao et al.[%d]" % mapping["Yao"]),
    ("Yao et al., bioRxiv 2023, PMID 36747806", "Yao et al.[%d]" % mapping["Yao"]),
    ("Overmyer et al. *Cell Systems* 2020", "Overmyer et al.[%d]" % mapping["Overmyer"]),
]
REPL_EN = REPL_ZH[:4] + [
    ("Kapp & Tibshirani 2007", "Kapp & Tibshirani 2007 [%d]" % mapping["Kapp"]),
    ("Yao et al. bioRxiv 2023, PMID 36747806", "Yao et al.[%d]" % mapping["Yao"]),
    ("Yao et al., bioRxiv 2023, PMID 36747806", "Yao et al.[%d]" % mapping["Yao"]),
    ("Overmyer et al. *Cell Systems* 2020", "Overmyer et al.[%d]" % mapping["Overmyer"]),
]

def apply_all(path, repls):
    t = open(path, encoding="utf-8", newline="").read()
    t = rewrite_numeric(t)  # 先做数字重排（避免作者-年份替换产生的新编号被二次映射）
    for old, new in repls:
        if old not in t:
            print(f"  !! 未找到替换目标: {old[:60]}")
        t = t.replace(old, new)
    open(path, "w", encoding="utf-8", newline="").write(t)
    print(f"  已重写 {path.split(chr(92))[-1]}")

print("标记替换：")
apply_all(ZH, REPL_ZH)
apply_all(EN, REPL_EN)

# ---------- 3. 生成 GB/T 7714-2025 条目 ----------
INST_HINTS = ("Consortium", "Task Force", "Network", "Group", "eQTLGen", "GTEx", "UniProt", "collaboration", "Collaboration")

def fmt_authors(authors):
    out = []
    for a in authors:
        if any(k in a for k in INST_HINTS) or " " not in a.strip():
            out.append(a)
        else:
            parts = a.split()
            out.append(" ".join(parts[:-1]).upper() + " " + parts[-1])
    return ", ".join(out[:3]) + ", et al" if len(out) > 3 else ", ".join(out)

def fmt_journal(m):
    au = fmt_authors(m["authors"])
    vi = f"{m['volume']}({m['issue']})" if m["issue"] else (m["volume"] or m["year"])
    return f"{au}. {m['title']}[J]. {m['journal']}, {m['year']}, {vi}: {m['pages']}."

GEO_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
geo_entries = {
    "18": ("GSE165659", None), "27": ("GSE212865", None), "65": ("GSE271370", None),
    "67": ("GSE310929", None), "69": ("GSE148871", None),
}
def fmt_geo(acc):
    g = META["geo"][acc]
    return (f"NCBI Gene Expression Omnibus. {acc}: {g['title'].rstrip('.')}[DS/OL]. "
            f"({g['pdat']})[{CITE_DATE}]. {GEO_URL}{acc}.")

def entry_for(slot):
    if slot in geo_entries:
        return fmt_geo(geo_entries[slot][0])
    if slot == "74":
        return (f"MetaboLights. MTBLS6844: serum lipidome of COVID-19 patients (FIA-MS positive "
                f"ion mode, 53 samples)[DS/OL]. EMBL-EBI[{CITE_DATE}]. https://www.ebi.ac.uk/metabolights/MTBLS6844.")
    if slot == "77":
        return (f"UniProt Consortium. UniProt: the Universal Protein Knowledgebase[DB/OL]. "
                f"[{CITE_DATE}]. https://www.uniprot.org.")
    if slot == "Yao":
        m = META["pubmed"]["N:Yao"]
        doi = (m.get("doi") or "").replace("doi: ", "").replace("DOI: ", "")
        url = f"https://doi.org/{doi}" if doi else f"https://pubmed.ncbi.nlm.nih.gov/{m['pmid']}/"
        return (f"{fmt_authors(m['authors'])}. {m['title']}[PP/OL]. bioRxiv({m['year']})[{CITE_DATE}]. {url}.")
    key = str(slot) if str(slot) in META["pubmed"] else "A:" + slot if "A:" + slot in META["pubmed"] else ("N:" + slot if "N:" + slot in META["pubmed"] else None)
    if key is None:
        return None
    return fmt_journal(META["pubmed"][key])

lines_out = []
for idx, slot in enumerate(order, 1):
    e = entry_for(slot)
    if e is None:
        print(f"  !! 缺条目: {slot}")
        continue
    lines_out.append(f"[{idx}] {e}")

ref_block = "\n".join(lines_out)
print("\n参考文献表（%d 条）首 3 条与末 2 条：" % len(lines_out))
for l in lines_out[:3] + lines_out[-2:]:
    print("  " + l[:160])

# ---------- 4. 插入双版 References 节 ----------
PLACEHOLDER_ZH = "【占位——参考文献全文见 `E:\\SCI\\REFERENCES_README.md`，待按目标期刊格式插回此处；正文 [1]–[77] 编号与该清单逐条对应。】"
PLACEHOLDER_EN = "【Placeholder—the full reference list is in `E:\\SCI\\REFERENCES_README.md`; to be inserted here in the target journal's format; the in-text [1]–[77] numbering corresponds one-to-one to that list.】"

for path, ph in ((ZH, PLACEHOLDER_ZH), (EN, PLACEHOLDER_EN)):
    t = open(path, encoding="utf-8", newline="").read()
    if ph not in t:
        print(f"  !! 占位符未找到: {path}")
        continue
    t = t.replace(ph, ref_block)
    open(path, "w", encoding="utf-8", newline="").write(t)
    print(f"  References 节已插入: {path.split(chr(92))[-1]}")

# ---------- 5. 一致性终检 ----------
print("\n终检：")
for path in (ZH, EN):
    t = open(path, encoding="utf-8", newline="").read()
    body = t.split("# References")[0]
    marks = [m for m in re.findall(r"\[[\d,\-\s]+\]", body) if m.strip() != "[0,1]"]
    nums = sorted({int(x) for m in marks for x in re.findall(r"\d+", m)})
    expand = set()
    for m in marks:
        for part in m[1:-1].split(","):
            part = part.strip()
            if "-" in part:
                a, b = part.split("-")
                expand |= set(range(int(a), int(b) + 1))
            elif part:
                expand.add(int(part))
    leftover = re.findall(r"\[Wang[^\]]*\]|\[Galluzzi[^\]]*\]|\[PMID \d+\]|PMID 36747806", t)
    n_entries = len(re.findall(r"^\[\d+\] ", t, re.M))
    print(f"  {path.split(chr(92))[-1]}: 标记 {len(marks)} 处 | 展开编号 {min(expand)}–{max(expand)} 共 {len(expand)} 个 | "
          f"覆盖 1..{len(order)}: {expand == set(range(1, len(order)+1))} | 条目行 {n_entries}")
zh_marks = re.findall(r"\[[\d,\-\s]+\]", open(ZH, encoding='utf-8').read().split("# References")[0])
en_marks = re.findall(r"\[[\d,\-\s]+\]", open(EN, encoding='utf-8').read().split("# References")[0])
print("  双版标记序列一致:", [m for m in zh_marks if m.strip() != '[0,1]'] == [m for m in en_marks if m.strip() != '[0,1]'])

# 保存映射留痕
with open(BASE + r"\_intermediate\P0_refs_gbt_mapping.json", "w", encoding="utf-8") as f:
    json.dump({"order": order, "mapping": mapping, "entries": dict(zip([str(i) for i in range(1, len(order)+1)], lines_out))},
              f, ensure_ascii=False, indent=1)
print("saved _intermediate/P0_refs_gbt_mapping.json")
