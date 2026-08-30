# -*- coding: utf-8 -*-
"""
P0_refs_gbt_fixup.py — 插入修复（2026-08-29）：
  1. EN 版残留的作者-年份括号（分号变体）
  2. 修正条目键名（Galluzzi→A:NCCD 等）与名缩写加空格（LAFFEY JG→LAFFEY J G）
  3. 以模糊正则替换被数字重排波及的占位行，插入 58 条 GB/T 7714-2025 列表
  4. 重跑终检
注意：不重跑数字重排（映射不可幂等）。
"""
import json, re, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
ZH, EN = BASE + r"\111文稿_v5.md", BASE + r"\111文稿_v5_英文版.md"
META = json.load(open(BASE + r"\_intermediate\P0_refs_gbt_metadata.json", encoding="utf-8"))
MAP = json.load(open(BASE + r"\_intermediate\P0_refs_gbt_mapping.json", encoding="utf-8"))
order = MAP["order"]; mapping = MAP["mapping"]
CITE_DATE = "2026-08-29"

# ---- 1. EN 残留作者-年份括号 ----
t = open(EN, encoding="utf-8", newline="").read()
old_en = "[Wang 2025; Casadio *Nat Cell Biol* 2026, PMID 41540269]"
if old_en in t:
    t = t.replace(old_en, "[%d,%d]" % (mapping["Wang"], mapping["Casadio"]))
    open(EN, "w", encoding="utf-8", newline="").write(t)
    print("EN 分号变体已替换")
else:
    print("EN 分号变体不存在（可能已处理）:", re.findall(r"\[Wang[^\]]*\]", t))

# ---- 2. 条目构建 ----
SLOT2KEY = {"Galluzzi": "A:NCCD", "Wang": "A:Wang", "Casadio": "A:Casadio",
            "Kapp": "A:Kapp", "Yao": "N:Yao", "Overmyer": "N:Overmyer"}
INST_HINTS = ("Consortium", "Task Force", "Network", "Group", "eQTLGen", "GTEx", "UniProt", "Collaboration")

def fmt_authors(authors):
    out = []
    for a in authors:
        if any(k in a for k in INST_HINTS) or " " not in a.strip():
            out.append(a); continue
        parts = a.split()
        ini = parts[-1]
        ini = " ".join(list(ini)) if len(ini) >= 2 and ini.isalpha() and ini.isupper() else ini
        out.append(" ".join(parts[:-1]).upper() + " " + ini)
    return (", ".join(out[:3]) + ", et al") if len(out) > 3 else ", ".join(out)

def fmt_journal(m):
    vi = f"{m['volume']}({m['issue']})" if m["issue"] else (m["volume"] or m["year"])
    return f"{fmt_authors(m['authors'])}. {m['title']}[J]. {m['journal']}, {m['year']}, {vi}: {m['pages']}."

GEO_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
SLOT2ACC = {"18": "GSE165659", "27": "GSE212865", "65": "GSE271370", "67": "GSE310929", "69": "GSE148871"}
def fmt_geo(acc):
    g = META["geo"][acc]
    return (f"NCBI Gene Expression Omnibus. {acc}: {g['title'].rstrip('.')}[DS/OL]. "
            f"({g['pdat']})[{CITE_DATE}]. {GEO_URL}{acc}.")

def entry_for(slot):
    if slot in SLOT2ACC:
        return fmt_geo(SLOT2ACC[slot])
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
        return f"{fmt_authors(m['authors'])}. {m['title']}[PP/OL]. bioRxiv({m['year']})[{CITE_DATE}]. {url}."
    key = SLOT2KEY.get(slot, str(slot))
    key = key if key in META["pubmed"] else str(slot)
    if key not in META["pubmed"]:
        return None
    return fmt_journal(META["pubmed"][key])

entries = []
for idx, slot in enumerate(order, 1):
    e = entry_for(slot)
    if e is None:
        print(f"!! 缺条目: {slot}"); continue
    entries.append(f"[{idx}] {e}")
print(f"条目构建完成：{len(entries)}/{len(order)}")
ref_block = "\n".join(entries)

# ---- 3. 替换占位行（模糊匹配：数字重排已把 [1]–[77] 改写） ----
PH_PAT_ZH = re.compile(r"【占位——参考文献全文见[^\n]*?】")
PH_PAT_EN = re.compile(r"【Placeholder[^\n]*?】")
for path, pat in ((ZH, PH_PAT_ZH), (EN, PH_PAT_EN)):
    t = open(path, encoding="utf-8", newline="").read()
    m = pat.search(t)
    if not m:
        print(f"!! 占位行未找到: {path}"); continue
    t = t[:m.start()] + ref_block + t[m.end():]
    open(path, "w", encoding="utf-8", newline="").write(t)
    print(f"References 节已插入: {path.split(chr(92))[-1]}（替换了 {len(m.group(0))} 字符占位行）")

# ---- 4. 终检 ----
print("\n终检：")
all_marks = {}
for path in (ZH, EN):
    t = open(path, encoding="utf-8", newline="").read()
    body = t.split("# References")[0]
    marks = [m for m in re.findall(r"\[[\d,\-\s]+\]", body) if m.strip() != "[0,1]"]
    expand = set()
    for m in marks:
        for part in m[1:-1].split(","):
            part = part.strip()
            if "-" in part:
                a, b = part.split("-"); expand |= set(range(int(a), int(b) + 1))
            elif part:
                expand.add(int(part))
    n_entries = len(re.findall(r"^\[\d+\] ", t, re.M))
    leftover = re.findall(r"\[Wang[^\]]*\]|\[Galluzzi[^\]]*\]|\[PMID \d+\]|PMID 36747806|Overmyer et al\. \*Cell", t)
    all_marks[path] = [m for m in marks]
    ok = expand == set(range(1, len(order) + 1))
    print(f"  {path.split(chr(92))[-1]}: 标记 {len(marks)} 处 | 编号覆盖 1..{len(order)}: {ok} | 条目行 {n_entries} | 残留作者-年份: {leftover[:3]}")
print("  双版标记序列一致:", all_marks[ZH] == all_marks[EN])

json.dump({"entries": dict(zip([str(i) for i in range(1, len(entries) + 1)], entries))},
          open(BASE + r"\_intermediate\P0_refs_gbt_final_entries.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("saved _intermediate/P0_refs_gbt_final_entries.json")
print("\n样例条目：")
for i in (6, 20, 27, 37, 46, 47, 48, 49):
    print("  " + entries[i - 1][:170])
