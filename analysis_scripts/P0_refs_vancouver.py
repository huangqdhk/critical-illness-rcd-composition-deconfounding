# -*- coding: utf-8 -*-
"""
P0_refs_vancouver.py — Vancouver/NLM 格式参考文献重排（2026-08-29）
编号不变（首次出现序=Vancouver 标准序），仅重排条目著录：
  N. 作者(≤6 全列，>6 前6+et al，姓+名缩写不加点). 句首式题名. 刊名缩写. 年;卷(期):页. doi:xxx
  网络/数据资源：NLM [Internet] 式；预印本：bioRxiv 年 + doi
"""
import json, re, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
ZH, EN = BASE + r"\111文稿_v5.md", BASE + r"\111文稿_v5_英文版.md"
META = json.load(open(BASE + r"\_intermediate\P0_refs_gbt_metadata.json", encoding="utf-8"))
MAP = json.load(open(BASE + r"\_intermediate\P0_refs_gbt_mapping.json", encoding="utf-8"))
order = MAP["order"]

SLOT2KEY = {"Galluzzi": "A:NCCD", "Wang": "A:Wang", "Casadio": "A:Casadio",
            "Kapp": "A:Kapp", "Yao": "N:Yao", "Overmyer": "N:Overmyer"}
SLOT2ACC = {"18": "GSE165659", "27": "GSE212865", "65": "GSE271370", "67": "GSE310929", "69": "GSE148871"}

# ---------- 句首大写转换 ----------
PROTECT_PHRASES = ["Nomenclature Committee on Cell Death"]
PROTECT_WORDS = {"Ca2+", "Berlin", "Meningitidis", "Perturb-seq"}

def lower_word(w):
    if w in PROTECT_WORDS:
        return w
    if not re.match(r"^[A-Z]", w):          # 非首字母大写：不动
        return w
    if re.search(r"[A-Z]", w[1:]):           # 词内还有大写（缩略词/混合）：不动
        return w
    if re.match(r"^[A-Z]\d", w):             # L1000 / IFI27 型
        return w
    if re.match(r"^[A-Z]\.$", w):            # 属名缩写 N.
        return w
    if "-" in w:                              # Large-Scale → large-scale（含内部大写的 SARS-CoV-2 已在上一步保留）
        segs = w.split("-")
        if all(re.match(r"^[A-Za-z0-9.]+$", s) for s in segs):
            fixed = [ (s[0].lower() + s[1:]) if re.match(r"^[A-Z][a-z]+$", s) else s for s in segs ]
            return "-".join(fixed)
        return w
    return w[0].lower() + w[1:]

def sentence_case(title):
    for i, ph in enumerate(PROTECT_PHRASES):          # 短语保护
        title = title.replace(ph, f"\x00{i}\x00")
    words = title.split(" ")
    out = []
    for i, w in enumerate(words):
        out.append(w if i == 0 else lower_word(w))    # 仅整题首词保大写
    s = " ".join(out)
    for i, ph in enumerate(PROTECT_PHRASES):
        s = s.replace(f"\x00{i}\x00", ph)
    return s

# ---------- 著录 ----------
INST = ("Consortium", "Task Force", "Network", "Group", "eQTLGen", "UniProt", "Collaboration")

def fmt_authors(authors):
    names = [a for a in authors if a.strip()]
    if len(names) > 6:
        return ", ".join(names[:6]) + ", et al"
    return ", ".join(names)

def doi_of(m):
    d = (m.get("doi") or "").strip()
    d = re.sub(r"^(doi|DOI)\s*:\s*", "", d).strip()
    return d or None

def fmt_journal(n, m):
    au = fmt_authors(m["authors"])
    vi = m["volume"] and f"{m['volume']}({m['issue']})" if m["issue"] else (m["volume"] or "")
    core = f"{au}. {sentence_case(m['title'])}. {m['journal']}. {m['year']}"
    if vi:
        core += f";{vi}:{m['pages']}"
    d = doi_of(m)
    return core + ("." if not d else f". doi:{d}")

def fmt_geo(n, acc):
    g = META["geo"][acc]
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}"
    return (f"NCBI Gene Expression Omnibus. {acc}: {sentence_case(g['title'].rstrip('.'))} [Internet]. "
            f"Bethesda, MD: National Center for Biotechnology Information; {g['pdat'][:4]} "
            f"[cited 2026 Aug 29]. Available from: {url}")

def entry_for(i, slot):
    if slot in SLOT2ACC:
        return fmt_geo(i, SLOT2ACC[slot])
    if slot == "74":
        return ("MetaboLights. MTBLS6844: serum lipidome of COVID-19 patients (FIA-MS positive ion mode, "
                "53 samples) [Internet]. Hinxton, UK: European Molecular Biology Laboratory, European "
                "Bioinformatics Institute [cited 2026 Aug 29]. Available from: https://www.ebi.ac.uk/metabolights/MTBLS6844")
    if slot == "77":
        return ("UniProt Consortium. UniProt: the Universal Protein Knowledgebase [Internet]. "
                "[cited 2026 Aug 29]. Available from: https://www.uniprot.org")
    if slot == "Yao":
        m = META["pubmed"]["N:Yao"]
        d = doi_of(m)
        tail = f"doi:{d}" if d else f"https://pubmed.ncbi.nlm.nih.gov/{m['pmid']}/"
        return f"{fmt_authors(m['authors'])}. {sentence_case(m['title'])}. bioRxiv. {m['year']}. {tail}"
    key = SLOT2KEY.get(slot, str(slot))
    key = key if key in META["pubmed"] else str(slot)
    return fmt_journal(i, META["pubmed"][key])

entries = []
for i, slot in enumerate(order, 1):
    entries.append(f"{i}. {entry_for(i, slot)}")

print("=== Vancouver 条目全文（供人工核对） ===")
for e in entries:
    print(e)

# ---------- 替换双版 References 块 ----------
block = "\n".join(entries)
pat = re.compile(r"(# References\n\n).*?(\n\n## Data availability)", re.S)
for path in (ZH, EN):
    t = open(path, encoding="utf-8", newline="").read()
    if not pat.search(t):
        print(f"!! References 块定位失败: {path}"); continue
    t2 = pat.sub(lambda m: m.group(1) + block + m.group(2), t, count=1)
    open(path, "w", encoding="utf-8", newline="").write(t2)
    n = len(re.findall(r"^\d+\. ", t2.split("# References")[1], re.M))
    print(f"\n已替换 {path.split(chr(92))[-1]}：Vancouver 条目 {n} 条")

json.dump({"format": "Vancouver/NLM", "entries": entries},
          open(BASE + r"\_intermediate\P0_refs_vancouver_entries.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)

# ---------- 终检 ----------
print("\n终检：")
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
    print(f"  {path.split(chr(92))[-1]}: 正文标记 {len(marks)} | 覆盖 1..{len(order)}: {expand == set(range(1, len(order)+1))}")
