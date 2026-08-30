# -*- coding: utf-8 -*-
"""
P0_refs_vancouver_fix.py — Vancouver 重排修正版（2026-08-29）
修正：行尾统一 CRLF（原文件主体 CRLF，前次插入块 LF 导致混合+正则失配）；
  [33] 无页码电子刊用 年;卷(期). doi 形式；[48] elocationid pii+doi 提取真 DOI；
  [49] 句首连字符词 Large-Scale→Large-scale；[39]/[40] 保护 UK Biobank / GTEx Consortium；
  [25] 去 gds 题名 " [Array]" 后缀；[47] 题名以 ?/! 结尾不补句点。
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

PROTECT_PHRASES = ["Nomenclature Committee on Cell Death", "UK Biobank", "GTEx Consortium"]
PROTECT_WORDS = {"Ca2+", "Berlin", "Perturb-seq"}

def lower_word(w):
    if w in PROTECT_WORDS or not re.match(r"^[A-Z]", w):
        return w
    if re.search(r"[A-Z]", w[1:]) or re.match(r"^[A-Z]\d", w) or re.match(r"^[A-Z]\.$", w):
        return w
    if "-" in w:
        segs = w.split("-")
        if all(re.match(r"^[A-Za-z0-9.]+$", s) for s in segs):
            return "-".join((s[0].lower() + s[1:]) if re.match(r"^[A-Z][a-z]+$", s) else s for s in segs)
        return w
    return w[0].lower() + w[1:]

def sentence_case(title):
    for i, ph in enumerate(PROTECT_PHRASES):
        title = title.replace(ph, f"\x00{i}\x00")
    words = title.split(" ")
    out = [words[0] if i == 0 else lower_word(w) for i, w in enumerate(words)]
    # 句首连字符词：首段保大写、后续段小写（Large-Scale→Large-scale）
    if "-" in out[0]:
        segs = out[0].split("-")
        if all(re.match(r"^[A-Za-z0-9.]+$", s) for s in segs):
            out[0] = "-".join(s if (k == 0 or not re.match(r"^[A-Z][a-z]+$", s)) else s[0].lower() + s[1:]
                               for k, s in enumerate(segs))
    s = " ".join(out)
    for i, ph in enumerate(PROTECT_PHRASES):
        s = s.replace(f"\x00{i}\x00", ph)
    return s

def fmt_authors(authors):
    names = [a for a in authors if a.strip()]
    return (", ".join(names[:6]) + ", et al") if len(names) > 6 else ", ".join(names)

def doi_of(m):
    eloc = (m.get("doi") or "").strip()
    mm = re.search(r"10\.\d{4,}/[^\s]+", eloc)          # elocationid 可能含 "pii: ... doi: 10.1101/..."
    d = mm.group(0).rstrip(".") if mm else ""
    return d or None

def fmt_journal(m):
    au = fmt_authors(m["authors"])
    ti = sentence_case(m["title"]).rstrip(".")
    head = f"{au}. {ti}. {m['journal']}. {m['year']}"
    if m["volume"]:
        vi = f"{m['volume']}" + (f"({m['issue']})" if m["issue"] else "")
        head += f";{vi}" + (f":{m['pages']}" if m["pages"] else "")
    d = doi_of(m)
    return head + ("?" if ti.endswith("?") else ".") + (f" doi:{d}" if d else "")

def fmt_geo(n, acc):
    g = META["geo"][acc]
    title = re.sub(r"\s*\[(Array|Series|GDS|Platform)\]$", "", g["title"].rstrip("."))
    return (f"NCBI Gene Expression Omnibus. {acc}: {sentence_case(title)} [Internet]. "
            f"Bethesda, MD: National Center for Biotechnology Information; {g['pdat'][:4]} "
            f"[cited 2026 Aug 29]. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}")

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
    return fmt_journal(META["pubmed"][key])

entries = [f"{i}. {entry_for(i, slot)}" for i, slot in enumerate(order, 1)]
print("=== 修正后关键条目复核 ===")
for i in (2, 6, 17, 20, 25, 33, 39, 40, 44, 47, 48, 49):
    print("  " + entries[i - 1][:165])

block = "\n".join(entries)
pat = re.compile(r"(# References\r?\n\r?\n).*?(\r?\n\r?\n## Data availability)", re.S)
for path in (ZH, EN):
    raw = open(path, encoding="utf-8", newline="").read()
    raw = raw.replace("\r\n", "\n")            # 统一：先归一 LF
    if not pat.search(raw):
        print(f"!! 定位失败: {path}"); continue
    raw = pat.sub(lambda m: m.group(1) + block + m.group(2), raw, count=1)
    raw = raw.replace("\n", "\r\n")            # 再统一回 CRLF（原文件主体行尾）
    open(path, "w", encoding="utf-8", newline="").write(raw)
    n = len(re.findall(r"^\d+\. ", raw.split("# References")[1], re.M))
    print(f"已替换 {path.split(chr(92))[-1]}：Vancouver 条目 {n} 条")

json.dump({"format": "Vancouver/NLM", "entries": entries},
          open(BASE + r"\_intermediate\P0_refs_vancouver_entries.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)

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
    crlf, lf = t.count("\r\n"), t.count("\n") - t.count("\r\n")
    print(f"  {path.split(chr(92))[-1]}: 标记 {len(marks)} | 覆盖 1..{len(order)}: {expand == set(range(1, len(order)+1))} | 行尾 CRLF={crlf} 纯LF={lf}")
