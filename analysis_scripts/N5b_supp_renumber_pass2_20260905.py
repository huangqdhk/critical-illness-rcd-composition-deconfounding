# -*- coding: utf-8 -*-
"""2026-09-05（下午 QC 批次）附图二次重排：按正文首次引用顺序严格升序重新编号 S1–S12。

背景：上午首排后补入 S1–S7 正文首引，附图首引顺序出现错位（S7@§3、S4@§8、S3@§16 等），
不符合 Nature Portfolio（含 CDD）"附图按引用顺序编号"的投稿要求。本脚本执行：
  1) Supplementary/ CSV 按 {旧号→新号} 字母保持式两阶段重命名；
  2) FIGURE_RENUMBER_MAP 台账追加二次迁移行；
  3) RESULTS_MANIFEST_v2.0.csv 文件名同步（sha256 不变、数值零改动）；
  4) lint_package.py 白名单键 Figure_S9I.csv → Figure_S6I.csv；
  5) 文稿内 "Figure(s) S…" 引用令牌整体重映射（Table S 引用不动；独立散文段与 causal note 由人工随后重写）。
映射（按首引位置）：S1→S1, S7→S2, S2→S3, S8→S4, S4→S5, S9→S6, S10→S7, S5→S8, S6→S9, S11→S10, S12→S11, S3→S12。
"""
import os, re, csv, json

ROOT = os.path.dirname(os.path.abspath(__file__))
SUP = os.path.join(ROOT, "01_FIGURE_DATA_CSV", "Supplementary")
MAP = {"S1": "S1", "S7": "S2", "S2": "S3", "S8": "S4", "S4": "S5", "S9": "S6",
       "S10": "S7", "S5": "S8", "S6": "S9", "S11": "S10", "S12": "S11", "S3": "S12"}

def rename_supplementary():
    files = [f for f in os.listdir(SUP) if re.match(r"Figure_S\d+[A-Z]+\.csv$", f)]
    tmp = []
    for f in files:
        m = re.match(r"Figure_(S\d+)([A-Z]+)\.csv$", f)
        new = f"Figure_{MAP[m.group(1)]}{m.group(2)}.csv"
        if new == f:
            continue
        os.rename(os.path.join(SUP, f), os.path.join(SUP, "__tmp2_" + f))
        tmp.append((f, new))
    for old, new in tmp:
        dst = os.path.join(SUP, new)
        assert not os.path.exists(dst), f"collision {new}"
        os.rename(os.path.join(SUP, "__tmp2_" + old), dst)
    return tmp

def append_ledger(moves):
    p = os.path.join(ROOT, "01_FIGURE_DATA_CSV", "FIGURE_RENUMBER_MAP_20260905.csv")
    with open(p, "a", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        for old, new in moves:
            w.writerow([old, new, "Supplementary", "2026-09-05(二次)",
                        "CDD 适配二次重排：附图按正文首次引用顺序严格升序重编号（S7→S2/S2→S3/S8→S4/S4→S5/S9→S6/S10→S7/S5→S8/S6→S9/S11→S10/S12→S11/S3→S12；数值零改动）"])

def update_manifest(moves):
    p = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "RESULTS_MANIFEST_v2.0.csv")
    rows = list(csv.reader(open(p, encoding="utf-8-sig")))
    hdr = rows[0]
    i_fn, i_nt = hdr.index("filename"), hdr.index("notes")
    mv = {old: new for old, new in moves}
    n = 0
    for r in rows[1:]:
        if len(r) > i_fn and r[i_fn] in mv:
            old = r[i_fn]
            r[i_fn] = mv[old]
            r[i_nt] = (r[i_nt] + "；" if r[i_nt] else "") + \
                f"2026-09-05 二次重排：原 {old}（首引顺序升序化，数值零改动）"
            n += 1
    with open(p, "w", newline="", encoding="utf-8-sig") as f:
        csv.writer(f).writerows(rows)
    return n

def update_lint():
    p = os.path.join(ROOT, "lint_package.py")
    s = open(p, encoding="utf-8").read()
    assert '"Figure_S9I.csv", "FDR"' in s
    s = s.replace('("Figure_S9I.csv", "FDR")', '("Figure_S6I.csv", "FDR")')
    s = s.replace("（原 Figure_5I，任务8 降入补充 S9）",
                  "（原 Figure_5I，任务8 降入补充；2026-09-05 二次重排 S9I→S6I）")
    open(p, "w", encoding="utf-8").write(s)

LIST_PAT = re.compile(
    r"(Figures? )((?:S?\d+[A-Z]*)(?:(?:\s*[–/,]\s*|\s+and\s+)(?:S?\d+[A-Z]*)+)*)")

def remap_tokens(m):
    head, lst = m.group(1), m.group(2)
    def tok(t):
        mm = re.match(r"S(\d+)([A-Z]*)$", t)
        if mm and ("S" + mm.group(1)) in MAP:
            return MAP["S" + mm.group(1)] + mm.group(2)
        return t
    parts = re.split(r"(\s*[–/,]\s*|\s+and\s+)", lst)
    out = "".join(tok(p) if not re.match(r"\s*[–/,]\s*$", p) and p.strip() and not re.match(r"\s+and\s+", p) else p for p in parts)
    return head + out

def remap_manuscript():
    p = os.path.join(ROOT, "111文稿_v5_英文版.md")
    s = open(p, encoding="utf-8").read()
    s2 = LIST_PAT.sub(remap_tokens, s)
    open(p, "w", encoding="utf-8").write(s2)

def main():
    moves = rename_supplementary()
    print("renamed:", len(moves))
    append_ledger(moves)
    print("manifest rows updated:", update_manifest(moves))
    update_lint()
    remap_manuscript()
    print("done; 散文段与 causal note 待人工重写")

if __name__ == "__main__":
    main()
