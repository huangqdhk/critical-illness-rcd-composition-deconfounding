# -*- coding: utf-8 -*-
"""2026-09-05 CDD 投稿适配：主图/附图面板重排（两阶段安全重命名）。

主图：11 图 64 面板（2 图 16 面板超载 + 9/10/11 图各仅 3 面板 + 引用乱序：§2 引图3、§3 引图2I–2P、§4 引图2A–2H）
  → 9 图 64 面板（6-6-8-8-10-6-7-7-6），全部按正文引用顺序编号，每图 6–10 面板。
附图：12 图 84 面板（S3 仅 4 面板、S12 达 10 且三主题混杂、S12 先于 S9/S10 被引）
  → 12 图 84 面板（5-10-8-9-6-8-7-5-9-7-5-5），S2+S3 合并（同属无监督结构发现），
    旧 S12 拆分为因果层（新 S8，§8 首引）与风险预测/深度学习/表观时钟（新 S12，§16 首引），
    全文引用严格升序。

数值零改动：仅文件名迁移；sha256 不变；RESULTS_MANIFEST/lint 同步更新（任务8 图号迁移同口径）。
"""
import os, csv, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(ROOT, "01_FIGURE_DATA_CSV")

# 主图映射：旧面板 → 新面板（同名省略）
MAIN_MAP = {
    # 旧 Figure 3（模块级解离）→ 新 Figure 2（§2 首引）
    "3A": "2A", "3B": "2B", "3C": "2C", "3D": "2D", "3E": "2E", "3F": "2F",
    # 旧 Figure 2I–2P（跨平台塌陷+GSE32707）→ 新 Figure 3（§3）
    "2I": "3A", "2J": "3B", "2K": "3C", "2L": "3D", "2M": "3E", "2N": "3F", "2O": "3G", "2P": "3H",
    # 旧 Figure 2A–2H（五维收敛）→ 新 Figure 4（§4）
    "2A": "4A", "2B": "4B", "2C": "4C", "2D": "4D", "2E": "4E", "2F": "4F", "2G": "4G", "2H": "4H",
    # 旧 Figure 4（MDI 临床化）→ 新 Figure 5A–5E（§6）
    "4A": "5A", "4B": "5B", "4C": "5C", "4D": "5D", "4E": "5E",
    # 旧 Figure 5（生化层）→ 新 Figure 5F–5J（§7，与 MDI 合并为 Part II 边界图）
    "5A": "5F", "5B": "5G", "5C": "5H", "5D": "5I", "5E": "5J",
    # 旧 Figure 9（干预响应）→ 新 Figure 8E–8G（§13，与时序合并）
    "9A": "8E", "9B": "8F", "9C": "8G",
    # 旧 Figure 10（开放工具）→ 新 Figure 9A–9C（§14）
    "10A": "9A", "10B": "9B", "10C": "9C",
    # 旧 Figure 11（力学边界）→ 新 Figure 9D–9F（§15，与工具合并为 Part V 推广图）
    "11A": "9D", "11B": "9E", "11C": "9F",
    # Figure 1A–F / 6A–F / 7A–G / 8A–D 不动
}

# 附图映射：旧面板（去 S 前缀）→ 新面板
SUPP_MAP = {
    # 旧 S3（EFA，4 面板）并入新 S2（WGCNA+EFA 无监督结构发现）
    "3A": "2G", "3B": "2H", "3C": "2I", "3D": "2J",
    # 旧 S4–S8 顺移为 S3–S7
    "4A": "3A", "4B": "3B", "4C": "3C", "4D": "3D", "4E": "3E", "4F": "3F", "4G": "3G", "4H": "3H",
    "5A": "4A", "5B": "4B", "5C": "4C", "5D": "4D", "5E": "4E", "5F": "4F", "5G": "4G", "5H": "4H", "5I": "4I",
    "6A": "5A", "6B": "5B", "6C": "5C", "6D": "5D", "6E": "5E", "6F": "5F",
    "7A": "6A", "7B": "6B", "7C": "6C", "7D": "6D", "7E": "6E", "7F": "6F", "7G": "6G", "7H": "6H",
    "8A": "7A", "8B": "7B", "8C": "7C", "8D": "7D", "8E": "7E", "8F": "7F", "8G": "7G",
    # 旧 S12 因果层面板 → 新 S8（§8 因果边界首引）
    "12C": "8A", "12D": "8B", "12E": "8C", "12I": "8D", "12J": "8E",
    # 旧 S12 风险预测/深度学习/表观时钟 → 新 S12（保留 S12 号，§16 引用）
    "12F": "12C", "12G": "12D", "12H": "12E",
    # S1A–E / S2A–F / S9A–I / S10A–G / S11A–E / S12A–B 不动
}

def panel_sort_key(p):
    import re
    m = re.match(r"(\d+)([A-Z]+)", p)
    return (int(m.group(1)), m.group(2))

def apply(folder, prefix, mapping):
    moved, skipped = [], []
    # 两阶段：先全部改为临时名，再落终名（避免新旧名碰撞）
    tmp = {}
    for old, new in sorted(mapping.items(), key=lambda kv: panel_sort_key(kv[0])):
        if old == new:
            continue
        src = os.path.join(BASE, folder, f"{prefix}{old}.csv")
        t = os.path.join(BASE, folder, f"__tmp_{prefix}{old}.csv")
        if not os.path.exists(src):
            print(f"[WARN] 源缺失: {src}")
            continue
        os.rename(src, t)
        tmp[old] = (t, new)
    for old, (t, new) in tmp.items():
        dst = os.path.join(BASE, folder, f"{prefix}{new}.csv")
        if os.path.exists(dst):
            print(f"[FATAL] 目标已存在（不应发生）: {dst}")
            sys.exit(1)
        os.rename(t, dst)
        moved.append((f"{prefix}{old}", f"{prefix}{new}"))
    return moved

def main():
    m1 = apply("Main", "Figure_", MAIN_MAP)
    m2 = apply("Supplementary", "Figure_S", SUPP_MAP)
    # 映射台账
    led = os.path.join(BASE, "FIGURE_RENUMBER_MAP_20260905.csv")
    with open(led, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["old_panel_file", "new_panel_file", "package", "change_date", "reason"])
        for old, new in m1:
            w.writerow([old + ".csv", new + ".csv", "Main", "2026-09-05",
                        "CDD 投稿适配：主图 11→9、面板均衡 6–10/图、按正文引用顺序重编号（数值零改动）"])
        for old, new in m2:
            w.writerow([old + ".csv", new + ".csv", "Supplementary", "2026-09-05",
                        "CDD 投稿适配：附图 S2+S3 合并、旧 S12 拆分（因果层→S8；ML/DL/表观→S12），引用升序（数值零改动）"])
    print(f"renamed: Main {len(m1)} + Supplementary {len(m2)} = {len(m1)+len(m2)}")
    print("ledger:", led)

if __name__ == "__main__":
    main()
