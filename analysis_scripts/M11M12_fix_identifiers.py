# -*- coding: utf-8 -*-
"""M11M12_fix_identifiers.py — 修复标识列与共享库 BH bug（一次性补丁，留痕）"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import pandas as pd
import numpy as np

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
INTER = ROOT + r"\_intermediate"

# 1) 修复 per_sample 中 GSE215865 的 sample_id（元组 repr -> subj|day）
S = pd.read_csv(INTER + r"\M11M12_step1_per_sample.csv", low_memory=False)
mask = S["cohort"] == "GSE215865"
S.loc[mask, "sample_id"] = S.loc[mask, "subject_norm"].astype(str) + "|" + S.loc[mask, "day_label"].astype(str)
S.to_csv(INTER + r"\M11M12_step1_per_sample.csv", index=False)
print("GSE215865 sample_id 修复:", S.loc[mask, "sample_id"].head(3).tolist())

# 2) 修复 mdi_lib.bh（顺序性 bug：最小 p 不在首位时 q 错位）
lib_path = ROOT + r"\mdi_lib.py"
src = open(lib_path, encoding="utf-8").read()
old_bh = '''def bh(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum(1, p[order] * n / (np.arange(n) + 1))
    q = np.minimum.accumulate(q[np.argsort(order)][::-1])[::-1]
    return q'''
new_bh = '''def bh(pvals):
    """BH 校正（2026-08-27 修复：原实现当最小 p 不在首位时 q 值错位，见 03_LOGS/M11M12_step4_log 审计注）。"""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q_sorted = np.minimum(1, p[order] * n / (np.arange(n) + 1))
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q = np.empty(n)
    q[order] = q_sorted
    return q'''
assert old_bh in src, "bh 源文本未匹配"
src = src.replace(old_bh, new_bh)
open(lib_path, "w", encoding="utf-8").write(src)
import importlib, mdi_lib
importlib.reload(mdi_lib)
print("bh 修复验证:", mdi_lib.bh(np.array([0.1231, 0.0003249])))  # 应为 [0.1231, 0.0006498]
print("bh 修复验证2:", np.round(mdi_lib.bh(np.array([0.0506, 0.0085, 0.539, 0.0075, 0.189, 0.113])), 4))
