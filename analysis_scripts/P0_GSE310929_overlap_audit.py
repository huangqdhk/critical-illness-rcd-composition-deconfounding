# -*- coding: utf-8 -*-
"""
P0_GSE310929_overlap_audit.py — Phase 0 任务3：GSE310929 与其余 M2 队列样本重叠审计
====================================================================================
依据：《黄裕荣创新提质方案 v2（重塑版）》2026-08-19 §3.1「meta 重复纳入（致命）」
问题：M2_step2_replication_meta.py 把 GSE185263、GSE32707 与 GSE310929（多数据集再处理
      合并队列）同时纳入同一个随机效应 meta（k=3）。若 GSE310929 内含前两者样本，
      则同一生物学样本被重复计数，meta 独立性假设被破坏。
判定源（唯一）：GSE310929_AllSampleMetadataSubmitted.xlsx!all_samples_meta_cluster
      的 Dataset / GEO Accession / Patient ID 列；并与 M2_per_sample_scores.csv 的
      Dataset 列交叉一致性核对。
输出：04_AUDIT_GOVERNANCE/P0_GSE310929_Overlap_Audit_Report.md
      03_LOGS/P0_gse310929_overlap_log.txt
"""
import os, sys
import pandas as pd

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
XLSX = ROOT + r"\00_RAW_DATA\GSE310929_Sepsis\GSE310929_AllSampleMetadataSubmitted.xlsx"
M2CSV = ROOT + r"\_intermediate\M2_per_sample_scores.csv"
REPORT = ROOT + r"\04_AUDIT_GOVERNANCE\P0_GSE310929_Overlap_Audit_Report.md"
LOGF = ROOT + r"\03_LOGS\P0_gse310929_overlap_log.txt"

log = []
def note(m=""):
    log.append(m); print(m, flush=True)

# ---------------- 1. 读 GSE310929 逐样本元数据（header 在说明行下一行） ----------------
note("== 1. GSE310929 逐样本元数据（唯一判定源）")
meta = pd.read_excel(XLSX, sheet_name="all_samples_meta_cluster", header=1)
note(f"   shape={meta.shape}")
note(f"   关键列存在性: Dataset={('Dataset' in meta.columns)}, "
     f"GEO Accession={('GEO Accession' in meta.columns)}, Patient ID={('Patient ID' in meta.columns)}")

ds_counts = meta["Dataset"].value_counts()
note(f"\n== 2. 原始来源数据集全景（n={ds_counts.size} 个）")
for ds, n in ds_counts.items():
    note(f"   {ds}\t{n}")

# ---------------- 3. 与 M2 其余队列的来源级重叠 ----------------
M2_OTHERS = ["GSE185263", "GSE32707", "GSE212865", "GSE148871", "GSE188309"]
note("\n== 3. 来源级重叠（meta 双计数风险判定）")
overlap_rows = {}
for t in M2_OTHERS:
    sub = meta[meta["Dataset"].astype(str).str.strip().str.upper() == t]
    overlap_rows[t] = sub
    if len(sub):
        note(f"   [重叠] {t}: n={len(sub)} 例进入 GSE310929")
        dis = sub["Disease Simplified"].value_counts().to_dict()
        note(f"          Disease Simplified 构成: {dis}")
    else:
        note(f"   [无重叠] {t}: 0 例")

# ---------------- 4. GSM 层交叉（对用 GSM 号的队列做独立复核） ----------------
note("\n== 4. GSM 层独立复核（M2 分数表实际用到的样本号）")
m2 = pd.read_csv(M2CSV)
meta_gsm = set(meta["GEO Accession"].dropna().astype(str))
for t in ["GSE32707", "GSE212865", "GSE148871", "GSE188309"]:
    used = set(m2[m2["cohort"] == t]["sample"].dropna().astype(str))
    inter = used & meta_gsm
    note(f"   {t}: M2 用样 {len(used)}，与 GSE310929 的 GSM 交集 = {len(inter)}"
         + (f" （示例 {sorted(inter)[:5]}）" if inter else ""))

# GSE185263 用 sepcol* 患者号，无 GSM；用 Patient ID 列做模式复核
pat = meta["Patient ID"].dropna().astype(str)
sepcol_like = pat[pat.str.lower().str.contains("sepcol", na=False)]
note(f"   GSE185263: M2 用样 {len(set(m2[m2['cohort']=='GSE185263']['sample']))}（sepcol* 患者号），"
     f"GSE310929 Patient ID 含 sepcol 模式 = {len(sepcol_like)}")

# ---------------- 5. M2 分数表 Dataset 列交叉一致性核对 ----------------
note("\n== 5. M2_per_sample_scores.csv 自带 Dataset 列的一致性")
if "Dataset" in m2.columns:
    m2_310 = m2[m2["cohort"] == "GSE310929"]
    m2_ds = m2_310["Dataset"].dropna().astype(str).value_counts()
    note(f"   GSE310929 行数={len(m2_310)}，Dataset 列取值 n={m2_ds.size}")
    for ds, n in m2_ds.head(40).items():
        note(f"   {ds}\t{n}")
    # 两个来源的来源级重叠是否一致
    for t in M2_OTHERS:
        n_xlsx = len(overlap_rows[t])
        n_m2 = int(m2_ds.get(t, 0))
        flag = "一致" if n_xlsx == n_m2 else "不一致"
        note(f"   核对 {t}: xlsx={n_xlsx} vs M2分数表={n_m2} → {flag}")
else:
    note("   M2 分数表无 Dataset 列（跳过一致性核对）")

# ---------------- 6. 结论与影响面 ----------------
note("\n== 6. 审计结论（自动判定）")
meta_over = {t: len(overlap_rows[t]) for t in M2_OTHERS}
meta_over_nz = {t: n for t, n in meta_over.items() if n > 0}
if not meta_over_nz:
    note("   判定：GSE310929 与 M2 其余队列无样本重叠；「meta 重复纳入」指控不成立（虚惊）。")
    verdict = "NO_OVERLAP"
else:
    note(f"   判定：重叠成立 → {meta_over_nz}")
    note("   受污染分析：M2_meta.csv（k=3 随机效应 meta）与 GSE310929 的全部单队列对比；")
    note("   修复路径（P0-4）：以 Dataset 列拆回原始数据集 → 数据集级效应 + 来源互斥 + leave-one-source-out。")
    verdict = "OVERLAP_CONFIRMED"

# ---------------- 输出报告 ----------------
os.makedirs(os.path.dirname(REPORT), exist_ok=True)
with open(REPORT, "w", encoding="utf-8") as f:
    f.write("# P0-3 GSE310929 样本重叠审计报告\n\n")
    f.write(f"- 日期：2026-08-20\n")
    f.write("- 依据：《黄裕荣创新提质方案 v2（重塑版）》§3.1 meta 重复纳入（致命）\n")
    f.write("- 判定源：GSE310929_AllSampleMetadataSubmitted.xlsx（Dataset/GEO Accession/Patient ID）\n")
    f.write(f"- 脚本：P0_GSE310929_overlap_audit.py\n\n")
    f.write(f"## 判定：{verdict}\n\n")
    f.write("## 与 M2 其余队列的来源级重叠\n\n| 来源队列 | GSE310929 内例数 |\n|---|---|\n")
    for t in M2_OTHERS:
        f.write(f"| {t} | {meta_over[t]} |\n")
    f.write(f"\n## 原始来源数据集全景（{ds_counts.size} 个）\n\n| Dataset | n |\n|---|---|\n")
    for ds, n in ds_counts.items():
        f.write(f"| {ds} | {n} |\n")
    f.write("\n## 审计日志\n\n```\n" + "\n".join(log) + "\n```\n")

with open(LOGF, "w", encoding="utf-8") as f:
    f.write("\n".join(log))
print("\nDONE P0-3 audit →", REPORT)
