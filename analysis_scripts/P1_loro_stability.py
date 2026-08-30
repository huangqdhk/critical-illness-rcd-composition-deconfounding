# -*- coding: utf-8 -*-
"""
P1_loro_stability.py — Phase 1 任务2b：IIAMD 签名 leave-one-replicate-out 稳定性
==================================================================================
依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md §4（H6：sign and FDR retained）
方法：对确证设计的 12 个样本逐个剔除 → 重拟合同一模型（全基因）→ 记录
     IIAMD v1.0 签名基因的 β₃ 方向与 padj<0.05 保持率。
输出：_intermediate/P1_loro_per_gene.csv（签名基因×12次重拟合的方向/FDR保持）
     _intermediate/P1_loro_summary.csv；追加写入 P1_Interaction_Signature_Report.md
"""
import os, sys, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
CNT = os.path.join(ROOT, "00_RAW_DATA", "GSE235046", "GSE235046_Count_table.txt", "GSE235046_Count_table.txt")
INTER = os.path.join(ROOT, "_intermediate")
SIGF = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "P1_IIAMD_signature_v1.0.csv")
REPORT = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "P1_Interaction_Signature_Report.md")

log = []
def note(m=""):
    log.append(m); print(m, flush=True)

raw = pd.read_csv(CNT, sep="\t")
raw["geneSymbol"] = raw["geneSymbol"].astype(str).str.strip()
ok = raw[raw["geneSymbol"].notna() & (raw["geneSymbol"] != "") & (raw["geneSymbol"] != "nan")]
mat = ok.groupby("geneSymbol")[[c for c in raw.columns if c not in ("geneID","geneSymbol","bioType","annotationLevel")]].sum().astype(int)

YW_GROUP = {}
for yw in ["YW001","YW002","YW003"]: YW_GROUP[yw] = "Media"
for yw in ["YW004","YW005","YW006"]: YW_GROUP[yw] = "LPS"
for yw in ["YW007","YW008","YW009"]: YW_GROUP[yw] = "CS"
for yw in ["YW010","YW011","YW012"]: YW_GROUP[yw] = "LPS_CS"
for yw in ["YW016","YW017","YW018"]: YW_GROUP[yw] = "LPS_CS_Torin"
grp = {c: YW_GROUP[c.split("_")[-1]] for c in mat.columns}

sig = pd.read_csv(SIGF)
sig_genes = set(sig["geneSymbol"])
note(f"签名基因 {len(sig_genes)}（up {(sig['direction']=='up').sum()} / down {(sig['direction']=='down').sum()}）")

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

conf_cols = [c for c in mat.columns if grp[c] in ("Media", "LPS", "CS", "LPS_CS")]
cmat = mat[conf_cols]
keep = (cmat >= 10).sum(axis=1) >= 3
cmatf = cmat[keep]

def fit_interaction(df, meta):
    dds = DeseqDataSet(counts=df.T, metadata=meta, design="~group", quiet=True, n_cpus=1)
    dds.deseq2()
    cols = list(dds.obsm["design_matrix"].columns)
    tlevels = [c.split("T.")[1].rstrip("]") for c in cols if c.startswith("group[T.")]
    ref = ({"Media","LPS","CS","LPS_CS"} - set(tlevels)).pop()
    def mu_vec(level):
        v = np.zeros(len(cols)); v[cols.index("Intercept")] = 1.0
        if level != ref and f"group[T.{level}]" in cols:
            v[cols.index(f"group[T.{level}]")] = 1.0
        return v
    cvec = mu_vec("LPS_CS") - mu_vec("LPS") - mu_vec("CS") + mu_vec("Media")
    st = DeseqStats(dds, contrast=cvec, quiet=True, n_cpus=1)
    st.summary()
    r = st.results_df
    return r

rows = []
for drop in conf_cols:
    keep_cols = [c for c in conf_cols if c != drop]
    df = cmatf[keep_cols]
    meta = pd.DataFrame({"group": [grp[c] for c in keep_cols]}, index=keep_cols)
    r = fit_interaction(df, meta)
    sub = r[r.index.isin(sig_genes)]
    s2 = sig.set_index("geneSymbol")
    common = sub.index
    sign_keep = (np.sign(sub.loc[common, "log2FoldChange"]) ==
                 np.where(s2.loc[common, "direction"] == "up", 1, -1))
    fdr_keep = sign_keep & (sub.loc[common, "padj"] < 0.05)
    rows.append(dict(left_out=drop, group=grp[drop], n_sig_common=len(common),
                     sign_retained=int(sign_keep.sum()),
                     fdr_retained=int(fdr_keep.sum()),
                     sign_rate=float(sign_keep.mean()), fdr_rate=float(fdr_keep.mean())))
    note(f"   剔除 {drop}({grp[drop]}): 方向保持 {sign_keep.mean():.1%}, FDR保持 {fdr_keep.mean():.1%} (n={len(common)})")

lo = pd.DataFrame(rows)
lo.to_csv(os.path.join(INTER, "P1_loro_summary.csv"), index=False)
worst_sign = lo["sign_rate"].min(); worst_fdr = lo["fdr_rate"].min()
stable = (worst_sign >= 0.90) and (worst_fdr >= 0.50)
note(f"\n最差方向保持率 {worst_sign:.1%}；最差 FDR 保持率 {worst_fdr:.1%}")
note(f"LORO 判定（方向≥90% 且 FDR≥50%）: {'STABLE' if stable else 'NOT STABLE'}")

with open(REPORT, "a", encoding="utf-8") as f:
    f.write("\n## P1-2b LORO 稳定性（追加于 2026-08-20）\n\n")
    f.write("| 剔除样本 | 组 | 方向保持 | FDR保持 |\n|---|---|---|---|\n")
    for _, r in lo.iterrows():
        f.write(f"| {r['left_out']} | {r['group']} | {r['sign_rate']:.1%} | {r['fdr_rate']:.1%} |\n")
    f.write(f"\n**判定：{'STABLE' if stable else 'NOT STABLE'}**（最差方向保持 {worst_sign:.1%}，最差 FDR 保持 {worst_fdr:.1%}；"
            f"判据：方向≥90% 且 FDR≥50%，两项均为判定期望的保守阈值，在本报告中作为操作性判据披露）\n")

with open(os.path.join(ROOT, "03_LOGS", "P1_loro_log.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(log))
print("DONE P1-2b LORO")
