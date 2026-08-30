# -*- coding: utf-8 -*-
"""
P1_interaction_signature.py — Phase 1 任务2：GSE235046 2×2 交互模型 → IIAMD 签名（H6，Gate 1）
============================================================================================
依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md（OSF https://osf.io/C7RYD/）§2/§4/§5
设计（Wang et al. Cell 2025, BMDM WT, 5组×3重复）：
  确证模型（12 样本，Torin 臂不入确证）：
      Y ~ group ∈ {Media, LPS, CS, LPS_CS}，交互检验用对比向量
      β₃ = μ_LPS_CS − μ_LPS − μ_CS + μ_Media（DESeq2 Wald，pydeseq2 ~group + ndarray contrast）
  映射（结构对齐 + 运行号块）：YW001-003=Media / YW004-006=LPS / YW007-009=CS /
      YW010-012=LPS_CS / YW016-018=LPS_CS_Torin；LPS 标记基因方向校验必须 PASS 才继续
  IIAMD-up = {β₃>0, BH padj<0.05}；IIAMD-down = {β₃<0, padj<0.05}；签名冻结为 v1.0 资产
  探索性（明确标注）：Torin 救援臂对比 LPS_CS_Torin vs LPS_CS（未预注册）
输入：00_RAW_DATA/GSE235046/GSE235046_Count_table.txt/GSE235046_Count_table.txt
输出：_intermediate/P1_GSE235046_interaction_DE.csv
     04_AUDIT_GOVERNANCE/P1_IIAMD_signature_v1.0.csv（冻结资产）/ P1_Interaction_Signature_Report.md
     03_LOGS/P1_interaction_signature_log.txt
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
LOGF = os.path.join(ROOT, "03_LOGS", "P1_interaction_signature_log.txt")

log = []
def note(m=""):
    log.append(m); print(m, flush=True)

# ---------------- 1. 读入与映射 ----------------
note("== 1. 读入 count 表与组映射")
raw = pd.read_csv(CNT, sep="\t")
note(f"   原始 shape={raw.shape}")
meta_cols = ["geneID", "geneSymbol", "bioType", "annotationLevel"]
samples = [c for c in raw.columns if c not in meta_cols]
note(f"   样本列 {len(samples)}: {samples[:3]} ... {samples[-3:]}")

YW_GROUP = {}
for i, yw in enumerate(["YW001","YW002","YW003"], 1): YW_GROUP[yw] = "Media"
for i, yw in enumerate(["YW004","YW005","YW006"], 1): YW_GROUP[yw] = "LPS"
for i, yw in enumerate(["YW007","YW008","YW009"], 1): YW_GROUP[yw] = "CS"
for i, yw in enumerate(["YW010","YW011","YW012"], 1): YW_GROUP[yw] = "LPS_CS"
for i, yw in enumerate(["YW016","YW017","YW018"], 1): YW_GROUP[yw] = "LPS_CS_Torin"

def col_group(c):
    yw = c.split("_")[-1]
    return YW_GROUP.get(yw)
grp_map = {c: col_group(c) for c in samples}
missing = [c for c, g in grp_map.items() if g is None]
if missing:
    note(f"   [ERROR] 未映射样本列: {missing}")
    sys.exit(1)
note("   映射: " + str({g: [c.split('_')[-1] for c in samples if grp_map[c] == g] for g in sorted(set(grp_map.values()))}))

# 按基因符号聚合（重复符号求和；空符号丢弃）
raw["geneSymbol"] = raw["geneSymbol"].astype(str).str.strip()
ok = raw[raw["geneSymbol"].notna() & (raw["geneSymbol"] != "") & (raw["geneSymbol"] != "nan")]
mat = ok.groupby("geneSymbol")[samples].sum()
mat = mat.astype(int)
note(f"   按符号聚合后: {mat.shape[0]} 基因 × {mat.shape[1]} 样本")

# ---------------- 2. 映射方向校验（LPS 标记基因） ----------------
note("\n== 2. LPS 标记基因方向校验（映射 PASS 才继续）")
LPS_MARKERS = ["Tnf", "Il1b", "Cxcl10", "Icam1", "Nfkbia", "Ifit3"]
CS_MARKERS = ["Trib3", "Ddit3", "Asns", "Slc7a5"]
gmean = mat.T.groupby([grp_map[c] for c in mat.columns]).mean().T  # 基因×组 均值
hits = 0
for g in LPS_MARKERS:
    if g in gmean.index:
        m, l, lc = gmean.loc[g, "Media"], gmean.loc[g, "LPS"], gmean.loc[g, "LPS_CS"]
        up = (l > m) and (lc > m)
        hits += up
        note(f"   {g}: Media={m:.1f} LPS={l:.1f} LPS_CS={lc:.1f} → {'诱导✓' if up else '未诱导✗'}")
    else:
        note(f"   {g}: 不在表内")
for g in CS_MARKERS:
    if g in gmean.index:
        note(f"   [CS参考] {g}: Media={gmean.loc[g,'Media']:.1f} CS={gmean.loc[g,'CS']:.1f}")
if hits < 3:
    note(f"   [ERROR] LPS 标记诱导仅 {hits}/{len([g for g in LPS_MARKERS if g in gmean.index])}——映射可疑，终止")
    sys.exit(1)
note(f"   映射校验 PASS（{hits} 个 LPS 标记在 LPS 与 LPS_CS 均诱导）")

# ---------------- 3. 确证交互模型 ----------------
note("\n== 3. 确证 2×2 交互模型（12 样本；pydeseq2 ~group + 交互对比向量）")
conf_cols = [c for c in samples if grp_map[c] in ("Media", "LPS", "CS", "LPS_CS")]
cmat = mat[conf_cols]
# 预过滤：≥10 counts 在 ≥3 样本（标准做法，披露）
keep = (cmat >= 10).sum(axis=1) >= 3
cmatf = cmat[keep]
note(f"   预过滤后 {cmatf.shape[0]} 基因（≥10 counts 于 ≥3 样本）")

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
counts_df = cmatf.T  # 样本×基因
meta = pd.DataFrame({"group": [grp_map[c] for c in cmatf.columns]}, index=cmatf.columns)
dds = DeseqDataSet(counts=counts_df, metadata=meta, design="~group", quiet=True, n_cpus=1)
dds.deseq2()
cols = list(dds.obsm["design_matrix"].columns)
note(f"   设计列: {cols}")
ref = "CS" if "group[T.Media]" in cols and "group[T.CS]" not in cols else None
# 通用：确定参考水平
tlevels = [c.split("T.")[1].rstrip("]") for c in cols if c.startswith("group[T.")]
ref = ({"Media", "LPS", "CS", "LPS_CS"} - set(tlevels)).pop()
def mu_vec(level):
    v = np.zeros(len(cols)); v[cols.index("Intercept")] = 1.0
    if level != ref and f"group[T.{level}]" in cols:
        v[cols.index(f"group[T.{level}]")] = 1.0
    return v
cvec = mu_vec("LPS_CS") - mu_vec("LPS") - mu_vec("CS") + mu_vec("Media")
note(f"   参考水平={ref}；交互对比向量={cvec.tolist()}")
st = DeseqStats(dds, contrast=cvec, quiet=True, n_cpus=1)
st.summary()
res = st.results_df.reset_index().rename(columns={"index": "geneSymbol"})
res = res.sort_values("padj")
res.to_csv(os.path.join(INTER, "P1_GSE235046_interaction_DE.csv"), index=False)

sig_up = res[(res["padj"] < 0.05) & (res["log2FoldChange"] > 0)].copy()
sig_dn = res[(res["padj"] < 0.05) & (res["log2FoldChange"] < 0)].copy()
note(f"   交互显著：up={len(sig_up)}, down={len(sig_dn)}（padj<0.05）")

# ---------------- 4. 签名冻结 ----------------
sig = pd.concat([
    sig_up.assign(direction="up"), sig_dn.assign(direction="down")
])[["geneSymbol", "direction", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"]]
sig.insert(0, "signature_version", "IIAMD_v1.0")
sig.insert(1, "frozen_date", "2026-08-20")
sig.insert(2, "source", "GSE235046 LPS×CS interaction (DESeq2 Wald, contrast vector)")
sig.to_csv(SIGF, index=False)
note(f"   冻结签名 → {SIGF}（up {len(sig_up)} / down {len(sig_dn)}）")

# 机制节点覆盖检查（描述性）
ANCHOR = ["Bax","Bak1","Bid","Rictor","Rhoa","Ninj1","Rptor","Mtor","Akt1","Tlr4","Myd88","Tnf"]
cov = {g: (g in gmean.index) for g in ANCHOR}
note(f"   锚定节点在表内: {cov}")

# ---------------- 5. 探索性：Torin 救援臂（未预注册，仅描述） ----------------
note("\n== 5. 探索性：LPS_CS_Torin vs LPS_CS（未预注册）")
try:
    exp_cols = [c for c in samples if grp_map[c] in ("LPS_CS", "LPS_CS_Torin")]
    emat = mat.loc[cmatf.index, exp_cols]
    emeta = pd.DataFrame({"group": [grp_map[c] for c in exp_cols]}, index=exp_cols)
    edds = DeseqDataSet(counts=emat.T, metadata=emeta, design="~group", quiet=True, n_cpus=1)
    edds.deseq2()
    torin_lvl = "LPS_CS_Torin"
    tcols = list(edds.obsm["design_matrix"].columns)
    ttlevels = [c.split("T.")[1].rstrip("]") for c in tcols if c.startswith("group[T.")]
    tref = ({"LPS_CS", "LPS_CS_Torin"} - set(ttlevels)).pop()
    st2 = DeseqStats(edds, contrast=["group", torin_lvl, tref], quiet=True, n_cpus=1)
    st2.summary()
    r2 = st2.results_df.reset_index().rename(columns={"index": "geneSymbol"})
    # IIAMD 签名基因在 Torin 下的方向（描述：救援=up 签名被 Torin 拉回 0/反向）
    merged = sig.merge(r2[["geneSymbol", "log2FoldChange"]], on="geneSymbol", suffixes=("", "_torin"))
    reverted = ((merged["direction"] == "up") & (merged["log2FoldChange_torin"] < 0)) | \
               ((merged["direction"] == "down") & (merged["log2FoldChange_torin"] > 0))
    note(f"   Torin 对签名基因：{int(reverted.sum())}/{len(merged)} 个方向逆转（探索性）")
    merged.to_csv(os.path.join(INTER, "P1_exploratory_Torin_signature_effect.csv"), index=False)
except Exception as e:
    note(f"   探索性 Torin 对比失败（{type(e).__name__}: {e}）——不影响确证结果")

# ---------------- 6. Gate 1 判定 + 报告 ----------------
gate1 = (len(sig_up) + len(sig_dn)) > 0
note(f"\n== 6. Gate 1 初判：{'存在非空 LPS×CS 交互签名（待 LORO 稳定性确认）' if gate1 else '无显著交互——不构建 Mitoxyperilysis 转录签名'}")

with open(REPORT, "w", encoding="utf-8") as f:
    f.write("# P1-2 GSE235046 交互签名报告（H6，Gate 1）\n\n")
    f.write("- 日期：2026-08-20；冻结依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md（OSF https://osf.io/C7RYD/）\n")
    f.write("- 确证模型：12 样本（Media/LPS/CS/LPS+CS 各 3），DESeq2（pydeseq2 0.5.4）~group + 交互对比向量\n")
    f.write(f"  β₃ = μ_LPS_CS − μ_LPS − μ_CS + μ_Media；参考水平={ref}\n")
    f.write(f"- 预过滤：≥10 counts 于 ≥3 样本，剩 {cmatf.shape[0]} 基因（标准做法，披露）\n")
    f.write(f"- 映射校验：LPS 标记 {hits} 个诱导 PASS\n\n")
    f.write(f"## 结果：IIAMD-up = {len(sig_up)}，IIAMD-down = {len(sig_dn)}（padj<0.05）\n\n")
    f.write("### IIAMD-up（前 25，按 padj）\n\n| gene | log2FC(β₃) | padj |\n|---|---|---|\n")
    for _, r in sig_up.head(25).iterrows():
        f.write(f"| {r['geneSymbol']} | {r['log2FoldChange']:+.2f} | {r['padj']:.2e} |\n")
    f.write("\n### IIAMD-down（前 25，按 padj）\n\n| gene | log2FC(β₃) | padj |\n|---|---|---|\n")
    for _, r in sig_dn.head(25).iterrows():
        f.write(f"| {r['geneSymbol']} | {r['log2FoldChange']:+.2f} | {r['padj']:.2e} |\n")
    f.write(f"\n## 机制锚定节点覆盖（描述性）\n\n{cov}\n")
    f.write(f"\n## Gate 1 初判\n\n{'非空签名存在；正式判定待 leave-one-replicate-out 稳定性（P1-2b）。' if gate1 else '空签名：按冻结计划，不构建 Mitoxyperilysis 转录签名，只检验原始机制节点。'}\n")
    f.write("\n## 待办（P2 前置）\n\n- 小鼠→人一对一 ortholog 映射（Ensembl BioMart 导出，需人工下载）\n")
    f.write("- ST002738 代谢组锚定（mwTab 注释解析；描述性）\n")

with open(LOGF, "w", encoding="utf-8") as f:
    f.write("\n".join(log))
print("\nDONE P1-2 ->", REPORT)
