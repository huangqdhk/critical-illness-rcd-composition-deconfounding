# -*- coding: utf-8 -*-
"""
M4 蛋白层：Nie2021 Cell 肺组织蛋白组 × 80 基因 manifest 方向一致性。
- mmc4 'A_COVID vs Non-COVID' 的 Lung log2FC（COVID vs 非COVID，正=重症更高）
- manifest 的 ARDS_vs_Control_log2FC（转录层，正=ARDS 更高）
- 口径：跨队列（COVID肺 vs ARDS肺）、跨分子层（蛋白 vs 转录）的方向一致性
"""
import pandas as pd
import io, sys, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
import numpy as np
from scipy.stats import spearmanr

BASE = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
DL = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\00_RAW_DATA\Nie2021_Cell_MultiOrgan_Proteomics"
man = pd.read_csv(BASE + r"\04_AUDIT_GOVERNANCE\Mitoxyperilysis_Gene_Manifest_v1.0.csv")

# ---- 1. mmc4: 多级表头，取 Lung 三列 ----
d = pd.read_excel(os.path.join(DL, "mmc4.xlsx"), sheet_name="A_COVID vs Non-COVID",
                  header=None, nrows=5400)
hdr0 = d.iloc[0].tolist()
hdr1 = d.iloc[1].tolist()
# 找 Lung 列的起止
lung_idx = [i for i, h in enumerate(hdr0) if str(h).strip() == "Lung"]
assert lung_idx, "未找到 Lung 列"
i0 = lung_idx[0]
cols = []
for j in range(i0, i0 + 3):
    cols.append(f"{hdr1[j]}_{j}")
lung = d.iloc[2:].iloc[:, i0:i0 + 3].copy()
lung.columns = ["log2FC", "p", "padj"]
lung.insert(0, "gene", d.iloc[2:, 1].astype(str))
lung["log2FC"] = pd.to_numeric(lung["log2FC"], errors="coerce")
lung["padj"] = pd.to_numeric(lung["padj"], errors="coerce")
lung["p"] = pd.to_numeric(lung["p"], errors="coerce")
print("mmc4 lung rows:", len(lung), "| 有 log2FC 的:", lung["log2FC"].notna().sum())

# 一个基因可能多行（多 Uniprot ID）→ 取 |log2FC| 最大的一行作代表，并记录行数
g = lung.dropna(subset=["log2FC"]).groupby("gene", as_index=False).apply(
    lambda x: x.loc[x["log2FC"].abs().idxmax()], include_groups=False).reset_index(drop=True)
g["n_rows"] = lung.dropna(subset=["log2FC"]).groupby("gene").size().reindex(g["gene"]).values
print("去重后基因数:", len(g))

# ---- 2. 与 manifest 对齐 ----
m = man[["gene_symbol", "arm", "ARDS_vs_Control_log2FC", "ARDS_vs_Control_padj"]].copy()
m.columns = ["gene", "arm", "tx_log2FC", "tx_padj"]
m["tx_log2FC"] = pd.to_numeric(m["tx_log2FC"], errors="coerce")
m["tx_padj"] = pd.to_numeric(m["tx_padj"], errors="coerce")
j = m.merge(g, on="gene", how="left")
print("\n80 基因中在肺蛋白组测到的:", j["log2FC"].notna().sum(), "/80")

# ---- 3. 方向一致率 ----
j = j.dropna(subset=["log2FC"]).copy()
j["dir_tx"] = np.sign(j["tx_log2FC"])
j["dir_px"] = np.sign(j["log2FC"])
j["concordant"] = j["dir_tx"] == j["dir_px"]
# 转录层不显著的基因不参与方向判定（避免噪音）
j_judge = j[j["tx_padj"] < 0.05].copy()
n = len(j_judge)
n_conc = j_judge["concordant"].sum()
print(f"\n可判定基因（转录层 padj<0.05 且有蛋白数据）: {n}")
print(f"方向一致: {n_conc} | 一致率: {n_conc/n*100:.1f}%")
rho, sp = spearmanr(j["tx_log2FC"], j["log2FC"])
print(f"全可测基因 Spearman rho = {rho:.3f} (p={sp:.3g}, n={len(j)})")

# 蛋白层显著的基因的一致率
j_sig = j[j["padj"] < 0.05]
if len(j_sig):
    c2 = (np.sign(j_sig["tx_log2FC"]) == np.sign(j_sig["log2FC"])).sum()
    print(f"蛋白层显著(padj<0.05)基因: {len(j_sig)} 个, 其中与转录同向: {c2} ({c2/len(j_sig)*100:.1f}%)")

# ---- 4. 明细表 ----
out = j[["gene", "arm", "tx_log2FC", "tx_padj", "log2FC", "p", "padj", "dir_tx", "dir_px", "concordant", "n_rows"]]
out.columns = ["gene", "arm", "transcript_log2FC", "transcript_padj",
               "protein_log2FC(Nie_lung)", "protein_p", "protein_padj",
               "dir_tx", "dir_protein", "concordant", "protein_rows"]
out = out.sort_values("protein_padj")
out.to_csv(BASE + r"\01_FIGURE_DATA_CSV\M4_Nie_lung_protein_vs_transcript.csv", index=False)
print("\n===== 按 arm 分组的明细 =====")
for arm in ["upstream_collapse", "execution_induction", "not_in_dissociation_arms"]:
    sub = out[out["arm"] == arm].dropna(subset=["protein_log2FC(Nie_lung)"])
    if len(sub) == 0:
        continue
    cc = (sub["dir_tx"] == sub["dir_protein"]).sum()
    print(f"\n[{arm}] n={len(sub)}, 一致={cc} ({cc/len(sub)*100:.0f}%)")
    print(sub.to_string(index=False, max_colwidth=18))
