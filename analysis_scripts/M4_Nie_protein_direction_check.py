# -*- coding: utf-8 -*-
"""
M4 蛋白层方向一致性（修正参照）：
参照1（文档口径，全血）: manifest ARDS_vs_Control_log2FC (GSE185263 全血 bulk)
参照2（组织匹配）: Table_S2b GSE145926 BALF scRNA d_mean（重症-健康）
蛋白层: Nie2021 尸检肺 mmc4 log2FC（缺测基因用 mmc3 C/N 比值补）
"""
import pandas as pd
import io, sys, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
import numpy as np
from scipy.stats import spearmanr

BASE = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
DL = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\00_RAW_DATA\Nie2021_Cell_MultiOrgan_Proteomics"

man = pd.read_csv(BASE + r"\04_AUDIT_GOVERNANCE\Mitoxyperilysis_Gene_Manifest_v1.0.csv")
s2b = pd.read_csv(BASE + r"\02_SUPPLEMENTARY_TABLES\SUPPLEMENTARY_Tables_CSV\Table_S2b_scRNA_DEG_80genes.csv")
mmc3 = pd.read_csv(BASE + r"\01_FIGURE_DATA_CSV\M4_Nie_lung_mmc3_key_genes.csv")

# 蛋白方向：优先 mmc4（论文统计），缺测用 mmc3 C/N（粗）
def px(gene):
    hit = mmc3[mmc3["gene"] == gene]
    return hit["C_N_ratio"].iloc[0] if len(hit) else np.nan

# 读 mmc4 lung（复用上一轮结果文件更稳）
mmc4_file = BASE + r"\01_FIGURE_DATA_CSV\M4_Nie_lung_protein_vs_transcript.csv"
if os.path.exists(mmc4_file):
    prev = pd.read_csv(mmc4_file)
    prev = prev[["gene", "protein_log2FC(Nie_lung)", "protein_padj"]]
    prev.columns = ["gene", "px_log2FC", "px_padj"]
else:
    prev = pd.DataFrame(columns=["gene", "px_log2FC", "px_padj"])

rows = []
for _, r in man.iterrows():
    g = r["gene_symbol"]
    p = prev[prev["gene"] == g]
    if len(p):
        px_fc, px_padj = p["px_log2FC"].iloc[0], p["px_padj"].iloc[0]
        px_src = "mmc4"
    else:
        ratio = px(g)
        px_fc = np.log2(ratio) if ratio and ratio > 0 else np.nan
        px_padj = np.nan
        px_src = "mmc3_ratio"
    b = s2b[s2b["gene_symbol"] == g]
    bal = b["GSE145926_BALF_d_mean"].iloc[0] if len(b) else np.nan
    bal_p = b["GSE145926_BALF_bh_padj"].iloc[0] if len(b) else np.nan
    tx = r["ARDS_vs_Control_log2FC"]
    tx_p = r["ARDS_vs_Control_padj"]
    rows.append({"gene": g, "arm": r["arm"],
                 "blood_tx_log2FC": tx, "blood_tx_padj": tx_p,
                 "balf_d": bal, "balf_padj": bal_p,
                 "px_log2FC": px_fc, "px_padj": px_padj, "px_src": px_src})

df = pd.DataFrame(rows)
df["blood_tx_log2FC"] = pd.to_numeric(df["blood_tx_log2FC"], errors="coerce")
df["balf_d"] = pd.to_numeric(df["balf_d"], errors="coerce")
df["px_log2FC"] = pd.to_numeric(df["px_log2FC"], errors="coerce")

def sign_series(x):
    return np.sign(x).astype(float)

print("=" * 66)
print("口径 A（文档原口径：全血转录 vs 尸检肺蛋白）")
a = df.dropna(subset=["px_log2FC", "blood_tx_log2FC"]).copy()
a = a[a["blood_tx_padj"] < 0.05]
a["conc"] = sign_series(a["blood_tx_log2FC"]) == sign_series(a["px_log2FC"])
print(f"  可判定 n={len(a)}, 一致 {a['conc'].sum()} ({a['conc'].mean()*100:.0f}%)")
rho, pv = spearmanr(a["blood_tx_log2FC"], a["px_log2FC"])
print(f"  Spearman rho={rho:.3f}, p={pv:.3g}")

print("\n口径 B（组织匹配：BALF scRNA vs 尸检肺蛋白）")
b = df.dropna(subset=["px_log2FC", "balf_d"]).copy()
b = b[b["balf_padj"] < 0.05]
b["conc"] = sign_series(b["balf_d"]) == sign_series(b["px_log2FC"])
print(f"  可判定 n={len(b)}, 一致 {b['conc'].sum()} ({b['conc'].mean()*100:.0f}%)")
rho2, pv2 = spearmanr(b["balf_d"], b["px_log2FC"])
print(f"  Spearman rho={rho2:.3f}, p={pv2:.3g}")

print("\n===== 口径 B 明细（按 arm）=====")
for arm in ["upstream_collapse", "execution_induction", "not_in_dissociation_arms"]:
    sub = b[b["arm"] == arm]
    if len(sub) == 0:
        continue
    c = sub["conc"].sum()
    print(f"\n[{arm}] n={len(sub)}, 一致={c} ({c/len(sub)*100:.0f}%)")
    print(sub[["gene", "balf_d", "balf_padj", "px_log2FC", "px_src", "conc"]]
          .to_string(index=False, max_colwidth=14))

df.to_csv(BASE + r"\01_FIGURE_DATA_CSV\M4_Nie_protein_direction_two_refs.csv", index=False)
print("\n已保存 M4_Nie_protein_direction_two_refs.csv")
