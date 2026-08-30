#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bridge_test_narrative.py — 桥接检验：为"Mitoxyperilysis 主题冲 IF10"定叙事骨架。
回答三个决定框架的问题：
  TEST 1  8 模块在 3 个平台（单细胞 BALF / 肺组织 bulk / 外周血）的转录方向是否一致？（重塑骨干是否扎实）
  TEST 2  MAM-铁-呼吸链节点的"多组学收敛打分" → 提名锚点基因（gnomAD约束+AM变异负荷+SCENIC中枢性+转录重塑+可成药）
  TEST 3  死亡【执行者】基因 vs 上游【重塑】基因 的方向对比 → 数据支持"正在执行死亡"还是"已布线但未执行"？

N1 修复（2026-08-15）：外部数据层与权威 80 基因清单的符号连接现经 HGNC 别名归一化
  （ITPR1→IP3R1、HSPD1→HSP60）。既往版本按精确符号 merge，导致 IP3R1 行 pLI/N_Pathogenic
  全空（gnomAD/AM 层实际以 ITPR1 记录：pLI=1.00、N_Pathogenic=10,614）、HSP60 行 N_Pathogenic
  全空（AM 层以 HSPD1 记录：2,003）。修复后模块排序翻转：MAM 0.458 > iron 0.439
  （修复前 iron 0.455 > MAM 0.441），IP3R1 收敛得分 0.304→0.656 跃居第 2。
  新增 TEST 2b 数据可得性敏感性（稀疏维度检验）与 alias_connected/n_real_dims/
  convergence_sens_available 三列。详见 04_AUDIT_GOVERNANCE/N1_Bridge_Alias_Reconnection_Report.md。
"""
import os, sys, re, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass

T = "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV"
OUT = "04_AUDIT_GOVERNANCE/Bridge_Test_Result.csv"

def L(*a): print(*a, flush=True)

# ---- 80 通路基因 + 模块（来自 S2） ----
# N3 修复（2026-08-15）：单细胞转录维度数据源改为 Table_S2b（真·单细胞 DEG，
#   在修正标签后的 combined_processed.h5ad 上重算，GSE145926 BALF severe vs Healthy）。
#   既往版本误将 Table_S2_DEGs_Analysis.csv 的 ARDS_vs_Control 列当作"单细胞BALF"，
#   经数值指纹核对（ACTB≈13.0、EPCAM≈0.5、样本设计 82/44）该表实为
#   GSE185263 全血 Bulk DESeq2 结果。d_mean 为 log1p 尺度均值差（自然对数），
#   与 bulk log2FC 尺度不同，但 ranknorm 仅用相对排序，不影响收敛打分口径。
s2b = pd.read_csv(f"{T}/Table_S2b_scRNA_DEG_80genes.csv")
s2b["gene"] = s2b["gene_symbol"].astype(str)
s2b["module"] = s2b["Mitoxyperilysis_module"].str.replace("mitoxy_", "", regex=False)
g = s2b[["gene", "module", "GSE145926_BALF_d_mean", "GSE145926_BALF_bh_padj"]].rename(
    columns={"GSE145926_BALF_d_mean": "sc_log2FC", "GSE145926_BALF_bh_padj": "sc_padj"})
L(f"[S2b] 通路基因（单细胞 BALF severe vs Healthy, GSE145926）: {len(g)} 个，模块: {sorted(g['module'].unique())}")

# ---- N1 修复（2026-08-15）：HGNC 别名归一化后再 merge ----
# 权威清单用 IP3R1/HSP60，外部层（gnomAD/AlphaMissense）以 HGNC 官方符号 ITPR1/HSPD1 记录
ALIAS = {"ITPR1": "IP3R1", "HSPD1": "HSP60"}   # 外部符号 -> 权威清单符号
alias_log = {}
def merge_ext(g, df, gene_col, val_cols, src):
    df = df.copy()
    raw = df[gene_col].astype(str)
    canon = raw.replace(ALIAS)
    hit = raw[raw != canon]
    for ext in hit.unique():
        alias_log.setdefault(canon[raw == ext].iloc[0], []).append(f"{ext}({src})")
    df["gene"] = canon
    return g.merge(df[["gene"] + val_cols], on="gene", how="left")

# ---- 全血 bulk（S46a: Sepsis_COVID vs Control，GSE185263；R1 正名） ----
a = pd.read_csv(f"{T}/Table_S46a_Circulating_Mitoxy_Genes_GSE185263.csv")
g = merge_ext(g, a, "Gene", ["COVID_vs_Control_log2FC"], "S46a").rename(columns={"COVID_vs_Control_log2FC": "GSE185263_log2FC"})

# ---- 外周血微阵列（S46b: SDRA vs Control，GSE212865） ----
b = pd.read_csv(f"{T}/Table_S46b_Circulating_Mitoxy_Genes_GSE212865.csv")
g = merge_ext(g, b, "Gene", ["SDRA_vs_Control_log2FC"], "S46b").rename(columns={"SDRA_vs_Control_log2FC": "GSE212865_log2FC"})

# ---- gnomAD 约束（S40） ----
s40 = pd.read_csv(f"{T}/Table_S40_gnomAD_Constraint.csv")
s40["pLI"] = pd.to_numeric(s40["pLI"], errors="coerce")
g = merge_ext(g, s40, "Gene_Symbol_HGNC", ["pLI"], "gnomAD")

# ---- AlphaMissense 变异负荷（S18） ----
s18 = pd.read_csv(f"{T}/Table_S18_AlphaMissense_Main.csv")
s18["N_Pathogenic"] = pd.to_numeric(s18["N_Pathogenic"], errors="coerce")
g = merge_ext(g, s18, "Gene_Symbol", ["N_Pathogenic"], "AlphaMissense")

# ---- SCENIC 调控中枢性（S12: 每个靶基因被多少 TF 调控） ----
s12 = pd.read_csv(f"{T}/Table_S12_SCENIC_Regulons.csv")
deg = s12[s12["Is_Mitoxyperilysis_Target"] == True].groupby("Target_Gene").size().reset_index(name="TF_degree")
g = merge_ext(g, deg, "Target_Gene", ["TF_degree"], "SCENIC")
g["TF_degree"] = g["TF_degree"].fillna(0)

# ---- 可成药（S36 对接靶蛋白） ----
s36 = pd.read_csv(f"{T}/Table_S36_Docking_Results.csv")
targets = set(s36["Target_Protein"].astype(str).unique())
g["druggable"] = g["gene"].isin(targets).astype(int)

# =====================================================================
L("\n" + "=" * 78); L("TEST 1 | 8 模块 × 3 平台 转录方向一致性"); L("=" * 78)
plat = {"sc_log2FC": "单细胞BALF", "GSE185263_log2FC": "脓毒症全血Bulk", "GSE212865_log2FC": "外周血微阵列"}
rows = []
for mod, sub in g.groupby("module"):
    r = {"module": mod, "n": len(sub)}
    dirs = []
    for col, lab in plat.items():
        m = sub[col].mean()
        r[lab] = round(m, 2)
        dirs.append(np.sign(m))
    r["方向一致"] = "是" if len(set(dirs)) == 1 else "否"
    rows.append(r)
t1 = pd.DataFrame(rows).sort_values("n", ascending=False)
L(t1.to_string(index=False))
consist_modules = t1[t1["方向一致"] == "是"]["module"].tolist()
L(f"\n→ 三平台方向一致的模块: {len(consist_modules)}/{len(t1)} = {consist_modules}")
L("→ 呼吸链/线粒体功能模块若三平台均为负（下调），即为最稳的重塑骨干。")

# =====================================================================
L("\n" + "=" * 78); L("TEST 2 | 多组学收敛打分 → 锚点基因提名"); L("=" * 78)
def ranknorm(s):
    s = s.rank(pct=True)
    return s.fillna(s.min() / 2)
g["r_transcript"] = ranknorm(g[["sc_log2FC", "GSE185263_log2FC", "GSE212865_log2FC"]].abs().mean(axis=1))
g["r_constraint"] = ranknorm(g["pLI"])
g["r_burden"]     = ranknorm(g["N_Pathogenic"])
g["r_central"]    = ranknorm(g["TF_degree"])
g["r_drug"]       = g["druggable"].astype(float)
g["convergence"] = g[["r_transcript", "r_constraint", "r_burden", "r_central", "r_drug"]].mean(axis=1)

# ---- N1 新增列：别名连接记录 + 真实维度计数 ----
g["n_real_dims"] = 3 + g["pLI"].notna().astype(int) + g["N_Pathogenic"].notna().astype(int)
g["alias_connected"] = g["gene"].map(lambda x: "; ".join(alias_log.get(x, []))).replace("", np.nan)
top = g.sort_values("convergence", ascending=False).head(15)[
    ["gene", "module", "convergence", "sc_log2FC", "GSE185263_log2FC", "GSE212865_log2FC", "pLI", "N_Pathogenic", "TF_degree", "druggable"]]
L("Top 15 收敛得分基因:"); L(top.to_string(index=False))
L("\n→ 收敛得分最高的模块分布（锚点候选）:")
L(g.groupby("module")["convergence"].mean().sort_values(ascending=False).round(3).to_string())
L(f"\n→ [N1] 别名重连: {alias_log}")
L(f"[N1] 五维均有真实数据的基因: {g.loc[g['n_real_dims']==5,'gene'].tolist()}")

# =====================================================================
L("\n" + "=" * 78); L("TEST 2b | N1b 数据可得性敏感性（稀疏维度检验）"); L("=" * 78)
avail = g.groupby("module").agg(n=("gene","size"), pLI有数据=("pLI", lambda x: x.notna().sum()),
    AM有数据=("N_Pathogenic", lambda x: x.notna().sum()), docking=("druggable", lambda x: (x==1).sum()))
L("模块×维度可得性:"); L(avail.to_string())
g["conv_fullcov"] = g[["r_transcript", "r_central", "r_drug"]].mean(axis=1)
L("\n敏感性 A（仅全覆盖维度 transcript+central+drug，无填充默认）模块均值:")
L(g.groupby("module")["conv_fullcov"].mean().sort_values(ascending=False).round(4).to_string())
def conv_avail(row):
    v = {"r_transcript": row["r_transcript"], "r_central": row["r_central"], "r_drug": row["r_drug"]}
    if pd.notna(row["pLI"]): v["r_constraint"] = row["r_constraint"]
    if pd.notna(row["N_Pathogenic"]): v["r_burden"] = row["r_burden"]
    return np.mean(list(v.values()))
g["convergence_sens_available"] = g.apply(conv_avail, axis=1)
sub = g[g["n_real_dims"] >= 4]
L(f"\n敏感性 B（逐基因仅用真实维度，≥4/5 维有数据，n={len(sub)}/80）模块均值:")
L(sub.groupby("module")["convergence_sens_available"].mean().sort_values(ascending=False).round(4).to_string())
L("\n→ 解读: 约束维度仅 4 个 MAM 基因有 gnomAD 数据、变异负荷 36/80——")
L("  若剔除稀疏维度后铁代谢反居首，则 MAM/铁层内先后由可得性驱动；")
L("  稳健结论应为'MAM–铁共同构成收敛首位层级、远高于执行端'。")

# =====================================================================
L("\n" + "=" * 78); L("TEST 3 | 死亡【执行者】 vs 上游【重塑】（决定'执行'还是'布线'框架）"); L("=" * 78)
EXEC = {"cell_death"}                      # CASP/GSDM/PYCARD/NLRP3/IL1B/CYCS/APAF1 等执行机器
UPSTREAM = {"MAM_integrity", "iron_metabolism", "mitochondrial_function", "oxidative_stress"}
def block_stats(sub):
    out = {}
    for col, lab in plat.items():
        v = sub[col].dropna()
        out[lab + "_mean"] = round(v.mean(), 2)
        out[lab + "_%up"] = round((v > 0).mean() * 100)
    return out
exec_g = g[g["module"].isin(EXEC)]
up_g = g[g["module"].isin(UPSTREAM)]
L(f"执行者模块(cell_death, n={len(exec_g)}): {block_stats(exec_g)}")
L(f"上游重塑(MAM+铁+线粒体+氧化, n={len(up_g)}): {block_stats(up_g)}")
# 判别指数：上游 |log2FC| 中位 vs 执行者 |log2FC| 中位（脓毒症全血bulk，R1 正名）
up_mag = up_g["GSE185263_log2FC"].abs().median()
ex_mag = exec_g["GSE185263_log2FC"].abs().median()
L(f"\n脓毒症全血bulk |log2FC| 中位数: 上游={up_mag:.2f}  执行者={ex_mag:.2f}")
ex_uprate = (exec_g["GSE185263_log2FC"] > 0).mean()
L(f"执行者基因在脓毒症全血bulk中上调比例: {ex_uprate*100:.0f}%")
L("\n→ 解读: 若上游 |log2FC| 明显大于执行者，且执行者未呈现协调强上调，")
L("  则数据支持'线路已被重塑/布线，但未在 bulk 级别主动执行死亡'——即可用")
L("  'primed death-competent state（已启动的死亡胜任态）'框架防守 IF10。")

# 保存（N1：新增 alias_connected / n_real_dims / convergence_sens_available / conv_fullcov 四列；
#   conv_fullcov 为敏感性 A 的逐基因值，供复查；主 convergence 口径与既往版本一致）
g = g.sort_values("convergence", ascending=False)
SAVE_COLS = ["gene", "module", "sc_log2FC", "sc_padj", "GSE185263_log2FC", "GSE212865_log2FC",
             "pLI", "N_Pathogenic", "TF_degree", "druggable",
             "r_transcript", "r_constraint", "r_burden", "r_central", "r_drug", "convergence",
             "alias_connected", "n_real_dims", "convergence_sens_available", "conv_fullcov"]
g[SAVE_COLS].to_csv(OUT, index=False, encoding="utf-8-sig")
L(f"\n[已保存] {OUT}")
