#!/usr/bin/env python3
"""
M3 步骤9：输出组装与判定
========================
1) 生成 Figure_8I.csv（pQTL-MR 森林图数据）与 Figure_8J.csv（SMR/TWAS 联合数据）
2) 全部新表 SHA256 冻结，追加 RESULTS_MANIFEST_v2.0.csv 与 FIGURE_DATA_MANIFEST.csv
3) 撰写 _intermediate/M3_judgment_report.md（预注册判定）
判定（预注册 M3_pre_registration_20260818.md）：
  R1 = pQTL-MR IVW BH<0.05 + PP.H4>=0.8，或 SMR BH<0.05 且 HEIDI 通过 -> 因果锚点升级
  R3 = 四层（eQTL/pQTL/SMR/TWAS）均不支持 -> 严谨排除写入 Discussion
"""
import hashlib
import sys
from pathlib import Path

import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = Path(r"E:\SCI\proj")
OUT = BASE / "02_SUPPLEMENTARY_TABLES" / "SUPPLEMENTARY_Tables_CSV"
FIG = BASE / "01_FIGURE_DATA_CSV" / "Main"
INT = BASE / "_intermediate"
GOV = BASE / "04_AUDIT_GOVERNANCE"
GENE_SET_VERSION = "Mitoxy-80_v1.0"


def sha256(f):
    h = hashlib.sha256()
    with open(f, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_fig8i():
    mr = pd.read_csv(OUT / "Table_S15e_pQTL_MR_Results.csv")
    cols = ["exposure", "outcome", "method", "nsnp", "b", "se", "pval", "OR",
            "OR_CI_low", "OR_CI_high", "bh_p", "bh_p_pooled", "Q", "Q_pval", "F_stat"]
    df = mr[cols].copy()
    df["layer"] = "pQTL-MR"
    df["gene_set_version"] = GENE_SET_VERSION
    df["score_version"] = "m3pqtl_v1.0"
    df.to_csv(FIG / "Figure_8I.csv", index=False)
    return df


def build_fig8j():
    smr = pd.read_csv(OUT / "Table_S15f_SMR_HEIDI_Results.csv")
    tw = pd.read_csv(OUT / "Table_S15g_TWAS_Results.csv")
    s = smr[["layer", "gene_symbol", "top_snp", "n_snps_matched", "top_p_eqtl",
             "b_smr", "se_smr", "p_smr", "bh_p_pooled", "bh_p_layer",
             "heidi_T", "heidi_df", "heidi_p", "heidi_pass"]].copy()
    s["method"] = "SMR"
    s["score_version"] = "m3smr_v1.0"
    t = tw[["tissue", "gene", "model", "method", "n_snps_weight",
            "n_snps_matched_gwas", "n_snps_ld", "twas_z", "twas_p",
            "bh_p_pooled", "bh_p_tissue"]].copy()
    t["score_version"] = "m3twas_v1.0"
    s["gene_set_version"] = GENE_SET_VERSION
    t["gene_set_version"] = GENE_SET_VERSION
    df = pd.concat([s, t], ignore_index=True, sort=False)
    df.to_csv(FIG / "Figure_8J.csv", index=False)
    return df


def register_manifest(files):
    """追加 RESULTS_MANIFEST_v2.0.csv。"""
    mf = GOV / "RESULTS_MANIFEST_v2.0.csv"
    man = pd.read_csv(mf, dtype={"s_number": str})
    for f, meta in files.items():
        if not f.exists():
            print("MISSING:", f)
            continue
        n = sum(1 for _ in open(f, encoding="utf-8")) - 1
        man.loc[len(man)] = {
            "s_number": meta["s_number"], "filename": f.name,
            "location": meta["location"], "class": meta["class"],
            "status": "canonical", "gene_universe": meta["gene_universe"],
            "gene_set_version": GENE_SET_VERSION, "score_version": meta["score_version"],
            "n_data_rows": n, "size_bytes": f.stat().st_size,
            "sha256": sha256(f), "analysis": meta["analysis"], "notes": meta["notes"],
        }
    man.to_csv(mf, index=False)
    print(f"RESULTS_MANIFEST_v2.0.csv updated: {len(man)} rows")


def register_figure_manifest():
    mf = FIG.parent / "FIGURE_DATA_MANIFEST.csv"
    man = pd.read_csv(mf)
    rows = [
        {"panel": "Figure_8I", "package_file": "Main/Figure_8I.csv",
         "source_file(01_RESULTS_TABLES)": "Table_S15e_pQTL_MR_Results.csv",
         "figure_title": "风险预测、因果线索与深度扰动（M3 增补）",
         "panel_content": "pQTL-MR（deCODE+UKB-PPP cis-pQTL 工具 × IEU 脓毒症/FinnGen ARDS）逐暴露-结局 IVW/Wald/加权中位数/Egger OR 与 BH（全阴：ITPR1 脓毒症 OR=1.142 p=0.406；TFRC ARDS UKB-PPP OR=1.207 p=0.188 为最小 p）",
         "suggested_plot": "森林图（暴露×结局，OR±95%CI）",
         "notes": "M3 (2026-08-18)；score_version=m3pqtl_v1.0；覆盖度 3/8 基因（ITPR1/SOD2/TFRC）"},
        {"panel": "Figure_8J", "package_file": "Main/Figure_8J.csv",
         "source_file(01_RESULTS_TABLES)": "Table_S15f_SMR_HEIDI_Results.csv",
         "figure_title": "风险预测、因果线索与深度扰动（M3 增补）",
         "panel_content": "SMR/HEIDI（eQTLGen 全血 72 基因 + GTEx v8 肺 31 基因 × FinnGen ARDS；103 检验全阴，min p=0.018）+ TWAS（GTEx v8 肺/全血权重 107 best-CV 检验；SLC7A11 全血 z=+5.13 BH=6.9e-5 但 top1 不印证、弥散模型驱动；MAP1LC3B 肺 z=-3.29、SLC40A1 全血 z=+3.22 名义未过 BH）",
         "suggested_plot": "双面板：SMR p 值排序点图（HEIDI 通过标记）+ TWAS 森林/点图（z±SE，BH 着色）",
         "notes": "M3 (2026-08-18)；score_version=m3smr_v1.0;m3twas_v1.0"},
    ]
    for r in rows:
        man.loc[len(man)] = r
    man.to_csv(mf, index=False)
    print(f"FIGURE_DATA_MANIFEST.csv updated: {len(man)} rows")


def write_judgment():
    mr = pd.read_csv(OUT / "Table_S15e_pQTL_MR_Results.csv")
    coloc = pd.read_csv(OUT / "Table_S15e_coloc_Results.csv")
    smr = pd.read_csv(OUT / "Table_S15f_SMR_HEIDI_Results.csv")
    tw = pd.read_csv(OUT / "Table_S15g_TWAS_Results.csv")
    best = tw[tw["method"] == "best"]
    top1 = tw[tw["method"] == "top1"]

    m = mr[mr["method"].isin(["Inverse variance weighted", "Wald ratio"])]
    min_p = m["pval"].min()
    min_row = m.loc[m["pval"].idxmin()]
    c1 = coloc[coloc["p12"] == 1e-5]
    r1 = (m["bh_p"] < 0.05).any() and (c1["PP.H4"] >= 0.8).any()
    r1b = ((smr["bh_p_pooled"] < 0.05) & (smr["heidi_pass"] == True)).any()
    judgment = "R3" if not (r1 or r1b) else "R1"
    slc = best[(best["gene"] == "ENSG00000151012") & (best["tissue"] == "GTExv8.EUR.Whole_Blood")]
    slc_top1 = top1[(top1["gene"] == "ENSG00000151012") & (top1["tissue"] == "GTExv8.EUR.Whole_Blood")]

    txt = """# M3 因果层翻盘——判定报告（预注册 R1–R3）

> 预注册：`M3_pre_registration_20260818.md`（判定规则先于结果）
> 执行日期：2026-08-18｜代码：`M3_step1..9` 系列（见 `M3_README_分析流程.md`）
> 数据源全部经逐源核实（见 `N4_资源核实报告.md` M3 小节）。

## 一、各腿结果（全部真实数值）

### 腿1 pQTL-MR（Table S15e；m3pqtl_v1.0）
- 覆盖度：8 目标基因中 3 个（ITPR1/SOD2/TFRC）在 deCODE（SomaScan）与 UKB-PPP（Olink）两套资源中均有 cis-pQTL 工具；VDAC1/CANX/SLC40A1/GPX4/FTH1 两源均无显著 cis-pQTL（不可仪器化，如实登记）。跨源 QC：同一变体（rs76604555/rs5746105）两源 BETA/SE/A1FREQ 完全一致。
- deCODE 层工具：ITPR1 1 个、SOD2 2 个（rs1551220/rs5746105，方向相反的双信号）、TFRC 6 个独立信号（rs2641379/rs546085354/rs114237070/rs28427561/rs56224756/rs113445087；rs28427561 F=1721 为主信号）；UKB-PPP 层 3 个单变体（敏感性）。
- IEU 脓毒症（ieu-b-69，主结局）：ITPR1 Wald OR=1.142（95%CI 0.835–1.562，p=0.406）、SOD2 IVW OR=1.028（p=0.753）、TFRC IVW OR=1.022（0.948–1.101，p=0.571）；全部 BH≥0.86。UKB-PPP 敏感性同向阴性（TFRC Wald OR=1.016，p=0.551）。
- FinnGen ARDS（次结局）：SOD2 IVW OR=0.701（p=0.442）、TFRC IVW OR=0.857（p=0.769）；UKB-PPP TFRC Wald OR=1.207（0.912–1.598，p=0.188，全表最小 p）。4 个工具变体（rs76604555/rs28427561/rs56224756/rs546085354）在 FinnGen R10 中缺席（rsid 列缺失或未通过 FinnGen QC，按 chr:pos 检索证实），ITPR1×ARDS 不可检，如实披露。
- 敏感性：多 SNP 层 Q 异质性全部不显著（Qp 0.25–0.58）；Egger/加权中位数（TFRC 层，n≥3）与 IVW 方向一致。
- 共定位（coloc.abf 5.2.3，lead±100kb，p12=1e-5 主口径）：6 组合 PP.H4 0.0072–0.0643，全部 <0.8；PP.H1 主导（0.86–0.95），即区域内仅有 pQTL 自身信号、无结局信号共享。敏感性 p12=1e-6/5e-5：PP.H4 最大 0.256（TFRC×FinnGen p12=5e-5），仍 <0.8。

### 腿2 SMR/HEIDI（Table S15f；m3smr_v1.0）
- eQTLGen 全血（n=31,470，Z→b 转换、AF 取 1000G EUR 近似）：72 基因可检；GTEx v8 肺（Portal API，NES 缩放口径）：31 基因可检（80 基因中仅 31 个为肺 eGene）。
- 103 个 gene×layer 检验全部 BH>0.70；最小 p=0.0176（CASP4，eQTLGen，HEIDI p=4.3e-45 未通过——连锁/多效异质性）；"SMR BH<0.05 且 HEIDI 通过"：0 个。
- 注：HEIDI 使用显著性 eQTL 对子集（非全 cis SNP），为文献通行口径，已在方法中披露。

### 腿3 TWAS（Table S15g；m3twas_v1.0）
- GTEx v8 肺+全血 FUSION 权重（best-CV 主口径，top1 单 SNP 敏感性）× FinnGen ARDS；107 个 best-CV 检验。
- 唯一 BH<0.05：SLC7A11 全血 susie（z=+5.13，p=2.9e-7，BH=6.9e-5）——但 (i) top1 敏感性 z=+0.47（p=0.64）完全不印证；(ii) 模型为 476 个 SNP 的弥散权重（|w|≤4.8e-4）；(iii) 该位点 FinnGen 最大 |z|=3.73（p=1.9e-4）为宽域弱信号；(iv) GTEx v8 全血该位点无任何显著 eQTL（coloc 不可行）。判为**模型依赖性信号，不作阳性**。
- 名义信号（未过 BH）：MAP1LC3B 肺 susie z=−3.29（p=9.9e-4，BH=0.10，top1 z=−1.72 p=0.085 同向弱印证）；SLC40A1 全血 enet z=+3.22（p=1.3e-3，BH=0.10，top1 z=+1.65 p=0.098 同向弱印证）。top1 敏感性全部 BH>0.9。

### 腿4 SeismicGWAS 空间链（Table S19i；m3seis_v1.0）
- FinnGen ARDS 96 个提示性位点（p<1e-5）内 847 个基因完成"位点→Visium 空间域（Banksy-lite k=5，逐切片域内 z+跨切片 Stouffer+域间 KW）→细胞类型（CosMx 27 类型均值，缺席回退 scRNA 8 类型）"链表；80 基因框架全域归因另存 `_intermediate/M3_seismic_80gene_spatial_attribution.csv`。链表为描述性定位，不产生新的因果主张。

## 二、判定（预注册规则）

- R1(i)：pQTL-MR IVW BH<0.05 且 PP.H4≥0.8 → **不成立**（全部 IVW BH≥0.77；PP.H4≤0.064 主口径）。
- R1(ii)：SMR BH<0.05 且 HEIDI 通过 → **不成立**（0/103）。
- TWAS 层：唯一 BH 显著结果为弥散模型驱动、top1 不印证、无 eQTL 共定位支持 → 不计阳性。
- **总判定：R3（阴性排除）**——四层工具变量证据（既有 eQTL-MR Table S15a/b + 本层 pQTL-MR/SMR/TWAS）均未支持任一候选节点与脓毒症/ARDS 的因果关联；SeismicGWAS 空间链为描述性定位。按预注册规则，以"四层工具变量均未支持因果"作为对概念边界的严谨排除写入 Discussion（阴性亦有发表价值，且排除审稿人质疑）。

## 三、与既有证据的关系
- 与 Table S15a（eQTL-MR：IP3R1 名义显著、未过多重校正）一致：pQTL 层同一方向无复制（ITPR1 pQTL Wald p=0.406）。
- 与 Table S15c（eQTL 共定位 PP.H1≈0.88–0.90 主导）一致：pQTL 共定位同样 PP.H1 主导、PP.H4 全低。
- 即：因果层翻盘未实现；M3 的贡献是**把"因果锚点"假设从 eQTL 单层排除升级为四层工具变量的严谨排除**。

## 四、产出清单
- 附表：Table_S15e_pQTL_MR_{Results,Sensitivity,Instruments,LeaveOneOut}.csv、Table_S15e_coloc_{Results,snp_pp}.csv、Table_S15f_SMR_HEIDI_{Results,GeneDetail}.csv、Table_S15g_TWAS_Results.csv、Table_S19i_SeismicGWAS_Spatial_CellType_Chain.csv
- 图面板数据：Figure_8I.csv（pQTL-MR 森林图数据）、Figure_8J.csv（SMR/TWAS 联合数据）
- 全部入 RESULTS_MANIFEST_v2.0（SHA256 冻结）。
"""
    (INT / "M3_judgment_report.md").write_text(txt, encoding="utf-8")
    print("judgment report written")
    print(f"判定: {judgment}（最小 pQTL-MR p={min_p:.4g}，{min_row['exposure']}×{min_row['outcome']}）")
    return judgment


def main():
    print("=" * 70)
    print("M3 step9: outputs & judgment")
    print("=" * 70)
    build_fig8i()
    build_fig8j()
    files = {
        OUT / "Table_S15e_pQTL_MR_Results.csv": {
            "s_number": "S15e", "location": "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
            "class": "supplementary_table", "gene_universe": "mitoxy_80",
            "score_version": "m3pqtl_v1.0", "analysis": "M3 pQTL-MR 主结果（deCODE+UKB-PPP × IEU脓毒症/FinnGen ARDS）",
            "notes": "M3 (2026-08-18); 全部 IVW/Wald BH>=0.77"},
        OUT / "Table_S15e_pQTL_MR_Sensitivity.csv": {
            "s_number": "S15e", "location": "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
            "class": "supplementary_table", "gene_universe": "mitoxy_80",
            "score_version": "m3pqtl_v1.0", "analysis": "M3 pQTL-MR 工具敏感性（F 统计等）",
            "notes": "M3 (2026-08-18)"},
        OUT / "Table_S15e_pQTL_MR_Instruments.csv": {
            "s_number": "S15e", "location": "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
            "class": "supplementary_table", "gene_universe": "mitoxy_80",
            "score_version": "m3pqtl_v1.0", "analysis": "M3 pQTL-MR harmonized 工具×结局表",
            "notes": "M3 (2026-08-18); 4 变体在 FinnGen R10 缺席已披露"},
        OUT / "Table_S15e_pQTL_MR_LeaveOneOut.csv": {
            "s_number": "S15e", "location": "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
            "class": "supplementary_table", "gene_universe": "mitoxy_80",
            "score_version": "m3pqtl_v1.0", "analysis": "M3 pQTL-MR leave-one-out",
            "notes": "M3 (2026-08-18)"},
        OUT / "Table_S15e_coloc_Results.csv": {
            "s_number": "S15e", "location": "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
            "class": "supplementary_table", "gene_universe": "mitoxy_80",
            "score_version": "m3pqtl_v1.0", "analysis": "M3 pQTL×结局 coloc.abf（p12=1e-5/1e-6/5e-5）",
            "notes": "M3 (2026-08-18); 主口径 PP.H4 全部 <0.065"},
        OUT / "Table_S15e_coloc_snp_pp.csv": {
            "s_number": "S15e", "location": "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
            "class": "supplementary_table", "gene_universe": "mitoxy_80",
            "score_version": "m3pqtl_v1.0", "analysis": "M3 pQTL×结局 coloc 逐 SNP PP.H4",
            "notes": "M3 (2026-08-18)"},
        OUT / "Table_S15f_SMR_HEIDI_Results.csv": {
            "s_number": "S15f", "location": "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
            "class": "supplementary_table", "gene_universe": "mitoxy_80",
            "score_version": "m3smr_v1.0", "analysis": "M3 SMR/HEIDI（eQTLGen 全血 + GTEx v8 肺 × FinnGen ARDS）",
            "notes": "M3 (2026-08-18); 103 检验全部 BH>0.70"},
        OUT / "Table_S15f_SMR_HEIDI_GeneDetail.csv": {
            "s_number": "S15f", "location": "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
            "class": "supplementary_table", "gene_universe": "mitoxy_80",
            "score_version": "m3smr_v1.0", "analysis": "M3 SMR/HEIDI 逐 SNP 明细",
            "notes": "M3 (2026-08-18)"},
        OUT / "Table_S15g_TWAS_Results.csv": {
            "s_number": "S15g", "location": "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
            "class": "supplementary_table", "gene_universe": "mitoxy_80",
            "score_version": "m3twas_v1.0", "analysis": "M3 TWAS（GTEx v8 肺+全血 FUSION 权重 × FinnGen ARDS；best-CV+top1）",
            "notes": "M3 (2026-08-18); SLC7A11 全血 BH<0.05 为弥散模型驱动、top1 不印证，不作阳性"},
        OUT / "Table_S19i_SeismicGWAS_Spatial_CellType_Chain.csv": {
            "s_number": "S19i", "location": "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
            "class": "supplementary_table", "gene_universe": "all_locus_genes",
            "score_version": "m3seis_v1.0", "analysis": "M3 SeismicGWAS 空间链：96 位点→Visium 域→细胞类型→基因",
            "notes": "M3 (2026-08-18); 描述性定位，不产生新因果主张"},
        FIG / "Figure_8I.csv": {
            "s_number": "", "location": "01_FIGURE_DATA_CSV/Main",
            "class": "figure_data", "gene_universe": "mitoxy_80",
            "score_version": "m3pqtl_v1.0", "analysis": "Figure 8I: pQTL-MR 森林图数据",
            "notes": "M3 (2026-08-18)"},
        FIG / "Figure_8J.csv": {
            "s_number": "", "location": "01_FIGURE_DATA_CSV/Main",
            "class": "figure_data", "gene_universe": "mitoxy_80",
            "score_version": "m3smr_v1.0;m3twas_v1.0", "analysis": "Figure 8J: SMR/HEIDI + TWAS 联合数据",
            "notes": "M3 (2026-08-18)"},
    }
    register_manifest(files)
    register_figure_manifest()
    write_judgment()
    print("\nDONE")


if __name__ == "__main__":
    main()
