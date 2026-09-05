# 图数据包说明（01_FIGURE_DATA_CSV）

> 生成日期：2026-08-16；2026-08-17 目录重组后路径更新；**2026-08-21 主稿重构（任务8）后图号体系更新**。本目录为**主图/附图可视化支撑数据包**：全部 CSV 复制自 `02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV/`（源文件受 `RESULTS_MANIFEST` SHA256 冻结与 `lint_package.py` 校验保护，**原件一律未改动/未改名**）。面板→源文件逐行映射见 `FIGURE_DATA_MANIFEST.csv`。

## 主图方案（6 个主图，39 个面板，每图 ≤16 面板）

| 图号 | 主题 | 面板 | 对应文稿主结果 |
|---|---|---|---|
| Figure 1 | 研究设计、队列流程与 80 基因通路框架 | 1A–1F | Results §1 |
| Figure 2 | 五维正交证据收敛于 MAM–铁节点（A–H）+ 呼吸链/MAM 跨平台塌陷与 GSE32707 外部验证（I–P） | 2A–2P | Results §2–3 |
| Figure 3 | 模块级解离（上游塌陷 vs 执行诱导）与复合评分方法依赖性 | 3A–3F | Results §4 |
| Figure 4 | MDI 临床化与跨疾病复制（判定 R2） | 4A–4E | Results §9 |
| Figure 5 | 生化层补证（判定阴性/部分支持：蛋白双口径一致率 33%/49%、铁轴签名、脂质/oxylipin、GEM、亚细胞定位） | 5A–5E | Results §10 |
| Figure 6 | 预注册统计单位重算与实验锚定验证（来源互斥 meta k=16、CosMx 患者级置换、IIAMD 签名稳定性、Gate 2 头对头） | 6A–6F | Results §11 |

注：Figure_2A/2B（同源 Bridge_Test_Result.csv 的模块/基因双视图）、Figure_2J/3A（同源 Table_S2 的基因/臂双视图）、Figure_2M/3F（同源 Table_S51 的模块/臂双视图）为同源文件的双面板复用，已各复制一份以保证"一面板一文件"。

## 附图方案（12 个附图，84 个面板，每图 ≤10 面板）

| 图号 | 主题 | 面板 |
|---|---|---|
| Figure S1 | 单细胞评分、QC 与细胞组成补充 | S1A–S1E |
| Figure S2 | WGCNA/hdWGCNA 共表达网络 | S2A–S2F |
| Figure S3 | 探索性因子分析（EFA） | S3A–S3D |
| Figure S4 | 染色质可及性、motif、SCENIC+、DTU | S4A–S4H |
| Figure S5 | MR 工具变量/敏感性、共定位逐 SNP、GWAS 位点 | S5A–S5I |
| Figure S6 | 细胞通讯与轨迹补充 | S6A–S6F |
| Figure S7 | 药物数据库、MD 状态、组合扰动、跨物种 | S7A–S7H |
| Figure S8 | 循环转录组与免疫浸润补充 | S8A–S8G |
| Figure S9 | 髓系细胞为主要效应载体（Scissor/MiloR/CIBERSORTx） | S9A–S9I |
| Figure S10 | 细胞通讯与髓系调控（CellPhoneDB/SCENIC/velocity） | S10A–S10G |
| Figure S11 | 可成药节点与药物重定位（scFOCAL/CMap/对接） | S11A–S11E |
| Figure S12 | 风险预测、MR/共定位、Geneformer、表观时钟、pQTL-MR、SMR/TWAS | S12A–S12J |

> 2026-08-21 任务8 图号体系变更记录：原主图 Figure 5–8（髓系载体/通讯/药物/风险因果）按主稿重构方案降入补充，编号顺延为 Figure S9–S12；原 Figure 9/10/11 前移为 Figure 4/5/6；原 Figure 3（跨平台塌陷+外部验证）并入新 Figure 2（面板 2I–2P）；原 Figure 4 顺延为 Figure 3。文件与 manifest 已同步重命名，FIGURE_DATA_MANIFEST.csv 面板映射已同步更新。

包内曾有历史图片 `Figure_S6_Dissociation_hexbin.png`（旧 S6 编号遗物，与本方案 Figure S6 编号含义不同；2026-08-16 已按 D4 归档流程移入归档区并在 RESULTS_MANIFEST 登记 archived；2026-09-04 起归档区统一为 `归档/`，原 archive/ 并入），绘图时勿混淆。

## 主表（置于正文）与附表方案

**正文主表（3 个）：**
- **Table 1**：Mitoxyperilysis 综合评分按细胞类型×临床状态 —— 数据源 `Main/Figure_1F.csv`
- **Table 2**：模块级五维收敛得分 —— 数据源 `Main/Figure_2A.csv`
- **Table 3**：GSE32707 外部验证模块级结果（A1 主对比 Δz/q 值 + 预注册判定 R4） —— 数据源 `Main/Figure_2M.csv`

**附表**：沿用包内现有 **Table S1–S60** 编号（文稿全文已按此引用，改号会破坏交叉引用），投稿时按期刊要求合并为单个 xlsx（每表一个 sheet）即可。*.txt/*.md 报告类文件作为 Supplementary Notes 一并提交。

## 未纳入图数据包的文件及原因（治理披露）

| 文件 | 原因 |
|---|---|
| Table_S22b_scFOCAL_CellType_IC50.csv | 源表自带 `score_version=legacy_scale_v0_not_for_text` 标记，按包内治理规则**不入正文/主图**；细胞类型维度已由 Figure_S11B 的 max_expression_celltype 列覆盖 |
| Table_S8_Dataset_DEG_Summary.csv | ⚠️ 旧标签遗物：将 GSE185263 误写为 "bulk, **lung**"（实为全血），建议回源修正后再用 |
| Table_S1_*_prev38.csv、Table_S2_DEGs_Analysis_pre_fix_backup.csv、Table_S18_*（非 CORRECTED 版）、Table_s12_*（小写重复） | 历史版本/备份/重复，已被 canonical 版取代 |
| Table_S5_PCD_Scores_RAW.csv | 未标准化原始尺度；Figure_3D 用 z 评分版 |
| GSE67530_beta_for_GrimAge.csv.gz（234 MB）、Table_S11_TF_Activity_Matrix_AllCells.csv（689 MB）、Table_S11_Mitoxyperilysis_TF_PerCell.csv（87.8 MB）、Table_S7k_Velocity_Pseudotime.csv（18.8 MB）、Table_S26d/S26i/S26k、Table_S26a/S26b_HLCA、Table_S26_EpiAgent、Table_S41a/c/d、Table_S4b/S4d、Table_S12_SCENIC_Regulons.csv、Table_S13 详细版、Table_S19_GWAS_Mitoxy_Integration.csv、Table_S20a/b/e/f、Table_S21c、Table_S6_Scissor_Marker_Genes.csv、Table_S8/S8_GSE212865/S8_Meta 全量版、Table_S29_scATAC_Peaks.csv、Table_S31_SCENICplus_GRN.csv、Table_S44_hdWGCNA_Modules.csv、Table_S46c/d 表达矩阵 等（2026-09-04：`WGCNA_intermediate_expr_block.csv` 已删除——R1 弃用模块中间产物，零现役引用；其 canonical 输出 Table S44/S44b/S45/S49 不受影响） | 原始矩阵/中间产物/全量明细，属数据存放层（建议存 GEO/figshare/Zenodo 并给 DOI），不适合作单面板绘图底表；其汇总层均已收入对应面板 |
| MR_bio_CRP_Sepsis.csv、MR_bio_Ferritin_Sepsis.csv | 生物标志物 MR 旁支分析，未进入文稿主线；如需可加 Figure_S5J–S5K |
| _pcd_score_summary.json | JSON 摘要（非 CSV），其内容已含于 Figure_3D 源表 |

## 绘图注意事项

1. 多数 CSV 首列带 UTF-8 BOM（`pandas.read_csv` 默认兼容；R 请用 `readr::read_csv` 或 `fileEncoding="UTF-8-BOM"`）。
2. Figure_3D 源表列名 `Cuprotosis` 为既有拼写（=Cuproptosis），绘图标签请改正。
3. Figure_S9A 的 Scissor 表型锚为"Sepsis_COVID-vs-对照"（非 ARDS 结局），图注须如实注明。
4. Figure_S12A 内部 AUC≈0.99 存在乐观偏倚（65/80 特征被选），图中应同时标注外部验证 AUC=0.705。
5. Figure_S11D 对接最优靶 SLC40A1 为 AlphaFold 预测结构、35 对 MD 全部 Pending，图注须声明。
6. 所有版本列（`gene_set_version`/`score_version`）请保留在图数据脚注，投稿审查时用于溯源。
