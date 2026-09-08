# 图数据包说明（01_FIGURE_DATA_CSV）

> 生成日期：2026-08-16；2026-08-17 目录重组后路径更新；2026-08-21 主稿重构（任务8）图号体系更新；**2026-09-05 CDD 投稿适配图号体系重排（本次）**。本目录为**主图/附图可视化支撑数据包**：全部 CSV 复制自 `02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV/`（源文件受 `RESULTS_MANIFEST` SHA256 冻结与 `lint_package.py` 校验保护，**原件一律未改动/未改名**）。面板→源文件逐行映射见 `FIGURE_DATA_MANIFEST.csv`；新旧图号逐行映射见 `FIGURE_RENUMBER_MAP_20260905.csv`。

## 主图方案（9 个主图，64 个面板，每图 6–10 面板，按正文引用顺序编号）

| 图号 | 主题 | 面板 | 对应文稿主结果 |
|---|---|---|---|
| Figure 1 | 研究设计、队列流程与 80 基因通路框架 | 1A–1F（6） | Results §1 |
| Figure 2 | 模块级解离（上游塌陷 vs 执行诱导）与复合评分方法依赖性 | 2A–2F（6） | Results §2 |
| Figure 3 | 呼吸链/MAM 塌陷跨平台复现与 GSE32707 外部验证 | 3A–3H（8） | Results §3 |
| Figure 4 | 五维正交证据收敛于 MAM–铁节点 | 4A–4H（8） | Results §4 |
| Figure 5 | 临床与生化边界：MDI 临床化与跨疾病复制（A–E）+ 生化层补证（F–J） | 5A–5J（10） | Results §6–7 |
| Figure 6 | 预注册统计单位重算与实验锚定验证（来源互斥 meta k=16、CosMx 患者级置换、IIAMD 签名稳定性、Gate 2 头对头） | 6A–6F（6） | Results §10 |
| Figure 7 | 组成枢纽检验、去混淆审计与扰动三角（组成是主要混杂，上游塌陷含髓系内在成分） | 7A–7G（7） | Results §11 |
| Figure 8 | 时序先后（A–D）与干预响应（E–G） | 8A–8G（7） | Results §12–13 |
| Figure 9 | 开放工具 mitoxdi（A–C）与力学边界（D–F） | 9A–9F（6） | Results §14–15 |

注：Figure_4A/4B（同源 Bridge_Test_Result.csv 的模块/基因双视图）、Figure_3B/2A（同源 Table_S2 的基因/臂双视图）、Figure_3E/2F（同源 Table_S51 的模块/臂双视图）为同源文件的双面板复用，已各复制一份以保证"一面板一文件"。

## 附图方案（12 个附图，84 个面板，每图 5–10 面板，**按正文首次引用顺序编号**）

| 图号 | 主题 | 面板 |
|---|---|---|
| Figure S1 | 单细胞评分、QC 与细胞组成补充（§1 首引） | S1A–S1E（5） |
| Figure S2 | 循环转录组与免疫浸润补充（§3 首引） | S2A–S2G（7） |
| Figure S3 | WGCNA/hdWGCNA 共表达网络与 EFA（§3 首引） | S3A–S3J（10） |
| Figure S4 | 因果层补充：MR 森林/共定位/pQTL-MR/SMR+HEIDI/TWAS（§8 首引） | S4A–S4E（5） |
| Figure S5 | MR 工具变量/敏感性、共定位逐 SNP、GWAS 位点（§8 首引） | S5A–S5I（9） |
| Figure S6 | 髓系细胞为主要效应载体（Scissor/MiloR/CIBERSORTx，§9 首引） | S6A–S6I（9） |
| Figure S7 | 细胞通讯与髓系调控（CellPhoneDB/SCENIC/velocity，§9 首引） | S7A–S7G（7） |
| Figure S8 | 细胞通讯与轨迹补充（§9 首引） | S8A–S8F（6） |
| Figure S9 | 药物数据库、MD 状态、组合扰动、跨物种（§11.7 首引） | S9A–S9H（8） |
| Figure S10 | 可成药节点与药物重定位（§16 首引） | S10A–S10E（5） |
| Figure S11 | 风险预测、深度扰动与表观时钟（§16 首引） | S11A–S11E（5） |
| Figure S12 | 染色质可及性、motif、SCENIC+、DTU（§16 首引） | S12A–S12H（8） |

> **2026-09-05（下午）二次重排记录**：上午首排补入 S1–S7 正文首引后，附图首引顺序出现错位（循环层@§3 早于 WGCNA、因果层先于 MR 细节、scATAC@§16 居末），不符合 Nature Portfolio"附图按引用顺序编号"惯例；遂按**正文首次引用位置**对附图二次重编号（79 个 CSV 迁移、数值零改动）：S7→S2、S2→S3、S8→S4、S4→S5、S9→S6、S10→S7、S5→S8、S6→S9、S11→S10、S12→S11、S3→S12（S1 不动）。主图编号不受影响。台账见 `FIGURE_RENUMBER_MAP_20260905.csv`（二次行），manifest/lint/文稿/图数据 manifest/Excel 源数据已全部同步。

> **2026-09-05 CDD 投稿适配图号重排记录**（匹配 Cell Death & Differentiation 惯例：主图每图 ≤10 面板且各图面板数均衡、全文引用严格升序）：
> **主图 11→9**：旧 Figure 3 → 新 Figure 2；旧 Figure 2I–2P → 新 Figure 3；旧 Figure 2A–2H → 新 Figure 4；旧 Figure 4 + 旧 Figure 5 → 新 Figure 5（Part II 临床+生化边界合并）；旧 Figure 6/7 不动；旧 Figure 8 + 旧 Figure 9 → 新 Figure 8（时序+干预）；旧 Figure 10 + 旧 Figure 11 → 新 Figure 9（工具+力学边界）。消除了旧 Figure 2 的 16 面板超载、旧 Figure 9/10/11 各仅 3 面板的过薄、以及"§2 引图 3、§3 引图 2I–2P、§4 引图 2A–2H"的乱序。
> **附图 12→12**：旧 S3（EFA，4 面板）并入新 S2（WGCNA+EFA 同属无监督结构发现）；旧 S4–S8 顺移为 S3–S7；旧 S12 三主题拆分——因果层 5 面板（C/D/E/I/J）独立为新 S8（§8 因果边界首引），风险预测/深度学习/表观时钟 5 面板（A/B/F/G/H）保留为新 S12；S9–S11 不动。全文附图引用由此严格升序（§8→S8、§9→S9/S10、§16→S11/S12）。
> 文件已两阶段安全重命名（91 个迁移），RESULTS_MANIFEST_v2.0.csv 文件名同步（sha256 不变、数值零改动），lint_package.py 同步，**lint 74/0 全绿**；FIGURE_DATA_MANIFEST.csv 已重建并补登 M15/M16 六面板（原 10A–10C/11A–11C → 9A–9F）。
>
> ⚠️ **历史脚本再生旧名警示**：根目录与 `05_RELEASE_GITHUB/repo/analysis_scripts/` 中的历史导出/登记脚本（M2_step4_outputs_judgment.py、M4_finalize_outputs.py、M11M12_step4_export.py、M15_step4_export_and_report.py、M16_step5_export_register.py、P13 系列等）仍按**旧图号**读写 `Figure_*.csv`（如 Figure_9A=干预、Figure_10A=工具、Figure_11A=力学）。按任务8 先例保留脚本原样（执行溯源优先，不追溯改写已运行代码）；**如需重跑，输出文件必须按 `FIGURE_RENUMBER_MAP_20260905.csv` 重命名后方可入包/登记**，否则将同时产生旧名幽灵文件并污染 manifest。

> 2026-08-21 任务8 图号体系变更记录：原主图 Figure 5–8（髓系载体/通讯/药物/风险因果）按主稿重构方案降入补充，编号顺延为 Figure S9–S12；原 Figure 9/10/11 前移为 Figure 4/5/6；原 Figure 3（跨平台塌陷+外部验证）并入新 Figure 2（面板 2I–2P）；原 Figure 4 顺延为 Figure 3。文件与 manifest 已同步重命名，FIGURE_DATA_MANIFEST.csv 面板映射已同步更新。

包内曾有历史图片 `Figure_S6_Dissociation_hexbin.png`（旧 S6 编号遗物，与本方案 Figure S6 编号含义不同；2026-08-16 已按 D4 归档流程移入 `archive/` 并在 RESULTS_MANIFEST 登记 archived），绘图时勿混淆。

## 主表（置于正文）与附表方案

**正文主表（4 个）：**
- **Table 1**：80 基因框架综合评分按细胞类型×临床状态 —— 数据源 `Main/Figure_1F.csv`
- **Table 2**：GSE32707 外部验证模块级结果（A1 主对比 Δz/q 值） —— 数据源 `Main/Figure_3E.csv`
- **Table 3**：模块级五维收敛得分 —— 数据源 `Main/Figure_4A.csv`
- **Table 4**：证伪级联预注册判定台账 —— 数据源 `Table4_criteria_verdicts.csv`

**附表**：沿用包内现有 **Table S1–S93** 编号（文稿全文已按此引用，改号会破坏交叉引用），投稿时按期刊要求合并为单个 xlsx（每表一个 sheet）即可。*.txt/*.md 报告类文件作为 Supplementary Notes 一并提交。

**图源数据（Source Data）**：按 Nature Portfolio（含 CDD）惯例，每图合并为单个 xlsx（每面板一个 sheet），见 `投稿CELL DEATH AND DIFFERENTIATION/Source_Data/Main|Supplementary/`（2026-09-05 由本目录 CSV 逐字节生成）。

## 未纳入图数据包的文件及原因（治理披露）

| 文件 | 原因 |
|---|---|
| Table_S22b_scFOCAL_CellType_IC50.csv | 源表自带 `score_version=legacy_scale_v0_not_for_text` 标记，按包内治理规则**不入正文/主图**；细胞类型维度已由 Figure_S11B 的 max_expression_celltype 列覆盖 |
| Table_S8_Dataset_DEG_Summary.csv | ⚠️ 旧标签遗物：将 GSE185263 误写为 "bulk, **lung**"（实为全血），建议回源修正后再用 |
| Table_S1_*_prev38.csv、Table_S2_DEGs_Analysis_pre_fix_backup.csv、Table_S18_*（非 CORRECTED 版）、Table_s12_*（小写重复） | 历史版本/备份/重复，已被 canonical 版取代 |
| Table_S5_PCD_Scores_RAW.csv | 未标准化原始尺度；Figure_2D 用 z 评分版 |
| GSE67530_beta_for_GrimAge.csv.gz（234 MB）、Table_S11_TF_Activity_Matrix_AllCells.csv（689 MB）、Table_S11_Mitoxyperilysis_TF_PerCell.csv（87.8 MB）、Table_S7k_Velocity_Pseudotime.csv（18.8 MB）、Table_S26d/S26i/S26k、Table_S26a/S26b_HLCA、Table_S26_EpiAgent、Table_S41a/c/d、Table_S4b/S4d、Table_S12_SCENIC_Regulons.csv、Table_S13 详细版、Table_S19_GWAS_Mitoxy_Integration.csv、Table_S20a/b/e/f、Table_S21c、Table_S6_Scissor_Marker_Genes.csv、Table_S8/S8_GSE212865/S8_Meta 全量版、Table_S29_scATAC_Peaks.csv、Table_S31_SCENICplus_GRN.csv、Table_S44_hdWGCNA_Modules.csv、Table_S46c/d 表达矩阵、WGCNA_intermediate_expr_block.csv 等 | 原始矩阵/中间产物/全量明细，属数据存放层（建议存 GEO/figshare/Zenodo 并给 DOI），不适合作单面板绘图底表；其汇总层均已收入对应面板 |
| MR_bio_CRP_Sepsis.csv、MR_bio_Ferritin_Sepsis.csv | 生物标志物 MR 旁支分析，未进入文稿主线；如需可加 Figure_S5J–S5K |
| _pcd_score_summary.json | JSON 摘要（非 CSV），其内容已含于 Figure_2D 源表 |

## 绘图注意事项

1. 多数 CSV 首列带 UTF-8 BOM（`pandas.read_csv` 默认兼容；R 请用 `readr::read_csv` 或 `fileEncoding="UTF-8-BOM"`）。
2. Figure_2D 源表列名 `Cuprotosis` 为既有拼写（=Cuproptosis），绘图标签请改正。
3. Figure_S6A 的 Scissor 表型锚为"Sepsis_COVID-vs-对照"（非 ARDS 结局），图注须如实注明。
4. Figure_S11A 内部 AUC≈0.99 存在乐观偏倚（65/80 特征被选），图中应同时标注外部验证 AUC=0.705。
5. Figure_S10D 对接最优靶 SLC40A1 为 AlphaFold 预测结构、35 对 MD 全部 Pending，图注须声明。
6. 所有版本列（`gene_set_version`/`score_version`）请保留在图数据脚注，投稿审查时用于溯源。
