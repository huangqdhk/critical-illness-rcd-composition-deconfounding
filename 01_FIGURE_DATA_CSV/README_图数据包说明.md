# 图数据包说明（01_FIGURE_DATA_CSV）

> **当前状态（本说明以正文图号为准）**：主图 **8 张（62 面板）** + 附图 **16 张（106 面板）**，合计 **168 面板**；全部 170 mm 宽，页高 20/21 ≤ 247 mm（唯一例外 Figure_S3 = 268 mm，见正文说明）。正文引用严格升序。
>
> **本目录是主图/附图可视化支撑数据包**：除湿实验三图（正文 Figure 8 / S15 / S16，见下文 2026-10-02 沿革）外，全部面板 CSV 复制自 `02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV/`（源文件受 `04_AUDIT_GOVERNANCE/RESULTS_MANIFEST_v2.0.csv` 的 SHA256 冻结与根目录 `lint_package.py` 校验保护，**原件一律未改动、未改名**）；湿实验三图的面板 CSV 取自图数据包外的湿实验源包（`Figure_8\Figure_8\`、`Figure_S15-16待可视化数据\`），并入时同样逐字节复制、数值零改动。
>
> **⚠️ 包内文件名不等于正文图号**：面板 CSV 的文件名沿用 **pkg 世代**编号（如 `Figure_5A.csv`、`Figure_8A.csv`、`Figure_S9A.csv`），正文用的是当前 **ms 世代**编号（Figure 4A、Figure 6A、Figure S11A）。**例外**：正文湿实验 Figure 8 的面板为 `Main/Figure_8wA–8wH.csv`（`w` = wet-lab，入包时命名，用以与包内既有同名不同内容的 `Figure_8A–8G.csv`〔正文 Figure 6 世代〕区隔）。逐行映射见随包 **`FIGURE_CROSSWALK.csv`**（列：`ms_figure` / `ms_panel` / `pkg_panel` / `csv_relpath` / `role` / `ms_title` / `suggested_plot` / `panel_content` / `pkg_source_file`），**这是本包唯一的"正文图号 ↔ 文件"桥接表**；历史（11→9 图）映射见 `FIGURE_RENUMBER_MAP_20260905.csv`。面板→源表映射见 `FIGURE_DATA_MANIFEST.csv`。

## 主图方案（8 张，62 面板，按正文引用顺序编号）

| 图号 | 主题 | 面板 | 对应文稿 |
|---|---|---|---|
| Figure 1 | 研究设计、队列流程与 80 基因通路框架 | 1A–1F（6） | §1 |
| Figure 2 | 模块级解离、跨平台复现与 GSE32707 外部验证 | 2A–2H（8） | §2–3 |
| Figure 3 | 五条正交证据线收敛于 MAM–铁节点 | 3A–3H（8） | §4 |
| Figure 4 | 解离状态的临床与生化边界：MDI 临床转化与跨疾病复制 | 4A–4J（10） | §6–7 |
| Figure 5 | 统计单位重算、组成审计与扰动三角验证 | 5A–5I（9） | §10–11 |
| Figure 6 | 时序先后（A–D）与干预响应（E–G） | 6A–6G（7） | §12–13 |
| Figure 7 | 泛化性：已发布工具 mitoxdi（A–C）与力学边界（D–F） | 7A–7F（6） | §14–15 |
| Figure 8 | 湿实验验证：两臂解离与 Deferasirox 干预（小鼠肺损伤 W1–W3 与 THP-1） | 8A–8H（8） | §16；Methods 4.29 |

注：正文 Figure 8 为湿实验图，8 个面板的源文件为 `Main/Figure_8wA–8wH.csv`（2026-10-02 随包；此前以外部路径登记，见 `WETLAB_FIG8_INTEGRATION_STATUS_20260928.md`）。

注：正文 Figure 2 的 8 个面板取自包内两个世代的文件——2A–2D ← `Main/Figure_2A–2D.csv`，2E/2F/2G ← `Main/Figure_3A/3B/3C.csv`，2H ← `Main/Figure_3E.csv`（09-16 主图压缩把旧 Figure 2 与旧 Figure 3 合并为新 Figure 2，文件未改名）。

## 附图方案（16 张，106 面板，按正文首次引用顺序编号）

| 图号 | 主题 | 面板 | 首次引用 |
|---|---|---|---|
| Figure S1 | 单细胞评分、QC 与细胞组成补充 | S1A–S1E（5） | §2 |
| Figure S2 | 方法依赖性与外部验证补充 | S2A–S2F（6） | §2–3 |
| Figure S3 | 循环转录组与免疫浸润补充 | S3A–S3G（7） | §2–3 |
| Figure S4 | WGCNA/hdWGCNA 共表达网络与探索性因子分析 | S4A–S4J（10） | §4 |
| Figure S5 | 因果层补充 | S5A–S5E（5） | §8 |
| Figure S6 | MR 工具变量、敏感性、逐 SNP 共定位与 GWAS 位点 | S6A–S6I（9） | §8 |
| Figure S7 | 髓系细胞为主要效应载体 | S7A–S7I（9） | §9 |
| Figure S8 | 细胞通讯与髓系调控 | S8A–S8G（7） | §9 |
| Figure S9 | 通讯与轨迹补充 | S9A–S9F（6） | §9 |
| Figure S10 | 统计单位与区室补充 | S10A–S10E（5） | §10–11 |
| Figure S11 | 可成药节点、药物库与跨物种补充 | S11A–S11G（7） | §13 |
| Figure S12 | 可成药节点与药物重定位 | S12A–S12E（5） | §13 |
| Figure S13 | 风险预测、深度扰动与表观遗传时钟 | S13A–S13E（5） | §12–13 |
| Figure S14 | scATAC、motif、SCENIC+ 与转录本使用补充 | S14A–S14H（8） | §12 |
| Figure S15 | W1 小鼠损伤/铁轴/解离读数补充（湿实验） | S15A–S15E（5） | §16 |
| Figure S16 | 髓系与上皮细胞实验 C1–C3（湿实验） | S16A–S16F（6） | §16、§19 |

> **面板映射说明（2026-10-02 更新）**：`FIGURE_CROSSWALK.csv` 已按出图管线台账（`可视化/00_对照表/FIGURE_CROSSWALK.csv`）重建——正文 Figure S10 的 5 个面板均有一行映射（其中 **S10D 与 S10E 同源于 `Main/Figure_7E.csv`**，为 Layer 1–2 与 Layer 3 两个视图）；正文 Figure S15 的 5 个面板中第 5 个（**S15E**）的源文件为 `Supplementary/Figure_S15F.csv`（2026-09-30 二次改版撤出两块面板后沿用；另两个历史面板源 `Figure_S15E.csv`、`Figure_S15G.csv` 一并随包保留，当前图未使用）。

## 图号世代沿革（历史留痕）

> **2026-10-02 湿实验三图随包（本次）**：正文 Figure 8（8 面板）、Figure S15（5 面板）、Figure S16（6 面板）的源数据并入本包——新增 `Main/Figure_8wA–8wH.csv`（图 8 专用，`w` = wet-lab）、`Supplementary/Figure_S15A–S15G.csv`、`Supplementary/Figure_S16A–S16F.csv`，共 21 个文件（逐字节复制、数值零改动）。`FIGURE_CROSSWALK.csv` 新增 19 行并同步至出图管线口径（含 S10/S11 的既有修正：S10 补 S10E 行、S11 由 8 行订正为 7 行）；`FIGURE_DATA_MANIFEST.csv` 同步重建（168 行）。湿实验源包（`Figure_8\`、`Figure_S15-16待可视化数据\`）仍为唯一原始来源；S15 历史面板源 `Figure_S15E/S15G.csv` 一并保留（当前图未用）。接入前状态见 `WETLAB_FIG8_INTEGRATION_STATUS_20260928.md`。

> **2026-09-16 主图压缩（9→7，迎接湿实验主图 Figure 8）**：主图 9→7、附图 12→14。旧 Figure 2 + 旧 Figure 3 → 新 Figure 2；旧 Figure 6 + 旧 Figure 7 → 新 Figure 5；旧 Figure 4→3、5→4、8→6、9→7 顺移。降级面板：旧 2E/2F/3D/3F/3G/3H → 新 S2A–S2F，旧 6B/6C/7G/7E → 新 S10A–S10D。附图级联：S2→S3、S3→S4、S4→S5、S5→S6、S6→S7、S7→S8、S8→S9、S9→S11、S10→S12、S11→S13、S12→S14。**本轮只改正文图号，数据文件名一概未动**（因此包内文件名与正文图号自本轮起分属两套编号，映射由 `FIGURE_CROSSWALK.csv` 承担）。

> **2026-09-05 CDD 投稿适配图号重排（11→9 主图）**：为匹配 Cell Death & Differentiation 惯例（主图每图 ≤10 面板且各图面板数均衡、全文引用严格升序）重排主图编号；附图同日下午按正文首次引用位置二次重编号（79 个 CSV 迁移，数值零改动）。逐行记录见 `FIGURE_RENUMBER_MAP_20260905.csv`。
>
> ⚠️ **历史脚本再生旧名警示**：根目录与 `05_RELEASE_GITHUB/repo/analysis_scripts/` 中的历史导出/登记脚本（`M2_step4_outputs_judgment.py`、`M4_finalize_outputs.py`、`M11M12_step4_export.py`、`M15_step4_export_and_report.py`、`M16_step5_export_register.py`、P13 系列等）仍按**更早世代的图号**读写 `Figure_*.csv`。按既有先例保留脚本原样（执行溯源优先，不追溯改写已运行代码）；**如需重跑，输出文件必须按 `FIGURE_CROSSWALK.csv`（当前）与 `FIGURE_RENUMBER_MAP_20260905.csv`（历史）核对重命名后方可入包/登记**，否则会产生旧名幽灵文件并污染 manifest。

> 2026-08-21 任务8 图号体系变更记录：原主图 Figure 5–8（髓系载体/通讯/药物/风险因果）按主稿重构方案降入补充；原 Figure 9/10/11 前移为 Figure 4/5/6；原 Figure 3 并入新 Figure 2（面板 2I–2P）；原 Figure 4 顺延为 Figure 3。

> 包内曾有历史图片 `Figure_S6_Dissociation_hexbin.png`（旧 S6 编号遗物，与本方案 Figure S6 编号含义不同；已按 D4 归档流程移入 `归档/` 并在 RESULTS_MANIFEST 登记 archived），绘图时勿混淆。

## 主表（置于正文）与附表方案

**正文主表（4 张）：**
- **Table 1**：80 基因框架综合评分按细胞类型×临床状态 —— 数据源 `Main/Figure_1F.csv`（正文 1F）
- **Table 2**：GSE32707 外部验证模块级结果（A1 主对比 Δz/q 值） —— 数据源 `Main/Figure_3E.csv`（正文 2H）
- **Table 3**：模块级五维收敛得分 —— 数据源 `Main/Figure_4A.csv`（正文 3A）
- **Table 4**：证伪级联预注册判定台账 —— 数据源 `Table4_criteria_verdicts.csv`

**附表**：沿用包内 **Table S1–S93** 编号（文稿全文已按此引用，改号会破坏交叉引用），投稿时按期刊要求合并为单个 xlsx（每表一个 sheet）。*.txt/*.md 报告类文件作为 Supplementary Notes 一并提交。

**图源数据（Source Data）**：按 Nature Portfolio（含 CDD）惯例，每图合并为单个 xlsx（每面板一个 sheet），见 `投稿CELL DEATH AND DIFFERENTIATION/Source_Data/`（**22 个**：主图 7 + Table 4 + 附图 14；由 `N5b_build_source_data_xlsx_crosswalk.py` 以 `FIGURE_CROSSWALK.csv` 为唯一真源、对本目录 CSV 逐字节生成）。每个 xlsx 的 `Index` sheet 逐面板列出 `package_panel` 与 `csv_relpath`，可据此回溯；面板的 `gene_set_version`/`score_version` 亦在该页登记。（2026-10-02 注：`投稿CELL DEATH AND DIFFERENTIATION/Source_Data/` 的 xlsx 为 2026-09 的投稿树交付形态；当前交付以本目录逐面板 CSV + 发布仓库〔GitHub/Zenodo，code DOI 10.5281/zenodo.22655285 的 concept DOI〕为准。）

## 未纳入图数据包的文件及原因（治理披露）

| 文件 | 原因 |
|---|---|
| Table_S22b_scFOCAL_CellType_IC50.csv | 源表自带 `score_version=legacy_scale_v0_not_for_text` 标记，按包内治理规则**不入正文/主图**；细胞类型维度已由 `Figure_S11B.csv`（正文 S13B）的 max_expression_celltype 列覆盖 |
| Table_S8_Dataset_DEG_Summary.csv | ⚠️ 旧标签遗物：将 GSE185263 误写为 "bulk, **lung**"（实为全血），建议回源修正后再用 |
| Table_S1_*_prev38.csv、Table_S2_DEGs_Analysis_pre_fix_backup.csv、Table_S18_*（非 CORRECTED 版）、Table_s12_*（小写重复） | 历史版本/备份/重复，已被 canonical 版取代 |
| Table_S5_PCD_Scores_RAW.csv | 未标准化原始尺度；`Figure_2D.csv` 用 z 评分版 |
| GSE67530_beta_for_GrimAge.csv.gz（234 MB）、Table_S11_TF_Activity_Matrix_AllCells.csv（689 MB）、Table_S11_Mitoxyperilysis_TF_PerCell.csv（87.8 MB）、Table_S7k_Velocity_Pseudotime.csv（18.8 MB）、Table_S26d/S26i/S26k、Table_S26a/S26b_HLCA、Table_S26_EpiAgent、Table_S41a/c/d、Table_S4b/S4d、Table_S12_SCENIC_Regulons.csv、Table_S13 详细版、Table_S19_GWAS_Mitoxy_Integration.csv、Table_S20a/b/e/f、Table_S21c、Table_S6_Scissor_Marker_Genes.csv、Table_S8/S8_GSE212865/S8_Meta 全量版、Table_S29_scATAC_Peaks.csv、Table_S31_SCENICplus_GRN.csv、Table_S44_hdWGCNA_Modules.csv、Table_S46c/d 表达矩阵、WGCNA_intermediate_expr_block.csv 等 | 原始矩阵/中间产物/全量明细，属数据存放层（建议存 GEO/figshare/Zenodo 并给 DOI），不适合作单面板绘图底表；其汇总层均已收入对应面板 |
| MR_bio_CRP_Sepsis.csv、MR_bio_Ferritin_Sepsis.csv | 生物标志物 MR 旁支分析，未进入文稿主线；如需可加为正文 Figure S6 的增补面板（包内 `Figure_S5J/S5K.csv`） |
| _pcd_score_summary.json | JSON 摘要（非 CSV），其内容已含于 `Figure_2D.csv` |

**包内非面板文件**：`FIGURE_CROSSWALK.csv`（正文图号 ↔ 文件映射，见文首）、`FIGURE_RENUMBER_MAP_20260905.csv`（历史映射）、`FIGURE_DATA_MANIFEST.csv`（面板→源表映射）、`Table4_criteria_verdicts.csv`（正文 Table 4 数据源）、`Figure_S2E_fit_stats.csv`（**非面板派生表**：原混在 `Figure_S2E.csv` 表尾的 3 行聚合统计，量纲与逐细胞类型占比不同，已拆为独立小表、数值零改动；正文面板 S2E 不使用这些行）、`WETLAB_FIG8_INTEGRATION_STATUS_20260928.md`（2026-09-28 湿实验接入台账；接入已于 2026-10-02 完成，保留为审计留痕，见沿革）。

## 绘图注意事项

1. 多数 CSV 首列带 UTF-8 BOM（`pandas.read_csv` 默认兼容；R 请用 `readr::read_csv` 或 `fileEncoding="UTF-8-BOM"`）。
2. `Figure_2D.csv`（正文 2D）源表列名 `Cuprotosis` 为既有拼写（=Cuproptosis），绘图标签请改正；该拼写按已登记的冻结溯源决定保留，并在 Source Data 的 `Index` 页加注说明。
3. `Figure_S6A.csv`（正文 S7A）的 Scissor 表型锚为"Sepsis_COVID-vs-对照"（非 ARDS 结局），图注须如实注明。
4. `Figure_S11A.csv`（正文 S13A）内部 AUC≈0.99 存在乐观偏倚（65/80 特征被选），图中应同时标注外部验证 AUC=0.705。
5. `Figure_S10D.csv`（正文 S12D）对接最优靶 SLC40A1 为 AlphaFold 预测结构、35 对 MD 全部 Pending，图注须声明。
6. 面板级 N 与误差棒定义已逐图写入正文图注（CDD 对 N 与误差棒描述均为 "must"），改图时勿删。
7. 所有版本列（`gene_set_version`/`score_version`）请保留在图数据脚注，投稿审查时用于溯源；M10–M16 模块的图数据按既有导出设计不含该两列（manifest 记 NA）。
