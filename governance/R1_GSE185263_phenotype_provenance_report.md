# R1 · GSE185263 表型 provenance 核查报告

**核查日期**：2026-08-15
**核查人**：自动核查管线（对应对外质控报告 R1 条目）
**结论**：GSE185263 为**全血 RNA-seq 脓毒症队列**（Baghela & Hancock, UBC）；既往 82 例 "ARDS" 标签**无任何元数据支持**，实为 `sepcv*` 前缀样本 = **脓毒症合并 COVID-19（Sepsis_COVID）**。已全包更正。

---

## 一、五重独立证据链

### 证据 1 · GEO 官方记录（NCGI accession GSE185263，2026-08-15 直连复核）

| 字段 | 官方内容 |
|---|---|
| Title | **Predicting sepsis severity at first clinical presentation: the role of endotypes and mechanistic signatures** |
| Organism / Platform | Homo sapiens / Illumina HiSeq 2500 (GPL16791) |
| Summary | "…we characterized the **blood immune profiles** of patients with early/pre-sepsis to identify signatures reflecting disease severity, organ dysfunction, mortality, and specific endotypes/mechanisms" |
| **Overall design** | "**Whole blood RNA-seq** and clinical data was collected from **348 patients** from four emergency rooms and one intensive care unit" |
| Contributors | Baghela A, Hancock RE（University of British Columbia） |
| Samples | 392（GSM5608946–GSM5609337），样本名与本地矩阵逐一对应 |

→ **官方记录通篇无 "ARDS" 字样；组织=全血；患者=脓毒症（急诊+ICU）。**
→ "348 患者 + 44 健康对照" 与本包方案 line 18/201 的拆分一致，与 groups.csv 的 "266+82+44" 拆分为**本包内部再分层**。

### 证据 2 · 原始矩阵样本前缀（`00_RAW_DATA/GSE185263_Lung_ARDS/GSE185263_raw_counts.csv`，392 样本）

| 前缀 | n | 含义（按本包管线一致映射） |
|---|---|---|
| sepcol | 67 | 脓毒症（Sepsis 组） |
| sepnet | 104 | 脓毒症（Sepsis 组） |
| sepwes | 84 | 脓毒症（Sepsis 组） |
| sepvh | 11 | 脓毒症（Sepsis 组） |
| **sepcv** | **82** | **脓毒症合并 COVID-19（'cv'=COVID），全血**；样本名后缀 T0/T1/T2 为采血时点（sepcv702W1/703W1 为周采血） |
| hccol / hchlD / hcwes / hcwimr | 6+9+24+5 = 44 | 健康对照（Control 组） |

合计：Sepsis 266 + Sepsis_COVID 82 + Control 44 = 392 ✓（与 groups.csv 数量一致）

### 证据 3 · 本包内早已存在的正确映射（"ARDS" 是后期覆盖）

- `CIBERSORTx_Immune_Infiltration_14.1.py` L522–534（classify 函数）：`sepcv → ('Sepsis_COVID', 'Blood')`，sepcol→Colon、sepnet→Blood_NETs、hchlD→Lung 等；
- `Table_S16e_Sample_Condition_Info.csv`（392 行）：Condition = {Sepsis 266, Healthy_Control 44, **Sepsis_COVID 82**}；Tissue = {Colon 73, Blood_NETs 104, Other 124, Lung 9, **Blood 82**}；
- `35.1_circulating_transcriptome_analysis.py` L156 注释："Blood samples: sepnet (Blood_NETs, n=104), **sepcv (Blood, n=82)**"，其输出报告同样写 "Sepsis_COVID: 82 (Blood)"；
- `33.1_transcript_analysis.py` L108：sepcv 归入 Sepsis（并按 T0/T1/T2 时点细分）。
- 注：S16e 的 Tissue 列为本包软注释（非 GEO 字段直出），**Condition 映射（Sepsis/Sepsis_COVID/Control）则由前缀确定性给出**。

### 证据 4 · "ARDS" 标签的发明点（无元数据来源）

- 归档脚本 `03_FIGURES/归档/Figure4_ModuleEigengenes_Heatmap_script.py` L29–31：`if s.startswith("sepcv"): return "ARDS"`——直接把 sepcv 等同 ARDS；
- `fix_s9_ml_model_v3_honest.py` L127：`is_ards = sids.str.startswith("sepcv")`——ML "Sepsis_vs_ARDS" 任务同源假设；
- 旧 `Table_S5_PCD_Scores.csv` 的 group 列沿用同一映射，`build_Table_S5.py` L91 又将该列导出为 `GSE185263_groups.csv`，标签自此扩散。
- **推断**：早期管线将"脓毒症+COVID-19"臆断为"COVID 相关 ARDS"，未经任何临床字段验证。GEO 该系列不提供每样本 ARDS 字段，provenance 链到此中断——**结论：ARDS 标签不成立**。

### 证据 5 · 方案文档自相矛盾（已另行修订）

- 方案 line 18 数据集表写"脓毒症患者**肺组织**…82 合并 ARDS/COVID"；
- 方案 line 201 脚本注释写 "(lung, Sepsis vs HC, n=392: 348 Sepsis + 44 HC)"；
- 方案 line 4539 却写 "**GSE185263 全血数据加载**"——三处互相冲突，其中仅"全血"与 GEO 官方一致。

---

## 二、已执行的更正（2026-08-15）

| 文件 | 操作 | 数值影响 |
|---|---|---|
| `GSE185263_groups.csv` | group 列 82 行 "ARDS"→"Sepsis_COVID"（备份 `*_pre_R1_backup.csv`） | 无（同一批样本，仅改名） |
| `Table_S5_PCD_Scores.csv` | 同上 | 无 |
| `Table_S5_PCD_Scores_RAW.csv` | 同上 | 无 |
| `Table_S5_Method_Note.txt` | 头部更正"肺组织/ARDS"→"全血/Sepsis_COVID"+本文附录+历史列名备注 | 无 |
| `文稿_初稿.md` | Methods 数据来源、§2/§3/§4/§6、Discussion、Limitations 的 GSE185263 相关表述更正（全血、Sepsis_COVID、"肺组织特异性"叙事改为梯度表述） | 无 |
| `SCI论文1_分析步骤详细方案.md` | line 18 数据集表、line 201/259 等处更正 | 无 |

## 三、遗留事项（并入 D4 manifest 治理，非本项范围）

1. **历史列名残留**：`Table_S5c`（Control_z/ARDS_z/delta_ARDS_Control）、`Table_S9_Diagnostic_Report.txt`（Task B "ARDS vs Sepsis"）、`Table_S16f/S16g`、`Bridge_Test_Result.csv`（lung_log2FC 列名）等结果表的"ARDS/lung"字样为历史命名，所指均为 Sepsis_COVID 组/ GSE185263——**待 D4 统一 manifest 重生成时更名**，数值全部不变；
2. **叙述层重锚定（D1 决策）**：ARDS 主张已下锚到 GSE212865 SDRA 亚组（34 例）+ COVID-19 BALF 单细胞层；是否进一步全文正名为 "Sepsis/COVID-19" 或补真实 ARDS 队列（D1 三选一）需导师定夺；
3. sepcv 样本含 T0/T1/T2 时点（纵向采血），作为组间比较时存在重复测量未建模的潜在问题——建议在 Limitations 提及或按 T0 子集做敏感性分析。

## 四、方法学备注

- 本核查未改任何统计数值：分组样本构成（同 82 例）与打分/检验输入完全一致，仅更正组名标签；
- 备份文件（`*_pre_R1_backup.csv`）保留于 01_RESULTS_TABLES，供审计比对；
- 若后续按 T0-only 重算，属新增敏感性分析，需另立输出表。

---

## 补记（2026-08-15，D4 治理）

遗留事项 1（历史列名残留）已于 D4 全部完成：`Table_S5b`（ARDS_z→Sepsis_COVID_z 等 3 列）、`Table_S5c`（2 列）、`Table_S9_Model_Performance.csv`（task 值 ARDS_vs_Control→Sepsis_COVID_vs_Control、Sepsis_vs_ARDS→Sepsis_vs_Sepsis_COVID）与 `Table_S9_Diagnostic_Report.txt`、`Bridge_Test_Result.csv`（lung_log2FC→GSE185263_log2FC、blood_log2FC→GSE212865_log2FC）均已正名（数值零改动，写入脚本同步）；S16f 经复核无 ARDS 列（无需改动）、S16g 此前已正确。`Table_S2`/基因清单的 `ARDS_vs_Control_*` 列经查为 **scRNA COVID_severe vs Healthy** 对比（语义正确的 ARDS 代理命名），保留并在 `RESULTS_MANIFEST.md` §五 登记。各备份文件已移入 `01_RESULTS_TABLES/archive/`。遗留事项 2（D1 叙述锚定）与 3（T0 纵向）仍开放。

**补记 2（2026-08-15，S5_RAW 恢复）**：本报告第二节所列 `Table_S5_PCD_Scores_RAW.csv` 的更正在后续治理过程中意外丢失（顶层仅存 `archive/Table_S5_PCD_Scores_RAW_pre_R1_backup.csv` 修正前备份）。已从该备份恢复：逐行施加与主表相同的 R1 标签更正（group 列 82 行 ARDS→Sepsis_COVID），数值列零改动、样本序与主表逐行一致（392 行，266/82/44）；已加版本戳并入 manifest 重冻结（220 行 = 211 canonical + 9 archived），lint 组计数不变量同步纳入 S5_RAW。
