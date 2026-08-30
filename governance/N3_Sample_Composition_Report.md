# N3 样本构成核查与修复报告（Sample Composition Report）

**日期**：2026-08-15（2026-08-15 二次更新：S2 溯源、scFOCAL 气道子集、QC 表标签）
**对象**：GSE145926 / GSE158055 两套单细胞数据集的样本构成、条件分组标签与组织来源；Table_S2 数据出处溯源；scFOCAL 输入子集组织构成
**结论**：主稿原 Methods 的两条样本构成声明均与 GEO 官方元数据不符；其中 GSE145926 存在**样本-条件标签错位**（3 例重症患者被标为健康对照），GSE158055 存在**组织来源错标**（51 例全部为 PBMC，非 BALF）。已按 GEO 官方元数据修复。二次核查进一步发现：Table_S2 的差异表达列实为 GSE185263 全血 Bulk 输出（曾被误标"单细胞BALF"）；scFOCAL 输入目录 BALF_filtered 的 34 个样本按 GEO source_name 为 **12 BALF + 22 sputum**（非全 BALF）。

---

## 一、GSE145926（scRNA-seq，BALF，12 例）

### 1.1 GEO 官方构成（`00_RAW_DATA/GSE145926_scRNA/GSE145926_family.soft.gz`，逐 GSM 核对）

| 样本 | GEO 官方分组 |
|---|---|
| C51、C52、C100 | healthy control |
| C141、C142、**C144** | **mild** |
| **C143**、C145、C146、**C148、C149、C152** | **severe** |

官方构成：**3 健康 + 3 轻型 + 6 重症**。

### 1.2 早期管线分组（`step1_load.py` 原分组字典，已实锤）

| 早期标签 | 样本 | 问题 |
|---|---|---|
| Healthy（6） | C51、C52、C100、**C148、C149、C152** | C148/C149/C152 官方为 **severe**，10,646 细胞被误归健康对照 |
| COVID_mild（2） | C141、C142 | C144（官方 mild，817 细胞）缺失 |
| COVID_severe（4） | C143、**C144**、C145、C146 | C144 官方为 **mild** |

### 1.3 修复（2026-08-15）

- `step1_load.py` 分组字典已按官方元数据订正（3/3/6）；
- `combined_processed.h5ad` 的 `obs['condition']` 原位修正 11,463 个细胞的标签（C148 1,999 + C149 2,349 + C152 6,298 → severe；C144 817 → mild），备份为 `combined_processed_pre_N3_backup.h5ad`；
- 修复后 GSE145926 各条件细胞数：Healthy **26,138**（原 36,784）/ mild **7,050**（原 6,233）/ severe **50,764**（原 40,935）；
- 修复后全合并矩阵（138,941 细胞）：Healthy 33,569 / mild 30,571 / severe 74,801（修复前 44,215 / 29,754 / 64,972）；
- 全部细胞级评分数值不变，仅条件分组归属更正；受影响结果表已重算：Table_S1（明细+汇总）、Table_S48（解离检验）。

## 二、GSE158055（scRNA-seq，51 例）

### 2.1 GEO 官方构成（`00_RAW_DATA/GSE158055_scRNA/GSE158055_family.soft.gz`，逐 GSM 核对）

- GSE158055（Ren et al., Cell 2021，COVID-19 单细胞图谱）共 **284 个 GSM**：PBMC 172 / Sputum 22 / PBMC-B 53 / PBMC-BT 11 / PBMC-T 13 / **BALF 12** / PFMC 1；
- severity 分布：mild/moderate 122、severe/critical 134、control 28；
- **BALF 子集仅 12 例**（S-M074-1/S-M075/S-M076-1 为 mild；S-S006/S-S008/S-S009/S-S085-1~S-S090-1 为 severe），**无健康对照**。

### 2.2 本包实际整合的 51 例

- 管线从 GSE158055 官方整合 counts 矩阵前 20 万 barcodes 均衡采样（每 subject ≤1,500 细胞），得 51 个样本、54,989 细胞；
- 逐样本核对 GEO 元数据：**51 例全部为 PBMC**（5 健康对照 S-HC008~012 + 21 mild + 25 severe）；BALF 子集的 12 个样本均不在本包样本清单内；
- 样本数与严重度前缀解析（HC/M/S）与 GEO 一致 ✓；细胞数 Healthy 7,431 / mild 23,521 / severe 24,037 ✓（与 Table_S1 一致）。

### 2.3 结论与修复

- 主稿原写"GSE158055：51 例 BALF（20 健康 + 11 轻型 + 20 重症）"**三重错误**：组织（BALF→PBMC）、样本构成（20/11/20 → 5/21/25）、以及由此产生的"两套 BALF 单细胞"叙事；
- 修复方式（R1 先例：数值零改动、标签全更正）：54,989 细胞的全部评分数值不变，仅组织归属更正为 **PBMC**；主稿 Methods/Abstract/正文所有"BALF 单细胞"表述已按"83,952 BALF + 54,989 PBMC"分层描述；
- BALF 主张现由 GSE145926（Liao et al. 2020，12 例官方 3/3/6）单独承载。

## 三、遗留与 D4 合并项

1. ✅ **已处理（2026-08-15）**：`Table_S2_QC_Metrics.csv` 条件标签已按 §1.3 修复（C148/C149/C152 → COVID_severe、C144 → COVID_mild；R1 式标签-only 改动，细胞数合计与修正后 h5ad 完全一致：7,050/50,764/26,138）；
2. ✅ **已处理（2026-08-16，R8 收口）**：scFOCAL 所用 57,225 细胞（气道子集：12 BALF + 22 sputum，11 类体系）与主整合 54,989 细胞（PBMC 51 例，8 类体系）的基数差异已在主稿 Methods scFOCAL 句注明（两套基数来自不同样本筛选与质控/分型流程、不可相互换算）；
3. `combined_processed.h5ad` 未加 tissue 列（anndata categorical 手工注入兼容性问题），组织信息以本报告 + `RESULTS_MANIFEST.md` 决策记录为准。

## 三b、N3 残留重算完成记录（2026-08-16，P2 批次）

按更正后 h5ad 完成全部条件分层结果表重算（脚本 `N3_residual_S3_S6_recompute.py` / `N3_residual_S17_postprocess.py` / `N3_residual_S7_label_fix.py` + `run_milo_correct.py` 全量重跑；备份 10 件入 archive/，`*_pre_N3res_backup.*`；manifest 重冻结 233 行，lint 全绿 50 项）：

| 表 | 处理 | 关键变化 |
|---|---|---|
| Table_S3_Cell_Type_Composition | obs groupby 重算（3,024 行全交叉结构不变、非 4 样本行块逐格一致） | GSE145926 条件计数 36,784/6,233/40,935 → **26,138/7,050/50,764** |
| Table_S6_Scissor_CellType_Composition / _Score_Summary | Scissor 打分不依赖 condition（已核 step5_2_scissor.py）→ 仅条件列重交叉表（72 行不变，边际逐格一致） | 条件计数 44,215/29,754/64,972 → **33,569/30,571/74,801**；**巨噬细胞 Scissor+ = 7/7 全部来自重症组**（原 5/7+2/7 中 2 个"健康"为 C148/C149/C152 误标细胞） |
| Table_S17 / S17b / S17c | `run_milo_correct.py` 全量重跑（3,000 邻域、5,000 次样本级置换；n=1,422 不变；组样本量 severe 31 / mild 24 / Healthy 8） | S17b FDR 显著：**T_cell 重症vs轻型 log2FC=−0.6449（FDR=0.0023）**、**新增 NK −0.8631（FDR=0.034）**；邻域层 27 个 FDR<0.05 全在重症vs健康（Mono_c14/T_cell 减少）；P<0.05 率 6.5%/31.4%/37.6%，reliable 格局不变；S17c 数值逐格不变 |
| Table_S7k / S7m / S7n（本批全面清查新发现） | S7k 逐细胞条件列字符串级更正（11,463 细胞，其余列逐字节保留）；S7m/S7n 从更正后 S7k 重新聚合 | 伪时间条件均值更新（Healthy 0.0890→0.0593、severe 0.0954→0.1067、mild 0.1550→0.1560）；计数 33,569/30,571/74,801 |
| 不受影响（已核） | S17d（无条件列）、S17e（GSE158055-only，其条件标签本就正确）、S6_Enrichment/Feature_Genes/Marker_Genes（无条件维度）、S11 逐细胞大表（无 condition 列）、S1/S1_prev38（N3 修复时已重算） | — |

Figure1–3 条件面板重出图维持暂缓 ⏭️（可视化类，作者决定跳过）。

## 四、Table_S2 数据出处溯源（2026-08-15，二次发现）

### 4.1 问题

`Table_S2_DEGs_Analysis.csv` 的 `ARDS_vs_Control_*` / `Sepsis_vs_Control_*` 列此前被主稿 §2/§3/§5、Bridge 收敛表（sc_log2FC）与旧 Table_S8 Fisher Meta 误标为"单细胞BALF 差异表达"。5.1 管线实际从未产出过单细胞 DEG 表。

### 4.2 数值指纹定案（GSE185263 全血 Bulk DESeq2）

| 指纹 | 证据 |
|---|---|
| 样本设计 | ARDS 对比 82 vs 44、Sepsis 对比 348 vs 44 → 与 GSE185263（Sepsis_COVID 82 + Control 44；Sepsis 348 + Control 44）精确匹配 |
| 表达尺度 | ACTB meanExpr ≈ 13.0（bulk 计数尺度，单细胞稀疏表达不可能） |
| 组织 | EPCAM ≈ 0.5（全血特征；BALF 上皮标志应高表达） |
| 归一化 | mean raw counts 452 ÷ DESeq2 中位数比率因子 3.3 → log1p 4.93 = 存储值 4.932（逐位吻合） |

结论：S2 差异表达列为 **GSE185263 全血 Bulk DESeq2 输出**，曾两重误标（"单细胞" + "BALF"）。

### 4.3 修复（R1 先例：数值零改动、标签全更正）

- S2 数值全部保留（真实 DESeq2 结果），主稿/方案/清单全部改标为 GSE185263 全血 Bulk；
- 新增 **Table_S2b_scRNA_DEG_80genes.csv**（脚本 `N3_scRNA_DEG_recompute.py`）：在条件标签修正后的合并 h5ad 上重算 80 通路基因真·单细胞 DEG——GSE145926 BALF / GSE158055 PBMC / 合并三口径，Mann-Whitney + BH；
- 新增 **Table_S8b_80gene_ThreePlatform_Meta.csv**（脚本 `N3_meta_80gene_rebuild.py`）：诚实三平台 Fisher Meta（单细胞合并矩阵 + GSE185263 Bulk + GSE212865 微阵列），33/80 共识（20 下调 + 13 上调）。旧版全基因组 Meta（4,536"共识基因"）曾把同一队列的两次分析（S2 ARDS + S2 Sepsis）当作两个独立数据集合并，显著性被系统性放大，该口径已弃用；
- **Bridge_Test_Result.csv** 重算（sc_log2FC 改为 S2b 的 GSE145926 BALF Δmean）：模块收敛排序更新为 iron 0.473 > MAM 0.443 > oxidative_stress 0.371 > mito 0.326 > TF 0.241 > cell_death 0.207；
- 主稿摘要/§1/§2/§3/§4/§6/讨论/方法已按新表重写（含旧 4,536 Meta 弃用披露）。

### 4.4 关键事实升级

- **TFRC**：真·单细胞层面两数据集一致**下调**（BALF −0.32、PBMC −0.29）——比旧"bulk +0.02 n.s."更有力地支持"铁输入下调"；bulk 层仍 n.s.（+0.007）；
- **呼吸链**：单细胞 BALF 8/10 线粒体基因下调（MT-CO1 −0.79、MT-CYB −0.70、MT-ATP6 −0.48），MT-ATP8 +0.64 与 HSP60 +0.24 例外上调；bulk 层 MT-ATP8 −2.74 最强下调。

## 五、scFOCAL 输入子集组织构成（2026-08-15，二次发现）

### 5.1 问题

§二十一 scFOCAL 的输入目录 `00_RAW_DATA/GSE158055_scRNA/BALF_filtered/`（57,225 细胞、34 个 sampleID）曾被整体称为"34 例 BALF 样本"。

### 5.2 GEO source_name 逐样本核对

34 个 sampleID 按 `GSE158055_family.soft.gz` 的 source_name：**12 BALF + 22 sputum**（无 HC）。

| 组织 | 样本 | 细胞数 |
|---|---|---|
| BALF（12） | S-M074-1、S-M075、S-M076-1（3 mild）；S-S006、S-S008、S-S009、S-S085-1~S-S090-1（9 severe） | 42,723 |
| sputum（22） | S-M002-1/2、S-M003-1/2/3（5 mild）；S-S001-1、S-S002、S-S003、S-S004-1/2/3、S-S005-1/2、S-S007-1/2、S-S010-1/2、S-S011-1/2、S-S012-1/2/3（17 severe） | 14,502 |

### 5.3 修复

- 方案 §21.1、`Table_S22_scFOCAL_Analysis_Report.txt`、`scFOCAL_S22_python.py` 头部注释已改标"气道样本（12 BALF + 22 sputum）"；数值零改动（评分全为真实计算）；
- 主稿 scFOCAL 句未含组织断言（"基于GSE158055 57,225个细胞"），无需改动；方案数据集表 GSE158055 行已注明该子集构成。

## 六、证据文件

- `00_RAW_DATA/GSE145926_scRNA/GSE145926_family.soft.gz`（12 GSM 逐样本 patient group）
- `00_RAW_DATA/GSE158055_scRNA/GSE158055_family.soft.gz`（284 GSM 逐样本 source_name/severity）
- `00_RAW_DATA/_scrna_work/step1_load.py`（原分组字典；已订正）
- `00_RAW_DATA/_scrna_work/step1_load_GSE158055.py`（均衡采样逻辑，51 subject 来源）
- `00_RAW_DATA/_scrna_work/combined_processed.h5ad`（obs condition 已修正；备份 `*_pre_N3_backup.h5ad`）
- `00_RAW_DATA/GSE158055_scRNA/BALF_filtered/cell_annotation.csv`（34 sampleID × 细胞，用于 §五 计数）
- `N3_scRNA_DEG_recompute.py` / `N3_meta_80gene_rebuild.py` / `N3_check_S2_deg_labels.py`（本次新增脚本）

## 四、证据文件

- `00_RAW_DATA/GSE145926_scRNA/GSE145926_family.soft.gz`（12 GSM 逐样本 patient group）
- `00_RAW_DATA/GSE158055_scRNA/GSE158055_family.soft.gz`（284 GSM 逐样本 source_name/severity）
- `00_RAW_DATA/_scrna_work/step1_load.py`（原分组字典；已订正）
- `00_RAW_DATA/_scrna_work/step1_load_GSE158055.py`（均衡采样逻辑，51 subject 来源）
- `00_RAW_DATA/_scrna_work/combined_processed.h5ad`（obs condition 已修正；备份 `*_pre_N3_backup.h5ad`）
