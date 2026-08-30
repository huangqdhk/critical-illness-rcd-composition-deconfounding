# GSE235046 GSM↔YW 样本对应关系核验报告

**日期**：2026-08-23｜**触发**：111提质v2_未完成项标注 §七"GEO 登记 15 样本 vs 本地计数表 18 列的对应关系未文档化"

## 一、结论（先说）

1. **本地计数表实为 15 个样本列，与 GEO 登记的 15 个样本 1:1 完全对应**。审计所记"18 列"系误计：把样本编号跨度（前缀运行号 2544278–2544295 跨 18 个位置）当成了列数。
2. 计数表总列数 19 = 4 元数据列（geneID/geneSymbol/bioType/annotationLevel）+ 15 样本列。样本编号 YW001–YW012 连续 + YW016–YW018；**YW013–YW015 不存在于任何本地文件**（_intermediate 全目录检索 0 命中）——编号缺口来自源数据发布本身，非本地丢失。
3. 文稿 v4 中英双版（§4.19 与 Gate 1 段）表述"Media/LPS/CS/LPS+CS/LPS+CS+Torin 各×3""12 次 leave-one-replicate-out""count 表列块映射经 LPS 标记基因 6/6 方向校验"均与本核验一致，**文稿无需改动**。

## 二、GSM↔YW 逐样本对应表

| 组 | GEO Sample Title | GSM 登记号 | 本地计数表列 | 前缀运行号 |
|---|---|---|---|---|
| Media | BMDM_Media_01–03 | GSM7493746–7493748 | YW001–YW003 | 2544278–2544280 |
| LPS | BMDM_LPS_01–03 | GSM7493749–7493751 | YW004–YW006 | 2544281–2544283 |
| CS | BMDM_CS_01–03 | GSM7493752–7493754 | YW007–YW009 | 2544284–2544286 |
| LPS+CS | BMDM_LPS+CS_01–03 | GSM7493755–7493757 | YW010–YW012 | 2544287–2544289 |
| LPS+CS+Torin | BMDM_LPS+CS+Torin_01–03 | GSM7493758–7493760 | YW016–YW018 | 2544293–2544295 |

## 三、映射依据与限定

- series_matrix 无逐样本补充文件（`!Sample_supplementary_file_1` 全为 NONE；count 表为 series 级补充文件），故 GSM↔YW 对应**不能**由文件名直接证明。现有依据三层：
  1. **组结构对齐**：GEO 5 组×3 重复 = 本地 5 个 YW 连续块×3，组序一致（Media/LPS/CS/LPS+CS/Torin）；
  2. **运行号块连续性**：YW001–012 前缀连续（2544278–89），缺号 YW013–015（2544290–92）后 Torin 块自 YW016（2544293）起；
  3. **功能校验**：`P1_interaction_signature.py` 设计要求 LPS 标记基因 6/6 方向校验 PASS 才继续（已 PASS；文稿 §4.19 亦记载）。
- BioSample（SAMN35767461–475）与 SRA（SRX20705960–974）逐样本 ID 已备案于 series_matrix；如需铁证可经 SRA run 表对测序文件名再核（未做，当前证据链已足够）。

## 四、分析使用口径（与文稿一致性）

- **确证模型（预注册）**：12 样本 = Media/LPS/CS/LPS+CS 四组，**Torin 臂不入确证**；交互对比 β₃=μ_LPS+CS−μ_LPS−μ_CS+μ_Media（DESeq2 Wald）。
- **Torin 臂（YW016–018）**：仅探索性对比 LPS_CS_Torin vs LPS_CS（未预注册，文稿已如实标注）。
- **leave-one-replicate-out 12 次** = 12 个确证样本逐个留一（文稿"12 次"与此一致，非 15 样本全量）。

## 五、溯源

- 数据文件：`00_RAW_DATA/GSE235046/`（SHA256 见 `04_AUDIT_GOVERNANCE/RAW_DATA_DOWNLOAD_SHA256_20260823.csv`，2026-08-23 固化）
- series_matrix 样本元数据行：`!Sample_title` / `!Sample_geo_accession` / `!Sample_characteristics_ch1`（treatment）
- 权威映射代码：`P1_interaction_signature.py` L47–56（YW_GROUP 字典，即本表第三层依据）
- 本报告核验操作：计数表表头 awk 逐列抽取、series_matrix 元数据行 grep 抽取、_intermediate 全目录 YW013–015 检索（0 命中）
