# D1b · GSE32707 元数据与平台兼容性核验报告

**核查日期**：2026-08-16
**性质**：D1b 补队列路线第 1 步核验（本地实测已下载文件，只读诊断，未改任何包内文件）
**结论**：GSE32707 平台、组定义、样本数与 GEO 登记一致，**可纳入**；80 基因可测 **74/80**（6 个 MT 基因为平台注释层真实缺席，按 N1 先例登记、不得虚构）；**d0/d7 为非配对横断面设计**（与"纵向队列"的通常理解不同）；对照组合 12 个技术重复样本需处置。

---

## 一、实测元数据（`GSE32707_series_matrix.txt.gz`，12.4MB，gzip 校验通过）

| 字段 | 实测值 |
|---|---|
| 平台 / 组织 | GPL10558（Illumina HumanHT-12 V4.0）；144 样本全为 whole blood |
| 组×时点 | untreated 对照 34；SIRS d0 21；**Sepsis d0 30 / d7 28；se/ARDS d0 18 / d7 13** |
| d0 子集 | 103（18 ARDS + 30 Sepsis + 21 SIRS + 34 对照，与文献再分析报道一致） |
| 配对结构 | Sepsis 58 样本=58 受试者；se/ARDS 31=31；**无双时点受试者 → 非配对** |
| 对照独立性 | untreated 34 样本仅 22 受试者；**subject 128115 占 12 个**（疑似技术重复） |
| 矩阵 | 47,221 探针 × 144 列纯 VALUE；无 Detection Pval 列（原始文件含，可选下载） |

## 二、80 基因可测性（`GPL10558.annot.gz`，Aug 09 2016 官方注释，EBI 镜像 SOFT 格式，48,107 行）

- **73/80 直接命中**（123 探针；塌陷规则 max-mean 照 `35.1_circulating_transcriptome_analysis.py`）；
- **GSDME 以旧名 DFNA5 存在**（探针 ILMN_1670145，Entrez 1687，矩阵 144/144 样本有值）→ 别名登记后 **74/80**；
- **6 个 MT 基因零条目**：MT-ATP6 / MT-ATP8 / MT-CO1 / MT-CYB / MT-ND1 / MT-ND2——该版官方注释完全无线粒体编码基因条目，属**平台注释层真实缺席**（同 S46a 等"源层基因真实缺席"登记口径）；
- 对比：GSE212865 Clariom D 旧 38 清单 18/38（47%）→ GSE32707 **74/80（92.5%）**，平台选择成立。

## 三、对分析设计的影响（已并入 D1b 清单）

1. **非配对 d0/d7**：取消配对建模；d0/d7 作为横断面两时点（时点作分层/协变量），主分析用 d0 子集（病例数更全）；
2. **对照技术重复**：主分析保留 34 例（发表口径），敏感性分析按受试者塌陷（22 例）；
3. **MT 基因缺席**：一致性矩阵按 74/80 计算、6 MT 独立披露；MAM/线粒体模块评分用 `filter_set` 在场基因模式（`build_Table_S5.py` 已有此模式）；文稿 EFA Factor 3 的 MT-ND1/CYB/ND2 负荷（GSE185263 层）**无法在本平台复验**；
4. **GSDME↔DFNA5**：待 D4 治理登记进 `Mitoxyperilysis_Gene_Manifest_v1.0.csv` aliases 列（同 IP3R1↔ITPR1 先例）；登记前分析脚本本地归一化即可。

## 四、遗留事项（并入第 2 步）

| # | 事项 | 处置 |
|---|---|---|
| 1 | GSDME↔DFNA5 别名登记 | D4 治理动作（改冻结 manifest + lint 归一化） |
| 2 | 6 个 MT 基因平台缺席 | SAMPLE_MANIFEST dataset 行登记 + 矩阵叙事披露 |
| 3 | subject 128115 × 12 | 主分析保留 + 按受试者塌陷敏感性 |
| 4 | Detection Pval 缺失 | 可选：下载原始文件做检测过滤；VALUE 为官方归一化矩阵，非必须 |

**证据文件**：`00_RAW_DATA/GSE32707_data/GSE32707_series_matrix.txt.gz`、`GPL10558.annot.gz`；诊断脚本 `00_RAW_DATA/GSE32707_data/_d1b_metadata_check.py`（只读，可复跑）。
