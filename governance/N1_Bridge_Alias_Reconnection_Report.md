# N1 核查报告：Bridge 收敛打分别名重连与重算（2026-08-15）

**对象**：`Bridge_Test_Result.csv`（主稿 §1 唯一数据源；生成脚本 `bridge_test_narrative.py`）
**性质**：P0 级数据连接 bug——既往版本按**精确符号** merge 外部数据层，权威 80 基因清单使用常用名（IP3R1/HSP60），而 gnomAD/AlphaMissense 层使用 HGNC 官方符号（ITPR1/HSPD1），导致该两行约束/变异负荷新数据漏连、由秩填充默认值支撑。
**修复**：merge 前按 HGNC 别名表归一化（ITPR1→IP3R1、HSPD1→HSP60），算法（秩归一化+等权五维平均）**零改动**；全部数值以修复后数据重算。

---

## 一、连接诊断（修复前，2026-08-15 Python 全量复核）

| 外部层 | 按精确符号连接率 | 说明 |
|---|---|---|
| S46a（GSE185263 bulk log2FC） | 30/80 | 源表基于旧 38 基因清单，缺基因为清单差异，非别名问题 |
| S46b（GSE212865 微阵列 log2FC） | 17/80 | 同上（18 行中 1 行为伪符号 "NADH"，见 R6，不涉别名） |
| **S40（gnomAD pLI）** | **3/80**（CANX、HSPA9、VDAC1） | **ITPR1↔IP3R1 断连**：S40 中 ITPR1 pLI=1.00（全表最强约束）未进 Bridge |
| **S18（AlphaMissense N_Pathogenic）** | **34/80** | **ITPR1↔IP3R1**（10,614 致病变异，全通路最多）、**HSPD1↔HSP60**（2,003）断连 |
| S12（SCENIC TF_degree） | 18/80 | HSPA9/VDAC1 等本就不是 SCENIC 靶基因，非别名问题 |
| S36（对接靶点） | 7/80 | 无别名问题 |

**别名断连仅 2 基因、3 条连接**（其余缺失均为源层基因真实缺席，不得虚构连接）：

| 权威清单符号 | 外部层符号 | 重连数据 | 对主稿的意义 |
|---|---|---|---|
| IP3R1 | ITPR1 | pLI=1.00（gnomAD）；N_Pathogenic=10,614（AM） | 主稿 §1 明文引用"ITPR1…10,614 个致病性变异（全通路最多）"，但该数据此前根本没连进 Bridge 的收敛打分 |
| HSP60 | HSPD1 | N_Pathogenic=2,003（AM） | 使 HSP60 进入收敛 Top10 |

注：S2/权威清单用 HSPA9（非 GRP75），该基因与 gnomAD/AM 的连接**原本就通**（pLI=0.97、N_Pathogenic=2,241），无需修复。

## 二、重算结果（before → after）

**模块均值排序翻转**（主收敛口径，算法不变）：

| 模块 | 修复前 | 修复后 | 变化 |
|---|---|---|---|
| **MAM_integrity** | 0.4405（第2） | **0.4578（第1）** | ↑ 接回 IP3R1 真实数据 |
| **iron_metabolism** | 0.4547（第1） | **0.4394（第2）** | ↓ 填充默认值随覆盖扩大而下调 |
| mitochondrial_function | 0.3225 | 0.3287 | |
| oxidative_stress | 0.3380 | 0.3264 | |
| ferroptosis_cuproptosis | 0.2713 | 0.2622 | |
| cell_death | 0.2576 | **0.2477** | 主稿 §1 引用值需 0.258→0.248 |
| transcription_factors | 0.2136 | 0.2045 | 0.214→0.205 |
| autophagy | 0.1768 | 0.1676 | 0.177→0.168 |

（所有模块值均受填充默认值下移影响——覆盖基因增多后维度内最低百分位秩变小，属算法内在行为。）

**基因排名**：

| 排名 | 修复前 | 修复后 |
|---|---|---|
| 1 | SLC40A1 0.7039 | SLC40A1 0.6861 |
| 2 | CANX 0.6767 | **IP3R1 0.6562**（原 0.3043，QC 预测 ≈0.656 精确兑现） |
| 3 | VDAC1 0.6438 | CANX 0.6512 |
| 4 | HSPA9 0.6270 | VDAC1 0.5963 |
| 5 | TFRC 0.5501 | HSPA9 0.5801 |
| 6–10 | SOD2、GPX4、PACS2、FTL、SLC25A37 | TFRC 0.5356、SOD2 0.5067、GPX4 0.4603、**HSP60 0.4496**（原 0.2893，新入 Top10）、PACS2 0.4465 |

- 修复后 Top10：**SLC40A1、IP3R1、CANX、VDAC1、HSPA9、TFRC、SOD2、GPX4、HSP60、PACS2**（FTL 0.4424、SLC25A37 0.4361 退至 11–12）。
- **IP3R1 与 CANX、VDAC1、HSPA9 为全表仅有的 4 个五维均有真实数据的基因**（全属 MAM 模块）。
- 排序翻转方向与论文自身叙事更自洽：MAM 居首 + IP3R1 跃居第 2，与 gnomAD Tier 3 最强约束（ITPR1 pLI=1.00）、§1 的 MAM 钙轴叙事、Discussion D2 的"约束二分性"主线一致。

## 三、N1b 敏感性：数据可得性检验（新增 TEST 2b）

**模块×维度可得性**（修复后）：

| 模块 | n | pLI 有数据 | AM 有数据 | 对接靶点 |
|---|---|---|---|---|
| MAM_integrity | 9 | **4** | 9 | 2 |
| iron_metabolism | 9 | 0 | 9 | 3 |
| mitochondrial_function | 10 | 0 | 10 | 0 |
| oxidative_stress | 11 | 0 | 7 | 2 |
| cell_death | 13 | 0 | 1 | 0 |
| ferroptosis_cuproptosis / transcription_factors / autophagy | 28 | 0 | 0 | 0 |

真实 pLI 仅 4 基因且**全部属 MAM**（IP3R1 1.00、VDAC1 0.97、HSPA9 0.97、CANX 0.87）；铁代谢模块 gnomAD 覆盖为 0。

**敏感性 A**（仅全覆盖维度 transcript+central+drug，无填充默认）：
**iron 0.4757 > MAM 0.4199** > oxidative 0.4097 > mito 0.4090 > ferro/cupro 0.3722 > cell_death 0.3277 > TF 0.2759 > autophagy 0.2146 —— **铁反超**。

**敏感性 B**（逐基因仅用真实维度，≥4/5 维有真实数据，n=36/80）：
**iron 0.5180 > MAM 0.4859** > cell_death 0.4668（n=1，APAF1，无意义）> oxidative 0.4401 > mito 0.3796 —— **铁反超**（cell_death 层内均值仅 1 基因支撑，不可比较）。

**结论（N1b）**：
1. 主口径下 MAM 居首**由约束维度驱动，而该维度只有 MAM 4 基因有数据**——MAM/铁层内先后部分由数据可得性驱动；
2. 稳健结论是 **"MAM–铁"共同构成收敛首位层级，二者在全部敏感性下均远高于死亡执行/TF/自噬模块**（敏感性 A 中首位层级与 cell_death 差 ≥0.09）；
3. 基因层面 SLC40A1（0.686）与 IP3R1（0.656）居前二，后者是五维全真实数据基因中得分最高者；
4. 上述敏感性已作为新列写入 `Bridge_Test_Result.csv`（`convergence_sens_available`、`conv_fullcov`、`n_real_dims`），并在主稿 Methods/§1/Discussion §1 如实披露。

## 四、已执行更正清单

| 文件 | 更正 |
|---|---|
| `Bridge_Test_Result.csv` | 重算重写（80 行×20 列；新增 alias_connected/n_real_dims/convergence_sens_available/conv_fullcov；备份 `Bridge_Test_Result_pre_N1_backup.csv`） |
| `bridge_test_narrative.py` | merge 前别名归一化（ALIAS 表）+ TEST 2b 敏感性 + 新列落盘（备份 `bridge_test_narrative_pre_N1_backup.py`） |
| `文稿_初稿.md` | ① Abstract Results 首句：收敛得分 MAM 0.458、铁 0.439、执行 0.248；② Methods 多组学收敛打分段：别名归一化+填充策略+两种敏感性说明；③ §1 首段：模块排序/Top10 基因清单/IP3R1 五维全数据/可得性敏感性披露；④ Discussion §1 首段：MAM 居首+可得性警示 |
| `SCI论文1_分析步骤详细方案.md` | 补 bridge 脚本登记行（含 N1 更正注） |
| `修复和优化提质20260814修改版.md` / `…待完成版.md` | N1/N1b 标记解决，优先级表/执行顺序/终审同步 |

## 五、遗留事项（并入 D4/N22 治理）

1. 权威清单其余常用名符号（如 HSP60）与 HGNC 官方符号的映射表应纳入 D4 canonical manifest（`gene_symbol_alias` 列），全表 lint 时统一校验——本次仅修复 Bridge 链路；
2. `Mitoxyperilysis_Pathway_Gene_List.csv` 的 gene_symbol 列仍用常用名（IP3R1/HSP60/HSPA9），是否整体切换为 HGNC 官方符号（ITPR1/HSPD1）属全局重命名决策，待 D4 一并定夺；
3. S46a/S46b 基于旧 38 基因清单生成（30/17 基因可连），是否按 80 基因清单重生成属 D1/D4 范畴，本次不动数值。

**证据文件**：`Bridge_Test_Result.csv`（修复后）、`Bridge_Test_Result_pre_N1_backup.csv`（修复前）、`Table_S40_gnomAD_Constraint.csv`（16 基因，ITPR1 pLI=1.00）、`Table_S18_AlphaMissense_Main_CORRECTED.csv`（55 基因，ITPR1 10,614/HSPD1 2,003）、`Mitoxyperilysis_Pathway_Gene_List.csv`（80 基因，IP3R1/HSP60 写法）。

---

## 补记（2026-08-15，D4/N22 治理）

1. 本报告遗留事项 1–3 已由 D4 治理承接：别名表（IP3R1/ITPR1、HSP60/HSPD1，并新增登记 HSPA9/GRP75、CYCS/CYTC、SLC11A2/DMT1）已入 `Mitoxyperilysis_Gene_Manifest_v1.0.csv` 的 hgnc_symbol/aliases 列，lint 符号校验按此归一化；符号策略决策为**维持常用名**（不整体切 HGNC）。
2. `Bridge_Test_Result.csv` 的 `lung_log2FC`/`blood_log2FC` 列已更名 `GSE185263_log2FC`/`GSE212865_log2FC`（数值零改动）；本报告所述备份已移入 `01_RESULTS_TABLES/archive/`。
3. 文中 `Table_S18_AlphaMissense_Main_CORRECTED.csv` 等文件现以晋升后的正式名 `Table_S18_AlphaMissense_Main.csv` 存在（内容一致，原版入 archive/）。
