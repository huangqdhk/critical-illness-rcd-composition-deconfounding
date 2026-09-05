# RESULTS_MANIFEST_v1.0 — 结果表包冻结说明（D4/N22 治理）

**2026-08-29 DOI 回填与注册 ID 纠错（投稿准备批次）**：三份 OSF 注册已于 2026-08-29 解禁公开并取得 DOI——**C7RYD**（冻结分析计划，10.17605/OSF.IO/C7RYD，2026-08-20 注册）、**ETVMJ**（M10/M13/M14 二次数据预注册，10.17605/OSF.IO/ETVMJ，2026-08-24 注册）、**98CM3**（M11/M12 预注册，10.17605/OSF.IO/98CM3，2026-08-27 13:45 注册，OSF 官方 date_registered=2026-08-27T05:45:15Z=北京时间 13:45:15）。DOI 已回填四处：注册文档头（M10_M13_M14_pre_registration_20260824.md / M11_M12_pre_registration_20260827_draft.md）、双版文稿注册表述（版本说明、4.22、4.25）、本文件变更行、冻结计划 DOI 占位行（P0_FROZEN_ANALYSIS_PLAN_v1.0.md）。本地文档误记的注册 ID（3y2rb、3e92j）已全局纠正为真实 ID（ETVMJ、C7RYD；39 文件 108 处，替换日志 _intermediate/DOI_replace_log.txt）；OSF API 核验 date_registered=2026-08-27T05:45:15Z = 北京时间 13:45:15，A1 时序裁定获官方记录二次确认。2026-08-29 终稿按作者指示：双版文稿内注册号（osf.io/xxx）全部移除、仅保留 DOI 表述。

**2026-08-29 B1 溯源修订（v4.7 终审未竟项批次）**：`M13_Confounding_Audit_Report.md` 数值修订（§1 补录 myeloid 零模型 bootstrap 数值对：null beta_myeloid=1.612±0.063、经验 p=0.0005；来源 `M13_step2_null_bootstrap.py` B=2000 seed=0，B=10,000 稳健核验 p=0.0001/零分布 1.613±0.062），manifest 行 SHA256/size/n_data_rows 同步更新；新增 `03_LOGS/M13_null_bootstrap_myeloid_log.txt`（核验日志，03_LOGS 不纳入 manifest 登记）。**事故记录**：本次修订过程中 manifest CSV 曾因 pandas 往返写回将字面量 "NA" 单元格误写为空串（并引入浮点化/换行符差异），已按 `P0_rebuild_manifest_v2.py`（v1.0 基底 + 磁盘实况）重建参考文件逐格比对修复（恢复 955 格），notes/analysis 治理文本保留；同时恢复了 M1x_step12 清理批次 64 行的登记约定（3 个更新行 location/status=archived + 61 个新增行 class/版本列空串），lint 复绿 74/0。

**2026-08-21 v2.1 变更记录（投稿前硬阻断项收尾批次）**：RESULTS_MANIFEST_v2.0.csv 以 v1.0（279 行，2026-08-17 冻结）为基底重建为 591 行（重建脚本 `P0_rebuild_manifest_v2.py`；v2.0 曾因修复脚本写回异常被截断，重建以磁盘为权威重算全部 canonical 文件 SHA256，v1.0 文件保留历史冻结链）。本批次变更：① 修正 M2/M3/M4 批次引入的版本戳漂移（manifest 空串 vs 冻结 CSV 'NA'，164 行）；② 5 个 CSV 补版本列（Figure_10D、S15e Instruments/LeaveOneOut/Sensitivity、S60）+ 4 个 Figure CSV（1B/1C/2E/S2F）；③ N22 s_number 归一（S15e/f/g→S15、S19i→S19）；④ archive 子目录重组后 location 修正（R1_table_backups/M1_check_debug/R1_artifacts）；⑤ 补登记 P0/P1/P2 系列与 2026-08-21 新增治理文件（P0-8、查新 v1.1/重跑日志、参考文献核实报告、manuscript-to-source 报告、软件快照、pip freeze×2、P1-4 报告等）；⑥ lint 同步更新（FDR_NULL_OK 增 M4 蛋白层 4 表与 Figure_5I；版本戳多值 ';' 集合语义；archive 递归比对）；**lint 全绿：74 项通过 / 0 失败**（2026-08-21）。

**2026-08-31 M16 力学边界模块批次**：新增第九条候选边界（力学边界）全套产物——5 个附表（S88 mech_v1.0 模块 23 基因三层+PMID 活验溯源 / S89 四套 VILI 同口径对比+合并 / S90 GSE2411 力学×LPS 2×2 互作 / S91 人体桥接+Visium 空间 / S92 同源映射覆盖审计）+ 3 个主图数据（Figure_11A/B/C）+ 2 个治理文件（GEO 核验、执行报告）；预注册 M16_pre_registration_20260831.md 按 M1/M2/M10 纪律不入册（根目录治理文档）；8 个 CSV 按 M10–M15 导出设计登记为无内嵌版本列（manifest 记 NA，lint 白名单 M1X_NO_VERSION_STAMP 同步补登）；双版文稿（v5 中文/英文版）同步更新：数据集表 +4 行、方法 4.28、结果第 17 节、讨论四向→五向边界、图 11 注、参考文献 59–66（Vancouver，PMID 全经 eutils 活验）、数据可得性。**lint 全绿：74 项通过 / 0 失败**（2026-08-31）。判据裁定：A 阴性、B 形式未达（直接增强显著）、C 阴性/混合、D 达提示级——阴性/部分结果如实界定，不升级措辞。

**变更记录（References 落盘批次，2026-08-29）**：本批次变更：① P0_References_1to64_Verification_Report.md 修订（11 条待确认条目裁定记录 + round5 系列检索留痕 + GB/T→Vancouver 落盘与 README 退役记录），manifest 同步 size/sha256（8861→16653）；② 补登记 2 个治理文件：P0_Manuscript_Source_Check_Report_v2.md（§四-2 扩展终检 13/13 PASS，含 References 完整性断言）、P0_Novelty_Search_Rerun_Log_20260829.md（§四-1 查新重跑，空白维持）；修订方式为行级精确写入（非 pandas 往返，规避 B1 批次误伤模式）。**lint 全绿：74 项通过 / 0 失败**（2026-08-29）。

**冻结日期**：2026-08-15
**机器可读清单**：`RESULTS_MANIFEST_v1.0.csv`（222 行 = 211 canonical + 11 archived；含 SHA256 校验和；2026-08-15 N3 修复后重冻结新增 S2b/S8b/N3 审计文件，同日补恢复 R1 期间丢失的 `Table_S5_PCD_Scores_RAW.csv`（82 行 ARDS→Sepsis_COVID 与主表同步，数值零改动）；**同日 N2 臂定义统一后再次重冻结**——S5b 增 `arm` 列、S5 方法注双臂归组更正，2 个 `*_pre_N2_backup.*` 入 archive/；**同日 R7 结构注释更正后再冻结为 223 行 = 211 canonical + 12 archived**（`Table_S36c` PDB 注释按 RCSB 官方订正，备份入 archive/，见 §七c）；**2026-08-16 P2 批次（N3 残留重算）后冻结 233 行 = 211 canonical + 22 archived**；**同日 P3 批次（N3 残留4）后再冻结 235 行 = 213 canonical + 22 archived**——新增 `SAMPLE_MANIFEST_v1.0.csv`（样本级 manifest，800 行 = 780 sample + 7 dataset_summary + 13 external_source，见 `N3_Sample_Manifest_Report.md`）与该报告；lint 同批新增样本级不变量检查）；**同日 D1b 批次（GSE32707 外部验证，S50–S54）后冻结 242 行 = 220 canonical + 22 archived**（同载 `03_LOGS/D4_build_log.txt`）；**同日 S6 旧编号 PNG 归档后再冻结 242 行 = 219 canonical + 23 archived**——`Figure_S6_Dissociation_hexbin.png`（旧 S6 编号遗物，解离检验图现编号 S48，见 `01_FIGURE_DATA_CSV/README_图数据包说明.md`）移入 archive/，数值零改动）；**同日目录重组（导师版结构）后再冻结**——`01_RESULTS_TABLES/` 拆分为 `02_SUPPLEMENTARY_TABLES/`（`SUPPLEMENTARY_Tables_CSV` + `Supplementary_Notes`）与 `04_AUDIT_GOVERNANCE/`，`archive/` 上移根目录，日志目录 `04_LOGS`→`03_LOGS`；manifest 新增 `location` 列；`Table_S5_PCD_Scores_RAW.csv` 原文件未随包迁移、由 canonical 脚本 `build_Table_S5.py` 重生成并按 §八 重冻结（z 值逐样本一致 392/392，见该行 notes）；`Table_S18` 三表自 archive 字节级恢复（SHA 核对一致）；`Table_S1` prev38 ×2 与 `Table_S11` 大矩阵 ×2 登记缺席留痕（与 lint `MISSING_OK` 白名单同源）
**校验入口**：项目根目录 `python lint_package.py` → `03_LOGS/lint_report.md`（当前 **全绿**：60 项通过 / 0 失败）
**构建脚本**：`build_D4_package.py`（幂等，可重跑；操作日志 `03_LOGS/D4_build_log.txt`）

对应质控报告条目：**D4 治理工程化**（canonical 基因 manifest + 全表版本列 + lint）与 **N22 补充表编号冲突与版本并存**。

---

## 一、N22 编号去重决策（每个 S 编号唯一对应一个分析族）

| 冲突 | 决策 | 理由 |
|---|---|---|
| S6 双占用（Scissor vs 单细胞解离检验） | **Scissor 保持 S6**；解离检验改 **S48**（`Table_S48_Dissociation_singlecell_correlation.csv` + `Table_S48_Dissociation_Method_Note.txt`） | 方案 §8.1 原始登记 S6=Scissor（步骤 6、主图 Figure 7）；解离检验为 2026-08-13 新增分析，顺延取新号。文稿 2 处引用同步改 S48（含 Figure S6→S48） |
| S3 双占用（WGCNA 模块 vs 细胞类型组成） | **细胞类型组成保持 S3**；WGCNA 模块改 **S49**（`Table_S49_WGCNA_Modules.csv`） | 文稿 line 139 的 "Table S3" 指细胞组成（引用稳定性优先）；WGCNA 表未被以编号引用，零引用风险 |
| `Table_s12_*`（小写 s）与 S12 并存 | 改 **`Table_S12b_Mitoxyperilysis_Upstream_TF.csv`** | 并入 S12 子表序列（main + Summary + b），消除大小写双命名空间 |
| S26b 双占用（Cluster_Statistics vs HLCA 注释） | **HLCA 注释改 S26a**（方案 §26 行内既定建议）；Cluster_Statistics 保持 S26b | S26a 空闲；b–k 序列不变 |
| S18 `_CORRECTED` 与原版并存 | **`_CORRECTED` 晋升为正式表名**（Main/S18b/S18c）；原非修正版（内容已逐字节一致，2026-08-15 修正时同步）移入 `archive/*_superseded.csv` | 方案 §8.1 既定"以 CORRECTED 为准"；消除投稿包内版本并存 |
| 备份文件与顶层混杂 | 6 个 `*_pre_R1/R2/N1/fix_backup*` 全部移入 `archive/` | 备份不进入投稿包命名空间；审计链保留 |

**编号总况**：S1–S9、S11–S24、S26–S33、S36–S41、S44–S49 共 43 个在用编号 + S10/S25（跳过：空间转录组/影像组学无数据）+ S34/S35/S42/S43（未启用）。多对一映射（合法）：`Table_S22_S23_*`、`Table_S29_S31_*`（报告）、`Table_S32_S33_*`（方法注）、`Table_S38_S39_*`（报告）。

## 二、canonical 基因 manifest（D4 第一交付物）

`Mitoxyperilysis_Gene_Manifest_v1.0.csv`：80 基因 × 12 列（version / gene_symbol / hgnc_symbol / aliases / ensembl_gene_id / gene_name / module / module_size / arm / is_mt_gene / ARDS_vs_Control_log2FC / ARDS_vs_Control_padj）。

- **符号策略决策**（N1 遗留 2）：canonical 符号**维持现用常用名**（IP3R1、HSP60），不整体切换 HGNC——全部管线与结果表已用此写法，切换无分析收益且有断链风险；HGNC 官方名以 `hgnc_symbol` 列提供，别名以 `aliases` 列登记。
- **包内别名登记**：IP3R1↔ITPR1（gnomAD/AlphaMissense 层）、HSP60↔HSPD1（同前）、HSPA9↔GRP75（ATAC/转录本层）、CYCS↔CYTC（ATAC/TF 结合层）、SLC11A2↔DMT1（MR 层）。lint 符号校验按此归一化。
- **双臂归属**：上游塌陷臂 30（MAM 9+线粒体功能 10+铁死/铜死 11）、执行诱导臂 33（铁代谢 9+cell_death 13+氧化 11）、不在臂内 17（TF 11+自噬 6）——与 `singlecell_dissociation_test.py`（Table_S48，上游臂经别名归一化 30/30 全部命中，早期 28/30 系 IP3R1/HSP60 别名不一致）/ bulk Table_S5b / 主稿 §3 bulk log2FC 叙事为同一划分（**2026-08-15 N2 臂定义统一**：此前主稿叙事曾用"MAM+铁+线粒体+氧化 39 基因"旧口径、S5b 方法注曾按 Δz 符号事后归组把 TF/自噬并入上游组，均已统一为本 manifest `arm` 列；S5b 已增 `arm` 列）。
- 模块尺寸不变量：9/10/9/11/13/11/11/6 = 80（lint 强校验）。

## 三、版本戳（gene_set_version / score_version）

全部 173 个 canonical CSV 追加两常量列；21 个 TXT 报告追加版本页脚。取值：

| 版本号 | 含义 |
|---|---|
| `Mitoxy-80_v1.0` | 依赖 80 基因清单的分析（本 manifest 冻结版） |
| `PCD-6panel_v1.0` | S5 六种程序性死亡评分基因面板 |
| `score_genes_v1` | scanpy score_genes 逐细胞评分（S1、S6 Score_Summary、S7q、S48） |
| `ssGSEA_v1` | running-sum ssGSEA（β 参数见表内列；S5/S5b/S5c/S16f） |
| `bridge_v2` | Bridge 五维收敛（v2 = N1 别名归一化修复后重算） |
| `NA` | 不依赖基因集/评分版本（分组定义、全基因组 DEG 附属表、注释表等） |

**一致性校验**：lint 逐行核对每个 CSV 的版本列值 == manifest 记录值；TXT 页脚版本 == manifest。

## 四、登记的例外（lint 白名单，均附理由）

1. **TOMM70A**（层内扩展基因）：ATAC（S26f/S26j）、GWAS（S19c）、转录本（S41*）层的基因面板自带，非 80 清单成员——真实基因、真实数据，保留并登记。
2. **NADH**（伪符号）：GSE212865 微阵列探针注释伪影（R6 已登记）+ 转录本层沿袭（S41a–d/f）——保留源数据原样，不参与符号合法性判定。
3. **S40 对照 12 基因**：BAD/BAK1/BAX/BBC3/BCL2/BCL2L1/BCL2L11/BCL2L2/BID/MCL1/MCU/PMAIP1——正文 §1 Tier 比较的凋亡/MAM 钙对照，设计内扩展。
4. **S18 自选面板**：AlphaMissense 55 基因面板 = 80 清单成员 + 历史扩展基因（C2CD5、PRDX4、TOP2A 等），登记为 custom panel，不做 80 子集强校验。
5. **FDR/p 空值白名单**：S47b `FDR`（6 个恒定比例细胞类型相关未定义 → NaN，`FDR_note` 列逐行说明）；S8 `GSE212865_padj`/`meta_Fisher_padj`（基因不在该平台/某数据集 p 缺失，数据集缺席 NA）；S15a `fdr_q`（逐 IV 敏感度行，无 FDR 义务）。
6. **不加版本列的文件**：`GSE185263_groups.csv`（分组定义，读取方按列 merge）、`MR_bio_*.csv`（探索性 MR 结局侧，状态待定）、`GSE67530_beta_for_GrimAge.csv.gz`（输入数据）。

## 五、历史列名统一（R1 遗留①，数值零改动）

| 文件 | 改动 |
|---|---|
| `Table_S5b` | `ARDS_z`→`Sepsis_COVID_z`、`ARDS_pct_pos`→`Sepsis_COVID_pct_pos`、`delta`→`delta_Sepsis_COVID_vs_Control` |
| `Table_S5c` | `ARDS_z`→`Sepsis_COVID_z`、`delta_ARDS_Control`→`delta_Sepsis_COVID_vs_Control` |
| `Table_S9_Model_Performance.csv` | task 值 `ARDS_vs_Control`→`Sepsis_COVID_vs_Control`、`Sepsis_vs_ARDS`→`Sepsis_vs_Sepsis_COVID`；note 同步 |
| `Table_S9_Diagnostic_Report.txt` | Task A/B 标题与任务定位文字按 R1 正名（Sepsis_COVID 全血，非"ARDS 肺"） |
| `Bridge_Test_Result.csv` | `lung_log2FC`→`GSE185263_log2FC`（S46a 全血，历史"肺组织"误标）、`blood_log2FC`→`GSE212865_log2FC`（按数据集锚定） |

**保留不改的历史命名（登记备查）**：`Table_S46a` 的 `COVID_*` 列 = GSE185263 sepcv（Sepsis_COVID）——命名风格差异已在此登记，改动将级联大量读取脚本，收益为零。

**N3 溯源更正（2026-08-15，撤销本 manifest 旧登记）**：`Table_S2`/基因清单的 `ARDS_vs_Control_*` 列**不是** scRNA COVID_severe vs Healthy（本文件此前登记有误）。经数值指纹核对（ACTB≈13.0 为 bulk 计数尺度、EPCAM≈0.5 符合全血、样本设计 82/44 与 348/44 匹配 GSE185263、meanExpr=4.932 = log1p(452/3.3) 为 DESeq2 中位数比率归一化），该表实为 **GSE185263 全血 Bulk DESeq2 输出**（ARDS_vs_Control = 82 sepcv vs 44 ctrl；Sepsis_vs_Control = 348 sepsis vs 44 ctrl）。数值零改动，标签层面按 R1 先例更正；真·单细胞 DEG 由新增 `Table_S2b` 承担，旧全基因组 Fisher Meta 由新增 `Table_S8b` 替代（详见 `N3_Sample_Composition_Report.md`）。

## 六、lint 校验的不变量（计数与组定义）

GSE185263 = Sepsis 266 / Sepsis_COVID 82 / Control 44（共 392；groups、S16e、S5、S5_RAW 四处一致）；scRNA 合计 138,941 = GSE145926 83,952 + GSE158055 54,989（S1_Summary 与 S48 pooled 双处核验）；GSE145926 条件细胞数 26,138/7,050/50,764（N3 标签修复后，S2_QC 合计与 h5ad 一致）；S6 Scissor 组成比例逐组合计 100%；基因 manifest 80 行/8 模块/双臂 30-33-17；Bridge 80 行；**S2b 80 行（符号列严格校验）；S8b 80 行（符号列严格校验）**；S12b 23 TF；S13 12 KO；S15d 15 基因 34 cis+13 trans；S18 55 行；S46e 38 行；全部 p 值列 ∈[0,1]、FDR 列非空（白名单除外）；210 个 canonical 文件 SHA256 与本 manifest 一致。**（2026-08-16 S6 旧编号 PNG 归档后为 219 个。）**

## 七、未纳入本次冻结范围（显式登记）

1. **S46a/S46b/S46e 基于旧 38 基因清单生成**（N1 遗留 3）——是否按 80 清单重生成属 D1 数据决策，本次未改数值。
2. **`03_FIGURES/Figure*.csv`（140 个）**：历史出图用工作副本（独立 Figure 编号体系），**非补充表包成员**；其中 S16/S17/S18 等对应副本早于 2026-08-14 修复，出图前须自 canonical 表重新导出。
3. 叙述层重锚定（D1）、S40 gnomAD 提升（D2）等见质控报告路线图。

## 七a、旧版评分尺度标注（R5 处置，2026-08-15）

- `Table_S19h_CellType_Mitoxy_Score.csv` 与 `Table_S22b_scFOCAL_CellType_IC50.csv` 为**旧版评分尺度**（S19h：AUCell 式 0.46–0.90 值域，排序与现行 S1 score_genes 不一致；S22b：~20 量级评分 + 11 类含 Neu/Plasma/Mast 的细胞类型体系，与全稿 8 类 BALF 注释不兼容）——两表 `score_version` 列已标 `legacy_scale_v0_not_for_text`，**仅供历史参考，不用于正文数值引用**（R5"标注"选项；如需"重算"选项属后续工作量）。`Table_S1_rescore_report.txt` 标题已正名（实现为 score_genes，原"AUCell"名不符实）。

## 七b、登记的脚本侧同步（重跑一致性）

以下写入/读取脚本已同步至 canonical 命名（仅路径与标签文字，无算法改动）：`singlecell_dissociation_test.py`（S48）、`bridge_test_narrative.py`（GSE185263/GSE212865 列名+平台标签）、`build_Table_S5.py`（Sepsis_COVID 组名+新列名+叙述）、`NicheNet_Ligand_Receptor_Target_Analysis_v2.py`（S12b）、`Mitoxyperilysis_ARDS_GSE185263_WGCNA_v2.py` / `WGCNA_v2_rerun.py`（S49）、`annotate_cell_types_lung_atlas.py`（S26a）、`AlphaMissense_S18_Fix.py`（正式名输出）、`fix_s9_ml_model_v3_honest.py`（任务标签）、`rescore_mitoxy_80gene_aucell.py`（AUCell 注记）。

## 七c、Table_S36c 结构注释更正（R7 处置，2026-08-15）

`Table_S36c_Target_Protein_Info.csv` 的 Resolution/CoLigand/Function 列经 **RCSB 官方 REST API 逐条复核**后更正（备份 `archive/Table_S36c_pre_R7_backup.csv`，SHA 与冻结前原版一致）：

| 靶蛋白 | PDB | 旧（错误）值 | 新（RCSB 官方）值 |
|---|---|---|---|
| GPX4 | 5H5Q | 2.0Å；CoLigand "GXP" | **1.1Å**（X-ray）；CoLigand **GXpep-1**；Function 注明"结晶构建体含 Sec/Cys 替换" |
| VDAC1 | 2JK4 | 3.0Å | **4.1Å**（X-ray，Bayrhuber PNAS 2008）——注：质控报告 R7 曾称 2JK4 为"NMR 溶液结构"，经 RCSB 复核**该说法有误**（2JK4 为 X-ray 4.1Å；VDAC1 的 NMR 结构为 2K1N，未用于对接） |
| CANX | 1JHN | 未注物种/片段 | 2.9Å（不变）；Function 注明"**犬源**（*Canis lupus familiaris*）**腔域片段**，无跨膜/胞质段" |
| FTL | 2FHA | 2.5Å | **1.9Å**（X-ray） |
| TFRC | 1SUV | 3.2Å | **EM 7.5Å**（cryo-EM，人 TfR–transferrin 复合物；低分辨率模型，Function 已注明） |
| SOD2 | 1N0J | 2.2Å | 2.2Å（核验无误，未改） |
| SLC40A1 | AlphaFold | pLDDT=80.25 | 不变（预测结构标识原本正确） |

行数（7）与对接数值（S36/S37）零改动；仅注释列更正。主稿 Methods 分子对接段已同步披露结构证据分层（见文稿"深度学习与分子模拟"节）。manifest 重冻结 223 行 = 211 canonical + 12 archived，lint 全绿（50 项）。

## 八、重跑与再冻结

任何结果表变更后：`python build_D4_package.py`（仅对新文件补版本戳并重冻结 manifest）→ `python lint_package.py`（必须全绿）。冻结清单自身变更应升版本号（v1.1）并保留 v1.0 审计链。
