# N3 残留4 · 样本级 manifest 构建报告（Sample-Level Manifest Report）

**日期**：2026-08-16
**交付物**：`01_RESULTS_TABLES/SAMPLE_MANIFEST_v1.0.csv`（800 行 = 780 sample + 7 dataset_summary + 13 external_source；16 列）
**构建脚本**：`build_sample_manifest.py`（幂等，全部断言内置，任一不变量失败即中止写出）
**任务来源**：质控报告第三轮 N3 残留4——"逐样本（GSM→样本→数据集→条件）manifest 未建，数据源 inventory（N10）与孤儿脚本治理（N16）依赖它"；本报告同时收口 **N16 残留**（孤儿脚本处置=保留+登记入册，见 §四）。

---

## 一、覆盖范围与逐数据集定案

| 数据集 | 样本行 | 组成（全部经断言核验） | GSM 映射来源 |
|---|---|---|---|
| GSE145926（scRNA BALF） | 12 | Healthy 3 / COVID_mild 3 / COVID_severe 6（GEO 官方，N3 修复基准）；细胞 83,952 | family.soft.gz 逐 GSM（仅登记 12 个 scRNA-seq GSM；另 9 个 TCR-seq GSM 未用） |
| GSE158055（scRNA） | 51+34 | 主整合 PBMC 51 = 5 HC + 21 mild + 25 severe（GEO source_name 51/51 PBMC）；细胞 54,989；scFOCAL 气道子集 34 = BALF 12（3m+9s）+ sputum 22（5m+17s）= 57,225（42,723+14,502） | family.soft.gz（title=sampleID 精确匹配） |
| GSE171668（参考） | 6 | 原始 h5 载入 0 细胞（barcode 不匹配）；S3 组成参考实为发表版 metadata（24 供体/106,792 细胞；106,043 入表，749 compartment-NA 剔除）；**物种=GEO 官方 Homo sapiens（COVID-19 尸检肺，Delore et al. Cell 2021 图谱）——方案旧登记"小鼠"有误已订正** | 目录名 GSM5229971–76 + GEO series 记录 |
| GSE185263（全血 Bulk 主队列） | 392 | Sepsis 266 / Sepsis_COVID 82 / Control 44；前缀亚队列 sepcol67+sepnet104+sepwes84+sepvh11+sepcv82+hccol6+hchl9+hcwes24+hcwimr5（=R1 证据2 精确复现）；groups.csv × S16e 双源 392/392 一致 | **per-sample GSM 已补全（2026-08-16 VPN 恢复后）**：GEO targ=gsm 官方列表 392 个 title 与本地样本名零缺失零多余，存档于 `00_RAW_DATA/GSE185263_Lung_ARDS/GSE185263_gsm_sample_list.txt`，title 精确匹配入册 |
| GSE212865（微阵列验证队列） | 137 | Healthy 51 / COVID-19 52 / COVID-19+SDRA 34（GEO disease state）；采血时点 D0 50 / D7 36 / 对照未标注 51；**编号层面非互异：96 个患者号中 37 个编号贡献 78 个样本（D0/D7 纵向重复，原文设计内行为），6 个编号（5/18/26/66/80/85）跨对照↔疾病组——2026-08-16 三源核验判定为编号撞号（见 §三 3）** | series_matrix 逐 GSM + GEO series 摘要 + 原文 PMC10216228 |
| GSE67530（450K 甲基化） | 144 | Healthy 30 / ICU 非 ARDS 75 / ARDS 39（GEO ards 字段） | series_matrix 顺序映射 KS001–KS144 ↔ GSM1648896–GSM1649039，与 S27 分组零错配（144/144） |
| GSE165659（scATAC 基线） | 4 | 健康供体 4（主稿 Methods 口径）；raw barcodes 5,126/5,307/5,502/7,339（合计 23,274）→ QC 后 21,790 | 目录名 GSM5047856–59 |

**scRNA 合计断言**：83,952 + 54,989 = 138,941（与 lint 不变量一致）。

## 二、external_source 行（13 条，非样本型数据源登记）

- **canonical（7）**：ieu-b-69（脓毒症 GWAS，n=486,484）、ieu-b-5086（28 天死亡 GWAS）、FinnGen R10 J10_ARDS（357/406,536）、eQTLGen（n=31,470，cis-only 后 34 工具/15 基因）、gnomAD v2.1.1（141,456）、AlphaMissense、CMap/L1000；
- **reference（3）**：1000G PLINK 参考面板（AFR.fam 含 HG* 个体 ID 核验）、HLCA 肺图谱注释、cCRE/EpiAgent scATAC 参考；
- **downloaded_unused（2）**：GSE171524（IPF 空间转录组，54 GSM 已下载；B1 暂缓 ⏭️）、GSE200042（小鼠 LPS ALI 空间，随 B1 暂缓）；
- **exploratory_asset（1）**：MosMedData_CT（radiomics 孤儿脚本输入；N16 处置见 §四）。

## 三、构建过程中的四项事实更正（主稿已同步修订；①–③ 见首版，④ 为 GSM 补全时新发现）

1. **GSE67530 分组数字错误（主稿 Methods 数据来源 (7)）**：主稿原写"144 例（24 健康 + 72 ICU 非 ARDS + 48 ARDS）"。经 GEO series matrix `ards` 字段逐样本核验（并与 Table_S27 分组零错配互证）：实为 **30 Healthy / 75 ICU 非 ARDS / 39 ARDS**（合计 144 不变）。主稿已按此更正。附加登记：KS123（Excluded_data_error）与 KS050（Notable_biological_extreme）的 S27b QC 注。
2. **sepcv"重复测量"担忧与样本名证据不符（主稿局限 (10)）**：局限原文"sepcv 样本含 T0/T1/T2 纵向采血时点，组间比较存在重复测量未建模的潜在混杂"。逐样本核验显示 82 个 sepcv 编号**互异、无同一患者多时点重复采样**（时点标注 T0 21 / T1 58（含 sepcv068T1a）/ T2 1 / W1 2 为各样本单次采血时点）。主稿局限 (10) 已按证据改写为"时点构成组内不均衡的混杂"表述，撤回"重复测量"措辞（R1 报告遗留事项 3 的推测性 hedge 一并修正）。
3. **GSE212865 存在编号层面纵向重复采样（新发现，主稿 Methods (4) + 局限 (3) 已补）；跨组 6 编号已于 2026-08-16 三源核验判定为编号撞号**：GSM title 去掉 `_D0/_D7` 后缀后，137 个样本仅对应 **96 个患者编号**——37 个编号贡献 78 个样本（D0/D7 纵向重复采样），其中 **6 个编号（5/18/26/66/80/85）同时在对照组与疾病组出现**。初版曾登记"配对设计 vs 编号撞号不可裁决"，当日经 VPN 恢复后完成三源核验，**判定为编号撞号（对照与患者为两套独立去标识化编号空间，配对随访无任何设计依据）**：
   - **GEO series 摘要**（存档 `00_RAW_DATA/GSE212865_data/GSE212865_series_self_record.txt`）："We recruited **60 hospitalized patients** with microbiology-confirmed COVID-19, among whom **19 developed ARDS**. Peripheral blood was collected... **within 24 hours of admission and at day 7**"——招募句只提患者；D0/D7 重复采样为设计内行为（入院 24h 内 + 第 7 天各一次）；
   - **原文全文**（Rombauts A et al. *Biomedicines* 2023;11:1348, PMID 37239019 / PMC10216228，摘要存档 `00_RAW_DATA/GSE212865_data/PMID37239019_abstract.txt`；SuperSeries GSE212866 记录存档 `GSE212866_superseries_record.txt`）：全文 **0 处**提及 healthy/microarray——51 例健康对照未在原文任何处描述，仅存在于微阵列子系列；
   - **结构与算术互证**：对照占 series 第一个连续 GSM 块（GSM6559856–9906）且全部无时点标注、对照编号空间稀疏（5–295），患者编号空间独立（1–1320，含 1304–1320 块），6 个重叠全为 ≤85 小编号；**SDRA 组 34 样本恰来自 19 个编号（15 对 D0/D7 + 4 单时点），与原文"19 developed ARDS"精确吻合**（已入 builder 断言）。若 6 跨组编号为同一人配对随访，则对照组需有 6 人后来住院确诊 COVID——如此戏剧性设计不可能在原文只字未提。
   **ML 外部验证子集（COVID-19 vs COVID-19+SDRA，n=86）中 70/86 样本来自 35 个编号的 D0/D7 重复**——外部 AUC=0.705 的独立受试者基础弱于名义样本量，主稿局限 (3) 已补重复测量未建模的乐观性风险告诫（该告诫与撞号判定无关、依然成立）。
4. **GSE171668 物种登记错误（方案层更正）**：方案数据集表原登记"ARDS 小鼠肺组织/动物模型"，经 GEO 官方 series 记录核验为 **Homo sapiens COVID-19 尸检组织图谱**（Delore et al. Cell 2021：23 肺+16 肾+15 肝+18 心尸检供体）——该数据集 0 细胞入主整合、主稿未提及，故仅订正方案数据集表与 manifest tissue/notes，主稿无需改动。

## 四、N16 残留收口（孤儿脚本处置）

**处置决定（2026-08-16，作者确认）**：7 个孤儿脚本（`radiomics_imaging_genomics_24.1.py`、`radiomics_v2_real_data.py`、`run_mofa_multiview_analysis.py`、`run_mofa_supplementary.py`、`GEARS_Combinatorial_Perturbation_Analysis.py`、`celloracle_grn_causal_perturbation.py`、`celloracle_grn_enhanced.py`）**保留为探索性资产、登记入册（MosMedData_CT 以 exploratory_asset 状态行入 SAMPLE_MANIFEST），不删除**；是否跑完纳入（C3 档位：CellOracle in silico KO / GEARS 组合扰动）留档位 C 启动时定夺——与质控报告"机会点：GEARS/CellOracle 恰可支撑纯干路线 C3 档位"一致。P2 层至此全部关闭。

## 五、包治理状态

- `SAMPLE_MANIFEST_v1.0.csv` 经 `build_D4_package.py` 重冻结纳入 `RESULTS_MANIFEST_v1.0.csv`（版本戳 NA/NA——样本清单不依赖基因集/评分版本），`lint_package.py` 全绿（51 项——lint 同批新增样本级不变量检查：record_type 结构、逐数据集样本行数、GSE67530/GSE212865/GSE185263/GSE145926 组计数、全部 780 样本行 GSM 非空——GSE185263 392 例已于 2026-08-16 按 GEO 官方存档补全）；
- 列名经 lint 正则安全核验（无 p/FDR 误匹配列）；行数与组计数不变量已由构建脚本内置断言 + lint 双层保障；
- GSE185263 per-sample GSM 已补全（GEO targ=gsm 列表存档 `00_RAW_DATA/GSE185263_Lung_ARDS/GSE185263_gsm_sample_list.txt`，392 title 零缺失零多余）；后续样本层变更应重跑 `build_sample_manifest.py` 并重冻结。

## 六、证据文件

- `00_RAW_DATA/GSE145926_scRNA/GSE145926_family.soft.gz`、`00_RAW_DATA/GSE158055_scRNA/GSE158055_family.soft.gz`（284 GSM）
- `00_RAW_DATA/GSE212865_data/GSE212865_series_matrix.txt.gz`、`00_RAW_DATA/GSE67530_series_matrix.txt.gz`
- `00_RAW_DATA/GSE212865_data/GSE212865_series_self_record.txt` + `GSE212866_superseries_record.txt` + `PMID37239019_abstract.txt`（2026-08-16 撞号判定三源证据）
- `01_RESULTS_TABLES/GSE185263_groups.csv` + `Table_S16e_Sample_Condition_Info.csv` + `Table_S27_Epigenetic_Clocks.csv` + `Table_S2_QC_Metrics.csv`
- `00_RAW_DATA/_scrna_work/load_summary.json`（GSE171668 0 细胞载入证据）
- `04_LOGS/sample_manifest_build.json`（构建摘要）
