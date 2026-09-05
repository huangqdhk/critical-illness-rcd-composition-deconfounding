# P0-7 空间 patient–slide–FOV 映射报告（H5 前置）

- 日期：2026-08-20；冻结依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md §4 H5
- 映射表：_intermediate/P0_spatial_patient_slide_fov_map.csv（139 行）

## 嵌套结构

Visium: 23 切片 / 23 患者（假设A下相等）；分组 {'ProliferativeDAD': 12, 'AcuteDAD': 7, 'Control': 4}
CosMx: 32 FOV / 18 Case / 4 TMA slide；每 Case FOV 数分布 {8: 10, 4: 4, 6: 3, 2: 1}
CosMx 多 FOV 共 Case：18/18 个 Case 有 >1 FOV（最大 8）——FOV 级检验构成患者内伪重复；患者级汇总可行（Case 级 n=18）

## H5 检验可行性判定

1. **Visium（GSE271370）**：23 切片在『1 切片=1 患者』假设下即患者级单元；该假设 GEO 元数据无法直接证实，投稿前须对照原文（Delorey 2021）核对是否存在共患者多切片。若假设成立，现有 23 切片合并统计即为患者级；若否，需按患者聚合后重算。
   **✅ 已核验关闭（2026-08-23）**：数据集无已发表文献（PubMed 按标题与登录号检索 0 命中；此前报告引"Delorey 2021"为错误归属，该文为另一 Nature 数据集），GEO series 登记即唯一权威来源。其 overall design 原文："a cohort of 23 lungs, including those with normal histology (n=4) and those that underwent DAD upon the course of fatal COVID-19 (n=19). … both acute (n=7) and proliferative (n=12) stages of DAD"——23 个肺对应 23 张切片、4/7/12 分组与本地完全一致，**假设成立**。旁证：逐样本字段仅 tissue/disease state、无供者标识；样本前缀 L14P（Acute DAD）与 L14C（Control）共存，证明数字前缀为编号习惯而非患者 ID。证据存档：`03_LOGS/GSE271370_gsm_brief_20260823.txt`（含 series 级与逐 GSM 级原文）。
2. **CosMx（GSE253474）**：FOV 嵌套于 Case 结构确定，原 FOV 级检验（如空间聚集置换）的有效 n 应以 Case 数而非 FOV 数计；按冻结计划，患者（Case）内汇总或患者阻断置换为正确口径，本轮先建立映射并披露，置换重算列于后续（依赖 M1_visium_scored.h5ad / CosMx 单细胞表）。
3. 按注册回退条款：在完成患者口径重算前，空间层数字以描述性身份呈现，不进入确证结论。

## P0-7b 患者级置换重算（2026-08-21 补做，完成 P0-7 遗留项）

口径：Case=患者（n=18，均死于 ARDS、无对照）；模块与单基因 z 于 Case 内标准化；kNN k=8（Case 内）；
患者级统计量 = 18 Case 的 Moran's I 等权均值；患者阻断置换 = 每 Case 内部打乱模块值后重算 I、取 Case 均值，999 次。

### 结果

| 统计 | 值 |
|---|---|
| 炎症小体模块 Case 级 I 均值（n=18） | +0.02276 |
| Wilcoxon 符号秩 vs 0 | stat=4.0, p=5.341e-05 |
| 患者阻断置换（999 次，双尾） | null mean=-0.00018 (sd=0.00190), **p < 0.001**（0/999 超过观测，按置换 p 最小值惯例报告 p < 2/1000） |
| Case 内模块×邻域髓系比例 Spearman ρ（18 Case 均值） | +0.0167（p 均值 0.334） |

### 判定

1. 炎症小体模块的空间聚集在**患者级统计单位**下仍成立（Case 级 I 均值>0，患者阻断置换 p 见上表）；
   原 FOV 级结论（4/4 TMA 聚集）的统计单位缺陷已按冻结计划补正，患者级结论取代 FOV 级结论。
2. 全部 18 Case 均为 ARDS 死亡患者（无对照）——本层结论身份仍为**病例内描述**（ARDS 肺内模块空间组织），
   不构成疾病-对照检验；与冻结计划的描述性定位一致。
3. 中间表：`_intermediate/P0_cosmx_patient_level_moran.csv`（逐 Case）、`P0_cosmx_patient_level_coloc.csv`、
   `P0_cosmx_patient_level_permutation.csv`。
