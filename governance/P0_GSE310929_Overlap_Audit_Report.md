# P0-3 GSE310929 样本重叠审计报告

- 日期：2026-08-20
- 依据：《黄裕荣创新提质方案 v2（重塑版）》§3.1 meta 重复纳入（致命）
- 判定源：GSE310929_AllSampleMetadataSubmitted.xlsx（Dataset/GEO Accession/Patient ID）
- 脚本：P0_GSE310929_overlap_audit.py

## 判定：OVERLAP_CONFIRMED

## 与 M2 其余队列的来源级重叠

| 来源队列 | GSE310929 内例数 |
|---|---|
| GSE185263 | 392 |
| GSE32707 | 144 |
| GSE212865 | 0 |
| GSE148871 | 0 |
| GSE188309 | 0 |

## 原始来源数据集全景（28 个）

| Dataset | n |
|---|---|
| GSE65682 | 479 |
| GSE236713 | 447 |
| GSE185263 | 392 |
| GSE134347 | 298 |
| GSE189400 | 201 |
| Faheem | 193 |
| GSE236892 | 189 |
| GSE69063 | 153 |
| GSE32707 | 144 |
| GSE131761 | 129 |
| GSE95233 | 120 |
| GSE74224 | 105 |
| GSE154918 | 105 |
| SRP198776 | 105 |
| GSE63311 | 83 |
| GSE57065 | 80 |
| GSE216902 | 74 |
| GSE13015 | 67 |
| GSE222393 | 58 |
| GSE137340 | 55 |
| GSE66890 | 54 |
| GSE100159 | 47 |
| GSE196117 | 40 |
| GSE232753 | 28 |
| GSE199816 | 27 |
| GSE33118 | 20 |
| GSE211210 | 10 |
| GSE232404 | 10 |

## 审计日志

```
== 1. GSE310929 逐样本元数据（唯一判定源）
   shape=(3713, 68)
   关键列存在性: Dataset=True, GEO Accession=True, Patient ID=True

== 2. 原始来源数据集全景（n=28 个）
   GSE65682	479
   GSE236713	447
   GSE185263	392
   GSE134347	298
   GSE189400	201
   Faheem	193
   GSE236892	189
   GSE69063	153
   GSE32707	144
   GSE131761	129
   GSE95233	120
   GSE74224	105
   GSE154918	105
   SRP198776	105
   GSE63311	83
   GSE57065	80
   GSE216902	74
   GSE13015	67
   GSE222393	58
   GSE137340	55
   GSE66890	54
   GSE100159	47
   GSE196117	40
   GSE232753	28
   GSE199816	27
   GSE33118	20
   GSE211210	10
   GSE232404	10

== 3. 来源级重叠（meta 双计数风险判定）
   [重叠] GSE185263: n=392 例进入 GSE310929
          Disease Simplified 构成: {'Sepsis': 348, 'Control': 44}
   [重叠] GSE32707: n=144 例进入 GSE310929
          Disease Simplified 构成: {'Sepsis': 58, 'Intubated subjects undergoing mechanical ventilation': 34, 'Sepsis - ARDs': 31, 'SIRS': 21}
   [无重叠] GSE212865: 0 例
   [无重叠] GSE148871: 0 例
   [无重叠] GSE188309: 0 例

== 4. GSM 层独立复核（M2 分数表实际用到的样本号）
   GSE32707: M2 用样 144，与 GSE310929 的 GSM 交集 = 144 （示例 ['GSM812609', 'GSM812610', 'GSM812611', 'GSM812612', 'GSM812613']）
   GSE212865: M2 用样 137，与 GSE310929 的 GSM 交集 = 0
   GSE148871: M2 用样 304，与 GSE310929 的 GSM 交集 = 0
   GSE188309: M2 用样 198，与 GSE310929 的 GSM 交集 = 0
   GSE185263: M2 用样 392（sepcol* 患者号），GSE310929 Patient ID 含 sepcol 模式 = 0

== 5. M2_per_sample_scores.csv 自带 Dataset 列的一致性
   GSE310929 行数=3713，Dataset 列取值 n=28
   GSE65682	479
   GSE236713	447
   GSE185263	392
   GSE134347	298
   GSE189400	201
   Faheem	193
   GSE236892	189
   GSE69063	153
   GSE32707	144
   GSE131761	129
   GSE95233	120
   GSE154918	105
   GSE74224	105
   SRP198776	105
   GSE63311	83
   GSE57065	80
   GSE216902	74
   GSE13015	67
   GSE222393	58
   GSE137340	55
   GSE66890	54
   GSE100159	47
   GSE196117	40
   GSE232753	28
   GSE199816	27
   GSE33118	20
   GSE211210	10
   GSE232404	10
   核对 GSE185263: xlsx=392 vs M2分数表=392 → 一致
   核对 GSE32707: xlsx=144 vs M2分数表=144 → 一致
   核对 GSE212865: xlsx=0 vs M2分数表=0 → 一致
   核对 GSE148871: xlsx=0 vs M2分数表=0 → 一致
   核对 GSE188309: xlsx=0 vs M2分数表=0 → 一致

== 6. 审计结论（自动判定）
   判定：重叠成立 → {'GSE185263': 392, 'GSE32707': 144}
   受污染分析：M2_meta.csv（k=3 随机效应 meta）与 GSE310929 的全部单队列对比；
   修复路径（P0-4）：以 Dataset 列拆回原始数据集 → 数据集级效应 + 来源互斥 + leave-one-source-out。
```
