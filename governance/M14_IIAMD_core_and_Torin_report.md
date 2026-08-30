# M14 腿1/腿2 报告（IIAMD core + Torin 正式检验）

- 日期：2026-08-26；依据：M10_M13_M14_pre_registration_20260824.md（osf.io/ETVMJ）§3.1/§3.2
- 确证模型：12 样本 DESeq2 ~group，交互对比向量；参考水平=CS
- Torin 对比：LPS_CS_Torin vs LPS_CS（6 样本，DESeq2 Wald）
- 腿2 预注册口径：FDR<0.05 且 |b3|>1 且 |b3| > |b_lps|+|b_cs|（因子 1x）；目标 200–500 基因

## 1. core 构建

- full 签名：up=4367 / down=4720（合计 9087）
- **core：up=687 / down=1212（合计 1899）**
- 人源一对一映射：1383/1899（73%）
- 敏感性口径登记：{'fdr_only_up': 4367, 'fdr_only_dn': 4720, 'core_1x_up': 687, 'core_1x_dn': 1212, 'core_15x_up': 489, 'core_15x_dn': 951, 'core_2x_up': 386, 'core_2x_dn': 761, 'core_3x_up': 250, 'core_3x_dn': 517}

## 2. Gate 2 重跑（core 版，三层，同 full 版口径）

- 层1 GSE32707（ARDS_d0 vs Control）、层2 GSE185263（Sepsis_COVID vs Control）、层3 scRNA 供者级
- 全部基因集行见 _intermediate/M14_gate2_core_scores_*.csv / M14_gate2_core_scRNA_localization.csv
- 预注册条件语句：core 仍不迁移 -> 原 Gate 2 结论稳健；core 迁移 -> 原结论修正。判定见下。

## 3. Torin 正式检验（腿1）

- 统计量：签名基因在 Torin 臂的方向逆转比例；H0 = 与 baseMean 十分位匹配随机集相同
- 主检验（注册口径，P1 同款）：全池 baseMean 匹配，统计量 = sign(b_torin)!=sign(b3)；
  敏感性：方向受限池（b3 同号）内统计量 = 方向定义逆转比例（池容量足够时）
- 置换：B=2000, seed=0；单侧 p = (1+#null>=obs)/(B+1)；负对照自校准 = 零分布中心应≈0.5

| 基因集 | 方向 | n | 观测逆转率 | 零分布 mean±sd [2.5%,97.5%] | 单侧 p | 中心偏移 | 零模型 |
|---|---|---|---|---|---|---|---|---|
| IIAMD_full_up | up | 4367 | 0.6721 | 0.5693±0.0040 [0.5617,0.5771] | 0.0004998 | 0.0693 | baseMean-matched, full pool (P1-style) |
| IIAMD_full_dn | dn | 4720 | 0.6210 | 0.5693±0.0035 [0.5625,0.5761] | 0.0004998 | 0.0693 | baseMean-matched, full pool (P1-style) |
| IIAMD_core_up | up | 687 | 0.7948 | 0.5687±0.0175 [0.5342,0.6026] | 0.0004998 | 0.0687 | baseMean-matched, full pool (P1-style) |
| IIAMD_core_up | up | 687 | 0.7948 | 0.5216±0.0164 [0.4891,0.5546] | 0.0004998 | 0.0216 | direction-restricted pool (b3 same sign) |
| IIAMD_core_dn | dn | 1212 | 0.7129 | 0.5688±0.0127 [0.5429,0.5932] | 0.0004998 | 0.0688 | baseMean-matched, full pool (P1-style) |
| IIAMD_core_dn | dn | 1212 | 0.7129 | 0.6068±0.0112 [0.5842,0.6287] | 0.0004998 | 0.1068 | direction-restricted pool (b3 same sign) |

## 4. 判定

- 层1 GSE32707：core_up g=-0.84, Welch q=0.0378, OLS p=0.0405；core_dn g=+0.17, Welch q=0.727, OLS p=0.312
- 层2 GSE185263：core_up g=-1.96, Welch q=9.94e-24, OLS p=0.76；core_dn g=-1.05, Welch q=1.82e-09, OLS p=0.357
- 层3 scRNA 供者级：core 相关检验 21 项，最小 q=0.215

（判定结论见正文 Results；本报告只记录量化事实。）