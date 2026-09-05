# P2-2 Gate 2 头对头报告（H7）

- 日期：2026-08-20；冻结依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md（OSF https://osf.io/C7RYD/）
- 签名：IIAMD v1.0 人源一对一映射 7,434（up 3,615 / down 3,819）
- 组成调整：bulk 髓系标志物均值 z（LST1/S100A8/S100A9/CD14/FCGR3A/MS4A7）作为髓系比例的操作化代理（披露）

## Phase A：GSE32707（真 ARDS 全血，ARDS_d0 vs Control）

| geneset | n_genes_avail | coverage | n_case | n_ctrl | diff | hedges_g | Welch_p | BH_q | OLS_p_adj |
|---|---|---|---|---|---|---|---|---|---|
| IIAMD_up | 3505 | 0.9696 | 18 | 34 | -0.8472 | -0.6455 | 0.07405 | 0.1234 | 0.01278 |
| IIAMD_down | 3696 | 0.9678 | 18 | 34 | -0.2019 | -0.1528 | 0.6676 | 0.7417 | 0.04915 |
| cell_cycle | 23 | 1 | 18 | 34 | -0.8247 | -0.5815 | 0.01174 | 0.04197 | 0.1192 |
| hypoxia | 17 | 0.9444 | 18 | 34 | -1.18 | -0.8517 | 0.0044 | 0.04197 | 0.01291 |
| inflammation | 20 | 1 | 18 | 34 | 0.8071 | 0.6119 | 0.02719 | 0.05438 | 0.8474 |
| interferon | 25 | 1 | 18 | 34 | 0.4625 | 0.3905 | 0.2826 | 0.4037 | 0.3924 |
| mito_stress_oxphos | 19 | 0.5278 | 18 | 34 | -1.097 | -0.83 | 0.01259 | 0.04197 | 0.002788 |
| pyroptosis_necroptosis_ferroptosis | 20 | 0.9524 | 18 | 34 | 0.321 | 0.2549 | 0.4466 | 0.5582 | 0.377 |
| technical_mt_ribo | 25 | 0.6579 | 18 | 34 | -1.128 | -0.8634 | 0.02225 | 0.05438 | 0.0005873 |
| anchor_mechanism | 19 | 1 | 18 | 34 | 0.08616 | 0.06478 | 0.8479 | 0.8479 | 0.05878 |

（Sepsis_d0_vs_Control 次对比见 _intermediate/P2_gate2_scores_GSE32707.csv）

## Phase B：scRNA 供者级 IIAMD 定位

| dataset | cell_type | score | n_case | n_ctrl | diff | p | BH_q |
|---|---|---|---|---|---|---|---|
| GSE145926 | DC | up_z | 3 | 3 | 0.7388 | 0.4347 | 0.5484 |
| GSE145926 | DC | dn_z | 3 | 3 | 0.8809 | 0.3416 | 0.5484 |
| GSE145926 | DC | IIAMD_diff | 3 | 3 | -0.1421 | 0.04164 | 0.3735 |
| GSE145926 | Mono_c14 | up_z | 6 | 3 | -0.794 | 0.4439 | 0.5484 |
| GSE145926 | Mono_c14 | dn_z | 6 | 3 | -0.5635 | 0.5918 | 0.6733 |
| GSE145926 | Mono_c14 | IIAMD_diff | 6 | 3 | -0.2305 | 0.0755 | 0.3735 |
| GSE145926 | T_cell | up_z | 6 | 3 | 1.22 | 0.1387 | 0.4031 |
| GSE145926 | T_cell | dn_z | 6 | 3 | 1.33 | 0.08108 | 0.3735 |
| GSE145926 | T_cell | IIAMD_diff | 6 | 3 | -0.1103 | 0.3692 | 0.5484 |
| GSE158055 | B_cell | up_z | 19 | 5 | 0.2538 | 0.4208 | 0.5484 |
| GSE158055 | B_cell | dn_z | 19 | 5 | 0.1961 | 0.6092 | 0.6733 |
| GSE158055 | B_cell | IIAMD_diff | 19 | 5 | 0.05774 | 0.9027 | 0.9027 |
| GSE158055 | Mono_c14 | up_z | 19 | 5 | -0.4486 | 0.0993 | 0.3735 |
| GSE158055 | Mono_c14 | dn_z | 19 | 5 | -0.2661 | 0.3885 | 0.5484 |
| GSE158055 | Mono_c14 | IIAMD_diff | 19 | 5 | -0.1825 | 0.4426 | 0.5484 |
| GSE158055 | NK | up_z | 13 | 5 | 0.4682 | 0.2815 | 0.5484 |
| GSE158055 | NK | dn_z | 13 | 5 | 1.365 | 0.01118 | 0.2349 |
| GSE158055 | NK | IIAMD_diff | 13 | 5 | -0.8972 | 0.1067 | 0.3735 |
| GSE158055 | T_cell | up_z | 21 | 5 | 0.07595 | 0.7808 | 0.8199 |
| GSE158055 | T_cell | dn_z | 21 | 5 | 0.5226 | 0.1536 | 0.4031 |
| GSE158055 | T_cell | IIAMD_diff | 21 | 5 | -0.4467 | 0.2302 | 0.5372 |
## Gate 2 初判

按冻结判据：强通过需 ≥2 个来源互斥真 ARDS 队列同向复制；本层可用真 ARDS bulk 队列仅 GSE32707 一个，因此最高只能给『部分通过』档（候选相关程序），除非后续队列核验扩充。


## Gate 2 正式判定（2026-08-20，三层证据合并）

| 判据（冻结计划 §5） | 结果 | 判定 |
|---|---|---|
| ≥2 个来源互斥真 ARDS 队列同向复制（供者级 FDR<0.05） | 唯一真 ARDS bulk 队列 GSE32707：IIAMD_up q=0.123（未达显著） | 不成立 |
| 髓系定位（供者级） | scRNA 两队列 7 细胞类型×3 评分全部 q≥0.235 | 不成立 |
| 组成调整后存活 | 发现层 GSE185263：IIAMD_up 原始 g=−1.98 (q=1.1e-23) 但髓系调整后 p=0.499；IIAMD_down 同 (p=0.557) | 不成立（关联完全由髓系组成解释） |
| 不劣于最优竞争面板 | GSE32707 最优面板 hypoxia (q=0.042)、mito_stress_oxphos (q=0.042) 均优于 IIAMD_up (q=0.123) | 不成立 |

**判定：Gate 2 未通过（FAIL）。** 实验锚定签名 IIAMD v1.0 在人类危重症中的迁移不成立：
其原始关联由髓系组成与泛炎症-代谢应激解释，不构成 Mitoxyperilysis 在人 ARDS 中的独立证据。
按冻结计划的预注册回退条款：相关发现标记为 pan-inflammatory metabolic stress；
**"Mitoxyperilysis-aligned program in ARDS" 的疾病映射主张撤回**，主稿不得写
"Mitoxyperilysis participates in ARDS"。

**阳性残余（如实报告）**：跨层最稳健的组成调整后信号为 mito_stress_oxphos 抑制
（发现层 g=−2.48, p_adj=4.6e-16；真 ARDS 层 g=−0.83, q=0.042）——与 P0-4 来源互斥 meta
（UCS g=−0.94）与 P0-6 对照表共同指向：**线粒体基础设施抑制是本项目唯一三层存活的稳健信号**，
执行臂信号属泛危重症/组成驱动。该结论与《重塑方案》§十二预测的回退分支一致：
描述性文稿（危重症线粒体-炎症解离图谱，IF 4-8 档），Mitoxyperilysis 机制归属留待湿实验（Gate 3+）。

*本判定由 P2_gate2_headtohead.py 自动+人工合议生成；全部中间数据见
_intermediate/P2_gate2_scores_GSE32707.csv、P2_gate2_scores_GSE185263.csv、P2_gate2_scRNA_localization.csv。*
