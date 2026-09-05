# P0-6 原结果 vs 正确统计单位结果 对照表

- 日期：2026-08-20；冻结依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md（OSF https://osf.io/C7RYD/）
- 判定规则：旧 k=3 meta 与 GSE310929 独立队列地位作废（P0-3 审计）；其余为重算或保留

| layer | item | old | new | verdict |
|---|---|---|---|---|
| meta-pooled | UCS | k=3, g=-0.96, p=0, I2=20% (无效:样本双计) | k=16, g=-0.94 [-1.33,-0.56], 单侧p=5.3e-05, I2=85% | 旧值作废；新值LOSO稳定 |
| meta-pooled | EIS | k=3, g=+0.85, p=0.002, I2=86% (无效:样本双计) | k=16, g=+1.02 [+0.64,+1.40], 单侧p=1.8e-05, I2=82% | 旧值作废；新值LOSO稳定 |
| meta-pooled | MDI | k=3, g=+1.85, p=0.0048, I2=96% (无效:样本双计) | k=16, g=+1.50 [+1.03,+1.98], 单侧p=3.1e-06, I2=89% | 旧值作废；新值LOSO稳定 |
| dataset-effect | GSE185263 SepsisCOVID_vs_Control UCS | g=-1.23 (n=82v44) | g=-1.23 (n=82v44)  | 保留（样本级检验，单位本就正确） |
| dataset-effect | GSE185263 SepsisCOVID_vs_Control EIS | g=+1.36 (n=82v44) | g=+1.36 (n=82v44)  | 保留（样本级检验，单位本就正确） |
| dataset-effect | GSE185263 SepsisCOVID_vs_Control MDI | g=+3.44 (n=82v44) | g=+3.44 (n=82v44)  | 保留（样本级检验，单位本就正确） |
| dataset-effect | GSE32707 ARDSd0_vs_Control UCS | g=-0.93 (n=18v34) | g=-0.93 (n=18v34)  | 保留（样本级检验，单位本就正确） |
| dataset-effect | GSE32707 ARDSd0_vs_Control EIS | g=+0.02 (n=18v34) | g=+0.02 (n=18v34)  | 保留（样本级检验，单位本就正确） |
| dataset-effect | GSE32707 ARDSd0_vs_Control MDI | g=+0.77 (n=18v34) | g=+0.77 (n=18v34)  | 保留（样本级检验，单位本就正确） |
| dataset-effect | GSE212865 SDRA_vs_Control UCS | g=+0.23 (n=34v51) | g=+0.04 (n=19v51) 新值=患者级基线化（纵向重复剔除） | 重算（纵向去重后效应缩小） |
| dataset-effect | GSE212865 SDRA_vs_Control EIS | g=-0.78 (n=34v51) | g=-0.46 (n=19v51) 新值=患者级基线化（纵向重复剔除） | 重算（纵向去重后效应缩小） |
| dataset-effect | GSE212865 SDRA_vs_Control MDI | g=-0.61 (n=34v51) | g=-0.30 (n=19v51) 新值=患者级基线化（纵向重复剔除） | 重算（纵向去重后效应缩小） |
| scrna-arm | 旧：细胞级检验（把细胞当独立重复，显著性无效） | Table S2b 细胞级 Mann-Whitney（80基因） | 供者级 pseudobulk：28 检验中 2 个 BH q<0.05；最强：GSE158055 Mono_c14/UCS diff=-0.99 q=0.0045; GSE158055 NK/MDI diff=+2.14 q=0.00203 | 单位修正后：髓系UCS塌陷与NK-MDI存活；执行臂供者级未存活 |
| scrna-gene | 旧：无单细胞基因级 DE | （S2 实为 bulk DESeq2，N3 已改标签） | P0_pseudobulk_DE_all.csv：7 个数据集×细胞类型组合 pydeseq2（供者级） | 新增，供者级 |
| spatial | 切片/FOV 级检验 | M1 slide/FOV 级（单位缺陷） | 待 P0-7（patient–slide–FOV 映射 + 患者阻断） | 待办 |

## 说明

1. 旧 M2_meta（pooled g=1.85 等）因样本双计数无效，仅列作对照；
2. 直并行（GSE185263/GSE32707）数据集级效应数字基本保留（同为样本级检验，单位本就正确）；
3. GSE212865 旧值含纵向伪重复，新值为患者级基线化结果；
4. 单细胞旧检验单位错误，其显著性一律不引用；新供者级结果为唯一确证口径。
