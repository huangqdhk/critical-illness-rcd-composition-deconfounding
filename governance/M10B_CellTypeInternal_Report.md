# M10 检验 B 报告:细胞类型内部供者级检验(髓系内 UCS)

- 日期:2026-08-26;依据:M10_M13_M14_pre_registration_20260824.md §1.3(osf.io/ETVMJ)
- 计算留痕:`03_LOGS/M10B_cross_cohort_log.txt`、`03_LOGS/M10B_gse216009_log.txt`、`03_LOGS/M10B_gse180578_log.txt`
- 方法(沿用 ETVMJ §1.3 / P0 口径):cell type × donor pseudobulk(donor×celltype ≥50 细胞、双侧各 ≥3 供者);臂评分 = 供者级 log1p 均值 → celltype 内 z → Welch t;臂基因 = 上游 30 / 执行 33(IIAMD v1.0);BH 族 = 数据集内 celltype × score
- 队列:GSE158055(PBMC scRNA,已持有)、GSE145926(BALF scRNA,已持有)、GSE216009(全血 Rhapsody,272,993 细胞,87 样本)、GSE180578(COVID PBMC/TA,h5ad 2,000 基因矩阵不用,pseudobulk 用 RAW 全基因组 33,538 基因;去后缀 barcode 碰撞 8.1% 歧义细胞剔除)——共 4 个独立 scRNA 队列

## 1. 逐队列髓系(单核/巨噬)内部 UCS 层效应

| 队列 | 层 | n(病 vs 对照) | UCS Hedges' g | 备注 |
|---|---|---|---|---|
| GSE158055 | Mono_c14 | 19 vs 5 | **−1.030**(se 0.524) | PBMC |
| GSE145926 | Mono_c14 | 6 vs 3 | **−0.699**(se 0.726) | BALF;Macrophage 层 0v0 不足跳过 |
| GSE216009 | Classical_monocytes | 23 vs 13 | **−0.780**(se 0.359) | 全血;Welch diff=−0.753,p=0.0201,q=0.142 |
| GSE180578 | Classical monocytes | 21 vs 10 | **−1.399**(se 0.423) | Welch diff=−1.200,p=0.000154,q=0.00103 |
| GSE216009 | Non-classical_monocytes | 5 vs 4 | −0.188 | 层薄,仅描述 |
| GSE216009 | Mature_neutrophils | 26 vs 13 | **+0.322** | 中性粒内部 UCS 不降(见 §4) |
| GSE216009 | S100A8-9_hi_neutrophils | 26 vs 13 | **+0.520** | 同上 |
| GSE180578 | Intermediate/Non-classical | 4v0 / 2v0 | 跳过 | 供者不足 |

## 2. 合并(每队列单核/巨噬代表层,k=4)

- **合并:g = −1.003,95%CI(HK)= [−1.510, −0.495],单侧 p(neg) = 0.0041,I² = 0%,预测区间 [−1.993, −0.012]**
- 方向一致(g<0):**4/4(100%)**;LOSO 同向 **4/4**

## 3. 敏感性:单核+中性粒全髓系层(k=7)

- 合并 g = −0.419 [−1.113, +0.275],单侧 p = 0.095,方向一致 71%(并入中性粒层后减弱,如实记录)

## 4. 附注(与 M10 检验 A 图景的一致性)

中性粒层内部 UCS 为正(+0.32/+0.52)——即全血 UCS 下降并非在粒细胞内部重现,而单核/巨噬内部存在一致的细胞内在下降。这与检验 A 的结论(组成解释全血信号的大部分,残差中保留跨队列一致的细胞内效应)构成互补证据链:支持 H_int(细胞内在成分)而非纯 H_comp。

## 5. 判定

**检验 B 成立:≥3 个独立 scRNA 队列中髓系(单核/巨噬)内部 UCS 一致下降(4/4 队列,合并显著,LOSO 稳定)→ 支持 H_int。**

与 M10 检验 A 的 A-PASS(见 `M10A_Composition_Test_Report.md` 第 45 行)联合,按 ETVMJ §4 定义:**Track B 判定成立(2026-08-26)** → 触发 §5 条件条款:M11/M12 启动并另立注册(即 2026-08-27 的 M11/M12 前瞻注册)。
