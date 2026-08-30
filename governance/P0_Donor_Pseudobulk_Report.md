# P0-5 供者级 Pseudobulk 重算报告（H4）

- 日期：2026-08-20；冻结依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md（OSF https://osf.io/C7RYD/）
- 单位：cell type × donor；供者=有效 n；细胞级统计仅描述
- 主对比：COVID_severe vs Healthy（mild 不入确证）；分数据集独立检验
- 阈值（未注册参数，披露）：donor×celltype ≥50 细胞；双侧各 ≥3 供者

## 覆盖矩阵

| dataset | cell_type | n_donors | n_case | n_ctrl | eligible |
|---|---|---|---|---|---|
| GSE145926 | B_cell | 3 | 3 | 0 | 0 |
| GSE145926 | Club | 4 | 4 | 0 | 0 |
| GSE145926 | DC | 6 | 3 | 3 | 1 |
| GSE145926 | Epithelial | 7 | 6 | 1 | 0 |
| GSE145926 | Macrophage | 2 | 2 | 0 | 0 |
| GSE145926 | Mono_c14 | 9 | 6 | 3 | 1 |
| GSE145926 | NK | 6 | 4 | 2 | 0 |
| GSE145926 | T_cell | 9 | 6 | 3 | 1 |
| GSE158055 | B_cell | 24 | 19 | 5 | 1 |
| GSE158055 | Club | 0 | 0 | 0 | 0 |
| GSE158055 | DC | 1 | 1 | 0 | 0 |
| GSE158055 | Epithelial | 3 | 3 | 0 | 0 |
| GSE158055 | Macrophage | 0 | 0 | 0 | 0 |
| GSE158055 | Mono_c14 | 24 | 19 | 5 | 1 |
| GSE158055 | NK | 18 | 13 | 5 | 1 |
| GSE158055 | T_cell | 26 | 21 | 5 | 1 |

## 供者级臂评分检验

| dataset | cell_type | score | n_case | n_ctrl | diff | Welch p | BH q |
|---|---|---|---|---|---|---|---|
| GSE145926 | DC | UCS | 3 | 3 | +1.086 | 0.214 | 0.464 |
| GSE145926 | DC | EIS | 3 | 3 | +1.301 | 0.112 | 0.464 |
| GSE145926 | DC | MDI | 3 | 3 | +0.215 | 0.159 | 0.464 |
| GSE145926 | DC | MDI_nomt | 3 | 3 | +0.338 | 0.232 | 0.464 |
| GSE145926 | Mono_c14 | UCS | 6 | 3 | -0.775 | 0.37 | 0.531 |
| GSE145926 | Mono_c14 | EIS | 6 | 3 | +0.177 | 0.867 | 0.867 |
| GSE145926 | Mono_c14 | MDI | 6 | 3 | +0.952 | 0.177 | 0.464 |
| GSE145926 | Mono_c14 | MDI_nomt | 6 | 3 | +1.037 | 0.0477 | 0.464 |
| GSE145926 | T_cell | UCS | 6 | 3 | +0.864 | 0.328 | 0.531 |
| GSE145926 | T_cell | EIS | 6 | 3 | +0.244 | 0.736 | 0.803 |
| GSE145926 | T_cell | MDI | 6 | 3 | -0.620 | 0.532 | 0.638 |
| GSE145926 | T_cell | MDI_nomt | 6 | 3 | -0.979 | 0.399 | 0.531 |
| GSE158055 | B_cell | UCS | 19 | 5 | -0.670 | 0.0282 | 0.106 |
| GSE158055 | B_cell | EIS | 19 | 5 | -0.221 | 0.55 | 0.629 |
| GSE158055 | B_cell | MDI | 19 | 5 | +0.448 | 0.402 | 0.59 |
| GSE158055 | B_cell | MDI_nomt | 19 | 5 | -0.450 | 0.525 | 0.629 |
| GSE158055 | Mono_c14 | UCS | 19 | 5 | -0.993 | 0.000562 | 0.0045 |
| GSE158055 | Mono_c14 | EIS | 19 | 5 | -0.603 | 0.0543 | 0.145 |
| GSE158055 | Mono_c14 | MDI | 19 | 5 | +0.390 | 0.312 | 0.59 |
| GSE158055 | Mono_c14 | MDI_nomt | 19 | 5 | -0.367 | 0.339 | 0.59 |
| GSE158055 | NK | UCS | 13 | 5 | -0.758 | 0.0332 | 0.106 |
| GSE158055 | NK | EIS | 13 | 5 | +1.386 | 0.00951 | 0.0507 |
| GSE158055 | NK | MDI | 13 | 5 | +2.144 | 0.000127 | 0.00203 |
| GSE158055 | NK | MDI_nomt | 13 | 5 | -0.160 | 0.522 | 0.629 |
| GSE158055 | T_cell | UCS | 21 | 5 | -0.317 | 0.221 | 0.505 |
| GSE158055 | T_cell | EIS | 21 | 5 | +0.037 | 0.885 | 0.885 |
| GSE158055 | T_cell | MDI | 21 | 5 | +0.354 | 0.406 | 0.59 |
| GSE158055 | T_cell | MDI_nomt | 21 | 5 | -0.066 | 0.764 | 0.815 |

## 基因级（pydeseq2）

- 已运行：见 _intermediate/P0_pseudobulk_DE_all.csv（全基因）与 P0_pseudobulk_DE_80genes.csv

## 与注册文本的偏差披露

1. limma-voom duplicateCorrelation 敏感性未运行（需 R；本环境仅 Python）——注册的回退条款针对 DESeq2 拟合失败，此处为环境限制，如实披露，不冒充已运行。
2. 合并（跨数据集 ~dataset+condition）敏感性本轮未跑（分数据集为主分析已覆盖注册主问题）。
3. MIN_CELLS/MIN_DONORS 阈值为未注册参数，取常规值并在报告与代码中披露。
