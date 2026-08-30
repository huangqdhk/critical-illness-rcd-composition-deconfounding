# _intermediate/ — WGCNA 中间产物（非结果表）

2026-08-16 包级清理（lint 全绿恢复）时自 `01_RESULTS_TABLES/` 顶层移出。

- 内容：Bulk WGCNA 的 5 个中间文件（ME_df / ME_var / SFT_curve / expr_block / metadata）
- 性质：`WGCNA_v2_rerun.py` / `Mitoxyperilysis_ARDS_GSE185263_WGCNA_v2.py` 的过程产物，
  非 RESULTS_MANIFEST 登记对象、不进入投稿包
- 重生成方式：重跑上述脚本即可（输入为 GSE185263 表达矩阵与 `GSE185263_groups.csv`）
