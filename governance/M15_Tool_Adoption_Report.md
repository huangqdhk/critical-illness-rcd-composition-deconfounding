# M15 工具封装与采纳演示治理报告（2026-08-27）

- 依据：111黄裕荣创新提质_20260822.md §五 M15（工具核心 = 两臂分解 + MDI + 组成预测值/组成残差 MDI；GitHub + Zenodo + 判定门文档化；2 个第三方数据集采纳演示）

## 1. 工具包

- 位置：M15_tool/（mitoxdi v1.0.0；score_version = mdi_v1.0；gene_set_version = Mitoxy-80_v1.0）
- 组成：mdi_core（两臂分解+评分，mdi_lib 冻结快照）/ composition（Monaco+ABIS 参考谱、NNLS 主口径 + OLS-CLS 交叉、合成谱评分、残差）/ report（对比统计 + ISED 解读）/ loaders（逐队列口径复刻）
- 判定门语言：M15_tool/GATE.md（M10 检验 A 判定门的通用化表述）
- 数据依赖（零下载）：Monaco 冻结缓存 _intermediate/M10A_monaco_ct_means.csv；ABIS sigmatrixRNAseq/sigmatrixMicro（00_RAW_DATA/ABIS/）；80 基因 manifest（内置副本与治理目录 SHA256 一致 346537be…e9965）；GPL570/GPL23159/GPL6947/GPL10295 探针映射（00_RAW_DATA）
- 发布件：README/LICENSE/requirements/CITATION.cff；GitHub + Zenodo DOI 随投稿前打包定稿（代码与文档，不含数据下载）

## 2. 内部回归测试（11 队列）

- 参照冻结表：M2（6 队列评分）、M10A（组成预测/残差）、M11M12（5 队列评分+组成+残差+比例）
- 结果：**254 族逐样本核对全部 PASS（254/254；队列级 17/17）**；族级 max|Δ| 最大值为 3.770×10⁻¹³（机器精度量级）；封装评分核心与项目共享 mdi_lib 同一矩阵 max|Δ|=0
- 表：Table S84/S84b；图 10C 面板数据；日志 03_LOGS/M15_regression_log.txt

## 3. 采纳演示（2 个第三方数据集，GEO 官方核验入列）

- **GSE157103**：n=126（COVID 100/非 COVID 26）；臂覆盖 up 24/30 ex 33/33；MDI R²(Monaco|NNLS)=0.809（OLS-CLS 0.475；ABIS 0.286）；观测 MDI g=−0.212（p=0.35，无疾病升高；UCS g=+0.751、EIS g=+0.575 两臂同升）→ 残差 MDI g=+0.322（p=0.094，方向翻转、衰减比例不适用）
- **GSE66099**：n=276（脓毒性休克 181/脓毒症 18/SIRS 30/对照 47）；臂覆盖 up 27/30 ex 32/33；MDI R²(Monaco|NNLS)=0.817（OLS-CLS 0.584；ABIS Micro 0.238 覆盖受限）；脓毒性休克观测 MDI g=+2.172（p=1.7×10⁻²²）→ 残差 MDI g=+0.261（p=0.105，**衰减 88.0%**）；中性粒比例 0.200→0.340、单核 0.102→0.126
- 核验留痕：03_LOGS/M15_demo_GEO_verification_20260827.md；报告：03_LOGS/M15_demo_*_report.md
- 解读纪律：仅陈述客观量；按 GATE.md 判定门语言解读；不评价原论文结论

## 4. 导出登记

- 主图面板：Figure_10A–10C（量化表，不出图——按项目'优先生成量化数据表格'输出规则）
- 补充表：Table S84/S84b、S85/S85b–S85d、S86/S86b–S86d（10 个文件）已写入 RESULTS_MANIFEST_v2.0.csv（sha256/行数/字节数）
- lint_package.py M1X_NO_VERSION_STAMP 已补登（13 个新文件名）

## 5. 最终状态（2026-08-27）

- 回归测试：11 队列 × 254 族核对 254/254 PASS（队列级 17/17；max|Δ|=3.77×10⁻¹³）
- 导出登记：Figure_10A–10C + Table S84–S86 全部登记入 RESULTS_MANIFEST_v2.0.csv（sha256/行数/字节数）
- 治理披露：M15_step4 重写 manifest 时因 pd.read_csv 默认 na 解析把字面量 'NA' 版本戳解析为 NaN 造成损伤，已由 M15_step4b_manifest_repair.py 按 lint 同口径从文件重建并修复 Figure_S12J 多值戳（留痕 03_LOGS/M15_manifest_repair_log.txt；M15_step4 已改以 keep_default_na=False 读取）
- lint 终态：**74/0 全绿**（03_LOGS/lint_report.md）
