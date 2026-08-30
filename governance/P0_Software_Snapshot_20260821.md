# 软件环境快照（2026-08-21 采集）

> 依据：《创新提质方案 v2》§8.3『固化 pip freeze、R sessionInfo()、随机种子和原始数据下载校验和』。
> 采集方式：2026-08-21 在本机从项目两个 venv 分别执行 `pip freeze`；R 环境不在本机，见 §3 处置。

## 1. Python 环境快照

| 环境 | 用途 | 快照文件 | 包数 |
|---|---|---|---|
| `.venv_spatial` | M1 空间/P0-P2 重算（scanpy、libpysal、esda、sklearn、pydeseq2） | `04_AUDIT_GOVERNANCE/pip_freeze_spatial_venv_20260821.txt` | 65 |
| `.venv_pyaging` | pyaging/表观时钟/GEM/COBRA 层 | `04_AUDIT_GOVERNANCE/pip_freeze_pyaging_venv_20260821.txt` | 167 |

关键版本（从快照摘录，投稿时可引用）：

| 包 | 版本 |
|---|---|
| numpy | 2.5.2（spatial）/ 2.5.1（pyaging） |
| scipy | 1.18.0 |
| pandas | 见快照（v4 文稿声明 3.0.5，以快照文件为准） |
| scikit-learn | 见快照（v4 声明 1.9.0，以快照文件为准） |
| pydeseq2 | 0.5.4（P1/P2 DESeq2 实现，P1 报告已登记） |
| libpysal / esda | 空间自相关（Moran's I 实现） |

## 2. 随机种子

- 全项目空间置换与抽样统一 `seed=0`（M1 脚本、P0-7b 患者阻断置换、P1-4 匹配基因集检验均使用 seed=0/999 次或 500 次置换）；P1 LORO 为确定性重拟合，无随机成分。
- 其余历史分析的种子以各 canonical 脚本内 `seed=` 为准（已在 lint/方法注登记）。

## 3. R 环境（sessionInfo）

- **本机无 R**（`Rscript` 不可用）；M2/M3 的 R 步骤（DESeq2、coloc 5.2.3、TwoSampleMR、TWAS FUSION 读取）在完成时未留存 `sessionInfo()`。
- 处置（如实登记，不伪造）：投稿前须在当初运行 R 步骤的机器上执行 `Rscript -e "sessionInfo()"` 并保存；R 包版本声明以那时输出为准。现有 R 中间产物（M3_coloc_*.csv、M3_smr_*、M3_twas_*、M2_*）已在 RESULTS_MANIFEST v2.0 冻结，版本快照缺失不影响结果可复现性审计，但影响软件可重复性声明完整性。

## 4. 原始数据下载校验和

- 结果表校验和：`RESULTS_MANIFEST_v1.0.csv` / `v2.0.csv`（SHA256，lint 全绿，2026-08-16 后 247 个 canonical 文件）。
- 原始下载文件校验和：`00_RAW_DATA/README.md` 与 `SAMPLE_MANIFEST_v1.0.csv` 已登记 accession/来源；逐文件 SHA256 清单未单独生成（投稿前如期刊要求，可对 `00_RAW_DATA/` 生成目录级 SHA256 清单）。
