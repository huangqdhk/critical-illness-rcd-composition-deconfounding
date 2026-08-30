# P0-8 80 基因框架用途更名交付文档

- 日期：2026-08-21；依据：`P0_FROZEN_ANALYSIS_PLAN_v1.0.md` §7（更名仅改用途描述，不动成员）与《创新提质方案 v2》Phase 0 任务 3
- 处置对象：`04_AUDIT_GOVERNANCE/Mitoxyperilysis_Gene_Manifest_v1.0.csv`（80 行 × 8 模块，成员与数值**零改动**）
- 性质：用途与术语归属声明（P0-8）；不改变任何统计结果、基因集成员、SHA256 校验对象

## 1. 更名决定

| 项目 | 旧用途（2026-08-19 前） | 新用途（自本文件起） |
|---|---|---|
| 框架名称 | Mitoxyperilysis canonical gene set（Mitoxyperilysis 规范基因集） | **80 基因线粒体基础设施–炎症执行解离框架**（80-gene mitochondrial-infrastructure / inflammatory-execution framework） |
| 框架用途 | 用于"界定 Mitoxyperilysis" | **描述危重症转录组中"上游基础设施塌陷（UCS）与炎症死亡执行诱导（EIS）解离"的组织层状态**；非任何已命名死亡方式的规范基因集 |
| 术语归属 | 本文命名/界定 | Mitoxyperilysis 术语与机制定义唯一归属于 Wang et al. *Cell* 2025（PMID 41317732）；本框架仅登记与之**表型对齐**的候选组织层状态的可操作测量（UCS/EIS/MDI） |
| 机制对应 | MAM-铁节点"界定"线粒体-质膜接触死亡 | 本框架的 MAM/线粒体模块属上游基础设施读出；**不构成**线粒体-质膜接触或任何亚细胞形态的直接证据 |

## 2. 文件层面处置（零数值改动）

1. CSV 文件名保持 `Mitoxyperilysis_Gene_Manifest_v1.0.csv` 不变——改名将破坏全部读取脚本与 RESULTS_MANIFEST SHA256 审计链，收益为零；用途纠偏以本文件+文稿 v4 术语口径为准（与 RESULTS_MANIFEST §二"符号策略决策"同逻辑：审计链优先）。
2. CSV 内 `gene_set_version=Mitoxy-80_v1.0` 保持为**版本标识符**（字符串），不解读为机制归属声明；本文件即该标识符的用途释义（登记于 RESULTS_MANIFEST 例外表）。
3. 模块名 `mitoxy_*` 前缀同理：历史标识符，自 v4 起所有文稿叙事已不再将其写作"Mitoxyperilysis 模块"，统一写作"上游臂/执行臂/调节模块"。
4. 任一依赖本框架的新增分析（含 P1/P2 下游）继续引用同一 manifest 文件与 `Mitoxy-80_v1.0` 版本戳，保证可复现链不断。

## 3. 术语使用禁令（同步登记）

- 不得写："我们命名 Mitoxyperilysis""本 80 基因集界定 Mitoxyperilysis""Mitoxyperilysis canonical gene set"。
- 不得把 UCS/EIS/MDI 评分写作死亡方式本身的证据；其身份是"候选组织层状态的可操作测量"。
- Deferasirox 相关表述仅限"铁稳态相关探索候选"。
- 与 Wang et al. 2025 机制对应的必要节点（BAX/BAK1/BID、RICTOR/mTORC2、RHOA、NINJ1、actin/lamellipodia）不在本框架内——引用本框架时必须披露该覆盖缺口，不得暗示覆盖完整。

## 4. 审计与验收

- 本文件为纯用途声明，不改 CSV → RESULTS_MANIFEST 与 lint 校验对象不变；`lint_package.py` 无需重跑（无文件内容变更）。
- 文稿 v4 已按此口径改写（标题、Highlights、Abstract、§4.19、Results §11、Discussion §12）；本文件补足治理层留痕（冻结计划 Phase 0 任务 3 的正式交付物）。
