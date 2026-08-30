# 分析计划执行完成状态总报告（2026-08-21）

> 依据：《111黄裕荣创新提质_v2.md》（2026-08-19 重塑版）逐项执行审计
> 审计人：分析作者 + AI 辅助执行；数据原则：全部结果可溯源、阴性不改口径

## 一、总判

**Phase 0 / Phase 1 / Phase 2 的计算任务已全部完成（含本轮补做的 5 项残留）；Phase 3 与 Phase 4 为湿实验阶段，无法以计算完成，须实验室资源投入——如实登记，未作任何替代性"完成"声明。**

计划 §8.3 投稿前硬阻断项中**可计算部分已全部执行**；两项需要作者本人操作：OSF 解除 embargo 公开（需登录 OSF 项目页手动 Make public）、作者/单位/伦理批件信息填写。

## 二、Phase 0（构念与统计急救）——完成

| # | 计划任务 | 状态 | 交付物 |
|---|---|---|---|
| 1 | 查新表 | ✅ | `P0_Novelty_Search_Table_v1.0.md`（08-20）+ **v1.1（08-21 PubMed 重跑更新）** + `P0_Novelty_Search_Rerun_Log_20260821.md` |
| 2 | Wang 2025 为唯一一级来源 | ✅ | 文稿 v4/v4.1 全文改写 |
| 3 | 80 基因 manifest 更名（不再称 canonical gene set） | ✅（本轮补做 P0-8） | `P0-8_Framework_Rename_and_Usage_Statement.md`（成员与数值零改动） |
| 4 | 供者级 pseudobulk + 患者阻断空间 + 来源互斥 meta | ✅ | P0-5 报告（08-20）+ **P0-7b 患者级置换（08-21 补做，置换 p<0.001）**；P0-4 报告 |
| 5 | 原结果 vs 正确统计单位对照表 | ✅ | `P0_Old_vs_New_Comparison.md` |
| 6 | OSF 冻结 | ✅（注册）/ ⚠️ 待作者公开 | osf.io/C7RYD，embargo 2030-06-30；**须作者手动 Make public** |

Go/No-Go：完整双臂在正确统计单位下存活（k=16 来源互斥 meta UCS/EIS/MDI 全部 LOSO 稳定）→ 按计划进入 Phase 1（已执行）。

## 三、Phase 1（实验锚定签名）——完成，Gate 1 通过

- GSE235046 2×2 交互模型：IIAMD v1.0 冻结（up 4,367 / down 4,720）；
- LORO 12 次：方向保持 100%、FDR 保持 93.0–95.6% → STABLE；
- **随机匹配基因集检验（本轮补做 P1-4，此前缺失）**：baseMean 十分位匹配随机集 500 次，mean|β₃|、|β₃|>1 计数、12 机制锚定节点命中三项统计量均显著超随机（单侧 p=0.002；锚定节点 up 臂命中 7/12）；
- ST002738 代谢组锚定：GSH 交互项全负，与原文"GSH 耗竭"方向一致（描述性）。
- **Gate 1：通过**（稳定交互签名 + 代谢方向一致）。

## 四、Phase 2（人类迁移验证）——完成，Gate 2 未通过（如实判定）

三层头对头：真 ARDS（GSE32707）IIAMD_up q=0.123 不显著且劣于 hypoxia/mito-stress 面板；发现层关联被髓系组成完全解释（调整后 p=0.499）；单细胞供者级无髓系定位。**按冻结条款撤回"Mitoxyperilysis-aligned program in ARDS"疾病映射**，唯一三层存活信号为线粒体基础设施抑制——与计划 §十二回退分支一致（描述性档位）。

## 五、本轮补做的 Phase 0/1 残留项（5 项）

1. **P0-8 manifest 用途更名**（计划任务 3 的正式交付物）；
2. **P1-4 随机匹配基因集检验**（计划 Phase 1 冻结三资产任务的第二项，此前只有 LORO）；
3. **P0-7b CosMx 患者级置换重算**（P0-7 报告"列于后续"的遗留项：18 Case 等权 Moran's I 均值 +0.0228，患者阻断置换 999 次 p<0.001，FOV 级结论在患者级成立）；
4. **查新重跑**（§8.3 硬阻断）：Q1=11、**Q2(Mitoxyperilysis×ARDS/肺)=0**、Q3(脓毒症)=1（评述）；2026 年新增 2 篇肿瘤"相关签名"论文与多篇综述已登记；**更正 v1.0 一处登记错误**（PMID 41856854 期刊应为 Trends Biochem Sci，非 Trends Pharmacol Sci）；
5. **空间层文稿升级**：CosMx 结论由"按注册回退条款降为描述性"升级为"患者级统计单位下成立（病例内描述）"。

## 六、§8.3 投稿前硬阻断项逐条

| 项 | 状态 |
|---|---|
| OSF embargo 解除并公开 | ⏳ **作者已确认将在投稿前 1–2 天手动 Make public**（无 OSF 凭据，需作者操作） |
| 重跑查新并更新查新表 | ✅（Q2=0 空白维持；v1.1 + 检索日志） |
| 补齐参考文献 | ✅（1-64 按正文引用上下文重建 + PubMed 逐条核实；65-77 全部核实 PMID；v1 原清单不在项目内，无法唯一锚定的 8 条标注【待作者确认】并附候选，绝不虚构） |
| 作者/单位/伦理 | ⚠️ 待作者填写（【TBD】占位保留） |
| 英文全文稿 | ✅ `111_manuscript_v4_EN.md`（完整翻译，数字与中文源完全一致，注明须领域专家核校） |
| 最终主图导出 | ✅ 新 P0/P1/P2 层 4 张 PNG（此前该层无图件）；既有 Figure 1-10 面板数据已在 `01_FIGURE_DATA_CSV/` |
| pip freeze / sessionInfo | ✅ pip freeze×2 快照；R sessionInfo 本机无 R，缺失已如实登记（`P0_Software_Snapshot_20260821.md`） |
| 代码公开 | ✅ 全部新脚本落盘项目根目录并可重跑 |
| manuscript-to-source 自动核对 | ✅ **17/17 PASS**（`P0_Manuscript_Source_Check_Report.md`） |
| lint 全绿 | ✅ **74 项通过 / 0 失败**（本轮从 194 失败修复：版本戳漂移 164 行、N22 编号归一、archive 子目录重组、P0-P2 系列登记、manifest v2.0 重建为 591 行） |

## 七、Phase 3 / Phase 4 —— 湿实验阶段，未执行（如实声明）

- **Phase 3（人源细胞最小机制闭环）**：需要 ≥3-5 名独立供者的人外周血单核细胞来源巨噬细胞、活细胞成像（线粒体-质膜接触时序）、MitoSOX/CellROX/4-HNE 局部氧化、GSH/TMRM/Seahorse、LDH/HMGB1/Sytox 破膜、p-AKT-S473/RICTOR/RhoA-GTP/F-actin 机制节点与救援矩阵（NAC/Torin-1/RICTOR KD/rapamycin/BAX-BAK1-BID 抑制/cytochalasin D + z-VAD/Nec-1s/Fer-1/MCC950 排除对照）。**属实验室工作，本计算执行无法完成。**
- **Phase 4（患者组织与动物因果闭环）**：需要独立 ARDS/ICU 非 ARDS 肺组织多重成像、TEM 定量线粒体-质膜距离、髓系特异 Rictor 敲除动物模型与肺损伤终点。**同上，属实验室/动物实验工作。**
- 计划给出的可操作下一步：Gate 2 已 FAIL 的事实下，Phase 3 的 2×2 设计与救援矩阵是**唯一能改变稿件档位的路径**；若实验室资源不可得，稿件按"描述性档位"（IF 4-8 区间）提交，结论维持 v4.1 现状（不得写"Mitoxyperilysis participates in ARDS"）。

## 八、数字与治理验证摘要

- 来源互斥 meta（k=16）：UCS g=−0.94、EIS g=+1.02、MDI g=+1.50，LOSO 全稳——与文稿 §11 逐位一致（17 项核对全 PASS）；
- 供者级 pseudobulk：28 项仅 2 项存活（GSE158055 单核 UCS q=0.0045；NK MDI 为 MT 驱动已披露）——与文稿一致；
- Gate 1 STABLE、Gate 2 FAIL——与文稿一致；
- lint 全绿；RESULTS_MANIFEST v2.0（591 行）v2.1 变更已登记于 RESULTS_MANIFEST.md 头部。
