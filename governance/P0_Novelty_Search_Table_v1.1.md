# P0-2 查新表：Mitoxyperilysis 文献登记 v1.1（投稿前重跑更新版）

> 检索日期：2026-08-21（重跑）；检索人：分析作者 + AI 辅助检索
> 检索接口：NCBI E-utilities（PubMed）；检索式与命中存档见 `P0_Novelty_Search_Rerun_Log_20260821.md`
> 基线：v1.0（2026-08-20，Web 聚合检索）；本版以 PubMed `[tiab]` 检索结果逐条核对后更新。
>
> **检索式（预注册口径）：**
> - Q1：`"mitoxyperilysis"[tiab] OR "mitoxyperiosis"[tiab]` → **11 条**
> - Q2：Q1 AND `("acute respiratory distress"[tiab] OR ARDS[tiab] OR lung[tiab])` → **0 条**
> - Q3：Q1 AND `sepsis[tiab]` → **1 条**（PMID 41856854，评述性展望，无数据）

## 核心结论（v1.1 更新）

**截至 2026-08-21，仍无任何以直接实验证据将 Mitoxyperilysis 与 ARDS/急性肺损伤联系的研究（0 条）。**
但相比 v1.0 出现两个必须登记的新动向：

1. **"相关签名"类论文已出现 2 篇**（均为肿瘤）：PMID 42491158（结直肠癌预后签名，Frontiers）、PMID 42170303（黑色素瘤单细胞+机器学习签名，Human Mutation）——按《创新提质方案》§二 Phase 0 任务 1/2 口径，此类论文只能作为**二级方法参考**，其基因集不得被误认为 Wang et al. 2025 的机制定义。
2. **v1.0 一处登记错误更正**：PMID 41856854 的期刊为 **Trends in Biochemical Sciences**（DOI 10.1016/j.tibs.2026.01.005），v1.0 误写为 Trends Pharmacol Sci 51(4):313-315；卷期页码亦相应更正为 TIBS 口径（卷期待出版信息以 PubMed 为准）。

## 逐条登记（v1.1 全量，PubMed 核实）

| # | 文献 | 类型 | 疾病语境 | 基因集来源 | 直接实验级别 | 链接 |
|---|---|---|---|---|---|---|
| 1 | Wang Y, Lu J, Carisey AF, et al. Innate immune and metabolic signals induce mitochondria-dependent membrane lysis via mitoxyperiosis. **Cell**. 2025;188(25):7155-7174.e25. PMID 41317732 | **原始研究（一级来源，术语与机制唯一权威定义）** | 体外 BMDM/肿瘤（黑色素瘤）；无 ARDS/肺 | 原文自有（BAX/BAK1/BID、RICTOR/mTORC2、RhoA-actin、NINJ1） | 细胞+遗传干预+活细胞成像+动物 | [PubMed](https://pubmed.ncbi.nlm.nih.gov/41317732/) |
| 2 | Casadio M. Mitoxyperilysis as a distinct cell death type. **Nat Cell Biol**. 2026. PMID 41540269 | 评论/亮点（二级） | 无疾病 | 引用 #1 | 无 | [PubMed](https://pubmed.ncbi.nlm.nih.gov/41540269/) |
| 3 | Mitoxyperilysis: fasting-induced cell death in immunometabolism and disease. **Trends Biochem Sci**. 2026. PMID 41856854（v1.0 误登记为 Trends Pharmacol Sci 51(4):313-315，已更正） | 评论（二级） | **脓毒症/癌症治疗展望（无数据）** | 引用 #1 | 无 | [PubMed](https://pubmed.ncbi.nlm.nih.gov/41856854/) |
| 4 | Mitochondria at the membrane provide a route to inflammatory cell death. **Cell Metab**. 2026. PMID 41638191 | 评论（二级） | 无疾病 | 引用 #1 | 无 | [PubMed](https://pubmed.ncbi.nlm.nih.gov/41638191/) |
| 5 | Mitochondria target the plasma membrane to cause mitoxyperiosis. **Cell Res**. 2026. PMID 41530348 | 评论（二级，新登记） | 无疾病 | 引用 #1 | 无 | [PubMed](https://pubmed.ncbi.nlm.nih.gov/41530348/) |
| 6 | Spatially gated oxidative killing: mitoxyperilysis redefines how ROS cause lytic cell death. **Apoptosis**. 2026. PMID 41721119 | 综述（二级，新登记） | 无疾病 | 引用 #1 | 无 | [PubMed](https://pubmed.ncbi.nlm.nih.gov/41721119/) |
| 7 | Mitoxyperilysis links organelle positioning to inflammatory lysis: A new lens for pharmacotherapy design. **Biomed Pharmacother**. 2026. PMID 41653904 | 综述（二级，新登记） | 无疾病 | 引用 #1 | 无 | [PubMed](https://pubmed.ncbi.nlm.nih.gov/41653904/) |
| 8 | Mitoxyperilysis: Rethinking oxidative stress as a spatially constrained lethal signal. **Int J Biol Sci**. 2026. PMID 42088416 | 综述（二级，新登记） | 无疾病 | 引用 #1 | 无 | [PubMed](https://pubmed.ncbi.nlm.nih.gov/42088416/) |
| 9 | Mitoxyperiosis: A Novel Mitochondria-Driven Cell Death Mechanism in Immunometabolic Stress. **Cell Biol Int**. 2026. PMID 42605933 | 综述（二级，新登记） | 无疾病 | 引用 #1 | 无 | [PubMed](https://pubmed.ncbi.nlm.nih.gov/42605933/) |
| 10 | A mitoxyperilysis-related signature stratifies prognosis and identifies an aggressive colorectal cancer ecosystem with immune remodeling. **Front Cell Dev Biol**. 2026. PMID 42491158 | **"相关签名"研究（二级方法参考，新登记）** | **结直肠癌（肿瘤预后签名）** | **自建"related signature"基因集（非 Wang 2025 机制清单）** | 队列转录组计算，无机制实验 | [PubMed](https://pubmed.ncbi.nlm.nih.gov/42491158/) |
| 11 | A Mitoxyperilysis-Related Single-Cell and Machine-Learning Framework Defines an Immune-Cold Melanoma Phenotype and a Robust Prognostic Signature. **Hum Mutat**. 2026. PMID 42170303 | **"相关签名"研究（二级方法参考，新登记）** | **黑色素瘤（肿瘤预后签名）** | **自建"related signature"基因集（非 Wang 2025 机制清单）** | 队列转录组计算，无机制实验 | [PubMed](https://pubmed.ncbi.nlm.nih.gov/42170303/) |
| 12 | Zenodo 预印本：对 Wang 2025 的批判性再评估 | 预印本评论（未同行评审，非 PubMed 索引） | 无疾病 | — | 无 | [Zenodo](https://zenodo.org/records/17797711) |
| 13 | 生物谷/《生物化学与生物物理进展》/知乎-MCE 等中文科普解读 | 科普/解读（三级，非 PubMed 索引） | 肿瘤免疫 | 引用 #1 | 无 | [生物谷](https://news.bioon.com/article/7c72913999aa.html) |

## 对本项目的含义（v1.1）

1. **优先权与空白**：Mitoxyperilysis×ARDS/肺 仍为 0 条；本项目"人 ARDS 中检验"的定位空白成立。但经 Gate 2（2026-08-20 判定 FAIL）后，本项目对该空白的贡献为**阴性边界**（实验锚定签名未迁移）而非确证。
2. **领域动向（相比 v1.0 显著变化）**：2026 年内新增 2 篇肿瘤"Mitoxyperilysis-related signature"论文与 5 篇综述/评论——领域正在快速靠近且已出现签名化滥用苗头。本项目的差异化定位（直接机制检验、阴性结果纪律、来源互斥统计）反而因该动向而更具必要性；投稿不应再拖延。
3. **引用更新义务**：v1.0 表第 3 条期刊名登记错误已在 v1.1 更正；文稿 References 中对应条目须同步更正（Trends Pharmacol Sci → Trends Biochem Sci）。
