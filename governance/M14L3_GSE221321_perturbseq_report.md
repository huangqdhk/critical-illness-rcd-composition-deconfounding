# M14 腿3 外部 Perturb-seq（GSE221321）GO 判定与三角化分析报告

> 生成日期：2026-08-26
> 前瞻注册：osf.io/ETVMJ（`M10_M13_M14_pre_registration_20260824.md` §3.3/§3.4），2026-08-24 冻结
> 主方案：`111黄裕荣创新提质_20260822.md` §四 M14
> 脚本：`M14_step2_perturbseq_gse221321.py`（项目根目录；seed=0 全流程可复现）
> 数据：`00_RAW_DATA/GSE221321_THP1_PerturbSeq/GSE221321_RAW/`（作者解压，18 文件；tar 入库校验 2026-08-26，18 条目，退出码 0）
> 输出表：`02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV/Table_S72–S75_M14L3_*.csv`

## 1. GO 判据核验（注册 §3.3 四项，全部 PASS）

| 判据 | 证据 | 结论 |
|---|---|---|
| ① 完整单细胞计数矩阵 | GSM6858447 h5ad 86,956 细胞 × 16,952 基因；GSM6858449 h5ad 66,283 细胞 × 18,017 基因（X 为作者 log2 归一化稀疏矩阵；gzip 解压自作者 RAW.tar） | PASS |
| ② gRNA 映射表 | 各 GSM `*_perturbations.txt.gz`：600 行（598 靶基因 + 2 NTC）× 细胞列（0–3 条 guide/细胞）；obs.Guides_collapsed_by_gene 与之口径一致 | PASS |
| ③ 非靶向对照（NTC） | `non-targeting` + `safe-targeting` 两行存在于 perturbations 表与 FRPerturb 表；KO 8,471 / KD 5,027 个纯 NTC 细胞 | PASS |
| ④ 髓系/LPS 刺激背景 | 人 THP-1 单核细胞系 + LPS 3h（Yao et al., bioRxiv 2023, PMID 36747806，Compressed Perturb-seq）；GEO 编号 2026-08-24 官方核验锁定，2026-08-26 复核 PMID | PASS |

**判定：腿3 = GO，进入分析。**（NO-GO 条款未触发）

## 2. 预指定分析口径（先于结果，脚本头部冻结）

- 评分：mdi_v1.0（臂均值 → 数据集内 z，全细胞；MDI = z(EIS) − z(UCS)；臂定义单一来源 = Mitoxyperilysis_Gene_Manifest_v1.0.csv）
- 细胞分类：NTC-only / 单扰动（恰 1 靶基因，主分析）/ 复合扰动（排除，计数披露）
- 逐基因效应：Δscore = mean(扰动细胞) − mean(NTC 细胞)；n_g ≥ 20；Cliff's δ + MWU p + BH q；NTC 匹配置换双侧 p（臂成员 B=10,000，其余 B=2,000，seed=0）
- 效应量级：官方 FRPerturb（Perturbed × Downstream LFC）投影 ΔUCS_FRP/ΔEIS_FRP/ΔMDI_FRP；臂成员集合检验 vs 表达量匹配随机等大小基因集（B=10,000，seed=0，NTC 均值表达十分位匹配），单侧（自抑方向）+ 双侧
- 负对照校准：NTC split-half（B=10,000）零中心检查；FRPerturb NTC 行投影 ≈0 检查
- 交叉复核：KO vs KD Spearman + 臂成员方向一致率
- 敏感性：MDI_nomt；臂成员 leave-self-out（排除被扰动基因自身后重算臂均值）
- 统计单位：逐基因效应为细胞级描述量；推断检验单位为基因集（置换在基因层面）

## 3. 结果

### 3.1 负对照校准（全部通过）
- NTC split-half 零分布中心 |mean| ≤ 0.0007（8 组，B=10,000）——估计器无系统偏倚
- FRPerturb NTC 行投影：|ΔUCS_FRP| ≤ 0.0023、|ΔEIS_FRP| ≤ 0.0102、|ΔMDI_FRP| ≤ 0.0124（≈0）

### 3.2 主检验（全阴）
- 扰动库内臂成员：7/63（UCS 臂 HSPD1；EIS 臂 CASP1/CASP4/GPX4/GSDMC/IL1B/NLRP3）
- 臂内自抑集合检验：两模态 8 项单侧 p = 0.20–0.94；最佳双侧 p = 0.081（KO UCS_arm_mdi）——无一显著

### 3.3 描述性结果（不入判定）
- HSPD1 KD → UCS −0.364（NTC 置换 p=0.0009；leave-self-out −0.248 方向保持）；KO 同向 −0.231（p=0.21，n=29）
- EIS 成员方向不一致（IL1B KD → EIS −0.241；GSDMC KD → EIS +0.322）
- 逐细胞 ρ(UCS,EIS) = 0.045（KO）/ 0.011（KD）——单细胞层近似不相关，与组织层"解耦"方向一致（零膨胀压低，描述性）
- KO vs KD 全屏 Spearman r = 0.11–0.18（p ≤ 0.009）；臂成员方向一致 4–6/7

## 4. 判定与措辞

- **腿3 = GO 且分析执行；预指定检验全阴 → 按注册条款如实入稿。** 本数据不构成对双臂结构的额外扰动支持，也不构成反驳（臂基因覆盖 7/63 的把握度局限如实披露，为向零偏倚）。
- M14 判定门不受影响：由腿1（Torin 正式检验，单侧 p=0.0005，负对照校准通过）承担，已通过。
- 定位不变：本层结果为预测层/三角化证据，不进中心结论。
- 文稿落点：Results §13.6 已按本报告改写（中文/英文双版）；Table S72–S75；局限与未来方向相应补记。

## 5. 留痕

- 中间表：`_intermediate/M14L3_celllevel_effects.csv`、`M14L3_frp_projections.csv`、`M14L3_arm_member_summary.csv`、`M14L3_cell_meta.csv`
- 附表：Table_S72–S75（见 `README_附表索引.md` 新模块段落）
- 数据入库：GSE221321_RAW.tar → `00_RAW_DATA/GSE221321_THP1_PerturbSeq/`（2026-08-26，tar 18 条目校验通过，gz CRC 校验通过）
- 复现：`M14_step2_perturbseq_gse221321.py`（seed=0；运行时约 90s；依赖 .venv_pyaging：numpy/pandas/scipy/anndata/h5py）
