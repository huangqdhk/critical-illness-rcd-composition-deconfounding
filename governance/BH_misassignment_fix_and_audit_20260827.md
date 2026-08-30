# mdi_lib.bh 错位赋值缺陷：历史影响审计与勘误再生记录（2026-08-27）

> 性质：共享统计库缺陷的全面影响审计 + 受累冻结输出的勘误再生。
> 触发：M11/M12 批次（2026-08-27）复用 `mdi_lib.bh` 时发现并当日修复
> （`M11M12_Analysis_Report_20260827.md` §6 将两个历史调用方标记为"待作者复核的历史项"）——
> 本记录即该复核的执行与闭环。

---

## 一、缺陷描述

- **原实现**（已不存在于当前 `mdi_lib.py`，修复版带 2026-08-27 docstring）：
  BH q 值按升序计算后**未反排回原输入位置**直接返回——当输入 p 数组非升序时，
  q 值被错位赋给错误条目（"最小 p 不在首位时 q 值错位"）。
- **修复版**：`q[order] = q_sorted`（`order = argsort(p)`），标准 BH + 反向累积最小，
  与 `P0_donor_pseudobulk.py` / `P2_gate2_headtohead.py` 内联的正确实现一致。

## 二、全仓调用方清点（2026-08-27 grep 复核）

| 调用方 | BH 来源 | 状态 |
|---|---|---|
| M2_step2_replication_meta.py:60 | `mdi_lib.bh`（缺陷版） | **受累**（本记录处理） |
| M14_step2_perturbseq_gse221321.py:281 | `mdi_lib.bh`（缺陷版） | **受累**（本记录处理） |
| M11M12_step1–4、M11M12_fix_identifiers | 修复版（08-27 当日） | 不受累 |
| P0_donor_pseudobulk.py:46 | 自有内联正确实现 | 不受累（单核 UCS q=0.0045 等供者级数值安全） |
| P2_gate2_headtohead.py:41 | 自有内联正确实现 | 不受累（Figure 6F/S15 族 BH 安全） |

## 三、受累输出与复算结果

### 1) `_intermediate/M2_contrasts.csv` → Table S57（M2 临床化对比表）

- BH 族 = 每 (cohort, contrast) 内 3 检验（UCS/EIS/MDI 的 MW_p）。
- 复算：**5/7 组存在 q 值错位**（GSE185263 Sepsis、GSE212865 COVID、GSE212865 SDRA、
  GSE32707 ARDSd0、GSE32707 Sepsisd0）。
- **文稿引用面（唯一）**：Results §6 判定段"GSE32707 在 BH 校正下 q=0.096"——
  该值实为 EIS 的 p 被错位赋给 MDI；**正确值 MDI BH q = 0.0225**（p=0.01497 × 3/2）。
  - 修正方向：q 从 0.096 → **0.0225**（变显著，R1(i) 判定只会更强，**R2 判定不变**）。
  - 双版文稿已同步勘误（v4.4 增量，就地括注勘误来源）。
- 其余错位值（GSE185263 UCS/MDI 互换、GSE212865 两对比等）未被文稿正文引用，
  随表格再生一并修正。
- meta 合并行（n_cohorts=3）无 BH 列值，不受影响；`check_m2_mdi` lint 锚点
  （S57 meta pooled_g≈+1.85、S56 4,888 行）不受影响。

### 2) `_intermediate/M14L3_celllevel_effects.csv` → Table S73（M14 腿3 逐基因扰动效应）

- BH 族 = 每 dataset（KO/KD）× 每 score 的全基因 mwu_p（~594–597 行/族）。
- 复算：10 个 (dataset × score) 族全部存在错位（每族 593–596/597 行 q 改变）。
- **文稿引用面：零**。§13.6 引用的均为 NTC 置换 p / 名义 p（`ntc_perm_p_*`、`mwu_p_*`），
  不含 bh_q_* 值；S75 集合检验与 S74 投影不走该函数。判定（预指定检验全阴）不变。

## 四、勘误再生（`03_LOGS/_m1x_closeout_20260827.py`，字符串级重写）

- 仅重写 BH 列（S57 `BH_q`、S73 `bh_q_UCS/EIS/MDI/UCS_nomt/MDI_nomt`），
  其余列文本原样保留；`_intermediate` 两份同步再生，保持与发表表一致。
- RESULTS_MANIFEST_v2.0 三行（S57/S73/S69c）sha256/size 已更新并加注勘误标签。
- S69c 同批追加 CRC 反向审计 not_executable 闭环行（独立事项，见
  `03_LOGS/M13_reverse_audit_CRC_closure_20260827.md`）。

## 五、结论

1. 缺陷影响限于两个冻结输出的 BH 列；**文稿唯一受累数字为 §6 的 q=0.096→0.0225**，
   已双版勘误；所有预注册判定（R2、M14 腿3 全阴）复核后均不变。
2. 修复版 `mdi_lib.bh` 自 2026-08-27 起为唯一在用实现；历史调用方已全部清点，
   无第三处受累。
3. 教训入库：共享库函数修复后必须 grep 全仓调用方并逐个复核冻结输出——本次
   M11/M12 报告的"待作者复核"标记即按此纪律执行完毕。
