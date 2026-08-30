# -*- coding: utf-8 -*-
"""
P0_old_vs_new.py — Phase 0 任务6：原结果 vs 正确统计单位结果 对照表
====================================================================
对照对象：_intermediate/M2_meta.csv + M2_contrasts.csv（旧，单位缺陷）
         vs P0_dataset_effects.csv / P0_meta_pooled.csv（新，来源互斥）
         + P0_pseudobulk_arm_tests.csv（新，供者级）vs 旧 N3 细胞级（Table S2b 仅定性标注）
输出：04_AUDIT_GOVERNANCE/P0_Old_vs_New_Comparison.md（+ _intermediate/P0_old_vs_new.csv）
"""
import os, sys
import pandas as pd
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
INTER = os.path.join(ROOT, "_intermediate")
OUTCSV = os.path.join(INTER, "P0_old_vs_new.csv")
OUTMD = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "P0_Old_vs_New_Comparison.md")

def df_md(d, cols=None):
    if d is None or len(d) == 0:
        return "(empty)"
    cols = cols or list(d.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, rr in d.iterrows():
        vals = [f"{rr[c]:.4g}" if isinstance(rr[c], (int, float)) and pd.notna(rr[c]) else str(rr[c]) for c in cols]
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)

rows = []

# ---- 1. meta 合成层 ----
old_meta = pd.read_csv(os.path.join(INTER, "M2_meta.csv")) if os.path.exists(os.path.join(INTER, "M2_meta.csv")) else None
new_pool = pd.read_csv(os.path.join(INTER, "P0_meta_pooled.csv"))
if old_meta is not None:
    for _, o in old_meta.iterrows():
        n = new_pool[(new_pool["family"] == o["score"]) & (new_pool["role"] == "confirmatory")]
        if len(n):
            n = n.iloc[0]
            rows.append(dict(layer="meta-pooled", item=o["score"],
                             old=f"k={o['n_cohorts']}, g={o['pooled_g']:+.2f}, p={o['meta_p']:.2g}, I2={o['I2']*100:.0f}% (无效:样本双计)",
                             new=f"k={n['k']}, g={n['pooled_g']:+.2f} [{n['ci_lo']:+.2f},{n['ci_hi']:+.2f}], 单侧p={min(n['p_pos'],n['p_neg']):.2g}, I2={n['I2']:.0f}%",
                             verdict="旧值作废；新值LOSO稳定" if n.get("loso_stable") else "旧值作废；新值LOSO不稳"))

# ---- 2. 直并行数据集级效应（旧 M2_contrasts vs 新 P0_dataset_effects）----
old_ct = pd.read_csv(os.path.join(INTER, "M2_contrasts.csv")) if os.path.exists(os.path.join(INTER, "M2_contrasts.csv")) else None
new_eff = pd.read_csv(os.path.join(INTER, "P0_dataset_effects.csv"))
PAIRS = [("GSE185263", "SepsisCOVID_vs_Control"), ("GSE32707", "ARDSd0_vs_Control"),
         ("GSE212865", "SDRA_vs_Control"), ("GSE310929", "Sepsis_vs_Control")]
if old_ct is not None:
    for coh, cid in PAIRS:
        for score in ("UCS", "EIS", "MDI"):
            o = old_ct[(old_ct["cohort"] == coh) & (old_ct["contrast"] == cid) & (old_ct["score"] == score)]
            nrow = new_eff[(new_eff["stratum"] == (coh if coh != "GSE212865" else "GSE212865_base")) &
                           (new_eff["family"] == score)]
            if len(o) and len(nrow):
                o = o.iloc[0]; nrow = nrow.iloc[0]
                notes = "GSE212865 新值=患者级基线化" if coh == "GSE212865" else ""
                if coh == "GSE310929":
                    notes = "旧队列解散→28来源拆分（单来源值不再对应）"
                    verdict = "结构作废"
                elif coh == "GSE212865":
                    notes = "新值=患者级基线化（纵向重复剔除）"
                    verdict = "重算（纵向去重后效应缩小）"
                else:
                    notes = ""
                    verdict = "保留（样本级检验，单位本就正确）"
                rows.append(dict(layer="dataset-effect", item=f"{coh} {cid} {score}",
                                 old=f"g={o['hedges_g']:+.2f} (n={o['n_case']}v{o['n_ctrl']})",
                                 new=(f"g={nrow['hedges_g']:+.2f} (n={nrow['n_case']}v{nrow['n_ctrl']}) {notes}"),
                                 verdict=verdict))

# ---- 3. 单细胞层：旧细胞级（不具显著性信息） vs 新供者级 ----
new_arm = pd.read_csv(os.path.join(INTER, "P0_pseudobulk_arm_tests.csv"))
sig = new_arm[new_arm["BH_q"] < 0.05]
top4 = "; ".join(
    "{} {}/{} diff={:+.2f} q={:.3g}".format(r["dataset"], r["cell_type"], r["score"], r["diff"], r["BH_q"])
    for _, r in sig.head(4).iterrows()
)
rows.append(dict(layer="scrna-arm", item="旧：细胞级检验（把细胞当独立重复，显著性无效）",
                 old="Table S2b 细胞级 Mann-Whitney（80基因）",
                 new=f"供者级 pseudobulk：{len(new_arm)} 检验中 {len(sig)} 个 BH q<0.05；最强：{top4}",
                 verdict="单位修正后：髓系UCS塌陷与NK-MDI存活；执行臂供者级未存活"))
rows.append(dict(layer="scrna-gene", item="旧：无单细胞基因级 DE",
                 old="（S2 实为 bulk DESeq2，N3 已改标签）",
                 new="P0_pseudobulk_DE_all.csv：7 个数据集×细胞类型组合 pydeseq2（供者级）",
                 verdict="新增，供者级"))

# ---- 4. 待补：空间层 ----
rows.append(dict(layer="spatial", item="切片/FOV 级检验", old="M1 slide/FOV 级（单位缺陷）",
                 new="待 P0-7（patient–slide–FOV 映射 + 患者阻断）", verdict="待办"))

out = pd.DataFrame(rows)
out.to_csv(OUTCSV, index=False)

with open(OUTMD, "w", encoding="utf-8") as f:
    f.write("# P0-6 原结果 vs 正确统计单位结果 对照表\n\n")
    f.write("- 日期：2026-08-20；冻结依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md（OSF https://osf.io/C7RYD/）\n")
    f.write("- 判定规则：旧 k=3 meta 与 GSE310929 独立队列地位作废（P0-3 审计）；其余为重算或保留\n\n")
    f.write(df_md(out, ["layer", "item", "old", "new", "verdict"]))
    f.write("\n\n## 说明\n\n1. 旧 M2_meta（pooled g=1.85 等）因样本双计数无效，仅列作对照；\n"
            "2. 直并行（GSE185263/GSE32707）数据集级效应数字基本保留（同为样本级检验，单位本就正确）；\n"
            "3. GSE212865 旧值含纵向伪重复，新值为患者级基线化结果；\n"
            "4. 单细胞旧检验单位错误，其显著性一律不引用；新供者级结果为唯一确证口径。\n")

print("DONE P0-6 ->", OUTMD)
