# -*- coding: utf-8 -*-
"""
P0_spatial_patient_mapping.py — Phase 0 任务7：空间 patient–slide–FOV 映射表
================================================================================
依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md §3.5 / §4 H5（切片/FOV 不作独立生物重复）
输入：04_AUDIT_GOVERNANCE/M1_spatial_sample_manifest.csv
规则：
  - GSE271370 Visium：unit=section；sample_id 即石蜡块/患者标识（L2P/L19P/…/CONTROL2/L3C）。
    GEO 元数据未披露多切片共患者 → 默认 1 section = 1 patient，标记为【假设A，须原文核对】；
  - GSE253474 CosMx：sample_id 形如 TMA_A_8_Case1_FOV27 → slide=TMA_A_8、patient=Case1、FOV=27，
    三级嵌套确定；多 FOV/同一 Case 存在 → FOV 级检验构成伪重复（H5 关注点坐实在 CosMx 层）。
输出：_intermediate/P0_spatial_patient_slide_fov_map.csv
     04_AUDIT_GOVERNANCE/P0_Spatial_Mapping_Report.md
"""
import os, sys, re
import pandas as pd
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
MAN = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "M1_spatial_sample_manifest.csv")
OUT = os.path.join(ROOT, "_intermediate", "P0_spatial_patient_slide_fov_map.csv")
REPORT = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "P0_Spatial_Mapping_Report.md")

m = pd.read_csv(MAN)
rows = []
for _, r in m.iterrows():
    sid = str(r["sample_id"])
    if r["dataset"] == "GSE253474":
        mt = re.match(r"(TMA_[A-Z]_\d+)_(Case\d+)_(FOV\d+)", sid)
        if mt:
            slide, patient, fov = mt.groups()
        else:
            slide, patient, fov = sid, sid, sid
        rows.append(dict(dataset=r["dataset"], gsm=r["gsm"], platform="CosMx_SMI",
                         patient_id=patient, slide_id=slide, fov_id=fov,
                         condition=r["condition"], unit_source="GSM title（三级嵌套确定）",
                         assumption="none"))
    else:
        rows.append(dict(dataset=r["dataset"], gsm=r["gsm"], platform="Visium_FFPE",
                         patient_id=sid, slide_id=sid, fov_id="NA",
                         condition=r["condition"], unit_source="sample_id=石蜡块",
                         assumption="1 section = 1 patient【假设A：GEO 未披露共患者；投稿前须对原文核对】"))

df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)

vis = df[df["platform"] == "Visium_FFPE"]
cos = df[df["platform"] == "CosMx_SMI"]
lines = []
lines.append(f"Visium: {vis['slide_id'].nunique()} 切片 / {vis['patient_id'].nunique()} 患者（假设A下相等）；"
             f"分组 {vis['condition'].value_counts().to_dict()}")
lines.append(f"CosMx: {cos['fov_id'].nunique()} FOV / {cos['patient_id'].nunique()} Case / "
             f"{cos['slide_id'].nunique()} TMA slide；每 Case FOV 数分布 "
             f"{cos.groupby('patient_id').size().value_counts().to_dict()}")
multi = cos.groupby("patient_id").size()
lines.append(f"CosMx 多 FOV 共 Case：{(multi > 1).sum()}/{len(multi)} 个 Case 有 >1 FOV（最大 {multi.max()}）——"
             f"FOV 级检验构成患者内伪重复；患者级汇总可行（Case 级 n={cos['patient_id'].nunique()}）")

txt = "\n".join(lines)
print(txt)

with open(REPORT, "w", encoding="utf-8") as f:
    f.write("# P0-7 空间 patient–slide–FOV 映射报告（H5 前置）\n\n")
    f.write("- 日期：2026-08-20；冻结依据：P0_FROZEN_ANALYSIS_PLAN_v1.0.md §4 H5\n")
    f.write(f"- 映射表：_intermediate/P0_spatial_patient_slide_fov_map.csv（{len(df)} 行）\n\n")
    f.write("## 嵌套结构\n\n" + txt + "\n\n")
    f.write("## H5 检验可行性判定\n\n")
    f.write("1. **Visium（GSE271370）**：23 切片在『1 切片=1 患者』假设下即患者级单元；该假设 GEO 元数据"
            "无法直接证实，投稿前须对照原文（Delorey 2021）核对是否存在共患者多切片。若假设成立，"
            "现有 23 切片合并统计即为患者级；若否，需按患者聚合后重算。\n")
    f.write("2. **CosMx（GSE253474）**：FOV 嵌套于 Case 结构确定，原 FOV 级检验（如空间聚集置换）"
            "的有效 n 应以 Case 数而非 FOV 数计；按冻结计划，患者（Case）内汇总或患者阻断置换为正确口径，"
            "本轮先建立映射并披露，置换重算列于后续（依赖 M1_visium_scored.h5ad / CosMx 单细胞表）。\n")
    f.write("3. 按注册回退条款：在完成患者口径重算前，空间层数字以描述性身份呈现，不进入确证结论。\n")
print("\nDONE P0-7 ->", REPORT)
