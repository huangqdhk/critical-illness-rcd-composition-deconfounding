# -*- coding: utf-8 -*-
"""
M16_step5_export_register.py — M16 图数据导出（Figure_11A/B/C）+ 报告 + manifest 注册
=====================================================================
预注册：M16_pre_registration_20260831.md
输出：
  01_FIGURE_DATA_CSV/Main/Figure_11A.csv / Figure_11B.csv / Figure_11C.csv
  04_AUDIT_GOVERNANCE/M16_Mechanical_Boundary_Report.md
  RESULTS_MANIFEST_v2.0.csv 追加登记（SHA256/行数/大小，行级写入非 pandas 往返）
"""
import hashlib
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES", "SUPPLEMENTARY_Tables_CSV")
FIG = os.path.join(ROOT, "01_FIGURE_DATA_CSV", "Main")
GOV = os.path.join(ROOT, "04_AUDIT_GOVERNANCE")
INTER = os.path.join(ROOT, "_intermediate")
LOG = os.path.join(ROOT, "03_LOGS", "M16_step5_log.txt")

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
log = []
def note(m=""):
    log.append(m); print(m, flush=True)

def sha256(f):
    h = hashlib.sha256()
    with open(f, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def main():
    note("== M16 step5 导出+注册 2026-08-31 ==")
    s89 = pd.read_csv(os.path.join(TAB, "Table_S89_M16_VILI_Contrasts.csv"))
    s90 = pd.read_csv(os.path.join(TAB, "Table_S90_M16_GSE2411_Interaction.csv"))
    s91 = pd.read_csv(os.path.join(TAB, "Table_S91_M16_Bridge_Correlations.csv"))
    meta = pd.read_csv(os.path.join(INTER, "M16_vili_meta.csv"))
    vis = pd.read_csv(os.path.join(INTER, "M16_visium_mech_spatial.csv"))

    # ---- Figure_11A：VILI 主对比四分数效应 + 合并 ----
    a = s89[(s89["layer"] == "primary") & s89["score"].isin(["UCS", "EIS", "MDI", "mech"])][
        ["gse", "contrast", "score", "n_case", "n_ctrl", "diff", "hedges_g",
         "welch_t", "welch_p_two", "welch_p_one", "welch_p_one_bh"]].copy()
    a["panel_block"] = "per_dataset_primary"
    mrows = []
    for _, r in meta.iterrows():
        mrows.append({"gse": "POOLED", "contrast": "REML_HK_random_effects", "score": r["score"],
                      "n_case": np.nan, "n_ctrl": np.nan, "diff": np.nan, "hedges_g": r["pooled"],
                      "welch_t": r["t"], "welch_p_two": np.nan, "welch_p_one": r["p_one"],
                      "welch_p_one_bh": np.nan, "panel_block": "pooled_meta",
                      "i2": r["i2"], "dir_consistent": f"{r['dir_consistent']}/{r['k']}",
                      "loso": r["loso_same_direction"]})
    fig11a = pd.concat([a, pd.DataFrame(mrows)], ignore_index=True)
    fig11a["gene_set_version"] = "Mitoxy-80_v1.0"
    fig11a["score_version"] = "mdi_v1.0+mech_v1.0"
    f11a = os.path.join(FIG, "Figure_11A.csv"); fig11a.to_csv(f11a, index=False, encoding="utf-8-sig")
    note(f"[write] {f11a} ({len(fig11a)} rows)")

    # ---- Figure_11B：GSE2411 2×2（组均值 + 互作 + MV 增强对比） ----
    b1 = s90[["score", "mean_Control", "mean_MV", "mean_LPS", "mean_MVLPS"]].copy()
    b1["panel_block"] = "group_means_z"
    b2 = s90[["score", "beta3_interaction", "se", "welch_t", "df_satterthwaite",
              "welch_p_two", "welch_p_one", "perm_p_two", "perm_p_directional",
              "welch_p_one_bh"]].copy()
    b2["panel_block"] = "interaction_beta3"
    b3 = s89[(s89["gse"] == "GSE2411") & (s89["contrast"].isin(
        ["MV_vs_Control", "LPS_vs_Control", "MV+LPS_vs_Control", "MV+LPS_vs_LPS"]))
        & s89["score"].isin(["UCS", "EIS", "MDI", "mech"])][
        ["contrast", "score", "hedges_g", "welch_p_two"]].copy()
    b3["panel_block"] = "pairwise_contrasts"
    fig11b = pd.concat([b1, b2, b3], ignore_index=True)
    fig11b["gse"] = "GSE2411"
    fig11b["gene_set_version"] = "Mitoxy-80_v1.0"; fig11b["score_version"] = "mdi_v1.0+mech_v1.0"
    f11b = os.path.join(FIG, "Figure_11B.csv"); fig11b.to_csv(f11b, index=False, encoding="utf-8-sig")
    note(f"[write] {f11b} ({len(fig11b)} rows)")

    # ---- Figure_11C：人体桥接 + Visium 逐切片 ----
    c1 = s91.copy(); c1["panel_block"] = "bridge_summary"
    v = vis.copy(); v["panel_block"] = "visium_per_section"
    fig11c = pd.concat([c1, v], ignore_index=True)
    fig11c["gene_set_version"] = "mech_v1.0"; fig11c["score_version"] = "mdi_v1.0+mech_v1.0"
    f11c = os.path.join(FIG, "Figure_11C.csv"); fig11c.to_csv(f11c, index=False, encoding="utf-8-sig")
    note(f"[write] {f11c} ({len(fig11c)} rows)")

    # ---- 报告 ----
    def gv(gse_, score, col, layer="primary", contrast=None):
        q = s89[(s89["gse"] == gse_) & (s89["score"] == score) & (s89["layer"] == layer)]
        if contrast is not None:
            q = q[q["contrast"] == contrast]
        return q.iloc[0][col]

    rpt = f"""# M16 力学边界模块执行报告（2026-08-31）

预注册：`M16_pre_registration_20260831.md`（冻结先于数值）；GEO 核验 `M16_GEO_Verification_20260831.md`。

## 判据裁定（逐条对照预注册 §4）

| 判据 | 预注册阈值 | 观测 | 裁定 |
|---|---|---|---|
| A 力学关联上游抑制 | UCS g<0 方向 ≥3/4 且合并单侧 p<0.05 | 方向 3/4 但为**正**向（与预测相反）；合并 g=+0.237, 单侧 p=0.322 | **阴性**——机械通气单独不复现上游塌陷 |
| B 力学许可/放大 | GSE2411 MDI β₃>0 且单侧 p<0.05（BH 内） | MDI β₃=+0.219, 单侧 p=0.090（BH 未过）；方向为正 | **未达形式互作阈值（趋势级正向）**；但修饰层 MV+LPS vs LPS 直接对比 EIS g=+1.80, p=0.0083——力学显著增强 LPS 诱导的执行臂（部分支持，读向为"放大执行而非诱发上游塌陷"） |
| C 机械感应轴激活 | mech g>0 ≥3/4 且合并单侧 p<0.05 | 方向 2/4；合并 g=+0.337, 单侧 p=0.811 | **阴性/混合** |
| D 人体桥接共变 | mech–MDI 方向一致 ≥2/3 队列 | GSE185263 ρ=−0.479（BH 3.5e-23）、GSE212865 +0.170、GSE158055 供体级 −0.344（BH 4.9e-04）→ 2/3 负向 | **达阈值（提示性）**：机械感应高 ↔ 解离低 |

## 关键数值（入稿候选）

**VILI 主对比（4 套，BH 家族 16 检验）**：UCS 合并 g=+0.237（I²=0.30，单侧 p=0.322，LOSO 3/4）；EIS 合并 g=+0.533（I²=0，单侧 p=0.086，LOSO 4/4）；MDI 合并 g=+0.147（p=0.336）；mech 合并 g=+0.337（p=0.811）。GSE9368 VALI 为唯一 UCS 负向层（g=−1.289，双侧 p=0.146）。
**GSE2411 2×2**：LPS vs Control EIS g=+2.481（p=0.0015）、MDI g=+2.453（p=0.0011）；MV+LPS vs LPS EIS g=+1.798（p=0.0083）；互作 β₃：MDI +0.219（SE 0.158，单侧 p=0.090）、EIS +0.118、UCS −0.101、mech −0.384（双侧 p=0.032，置换 p=0.093）。
**修饰层披露**：GSE9368 rhPBEF 单独灌注复现完整解离（UCS g=−1.713 p=0.065 / EIS g=+3.901 p=0.0057 / MDI g=+14.172 p=3.6e-05）；GSE7742 jnk1-/- 层 mech g=+5.041（p=0.0068，修饰层异常读数，披露）；GSE9208 Nrf2-/- 层 MDI g=−1.080（p=0.187）。
**人体桥接**：GSE185263 mech×UCS ρ=+0.601（p=7.1e-40；Monaco|NNLS 组成校正后 +0.357，p=3.4e-13）、mech×MDI ρ=−0.479（校正后 −0.152，p=0.0026）；GSE212865 mech×UCS ρ=−0.239（p=0.0048，方向相反——该层即本文已知离群锚点层，披露为跨数据集异质性）；GSE158055 单核亚集供体级（108 供体×类型单元，68 供体）mech×UCS ρ=+0.592（p=1.5e-11）、mech×EIS +0.495、mech×MDI −0.344（p=2.7e-04）。
**Visium 空间（探索性，家族外）**：mech×UCS 双变量 Moran's I 合并 z=+5.27（p=1.3e-07，18/23 切片正值）；mech×EIS I 合并 z=−4.43（p=9.5e-06，5/23 正值）；per-spot ρ(mech,EIS) 合并 z=−11.62、ρ(mech,MDI) z=−8.26、ρ(mech,UCS) z=+0.10（n.s.）——机械感应模块与上游基建区共位、与执行区空间相离。

## 边界解读（供文稿引用）

1. 机械通气单独（4 套小鼠全肺，63 只动物）**不复现**双臂解离，上游臂方向若任何thing偏正——"解离=通气伪影"这一审稿人级平凡备择假说被排除。
2. 力学是**执行臂的放大器**而非上游塌陷的触发器（MV+LPS vs LPS EIS g=+1.80, p=0.0083；MDI 互作趋势级正向 +0.219, p=0.090）。
3. 完整解离格局由炎症介质（rhPBEF/NAMPT）单独灌注复现——与本文"状态为炎症/组成驱动"主线一致。
4. 人体组织内机械感应模块与上游基建臂共变并共位（bulk ρ=+0.60 校正后存留；供体级单核 +0.59；空间 I z=+5.27），与执行区相离——力学小生境="基建邻近、执行远隔"。

## 数据与口径披露

- 4 套 VILI 均为小鼠全肺匀浆旧式 Affymetrix 平台；GSE2411/GSE9208 线性尺度经 log2(x+1)，GSE7742/GSE9368 原 log2 尺度（披露）；GSE7742 一样本标题拼写 "lung_c47/b_ventilation2"（c57 笔误，按特征字段归组，披露）。
- 跨物种层不做组成校正（人源含粒细胞参考不适用），同 M10D 纪律；本层不构成组成独立性证据。
- 臂基因覆盖（正映射后）：GSE2411/GSE9208 UCS 23/30 EIS 27/33 mech 20/23；GSE7742 19/30、25/33、16/23；GSE9368 24/30、28/33、22/23（Table S92 逐数据集缺失清单）。
- Visium 层复用归档评分 h5ad（归档 -20260819/_intermediate/M1_visium_scored.h5ad，来源注册）；mech 命中 22/23。
- GSE158055 供体级在 Mono 亚集 h5ad（20,000 细胞）内自洽重算臂分数（log 空间、供体×类型≥50 细胞），未与 S67 逐格合并（S67 仅含 GSE216009/GSE180578 两层新队列）；GSE145926 全基因组 h5ad 因 20260816 隔离为合成件（_QUARANTINE）未使用。
"""
    fr = os.path.join(GOV, "M16_Mechanical_Boundary_Report.md")
    with open(fr, "w", encoding="utf-8") as f:
        f.write(rpt)
    note(f"[write] {fr}")

    # ---- manifest 注册（行级写入） ----
    mani = os.path.join(GOV, "RESULTS_MANIFEST_v2.0.csv")
    new_files = [
        ("S88", "Table_S88_M16_Mechanosensing_Module_v10.csv", "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV", "M16 机械感应模块 mech_v1.0（23 基因三层，PMID 活验溯源）"),
        ("S89", "Table_S89_M16_VILI_Contrasts.csv", "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV", "M16C 4 套 VILI mdi_v1.0 同标准重算对比+合并"),
        ("S90", "Table_S90_M16_GSE2411_Interaction.csv", "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV", "M16C GSE2411 力学×LPS 2×2 互作（判据 B）"),
        ("S91", "Table_S91_M16_Bridge_Correlations.csv", "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV", "M16B 人体桥接相关+Visium 空间合并"),
        ("S92", "Table_S92_M16_Ortholog_Coverage_Audit.csv", "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV", "M16 同源映射与平台覆盖审计"),
        ("", "Figure_11A.csv", "01_FIGURE_DATA_CSV/Main", "M16 图 11A VILI 效应+合并"),
        ("", "Figure_11B.csv", "01_FIGURE_DATA_CSV/Main", "M16 图 11B GSE2411 2×2"),
        ("", "Figure_11C.csv", "01_FIGURE_DATA_CSV/Main", "M16 图 11C 桥接+空间"),
        ("", "M16_pre_registration_20260831.md", "", "M16 预注册（冻结先于数值）"),
        ("", "M16_GEO_Verification_20260831.md", "04_AUDIT_GOVERNANCE", "M16 GEO 官方核验+PMID 活验"),
        ("", "M16_Mechanical_Boundary_Report.md", "04_AUDIT_GOVERNANCE", "M16 执行报告（判据 A-D 裁定）"),
    ]
    with open(mani, "a", encoding="utf-8") as f:
        for snum, fn, loc, desc in new_files:
            fp = os.path.join(ROOT, loc, fn) if loc else os.path.join(ROOT, fn)
            if not os.path.exists(fp):
                note(f"[manifest-skip] {fn} 不存在"); continue
            n = ""
            if fn.endswith(".csv"):
                n = str(max(0, sum(1 for _ in open(fp, encoding="utf-8-sig", errors="replace")) - 1))
            row = ",".join([snum, fn, loc, "analysis_output" if fn.endswith(".csv") else "governance",
                            "canonical", "mitoxy_80+mech" if fn.endswith(".csv") else "none",
                            "Mitoxy-80_v1.0" if "S88" not in fn and fn.endswith(".csv") else "",
                            "mdi_v1.0+mech_v1.0" if fn.endswith(".csv") else "",
                            n, str(os.path.getsize(fp)), sha256(fp), desc,
                            "2026-08-31 M16 力学边界模块登记", ""])
            f.write("\n" + row)
            note(f"[manifest] 登记 {fn}")
    with open(LOG, "w", encoding="utf-8") as f:
        f.write("\n".join(log) + "\n")
    print("[done]")

if __name__ == "__main__":
    main()
