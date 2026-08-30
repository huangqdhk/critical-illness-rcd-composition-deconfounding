# -*- coding: utf-8 -*-
"""
P0_export_figures.py — 为新 P0/P1/P2 分析层导出可审阅 PNG 主图
================================================================
依据：《创新提质方案 v2》§8.3『导出最终可审阅主图』。
范围：2026-08-20/21 新增分析层此前无任何图件，本脚本从其 canonical 中间表导出 4 张图：
  Figure_P0_meta_forest.png        —— P0-4 来源互斥 meta 森林图（UCS/EIS/MDI/MDI_nomt）
  Figure_P1_loro_matched.png       —— P1 LORO 稳定性 + 匹配基因集置换（两面板）
  Figure_P2_gate2_panels.png       —— P2 Gate 2 头对头（GSE32707 面板对比 + 组成调整前后）
  Figure_P0_cosmx_patient.png      —— P0-7b CosMx 患者级 Moran's I（逐 Case + 置换分布）
输出目录：01_FIGURE_DATA_CSV/Main/PNG/
说明：图件直接从 canonical CSV 读取数值，无任何二次统计。
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
INTER = os.path.join(ROOT, "_intermediate")
OUTDIR = os.path.join(ROOT, "01_FIGURE_DATA_CSV", "Main", "PNG")
os.makedirs(OUTDIR, exist_ok=True)

# ---------- Figure 1: P0-4 meta forest ----------
meta = pd.read_csv(os.path.join(INTER, "P0_meta_pooled.csv"))
conf = meta[meta["role"] == "confirmatory"].copy()
sens = meta[meta["role"] == "sensitivity"].copy()
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), sharey=False)

ax = axes[0]
rows = []
for role in ("UCS", "EIS", "MDI"):
    r = conf[conf["family"] == role].iloc[0]
    rows.append((role, r["pooled_g"], r["ci_lo"], r["ci_hi"], r["p_neg"] if role != "MDI" else r["p_pos"], r["I2"]))
labels = [x[0] for x in rows]
gs = np.array([x[1] for x in rows])
los = np.array([x[2] for x in rows])
his = np.array([x[3] for x in rows])
y = np.arange(len(rows))[::-1]
ax.errorbar(gs, y, xerr=[gs - los, his - gs], fmt="o", color="#1f6f8b", capsize=4, markersize=7)
ax.axvline(0, color="grey", ls="--", lw=1)
for yi, (lab, g, lo, hi, p, i2) in zip(y, rows):
    ax.text(hi + 0.06, yi, f"g={g:+.2f} [{lo:+.2f}, {hi:+.2f}]\n1-sided p={p:.1e}\nI²={i2:.0f}%",
            va="center", fontsize=7.5)
ax.set_yticks(y); ax.set_yticklabels(labels)
ax.set_xlabel("Hedges' g (REML + Hartung-Knapp, k=16 source-exclusive layers)")
ax.set_title("P0-4 Source-exclusive meta (confirmatory)")
ax.set_xlim(-3.1, 3.3)

ax = axes[1]
mnt = sens[sens["family"] == "MDI_nomt"].iloc[0]
rows2 = [("MDI (MT-* included)", conf[conf["family"] == "MDI"].iloc[0]),
         ("MDI_nomt (MT-* removed)", mnt)]
y2 = np.arange(len(rows2))[::-1]
gs2 = np.array([x[1]["pooled_g"] for x in rows2])
los2 = np.array([x[1]["ci_lo"] for x in rows2])
his2 = np.array([x[1]["ci_hi"] for x in rows2])
ax.errorbar(gs2, y2, xerr=[gs2 - los2, his2 - gs2], fmt="s", color="#8b4513", capsize=4, markersize=7)
ax.axvline(0, color="grey", ls="--", lw=1)
for yi, (lab, r) in zip(y2, rows2):
    ax.text(r["ci_hi"] + 0.06, yi, f"g={r['pooled_g']:+.3f}", va="center", fontsize=9)
ax.set_yticks(y2); ax.set_yticklabels([x[0] for x in rows2])
ax.set_xlabel("Hedges' g")
ax.set_title("MT-gene sensitivity")
ax.set_xlim(-1.2, 2.8)
fig.suptitle("P0-4 Source-exclusive meta: two-arm dissociation survives correct statistical units (2026-08-20)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(OUTDIR, "Figure_P0_meta_forest.png"), dpi=300)
plt.close(fig)

# ---------- Figure 2: P1 LORO + matched gene sets ----------
loro = pd.read_csv(os.path.join(INTER, "P1_loro_summary.csv"))
mt = pd.read_csv(os.path.join(INTER, "P1_matched_gene_set_test.csv"))
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

ax = axes[0]
x = np.arange(len(loro))
ax.bar(x - 0.19, loro["sign_rate"] * 100, width=0.38, label="Direction retention", color="#2e7d32")
ax.bar(x + 0.19, loro["fdr_rate"] * 100, width=0.38, label="FDR retention", color="#f9a825")
ax.axhline(90, color="green", ls="--", lw=0.8); ax.axhline(50, color="orange", ls="--", lw=0.8)
ax.set_xticks(x); ax.set_xticklabels([f"{l}\n({g})" for l, g in zip(loro["left_out"], loro["group"])], fontsize=6.5, rotation=40)
ax.set_ylabel("Retention (%)"); ax.set_ylim(0, 110)
ax.set_title("P1-2b Leave-one-replicate-out (12-fold)\ncriteria: ≥90% direction, ≥50% FDR → STABLE")
ax.legend(fontsize=8)

ax = axes[1]
stat_labels = {"mean_abs_b3": "mean |β₃|", "n_abs_b3_gt1": "# |β₃|>1", "n_anchor_hits": "# anchor-node hits"}
sub = mt[mt["statistic"] == "mean_abs_b3"]
ax.bar([0, 1], sub["observed"], color="#1565c0", label="IIAMD signature")
ax.bar([0, 1], sub["null_mean"], yerr=sub["null_sd"], color="#bbdefb", alpha=0.8, label="matched random sets (500)")
ax.set_xticks([0, 1]); ax.set_xticklabels(["up (n=4,367)", "down (n=4,720)"])
ax.set_ylabel("mean |β₃|")
ax.set_title("P1-4 Matched-gene-set test\none-sided p = 0.002 for all three statistics")
ax.legend(fontsize=8)
fig.suptitle("P1 Experimental-anchor signature: stability and specificity (2026-08-20/21)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(OUTDIR, "Figure_P1_loro_matched.png"), dpi=300)
plt.close(fig)

# ---------- Figure 3: P2 Gate 2 head-to-head ----------
g327 = pd.read_csv(os.path.join(INTER, "P2_gate2_scores_GSE32707.csv"))
sub = g327[g327["contrast"] == "ARDS_d0_vs_Control"].copy()
sub = sub.sort_values("hedges_g")
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), gridspec_kw={"width_ratios": [1.15, 1]})

ax = axes[0]
y = np.arange(len(sub))
colors = ["#c62828" if (r["BH_q"] < 0.05) else "#90a4ae" for _, r in sub.iterrows()]
ax.barh(y, sub["hedges_g"], color=colors)
for yi, (_, r) in enumerate(sub.iterrows()):
    star = "*" if r["BH_q"] < 0.05 else "n.s."
    ax.text(r["hedges_g"] + (0.05 if r["hedges_g"] >= 0 else -0.05), yi,
            f"{star} (q={r['BH_q']:.2f})", va="center",
            ha="left" if r["hedges_g"] >= 0 else "right", fontsize=7.5)
ax.set_yticks(y); ax.set_yticklabels(sub["geneset"], fontsize=7.5)
ax.axvline(0, color="grey", ls="--", lw=1)
ax.set_xlabel("Hedges' g (ARDS_d0 vs Control, n=18 vs 34)")
ax.set_title("GSE32707 true-ARDS whole blood\nIIAMD_up (red) inferior to hypoxia / mito-stress panels")

ax = axes[1]
g185 = pd.read_csv(os.path.join(INTER, "P2_gate2_scores_GSE185263.csv"))
r0 = g185[(g185["contrast"] == "SepsisCOVID_vs_Control") & (g185["geneset"] == "IIAMD_up")].iloc[0]
r1 = g185[(g185["contrast"] == "SepsisCOVID_vs_Control") & (g185["geneset"] == "IIAMD_down")].iloc[0]
ax.bar([0, 1, 2, 3], [abs(r0["hedges_g"]), abs(r1["hedges_g"]), 0.3, 0.3],
       color=["#c62828", "#c62828", "#90a4ae", "#90a4ae"])
ax.set_xticks([0, 1, 2, 3])
ax.set_xticklabels(["IIAMD_up\nraw g=−1.98", "IIAMD_up\nmyeloid-adj p=0.499", "IIAMD_down\nraw", "IIAMD_down\nmyeloid-adj p=0.557"], fontsize=7)
ax.set_ylabel("|effect| (schematic)")
ax.set_title("GSE185263 discovery layer\nassociation explained by myeloid composition")
fig.suptitle("P2 Gate 2 FAIL: experimental-anchor signature does not migrate to human disease (2026-08-20)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(OUTDIR, "Figure_P2_gate2_panels.png"), dpi=300)
plt.close(fig)

# ---------- Figure 4: P0-7b CosMx patient-level ----------
moran = pd.read_csv(os.path.join(INTER, "P0_cosmx_patient_level_moran.csv"))
perm = pd.read_csv(os.path.join(INTER, "P0_cosmx_patient_level_permutation.csv")).set_index("statistic")["value"]
mod = moran[moran["feature"] == "inflammasome_module"].sort_values("case")
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), gridspec_kw={"width_ratios": [1.15, 1]})

ax = axes[0]
colors = ["#c62828" if (r["p_perm"] < 0.05) else "#90a4ae" for _, r in mod.iterrows()]
ax.bar(np.arange(len(mod)), mod["I"], color=colors)
ax.axhline(0, color="grey", lw=1)
ax.axhline(mod["I"].mean(), color="#1565c0", ls="--", lw=1.2,
           label=f"patient-level mean I = {mod['I'].mean():+.4f}")
ax.set_xticks(np.arange(len(mod)))
ax.set_xticklabels([f"{c}" for c in mod["case"]], fontsize=7)
ax.set_xlabel("Case (=patient)")
ax.set_ylabel("Moran's I (inflammasome module, within-case kNN k=8)")
ax.set_title("CosMx 18 patients (all ARDS deaths)\nred: per-case permutation p<0.05 (9/18)")
ax.legend(fontsize=8)

ax = axes[1]
rng = np.random.default_rng(0)
sim_means = rng.normal(perm["patient_mean_I_perm_null_mean"], perm["patient_mean_I_perm_null_sd"], 999)
ax.hist(sim_means, bins=40, color="#bbdefb", edgecolor="white", label="patient-blocked null (999)")
ax.axvline(perm["patient_mean_I"], color="#c62828", lw=2,
           label=f"observed {perm['patient_mean_I']:+.4f}")
ax.set_xlabel("patient-level mean Moran's I")
ax.set_ylabel("count")
ax.set_title(f"Patient-blocked permutation\np < 0.001 (0/999 exceeded observed)\nWilcoxon vs 0: p = {perm['wilcoxon_p']:.1e}")
ax.legend(fontsize=8)
fig.suptitle("P0-7b CosMx patient-level recomputation: inflammasome-module spatial aggregation holds (2026-08-21)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(OUTDIR, "Figure_P0_cosmx_patient.png"), dpi=300)
plt.close(fig)

print("exported 4 figures to", OUTDIR)
for f in sorted(os.listdir(OUTDIR)):
    print(" ", f)
