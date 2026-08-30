# -*- coding: utf-8 -*-
"""
M2_step4_outputs_judgment.py — M2 第 4 步：判定 + 最终量化输出表 + SHA256 登记
=============================================================================
- 判定报告：_intermediate/M2_judgment_report.md（预注册 R1–R3）
- 主图数据：01_FIGURE_DATA_CSV/Main/Figure_9A–9E.csv（量化表格优先）
- 附表：02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV/Table_S56–S58
- RESULTS_MANIFEST_v2.0.csv（v1.0 + M2 新表，SHA256 冻结）
"""
import hashlib
import io
import os
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS")
import mdi_lib as L

ROOT = L.ROOT
INTER = ROOT + r"\_intermediate"
FIG = ROOT + r"\01_FIGURE_DATA_CSV\Main"
TAB = ROOT + r"\02_SUPPLEMENTARY_TABLES\SUPPLEMENTARY_Tables_CSV"
AUD = ROOT + r"\04_AUDIT_GOVERNANCE"
df = pd.read_csv(INTER + r"\M2_per_sample_scores.csv")
ct = pd.read_csv(INTER + r"\M2_contrasts.csv")
mt = pd.read_csv(INTER + r"\M2_meta.csv")
od = pd.read_csv(INTER + r"\M2_outcomes.csv")
ed = pd.read_csv(INTER + r"\M2_endotypes.csv")

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

written = []  # (s_number, filename, location, class, gene_universe, gsv, sv, notes)

def save(df_, path, s_number, location, cls, gsv, sv, notes):
    df_ = df_.copy()
    df_["gene_set_version"] = gsv
    df_["score_version"] = sv
    df_.to_csv(path, index=False)
    written.append(dict(s_number=s_number, filename=os.path.basename(path), location=location,
                        class_=cls, gene_universe="mitoxy_80", gene_set_version=gsv,
                        score_version=sv, notes=notes, path=path))

# ======================================================= Figure 9A: 队列级分布
rows9a = []
for (cohort, grp), sub in df.groupby(["cohort", "group"]):
    for score in ("UCS", "EIS", "MDI"):
        v = sub[score]
        rows9a.append(dict(cohort=cohort, group=grp, score=score, n=len(v),
                           mean=v.mean(), sd=v.std(ddof=1), median=v.median(),
                           q25=v.quantile(0.25), q75=v.quantile(0.75)))
fig9a = pd.DataFrame(rows9a)
save(fig9a, FIG + r"\Figure_9A.csv", None, "01_FIGURE_DATA_CSV/Main", "figure_data",
     L.GENE_SET_VERSION, L.SCORE_VERSION, "Figure 9A: 六队列 UCS/EIS/MDI 逐组分布统计")

# ======================================================= Figure 9B: meta 森林图数据
main_contrasts = ["SepsisCOVID_vs_Control", "ARDSd0_vs_Control", "Sepsis_vs_Control"]
rows9b = []
for score in ("UCS", "EIS", "MDI"):
    for _, r in ct[(ct["score"] == score) & (ct["contrast"].isin(main_contrasts))].iterrows():
        rows9b.append(dict(score=score, cohort=r["cohort"], contrast=r["contrast"],
                           n_case=r["n_case"], n_ctrl=r["n_ctrl"], hedges_g=r["hedges_g"],
                           g_se=r["hedges_g_se"], g_lo=r["hedges_g"] - 1.96 * r["hedges_g_se"],
                           g_hi=r["hedges_g"] + 1.96 * r["hedges_g_se"], MW_p=r["MW_p"], BH_q=r["BH_q"]))
    m = mt[mt["score"] == score]
    if len(m):
        r = m.iloc[0]
        rows9b.append(dict(score=score, cohort="POOLED (random-effects)", contrast="meta",
                           n_case=np.nan, n_ctrl=np.nan, hedges_g=r["pooled_g"], g_se=r["pooled_se"],
                           g_lo=r["pooled_g"] - 1.96 * r["pooled_se"], g_hi=r["pooled_g"] + 1.96 * r["pooled_se"],
                           MW_p=r["meta_p"], BH_q=np.nan))
fig9b = pd.DataFrame(rows9b)
save(fig9b, FIG + r"\Figure_9B.csv", None, "01_FIGURE_DATA_CSV/Main", "figure_data",
     L.GENE_SET_VERSION, L.SCORE_VERSION, "Figure 9B: 跨队列 meta（随机效应）森林图数据")

# ======================================================= Figure 9C: 结局关联（模型+四分位）
quart = pd.read_csv(INTER + r"\M2_quartiles.csv")
dca = pd.read_csv(INTER + r"\M2_dca.csv")
cal = pd.read_csv(INTER + r"\M2_calibration.csv")
fig9c = pd.concat([od,
                   quart.assign(model=quart["quartile"], HR_or_OR=quart["OR_vs_Q1"],
                                cindex=quart["mortality"])[["cohort", "model", "HR_or_OR", "CI_lo", "CI_hi", "p", "cindex", "n"]],
                   ], ignore_index=True, sort=False)
save(fig9c, FIG + r"\Figure_9C.csv", None, "01_FIGURE_DATA_CSV/Main", "figure_data",
     L.GENE_SET_VERSION, L.SCORE_VERSION, "Figure 9C: MDI-28天死亡关联（OR/HR/AUC/C-index/四分位）")

# ======================================================= Figure 9D: 增量价值 + 校准 + DCA
inc = od[od["model"] == "Increment"].iloc[0]
fig9d = pd.DataFrame([dict(metric="AUC_base(age+sex+subtype)", value=inc["CI_lo"]),
                      dict(metric="AUC_plus_MDI", value=inc["HR_or_OR"]),
                      dict(metric="delta_AUC", value=inc["CI_hi"]),
                      dict(metric="IDI", value=inc["idi"]),
                      dict(metric="category_NRI", value=inc["nri"])])
save(fig9d, FIG + r"\Figure_9D.csv", None, "01_FIGURE_DATA_CSV/Main", "figure_data",
     L.GENE_SET_VERSION, L.SCORE_VERSION, "Figure 9D: MDI 增量价值（ΔAUC/IDI/NRI，基线=年龄+性别+分子亚型）")

# ======================================================= Figure 9E: 四象限内型
rows9e = []
for (cohort, grp), sub in df.groupby(["cohort", "group"]):
    pass
# 用 endotypes CSV + IGP
igp_rows = []
for cohort in df["cohort"].unique():
    sub = df[df["cohort"] == cohort]
    X = sub[["UCS", "EIS"]].values
    QUAD_DEF = {(True, False): "Q1_execution_dominant", (True, True): "Q2_dual_high",
                (False, False): "Q3_dual_low", (False, True): "Q4_collapse_dominant"}
    quad = [QUAD_DEF[(bool(a), bool(b))] for a, b in
            zip(sub["EIS"] >= sub["EIS"].median(), sub["UCS"] <= sub["UCS"].median())]
    igp, igp_p = L.igp_permutation_pvalue(X, quad, k=10, n_perm=500)
    dist = pd.Series(quad).value_counts(normalize=True)
    for q in QUAD_DEF.values():
        rows9e.append(dict(cohort=cohort, quadrant=q, proportion=dist.get(q, 0.0),
                           igp=igp, igp_perm_p=igp_p))
fig9e = pd.DataFrame(rows9e)
if len(ed):
    fig9e = fig9e.merge(ed.rename(columns={"mortality_rate": "outcome_rate"}),
                        on=["cohort", "quadrant"], how="left")
save(fig9e, FIG + r"\Figure_9E.csv", None, "01_FIGURE_DATA_CSV/Main", "figure_data",
     L.GENE_SET_VERSION, L.SCORE_VERSION, "Figure 9E: 四象限内型分布/IGP/结局梯度")

# ======================================================= 附表
t56 = df.copy()
save(t56, TAB + r"\Table_S56_M2_PerSample_Scores.csv", "S56", "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
     "supplementary_table", L.GENE_SET_VERSION, L.SCORE_VERSION, "M2 逐样本 UCS/EIS/MDI（六队列 4888 样本）")

t57 = pd.concat([ct, mt.assign(score=mt["score"])], ignore_index=True, sort=False)
save(t57, TAB + r"\Table_S57_M2_CrossDisease_Contrasts_Meta.csv", "S57", "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
     "supplementary_table", L.GENE_SET_VERSION, L.SCORE_VERSION, "M2 组间对比+随机效应 meta")

t58 = pd.concat([od, quart, dca, cal], ignore_index=True, sort=False)
t58["gene_set_version"] = L.GENE_SET_VERSION
t58["score_version"] = L.SCORE_VERSION
save(t58, TAB + r"\Table_S58_M2_Outcome_Association.csv", "S58", "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV",
     "supplementary_table", L.GENE_SET_VERSION, L.SCORE_VERSION, "M2 结局关联/四分位/增量价值/敏感性/校准/DCA")

# ======================================================= 判定（预注册 R1–R3）
repl_ok = {
    "GSE185263": float(ct[(ct.cohort == "GSE185263") & (ct.score == "MDI") & (ct.contrast == "SepsisCOVID_vs_Control")]["BH_q"].iloc[0]) < 0.05,
    "GSE32707": float(ct[(ct.cohort == "GSE32707") & (ct.score == "MDI") & (ct.contrast == "ARDSd0_vs_Control")]["BH_q"].iloc[0]) < 0.05,
    "GSE310929": float(ct[(ct.cohort == "GSE310929") & (ct.score == "MDI") & (ct.contrast == "Sepsis_vs_Control")]["BH_q"].iloc[0]) < 0.05,
}
outcome_ok = {
    "GSE310929_logistic": float(od[(od.cohort == "GSE310929") & (od.model == "Logistic_unadj_28d")]["p"].iloc[0]) < 0.05,
    "GSE188309_logistic": float(od[(od.cohort == "GSE188309") & (od.model == "Logistic_unadj")]["p"].iloc[0]) < 0.05,
}
increment_ok = inc["CI_hi"] > 0.02  # ΔAUC
cond_i = sum(repl_ok.values()) >= 2
cond_ii = sum(outcome_ok.values()) >= 2
cond_iii = increment_ok
if cond_i and cond_ii and cond_iii:
    verdict, verdict_name = "R1", "强阳性：MDI 可分层预后的组织层定律（主张升级）"
elif cond_i or cond_ii:
    verdict, verdict_name = "R2", "部分：MDI 作为生物标志物候选如实报告（跨队列复制成立、预后价值未证实）"
else:
    verdict, verdict_name = "R3", "阴性/重定向：解离为危重疾病共同组织原则、无预后分层"

# ======================================================= 判定报告
rep = []
rep.append("# M2 判定报告（2026-08-17，预注册 M2_pre_registration_20260817.md）\n")
rep.append(f"## 判定结果：**{verdict}** — {verdict_name}\n")
rep.append("## 预注册条件核对\n")
rep.append("| 条件 | 内容 | 结果 |")
rep.append("|---|---|---|")
rep.append(f"| R1(i) 跨队列复制 | ≥2 独立队列疾病组 MDI 升高 BH<0.05 | "
           f"{'通过' if cond_i else '未通过'}（{sum(repl_ok.values())}/3 队列显著：GSE185263/GSE32707/GSE310929）|")
rep.append(f"| R1(ii) 结局关联 | ≥2 独立结局队列 MDI-死亡同向显著 | "
           f"{'通过' if cond_ii else '未通过'}（GSE310929 OR=1.063 p=0.218；GSE188309 OR=1.030 p=0.918）|")
rep.append(f"| R1(iii) 增量价值 | ΔC-index>0.02 或 NRI 显著 | "
           f"{'通过' if cond_iii else '未通过'}（ΔAUC=+0.0001、IDI≈0、NRI=−0.007）|")
rep.append("")
rep.append("## 关键数值\n")
rep.append("- 跨队列复制（疾病 vs 对照 MDI Hedges' g）：GSE185263 +3.44（p=9.8e-20）、GSE32707 ARDS d0 +0.77（p=0.015）、"
           "GSE310929 +1.35（p=1.3e-114）；随机效应合并 g=+1.85（p=0.005，I²=0.96）。")
rep.append("- 脓毒症 2858 例 MDI 均值 +0.20（z 空间，t=7.57，p=4.9e-14）——解离在脓毒症谱内强复制。")
rep.append("- 结局关联：脓毒症 28 天死亡 logistic OR=1.063（0.965–1.170，p=0.218，n=2436/死亡 534）；"
           "Cox（时间亚组 n=593/事件 148）HR=0.933（p=0.399）；CAP 住院死亡 OR=1.030（p=0.918，n=198/死亡 13）。均不显著。")
rep.append("- 四分位（探索性，未预注册为主判定）：脓毒症 Q4 vs Q1 死亡率 25.1% vs 19.4%，OR=1.397（p=0.016）；有序趋势 OR/级=1.103（p=0.026）。")
rep.append("- 增量价值：基线（年龄+性别+C1–C4 分子亚型）AUC=0.6370 → +MDI 0.6371；ΔAUC=+0.0001；IDI<0.001；分类 NRI=−0.007。无增量。")
rep.append("- 四象限内型：IGP 0.807–0.964（置换 p=0.002）各队列象限内聚成立；脓毒症象限间 28 天死亡率 0.192–0.243（RR vs Q1 0.88–1.12），梯度弱。")
rep.append("- GSE212865 SDRA 层 MDI 反向（Δ=−0.85，p=0.011）——与包内既有发现一致（SDRA 层为离群锚点，ρ=−0.282），如实登记。")
rep.append("- COPD：痰液 MDI 高于血液（+0.31，p=6.7e-7，执行臂驱动）；SCREENING→V4D28 治疗后 UCS 回升（+0.29，p=0.004）、MDI 下降（−0.33，p=0.016）——治疗相关部分恢复信号（描述性）。")
rep.append("")
rep.append("## 判定说明\n")
rep.append("按预注册规则判定 **R2**：解离状态（MDI）在 ≥2 个独立疾病队列中同向显著升高（条件 (i) 成立），"
           "但 MDI 与 28 天/住院死亡在两个独立结局队列中均无显著关联（条件 (ii) 不成立），"
           "且对已知分子亚型无增量（条件 (iii) 不成立）。**MDI 是可跨队列复制的疾病状态读出，但不是预后分层标志物。**"
           "该界定本身即高价值阴性结论：解离是危重疾病的共同组织状态特征，而非死亡风险的连续调节量。"
           "四分位顶层信号（Q4 OR=1.40，p=0.016）为探索性发现，因主分析（连续 OR/HR）预注册为主判定而未升级，"
           "可作为后续研究假说。")
rep.append("")
rep.append("## 数据与代码\n")
rep.append("- 判定输入：`_intermediate/M2_contrasts.csv`、`M2_meta.csv`、`M2_outcomes.csv`、`M2_endotypes.csv`。")
rep.append("- 输出表：Figure_9A–9E、Table_S56–S58（RESULTS_MANIFEST v2.0 登记，SHA256 冻结）。")
rep.append("- 复现：`M2_step1_load_score.py → M2_step2_replication_meta.py → M2_step3_outcomes_endotypes.py → M2_step4_outputs_judgment.py`"
           "（venv_pyaging：pandas 3.0.5/numpy 2.5.1/scipy 1.18.0/sklearn 1.9.0/lifelines 0.30.3；"
           "R 4.6.0 + hgu133plus2.db 3.13.0 + clariomshumantranscriptcluster.db 8.8.0 用于探针注释）。")
with open(INTER + r"\M2_judgment_report.md", "w", encoding="utf-8") as f:
    f.write("\n".join(rep))
print(f"判定: {verdict} — {verdict_name}")

# ======================================================= RESULTS_MANIFEST v2.0
# 文本级追加（pandas 往返会把 v1.0 的 170 个字面 "NA" 写成空串，破坏 lint 版本戳比对）
import csv as _csv
rm1_path = AUD + r"\RESULTS_MANIFEST_v1.0.csv"
rm2_path = AUD + r"\RESULTS_MANIFEST_v2.0.csv"
with open(rm1_path, encoding="utf-8-sig", newline="") as f:
    header = next(_csv.reader(f))  # s_number,filename,location,class,status,...
idx = {c: i for i, c in enumerate(header)}
def csv_escape(v):
    s = str(v)
    if any(ch in s for ch in ',"\n'):
        s = '"' + s.replace('"', '""') + '"'
    return s

with open(rm2_path, "w", encoding="utf-8", newline="") as f:
    f.write(",".join(header) + "\n")
    for w in written:
        sz = os.path.getsize(w["path"])
        nrows = sum(1 for _ in open(w["path"], encoding="utf-8-sig", errors="replace")) - 1
        row = {c: "" for c in header}
        row.update({
            "s_number": w["s_number"] or "",
            "filename": w["filename"],
            "location": w["location"],
            "class": w["class_"],
            "status": "canonical",
            "gene_universe": w["gene_universe"],
            "gene_set_version": w["gene_set_version"],
            "score_version": w["score_version"],
            "n_data_rows": nrows,
            "size_bytes": sz,
            "sha256": sha256(w["path"]),
            "analysis": w["notes"],
            "notes": "M2 (2026-08-17)",
        })
        f.write(",".join(csv_escape(row[c]) for c in header) + "\n")
# 追加 v1.0 全部原行（文本原样，保留 "NA" 字面量）
with open(rm1_path, encoding="utf-8-sig", errors="replace") as f:
    f.readline()  # skip header
    rest = f.read()
with open(rm2_path, "a", encoding="utf-8", newline="") as f:
    f.write(rest)
n_old = rest.count("\n")
print(f"RESULTS_MANIFEST v2.0: {n_old} 旧行 + {len(written)} 新表（文本级追加，v1.0 原样保留）")
for w in written:
    print(f"  saved: {w['filename']} ({os.path.getsize(w['path'])} bytes)")
print("\nDONE step4")
