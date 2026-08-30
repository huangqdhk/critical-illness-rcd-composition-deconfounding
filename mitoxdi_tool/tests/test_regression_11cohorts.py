# -*- coding: utf-8 -*-
"""
test_regression_11cohorts.py — M15 工具内部回归测试（11 队列）
================================================================================
测试对象：M15_tool/mitoxdi 工具核心（两臂分解 + mdi_v1.0 评分 + 组成预测值/组成残差 MDI）。
参照：本项目冻结的逐样本评分表（前视注册批次产物，作工具内部回归测试集）：
  - M2_per_sample_scores.csv            （6 队列：UCS/EIS/MDI/nomt/ssgsea 逐样本）
  - M10A_synthetic_scores.csv + M10A_residual_per_sample.csv（6 队列：MDI_comp/MDI_resid）
  - M11M12_step1_per_sample.csv         （5 队列：评分 + 组成预测 + 组成残差 + 细胞比例）
共 11 个队列 × 3 类核对。全部加载器精确复刻原流水线口径（见 mitoxdi/loaders.py）。
判定：每族 max|Δ| ≤ tol（1e-4；GSE215865 为 float32 口径，tol=5e-3）。
输出：_intermediate/M15_regression_summary.csv / M15_regression_cohort_summary.csv；
      03_LOGS/M15_regression_log.txt。
"""
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.dirname(HERE)
sys.path.insert(0, TOOL)

import numpy as np
import pandas as pd

import mitoxdi
from mitoxdi import loaders as LD
from mitoxdi.mdi_core import load_manifest, arm_scores, sample_characteristics
from mitoxdi.composition import (load_monaco_ref, load_abis_refs, deconvolve_nnls,
                                 synthetic_scores, score_and_compose)

ROOT = mitoxdi.ROOT
INTER = os.path.join(ROOT, "_intermediate")
LOG = os.path.join(ROOT, "03_LOGS")
CACHE = os.path.join(INTER, "M15_regression_cache")

log = []
def note(m=""):
    log.append(m)
    print(m, flush=True)

_, arms, _ = load_manifest()
UP, EX = arms["upstream_collapse"], arms["execution_induction"]

M2CSV = os.path.join(INTER, "M2_per_sample_scores.csv")
SYNCSV = os.path.join(INTER, "M10A_synthetic_scores.csv")
RESCSV = os.path.join(INTER, "M10A_residual_per_sample.csv")
M11CSV = os.path.join(INTER, "M11M12_step1_per_sample.csv")

summary_rows = []
cohort_rows = []


def compare(new_df, ref_df, families, id_new, id_ref, tag, tol=1e-4):
    """对齐样本后逐族比较（位置对齐，iloc 索引）。"""
    ref_ids = [str(v) for v in ref_df[id_ref]]
    new_ids = [str(v) for v in new_df[id_new]]
    ref_map = {s: i for i, s in enumerate(ref_ids)}
    new_map = {s: i for i, s in enumerate(new_ids)}
    cmn = [s for s in new_ids if s in ref_map]
    if not cmn:
        note(f"   [{tag}] 无共同样本！")
        return
    new_idx = [new_map[s] for s in cmn]
    ref_idx = [ref_map[s] for s in cmn]
    for fam in families:
        if fam not in new_df.columns or fam not in ref_df.columns:
            continue
        a = new_df[fam].iloc[new_idx].astype(float).values
        b = ref_df[fam].iloc[ref_idx].astype(float).values
        if np.isnan(a).all() and np.isnan(b).all():
            # 双全 NaN（如 Mono_monaco：冻结口径名称匹配不到 Monaco 单核列）→ 视同一致
            summary_rows.append(dict(cohort=tag.split("|")[0], check=tag.split("|")[1],
                                     family=fam, n=len(cmn), r=1.0,
                                     max_abs_diff=0.0, mean_abs_diff=0.0,
                                     tol=tol, PASS=True))
            note(f"   {tag:<34} {fam:<16} n={len(cmn):<5} 双全 NaN（口径一致） -> PASS")
            continue
        d = np.abs(a - b)
        if a.std(ddof=1) > 0 and b.std(ddof=1) > 0:
            r = float(np.corrcoef(a, b)[0, 1])
        else:  # 常量列（如 n_up/n_ex/b0/b1）：直接以最大差判定
            r = 1.0 if np.nanmax(d) <= tol else 0.0
        ok = (np.nanmax(d) <= tol) and np.isfinite(r)
        summary_rows.append(dict(cohort=tag.split("|")[0], check=tag.split("|")[1],
                                 family=fam, n=len(cmn), r=r,
                                 max_abs_diff=float(np.nanmax(d)),
                                 mean_abs_diff=float(np.nanmean(d)),
                                 tol=tol, PASS=bool(ok)))
        note(f"   {tag:<34} {fam:<16} n={len(cmn):<5} r={r:+.8f} max|Δ|={np.nanmax(d):.3e} "
             f"mean|Δ|={np.nanmean(d):.3e} -> {'PASS' if ok else 'FAIL'}")
    return cmn


def close_cohort(cohort, check, n_fam, n_pass):
    cohort_rows.append(dict(cohort=cohort, check=check, n_families=n_fam, n_pass=n_pass,
                            PASS=bool(n_pass > 0 and n_pass == n_fam)))


t0 = time.time()

# ======================================================================
# A. M2 评分回归（6 队列）：UCS/EIS/MDI(+nomt/ssgsea)
# ======================================================================
note("== A. M2 评分回归（mdi_v1.0 冻结表）")
m2 = pd.read_csv(M2CSV, low_memory=False)
m2["sample"] = m2["sample"].astype(str)
FAM_SCORE = ["UCS_raw", "EIS_raw", "UCS", "EIS", "MDI",
             "UCS_nomt", "MDI_nomt", "UCS_ssgsea", "EIS_ssgsea"]

m2_checks = [
    ("GSE185263", LD.load_gse185263_m2()),
    ("GSE32707", LD.load_gse32707()),
    ("GSE212865", LD.load_gse212865_m2(arms)),
    ("GSE188309", LD.load_gse188309_m2(arms)),
    ("GSE310929", LD.load_gse310929()),
    ("GSE148871", LD.load_gse148871_m2(arms)),
]
for cname, mat in m2_checks:
    note(f"\n-- {cname}")
    sc = arm_scores(mat, arms)
    sc["sample"] = [str(i) for i in mat.index]
    ref = m2[m2["cohort"] == cname]
    tag = f"{cname}|M2_score"
    before = len(summary_rows)
    compare(sc, ref, FAM_SCORE, "sample", "sample", tag, tol=1e-4)
    n_new = len(summary_rows) - before
    n_pass = sum(r["PASS"] for r in summary_rows[-n_new:])
    close_cohort(cname, "M2_score", n_new, n_pass)
    note(f"   n={len(ref)} 行参照；核对族 {n_new}，通过 {n_pass}")

# 工具包 vs 项目根 mdi_lib 的直接等价性（同一矩阵双实现）
note("\n== A2. mitoxdi.mdi_core vs 项目根 mdi_lib 等价性（GSE32707 矩阵）")
sys.path.insert(0, ROOT)
import mdi_lib as Lroot  # noqa: E402
mat7 = LD.load_gse32707()
sc_tool = arm_scores(mat7, arms)
sc_root = Lroot.arm_scores(mat7, arms)
d = np.abs(sc_tool[["UCS", "EIS", "MDI"]].values - sc_root[["UCS", "EIS", "MDI"]].values).max()
note(f"   max|Δ|(UCS/EIS/MDI) = {d:.3e} -> {'PASS' if d < 1e-12 else 'FAIL'}")

# ======================================================================
# B. M10A 组成预测/残差回归（6 队列）：MDI_comp / MDI_resid
# ======================================================================
note("\n== B. M10A 组成预测与残差回归（Monaco|NNLS 主口径）")
monaco = load_monaco_ref()
syn = pd.read_csv(SYNCSV)
syn["sample"] = syn["sample"].astype(str)
res = pd.read_csv(RESCSV)
res["sample"] = res["sample"].astype(str)

m10a_checks = [
    ("GSE185263", LD.load_gse185263_m10a(), "GSE185263"),
    ("GSE32707", LD.load_gse32707(), "GSE32707"),
    ("GSE212865", LD.load_gse212865_m10a(monaco.columns), "GSE212865"),
    ("GSE188309", LD.load_gse188309_m10a(monaco.columns), "GSE188309"),
    ("GSE310929", LD.load_gse310929(), "GSE310929"),
    ("GSE148871_blood", LD.load_gse148871_m10a(monaco.columns), "GSE148871"),
]
for cname, mat, m2name in m10a_checks:
    note(f"\n-- {cname}")
    prop, common = deconvolve_nnls(mat, monaco)
    sy = synthetic_scores(prop, monaco, UP, EX)
    sy["sample"] = [str(i) for i in mat.index]
    ref = syn[(syn["cohort"] == cname) & (syn["reference"] == "Monaco")
              & (syn["method"] == "NNLS")]
    tag = f"{cname}|M10A_comp"
    before = len(summary_rows)
    cmn = compare(sy, ref, ["UCS_comp_raw", "EIS_comp_raw", "UCS_comp", "EIS_comp",
                            "MDI_comp"], "sample", "sample", tag, tol=1e-4)
    n_new = len(summary_rows) - before
    n_pass = sum(r["PASS"] for r in summary_rows[-n_new:])
    # 残差核对（观测值取自冻结 M2 表，与 M10A step6 完全同口径）。
    # 注意：残差表无 cohort 列，且 GSE310929 与 GSE185263/GSE32707 共用 GSM 标识——
    # 以 MDI_obs 与冻结 M2 表逐样本相等（1e-9 容差）消歧选出本队列行。
    o = m2[m2["cohort"] == m2name][["sample", "MDI"]]
    o = o.set_index("sample")
    syc = sy.set_index("sample")
    cmn2 = sorted(set(o.index) & set(syc.index) & set(res["sample"]))
    if len(cmn2) >= 10:
        o_map = dict(o["MDI"])
        sel = res[res["sample"].isin(cmn2)]
        sel = sel[[abs(float(r["MDI_obs"]) - o_map.get(r["sample"], np.nan)) < 1e-9
                   for _, r in sel.iterrows()]]
        A = np.column_stack([np.ones(len(cmn2)), syc.loc[cmn2, "MDI_comp"].values])
        beta, *_ = np.linalg.lstsq(A, o.loc[cmn2, "MDI"].values, rcond=None)
        rdf = pd.DataFrame({"sample": cmn2,
                            "MDI_obs": o.loc[cmn2, "MDI"].values,
                            "MDI_comp": syc.loc[cmn2, "MDI_comp"].values,
                            "MDI_resid": o.loc[cmn2, "MDI"].values - A @ beta,
                            "b0": beta[0], "b1": beta[1]})
        tag2 = f"{cname}|M10A_resid"
        before = len(summary_rows)
        compare(rdf, sel, ["MDI_obs", "MDI_comp", "MDI_resid", "b0", "b1"],
                "sample", "sample", tag2, tol=1e-4)
        n2 = len(summary_rows) - before
        n_pass += sum(r["PASS"] for r in summary_rows[-n2:])
        n_new += n2
        note(f"   残差核对 n={len(cmn2)}（b1={beta[1]:+.3f}）")
    close_cohort(cname, "M10A_comp+resid", n_new, n_pass)

# ======================================================================
# C. M11M12 回归（5 队列）：评分 + 组成预测 + 残差 + 细胞比例
# ======================================================================
note("\n== C. M11M12 回归（mdi_v1.0 + Monaco/ABIS 组成口径）")
abis_rna, abis_micro = load_abis_refs()
t11 = pd.read_csv(M11CSV, low_memory=False)
t11["sample_id"] = t11["sample_id"].astype(str)

FAM_11 = ["UCS_raw", "EIS_raw", "UCS", "EIS", "MDI", "UCS_nomt", "MDI_nomt",
          "UCS_ssgsea", "EIS_ssgsea", "n_up", "n_ex",
          "Neut", "Mono", "comp", "Neut_monaco", "Mono_monaco", "comp_monaco",
          "UCS_comp_raw", "EIS_comp_raw", "UCS_comp", "EIS_comp", "MDI_comp",
          "n_up_ref", "n_ex_ref",
          "UCS_comp_raw_abis", "EIS_comp_raw_abis", "UCS_comp_abis", "EIS_comp_abis",
          "MDI_comp_abis", "MDI_resid", "b0", "b1"]

# --- GSE215865（float32 口径） ---
note("\n-- GSE215865")
mat215 = LD.load_gse215865(cache_dir=CACHE)
out215, _, _ = score_and_compose(mat215, arms, monaco, ref_main=abis_rna)
if isinstance(mat215.index, pd.MultiIndex):
    out215["sample_id"] = [f"{i[0]}|{i[1]}" for i in mat215.index]
else:  # 缓存重建路径：索引即 "subject|day_label" 字符串
    out215["sample_id"] = [str(i) for i in mat215.index]
before = len(summary_rows)
compare(out215, t11[t11["cohort"] == "GSE215865"], FAM_11, "sample_id", "sample_id",
        "GSE215865|M11M12", tol=5e-3)
n_new = len(summary_rows) - before
n_pass = sum(r["PASS"] for r in summary_rows[-n_new:])
close_cohort("GSE215865", "M11M12", n_new, n_pass)

# --- GSE54514 ---
note("\n-- GSE54514")
mat545 = LD.load_gse54514()
out545, _, _ = score_and_compose(mat545, arms, monaco, ref_main=abis_micro)
out545["sample_id"] = [str(i) for i in mat545.index]
before = len(summary_rows)
compare(out545, t11[t11["cohort"] == "GSE54514"], FAM_11, "sample_id", "sample_id",
        "GSE54514|M11M12", tol=1e-4)
n_new = len(summary_rows) - before
n_pass = sum(r["PASS"] for r in summary_rows[-n_new:])
close_cohort("GSE54514", "M11M12", n_new, n_pass)

# --- GSE148871（304 全样本评分，血样过滤后比对） ---
note("\n-- GSE148871")
mat148, tis148 = LD.load_gse148871_m11()
out148, _, _ = score_and_compose(mat148, arms, monaco, ref_main=abis_micro)
out148["sample_id"] = [str(i) for i in mat148.index]
out148 = out148[[("blood" in str(t).lower()) for t in tis148]]
before = len(summary_rows)
compare(out148, t11[t11["cohort"] == "GSE148871"], FAM_11, "sample_id", "sample_id",
        "GSE148871|M11M12", tol=1e-4)
n_new = len(summary_rows) - before
n_pass = sum(r["PASS"] for r in summary_rows[-n_new:])
close_cohort("GSE148871", "M11M12", n_new, n_pass)

# --- GSE106878 ---
note("\n-- GSE106878")
mat878 = LD.load_gse106878()
out878, _, _ = score_and_compose(mat878, arms, monaco, ref_main=abis_micro)
out878["sample_id"] = [str(i) for i in mat878.index]
before = len(summary_rows)
compare(out878, t11[t11["cohort"] == "GSE106878"], FAM_11, "sample_id", "sample_id",
        "GSE106878|M11M12", tol=1e-4)
n_new = len(summary_rows) - before
n_pass = sum(r["PASS"] for r in summary_rows[-n_new:])
close_cohort("GSE106878", "M11M12", n_new, n_pass)

# --- GSE212865 ---
note("\n-- GSE212865")
mat212 = LD.load_gse212865_m11()
out212, _, _ = score_and_compose(mat212, arms, monaco, ref_main=abis_micro)
out212["sample_id"] = [str(i) for i in mat212.index]
before = len(summary_rows)
compare(out212, t11[t11["cohort"] == "GSE212865"], FAM_11, "sample_id", "sample_id",
        "GSE212865|M11M12", tol=1e-4)
n_new = len(summary_rows) - before
n_pass = sum(r["PASS"] for r in summary_rows[-n_new:])
close_cohort("GSE212865", "M11M12", n_new, n_pass)

# ======================================================================
# 汇总
# ======================================================================
summary = pd.DataFrame(summary_rows)
cohort_sum = pd.DataFrame(cohort_rows)
summary.to_csv(os.path.join(INTER, "M15_regression_summary.csv"), index=False)
cohort_sum.to_csv(os.path.join(INTER, "M15_regression_cohort_summary.csv"), index=False)

n_fam = len(summary)
n_pass = int(summary["PASS"].sum())
n_cohort_pass = int(cohort_sum["PASS"].sum())
note(f"\n== 汇总 == 11 队列 × {n_fam} 族核对；通过 {n_pass}/{n_fam}；队列级全通过 "
     f"{n_cohort_pass}/{len(cohort_sum)}")
verdict = "ALL PASS" if n_pass == n_fam and n_cohort_pass == len(cohort_sum) else "HAS FAILURES"
note(f"判定：{verdict}（用时 {time.time()-t0:.0f}s）")
with open(os.path.join(LOG, "M15_regression_log.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(log))
sys.exit(0 if verdict == "ALL PASS" else 1)
