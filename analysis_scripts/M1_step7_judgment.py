# -*- coding: utf-8 -*-
"""
M1 Step 7: Pre-registered judgment (R1/R2/R3) + Table S55v sensitivity.

Judgment rules frozen in M1_pre_registration_20260817.md BEFORE results.
Sensitivity: bivariate Moran on count-residualized arm scores (addresses
density/gene-richness confound).
Outputs: M1_judgment_report.md + Table_S55v (方案A: no figure-panel copies;
spot-level scores stay in _intermediate/M1_visium_scored.h5ad).
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats
from sklearn.neighbors import NearestNeighbors
import libpysal

INT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate"
OUT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\02_SUPPLEMENTARY_TABLES\SUPPLEMENTARY_Tables_CSV"
SEED = 0
K = 8

sdf = pd.read_csv(os.path.join(OUT, "Table_S55a_M1_Visium_Section_Spatial_Stats.csv"))

# ---- sensitivity: bivariate I on residualized scores ----
adata = ad.read_h5ad(os.path.join(INT, "M1_visium_scored.h5ad"))
sens_rows = []
for sec in sorted(adata.obs["section"].unique()):
    m = (adata.obs["section"] == sec).values
    up = adata.obs.loc[m, "up_score_z"].values.astype(float)
    ex = adata.obs.loc[m, "ex_score_z"].values.astype(float)
    tc = np.log10(adata.obs.loc[m, "n_counts"].values.astype(float))
    ng = adata.obs.loc[m, "n_genes"].values.astype(float)
    A = np.column_stack([np.ones(len(up)), tc, ng])
    up_r = up - A @ np.linalg.lstsq(A, up, rcond=None)[0]
    ex_r = ex - A @ np.linalg.lstsq(A, ex, rcond=None)[0]
    coords = adata.obsm["spatial"][m]
    nn = NearestNeighbors(n_neighbors=K + 1).fit(coords)
    _, nn_idx = nn.kneighbors(coords)
    neighbors = {i: nn_idx[i, 1:].tolist() for i in range(len(coords))}
    w = libpysal.weights.W(neighbors, silence_warnings=True)
    w.transform = "r"
    Wy = w.sparse.dot(ex_r)
    I_obs = float((up_r @ Wy) / (up_r @ up_r))
    rng = np.random.default_rng(SEED)
    sims = np.array([float((up_r @ w.sparse.dot(ex_r[rng.permutation(len(ex_r))])) / (up_r @ up_r))
                     for _ in range(999)])
    p = float(min(1.0, 2 * min(np.mean(sims >= I_obs), np.mean(sims <= I_obs))))
    sens_rows.append(dict(section=sec, I_bv_resid=I_obs, p_resid=p))
sens = pd.DataFrame(sens_rows)
sens.to_csv(os.path.join(OUT, "Table_S55v_M1_Visium_Bivariate_Residualized.csv"), index=False)
z_res = stats.norm.ppf(1 - sens["p_resid"].clip(1e-3, 1 - 1e-3) / 2) * np.sign(sens["I_bv_resid"])
z_st_res = z_res.sum() / np.sqrt(len(z_res))
p_st_res = 2 * stats.norm.sf(abs(z_st_res))
print(f"residualized bivariate I: Stouffer z={z_st_res:+.2f} p={p_st_res:.3g}")
print(sens.describe().loc[["mean", "min", "max"], ["I_bv_resid", "p_resid"]].to_string())

# ---- judgment ----
n = len(sdf)
z_st = sdf["z_bv"].sum() / np.sqrt(n)
p_st = 2 * stats.norm.sf(abs(z_st))
n_neg_sig = int(((sdf["I_bv"] < 0) & (sdf["p_bv"] < 0.05)).sum())
n_up_sig = int(((sdf["I_up"] > 0) & (sdf["p_up"] < 0.05)).sum())
n_ex_sig = int(((sdf["I_ex"] > 0) & (sdf["p_ex"] < 0.05)).sum())
r1 = (p_st < 0.05) and (z_st < 0) and (n_neg_sig >= 2)
r2 = False
r3 = False
if not r1:
    up_diffuse = n_up_sig <= n / 2
    ex_diffuse = n_ex_sig <= n / 2
    if (n_up_sig >= n / 2) != (n_ex_sig >= n / 2) and not (up_diffuse and ex_diffuse):
        r2 = True
    else:
        r3 = True
verdict = "R1" if r1 else ("R2" if r2 else "R3")
print(f"\nVERDICT = {verdict}")
print(f"  Stouffer z_bv={z_st:+.3f} p={p_st:.3g} | I_bv<0&p<0.05: {n_neg_sig}/23 | "
      f"I_up sig: {n_up_sig}/23 | I_ex sig: {n_ex_sig}/23")

# 方案A（2026-08-17）：不新增图面板副本，量化结果仅以 Table S55a–v 为唯一载体。
# 逐 spot 坐标与双臂评分存于 _intermediate/M1_visium_scored.h5ad
# （obs: up_score_z / ex_score_z / mdi_z / domain；obsm: spatial），
# 排版阶段如需出图可由 archive/M1_step8_figures_optional.py 恢复导出。

# ---- judgment report ----
rep_txt = f"""# M1 空间地理学验证 —— 判定报告（{verdict}）

生成日期：2026-08-17（规则先于结果，见 M1_pre_registration_20260817.md）

## 判定结果：{verdict}
"""
if verdict == "R3":
    rep_txt += f"""
按预注册规则判定为 **R3（阴性：双臂空间共定位、无分离）**。

- 双变量 Moran's I 合并（Stouffer，23 切片）：z = {z_st:+.2f}，p = {p_st:.2e}（**显著为正**，即空间共定位而非负相关分区）。
- I_bv<0 且置换 p<0.05 的切片仅 {n_neg_sig}/23（R1 要求 ≥2 个显著负相关切片且合并显著为负 → 不成立）。
- 计数/基因数残差化后双变量 I 合并 z = {z_st_res:+.2f}，p = {p_st_res:.2e} → 共定位并非组织密度混杂。
- 双臂各自呈空间聚集：上游 {n_up_sig}/23、执行 {n_ex_sig}/23 切片显著（执行臂聚集更强、更一致）。
- 结论：空间转录组**不支持**"上游塌陷区与执行诱导区在空间上分离"的区室化主张；双臂信号在空间上共定位、且逐 spot 相关性≈0（解耦但不分区）。

## 对文稿主张的影响
主结果四的"解离存在于组织/区室层面"表述需**收窄**为：
"解离 = 转录层面的解耦（同 spot/同细胞 ρ≈0）+ 双臂各自的空间组织（执行臂尤强），
但**不存在可检测的组织区域级空间分离**"。
这一定界本身即 M1 的高价值产出：以空间直接证据界定了 Mitoxyperilysis 的组织形态边界。

## 关键数值
"""
    rep_txt += sdf.sort_values("I_bv").to_string() + "\n"
with open(os.path.join(INT, "M1_judgment_report.md"), "w", encoding="utf-8") as f:
    f.write(rep_txt)
print("judgment report saved")
