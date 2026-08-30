# -*- coding: utf-8 -*-
"""
composition.py — 组成预测值 / 组成残差 MDI（M15 工具核心之二）
================================================================
冻结口径（与 M10A_step1 / M11M12_step1 完全一致）：
  - 参考谱：
      Monaco = GSE107011 29 种纯化免疫细胞 RNA-seq（TPM，含中性粒细胞），
               log2(TPM+1) 细胞型均值（本地冻结缓存 _intermediate/M10A_monaco_ct_means.csv，
               由 00_RAW_DATA/GSE107011 原始 TPM 生成；发布版经 Zenodo DOI 分发）。
      ABIS   = sigmatrixRNAseq.txt（17 细胞型，全血主用）与 sigmatrixMicro.txt（11 细胞型，
               微阵列跨平台敏感性口径），log2(+1)，基因符号去引号、重复符号取均值。
  - 反卷积：NNLS（scipy.optimize.nnls，经 cholesky 预分解；主口径）；
            方法交叉 = OLS 后置零截断（OLS-CLS；离线环境无法运行 CIBERSORTx，如实披露）。
  - 合成谱：X̂_comp(s) = Σ_c π_c(s)·μ_c（log2 空间）；臂均值 -> 队列内全样本 z（不分组）
            -> UCS_comp / EIS_comp / MDI_comp（mdi_v1.0 同口径）。
  - 组成残差：逐队列 OLS MDI ~ MDI_comp -> MDI_resid（b0/b1 一并输出）。
"""
import os

import numpy as np
import pandas as pd
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import nnls

from . import ROOT

RAW = os.path.join(ROOT, "00_RAW_DATA")
INTER = os.path.join(ROOT, "_intermediate")

MONACO_TPM = os.path.join(RAW, "GSE107011_Monaco_ImmuneRef", "GSE107011_Processed_data_TPM.txt.gz")
MONACO_CACHE = os.path.join(INTER, "M10A_monaco_ct_means.csv")
ABIS_RNA = os.path.join(RAW, "ABIS", "sigmatrixRNAseq.txt")
ABIS_MICRO = os.path.join(RAW, "ABIS", "sigmatrixMicro.txt")


def zc(s):
    s = np.asarray(s, float)
    sd = s.std(ddof=1)
    return (s - s.mean()) / sd if sd > 0 else np.zeros_like(s)


# ------------------------------------------------------------------ references
def load_monaco_ref(cache_path=None):
    """Monaco 29 细胞型参考谱（log2(TPM+1)，细胞型 x 基因）。主用冻结缓存。"""
    path = cache_path or MONACO_CACHE
    mu = pd.read_csv(path, index_col=0)
    mu.index = [str(i).strip().upper() for i in mu.index]
    mu.columns = [str(c).strip().upper() for c in mu.columns]
    return mu


def load_abis_refs(raw_dir=None):
    """返回 (abis_rna, abis_micro)：log2(+1) 细胞型 x 基因。基因符号去引号、重复取均值。"""
    base = raw_dir or RAW
    out = {}
    for fname, tag in [(r"ABIS\sigmatrixRNAseq.txt", "rna"),
                       (r"ABIS\sigmatrixMicro.txt", "micro")]:
        df = pd.read_csv(os.path.join(base, fname), sep="\t", index_col=0)
        df.index = [str(i).strip().strip('"\'').upper() for i in df.index]
        df.columns = [str(c).strip().strip('"\'') for c in df.columns]
        df = df.groupby(level=0).mean()
        out[tag] = np.log2(df + 1.0).T
    return out["rna"], out["micro"]


def ref_arm_coverage(ref, up, ex):
    return (len([g for g in up if g in ref.columns]),
            len([g for g in ex if g in ref.columns]))


# ------------------------------------------------------------------ deconvolution
def deconvolve_nnls(X, ref_log):
    """X: samples x genes(log2)；ref_log: celltypes x genes(log2)。返回 (proportions, common)。"""
    common = [g for g in X.columns if g in ref_log.columns]
    A = ref_log[common].values.T        # genes x celltypes
    Xc = X[common].values
    AtA = A.T @ A
    n_ct = A.shape[1]
    Lc = cholesky(AtA + np.eye(n_ct) * 1e-10, lower=True)
    rows = []
    for i in range(Xc.shape[0]):
        rhs = solve_triangular(Lc, A.T @ Xc[i], lower=True)
        pi, _ = nnls(Lc.T, rhs)
        s = pi.sum()
        rows.append((pi / s) if s > 0 else np.full(n_ct, 1.0 / n_ct))
    return pd.DataFrame(rows, index=X.index, columns=ref_log.index), common


def deconvolve_olcls(X, ref_log):
    """OLS 后置零截断（OLS-CLS）交叉口径。返回 (proportions, common)。"""
    common = [g for g in X.columns if g in ref_log.columns]
    A = ref_log[common].values.T
    Xc = X[common].values
    AtA = A.T @ A
    n_ct = A.shape[1]
    rows = []
    for i in range(Xc.shape[0]):
        pi = np.linalg.lstsq(AtA + np.eye(n_ct) * 1e-6, A.T @ Xc[i], rcond=None)[0]
        pi = np.clip(pi, 0, None)
        s = pi.sum()
        rows.append((pi / s) if s > 0 else np.full(n_ct, 1.0 / n_ct))
    return pd.DataFrame(rows, index=X.index, columns=ref_log.index), common


# ------------------------------------------------------------------ synthetic scores
def synthetic_scores(prop, ref_log, up_genes, ex_genes):
    """X_comp = pi @ mu；臂均值（log2 空间）-> 队列内 z -> UCS_comp/EIS_comp/MDI_comp。"""
    up_in = [g for g in up_genes if g in ref_log.columns]
    ex_in = [g for g in ex_genes if g in ref_log.columns]
    mu_up = ref_log[up_in].values          # celltypes x genes
    mu_ex = ref_log[ex_in].values
    P = prop.values                        # samples x celltypes
    uc = P @ mu_up.mean(axis=1)
    ec = P @ mu_ex.mean(axis=1)
    out = pd.DataFrame(index=prop.index)
    out["UCS_comp_raw"] = uc
    out["EIS_comp_raw"] = ec
    out["UCS_comp"] = zc(uc)
    out["EIS_comp"] = zc(ec)
    out["MDI_comp"] = out["EIS_comp"] - out["UCS_comp"]
    out["n_up_ref"] = len(up_in)
    out["n_ex_ref"] = len(ex_in)
    return out


def residual_mdi(mdi_obs, mdi_comp):
    """逐队列 OLS MDI ~ MDI_comp 的残差。返回 (MDI_resid, b0, b1)。"""
    A = np.column_stack([np.ones(len(mdi_obs)), np.asarray(mdi_comp, float)])
    beta, *_ = np.linalg.lstsq(A, np.asarray(mdi_obs, float), rcond=None)
    resid = np.asarray(mdi_obs, float) - A @ beta
    return resid, beta[0], beta[1]


def score_and_compose(mat_gene, arms, monaco_ref, ref_main=None):
    """工具一站式核心：mdi_v1.0 评分 + 组成预测 + 组成残差（Monaco 主口径）。
    mat_gene: samples x genes(log2, 大写符号)。
    - 评分：arm_scores（mdi_v1.0）
    - 比例：Monaco NNLS（主口径）；ref_main 给定时额外 ABIS NNLS 交叉
    - MDI_comp：Monaco 合成谱（mdi_v1.0 同口径）；ref_main 给定时同报 ABIS 合成（_abis 后缀）
    - 细胞比例汇总：Neut/Mono/comp（Monaco 与 ABIS 各一套，后者 _abis 后缀）
    - MDI_resid：队列内 OLS MDI ~ MDI_comp（b0/b1 一并输出）
    返回 (scores_df, prop_monaco, prop_main_or_None)。"""
    from .mdi_core import arm_scores
    sc = arm_scores(mat_gene, arms)
    prop_mon, common_mon = deconvolve_nnls(mat_gene, monaco_ref)
    sy = synthetic_scores(prop_mon, monaco_ref,
                          arms["upstream_collapse"], arms["execution_induction"])
    out = pd.concat([sc, sy], axis=1)
    nm_mon = neutrophil_monocyte(prop_mon).rename(columns=lambda c: c + "_monaco")
    out = pd.concat([out, nm_mon], axis=1)
    prop_main = None
    if ref_main is not None:
        prop_main, common_main = deconvolve_nnls(mat_gene, ref_main)
        nm = neutrophil_monocyte(prop_main)
        out = pd.concat([out, nm], axis=1)
        sy_abis = synthetic_scores(prop_main, ref_main,
                                   arms["upstream_collapse"], arms["execution_induction"])
        sy_abis = sy_abis.rename(columns=lambda c: c + "_abis")
        out = pd.concat([out, sy_abis], axis=1)
        out["n_common_abis"] = len(common_main)
    resid, b0, b1 = residual_mdi(out["MDI"].values, out["MDI_comp"].values)
    out["MDI_resid"] = resid
    out["b0"] = b0
    out["b1"] = b1
    out["n_common_monaco"] = len(common_mon)
    return out, prop_mon, prop_main


def neutrophil_monocyte(prop):
    """返回 DataFrame: Neut / Mono / comp（Neut+Mono 合计比例）。"""
    out = pd.DataFrame(index=prop.index)
    neut_cols = [c for c in prop.columns if "neutrophil" in c.lower()]
    mono_cols = [c for c in prop.columns if c.lower() in ("monocytes", "monocytes c", "monocytes nc+i")]
    out["Neut"] = prop[neut_cols].sum(axis=1) if neut_cols else np.nan
    out["Mono"] = prop[mono_cols].sum(axis=1) if mono_cols else np.nan
    out["comp"] = out["Neut"] + out["Mono"]
    return out
