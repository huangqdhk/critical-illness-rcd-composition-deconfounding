# -*- coding: utf-8 -*-
"""
mdi_lib.py — M2 解离指数（MDI）分析共享库（score_version = mdi_v1.0）
====================================================================
预注册：M2_pre_registration_20260817.md（判定规则先于结果）
- UCS / EIS / MDI 定义冻结：MDI = z(执行臂均值) - z(上游臂均值)
- 双臂基因唯一来源：Mitoxyperilysis_Gene_Manifest_v1.0.csv 的 arm 列
- 平台探针→基因塌陷：GEO 官方 GPL annot 文件，max-mean（每基因最多取整体均值最高的 3 个探针，逐样本取均值）
- 全部 z 标准化在队列内用全部样本（不分组），防止组内泄漏
"""
import gzip
import io
import math
import os
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from scipy import stats

try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
MANIFEST = ROOT + r"\04_AUDIT_GOVERNANCE\Mitoxyperilysis_Gene_Manifest_v1.0.csv"
GENE_SET_VERSION = "Mitoxy-80_v1.0"
SCORE_VERSION = "mdi_v1.0"
MT_GENES = {"MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CYB", "MT-ND1", "MT-ND2"}

# ---------------------------------------------------------------- manifest
def load_manifest():
    df = pd.read_csv(MANIFEST)
    arms = {}
    for arm in ("upstream_collapse", "execution_induction"):
        sub = df[df["arm"] == arm]
        genes = sorted(sub["hgnc_symbol"].astype(str).str.strip().str.upper().tolist())
        arms[arm] = genes
    # 别名表：gene -> 全部可匹配写法（大写）
    alias = defaultdict(set)
    for _, r in df.iterrows():
        key = str(r["hgnc_symbol"]).strip().upper()
        s = {str(r["gene_symbol"]).strip().upper(), key}
        for a in (str(r.get("aliases") or "").split(";")):
            a = a.strip()
            if a:
                s.add(a.upper())
        alias[key] = s
    return df, arms, alias


# ------------------------------------------------------- series matrix parser
def parse_series_matrix(path):
    """返回 (meta: dict tag->list[list], samples: list[str], mat: DataFrame 行=samples 列=probe)。
    用 csv.reader 处理引号与不等长行；非数值格置 NaN。"""
    import csv
    meta = defaultdict(list)
    samples = None
    rows = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        rd = csv.reader(f, delimiter="\t")
        in_table = False
        header_pending = False
        n_expect = None
        for parts in rd:
            if not parts:
                continue
            if not in_table:
                tag = parts[0]
                if tag == "!series_matrix_table_begin":
                    in_table = True
                    header_pending = True
                    continue
                if tag.startswith("!"):
                    meta[tag].append([p.strip('"') for p in parts[1:]])
                    continue
                if tag.strip('"') == "ID_REF":  # 无 marker 的兼容路径
                    samples = [c.strip('"') for c in parts[1:]]
                    n_expect = len(parts) - 1
                    in_table = True
                continue
            if header_pending:  # marker 之后的首行 = 表头
                header_pending = False
                samples = [c.strip('"') for c in parts[1:]]
                n_expect = len(parts) - 1
                continue
            if parts[0].strip('"').startswith("!series_matrix_table_end"):
                break
            vals = parts[1:]
            if n_expect and len(vals) < n_expect:
                vals = vals + [""] * (n_expect - len(vals))
            elif n_expect is None:
                n_expect = len(vals)
            else:
                vals = vals[:n_expect]
            rows.append((parts[0].strip('"'), vals))
    if samples is None:
        raise RuntimeError(f"series matrix 未找到 ID_REF 表头: {path}")
    data = {}
    for pid, vals in rows:
        data[pid] = pd.to_numeric(pd.Series([v.strip('"') for v in vals]), errors="coerce").values
    mat = pd.DataFrame(data)          # 行=sample, 列=probe
    mat.index = samples
    mat.index.name = "sample"
    mat.columns.name = "ID_REF"
    mat = mat.dropna(axis=1, how="all")
    return meta, samples, mat


def sample_characteristics(meta, key="!Sample_characteristics_ch1", label=None):
    """提取 characteristics 块中指定 label 的字段（如 'disease state'），返回 [samples]."""
    out = [None] * 0
    for block in meta.get(key, []):
        if label is None or (block and block[0].lower().startswith(label.lower())):
            out = [b.split(":", 1)[1].strip() if ":" in b else b.strip() for b in block]
            break
    return out


def sample_meta_dict(meta, tag="!Sample_title"):
    """返回 tag 的第一个块（逐样本字符串列表）。"""
    for block in meta.get(tag, []):
        if block:
            return block
    return []


# ------------------------------------------------- probe collapse via GPL annot
def load_probe_symbol_csv(path, wanted_aliases):
    """从 CSV(probe_id, symbol) 加载探针→canonical 基因映射（Bioconductor db 导出）。"""
    df = pd.read_csv(path)
    cols = [c for c in df.columns if c.lower() in ("probe_id", "probeid", "id")]
    syms = [c for c in df.columns if c.lower() in ("symbol", "gene_symbol")]
    if not cols or not syms:
        raise RuntimeError(f"probe-symbol CSV 列异常: {df.columns.tolist()}")
    probe2gene = defaultdict(set)
    for pid, sym in zip(df[cols[0]].astype(str), df[syms[0]].astype(str)):
        s = sym.strip().upper()
        for gene, aliases in wanted_aliases.items():
            if s in aliases:
                probe2gene[pid.strip()].add(gene)
    return probe2gene


def load_probe_map(annot_path, csv_path, wanted_aliases):
    """优先 GEO annot.gz，其次 Bioconductor CSV。返回 probe->genes。"""
    if annot_path and os.path.exists(annot_path):
        return load_gpl_annot(annot_path, wanted_aliases)
    if csv_path and os.path.exists(csv_path):
        return load_probe_symbol_csv(csv_path, wanted_aliases)
    return None


def load_gpl_annot(path, wanted_aliases):
    """解析 GEO GPL annot 文件，返回 probe -> 命中的 canonical 基因集合。
    wanted_aliases: canonical_gene -> set(大写别名)。"""
    probe2gene = defaultdict(set)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        in_table = False
        header = None
        for line in f:
            line = line.rstrip("\n")
            if not in_table:
                if line.startswith("!platform_table_begin"):
                    in_table = True
                continue
            if header is None:
                header = line.split("\t")
                # 找符号列
                sym_cols = []
                for i, c in enumerate(header):
                    cl = c.lower().strip('"')
                    if cl in ("gene symbol", "gene_symbol", "symbol", "genesymbol", "gene symbol(s)"):
                        sym_cols.append(i)
                if not sym_cols:
                    raise RuntimeError(f"annot 无 Gene Symbol 列: {header[:15]}")
                continue
            parts = line.split("\t")
            if len(parts) <= max(sym_cols):
                continue
            pid = parts[0].strip('"')
            if not pid or pid.startswith("!"):
                continue
            for i in sym_cols:
                for tok in parts[i].replace("///", ";").split(";"):
                    tok = tok.strip().strip('"').upper()
                    if not tok or tok == "---":
                        continue
                    for gene, aliases in wanted_aliases.items():
                        if tok in aliases:
                            probe2gene[pid].add(gene)
    return probe2gene


def collapse_probes_to_genes(mat, probe2gene, genes, max_probes=3):
    """mat: samples × probes(列名=probe id)。max-mean 塌陷：每基因取整体均值最高的 ≤3 个探针，逐样本均值。"""
    out = pd.DataFrame(index=mat.index)
    for g in genes:
        probes = [p for p in mat.columns if g in probe2gene.get(p, set())]
        if not probes:
            continue
        sub = mat[probes]
        if len(probes) > max_probes:
            best = sub.mean(axis=0).nlargest(max_probes).index
            sub = sub[best]
        out[g] = sub.mean(axis=1)
    return out


# ------------------------------------------------------------------ scoring
def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = x.std(ddof=1)
    return (x - x.mean()) / sd if sd > 0 else np.zeros_like(x)


def arm_scores(mat_gene, arms):
    """mat_gene: samples × genes(log2, 列名=大写符号)。返回 DataFrame[sample] UCS/EIS/MDI(z) + n_up/n_ex + UCS_nomt/EIS_nomt/MDI_nomt + ssGSEA 敏感性。"""
    up, ex = arms["upstream_collapse"], arms["execution_induction"]
    have = {g.upper() for g in mat_gene.columns}
    up_in = [g for g in up if g in have]
    ex_in = [g for g in ex if g in have]
    up_nomt = [g for g in up_in if g not in MT_GENES]
    res = pd.DataFrame(index=mat_gene.index)
    res["n_up"] = len(up_in)
    res["n_ex"] = len(ex_in)
    res["UCS_raw"] = mat_gene[up_in].mean(axis=1)
    res["EIS_raw"] = mat_gene[ex_in].mean(axis=1)
    res["UCS"] = zscore(res["UCS_raw"])
    res["EIS"] = zscore(res["EIS_raw"])
    res["MDI"] = res["EIS"] - res["UCS"]
    if up_nomt:
        res["UCS_raw_nomt"] = mat_gene[up_nomt].mean(axis=1)
        res["UCS_nomt"] = zscore(res["UCS_raw_nomt"])
        res["MDI_nomt"] = res["EIS"] - res["UCS_nomt"]
    # running-sum ssGSEA（Barbie 2009，β=0.25，逐臂）
    for nm, gs in (("UCS_ssgsea", up_in), ("EIS_ssgsea", ex_in)):
        if len(gs) >= 2:
            res[nm] = ssgsea_per_sample(mat_gene, gs, beta=0.25)
    return res


def ssgsea_per_sample(mat_gene, gene_set, beta=0.25):
    """经典 running-sum ssGSEA：逐样本基因排序，ES = 加权 KS 差。返回逐样本 ES 向量。"""
    gs = set(gene_set)
    scores = np.zeros(mat_gene.shape[0])
    for i, (_, row) in enumerate(mat_gene.iterrows()):
        r = row.values
        order = np.argsort(-r)  # 高表达在前
        n = len(order)
        p_hit = np.zeros(n)
        p_miss = np.zeros(n)
        hits = np.array([1.0 if (mat_gene.columns[j].upper() in gs) else 0.0 for j in order])
        w = np.abs(r[order]) ** beta
        wsum = (hits * w).sum()
        if wsum == 0:
            scores[i] = 0.0
            continue
        miss_w = (1.0 - hits)
        p_hit = np.cumsum(hits * w) / wsum
        denom = miss_w.sum()
        p_miss = np.cumsum(miss_w) / denom if denom > 0 else np.zeros(n)
        scores[i] = (p_hit - p_miss).max() - (p_hit - p_miss).min() * 0.0
        # 标准定义 ES = max(hit - miss)，对称版取 max(abs)。这里用 max(0, max) - min? 采用 Barbie: ES = max_i (P_hit_i - P_miss_i)
        scores[i] = (p_hit - p_miss).max()
    return scores


# ------------------------------------------------------------------- stats
def mwu_cliff(a, b):
    """Mann-Whitney U 双侧 p + Cliff's delta (b 相对 a)。"""
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    n, m = len(a), len(b)
    # Cliff's delta = P(b > a) - P(b < a)
    d = 0.0
    for x in a:
        d += (np.sum(b > x) - np.sum(b < x))
    delta = d / (n * m)
    return p, delta


def hedges_g(a, b):
    n1, n2 = len(a), len(b)
    s1, s2 = np.var(a, ddof=1), np.var(b, ddof=1)
    sp = math.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    g = (np.mean(b) - np.mean(a)) / sp
    # Hedges 校正
    corr = 1 - 3 / (4 * (n1 + n2) - 9)
    g = g * corr
    se = math.sqrt(1 / n1 + 1 / n2 + g * g / (2 * (n1 + n2)))
    return g, se


def der_simonian_laird(gs, ses):
    """随机效应 meta（Hedges' g + SE）。返回 pooled g, se, p, tau2, I2, Q, Q_df, Q_p。"""
    gs = np.asarray(gs, dtype=float)
    ses = np.asarray(ses, dtype=float)
    w = 1.0 / ses ** 2
    g_fe = (w * gs).sum() / w.sum()
    Q = (w * (gs - g_fe) ** 2).sum()
    k = len(gs)
    df_q = k - 1
    c = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (Q - df_q) / c) if df_q > 0 else 0.0
    w_re = 1.0 / (ses ** 2 + tau2)
    g_re = (w_re * gs).sum() / w_re.sum()
    se_re = math.sqrt(1.0 / w_re.sum())
    z = g_re / se_re
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    I2 = max(0.0, (Q - df_q) / Q) if Q > 0 else 0.0
    Qp = stats.chi2.sf(Q, df_q) if df_q > 0 else 1.0
    return g_re, se_re, p, tau2, I2, Q, df_q, Qp


def stouffer_z(ps):
    """Stouffer 合并 p（单侧方向由调用方保证）。"""
    ps = np.asarray(ps, dtype=float)
    z = np.sum(stats.norm.ppf(1 - ps)) / math.sqrt(len(ps))
    return z, 1 - stats.norm.cdf(z)


def bh(pvals):
    """BH 校正（2026-08-27 修复：原实现当最小 p 不在首位时 q 值错位，见 03_LOGS/M11M12_step4_log 审计注）。"""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q_sorted = np.minimum(1, p[order] * n / (np.arange(n) + 1))
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q = np.empty(n)
    q[order] = q_sorted
    return q


def clf_delta_ci(a, b, n_boot=2000, seed=0):
    rng = np.random.default_rng(seed)
    _, obs = mwu_cliff(a, b)
    boots = []
    for _ in range(n_boot):
        a1 = rng.choice(a, len(a), replace=True)
        b1 = rng.choice(b, len(b), replace=True)
        _, d = mwu_cliff(a1, b1)
        boots.append(d)
    return obs, np.percentile(boots, 2.5), np.percentile(boots, 97.5)


def cochran_armitage(table):
    """table: 2×k（行=结局 0/1，列=有序分组）。返回 z, p（双侧）。"""
    table = np.asarray(table, dtype=float)
    n = table.sum()
    k = table.shape[1]
    scores = np.arange(k)
    p1 = table[1].sum() / n
    n_col = table.sum(axis=0)
    p_col = table[1] / n_col
    num = (n_col * (p_col - p1) * scores).sum()
    den = p1 * (1 - p1) * (n_col * (scores - (n_col * scores).sum() / n) ** 2).sum()
    z = num / math.sqrt(den)
    return z, 2 * (1 - stats.norm.cdf(abs(z)))


# ------------------------------------------------------------ 校准与 DCA
def hosmer_lemeshow(y, phat, g=10):
    """Hosmer-Lemeshow C 统计量与 p（g=10 组）。"""
    order = np.argsort(phat)
    bins = np.array_split(order, g)
    obs = np.array([y[b].sum() for b in bins], dtype=float)
    exp = np.array([phat[b].sum() for b in bins], dtype=float)
    n = np.array([len(b) for b in bins], dtype=float)
    c = ((obs - exp) ** 2 / (exp * (1 - exp / n) + 1e-12)).sum()
    p = stats.chi2.sf(c, g - 2)
    return c, p


def calibration_curve_stats(y, phat):
    """logistic 校准：logit(phat) ~ y 的斜率与截距（预测模型校准用）。"""
    eps = 1e-6
    lo = np.log(np.clip(phat, eps, 1 - eps) / np.clip(1 - phat, eps, 1 - eps))
    X = np.column_stack([np.ones_like(lo), lo])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta[0], beta[1]  # 截距, 斜率


def dca_net_benefit(y, phat, thresholds=(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)):
    """决策曲线净获益。返回 {threshold: (NB_model, NB_all, NB_none)}。"""
    out = {}
    y = np.asarray(y)
    for t in thresholds:
        pred_pos = phat >= t
        tp = ((pred_pos) & (y == 1)).sum()
        fp = ((pred_pos) & (y == 0)).sum()
        n = len(y)
        nb_model = tp / n - fp / n * (t / (1 - t))
        nb_all = y.sum() / n - (1 - y).sum() / n * (t / (1 - t))
        nb_none = 0.0
        out[t] = (nb_model, nb_all, nb_none)
    return out


def brier_score(y, phat):
    return float(np.mean((y - phat) ** 2))


# --------------------------------------------------------------------- IGP
def in_group_proportion(X, labels, k=10):
    """IGP (Kapp & Tibshirani 2007)：每样本 k 近邻中同类比例的平均。X: n×2 特征矩阵。"""
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=k + 1).fit(X)
    _, idx = nn.kneighbors(X)
    labs = np.asarray(labels)
    fracs = []
    for i in range(len(labs)):
        neigh = labs[idx[i, 1:]]
        fracs.append(np.mean(neigh == labs[i]))
    return float(np.mean(fracs))


def igp_permutation_pvalue(X, labels, k=10, n_perm=500, seed=0):
    rng = np.random.default_rng(seed)
    obs = in_group_proportion(X, labels, k)
    cnt = 0
    labs = np.asarray(labels)
    for _ in range(n_perm):
        perm = rng.permutation(labs)
        if in_group_proportion(X, perm, k) >= obs:
            cnt += 1
    return obs, (cnt + 1) / (n_perm + 1)
