# -*- coding: utf-8 -*-
"""
mdi_core.py — 两臂分解 + mdi_v1.0 评分核心（自包含冻结副本）
================================================================
本模块为项目根 mdi_lib.py 的发布快照（同版本、同数值口径），供 M15 工具包
独立发布使用。项目内运行时与 mdi_lib.py 的等价性由 tests/test_regression_11cohorts.py
对 11 个队列的逐样本核对证明（max|Δ| 报告于回归表）。
冻结定义（mdi_v1.0）：
  - MDI = z(执行臂均值) - z(上游臂均值)；z 在队列内全样本（不分组）标准化
  - 微阵列探针→基因塌陷：GEO 官方 GPL annot 或 Bioconductor db 导出的
    probe→symbol CSV，别名匹配（manifest 别名表），max-mean（每基因取整体均值
    最高的 ≤3 个探针，逐样本取均值）
  - 敏感性：UCS_nomt/MDI_nomt（去 6 个 MT-* 基因）、UCS_ssgsea/EIS_ssgsea（β=0.25）
"""
import gzip
import math
import os
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

from . import GENE_SET_VERSION, SCORE_VERSION, ROOT

# 冻结 manifest 路径（唯一基因集来源）：优先项目治理目录，其次工具包内置副本
CANONICAL_MANIFEST = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "Mitoxyperilysis_Gene_Manifest_v1.0.csv")
BUNDLED_MANIFEST = os.path.join(os.path.dirname(__file__), "..", "data",
                                "Mitoxyperilysis_Gene_Manifest_v1.0.csv")

MT_GENES = {"MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CYB", "MT-ND1", "MT-ND2"}


# ------------------------------------------------------------------ manifest
def load_manifest(manifest_path=None):
    """返回 (df, arms, alias)。arms: {'upstream_collapse': [genes], 'execution_induction': [genes]}；
    alias: canonical_gene -> set(大写别名写法)。"""
    path = manifest_path or (CANONICAL_MANIFEST if os.path.exists(CANONICAL_MANIFEST)
                             else os.path.abspath(BUNDLED_MANIFEST))
    df = pd.read_csv(path)
    arms = {}
    for arm in ("upstream_collapse", "execution_induction"):
        sub = df[df["arm"] == arm]
        genes = sorted(sub["hgnc_symbol"].astype(str).str.strip().str.upper().tolist())
        arms[arm] = genes
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
    """提取 characteristics 块中指定 label 的字段（如 'disease'），返回 [samples]."""
    out = []
    for block in meta.get(key, []):
        if label is None or (block and block[0].lower().startswith(label.lower())):
            out = [b.split(":", 1)[1].strip() if ":" in b else b.strip() for b in block]
            break
    return out


def sample_meta_dict(meta, tag="!Sample_title"):
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
    """解析 GEO GPL annot 文件，返回 probe -> 命中的 canonical 基因集合。"""
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


def collapse_probes_all(probe_mat, probe2sym, max_probes=3):
    """全探针塌陷：probe2sym: dict probe->symbol(大写)。返回 samples × genes。"""
    g2p = defaultdict(list)
    for p, s in probe2sym.items():
        if p in probe_mat.columns:
            g2p[s].append(p)
    out = pd.DataFrame(index=probe_mat.index)
    for g, ps in g2p.items():
        sub = probe_mat[ps]
        if len(ps) > max_probes:
            best = sub.mean(axis=0).nlargest(max_probes).index
            sub = sub[best]
        out[g] = sub.mean(axis=1)
    return out


def gene_matrix_from_symbols(mat, alias, max_dups="mean"):
    """符号矩阵→canonical 基因矩阵（RNA-seq 口径）：
    mat: samples × genes(列名=任意符号写法)；按 manifest 别名表解析到 canonical 臂基因，
    重复 canonical 取逐样本均值。返回 samples × canonical genes。"""
    canon = defaultdict(list)
    for c in mat.columns:
        cu = str(c).strip().upper()
        for gene, aliases in alias.items():
            if cu in aliases:
                canon[gene].append(c)
                break
    out = pd.DataFrame(index=mat.index)
    for g, cols in canon.items():
        out[g] = mat[cols].mean(axis=1)
    return out


# ------------------------------------------------------------------ scoring
def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = x.std(ddof=1)
    return (x - x.mean()) / sd if sd > 0 else np.zeros_like(x)


def arm_scores(mat_gene, arms):
    """mat_gene: samples × genes(log2, 列名=大写符号)。返回 DataFrame[sample] UCS/EIS/MDI(z)
    + n_up/n_ex + UCS_nomt/EIS_nomt/MDI_nomt + ssGSEA 敏感性（mdi_v1.0 冻结口径）。"""
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
    for nm, gs in (("UCS_ssgsea", up_in), ("EIS_ssgsea", ex_in)):
        if len(gs) >= 2:
            res[nm] = ssgsea_per_sample(mat_gene, gs, beta=0.25)
    return res


def ssgsea_per_sample(mat_gene, gene_set, beta=0.25):
    """经典 running-sum ssGSEA：逐样本基因排序，ES = 加权 KS 差（Barbie 2009，β=0.25）。"""
    gs = set(gene_set)
    scores = np.zeros(mat_gene.shape[0])
    for i, (_, row) in enumerate(mat_gene.iterrows()):
        r = row.values
        order = np.argsort(-r)
        hits = np.array([1.0 if (mat_gene.columns[j].upper() in gs) else 0.0 for j in order])
        w = np.abs(r[order]) ** beta
        wsum = (hits * w).sum()
        if wsum == 0:
            scores[i] = 0.0
            continue
        miss_w = (1.0 - hits)
        p_hit = np.cumsum(hits * w) / wsum
        denom = miss_w.sum()
        p_miss = np.cumsum(miss_w) / denom if denom > 0 else np.zeros(len(order))
        scores[i] = (p_hit - p_miss).max()
    return scores


# ------------------------------------------------------------------- stats
def mwu_cliff(a, b):
    """Mann-Whitney U 双侧 p + Cliff's delta (b 相对 a)。"""
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    n, m = len(a), len(b)
    d = 0.0
    for x in a:
        d += (np.sum(b > x) - np.sum(b < x))
    delta = d / (n * m)
    return p, delta


def hedges_g(a, b):
    """Hedges' g（b 相对 a）与 SE。"""
    n1, n2 = len(a), len(b)
    s1, s2 = np.var(a, ddof=1), np.var(b, ddof=1)
    sp = math.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    g = (np.mean(b) - np.mean(a)) / sp
    corr = 1 - 3 / (4 * (n1 + n2) - 9)
    g = g * corr
    se = math.sqrt(1 / n1 + 1 / n2 + g * g / (2 * (n1 + n2)))
    return g, se


def bh(pvals):
    """BH 校正（2026-08-27 修复版：与项目根 mdi_lib 一致）。"""
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
