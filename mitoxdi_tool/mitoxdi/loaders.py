# -*- coding: utf-8 -*-
"""
loaders.py — 队列加载器（回归测试与采纳演示共用）
================================================================================
每个加载器精确复刻原流水线（M2_step1 / M10A_step1 / M11M12_step1）的同一口径，
供 tests/test_regression_11cohorts.py 与 demos/ 复用：
  - M2 口径：微阵列探针→基因用别名匹配（manifest 别名表）+ max-mean ≤3 塌陷
  - M10A/M11M12 口径：直接符号匹配全探针塌陷（GPL CSV/Bioconductor db 导出）
  - GSE215865：官方 logCPM（float32）→ ENSG(去版本)→HUGO（本地缓存 + mygene 补缺，
    缓存 = M11M12_step1_ensg2sym_union.csv）→ 受试者×日标签重复取均值
  - GSE54514/GSE106878：non-normalized → quantile（lumiN 等价）+ log2(x+1)
  - GSE157103（演示）：TPM → log2(TPM+1)，符号矩阵
  - GSE66099（演示）：GPL570 系列矩阵 → 直接符号全塌陷
"""
import gzip
import os
import re
import warnings
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

from . import ROOT
from .mdi_core import (parse_series_matrix, sample_characteristics,
                       load_probe_symbol_csv, collapse_probes_to_genes,
                       gene_matrix_from_symbols, load_manifest)

RAW = os.path.join(ROOT, "00_RAW_DATA")
INTER = os.path.join(ROOT, "_intermediate")
GOV = os.path.join(ROOT, "04_AUDIT_GOVERNANCE")

GPL23159_CSV = os.path.join(RAW, "GPL23159_probe2symbol_clariomsdb.csv")
GPL570_CSV = os.path.join(RAW, "GPL570_probe2symbol_hgu133plus2db.csv")
GPL6947_TXT = os.path.join(RAW, "GPL6947_table.txt")
GPL10295_TXT = os.path.join(RAW, "GPL10295_table.txt")


# ---------------------------------------------------------------- utilities
def read_probe2sym_csv(path):
    """Bioconductor CSV → dict probe->symbol(大写)。"""
    gpl = pd.read_csv(path)
    pc = [c for c in gpl.columns if c.lower() in ("probe_id", "probeid", "id")][0]
    sc = [c for c in gpl.columns if c.lower() in ("symbol", "gene_symbol")][0]
    return dict(zip(gpl[pc].astype(str).str.strip(),
                    gpl[sc].astype(str).str.strip().str.upper()))


def read_gpl_table(path):
    """GEO 平台表（!platform_table_begin 之后）→ dict probe->symbol(大写)。"""
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = f.read().split("\n")
    hdr_idx = next(i for i, ln in enumerate(lines) if ln.startswith("!platform_table_begin"))
    header = lines[hdr_idx + 1].split("\t")
    sym_col = next(i for i, h in enumerate(header) if h.strip().strip('"').lower() == "symbol")
    id_col = next(i for i, h in enumerate(header) if h.strip().strip('"').lower() in ("id",))
    p2s = {}
    for ln in lines[hdr_idx + 2:]:
        if ln.startswith("!platform_table_end") or not ln.strip():
            break
        parts = ln.split("\t")
        if len(parts) > max(sym_col, id_col):
            s = parts[sym_col].strip().strip('"').upper()
            pid = parts[id_col].strip().strip('"')
            if s and s != "---":
                p2s[pid] = s
    return p2s


def collapse_all(probe_mat, p2s, max_probes=3):
    """全探针 max-mean ≤3 塌陷（M11M12 collapse_maxmean 同口径）。"""
    g2p = defaultdict(list)
    for p, s in p2s.items():
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


def fast_collapse_wanted(probe_mat, probe2sym_csv, wanted_symbols, max_probes=3):
    """M10A fast_collapse 同口径：仅塌陷 wanted_symbols 内的符号。"""
    pmap = pd.read_csv(probe2sym_csv)
    pmap.columns = [str(c).upper().strip() for c in pmap.columns]
    pc = [c for c in pmap.columns if c in ("PROBEID", "ID", "PROBE_ID")][0]
    sc = [c for c in pmap.columns if c in ("SYMBOL", "GENE_SYMBOL", "GENE.SYMBOL")][0]
    pmap = pmap[[pc, sc]]
    pmap[pc] = pmap[pc].astype(str).str.strip()
    pmap[sc] = pmap[sc].astype(str).str.strip().str.upper()
    g2p = {}
    for pid, sym in pmap.itertuples(index=False, name=None):
        if sym in wanted_symbols:
            g2p.setdefault(sym, []).append(pid)
    out = pd.DataFrame(index=probe_mat.index)
    for g, probes in g2p.items():
        ps = [p for p in probes if p in probe_mat.columns]
        if not ps:
            continue
        sub = probe_mat[ps]
        if len(ps) > max_probes:
            best = sub.mean(axis=0).nlargest(max_probes).index
            sub = sub[best]
        out[g] = sub.mean(axis=1)
    return out


def quantile_normalize(mat):
    """列分位数归一化（等价 lumiN）。mat: samples x features。"""
    M = mat.values.astype(float).T   # features x samples
    order = np.argsort(M, axis=1)
    sortedM = np.take_along_axis(M, order, axis=1)
    means = np.nanmean(sortedM, axis=0)
    Q = np.broadcast_to(means, M.shape)
    out = np.empty_like(M)
    for i in range(M.shape[0]):
        out[i] = Q[i][np.argsort(order[i])]
    return pd.DataFrame(out.T, index=mat.index, columns=mat.columns)


def parse_gsm_metadata(path):
    txt = open(path, encoding="utf-8", errors="replace").read()
    recs = txt.split("^SAMPLE")
    out = []
    for r in recs[1:]:
        d = {}
        for ln in r.split("\n"):
            ln = ln.strip()
            if ln.startswith("!Sample_geo_accession"):
                d["gsm"] = ln.split("=", 1)[1].strip()
            elif ln.startswith("!Sample_title"):
                d["title"] = ln.split("=", 1)[1].strip()
            elif ln.startswith("!Sample_description"):
                d.setdefault("descriptions", []).append(ln.split("=", 1)[1].strip())
            elif ln.startswith("!Sample_characteristics_ch1"):
                kv = ln.split("=", 1)[1].strip()
                if ":" in kv:
                    k, v = kv.split(":", 1)
                    d[k.strip().lower()] = v.strip()
        out.append(d)
    return pd.DataFrame(out)


# ------------------------------------------------------------ M2 / M10A 队列
def load_gse185263_m2():
    """M2 口径：manifest ENSG→symbol（80 基因）→ CPM → log2(CPM+1)。"""
    man, _, _ = load_manifest()
    ensg2sym = {}
    for _, r in man.iterrows():
        e = str(r.get("ensembl_gene_id") or "").strip()
        if e and e.upper().startswith("ENSG"):
            ensg2sym[e.upper()] = str(r["hgnc_symbol"]).strip().upper()
    cnt = pd.read_csv(os.path.join(RAW, "GSE185263_Lung_ARDS", "GSE185263_raw_counts.csv"), index_col=0)
    cnt.index = [str(i).upper() for i in cnt.index]
    cnt.columns = [str(c) for c in cnt.columns]
    mapped = {e: s for e, s in ensg2sym.items() if e in cnt.index}
    sub = cnt.loc[list(mapped.keys())].copy()
    sub.index = [mapped[e] for e in sub.index]
    sub = sub.groupby(level=0).mean()
    cpm = sub / sub.sum(axis=0) * 1e6
    return np.log2(cpm + 1).T


def load_gse185263_m10a():
    """M10A 口径：P2 缓存 ENSG→symbol（全基因组）→ CPM → log2(CPM+1)。"""
    cnt = pd.read_csv(os.path.join(RAW, "GSE185263_Lung_ARDS", "GSE185263_raw_counts.csv"), index_col=0)
    cnt.index = [str(i).upper() for i in cnt.index]
    cnt.columns = [str(c) for c in cnt.columns]
    e2s = dict(pd.read_csv(os.path.join(INTER, "P2_gse185263_en2sym.csv")).values)
    e2s = {str(k).upper(): str(v).upper() for k, v in e2s.items()}
    sub = cnt.loc[[e for e in e2s if e in cnt.index]].copy()
    sub.index = [e2s[e] for e in sub.index]
    sub = sub.groupby(level=0).mean()
    cpm = sub / sub.sum(axis=0) * 1e6
    return np.log2(cpm + 1).T


def load_gse32707():
    m7 = pd.read_csv(os.path.join(RAW, "GSE32707_data", "GSE32707_gene_symbol_log2_matrix.csv.gz"),
                     index_col=0)
    m7.index = [str(i).upper() for i in m7.index]
    m7.columns = [str(c) for c in m7.columns]
    return m7.groupby(level=0).mean().T


def load_gse212865_m2(arms):
    """M2 口径：GPL23159 别名匹配臂基因塌陷。"""
    meta, s, mat = parse_series_matrix(
        os.path.join(RAW, "GSE212865_data", "GSE212865_series_matrix.txt.gz"))
    man, _, alias = load_manifest()
    probe2g = load_probe_symbol_csv(GPL23159_CSV, alias)
    return collapse_probes_to_genes(mat, probe2g, arms["upstream_collapse"] + arms["execution_induction"])


def load_gse212865_m11():
    """M11M12 口径：GPL23159 直接符号全塌陷。"""
    meta, s, mat = parse_series_matrix(
        os.path.join(RAW, "GSE212865_data", "GSE212865_series_matrix.txt.gz"))
    p2s = read_probe2sym_csv(GPL23159_CSV)
    return collapse_all(mat, p2s)


def load_gse212865_m10a(monaco_symbols):
    """M10A 口径：fast_collapse（仅 Monaco 符号；CSV 行迭代、保留重复探针行）。"""
    meta, s, mat = parse_series_matrix(
        os.path.join(RAW, "GSE212865_data", "GSE212865_series_matrix.txt.gz"))
    return fast_collapse_wanted(mat, GPL23159_CSV, set(monaco_symbols))


def load_gse188309_m2(arms):
    meta, s, mat = parse_series_matrix(
        os.path.join(RAW, "GSE188309_CAP_WholeBlood", "GSE188309_series_matrix.txt.gz"))
    man, _, alias = load_manifest()
    probe2g = load_probe_symbol_csv(GPL23159_CSV, alias)
    return collapse_probes_to_genes(mat, probe2g, arms["upstream_collapse"] + arms["execution_induction"])


def load_gse188309_m10a(monaco_symbols):
    meta, s, mat = parse_series_matrix(
        os.path.join(RAW, "GSE188309_CAP_WholeBlood", "GSE188309_series_matrix.txt.gz"))
    return fast_collapse_wanted(mat, GPL23159_CSV, set(monaco_symbols))


def load_gse310929():
    tsv = os.path.join(RAW, "GSE310929_Sepsis", "GSE310929_AllSampleExpressionSubmitted.tsv",
                       "GSE310929_AllSampleExpressionSubmitted.tsv")
    m3 = pd.read_csv(tsv, sep="\t", index_col=0)
    m3.index = [str(i).strip('"').upper() for i in m3.index]
    m3.columns = [str(c).strip('"') for c in m3.columns]
    return m3.groupby(level=0).mean().T


def load_gse148871_m2(arms):
    meta, s, mat = parse_series_matrix(
        os.path.join(RAW, "GSE148871_COPD_AE", "GSE148871_series_matrix.txt.gz"))
    man, _, alias = load_manifest()
    probe2g = load_probe_symbol_csv(GPL570_CSV, alias)
    return collapse_probes_to_genes(mat, probe2g, arms["upstream_collapse"] + arms["execution_induction"])


def load_gse148871_m10a(monaco_symbols):
    meta, s, mat = parse_series_matrix(
        os.path.join(RAW, "GSE148871_COPD_AE", "GSE148871_series_matrix.txt.gz"))
    g = fast_collapse_wanted(mat, GPL570_CSV, set(monaco_symbols))
    tis = sample_characteristics(meta, label="tissue")
    blood = [str(t).lower() == "whole blood" for t in tis]
    return g[blood]


def load_gse148871_m11():
    """M11M12 口径：全塌陷（304 样本，含痰液）；返回 (gene_matrix, tissue_list)；
    血样过滤由调用方在评分后执行（冻结表的 z 在 304 全样本上计算）。"""
    meta, s, mat = parse_series_matrix(
        os.path.join(RAW, "GSE148871_COPD_AE", "GSE148871_series_matrix.txt.gz"))
    p2s = read_probe2sym_csv(GPL570_CSV)
    g = collapse_all(mat, p2s)
    tis = sample_characteristics(meta, label="tissue")
    return g, tis


# ------------------------------------------------------------- M11M12 队列
def load_gse215865(cache_dir=None):
    """M11M12 step1 同口径：logCPM(float32) → 符号化 → 受试者×日标签均值（1326 行）。"""
    if cache_dir:
        cp = os.path.join(cache_dir, "GSE215865_merged.npz")
        if os.path.exists(cp):
            z = np.load(cp, allow_pickle=True)
            return pd.DataFrame(z["M"], index=z["samples"].tolist(), columns=z["genes"].tolist())
    cols_file = os.path.join(RAW, "GSE215865_COVID19_WholeBlood_Longitudinal",
                             "GSE215865_rnaseq_logCPM_matrix.csv.gz")
    with gzip.open(cols_file, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n")
    colnames = header.split(",")
    cache_union = {}
    for c in (os.path.join(INTER, "M10A_monaco_ensg2sym.csv"),
              os.path.join(INTER, "P2_gse185263_en2sym.csv")):
        for k, v in pd.read_csv(c).values:
            cache_union[str(k).strip().upper().split(".")[0]] = str(v).strip().upper()
    union_csv = os.path.join(INTER, "M11M12_step1_ensg2sym_union.csv")
    if os.path.exists(union_csv):
        for k, v in pd.read_csv(union_csv).values:
            cache_union[str(k).strip().upper()] = str(v).strip().upper()

    df215 = pd.read_csv(cols_file, compression="gzip", index_col=0, low_memory=False,
                        na_values=["NA"], dtype={c: np.float32 for c in colnames[1:]})
    df215 = df215.astype(np.float32)
    sym_list = [cache_union.get(str(i).split(".")[0].upper(), "") for i in df215.index]
    df215["sym"] = sym_list
    df215 = df215[df215["sym"] != ""]
    df215 = df215.groupby("sym").mean(numeric_only=True)
    df215 = df215.fillna(0.0)
    pat = re.compile(r"^Subj_([0-9a-fA-F]+)T(\d+[A-Za-z]?)_Plate_(\d+)$")
    meta = []
    for c in df215.columns:
        m = pat.match(str(c))
        if m:
            meta.append((str(c), m.group(1), "T" + re.sub(r"[A-Za-z]+$", "", m.group(2))))
        else:
            meta.append((str(c), "", ""))
    mdf = pd.DataFrame(meta, columns=["col", "subject", "day_label"])
    mdf["key"] = mdf["subject"] + "|" + mdf["day_label"]
    mdf = mdf.set_index("col")
    df215_t = df215.T
    merged = df215_t.groupby(mdf["key"]).mean()
    merged.index = merged.index.str.split("|", expand=True)
    merged.index.names = ["subject", "day_label"]
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        np.savez_compressed(cp, M=merged.values.astype(np.float32),
                            samples=np.array([f"{i[0]}|{i[1]}" for i in merged.index], dtype=object),
                            genes=np.array(merged.columns, dtype=object))
    return merged


def load_gse54514():
    f = os.path.join(RAW, "GSE54514_Sepsis_PAXgene_WholeBlood", "GSE54514_non-normalized.txt.gz")
    raw54 = pd.read_csv(f, sep="\t", index_col=0, low_memory=False,
                        dtype={c: np.float32 for c in pd.read_csv(f, sep="\t", nrows=0).columns[1:]})
    expr_cols = [c for c in raw54.columns if not str(c).endswith("Detection Pval")]
    g545_raw = raw54[expr_cols].T
    gsm54 = parse_gsm_metadata(os.path.join(RAW, "GSE54514_Sepsis_PAXgene_WholeBlood",
                                            "GSE54514_gsm_metadata.txt"))
    gsm54["sentrix"] = gsm54["descriptions"].apply(lambda ds: ds[0] if ds else "")
    sentrix2gsm = dict(zip(gsm54["sentrix"], gsm54["gsm"]))
    g545_raw["gsm"] = g545_raw.index.map(lambda c: sentrix2gsm.get(str(c), ""))
    g545_raw = g545_raw[g545_raw["gsm"] != ""].set_index("gsm", drop=True)
    q545 = quantile_normalize(g545_raw)
    g545_log = np.log2(q545 + 1.0)
    p2s = read_gpl_table(GPL6947_TXT)
    return collapse_all(g545_log, p2s)


def load_gse106878():
    f = os.path.join(RAW, "GSE106878_Sepsis_Hydrocortisone_CORTICUS", "GSE106878_non-normalized_data.txt.gz")
    raw878 = pd.read_csv(f, sep="\t", index_col=0, low_memory=False,
                         dtype={c: np.float32 for c in pd.read_csv(f, sep="\t", nrows=0).columns[1:]})
    expr_cols = [c for c in raw878.columns if not str(c).endswith("Detection Pval")]
    g878_raw = raw878[expr_cols].T
    gsm878 = parse_gsm_metadata(os.path.join(RAW, "GSE106878_Sepsis_Hydrocortisone_CORTICUS",
                                             "GSE106878_gsm_metadata.txt"))
    gsm878["sentrix"] = gsm878["descriptions"].apply(
        lambda ds: ds[1] if len(ds) > 1 else (ds[0] if ds else ""))
    sentrix2gsm = dict(zip(gsm878["sentrix"], gsm878["gsm"]))
    g878_raw["gsm"] = g878_raw.index.map(lambda c: sentrix2gsm.get(str(c), ""))
    g878_raw = g878_raw[g878_raw["gsm"] != ""].set_index("gsm", drop=True)
    q878 = quantile_normalize(g878_raw)
    g878_log = np.log2(q878 + 1.0)
    p2s = read_gpl_table(GPL10295_TXT)
    return collapse_all(g878_log, p2s)


# ------------------------------------------------------------------ 演示数据
def load_gse157103():
    """GSE157103 TPM → log2(TPM+1)（126 样本 × 基因符号）。返回 (mat, groups_series)。"""
    f = os.path.join(RAW, "GSE157103_COVID19_WholeBlood_RNAseq", "GSE157103_genes.tpm.tsv.gz")
    tpm = pd.read_csv(f, sep="\t", index_col=0)
    tpm.columns = [str(c).strip('"') for c in tpm.columns]
    tpm.index = [str(i).strip().strip('"').upper() for i in tpm.index]
    tpm = tpm.groupby(level=0).mean()
    mat = np.log2(tpm + 1.0).T
    groups = pd.Series(
        [("COVID" if re.match(r"^C\d+$", s) else ("non-COVID" if re.match(r"^NC\d+$", s) else "?"))
         for s in mat.index], index=mat.index)
    return mat, groups


def load_gse66099():
    """GSE66099 系列矩阵（GPL570）→ 直接符号全塌陷。返回 (mat, groups_series)。"""
    meta, samples, mat = parse_series_matrix(
        os.path.join(RAW, "GSE66099_Pediatric_Sepsis_WholeBlood", "GSE66099_series_matrix.txt.gz"))
    p2s = read_probe2sym_csv(GPL570_CSV)
    g = collapse_all(mat, p2s)
    disease = sample_characteristics(meta, label="disease")
    groups = pd.Series(disease, index=samples)
    return g, groups
