# -*- coding: utf-8 -*-
"""
M2_step1_load_score.py — M2 第 1 步：六队列数据加载 + 逐样本 UCS/EIS/MDI 评分
========================================================================
预注册：M2_pre_registration_20260817.md；定义冻结：mdi_lib.py（mdi_v1.0）
- MDI = z(执行臂均值) - z(上游臂均值)；队列内全样本 z 标准化（不分组，防泄漏）
- 微阵列探针→基因：GEO 官方 GPL annot（GPL23159/GPL570/GPL10558），max-mean 塌陷（每基因 ≤3 个最高均值探针）
- 输出：_intermediate/M2_per_sample_scores.csv（长表，含临床字段）
"""
import gzip
import os
import io
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS")
import mdi_lib as L

ROOT = L.ROOT
RAW = ROOT + r"\00_RAW_DATA"
INTER = ROOT + r"\_intermediate"
os.makedirs(INTER, exist_ok=True)

man, arms, alias = L.load_manifest()
ENSG2SYM = {}
for _, r in man.iterrows():
    e = str(r.get("ensembl_gene_id") or "").strip()
    if e and e.upper().startswith("ENSG"):
        ENSG2SYM[e.upper()] = str(r["hgnc_symbol"]).strip().upper()

log = []
def note(msg):
    log.append(msg)
    print(msg, flush=True)

def score_and_append(df_gene, cohort, meta_df):
    """df_gene: 行=sample, 列=大写基因符号(log2)。meta_df: 同 index 的元数据。返回长表。"""
    sc = L.arm_scores(df_gene, arms)
    out = meta_df.copy()
    for c in sc.columns:
        out[c] = sc[c].values
    out["cohort"] = cohort
    out["platform_coverage_genes"] = len(df_gene.columns)
    return out

frames = []

# ================================================================ GSE185263
note("== GSE185263 (bulk RNA-seq, counts)")
cnt = pd.read_csv(RAW + r"\GSE185263_Lung_ARDS\GSE185263_raw_counts.csv", index_col=0)
cnt.index = [str(i).upper() for i in cnt.index]
cnt.columns = [str(c) for c in cnt.columns]
mapped = {e: s for e, s in ENSG2SYM.items() if e in cnt.index}
note(f"  manifest ENSG mapped in matrix: {len(mapped)}/80")
sub = cnt.loc[[e for e, s in mapped.items()]].copy()
sub.index = [mapped[e] for e in sub.index]
sub = sub.groupby(level=0).mean()          # 重复符号取均值
cpm = sub / sub.sum(axis=0) * 1e6
logc = np.log2(cpm + 1).T                  # samples × genes
groups = pd.read_csv(ROOT + r"\04_AUDIT_GOVERNANCE\GSE185263_groups.csv")
groups.columns = ["sample_id", "group"]
groups["sample_id"] = groups["sample_id"].astype(str)
meta = pd.DataFrame(index=logc.index)
meta["sample_id"] = meta.index
meta = meta.merge(groups, on="sample_id", how="left").set_index(meta.index)
meta["group"] = meta["group"].fillna("unknown")
note(f"  groups: {meta['group'].value_counts().to_dict()}")
frames.append(score_and_append(logc, "GSE185263", meta))

# ================================================================ GSE32707
note("== GSE32707 (Illumina HT-12 v4, gene-level log2 matrix)")
m7 = pd.read_csv(RAW + r"\GSE32707_data\GSE32707_gene_symbol_log2_matrix.csv.gz", index_col=0)
m7.index = [str(i).upper() for i in m7.index]
m7.columns = [str(c) for c in m7.columns]
m7 = m7.groupby(level=0).mean().T
meta7, _, _ = L.parse_series_matrix(RAW + r"\GSE32707_data\GSE32707_series_matrix.txt.gz")
src7 = meta7["!Sample_source_name_ch1"][0]
title7 = meta7["!Sample_title"][0]
gsm7 = meta7["!Sample_geo_accession"][0]
gmap = {
    "untreated": "Control",
    "SIRS Day 0": "SIRS_d0",
    "Sepsis Day 0": "Sepsis_d0", "Sepsis Day 7": "Sepsis_d7",
    "se/ARDS Day 0": "ARDS_d0", "se/ARDS Day 7": "ARDS_d7",
}
meta7df = pd.DataFrame({"gsm": gsm7, "source": src7, "title": title7})
meta7df["group"] = meta7df["source"].map(gmap)
unmapped = meta7df[meta7df["group"].isna()]
if len(unmapped):
    note(f"  WARNING unmapped GSE32707 groups: {unmapped['source'].value_counts().to_dict()}")
meta7df["sample_id"] = meta7df["gsm"]
note(f"  groups: {meta7df['group'].value_counts().to_dict()}")
assert meta7df["group"].notna().all(), "GSE32707 分组映射不完整"
meta7df = meta7df.set_index("gsm")
frames.append(score_and_append(m7, "GSE32707", meta7df[["sample_id", "group"]]))

# ================================================================ GSE212865
note("== GSE212865 (Clariom S probes)")
gpl23159_annot = RAW + r"\GPL23159.annot.gz"
gpl23159_csv = RAW + r"\GPL23159_probe2symbol_clariomsdb.csv"
if not (os.path.exists(gpl23159_annot) or os.path.exists(gpl23159_csv)):
    note("  SKIP: GPL23159 注释（annot.gz 或 Bioconductor CSV）不存在，跳过 GSE212865")
else:
    meta2, s2, m2_ = L.parse_series_matrix(RAW + r"\GSE212865_data\GSE212865_series_matrix.txt.gz")
    probe2g = L.load_probe_map(gpl23159_annot, gpl23159_csv, alias)
    note(f"  GPL23159 probes mapped to 80-gene set: {len(probe2g)}")
    genemat = L.collapse_probes_to_genes(m2_, probe2g, arms["upstream_collapse"] + arms["execution_induction"])
    note(f"  collapsed genes: {list(genemat.columns)}")
    chars = L.sample_characteristics(meta2, label="disease state")
    meta2df = pd.DataFrame({"sample_id": s2, "group": chars}).set_index("sample_id")
    note(f"  groups: {meta2df['group'].value_counts().to_dict()}")
    frames.append(score_and_append(genemat, "GSE212865", meta2df))

# ================================================================ GSE310929
note("== GSE310929 (merged sepsis atlas, gene symbols)")
tsv = RAW + r"\GSE310929_Sepsis\GSE310929_AllSampleExpressionSubmitted.tsv\GSE310929_AllSampleExpressionSubmitted.tsv"
m3 = pd.read_csv(tsv, sep="\t", index_col=0, nrows=None)
m3.index = [str(i).strip('"').upper() for i in m3.index]
m3.columns = [str(c).strip('"') for c in m3.columns]
m3 = m3.groupby(level=0).mean().T
meta3 = pd.read_excel(RAW + r"\GSE310929_Sepsis\GSE310929_AllSampleMetadataSubmitted.xlsx",
                      sheet_name=0, header=1)
meta3.columns = [str(c).strip() for c in meta3.columns]
meta3["Sample ID"] = meta3["Sample ID"].astype(str)
meta3 = meta3.set_index("Sample ID")
meta3 = meta3.loc[m3.index.intersection(meta3.index)]
note(f"  meta matched {len(meta3)}/{len(m3)} samples")
m3 = m3.loc[meta3.index]
keep = meta3[["Dataset", "Disease", "Disease Simplified", "Survival", "DaySurvivalEdited",
              "Time to Event", "TimetoEventEdited", "Age", "Gender", "MolecularSubtype",
              "PreviouslyDefinedSubtypes", "Sepsis/Septic Shock", "Tissue", "Timepoint"]].copy()
keep["sample_id"] = keep.index
keep["group"] = keep["Disease Simplified"]
note(f"  Disease Simplified: {keep['Disease Simplified'].value_counts().head(10).to_dict()}")
note(f"  DaySurvivalEdited: {keep['DaySurvivalEdited'].value_counts(dropna=False).to_dict()}")
frames.append(score_and_append(m3, "GSE310929", keep))

# ================================================================ GSE188309
note("== GSE188309 (CAP, Clariom S probes)")
if not (os.path.exists(gpl23159_annot) or os.path.exists(gpl23159_csv)):
    note("  SKIP: GPL23159 注释不存在，跳过 GSE188309")
else:
    meta8, s8, m8_ = L.parse_series_matrix(RAW + r"\GSE188309_CAP_WholeBlood\GSE188309_series_matrix.txt.gz")
    probe2g = L.load_probe_map(gpl23159_annot, gpl23159_csv, alias)
    genemat8 = L.collapse_probes_to_genes(m8_, probe2g, arms["upstream_collapse"] + arms["execution_induction"])
    note(f"  collapsed genes: {list(genemat8.columns)}")
    chars8 = {}
    for label in ("mortality", "age", "Sex", "batch", "hospital"):
        chars8[label] = L.sample_characteristics(meta8, label=label)
    meta8df = pd.DataFrame({"sample_id": s8, **chars8}).set_index("sample_id")
    meta8df["mortality"] = pd.to_numeric(meta8df["mortality"], errors="coerce")
    meta8df["age"] = pd.to_numeric(meta8df["age"], errors="coerce")
    meta8df["Sex"] = pd.to_numeric(meta8df["Sex"], errors="coerce")
    note(f"  mortality: {meta8df['mortality'].value_counts(dropna=False).to_dict()}")
    note(f"  age: mean={meta8df['age'].mean():.1f}, n={meta8df['age'].notna().sum()}")
    frames.append(score_and_append(genemat8, "GSE188309", meta8df))

# ================================================================ GSE148871
note("== GSE148871 (COPD AE, HG-U133 Plus 2.0 probes)")
gpl570_annot = RAW + r"\GPL570.annot.gz"
gpl570_csv = RAW + r"\GPL570_probe2symbol_hgu133plus2db.csv"
if not (os.path.exists(gpl570_annot) or os.path.exists(gpl570_csv)):
    note("  SKIP: GPL570 注释不存在，跳过 GSE148871")
else:
    meta4, s4, m4_ = L.parse_series_matrix(RAW + r"\GSE148871_COPD_AE\GSE148871_series_matrix.txt.gz")
    probe2g = L.load_probe_map(gpl570_annot, gpl570_csv, alias)
    genemat4 = L.collapse_probes_to_genes(m4_, probe2g, arms["upstream_collapse"] + arms["execution_induction"])
    note(f"  collapsed genes: {list(genemat4.columns)}")
    chars4 = {}
    for label in ("tissue", "treatment", "visit", "Sex", "subject"):
        chars4[label] = L.sample_characteristics(meta4, label=label)
    meta4df = pd.DataFrame({"sample_id": s4, **chars4}).set_index("sample_id")
    note(f"  tissue: {meta4df['tissue'].value_counts().to_dict()}")
    note(f"  treatment: {meta4df['treatment'].value_counts().to_dict()}")
    note(f"  visit: {meta4df['visit'].value_counts().to_dict()}")
    frames.append(score_and_append(genemat4, "GSE148871", meta4df))

# ================================================================ 汇总输出
allf = pd.concat(frames, ignore_index=False)
allf["gene_set_version"] = L.GENE_SET_VERSION
allf["score_version"] = L.SCORE_VERSION
outp = INTER + r"\M2_per_sample_scores.csv"
allf.to_csv(outp, index_label="sample")
note(f"\nTOTAL samples: {len(allf)}; saved to {outp}")
note("per-cohort rows: " + str(allf.groupby('cohort').size().to_dict()))
with open(INTER + r"\M2_step1_log.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(log))
