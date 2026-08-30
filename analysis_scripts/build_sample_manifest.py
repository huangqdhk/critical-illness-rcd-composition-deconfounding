# -*- coding: utf-8 -*-
"""
build_sample_manifest.py — N3 残留4：样本级 manifest 构建（GSM→样本→数据集→条件）
=================================================================================
构建日期：2026-08-16
输出：
  04_AUDIT_GOVERNANCE/SAMPLE_MANIFEST_v1.0.csv   样本级清单（sample 行 + dataset_summary + external_source）
  04_AUDIT_GOVERNANCE/N3_Sample_Manifest_Report.md  构建报告（含全部不变量断言结果）

证据源（全部本地，逐样本可追溯）：
  - GSE145926_family.soft.gz        12 个 scRNA GSM（title=C 样本名；patient group 官方字段）
  - GSE158055_family.soft.gz        284 GSM（title=S-* 样本名；source_name=组织）
  - GSE212865_series_matrix.txt.gz  137 GSM（disease state / time 字段）
  - GSE67530_series_matrix.txt.gz   144 GSM（ards 字段 0/1/NA；顺序与 KS001–KS144 对应，与 S27 零错配）
  - GSE185263_gsm_sample_list.txt   392 GSM（GEO 官方 targ=gsm 列表，2026-08-16 经代理抓取存档于
                                     00_RAW_DATA/GSE185263_Lung_ARDS/；title 与本地样本名零缺失零多余）+ groups.csv×S16e 双源
  - Table_S2_QC_Metrics.csv         63 个 scRNA 样本 per-sample 细胞数（canonical）
  - GSE158055 BALF_filtered/cell_annotation.csv（scFOCAL 34 气道样本细胞计数）
  - GSE171668/GSE165659 目录名 GSM 与 barcodes 行数（GSE171668 物种=GEO 官方 Homo sapiens 尸检肺）
设计：全部条件标签取 canonical 结果表 / GEO 官方字段，不引入任何新推断；
      每个数据集做组计数断言（与 lint/质控报告既定不变量一致），任一断言失败即中止。
"""
import csv
import gzip
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "00_RAW_DATA"
GOV = ROOT / "04_AUDIT_GOVERNANCE"
T = ROOT / "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV"

ROWS = []          # manifest 行累积
FAILS = []


def check(name, cond, detail=""):
    if cond:
        print(f"  [OK] {name}: {detail}")
    else:
        FAILS.append(f"{name}: {detail}")
        print(f"  [FAIL] {name}: {detail}")


def add_sample(dataset_id, gsm, sample_id, condition, tissue, assay, timepoint,
               subcohort, n_cells, role, status, evidence, notes=""):
    ROWS.append({
        "record_type": "sample", "dataset_id": dataset_id, "gsm": gsm, "sample_id": sample_id,
        "condition": condition, "tissue": tissue, "assay": assay, "timepoint": timepoint,
        "subcohort": subcohort, "n_cells": n_cells, "n_units": "", "group_summary": "",
        "role": role, "status": status, "evidence_source": evidence, "notes": notes,
    })


def add_dataset(dataset_id, n_units, group_summary, role, status, evidence, notes=""):
    ROWS.append({
        "record_type": "dataset_summary", "dataset_id": dataset_id, "gsm": "", "sample_id": "",
        "condition": "", "tissue": "", "assay": "", "timepoint": "", "subcohort": "",
        "n_cells": "", "n_units": n_units, "group_summary": group_summary,
        "role": role, "status": status, "evidence_source": evidence, "notes": notes,
    })


def add_external(dataset_id, n_units, group_summary, role, status, evidence, notes=""):
    ROWS.append({
        "record_type": "external_source", "dataset_id": dataset_id, "gsm": "", "sample_id": "",
        "condition": "", "tissue": "", "assay": "", "timepoint": "", "subcohort": "",
        "n_cells": "", "n_units": n_units, "group_summary": group_summary,
        "role": role, "status": status, "evidence_source": evidence, "notes": notes,
    })


# --------------------------------------------------------------------------
# 通用：GEO soft / series matrix 解析
# --------------------------------------------------------------------------
def parse_soft_samples(path):
    """返回 [(gsm, title, source_name, characteristics_str)]"""
    txt = gzip.open(path, "rt", encoding="utf-8", errors="replace").read()
    out = []
    for b in txt.split("^SAMPLE")[1:]:
        def grab(tag):
            m = re.search(tag + r" = \"?(.+?)\"?\s*$", b, re.M)
            return m.group(1).strip() if m else ""
        gsm = re.search(r"!Sample_geo_accession = (GSM\d+)", b).group(1)
        chars = " ;; ".join(re.findall(r"!Sample_characteristics_ch1 = \"?(.+?)\"?\s*$", b, re.M))
        out.append((gsm, grab(r"!Sample_title"), grab(r"!Sample_source_name_ch1"), chars))
    return out


def parse_series_matrix(path):
    """返回 (acc, title, src, per_sample_chars)：per_sample_chars[i] = {tag: value}"""
    gz = gzip.open(path, "rt", encoding="utf-8", errors="replace")
    acc, title, src = [], [], []
    chars_by_sample = defaultdict(list)
    for line in gz:
        if line.startswith("!Sample_geo_accession"):
            acc = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_title"):
            title = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_source_name_ch1"):
            src = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            for i, v in enumerate(line.rstrip("\n").split("\t")[1:]):
                chars_by_sample[i].append(v.strip('"'))
        elif line.startswith("!series_matrix_table_begin"):
            break
    per_sample = []
    for i in range(len(acc)):
        d = {}
        for v in chars_by_sample.get(i, []):
            if ":" in v:
                k2, v2 = v.split(":", 1)
                d[k2.strip()] = v2.strip()
        per_sample.append(d)
    return acc, title, src, per_sample


# ==========================================================================
# 1) GSE145926 — scRNA BALF，12 例（3 HC + 3 mild + 6 severe）
# ==========================================================================
print("[1] GSE145926 (scRNA BALF, 12)")
qc = pd.read_csv(T / "Table_S2_QC_Metrics.csv")
qc145 = qc[qc.dataset == "GSE145926"].set_index("sampleID")
soft = parse_soft_samples(RAW / "GSE145926_scRNA" / "GSE145926_family.soft.gz")
gsm145 = {}
for gsm, title, src, chars in soft:
    if "(scRNA-seq)" not in title:
        continue  # 排除 TCR-seq GSM（本包未用）
    sid = title.split(",")[-1].strip().split(" ")[0]
    grp = re.search(r"patient group: ([^;]+)", chars).group(1).strip()
    gsm145[sid] = (gsm, grp)
cond_map = {"healthy control": "Healthy", "mild": "COVID_mild",
            "severe": "COVID_severe", "severe COVID-19 patient": "COVID_severe"}
g145_cond = Counter()
for sid in sorted(qc145.index, key=lambda s: (s[:1], int(s[1:]))):
    gsm, grp_raw = gsm145[sid]
    cond = cond_map[grp_raw]
    g145_cond[cond] += 1
    add_sample("GSE145926", gsm, sid, cond, "BALF", "scRNA-seq", "", "",
               int(qc145.loc[sid, "n_cells_post_filter"]),
               "主整合 scRNA 队列（BALF）", "canonical",
               "GSE145926_family.soft.gz (patient group) + Table_S2_QC_Metrics",
               "GEO 官方分组（N3 修复基准 3/3/6）")
check("GSE145926 样本数", len(qc145) == 12, f"n={len(qc145)}")
check("GSE145926 分组 3/3/6", dict(g145_cond) == {"Healthy": 3, "COVID_mild": 3, "COVID_severe": 6},
      str(dict(g145_cond)))
check("GSE145926 细胞合计 83,952", int(qc145.n_cells_post_filter.sum()) == 83952,
      f"sum={int(qc145.n_cells_post_filter.sum())}")
add_dataset("GSE145926", "12", "Healthy 3 / COVID_mild 3 / COVID_severe 6（GEO 官方；N3 修复基准）",
            "主整合 scRNA 队列（BALF，Liao et al. 2020）；83,952 细胞", "canonical",
            "family.soft.gz 逐 GSM + S2_QC",
            "GSM 范围 GSM4339769–GSM4475056（本表仅登记 12 个 scRNA-seq GSM；另 9 个 TCR-seq GSM 未用）；N3 曾修正 C148/C149/C152→severe、C144→mild")

# ==========================================================================
# 2) GSE158055 — scRNA：主整合 51 PBMC + scFOCAL 气道子集 34
# ==========================================================================
print("[2] GSE158055 (scRNA, 51 PBMC + 34 airway)")
soft158 = parse_soft_samples(RAW / "GSE158055_scRNA" / "GSE158055_family.soft.gz")
gsm158 = {title: (gsm, src) for gsm, title, src, _ in soft158}
check("GSE158055 GEO GSM 总数 284", len(soft158) == 284, f"n={len(soft158)}")

qc158 = qc[qc.dataset == "GSE158055"].set_index("sampleID")
for sid in sorted(qc158.index):
    gsm, src = gsm158[sid]
    sev = "Healthy" if sid.startswith("S-HC") else ("COVID_mild" if sid.startswith("S-M") else "COVID_severe")
    add_sample("GSE158055", gsm, sid, sev, src + "（GEO source_name）", "scRNA-seq", "", "",
               int(qc158.loc[sid, "n_cells_post_filter"]),
               "主整合 scRNA 队列（PBMC）", "canonical",
               "GSE158055_family.soft.gz (source_name) + Table_S2_QC_Metrics",
               "GEO 组织字段核验 51 例全为 PBMC（N3 更正）")
c158 = Counter(("Healthy" if s.startswith("S-HC") else "COVID_mild" if s.startswith("S-M") else "COVID_severe")
               for s in qc158.index)
check("GSE158055 主整合 51 = 5/21/25", len(qc158) == 51 and dict(c158) == {"Healthy": 5, "COVID_mild": 21, "COVID_severe": 25},
      f"n={len(qc158)}, {dict(c158)}")
check("GSE158055 主整合组织全 PBMC", all(gsm158[s][1] == "PBMC" for s in qc158.index), "51/51 PBMC")
check("GSE158055 细胞合计 54,989", int(qc158.n_cells_post_filter.sum()) == 54989,
      f"sum={int(qc158.n_cells_post_filter.sum())}")

ann = pd.read_csv(RAW / "GSE158055_scRNA" / "BALF_filtered" / "cell_annotation.csv")
cells34 = ann.groupby("sampleID").size()
for sid in sorted(cells34.index):
    gsm, src = gsm158[sid]
    sev = "COVID_mild" if sid.startswith("S-M") else "COVID_severe"
    add_sample("GSE158055", gsm, sid, sev, src + "（GEO source_name）", "scRNA-seq", "", "",
               int(cells34[sid]), "scFOCAL 气道子集输入（S22）", "canonical",
               "GSE158055_family.soft.gz (source_name) + BALF_filtered/cell_annotation.csv",
               "12 BALF + 22 sputum（N3 更正，无健康对照）")
tis34 = Counter(gsm158[s][1] for s in cells34.index)
sev34 = Counter("COVID_mild" if s.startswith("S-M") else "COVID_severe" for s in cells34.index)
check("GSE158055 气道子集 34 = BALF 12 + Sputum 22", len(cells34) == 34 and dict(tis34) == {"BALF": 12, "Sputum": 22},
      f"n={len(cells34)}, {dict(tis34)}")
balf_cells = int(sum(v for s, v in cells34.items() if gsm158[s][1] == "BALF"))
sput_cells = int(sum(v for s, v in cells34.items() if gsm158[s][1] == "Sputum"))
check("scFOCAL 57,225 = 42,723 + 14,502", len(ann) == 57225 and balf_cells == 42723 and sput_cells == 14502,
      f"{balf_cells}+{sput_cells}={len(ann)}")
add_dataset("GSE158055", "51+34（本包使用；图谱共 284 GSM）",
            "主整合 PBMC 51 = 5 HC + 21 mild + 25 severe；气道子集 34 = BALF 12（3m+9s）+ sputum 22（5m+17s）",
            "主整合 scRNA 队列（PBMC）+ scFOCAL 气道子集（S22）", "canonical",
            "family.soft.gz 逐 GSM + S2_QC + cell_annotation.csv",
            "图谱其余 199 个 GSM（PBMC-B/BT/T、PFMC 等）本包未直接使用；N3 曾更正组织标注（旧 BALF→PBMC）")

# ==========================================================================
# 3) GSE171668 — 人 COVID-19 尸检肺 scRNA（Delore et al. Cell 2021 图谱子集），6 GSM 原始 h5 载入 0 细胞
# ==========================================================================
print("[3] GSE171668 (参考)")
h5dir = sorted((RAW / "GSE171668_Lung_ARDS_scRNA").glob("GSM*_raw_feature_bc_matrix.h5"))
g171 = []
for d in h5dir:
    m = re.match(r"(GSM\d+)_(.+)_raw_feature_bc_matrix", d.name)
    g171.append((m.group(1), m.group(2)))
for gsm, name in g171:
    add_sample("GSE171668", gsm, name, "未标注（本地元数据无分组字段）", "人 COVID-19 尸检肺组织（GEO 官方：Homo sapiens）",
               "scRNA-seq", "", "", "", "未参与主整合（S3 组成参考）", "reference",
               "目录名 GSM + GSE171668_bulk_metadata.csv（无疾病列）+ GEO 官方 series 记录",
               "原始 h5 载入 0 细胞（barcode 与 metadata 不匹配，load_summary.json）")
check("GSE171668 原始 GSM 6 个", len(g171) == 6, f"n={len(g171)}")
add_dataset("GSE171668", "6 GSM 原始（发表版 metadata 24 供体 / 106,792 细胞）",
            "无分组字段（COVID-19 尸检供体，非病例-对照设计）",
            "S3 组成参考表使用发表版 lung_metadata（106,043 细胞入表，749 compartment-NA 细胞剔除）；未参与主整合",
            "reference", "load_summary.json + S3 companion 核验 + GEO series 记录",
            "GEO 官方为 Homo sapiens COVID-19 尸检组织图谱（Delore et al. Cell 2021：23 肺+16 肾+15 肝+18 心供体）——"
            "方案旧登记“ARDS 小鼠肺组织/动物模型”有误，2026-08-16 按 GEO 订正为人尸检肺组织")

# ==========================================================================
# 4) GSE185263 — 全血 Bulk 主队列，392（266/82/44）
# ==========================================================================
print("[4] GSE185263 (bulk whole blood, 392)")
g185 = pd.read_csv(GOV / "GSE185263_groups.csv", index_col=0)
s16e = pd.read_csv(T / "Table_S16e_Sample_Condition_Info.csv")
s16e_map = dict(zip(s16e.iloc[:, 0], s16e["Condition"])) if "Condition" in s16e.columns else {}
s16e_map = {str(k): v for k, v in s16e_map.items()}
canon185 = {"Sepsis": "Sepsis", "Sepsis_COVID": "Sepsis_COVID", "Control": "Control"}
mis = [s for s, grp in g185.group.items()
       if s in s16e_map and canon185.get(grp) != {"Healthy_Control": "Control"}.get(s16e_map[s], s16e_map[s])]
check("GSE185263 groups.csv × S16e 双源一致", not mis, f"错配 {mis[:5]}" if mis else "392/392 一致")
# per-sample GSM：GEO 官方样本名列表（2026-08-16 经代理抓取存档；392 个 GSM title 与本地样本名零缺失零多余）
gsm185_file = RAW / "GSE185263_Lung_ARDS" / "GSE185263_gsm_sample_list.txt"
gsm185 = {}
if gsm185_file.exists():
    for b in gsm185_file.read_text(encoding="utf-8", errors="replace").split("^SAMPLE = ")[1:]:
        gsm = b.split("\n", 1)[0].strip()
        m = re.search(r"!Sample_title = (.+)", b)
        if m:
            gsm185[m.group(1).strip()] = gsm
    unmatched = [s for s in g185.index if s not in gsm185]
    check("GSE185263 per-sample GSM 映射（GEO title == 本地样本名）",
          len(gsm185) == 392 and not unmatched,
          f"{len(gsm185)} GSM；本地未匹配 {len(unmatched)}")
else:
    check("GSE185263 per-sample GSM 映射", False,
          "GEO 样本列表存档缺失（00_RAW_DATA/GSE185263_Lung_ARDS/GSE185263_gsm_sample_list.txt）")
sub_census = Counter()
for sid, grp in g185.group.items():
    pref = re.match(r"[a-z]+", sid).group(0)
    sub_census[pref] += 1
    if sid.startswith("sepcv"):
        m = re.search(r"(T\d|W\d)[a-z]?$", sid)
        tp = m.group(1) if m else ""
    else:
        tp = ""
    add_sample("GSE185263", gsm185.get(sid, ""), sid, grp, "全血（R1 核验）", "Bulk RNA-seq", tp, pref, "",
               "Bulk 主队列（DEG S2 / WGCNA S49 / S5 / ML S9 / CIBERSORTx S47 / Scissor S6）", "canonical",
               "GSE185263_groups.csv（R1 正名）+ Table_S16e 双源 + GEO gsm_sample_list（title 精确匹配）", "")
check("GSE185263 n=392 = 266/82/44", len(g185) == 392 and dict(Counter(g185.group)) ==
      {"Sepsis": 266, "Sepsis_COVID": 82, "Control": 44}, str(dict(Counter(g185.group))))
check("GSE185263 前缀亚队列 = R1 证据2", dict(sub_census) ==
      {"sepcol": 67, "sepnet": 104, "sepwes": 84, "sepvh": 11, "sepcv": 82,
       "hccol": 6, "hchl": 9, "hcwes": 24, "hcwimr": 5}, str(dict(sub_census)))
sepcv = [s for s in g185.index if s.startswith("sepcv")]
bases = [re.sub(r"(T\d|W\d)[a-z]?$", "", s) for s in sepcv]
tp_census = Counter((re.search(r"(T\d|W\d)[a-z]?$", s).group(1) if re.search(r"(T\d|W\d)[a-z]?$", s) else "未标注") for s in sepcv)
check("sepcv 无同患者多时点重复", len(bases) == len(set(bases)),
      f"82 个互异编号；时点 {dict(tp_census)}")
add_dataset("GSE185263", "392", "Sepsis 266 / Sepsis_COVID 82 / Control 44（前缀亚队列 sepcol67+sepnet104+sepwes84+sepvh11+sepcv82+hc*44）",
            "全血 Bulk 主队列（Baghela & Hancock 脓毒症队列；R1 正名：82 例 sepcv=脓毒症合并 COVID-19，非 ARDS）",
            "canonical", "groups.csv × S16e 双源 + R1 报告",
            "sepcv 采血时点 T0 21/T1 58（含 sepcv068T1a）/T2 1/W1 2；样本命名层面无同患者重复采样（N3 残留4 逐样本核验）")

# ==========================================================================
# 5) GSE212865 — 外周血微阵列验证队列，137（51/52/34）
# ==========================================================================
print("[5] GSE212865 (microarray blood, 137)")
acc, title, src, per_sample_chars = parse_series_matrix(RAW / "GSE212865_data" / "GSE212865_series_matrix.txt.gz")
g212_cond = {"Control": "Healthy_Control", "Covid19": "COVID19", "Covid19_SDRA": "COVID19_SDRA"}
c212 = Counter()
by_base = defaultdict(list)   # 患者编号（去 _D0/_D7 后缀）-> [(gsm, cond, time)]
for i, gsm in enumerate(acc):
    dis = per_sample_chars[i]["disease state"]
    cond = g212_cond[dis]
    c212[cond] += 1
    base = re.sub(r"_D\d$", "", title[i])
    tp = per_sample_chars[i].get("time", "").replace("NA", "")
    by_base[base].append((gsm, cond, tp))
cross_grp = {b for b, v in by_base.items() if len({x[1] for x in v}) > 1}
n_rep = sum(1 for v in by_base.values() if len(v) > 1)
n_rep_samples = sum(len(v) for v in by_base.values() if len(v) > 1)
for i, gsm in enumerate(acc):
    dis = per_sample_chars[i]["disease state"]
    cond = g212_cond[dis]
    base = re.sub(r"_D\d$", "", title[i])
    tp = per_sample_chars[i].get("time", "").replace("NA", "")
    note = "SDRA=COVID-19 相关 ARDS"
    if len(by_base[base]) > 1:
        grp_flag = ("，**跨疾病组（对照↔疾病）——判定为编号撞号（对照与患者为两套独立去标识化编号空间：原文仅招募 60 例住院"
                    "COVID 患者、51 例对照未在原文任何处出现；对照 GSM 连续成块且无时点标注；2026-08-16 经 GEO series 摘要"
                    "+原文 PMC10216228 核验）**") if base in cross_grp else ""
        note += f"；患者编号 {base} 共 {len(by_base[base])} 个样本（时点 {'/'.join(x[2] or 'NA' for x in by_base[base])}{grp_flag}）"
    add_sample("GSE212865", gsm, f"patient {title[i]}", cond,
               "外周血（全血 RNA）", "微阵列（Bulk）", tp,
               base, "", "外周血微阵列验证队列（S8 / S9 外部验证 / S46b / Bridge）", "canonical",
               "GSE212865_series_matrix.txt.gz (disease state/time)", note)
check("GSE212865 n=137 = 51/52/34", len(acc) == 137 and dict(c212) ==
      {"Healthy_Control": 51, "COVID19": 52, "COVID19_SDRA": 34}, str(dict(c212)))
check("GSE212865 患者编号层面重复已量化（非 137 互异患者）",
      len(by_base) == 96 and n_rep == 37 and n_rep_samples == 78 and len(cross_grp) == 6,
      f"96 个互异编号；{n_rep} 个编号共 {n_rep_samples} 个样本（D0/D7 纵向）；{len(cross_grp)} 个编号跨疾病组")
sdra_bases = {b for b, v in by_base.items() if any(x[1] == "COVID19_SDRA" for x in v)}
check("GSE212865 SDRA 患者号 19 与原文 19 例 ARDS 精确吻合（撞号判定的算术互证）",
      len(sdra_bases) == 19,
      f"SDRA 34 样本来自 {len(sdra_bases)} 个编号（15 对 D0/D7 + 4 单时点）")
add_dataset("GSE212865", "137 样本（编号层面 96 个患者号）",
            "Healthy 51 / COVID-19 52 / COVID-19+SDRA 34（GEO disease state 字段）",
            "外周血微阵列验证队列（唯一真 ARDS 表型层=SDRA 34 例）", "canonical",
            "series_matrix 逐 GSM + GEO series 摘要 + 原文 PMC10216228（撞号核验三源）",
            "采血时点 D0 50 / D7 36 / 对照未标注 51；编号层面非互异——37 个编号贡献 78 个样本（D0/D7 纵向重复采样，为原文设计内行为：入院 24h 内与第 7 天各采一次），"
            "其中 6 个编号（5/18/26/66/80/85）同时在对照组与疾病组出现——**2026-08-16 经 GEO series 摘要+原文（Rombauts et al. Biomedicines 2023, PMC10216228）核验判定为编号撞号**："
            "队列仅招募 60 例住院 COVID 患者（19 例 ARDS；SDRA 组 34 样本恰为 19 编号=15 对 D0/D7+4 单时点，与原文 19 例精确吻合）、"
            "51 例健康对照未在原文任何处描述（仅见于微阵列子系列，GSM 连续成块且无时点标注），对照编号空间（5–295 稀疏）与患者编号空间（1–1320，含 1304–1320 块）独立编制、小编号相撞；"
            "ML 外部验证子集（COVID+SDRA 86 例）中 70 例来自 35 个编号的 D0/D7 重复——外部 AUC 的独立样本基础弱于名义样本量（主稿局限 (3) 已补告诫）")

# ==========================================================================
# 6) GSE67530 — 全血 450K 甲基化，144（30 HC + 75 ICU + 39 ARDS）
# ==========================================================================
print("[6] GSE67530 (450K methylation, 144)")
acc6, title6, src6, flat6 = parse_series_matrix(RAW / "GSE67530_series_matrix.txt.gz")
s27 = pd.read_csv(T / "Table_S27_Epigenetic_Clocks.csv")
s27_map = dict(zip(s27["Sample_ID"], s27["Group"]))
mismatch = 0
c67530 = Counter()
for i, gsm in enumerate(acc6):
    ks = "KS%03d" % (i + 1)
    a = flat6[i]["ards"]
    geo_grp = "Healthy_Control" if a in ("NA", "") else ("ARDS" if a == "1" else "ICU_Control")
    if s27_map.get(ks) != geo_grp:
        mismatch += 1
    c67530[geo_grp] += 1
    note = ""
    if ks == "KS123":
        note = "S27b QC：Excluded_data_error（表内保留，时钟分析按脚本口径含全 144 例）"
    elif ks == "KS050":
        note = "S27b QC：Notable_biological_extreme（保留）"
    add_sample("GSE67530", gsm, ks, geo_grp, "whole blood（GEO tissue 字段）", "DNA 甲基化 450K", "", "", "",
               "表观遗传时钟队列（S27/S27b/S27c）", "canonical",
               "GSE67530_series_matrix.txt.gz (ards 字段) + Table_S27 零错配核验", note)
check("GSE67530 n=144 = 30/75/39（GEO ards 字段）", len(acc6) == 144 and dict(c67530) ==
      {"Healthy_Control": 30, "ICU_Control": 75, "ARDS": 39}, str(dict(c67530)))
check("GSE67530 顺序 GSM 映射 × S27 零错配", mismatch == 0, f"mismatch={mismatch}/144")
add_dataset("GSE67530", "144", "Healthy 30 / ICU 非 ARDS 75 / ARDS 39（GEO ards 字段；S27 同口径）",
            "表观遗传时钟队列（全血 450K）", "canonical",
            "series_matrix 逐 GSM（ards 字段）× S27 零错配",
            "主稿 Methods 原写 24/72/48 系误记，2026-08-16 按本 manifest 更正为 30/75/39；KS123 数据错误样本、KS050 极端值样本（S27b 注）")

# ==========================================================================
# 6b) GSE32707 — 全血 HumanHT-12 V4.0 微阵列，144（Davenport et al. 2016；D1b 2026-08-16）
#     Control 34 / SIRS d0 21 / Sepsis d0 30 / Sepsis d7 28 / se-ARDS d0 18 / se-ARDS d7 13
# ==========================================================================
print("[6b] GSE32707 (whole blood HumanHT-12, 144, D1b)")
acc32, title32, src32, flat32 = parse_series_matrix(RAW / "GSE32707_data" / "GSE32707_series_matrix.txt.gz")
check("GSE32707 n=144", len(acc32) == 144, f"n={len(acc32)}")
subj32 = [d.get("subject id", "") for d in flat32]
check("GSE32707 subject id 字段全覆盖", all(s for s in subj32), "144/144")


def _g32707(s):
    s = s.strip()
    if s.lower().startswith("untreated"):
        return "Control", "NA"
    up = s.upper()
    day = "d0" if "DAY 0" in up else ("d7" if "DAY 7" in up else "NA")
    if "SE/ARDS" in up or "SE-ARDS" in up:
        return "ARDS", day
    if "SIRS" in up:
        return "SIRS", day
    return "Sepsis", day


c32707 = Counter()
dup_subj = Counter()
for gsm, s, sb in zip(acc32, src32, subj32):
    cond32, tp32 = _g32707(s)
    c32707[(cond32, tp32)] += 1
    if cond32 == "Control":
        dup_subj[sb] += 1
    note = ""
    if cond32 == "Control" and dup_subj[sb] > 1:
        note = "对照技术重复（subject %s 共 %d 管，主分析保留、S54 受试者塌陷敏感性）" % (sb, dup_subj[sb])
    add_sample("GSE32707", gsm, gsm, cond32, "whole blood（GEO source_name）",
               "表达微阵列（GPL10558 HumanHT-12 V4.0）", tp32, sb, "",
               "D1b 全血 ARDS 外部验证队列（S50–S54）", "canonical",
               "GSE32707_series_matrix.txt.gz (source_name + subject id) + D1b 元数据核验报告",
               note)
check("GSE32707 组×时点 = 34/21/30/28/18/13",
      dict(c32707) == {("Control", "NA"): 34, ("SIRS", "d0"): 21, ("Sepsis", "d0"): 30,
                       ("Sepsis", "d7"): 28, ("ARDS", "d0"): 18, ("ARDS", "d7"): 13},
      str(dict(c32707)))
pair_chk = defaultdict(set)
for s, sb in zip(src32, subj32):
    c32, d32 = _g32707(s)
    if c32 in ("Sepsis", "ARDS"):
        pair_chk[sb].add(d32)
check("GSE32707 病例组无双时点受试者（非配对横断面）", all(len(v) == 1 for v in pair_chk.values()),
      f"{sum(1 for v in pair_chk.values() if len(v) > 1)} 个双时点")
check("GSE32707 对照 34 样本=22 受试者（subject 128115 ×12 技术重复）",
      len(dup_subj) == 22 and max(dup_subj.values()) == 12 and dup_subj.get("128115") == 12,
      f"n_subj={len(dup_subj)}, max={dup_subj.most_common(1)}")
add_dataset("GSE32707", "144",
            "Control 34 / SIRS d0 21 / Sepsis d0 30 / Sepsis d7 28 / se-ARDS d0 18 / se-ARDS d7 13（GEO source_name 逐样本）",
            "D1b 全血真 ARDS 表型外部验证队列（S50–S54；Davenport et al. Sci Transl Med 2016）", "canonical",
            "series_matrix 逐 GSM + D1b_GSE32707_Metadata_Check_Report.md + Table_S50_D1b_GSE32707_Groups.csv",
            "非配对横断面两时点（无双时点受试者；主分析 d0 子集）；对照 34 样本=22 受试者（subject 128115 占 12 管技术重复，主分析保留+受试者塌陷敏感性 S54）；"
            "80 基因可测 74/80——6 个 MT 基因（MT-ATP6/ATP8/CO1/CYB/ND1/ND2）为 GPL10558 官方注释层真实缺席（S52 登记）；GSDME 以旧名 DFNA5 检出（manifest aliases 2026-08-16 登记）")

# ==========================================================================
# 7) GSE165659 — scATAC 肺基线，4 健康供体
# ==========================================================================
print("[7] GSE165659 (sci-ATAC-seq lung, 4)")
atac_files = sorted((RAW / "GSE165659_Lung_scATAC").glob("GSM*_lung_*_barcodes.txt.gz"))
atac_n = 0
for f in atac_files:
    m = re.match(r"(GSM\d+)_lung_(SM-[A-Za-z0-9]+)_barcodes", f.name)
    n = sum(1 for _ in gzip.open(f, "rt"))
    atac_n += n
    add_sample("GSE165659", m.group(1), m.group(2), "健康供体（主稿 Methods 口径：图谱基线）",
               "肺组织", "sci-ATAC-seq", "", "", n,
               "scATAC 染色质可及性基线（S26）", "canonical",
               "目录名 GSM + barcodes 行数", "raw barcodes；SnapATAC2 QC 后全队列 21,790 细胞（主稿 Methods）")
check("GSE165659 4 样本", len(atac_files) == 4, f"n={len(atac_files)}；raw barcodes 合计 {atac_n}（≥QC 后 21,790）")
add_dataset("GSE165659", "4", "健康供体 4（无病例分组；主稿定位为基线非 ARDS 重塑）",
            "scATAC-seq 染色质可及性基线（S26）", "canonical", "目录名 + 主稿 Methods",
            "raw barcodes 23,274 → SnapATAC2 QC 后 21,790 细胞 / 562,444 bins")

# ==========================================================================
# 8) external_source 行（非样本型数据源；N10 inventory + N16 孤儿脚本输入登记）
# ==========================================================================
print("[8] external sources")
add_external("ieu-b-69", "486,484", "脓毒症 GWAS 汇总统计", "cis-MR 主分析结局（S15a）", "canonical",
             "主稿 Methods（IEU OpenGWAS）", "")
add_external("ieu-b-5086", "486,484", "脓毒症 28 天死亡 GWAS", "cis-MR 主分析结局（S15a）", "canonical",
             "主稿 Methods（IEU OpenGWAS）", "")
add_external("FinnGen_R10_J10_ARDS", "357 病例 / 406,536 对照（EUR）", "真 ARDS 结局 GWAS",
             "S19 疾病关联检索 + cis-MR 敏感性结局 + coloc（S15b/S15c）", "canonical",
             "00_RAW_DATA/GWAS/finngen_R10_J10_ARDS + 方案 §13.2", "")
add_external("eQTLGen", "31,470", "全血 cis-eQTL 汇总统计", "MR/coloc 暴露工具变量来源（S15a/S15b/S15d）", "canonical",
             "主稿 Methods", "cis-only 过滤后 34 工具/15 基因")
add_external("gnomAD_v2.1.1", "141,456 个体", "LoF/错义约束统计", "S40 纯化选择约束（Bridge 维度）", "canonical",
             "主稿 Methods", "")
add_external("AlphaMissense", "全库 182,567 错义（55 基因面板）", "错义致病性评分", "S18 致病变异负荷（Bridge 维度）", "canonical",
             "主稿 Methods", "")
add_external("CMap_L1000", "19 种通路相关药物", "转录连接性扰动签名", "S23 药物重定位评分", "canonical", "主稿 Methods", "")
add_external("1000G_PLINK_参考面板", "AFR 661 + AMR 等", "LD 参考基因型", "MR LD clumping（r²<0.001）", "reference",
             "00_RAW_DATA/AFR.fam 等（HG* 个体 ID 核验）", "")
add_external("HLCA_肺细胞图谱", "参考注释", "肺细胞类型参考映射", "单细胞注释（S26a）", "reference",
             "annotate_cell_types_lung_atlas.py", "")
add_external("cCRE_EpiAgent 参考", "cCRE.bed + document_frequency.npy", "染色质元件参考",
             "scATAC 分析参考（run_scATAC_SCENIC_joint_27.1.py）", "reference", "00_RAW_DATA/EpiAgent_reference/", "")
add_external("GSE171524", "54 GSM（已下载）", "IPF 肺组织空间转录组",
             "空间验证（B1）——2026-08-15 暂缓 ⏭️", "downloaded_unused",
             "00_RAW_DATA/GSE171524_RAW（C51ctr–C57ctr 等样本文件）", "未产出任何结果表；可视化类任务暂缓决策")
add_external("GSE200042", "已下载（tar+RAW）", "小鼠 LPS ALI 空间转录组",
             "空间验证——随 B1 暂缓 ⏭️", "downloaded_unused", "00_RAW_DATA/GSE200042_Mouse_LPS_ALI_Spatial/", "")
add_external("MosMedData_CT", "CT-0/CT-23 等 zip", "COVID-19 胸部 CT 影像",
             "影像组学×基因组（radiomics_imaging_genomics_24.1.py / radiomics_v2_real_data.py 输入）——N16 保留在册",
             "exploratory_asset", "00_RAW_DATA/MosMedData_CT + 孤儿脚本头注",
             "N16 处置（2026-08-16）：7 个孤儿脚本（radiomics×2、run_mofa×2、GEARS、celloracle×2）保留为探索性资产并登记入册，不删除；是否纳入留档位 C 启动时定夺")

# ==========================================================================
# 汇总断言 + 写出
# ==========================================================================
sc_total = int(qc145.n_cells_post_filter.sum()) + int(qc158.n_cells_post_filter.sum())
check("scRNA 合计 138,941", sc_total == 138941, f"{sc_total}")

if FAILS:
    print(f"\n!! {len(FAILS)} 项断言失败，中止写出：")
    for f in FAILS:
        print("   -", f)
    sys.exit(1)

df = pd.DataFrame(ROWS)
cols = ["record_type", "dataset_id", "gsm", "sample_id", "condition", "tissue", "assay",
        "timepoint", "subcohort", "n_cells", "n_units", "group_summary", "role", "status",
        "evidence_source", "notes"]
df = df[cols]
out = GOV / "SAMPLE_MANIFEST_v1.0.csv"
df.to_csv(out, index=False, encoding="utf-8-sig", lineterminator="\r\n")
print(f"\n写出 {out.name}: {len(df)} 行 = sample {sum(df.record_type == 'sample')} + "
      f"dataset_summary {sum(df.record_type == 'dataset_summary')} + "
      f"external_source {sum(df.record_type == 'external_source')}")
print("样本行按数据集：", dict(Counter(df[df.record_type == 'sample'].dataset_id)))
json.dump({"n_rows": len(df),
           "by_type": dict(Counter(df.record_type)),
           "by_dataset": dict(Counter(df[df.record_type == 'sample'].dataset_id))},
          open(ROOT / "03_LOGS" / "sample_manifest_build.json", "w", encoding="utf-8"), ensure_ascii=False)
print("全部断言通过（ALL CHECKS PASSED）")
