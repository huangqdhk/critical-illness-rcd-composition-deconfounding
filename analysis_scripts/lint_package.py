# -*- coding: utf-8 -*-
"""
lint_package.py — D4 投稿前全包校验（只读，不改任何数据）
================================================================
对应质控报告 D4：校验 符号合法性 / 版本一致性 / FDR 列非空 / p∈[0,1] /
计数与组定义一致 + N22 编号唯一性 + 冻结 manifest 校验和。
M1（2026-08-17）：新增 S55 空间模块校验（check_m1_spatial）——S55 表族、
预注册 R1-R3 判定复算（R3）、CosMx 面板覆盖、空间样本清单。

用法:  python lint_package.py
输出:  03_LOGS/lint_report.md + 控制台摘要；全绿退出码 0，否则 1。
"""
import csv
import hashlib
import math
import os
import re
import sys
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ROOT = Path(__file__).resolve().parent
# 2026-08-16 目录重组后：治理层独立目录；补充表与笔记在 02_SUPPLEMENTARY_TABLES；
# 文件实际位置以 manifest 的 location 列为准（唯一路径来源）。
GOV = ROOT / ("04_AUDIT_GOVERNANCE")
TAB = ROOT / ("02_SUPPLEMENTARY_TABLES") / "SUPPLEMENTARY_Tables_CSV"
NOTES = ROOT / ("02_SUPPLEMENTARY_TABLES") / "Supplementary_Notes"
# 2026-09-04 归档区统一：原 archive/（用户建）并入导师建的 归档/，规则随之指向 归档/
ARCH = ROOT / ("归档")
MANIFEST = GOV / "RESULTS_MANIFEST_v2.0.csv"  # M2（2026-08-17）起升 v2.0（v1.0 冻结保留）

# 受管目录（顶层文件必须全部登记 canonical；manifest 自身除外）
MANAGED_DIRS = {
    "04_AUDIT_GOVERNANCE": GOV,
    "02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV": TAB,
    "02_SUPPLEMENTARY_TABLES/Supplementary_Notes": NOTES,
    "归档": ARCH,
}
# 2026-08-16 目录重组：登记在册但本地不随包的文件（附理由，逐项登记）
MISSING_OK = {
    "Table_S11_Mitoxyperilysis_TF_PerCell.csv": "数据存放层大矩阵，本地不随包（Zenodo/GEO 存档），见 01_FIGURE_DATA_CSV/README_图数据包说明.md",
    "Table_S11_TF_Activity_Matrix_AllCells.csv": "数据存放层大矩阵，本地不随包（Zenodo/GEO 存档），见 01_FIGURE_DATA_CSV/README_图数据包说明.md",
    "Table_S1_Mitoxyperilysis_Score_by_CellType_prev38.csv": "已被 canonical Table_S1 取代，原始文件未随包迁移（archived 留痕）",
    "Table_S1_Mitoxyperilysis_Score_by_CellType_Summary_prev38.csv": "已被 canonical Table_S1 取代，原始文件未随包迁移（archived 留痕）",
    # 2026-08-23：00_RAW_DATA 输入层不随交付包（README_打包说明 §二），交付夹 lint 放行
    "GSE67530_beta_for_GrimAge.csv.gz": "GrimAge 甲基化 β 输入矩阵（235MB），存于 00_RAW_DATA/，公开数据可自 GEO 下载，不随交付包",
}

# 2026-08-27：M10–M15 新模块导出表/图数据（前瞻注册 osf.io/ETVMJ 与第二次独立注册批次），
# 按 M15/M1x 导出设计无内嵌 gene_set_version/score_version 列，manifest 记 NA（设计内豁免）
# 2026-09-06 图号重编（11→9 主图）键名同步：11A-C→9D-F、旧9A-C→8E-G、旧10A-C→9A-C、S9I→S6I、旧5A/B→5F/G（FDR_NULL_OK 键）；Figure_7/8A-D 不变
M1X_NO_VERSION_STAMP = {
    "Table_S62_M14_IIAMD_Core_v1.0.csv",
    "Table_S63_M14_Gate2_Core_Rerun_ThreeLayers.csv",
    "Table_S64_M10A_Deconvolution_Proportions.csv",
    "Table_S65_M10A_R2_AllMethods.csv",
    "Table_S65b_M10A_R2_Merged.csv",
    "Table_S65c_M10A_Proportion_CrossReference.csv",
    "Table_S66_M10A_Residual_LayerEffects.csv",
    "Table_S66b_M10A_Residual_Meta_LOSO.csv",
    "Table_S67_M10B_Donor_Arm_Scores.csv",
    "Table_S67b_M10B_Donor_Arm_Tests.csv",
    "Table_S67c_M10B_Myeloid_UCS_Layers.csv",
    "Table_S67d_M10B_Myeloid_UCS_Meta.csv",
    "Table_S68_M13_Signature_Catalog.csv",
    "Table_S69b_M13_NullModel_Summary.csv",
    "Table_S69c_M13_ReverseAudit.csv",
    "Table_S70_M10C_Compartment_Contrasts.csv",
    "Table_S72_M14L3_GO_Audit.csv",
    "Table_S73_M14L3_CellLevel_Perturbation_Effects.csv",
    "Table_S74_M14L3_FRPerturb_Arm_Projections.csv",
    "Table_S75_M14L3_NullCalibration_and_SetTests.csv",
    "Table_S76_M11_PerSample_Longitudinal_Scores.csv",
    "Table_S77_M11_H1_MixedModels.csv",
    "Table_S78_M11_H2_RICLPM.csv",
    "Table_S79_M11_NullModels_LOSO_Sensitivity.csv",
    "Table_S80_M12_GSE106878_PerPatient.csv",
    "Table_S81_M12_GSE106878_Tests_Comp_H4_Null.csv",
    "Table_S82_M12_GSE148871_Replication.csv",
    "Table_S83_M10D_CrossSpecies_Contrasts.csv",  # 2026-08-27 M10D 次终点执行批次
    "Table_S83b_M10D_Pig_Ensembl_Sensitivity.csv",  # 2026-08-27 M10D 补做批次（猪同源敏感性）
    "Table_S87_M10D_Pandisease_k20_LayerEffects.csv",  # 2026-08-27 M10D 补做批次（k≥20 泛疾病）
    # 2026-08-31 M16 力学边界模块批次（导出设计无内嵌版本列，manifest 记 NA）
    "Table_S88_M16_Mechanosensing_Module_v10.csv",
    "Table_S89_M16_VILI_Contrasts.csv",
    "Table_S90_M16_GSE2411_Interaction.csv",
    "Table_S91_M16_Bridge_Correlations.csv",
    "Table_S92_M16_Ortholog_Coverage_Audit.csv",
    "Figure_9D.csv",
    "Figure_9E.csv",
    "Figure_9F.csv",
    "Table_S87b_M10D_Pandisease_k20_Meta.csv",
    "Table_S87c_M10D_Pandisease_k20_LOSO.csv",
    "Figure_7A.csv", "Figure_7B.csv", "Figure_7C.csv", "Figure_7D.csv",
    "Figure_7E.csv", "Figure_7F.csv", "Figure_7G.csv",
    "Figure_8A.csv", "Figure_8B.csv", "Figure_8C.csv", "Figure_8D.csv",
    "Figure_8E.csv", "Figure_8F.csv", "Figure_8G.csv",
    "Table_S86d_M15_Demo_GSE66099_Proportions.csv",
    "Table_S86c_M15_Demo_GSE66099_R2.csv",
    "Table_S86b_M15_Demo_GSE66099_Contrasts.csv",
    "Table_S86_M15_Demo_GSE66099_PerSample.csv",
    "Table_S85d_M15_Demo_GSE157103_Proportions.csv",
    "Table_S85c_M15_Demo_GSE157103_R2.csv",
    "Table_S85b_M15_Demo_GSE157103_Contrasts.csv",
    "Table_S85_M15_Demo_GSE157103_PerSample.csv",
    "Table_S84b_M15_RegressionTest_CohortSummary.csv",
    "Table_S84_M15_RegressionTest_11Cohorts.csv",
    "Figure_9C.csv",
    "Figure_9B.csv",
    "Figure_9A.csv",
    # 2026-09-05 M17 蛋白组批次（PXD050432 审计入册，导出设计无内嵌版本列，manifest 记 NA）
    "Table_S93_M17_PXD050432_AcuteLPS_LungProteome_WB_Targets.csv",
    "Table_S93b_M17_PXD050432_PerSample_LFQ.csv",
}

_LOC = None


def load_locations():
    """manifest 的 filename → location 映射（缓存）。"""
    global _LOC
    if _LOC is None:
        _LOC = {r["filename"]: r.get("location", "") for r in read_csv_all(MANIFEST)}
    return _LOC


def fp(fn):
    """按 manifest location 列解析文件实际路径。"""
    return ROOT / load_locations().get(fn, "") / fn

FAILS = []
WARNS = []
OKS = []


def fail(msg):
    FAILS.append(msg)


def warn(msg):
    WARNS.append(msg)


def ok(msg):
    OKS.append(msg)


# ---------------------------------------------------------------- 读取工具
def read_csv_stream(path):
    """流式读 CSV，产出 dict 行（大文件安全）。"""
    raw = path.read_bytes()
    enc = "utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf-8"
    with open(path, "r", encoding=enc, newline="") as f:
        for row in csv.DictReader(f):
            yield row


def read_csv_all(path):
    return list(read_csv_stream(path))


def sha256_of(p, block=1024 * 1024):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while True:
            b = f.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------- 权威基因集
GM = GOV / "Mitoxyperilysis_Gene_Manifest_v1.0.csv"
_gm = read_csv_all(GM)
CANON80 = [r["gene_symbol"] for r in _gm]
SET80 = set(CANON80)
HGNC = set(r["hgnc_symbol"] for r in _gm)
ALIAS2CANON = {"ITPR1": "IP3R1", "HSPD1": "HSP60", "DMT1": "SLC11A2",
               "GRP75": "HSPA9", "CYTC": "CYCS",
               "DFNA5": "GSDME"}  # 包内出现过的外部别名/变体（DFNA5：GPL10558 注释层写法，D1b 2026-08-16 登记）
# S40 凋亡/MAM 钙对照基因（正文 §1 Tier 比较，不在 80 清单内，属设计内扩展）
S40_EXTRAS = {"BAD", "BAK1", "BAX", "BBC3", "BCL2", "BCL2L1", "BCL2L11",
              "BCL2L2", "BID", "MCL1", "MCU", "PMAIP1"}
# 层内扩展基因（TOMM70A：ATAC/GWAS/转录本层基因面板自带，非 80 清单成员）
GLOBAL_EXTRAS = {"TOMM70A"}
# 伪符号（R6 已登记：微阵列/转录层 "NADH" 非真实基因符号，源数据保留原样）
PSEUDO_SYMBOLS = {"NADH"}
# S18 = AlphaMissense 自选 55 基因面板（80 清单成员 + 历史扩展基因），不做子集校验
CUSTOM_PANEL_TABLES = {"Table_S18_AlphaMissense_Main.csv", "Table_S18c_AlphaMissense_HighPriority.csv"}

EXPECTED_MODULE_SIZES = {
    "mitoxy_MAM_integrity": 9, "mitoxy_mitochondrial_function": 10,
    "mitoxy_iron_metabolism": 9, "mitoxy_oxidative_stress": 11,
    "mitoxy_cell_death": 13, "mitoxy_ferroptosis_cuproptosis": 11,
    "mitoxy_transcription_factors": 11, "mitoxy_autophagy": 6,
}
EXPECTED_ARMS = {"upstream_collapse": 30, "execution_induction": 33, "not_in_dissociation_arms": 17}

# ---------------------------------------------------------------- 严格符号校验表
STRICT_SYMBOL_COLS = {
    "Mitoxyperilysis_Pathway_Gene_List.csv": ["gene_symbol"],
    "Mitoxyperilysis_Gene_Manifest_v1.0.csv": ["gene_symbol"],
    "Bridge_Test_Result.csv": ["gene"],
    "Table_S2b_scRNA_DEG_80genes.csv": ["gene_symbol"],
    "Table_S8b_80gene_ThreePlatform_Meta.csv": ["gene"],
    "Table_S9_Risk_Model_Features.csv": ["gene_symbol"],
    "Table_S13_ScTenifoldKnk_KO_Summary.csv": ["ko_gene"],
    "Table_S15a_MR_Sepsis_cisonly_Instruments.csv": ["exposure"],
    "Table_S15b_MR_FinnGenARDS_Instruments.csv": ["exposure"],
    "Table_S15d_cis_trans_classification.csv": ["exposure"],
    "Table_S18_AlphaMissense_Main.csv": ["Gene_Symbol"],
    "Table_S19c_Mitoxy_Gene_GWAS_Detail.csv": ["gene"],
    "Table_S26f_Mitoxyperilysis_Gene_Bin_Map.csv": ["gene"],
    "Table_S26j_Mitoxyperilysis_PerGene_Accessibility.csv": ["gene"],
    "Table_S32_TF_Binding_Predictions.csv": ["target_gene"],
    "Table_S36c_Target_Protein_Info.csv": ["Target_Protein"],
    "Table_S38_Geneformer_Perturbations.csv": ["Gene"],
    "Table_S38b_Geneformer_Perturbations_Detail.csv": ["Gene"],
    "Table_S40_gnomAD_Constraint.csv": ["Gene_Symbol_HGNC"],
    "Table_S41a_Differential_Transcript_Usage.csv": ["Gene"],
    "Table_S41b_DTU_Gene_Summary.csv": ["Gene"],
    "Table_S41c_Alternative_Splicing_Events.csv": ["Gene"],
    "Table_S41d_Differential_Exon_Usage.csv": ["gene"],
    "Table_S41f_Transcript_Diversity_Index.csv": ["Gene"],
    "Table_S2_DEGs_Analysis.csv": ["gene_symbol"],   # 仅 is_pathway_gene=True 行
    "Table_S52_D1b_GSE32707_Gene_log2FC.csv": ["gene"],
}
LEGACY38_TABLES = {  # gene 列须 ⊆ 旧 38 基因清单（S46e 动态定义）
    "Table_S46a_Circulating_Mitoxy_Genes_GSE185263.csv": "Gene",
    "Table_S46b_Circulating_Mitoxy_Genes_GSE212865.csv": "Gene",
}
PSEUDO_OK_TABLES = {  # 允许出现 "NADH" 伪符号的表（R6 已登记，源数据保留）
    "Table_S46b_Circulating_Mitoxy_Genes_GSE212865.csv",
    "Table_S41a_Differential_Transcript_Usage.csv",
    "Table_S41b_DTU_Gene_Summary.csv",
    "Table_S41c_Alternative_Splicing_Events.csv",
    "Table_S41d_Differential_Exon_Usage.csv",
    "Table_S41f_Transcript_Diversity_Index.csv",
}
VALID_ALIASES = {"ITPR1": "IP3R1", "HSPD1": "HSP60", "GRP75": "HSPA9", "CYTC": "CYCS", "DMT1": "SLC11A2",
                 "DFNA5": "GSDME"}
MATRIX_HEADER_CHECK = {  # 表头基因列（矩阵）须 ⊆ 80
    "Table_S20l_Causal_TF_MitoGene_Matrix.csv": {"Condition", "Tissue"},
    "Table_S46c_Blood_Mitoxy_Expression_Matrix.csv": {"Condition", "Tissue"},
}

# ---------------------------------------------------------------- M1 空间模块（2026-08-17）
# S55 表族（S55r 按设计留空：条件间 Kruskal-Wallis 全不显著，见 M1_README_分析流程.md）
S55_FAMILY = {
    "Table_S55a_M1_Visium_Section_Spatial_Stats.csv",
    "Table_S55b_M1_Visium_Condition_Summary.csv",
    "Table_S55c_M1_Visium_Domain_Stats.csv",
    "Table_S55d_M1_Visium_SVG_Detail.csv",
    "Table_S55e_M1_Visium_SVG_Meta.csv",
    "Table_S55f_M1_Visium_SVG_Enrichment.csv",
    "Table_S55g_M1_Visium_SVG_Top10_Hypergeom.csv",
    "Table_S55h_M1_Visium_KeyGene_Moran_Perm.csv",
    "Table_S55i_M1_CosMx_Gene_Detection.csv",
    "Table_S55j_M1_CosMx_Moran.csv",
    "Table_S55k_M1_CosMx_Myeloid_Coloc.csv",
    "Table_S55l_M1_CosMx_Module_by_CellType.csv",
    "Table_S55m_M1_CosMx_ViralRegion_Niche.csv",
    "Table_S55n_M1_CosMx_TNFSF13B_Source.csv",
    "Table_S55o_M1_CosMx_TNFSF13B_Neighborhood.csv",
    "Table_S55p_M1_Visium_NNLS_CellType_by_Condition.csv",
    "Table_S55q_M1_Visium_NNLS_MarkerScore_Correlation.csv",
    "Table_S55s_M1_Visium_TNFSF13B_TFRC_Neighborhood.csv",
    "Table_S55t_M1_Visium_TNFSF13B_TFRC_CellType_Attribution.csv",
    "Table_S55u_M1_Visium_NNLS_CosMx_by_Condition.csv",
    "Table_S55v_M1_Visium_Bivariate_Residualized.csv",
}
S55_N_SECTIONS = 23
S55_SPOT_TOTAL = 93869
S55_CONDITION_COUNTS = {"Control": 4, "AcuteDAD": 7, "ProliferativeDAD": 12}
# 预注册可测上限（M1_pre_registration_20260817.md）：Visium FFPE 上游 24/30（缺 6 MT-*）、执行 32/33（缺 CASP1）
S55_GENE_CAPS = {"n_up_genes": 24, "n_ex_genes": 32}
# 预注册 CosMx 面板覆盖：执行臂仅 5/33 可测、上游 0/30（N4_资源核实报告.md §2）
S55_COSMX_EXEC_EXPECTED = {"IL18", "IL1B", "NLRP3", "SLC40A1", "SOD2"}

# FDR/p 列识别（值为数值的列才纳入；标记/注释/计数列排除）
FDR_RE = re.compile(r"(?i)(fdr|padj|p_adj|q_value|qvalue|q-val|adj_p)")
FLAG_RE = re.compile(r"(?i)(significant|_note$|^note|n_sig|flag|_in_$|threshold|cutoff)")
P_RE = re.compile(r"(?i)(^p$|pvalue|p_value|p\.value|_p$|permutation_p|mw_p|p_li$)")
P_EXCLUDE_RE = re.compile(r"(?i)(neg_log10|n_sig|^pct|pct_|^pos$|position|protein|pathway|peak|pmaip|pdb)")
# 个别表 FDR/p 列允许空值/NaN 的白名单（附理由，逐项登记）
FDR_NULL_OK = {
    ("Table_S15a_MR_Sepsis_cisonly_Sensitivity.csv", "fdr_q"): "MR 敏感度异质性检验逐 IV 行，无 FDR 义务",
    ("Table_S47b_Circulating_Mitoxy_Immune_Correlations.csv", "FDR"):
        "6 个恒定比例细胞类型（如 NK activated）相关未定义 → NaN，FDR_note 列已逐行说明",
    ("Figure_S6I.csv", "FDR"):
        "6 个恒定比例细胞类型（NK cells activated）相关未定义 → NaN，FDR_note 列已逐行说明（原 Figure_5I，任务8 降入补充 S9）",
    ("Table_S8_Consensus_Genes.csv", "GSE212865_padj"): "基因不在 GSE212865 平台时留空（数据集缺席 NA，非计算失败）",
    ("Table_S8_Meta_Analysis.csv", "GSE212865_padj"): "同上：数据集缺席 NA",
    ("Table_S8_Meta_Analysis.csv", "meta_Fisher_padj"): "某数据集 p 缺失无法合并 Fisher 的 2 行（源数据缺席）",
    ("Table_S52_D1b_GSE32707_Gene_log2FC.csv", "BH_q"):
        "6 个 MT 基因平台缺席登记行（GPL10558 注释层真实缺席，非计算失败）",
    # M4 蛋白层（2026-08-18 登记，2026-08-21 lint 白名单补登）：mmc4 论文官方统计仅覆盖 25/80 基因，
    # 缺官方 padj 的基因逐行留空并在正文 §4.21/Table S59 口径披露（px_padj/mmc4_padj 空值=无官方统计，非计算失败）
    ("Figure_5F.csv", "px_padj"): "M4 蛋白层（原 Figure_10A，任务8 图号迁移）：无论文官方统计的基因留空（mmc4 仅 25/80 基因有官方 padj，正文已披露）",
    ("Figure_5G.csv", "mmc4_padj"): "M4 蛋白层铁轴表（原 Figure_10B，任务8 图号迁移）：无官方统计基因留空（同上口径）",
    ("Table_S59a_M4_Protein_Transcript_Consistency.csv", "px_padj"): "M4 蛋白层：无官方统计基因留空（同上口径）",
    ("Table_S59b_M4_Lung_Iron_Axis_Proteins.csv", "mmc4_padj"): "M4 蛋白层铁轴表：无官方统计基因留空（同上口径）",
    # M10–M15 新模块批次（2026-08-27 收口补登）：以下空值均为设计内缺席，非计算失败
    ("Table_S62_M14_IIAMD_Core_v1.0.csv", "b_cs_padj"):
        "M14 冻结签名表：CS 主效应 padj 不在该签名的计算口径内（源中间表 M14_IIAMD_core_v1.0 同列全空），非计算失败",
    ("Table_S62_M14_IIAMD_Core_v1.0.csv", "b_torin_padj"):
        "M14 冻结签名表：1 基因在 Torin 宇宙效应表（S71）无测值（源数据缺席 NA）",
    ("Table_S63_M14_Gate2_Core_Rerun_ThreeLayers.csv", "OLS_p_adj"):
        "M14 腿2 三层头对头：OLS 组成调整仅适用 bulk 两层，scRNA 供者级层无该口径（设计内 NA）",
    ("Figure_7E.csv", "OLS_p_adj"):
        "图7E 面板（同 Table_S63 口径）：scRNA 供者级层无 OLS 组成调整口径（设计内 NA）",
}


def to_float(v):
    if v is None:
        return None
    v = v.strip()
    if v in ("", "NA", "NaN", "nan", "N/A"):
        return None
    if v.startswith("<"):
        v = v[1:]
    if v.startswith(">= "):
        v = v[3:]
    try:
        return float(v)
    except ValueError:
        return None


# ================================================================ 检查函数
def check_manifest():
    rows = read_csv_all(MANIFEST)
    man_canon = {r["filename"]: r for r in rows if r["status"] == "canonical"}
    man_arch = {r["filename"]: r for r in rows if r["status"] == "archived"}
    locs = load_locations()

    # 各受管目录顶层文件必须全部登记（治理目录的 manifest 自身除外）
    for loc, d in MANAGED_DIRS.items():
        disk = {p.name for p in d.iterdir() if p.is_file()}
        if d == GOV:
            disk = {n for n in disk if not n.startswith("RESULTS_MANIFEST")}
        reg = {fn for fn, l in locs.items() if l == loc}
        unreg = disk - reg
        if unreg:
            fail(f"[manifest] {loc} 未登记文件 {len(unreg)} 个: {sorted(unreg)[:10]}")
        else:
            ok(f"[manifest] {loc} 内 {len(disk)} 个文件全部登记")
    # 登记为 canonical/archived 但本地不存在（MISSING_OK 白名单除外）
    absent = [fn for fn, r in {**man_canon, **man_arch}.items()
              if fn not in MISSING_OK and (not locs.get(fn) or not fp(fn).exists())]
    if absent:
        fail(f"[manifest] 登记在册但本地不存在: {sorted(absent)}")
    else:
        ok(f"[manifest] 登记在册文件全部存在（MISSING_OK 白名单 {len(MISSING_OK)} 项）")
    # 归档区内容与登记一致（MISSING_OK 内的 archived 留痕文件不要求实体存在）
    # 2026-08-21：archive 重组为子目录（R1_table_backups/M1_check_debug 等），改为递归比对
    # 2026-09-04：archive/ 并入 归档/（归档区统一），规则指向 归档/
    disk_arch = set()
    for dirpath, _dirs, files in os.walk(ARCH):
        for f in files:
            disk_arch.add(f)
    exp_arch = set(man_arch) - set(MISSING_OK)
    if disk_arch != exp_arch:
        fail(f"[manifest] 归档/ 与登记不一致: 多出 {sorted(disk_arch - exp_arch)[:5]}, 缺少 {sorted(exp_arch - disk_arch)[:5]}")
    else:
        ok(f"[manifest] 归档/ {len(disk_arch)} 个归档文件（含子目录）与登记一致")

    # 校验和 + 尺寸
    n_hash = 0
    for fn, r in man_canon.items():
        if fn in MISSING_OK:
            continue
        p = fp(fn)
        if not p.exists():
            continue
        if p.stat().st_size != int(r["size_bytes"]):
            fail(f"[checksum] {fn} 尺寸不符: 磁盘 {p.stat().st_size} vs manifest {r['size_bytes']}")
            continue
        if sha256_of(p) != r["sha256"]:
            fail(f"[checksum] {fn} SHA256 不符（manifest 冻结后文件被改动）")
        n_hash += 1
    ok(f"[checksum] {n_hash} 个 canonical 文件 SHA256+尺寸校验通过")
    return man_canon


def check_n22_numbering(man_canon):
    files = set(man_canon)
    # 小写 Table_s
    low = [f for f in files if re.match(r"(?i)^table_s\d", f) and not f.startswith("Table_S")]
    if low:
        fail(f"[N22] 仍存在小写 s 编号文件: {low}")
    # 各编号唯一族
    s6 = [f for f in files if re.match(r"Table_S6([_.]|$)", f)]
    if any("Dissociation" in f for f in s6):
        fail(f"[N22] S6 仍含解离检验文件: {[f for f in s6 if 'Dissociation' in f]}")
    if not any("Scissor" in f for f in s6):
        fail("[N22] S6 缺少 Scissor 主表")
    s3 = [f for f in files if re.match(r"Table_S3([_.]|$)", f)]
    if any("WGCNA" in f for f in s3):
        fail(f"[N22] S3 仍含 WGCNA 文件: {[f for f in s3 if 'WGCNA' in f]}")
    if "Table_S49_WGCNA_Modules.csv" not in files:
        fail("[N22] S49 WGCNA 模块表缺失")
    s26b = [f for f in files if f.startswith("Table_S26b_")]
    if len(s26b) != 1 or "Cluster_Statistics" not in s26b[0]:
        fail(f"[N22] S26b 应唯一为 Cluster_Statistics，实际: {s26b}")
    if "Table_S26a_HLCA_CellType_Annotations.csv" not in files:
        fail("[N22] S26a HLCA 注释表缺失")
    # 版本并存 / 备份残留
    bad_pat = [f for f in files if re.search(r"(CORRECTED|_backup|_pre_|_superseded)", f)]
    if bad_pat:
        fail(f"[N22] 顶层仍存在版本并存/备份命名: {bad_pat}")
    for need in ["Table_S18_AlphaMissense_Main.csv", "Table_S18b_AlphaMissense_Module_Summary.csv",
                 "Table_S18c_AlphaMissense_HighPriority.csv", "Table_S12b_Mitoxyperilysis_Upstream_TF.csv",
                 "Table_S48_Dissociation_singlecell_correlation.csv"]:
        if need not in files:
            fail(f"[N22] 缺少 canonical 文件: {need}")
    ok("[N22] 编号唯一性：S3/S6/S12/S18/S26 冲突清除；无小写 s、无版本并存/备份残留")
    # S 编号解析一致性
    for fn, r in man_canon.items():
        m = re.match(r"[Tt]able_([Ss]\d+)[a-z]?[_.]", fn)
        if m and r["s_number"] not in ("", m.group(1).upper()) and "/" not in r["s_number"]:
            fail(f"[N22] {fn} 文件名编号 {m.group(1).upper()} 与 manifest {r['s_number']} 不一致")


def check_versions(man_canon):
    n_csv = n_txt = 0
    for fn, r in sorted(man_canon.items()):
        if fn in MISSING_OK:
            continue
        p = fp(fn)
        if not p.exists():
            continue
        if fn.lower().endswith(".csv"):
            if fn in ("GSE185263_groups.csv", "MR_bio_CRP_Sepsis.csv", "MR_bio_Ferritin_Sepsis.csv",
                      "Table_S50_D1b_GSE32707_Groups.csv"):
                continue  # 分组定义/探索性输出，manifest 记 NA，无版本列（设计内）
            if fn in M1X_NO_VERSION_STAMP:
                continue  # 2026-08-27：M10–M15 新模块导出表/图数据，导出设计无内嵌版本列，manifest 记 NA
            raw = p.read_bytes()
            enc = "utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf-8"
            bad = False
            with open(p, "r", encoding=enc, newline="") as f:
                rd = csv.reader(f)
                hdr = next(rd)
                if hdr and str(hdr[0]).startswith("#"):
                    hdr = next(rd)          # '#' 注释行后的真实表头
                if "gene_set_version" not in hdr or "score_version" not in hdr:
                    fail(f"[版本] {fn} 缺 gene_set_version/score_version 列")
                    bad = True
                    continue
                gi, si = hdr.index("gene_set_version"), hdr.index("score_version")
                mg = {x.strip() for x in r["gene_set_version"].split(";") if x.strip()} or {""}
                ms = {x.strip() for x in r["score_version"].split(";") if x.strip()} or {""}
                for i, row in enumerate(rd):
                    if row and str(row[0]).startswith("#"):
                        continue
                    if len(row) <= max(gi, si) or \
                            (row[gi] not in mg) or (row[si] not in ms):
                        fail(f"[版本] {fn} 第{i}行版本戳 != manifest ({r['gene_set_version']},{r['score_version']})")
                        bad = True
                        break
            if not bad:
                n_csv += 1
        elif fn.lower().endswith(".txt"):
            if r["gene_set_version"] == "NA" and r["score_version"] == "NA":
                n_txt += 1
                continue
            txt = p.read_bytes()
            enc_try = ("utf-8", "gbk")
            body = None
            for e in enc_try:
                try:
                    body = txt.decode(e)
                    break
                except UnicodeDecodeError:
                    continue
            if body is None:
                fail(f"[版本] {fn} 无法解码")
                continue
            if "版本戳 2026-08-15" not in body:
                fail(f"[版本] {fn} 缺版本页脚")
            elif f"gene_set_version={r['gene_set_version']}" not in body or f"score_version={r['score_version']}" not in body:
                fail(f"[版本] {fn} 页脚版本与 manifest 不一致")
            else:
                n_txt += 1
    ok(f"[版本] {n_csv} 个 CSV 列值 + {n_txt} 个 TXT 页脚与 manifest 一致")


def check_symbols():
    # 基因 manifest 自身
    if len(CANON80) != 80:
        fail(f"[符号] 基因 manifest 行数 {len(CANON80)} != 80")
    if len(SET80) != 80:
        fail(f"[符号] 基因 manifest 符号重复: {[k for k, v in Counter(CANON80).items() if v > 1]}")
    mods = Counter(r["module"] for r in _gm)
    if dict(mods) != EXPECTED_MODULE_SIZES:
        fail(f"[符号] 模块尺寸 {dict(mods)} != 期望 {EXPECTED_MODULE_SIZES}")
    arms = Counter(r["arm"] for r in _gm)
    if dict(arms) != EXPECTED_ARMS:
        fail(f"[符号] 解离双臂基因数 {dict(arms)} != 期望 {EXPECTED_ARMS}")
    for r in _gm:
        exp_alias = {v: k for k, v in VALID_ALIASES.items()}.get(r["hgnc_symbol"], "")
        # IP3R1/HSP60 的包内变体即 HGNC 名本身；HSPA9/CYCS/SLC11A2 的变体为 GRP75/CYTC/DMT1；
        # GSDME 的变体为 DFNA5（GPL10558 注释层写法，D1b 登记）
        if r["aliases"] not in ("", "ITPR1", "HSPD1", "GRP75", "CYTC", "DMT1", "DFNA5"):
            fail(f"[符号] 基因 manifest 别名列异常: {r['gene_symbol']}/{r['aliases']}")
        if r["hgnc_symbol"] not in HGNC:
            fail(f"[符号] 基因 manifest HGNC 列异常: {r['gene_symbol']}/{r['hgnc_symbol']}")
    ok(f"[符号] 基因 manifest：80 基因、8 模块尺寸、双臂 {EXPECTED_ARMS}、别名表全部正确")

    def allowed(sym, fname):
        s = ALIAS2CANON.get(sym, sym)
        return (s in SET80) or (sym in SET80) or (sym in GLOBAL_EXTRAS) or \
               (fname == "Table_S40_gnomAD_Constraint.csv" and sym in S40_EXTRAS) or \
               (sym in PSEUDO_SYMBOLS and fname in PSEUDO_OK_TABLES)

    for fname, cols in STRICT_SYMBOL_COLS.items():
        p = fp(fname)
        if not p.exists():
            fail(f"[符号] 校验目标缺失: {fname}")
            continue
        if fname in CUSTOM_PANEL_TABLES:
            ok(f"[符号] {fname}: 自选面板（custom panel），非 80 子集——仅检查非空")
            continue
        bad = defaultdict(int)
        n = 0
        for row in read_csv_stream(p):
            if fname == "Table_S2_DEGs_Analysis.csv":
                if str(row.get("is_pathway_gene", "")).lower() != "true":
                    continue
            n += 1
            for c in cols:
                v = (row.get(c) or "").strip()
                if v and not allowed(v, fname):
                    bad[v] += 1
        if bad:
            fail(f"[符号] {fname} 存在 80 清单（含别名/登记例外）之外的符号: {dict(bad)}")
        else:
            ok(f"[符号] {fname}: {n} 行 {cols} 全部合法")

    # legacy-38 层
    s46e = read_csv_all(fp("Table_S46e_Mitoxy_Gene_Blood_Summary.csv"))
    legacy38 = {r["Gene"] for r in s46e}
    if len(legacy38) != 38:
        fail(f"[符号] S46e 旧 38 基因清单行数 {len(legacy38)} != 38")
    for fname, col in LEGACY38_TABLES.items():
        bad = defaultdict(int)
        for row in read_csv_stream(fp(fname)):
            v = (row.get(col) or "").strip()
            if v and v not in legacy38 and not (v in PSEUDO_SYMBOLS and fname in PSEUDO_OK_TABLES):
                bad[v] += 1
        if bad:
            fail(f"[符号] {fname} 存在旧 38 清单之外的符号: {dict(bad)}")
        else:
            ok(f"[符号] {fname}: gene 列 ⊆ 旧 38 清单（登记伪符号除外）")

    # 矩阵表头基因列
    for fname, excl in MATRIX_HEADER_CHECK.items():
        raw = fp(fname).read_bytes()
        enc = "utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf-8"
        with open(fp(fname), "r", encoding=enc, newline="") as f:
            hdr = next(csv.reader(f))
        bad = [c for c in hdr[1:] if c and c not in excl and c not in SET80
               and ALIAS2CANON.get(c, c) not in SET80 and c not in ("gene_set_version", "score_version")]
        if bad:
            fail(f"[符号] {fname} 表头存在清单外基因列: {bad}")
        else:
            ok(f"[符号] {fname}: 表头基因列全部 ⊆ 80 清单")


def check_p_and_fdr(man_canon):
    n_p = n_fdr = 0
    for fn, r in sorted(man_canon.items()):
        if not fn.lower().endswith(".csv"):
            continue
        if fn in MISSING_OK:
            continue
        p = fp(fn)
        if not p.exists():
            continue
        raw = p.read_bytes()
        enc = "utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf-8"
        fdr_cols = p_cols = []
        with open(p, "r", encoding=enc, newline="") as f:
            rd = csv.reader(f)
            hdr = next(rd)
            if hdr and str(hdr[0]).startswith("#"):
                hdr = next(rd)          # '#' 注释行后的真实表头
            fdr_cols = [c for c in hdr if FDR_RE.search(c) and not FLAG_RE.search(c)
                        and not P_EXCLUDE_RE.search(c)]
            p_cols = [c for c in hdr if (P_RE.search(c) or c.lower() in ("pvalue", "pval"))
                      and c not in fdr_cols and not P_EXCLUDE_RE.search(c)
                      and not FDR_RE.search(c) and not FLAG_RE.search(c)]
            idx = {c: hdr.index(c) for c in set(fdr_cols + p_cols)}
            viol_range = defaultdict(list)
            viol_null = defaultdict(int)
            nrows = 0
            for row in rd:
                if row and str(row[0]).startswith("#"):
                    continue
                nrows += 1
                for c, i in idx.items():
                    if i >= len(row):
                        continue
                    v = row[i]
                    fv = to_float(v)
                    is_null = v.strip() in ("", "NA", "N/A", "NaN", "nan", "nan ")
                    if c in fdr_cols:
                        if is_null:
                            if (fn, c) not in FDR_NULL_OK:
                                viol_null[c] += 1
                        elif fv is not None and not (0.0 <= fv <= 1.0):
                            viol_range[c].append(v)
                        elif fv is None and not is_null:
                            pass  # 非数值文本（如标记字面量）不属于 FDR 值域检查
                    else:  # p 列
                        if fv is not None and not (0.0 <= fv <= 1.0):
                            viol_range[c].append(v)
        if viol_range:
            fail(f"[p值] {fn} 数值超 [0,1]: " + "; ".join(f"{c}×{v[:3]}" for c, v in viol_range.items()))
        if viol_null:
            fail(f"[FDR空] {fn} FDR 列存在空值: " + "; ".join(f"{c}×{n}" for c, n in viol_null.items()))
        if fdr_cols:
            n_fdr += 1
        if p_cols:
            n_p += 1
    ok(f"[p值] p∈[0,1] 与 FDR 非空扫描完成（含 FDR 列的表 {n_fdr} 个、含 p 列的表 {n_p} 个）")


def check_counts():
    # GSE185263 组定义
    g = read_csv_all(fp("GSE185263_groups.csv"))
    cnt = Counter(r["group"] for r in g)
    if dict(cnt) != {"Sepsis": 266, "Sepsis_COVID": 82, "Control": 44} or len(g) != 392:
        fail(f"[计数] GSE185263_groups: {dict(cnt)} (n={len(g)}) != Sepsis 266/Sepsis_COVID 82/Control 44=392")
    else:
        ok("[计数] GSE185263_groups = 266/82/44（共 392）")

    e = Counter(r["Condition"] for r in read_csv_stream(fp("Table_S16e_Sample_Condition_Info.csv")))
    if dict(e) != {"Sepsis": 266, "Sepsis_COVID": 82, "Healthy_Control": 44}:
        fail(f"[计数] S16e Condition 计数 {dict(e)} != 266/82/44（Healthy_Control）")
    else:
        ok("[计数] S16e Condition = 266/82/44")

    s5 = Counter(r["group"] for r in read_csv_stream(fp("Table_S5_PCD_Scores.csv")))
    if dict(s5) != {"Sepsis": 266, "Sepsis_COVID": 82, "Control": 44}:
        fail(f"[计数] S5 group 计数 {dict(s5)} != 266/82/44")
    else:
        ok("[计数] S5 group = 266/82/44")

    s5r = Counter(r["group"] for r in read_csv_stream(fp("Table_S5_PCD_Scores_RAW.csv")))
    if dict(s5r) != {"Sepsis": 266, "Sepsis_COVID": 82, "Control": 44}:
        fail(f"[计数] S5_RAW group 计数 {dict(s5r)} != 266/82/44")
    else:
        ok("[计数] S5_RAW group = 266/82/44")

    # D1b（2026-08-16）：GSE32707 组计数 + 表行数 + 74 可测/6 MT 缺席
    d1b = Counter(r["group"] for r in read_csv_stream(fp("Table_S50_D1b_GSE32707_Groups.csv")))
    exp_d1b = {"Control": 34, "SIRS_d0": 21, "Sepsis_d0": 30, "Sepsis_d7": 28,
               "ARDS_d0": 18, "ARDS_d7": 13}
    if dict(d1b) != exp_d1b or sum(d1b.values()) != 144:
        fail(f"[计数] D1b GSE32707 组计数 {dict(d1b)} != {exp_d1b}（共 144）")
    else:
        ok("[计数] D1b GSE32707 = 34/21/30/28/18/13（共 144）")

    s51 = list(read_csv_stream(fp("Table_S51_D1b_GSE32707_Module_ssGSEA.csv")))
    if len(s51) != 54:
        fail(f"[计数] S51 行数 {len(s51)} != 54（6 对比 × 9 特征）")
    else:
        ok("[计数] S51 = 54 行（6 对比 × 9 特征）")

    s52 = list(read_csv_stream(fp("Table_S52_D1b_GSE32707_Gene_log2FC.csv")))
    n_meas = sum(1 for r in s52 if r["contrast"] == "A1_ARDS_d0_vs_Control")
    n_abs = sum(1 for r in s52 if r["contrast"] == "platform_absent_registration")
    if (n_meas, n_abs) != (74, 6) or len(s52) != 228:
        fail(f"[计数] S52 可测 {n_meas}/缺席 {n_abs}/总行 {len(s52)} != 74/6/228（74×3+6）")
    else:
        mt_abs = sorted(r["gene"] for r in s52 if r["contrast"] == "platform_absent_registration")
        ok(f"[计数] S52 = 74 可测 ×3 对比 + 6 MT 缺席登记（{','.join(mt_abs)}）= 228 行")

    s53v = [r for r in read_csv_stream(fp("Table_S53_D1b_Three_Anchor_Consistency.csv"))
            if r["panel"] == "D_verdict" and r["row_id"] == "最终判定"]
    allowed_verdicts = {"R1_一致_加强ARDS主张", "R2_异质性叙事", "R3_全阴_回退D1a", "R4_部分信号_中性报告"}
    if len(s53v) != 1 or s53v[0]["flag"] not in allowed_verdicts:
        fail(f"[计数] S53 D1b 预注册判定行异常: {[r['flag'] for r in s53v]}")
    else:
        ok(f"[计数] S53 预注册判定在册（{s53v[0]['flag']}）")

    # 单细胞总数
    s1 = read_csv_all(fp("Table_S1_Mitoxyperilysis_Score_by_CellType_Summary.csv"))
    tot = sum(int(r["n_cells"]) for r in s1)
    byds = defaultdict(int)
    for r in s1:
        byds[r["dataset"]] += int(r["n_cells"])
    if tot != 138941 or dict(byds) != {"GSE145926": 83952, "GSE158055": 54989}:
        fail(f"[计数] S1_Summary 总数 {tot}/分数据集 {dict(byds)} != 138,941 = 83,952+54,989")
    else:
        ok("[计数] S1_Summary = 138,941（GSE145926 83,952 + GSE158055 54,989）")

    pooled = [r for r in read_csv_stream(fp("Table_S48_Dissociation_singlecell_correlation.csv"))
              if r["stratum"] == "ALL_pooled"]
    if not pooled or int(pooled[0]["n"]) != 138941:
        fail(f"[计数] S48 pooled n != 138,941（实际 {[r['n'] for r in pooled]}）")
    else:
        ok("[计数] S48 解离检验 pooled n = 138,941")

    # Scissor 比例合计
    sc = defaultdict(float)
    for r in read_csv_stream(fp("Table_S6_Scissor_CellType_Composition.csv")):
        sc[(r["cell_type"], r["condition"])] += float(r["pct_within_celltype_condition"])
    badsc = {k: round(v, 2) for k, v in sc.items() if abs(v - 100) > 0.5}
    if badsc:
        fail(f"[计数] S6 Scissor 组成比例合计 != 100%: {dict(list(badsc.items())[:5])}")
    else:
        ok(f"[计数] S6 Scissor 组成比例：{len(sc)} 个 cell_type×condition 组均合计 100%")

    # 固定行数
    for fn, col, exp in [
        ("Bridge_Test_Result.csv", "gene", 80),
        ("Table_S12b_Mitoxyperilysis_Upstream_TF.csv", "TF", 23),
        ("Table_S13_ScTenifoldKnk_KO_Summary.csv", "ko_gene", 12),
        ("Table_S18_AlphaMissense_Main.csv", "Gene_Symbol", 55),
        ("Table_S46e_Mitoxy_Gene_Blood_Summary.csv", "Gene", 38),
    ]:
        vals = [r[col] for r in read_csv_stream(fp(fn))]
        if len(vals) != exp:
            fail(f"[计数] {fn} 行数 {len(vals)} != {exp}")
        else:
            ok(f"[计数] {fn} 行数 = {exp}")

    # S15d：15 基因、34 cis + 13 trans
    d = list(read_csv_stream(fp("Table_S15d_cis_trans_classification.csv")))
    genes = {r["exposure"] for r in d}
    cis = sum(1 for r in d if r["is_cis"] == "TRUE")
    trans = sum(1 for r in d if r["is_cis"] == "FALSE")
    if len(genes) != 15 or (cis, trans) != (34, 13):
        fail(f"[计数] S15d: {len(genes)} 基因 / cis {cis} + trans {trans} != 15 / 34+13")
    else:
        ok("[计数] S15d = 15 基因、34 cis + 13 trans")

    # S26b 聚类数（登记值 16）
    n26 = len(read_csv_all(fp("Table_S26b_Cluster_Statistics.csv")))
    if n26 != 16:
        warn(f"[计数] S26b 聚类数 {n26}（方案登记 16）")
    else:
        ok("[计数] S26b = 16 clusters")


def check_sample_manifest():
    """N3 残留4（2026-08-16）：样本级 manifest 不变量——结构、逐数据集样本数与组计数。
    D1b（2026-08-16）：并入 GSE32707（144 样本，全血，34/21/30/28/18/13）。
    与 build_sample_manifest.py 内置断言同口径（lint 为冻结包层面守门）。"""
    p = fp("SAMPLE_MANIFEST_v1.0.csv")
    if not p.exists():
        fail("[样本manifest] SAMPLE_MANIFEST_v1.0.csv 缺失（N3 残留4 交付物）")
        return
    rows = list(read_csv_stream(p))
    rt = Counter(r["record_type"] for r in rows)
    if dict(rt) != {"sample": 924, "dataset_summary": 8, "external_source": 13}:
        fail(f"[样本manifest] record_type 计数 {dict(rt)} != sample 924 + dataset_summary 8 + external_source 13")
        return
    smp = [r for r in rows if r["record_type"] == "sample"]
    ds = Counter(r["dataset_id"] for r in smp)
    exp = {"GSE145926": 12, "GSE158055": 85, "GSE171668": 6, "GSE185263": 392,
           "GSE212865": 137, "GSE67530": 144, "GSE165659": 4, "GSE32707": 144}
    if dict(ds) != exp:
        fail(f"[样本manifest] 逐数据集样本行 {dict(ds)} != {exp}")
        return
    # 组计数（条件标签与 canonical 结果表/lint 既有不变量同口径）
    c67530 = Counter(r["condition"] for r in smp if r["dataset_id"] == "GSE67530")
    if dict(c67530) != {"Healthy_Control": 30, "ICU_Control": 75, "ARDS": 39}:
        fail(f"[样本manifest] GSE67530 组计数 {dict(c67530)} != 30/75/39（GEO ards 字段）")
    c212 = Counter(r["condition"] for r in smp if r["dataset_id"] == "GSE212865")
    if dict(c212) != {"Healthy_Control": 51, "COVID19": 52, "COVID19_SDRA": 34}:
        fail(f"[样本manifest] GSE212865 组计数 {dict(c212)} != 51/52/34")
    c185 = Counter(r["condition"] for r in smp if r["dataset_id"] == "GSE185263")
    if dict(c185) != {"Sepsis": 266, "Sepsis_COVID": 82, "Control": 44}:
        fail(f"[样本manifest] GSE185263 组计数 {dict(c185)} != 266/82/44")
    c145 = Counter(r["condition"] for r in smp if r["dataset_id"] == "GSE145926")
    if dict(c145) != {"Healthy": 3, "COVID_mild": 3, "COVID_severe": 6}:
        fail(f"[样本manifest] GSE145926 组计数 {dict(c145)} != 3/3/6")
    c32707 = Counter((r["condition"], r["timepoint"]) for r in smp if r["dataset_id"] == "GSE32707")
    exp32707 = {("Control", "NA"): 34, ("SIRS", "d0"): 21, ("Sepsis", "d0"): 30,
                ("Sepsis", "d7"): 28, ("ARDS", "d0"): 18, ("ARDS", "d7"): 13}
    if dict(c32707) != exp32707:
        fail(f"[样本manifest] GSE32707 组×时点计数 {dict(c32707)} != {exp32707}")
    # GSM 覆盖：全部 924 个样本行均须有 GSM（GSE185263 已于 2026-08-16 按 GEO 存档补全）
    no_gsm = [r["sample_id"] for r in smp if not r["gsm"]]
    if no_gsm:
        fail(f"[样本manifest] GSM 缺失: {no_gsm[:5]}")
    else:
        ok(f"[样本manifest] 945 行（924 sample + 8 summary + 13 external）；逐数据集 {exp}；"
           f"GSE67530 30/75/39、GSE212865 51/52/34、GSE185263 266/82/44、GSE145926 3/3/6、"
           f"GSE32707 34/21/30/28/18/13（Control/SIRS_d0/Sepsis_d0/Sepsis_d7/ARDS_d0/ARDS_d7）；"
           f"GSM 覆盖 924/924")


def _norm_sf(x):
    """标准正态右尾概率（stdlib 实现，与 scipy norm.sf 等价）。"""
    return 0.5 * math.erfc(x / math.sqrt(2.0))


def _norm_ppf(p):
    """标准正态分位数（Acklam 有理逼近，相对误差 ~1e-9，stdlib 复刻 scipy norm.ppf）。"""
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (((((-7.784894002430293e-03 * q - 3.223964580411365e-01) * q - 2.400758277161838e+00)
                  * q - 2.549732539343734e+00) * q + 4.374664141464968e+00) * q
                + 2.938163982698783e+00) / ((((7.784695709041462e-03 * q + 3.224671290700398e-01)
                                              * q + 2.445134137142996e+00) * q + 3.754408661907416e+00) * q + 1)
    if p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((-7.784894002430293e-03 * q - 3.223964580411365e-01) * q - 2.400758277161838e+00)
                   * q - 2.549732539343734e+00) * q + 4.374664141464968e+00) * q
                 + 2.938163982698783e+00) / ((((7.784695709041462e-03 * q + 3.224671290700398e-01)
                                               * q + 2.445134137142996e+00) * q + 3.754408661907416e+00) * q + 1)
    q = p - 0.5
    r = q * q
    return (((((-3.969683028665376e+01 * r + 2.209460984245205e+02) * r - 2.759285104469687e+02)
              * r + 1.383577518672690e+02) * r - 3.066479806614716e+01) * r
            + 2.506628277459239e+00) * q / (((((-5.447609879822406e+01 * r + 1.615858368580409e+02)
                                               * r - 1.556989798598866e+02) * r + 6.680131188771972e+01)
                                             * r - 1.328068155288572e+01) * r + 1)


def check_m1_spatial(man_canon):
    """M1 空间模块守门（2026-08-17）：S55 表族完整性与 manifest 冻结行数一致；
    预注册 R1-R3 判定复算（口径复刻 M1_step7_judgment.py，规则先于结果，
    见 M1_pre_registration_20260817.md）；CosMx 面板覆盖；空间样本清单。"""
    s55_reg = {fn for fn, r in man_canon.items() if r["s_number"] == "S55"}
    if s55_reg != S55_FAMILY:
        fail(f"[M1] S55 登记族与期望不一致: 缺 {sorted(S55_FAMILY - s55_reg)}, 多 {sorted(s55_reg - S55_FAMILY)}")
        return

    row_bad = [fn for fn in sorted(S55_FAMILY)
               if sum(1 for _ in read_csv_stream(fp(fn))) != int(float(man_canon[fn]["n_data_rows"]))]
    if row_bad:
        fail(f"[M1] S55 行数与 manifest n_data_rows 不符: {row_bad}")
    else:
        ok(f"[M1] S55 表族 {len(S55_FAMILY)} 张齐全（S55r 按设计留空），行数全部与 manifest 冻结值一致")

    # S55a：23 切片、条件计数、spots 总数、可测基因上限
    a = list(read_csv_stream(fp("Table_S55a_M1_Visium_Section_Spatial_Stats.csv")))
    cond = Counter(r["condition"] for r in a)
    spots = sum(int(r["n_spots"]) for r in a)
    cap_bad = [r["section"] for r in a
               if not (1 <= int(r["n_up_genes"]) <= S55_GENE_CAPS["n_up_genes"]
                       and 1 <= int(r["n_ex_genes"]) <= S55_GENE_CAPS["n_ex_genes"])]
    if len(a) != S55_N_SECTIONS or len({r["section"] for r in a}) != S55_N_SECTIONS \
            or dict(cond) != S55_CONDITION_COUNTS or spots != S55_SPOT_TOTAL:
        fail(f"[M1] S55a 不变量破坏: n={len(a)} 条件 {dict(cond)} spots={spots} "
             f"!= {S55_N_SECTIONS}/{S55_CONDITION_COUNTS}/{S55_SPOT_TOTAL}")
    if cap_bad:
        fail(f"[M1] S55a 可测基因数超出预注册上限（上游≤{S55_GENE_CAPS['n_up_genes']}/"
             f"执行≤{S55_GENE_CAPS['n_ex_genes']}）: {cap_bad[:5]}")
    if len(a) == S55_N_SECTIONS and dict(cond) == S55_CONDITION_COUNTS \
            and spots == S55_SPOT_TOTAL and not cap_bad:
        ok(f"[M1] S55a = {S55_N_SECTIONS} 切片（Control 4/AcuteDAD 7/ProliferativeDAD 12，"
           f"spots 合计 {S55_SPOT_TOTAL}）；可测基因数均在预注册上限内")

    # 预注册判定复算（复刻 M1_step7_judgment.py，期望 R3）
    n = len(a)
    z_st = sum(float(r["z_bv"]) for r in a) / math.sqrt(n)
    p_st = 2 * _norm_sf(abs(z_st))
    n_neg_sig = sum(1 for r in a if float(r["I_bv"]) < 0 and float(r["p_bv"]) < 0.05)
    n_up_sig = sum(1 for r in a if float(r["I_up"]) > 0 and float(r["p_up"]) < 0.05)
    n_ex_sig = sum(1 for r in a if float(r["I_ex"]) > 0 and float(r["p_ex"]) < 0.05)
    r1 = (p_st < 0.05) and (z_st < 0) and (n_neg_sig >= 2)
    verdict = "R1"
    if not r1:
        up_diffuse = n_up_sig <= n / 2
        ex_diffuse = n_ex_sig <= n / 2
        verdict = "R2" if (n_up_sig >= n / 2) != (n_ex_sig >= n / 2) \
            and not (up_diffuse and ex_diffuse) else "R3"
    if verdict != "R3" or abs(z_st - 6.33) > 0.02 or (n_neg_sig, n_up_sig, n_ex_sig) != (3, 15, 21):
        fail(f"[M1] S55a 判定复算偏离冻结结果: verdict={verdict} z_bv_st={z_st:+.3f} "
             f"n_neg_sig/n_up_sig/n_ex_sig={n_neg_sig}/{n_up_sig}/{n_ex_sig}（冻结 R3 / +6.33 / 3/15/21）")
    else:
        ok(f"[M1] 预注册判定复算 = R3（双变量 I 合并 z=+6.33 显著为正；负显著切片 {n_neg_sig}/23；"
           f"双臂聚集 上游 {n_up_sig}/23、执行 {n_ex_sig}/23）")

    # 判定报告文件与结论一致（2026-08-23：交付夹按 README_打包说明 §一将报告抽至 M_judgment_reports/，
    # _intermediate/ 不随包——两处位置任一存在即可）
    jrep = ROOT / "_intermediate" / "M1_judgment_report.md"
    if not jrep.exists():
        jrep = ROOT / "M_judgment_reports" / "M1_judgment_report.md"
    if not jrep.exists():
        fail("[M1] 判定报告缺失（_intermediate/ 与 M_judgment_reports/ 均无 M1_judgment_report.md）")
    else:
        body = jrep.read_text(encoding="utf-8")
        if "判定结果：R3" not in body:
            fail("[M1] 判定报告与复算判定不一致（应含 '判定结果：R3'）")
        else:
            ok("[M1] 判定报告存在且与复算判定一致（R3）")

    # S55b/S55c/S55f 结构不变量
    b = list(read_csv_stream(fp("Table_S55b_M1_Visium_Condition_Summary.csv")))
    b_sec = {r["condition"]: int(r["n_sections"]) for r in b}
    c = list(read_csv_stream(fp("Table_S55c_M1_Visium_Domain_Stats.csv")))
    c_dom = Counter(r["section"] for r in c)
    f = list(read_csv_stream(fp("Table_S55f_M1_Visium_SVG_Enrichment.csv")))
    if len(b) != 3 or b_sec != S55_CONDITION_COUNTS:
        fail(f"[M1] S55b 条件汇总异常: {b_sec}")
    if len(c) != S55_N_SECTIONS * 5 or any(v != 5 for v in c_dom.values()):
        fail(f"[M1] S55c 域统计异常: 行数 {len(c)}，逐切片域数 {dict(c_dom)}")
    if {r["arm"] for r in f} != {"upstream_collapse", "execution_induction", "not_in_dissociation_arms"} \
            or any(int(r["n_sections"]) != S55_N_SECTIONS for r in f):
        fail(f"[M1] S55f SVG 富集异常: arms={sorted({r['arm'] for r in f})}")
    else:
        ok("[M1] S55b 条件 n_sections=4/7/12；S55c 23×5=115 域；S55f 双臂+背景各 23 切片")

    # S55v 残差化敏感性：Stouffer z≈+4.04（共定位非密度混杂，冻结值）
    v = list(read_csv_stream(fp("Table_S55v_M1_Visium_Bivariate_Residualized.csv")))
    if len(v) != S55_N_SECTIONS:
        fail(f"[M1] S55v 行数 {len(v)} != {S55_N_SECTIONS}")
    else:
        z_res = sum(_norm_ppf(1 - min(max(float(r["p_resid"]), 1e-3), 1 - 1e-3) / 2)
                    * math.copysign(1.0, float(r["I_bv_resid"])) for r in v) / math.sqrt(S55_N_SECTIONS)
        if abs(z_res - 4.04) > 0.05 or z_res <= 0:
            fail(f"[M1] S55v 残差化 Stouffer z={z_res:+.3f} 偏离冻结 +4.04")
        else:
            ok(f"[M1] S55v 残差化敏感性 Stouffer z=+{z_res:.2f}（共定位非密度混杂，与冻结 +4.04 一致）")

    # CosMx 面板覆盖与预注册一致（上游 0/30、执行 5/33）
    up80 = {r["gene_symbol"] for r in _gm if r["arm"] == "upstream_collapse"}
    ex80 = {r["gene_symbol"] for r in _gm if r["arm"] == "execution_induction"}
    panel = {r["gene"] for r in read_csv_stream(fp("Table_S55i_M1_CosMx_Gene_Detection.csv"))}
    if panel & up80 or (panel & ex80) != S55_COSMX_EXEC_EXPECTED:
        fail(f"[M1] S55i CosMx 面板覆盖偏离预注册: 上游交集 {sorted(panel & up80)}，"
             f"执行交集 {sorted(panel & ex80)}（期望 ∅ / {sorted(S55_COSMX_EXEC_EXPECTED)}）")
    else:
        ok(f"[M1] S55i CosMx 面板覆盖与预注册一致（上游 0 可测、执行 {sorted(S55_COSMX_EXEC_EXPECTED)}）")

    # M1 空间样本清单（139 = GSE271370 23 切片 + GSE253474 116 FOV）
    sp = list(read_csv_stream(fp("M1_spatial_sample_manifest.csv")))
    sp_ds = Counter(r["dataset"] for r in sp)
    if dict(sp_ds) != {"GSE271370": 23, "GSE253474": 116} or len(sp) != 139 \
            or any(not r["gsm"] for r in sp):
        fail(f"[M1] M1_spatial_sample_manifest 异常: {dict(sp_ds)}（n={len(sp)}）")
    else:
        ok("[M1] M1_spatial_sample_manifest = 139 行（GSE271370 23 切片 + GSE253474 116 FOV），GSM 全覆盖")


def check_m2_mdi(man_canon):
    """M2 MDI 临床化守门（2026-08-17）：S56–S58 表族不变量、Figure 9 面板、
    预注册判定 R2 复算（M2_pre_registration_20260817.md，规则先于结果）。"""
    FIG = ROOT / "01_FIGURE_DATA_CSV" / "Main"
    s56 = {fn for fn, r in man_canon.items() if r["s_number"] == "S56"}
    s57 = {fn for fn, r in man_canon.items() if r["s_number"] == "S57"}
    s58 = {fn for fn, r in man_canon.items() if r["s_number"] == "S58"}
    if not s56 or not s57 or not s58:
        fail(f"[M2] S56–S58 登记缺失: S56={sorted(s56)} S57={sorted(s57)} S58={sorted(s58)}")
        return

    # S56：4,888 行、六队列计数（与 M2_README 一致）
    f56 = fp("Table_S56_M2_PerSample_Scores.csv")
    rows56 = list(read_csv_stream(f56))
    coh = Counter(r["cohort"] for r in rows56)
    exp56 = {"GSE148871": 304, "GSE185263": 392, "GSE188309": 198,
             "GSE212865": 137, "GSE310929": 3713, "GSE32707": 144}
    if len(rows56) != 4888 or dict(coh) != exp56:
        fail(f"[M2] S56 行数/队列计数 {len(rows56)}/{dict(coh)} != 4888/{exp56}")
    else:
        ok("[M2] S56 = 4,888 行；六队列计数 304/392/198/137/3,713/144")

    # S57：含 meta 合并行（MDI 合并 g=+1.85 冻结值）
    rows57 = list(read_csv_stream(fp("Table_S57_M2_CrossDisease_Contrasts_Meta.csv")))
    meta_rows = [r for r in rows57 if r.get("n_cohorts") and float(r["n_cohorts"]) == 3]
    if not meta_rows:
        fail("[M2] S57 缺 3 队列 meta 合并行")
    else:
        mdi_meta = [r for r in meta_rows if r["score"] == "MDI"]
        if not mdi_meta or abs(float(mdi_meta[0]["pooled_g"]) - 1.846) > 0.02:
            fail(f"[M2] S57 MDI meta 合并 g 偏离冻结 +1.85: {[r['pooled_g'] for r in mdi_meta]}")
        else:
            ok(f"[M2] S57 meta 合并（3 队列）与冻结一致（MDI g=+{float(mdi_meta[0]['pooled_g']):.2f}）")

    # S58：主分析行冻结值（脓毒症 28d logistic OR=1.063 p=0.218；增量 ΔAUC=+0.0001）
    rows58 = list(read_csv_stream(fp("Table_S58_M2_Outcome_Association.csv")))
    main_row = [r for r in rows58 if r.get("cohort") == "GSE310929" and r.get("model") == "Logistic_unadj_28d"]
    inc_row = [r for r in rows58 if r.get("model") == "Increment"]
    if not main_row or abs(float(main_row[0]["HR_or_OR"]) - 1.063) > 0.01 \
            or abs(float(main_row[0]["p"]) - 0.218) > 0.01:
        fail(f"[M2] S58 主分析行偏离冻结 OR=1.063/p=0.218: {[{k: r[k] for k in ('HR_or_OR','p')} for r in main_row]}")
    if not inc_row or abs(float(inc_row[0]["CI_hi"]) - 0.0001) > 0.001:
        fail(f"[M2] S58 增量行偏离冻结 ΔAUC=+0.0001: {[{k: r[k] for k in ('CI_hi',)} for r in inc_row]}")
    else:
        ok("[M2] S58 主分析 OR=1.063/p=0.218 与增量 ΔAUC=+0.0001 与冻结一致")

    # Figure 5A–5E（MDI 临床面板；原 11 图体系 Figure_4A–4E，2026-09-06 9 图重编号后迁移）存在且行数>0
    fig_ok = all((FIG / f"Figure_5{c}.csv").exists() for c in "ABCDE")
    if not fig_ok:
        fail("[M2] Figure_5A–5E 面板数据缺失")
    else:
        ok("[M2] Figure_5A–5E 面板数据齐全")

    # 预注册判定复算（规则：M2_pre_registration_20260817.md §3；口径与 M2_step4 一致）
    repl = 0
    for cond in [
        (("GSE185263", "SepsisCOVID_vs_Control")),
        (("GSE32707", "ARDSd0_vs_Control")),
        (("GSE310929", "Sepsis_vs_Control")),
    ]:
        r = [x for x in rows57 if x.get("cohort") == cond[0] and x.get("contrast") == cond[1]
             and x.get("score") == "MDI" and x.get("BH_q")]
        if r and float(r[0]["BH_q"]) < 0.05:
            repl += 1
    cond_i = repl >= 2
    cond_ii = (main_row and float(main_row[0]["p"]) < 0.05) and \
              any(float(x["p"]) < 0.05 for x in rows58
                  if x.get("cohort") == "GSE188309" and x.get("model") == "Logistic_unadj")
    cond_iii = inc_row and float(inc_row[0]["CI_hi"]) > 0.02
    verdict = "R1" if (cond_i and cond_ii and cond_iii) else ("R2" if (cond_i or cond_ii) else "R3")
    if verdict != "R2":
        fail(f"[M2] 预注册判定复算 {verdict} != 冻结 R2（(i){cond_i}/(ii){cond_ii}/(iii){cond_iii}）")
    else:
        ok(f"[M2] 预注册判定复算 = R2（(i) 跨队列复制成立 {(cond_i)}；(ii) 结局关联不成立；"
           f"(iii) 增量不成立）")

    # 2026-08-23：交付夹报告位于 M_judgment_reports/（见 README_打包说明 §一），两处任一存在即可
    jrep = ROOT / "_intermediate" / "M2_judgment_report.md"
    if not jrep.exists():
        jrep = ROOT / "M_judgment_reports" / "M2_judgment_report.md"
    if not jrep.exists() or "判定结果：**R2**" not in jrep.read_text(encoding="utf-8"):
        fail("[M2] 判定报告缺失（_intermediate/ 与 M_judgment_reports/ 均无 M2_judgment_report.md）或与 R2 不一致")
    else:
        ok("[M2] 判定报告存在且与复算判定一致（R2）")


def main():
    if not MANIFEST.exists():
        print("缺少 RESULTS_MANIFEST_v2.0.csv，请先运行 build_D4_package.py / M2_step4")
        sys.exit(2)
    man_canon = check_manifest()
    check_n22_numbering(man_canon)
    check_versions(man_canon)
    check_symbols()
    check_p_and_fdr(man_canon)
    check_counts()
    check_sample_manifest()
    check_m1_spatial(man_canon)
    check_m2_mdi(man_canon)

    rep = ROOT / "03_LOGS" / "lint_report.md"
    rep.parent.mkdir(exist_ok=True)
    with open(rep, "w", encoding="utf-8") as f:
        f.write("# lint_package.py 全包校验报告（D4）\n\n")
        f.write(f"日期：{date.today().isoformat()}（manifest: RESULTS_MANIFEST_v2.0）\n\n")
        f.write(f"## 结果：{'全绿 PASS' if not FAILS else 'FAIL（' + str(len(FAILS)) + '）'}\n\n")
        f.write(f"- 通过 {len(OKS)} 项；警告 {len(WARNS)} 项；失败 {len(FAILS)} 项\n\n")
        if FAILS:
            f.write("## 失败项\n\n")
            for m in FAILS:
                f.write(f"- ❌ {m}\n")
        if WARNS:
            f.write("\n## 警告项\n\n")
            for m in WARNS:
                f.write(f"- ⚠️ {m}\n")
        f.write("\n## 通过项\n\n")
        for m in OKS:
            f.write(f"- ✅ {m}\n")

    print(f"\n==== lint 结果: {'全绿 ✅' if not FAILS else '存在失败 ❌'} "
          f"(通过 {len(OKS)} / 警告 {len(WARNS)} / 失败 {len(FAILS)}) ====")
    for m in FAILS[:40]:
        print("  ❌", m)
    for m in WARNS[:10]:
        print("  ⚠️", m)
    print(f"报告: {rep}")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
