# -*- coding: utf-8 -*-
"""
D4/N22 治理工程化构建脚本（2026-08-15）
=========================================
对应质控报告 D4（治理工程化）+ N22（补充表编号冲突与版本并存）条目。

执行内容（幂等，可重复运行，已应用的步骤自动跳过）：
  1. N22 编号去重：
     - S6 双占用：Scissor 保持 S6（§8.1 原始登记/主图 Figure 7），单细胞解离检验 → S48
       （Table_S6_Dissociation_*.csv → Table_S48_Dissociation_*.csv）
     - S3 双占用：文稿引用（line 139）的细胞类型组成保持 S3，WGCNA 模块 → S49
       （Table_S3_WGCNA_Modules.csv → Table_S49_WGCNA_Modules.csv）
     - Table_s12_Mitoxyperilysis_Upstream_TF.csv（小写 s）→ Table_S12b_*
     - S26b 双占用（方案 §26 建议）：HLCA 注释 → S26a
     - S18 三对 _CORRECTED 与原版内容已逐字节一致（2026-08-15 修正时同步）：
       _CORRECTED 晋升为正式表名，原非修正版移入 archive/
  2. 历史备份/被取代文件 → archive/
  3. R1 遗留① 历史列名统一（数值零改动）：
     - Table_S5b/S5c：ARDS_z → Sepsis_COVID_z 等（GSE185263 sepcv 组 R1 正名）
     - Table_S9_Model_Performance.csv task 值 + Table_S9_Diagnostic_Report.txt
     - Bridge_Test_Result.csv：lung_log2FC → GSE185263_log2FC、blood_log2FC → GSE212865_log2FC
  4. canonical 基因 manifest：Mitoxyperilysis_Gene_Manifest_v1.0.csv
     （80 基因 + HGNC 别名 + 模块 + 解离双臂归属 + 版本号）
  5. 全部 canonical 结果表加 gene_set_version / score_version 列（常量列，值按冻结配置）
  6. 冻结结果表 manifest：RESULTS_MANIFEST_v1.0.csv（S 编号↔文件唯一映射 + sha256 校验和）
     + RESULTS_MANIFEST.md（人读版：决策记录 + 审计轨迹）

不改动任何统计数值。运行后请运行 lint_package.py 校验至全绿。
"""
import csv
import hashlib
import io
import json
import os
import shutil
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ROOT = Path(__file__).resolve().parent
GOV = ROOT / "04_AUDIT_GOVERNANCE"
TAB = ROOT / "02_SUPPLEMENTARY_TABLES" / "SUPPLEMENTARY_Tables_CSV"
NOTES = ROOT / "02_SUPPLEMENTARY_TABLES" / "Supplementary_Notes"
ARCH = ROOT / "archive"

# 数据存放层输入（00_RAW_DATA 内，显式登记，不整目录扫描）
RAW_REGISTRY = {"GSE67530_beta_for_GrimAge.csv.gz": "00_RAW_DATA"}
# 登记在册但本地不随包的文件（留痕；与 lint_package.py MISSING_OK 白名单同源）
ABSENT_REGISTRY = {
    "Table_S11_Mitoxyperilysis_TF_PerCell.csv":
        ("canonical", "数据存放层大矩阵，本地不随包（Zenodo/GEO 存档），见 01_FIGURE_DATA_CSV/README_图数据包说明.md"),
    "Table_S11_TF_Activity_Matrix_AllCells.csv":
        ("canonical", "数据存放层大矩阵，本地不随包（Zenodo/GEO 存档），见 01_FIGURE_DATA_CSV/README_图数据包说明.md"),
    "Table_S1_Mitoxyperilysis_Score_by_CellType_prev38.csv":
        ("archived", "历史版本，已被 canonical Table_S1 取代；2026-08-16 目录重组确认原始文件未随包迁移（见 01_FIGURE_DATA_CSV/README_图数据包说明.md 治理披露）"),
    "Table_S1_Mitoxyperilysis_Score_by_CellType_Summary_prev38.csv":
        ("archived", "历史版本，已被 canonical Table_S1 取代；2026-08-16 目录重组确认原始文件未随包迁移（见 01_FIGURE_DATA_CSV/README_图数据包说明.md 治理披露）"),
}


def d4_path(fname):
    """按 manifest 的 location 列解析文件位置（治理层是唯一路径来源；v2.0 为现役 manifest，v1.0 留作回退）。"""
    for mf in ("RESULTS_MANIFEST_v2.0.csv", "RESULTS_MANIFEST_v1.0.csv"):
        m = GOV / mf
        if not m.exists():
            continue
        try:
            with open(m, encoding="utf-8-sig", newline="") as f:
                for r in csv.DictReader(f):
                    if r["filename"] == fname:
                        loc = r.get("location", "")
                        return (ROOT / loc / fname) if loc else (ROOT / fname)
        except Exception:
            continue
    return GOV / fname

GS_V = "Mitoxy-80_v1.0"      # 80 基因通路清单（canonical，常用名符号 IP3R1/HSP60）
PCD_V = "PCD-6panel_v1.0"    # S5 六种程序性死亡评分基因面板
SV_SG = "score_genes_v1"     # scanpy score_genes 逐细胞评分
SV_SSG = "ssGSEA_v1"         # running-sum ssGSEA（β 参数见表内列）
SV_BR = "bridge_v2"          # Bridge 五维收敛（v2 = N1 别名归一化修复后）

LOG = []


def L(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s, flush=True)


# --------------------------------------------------------------------------
# 1) 归档移动（备份 + 被取代版本）—— old -> archive 内文件名
# --------------------------------------------------------------------------
ARCHIVE_MOVES = {
    "Bridge_Test_Result_pre_N1_backup.csv": "Bridge_Test_Result_pre_N1_backup.csv",
    "GSE185263_groups_pre_R1_backup.csv": "GSE185263_groups_pre_R1_backup.csv",
    "Table_S2_DEGs_Analysis_pre_fix_backup.csv": "Table_S2_DEGs_Analysis_pre_fix_backup.csv",
    "Table_S5_PCD_Scores_RAW_pre_R1_backup.csv": "Table_S5_PCD_Scores_RAW_pre_R1_backup.csv",
    "Table_S5_PCD_Scores_pre_R1_backup.csv": "Table_S5_PCD_Scores_pre_R1_backup.csv",
    "Table_S47b_Circulating_Mitoxy_Immune_Correlations_pre_R2_backup.csv": "Table_S47b_Circulating_Mitoxy_Immune_Correlations_pre_R2_backup.csv",
    # S18 原非修正版（与 _CORRECTED 逐字节一致，2026-08-15 已核对；加 _superseded 后缀避免与晋升后的正式名同名）
    "Table_S18_AlphaMissense_Main.csv": "Table_S18_AlphaMissense_Main_superseded.csv",
    "Table_S18b_AlphaMissense_Module_Summary.csv": "Table_S18b_AlphaMissense_Module_Summary_superseded.csv",
    "Table_S18c_AlphaMissense_HighPriority.csv": "Table_S18c_AlphaMissense_HighPriority_superseded.csv",
}

# --------------------------------------------------------------------------
# 2) 编号去重重命名
# --------------------------------------------------------------------------
RENAMES = {
    "Table_S6_Dissociation_singlecell_correlation.csv": "Table_S48_Dissociation_singlecell_correlation.csv",
    "Table_S6_Method_Note.txt": "Table_S48_Dissociation_Method_Note.txt",
    "Table_S3_WGCNA_Modules.csv": "Table_S49_WGCNA_Modules.csv",
    "Table_s12_Mitoxyperilysis_Upstream_TF.csv": "Table_S12b_Mitoxyperilysis_Upstream_TF.csv",
    "Table_S26b_HLCA_CellType_Annotations.csv": "Table_S26a_HLCA_CellType_Annotations.csv",
    # S18 _CORRECTED 晋升为正式名（原版已先行移入 archive/）
    "Table_S18_AlphaMissense_Main_CORRECTED.csv": "Table_S18_AlphaMissense_Main.csv",
    "Table_S18b_AlphaMissense_Module_Summary_CORRECTED.csv": "Table_S18b_AlphaMissense_Module_Summary.csv",
    "Table_S18c_AlphaMissense_HighPriority_CORRECTED.csv": "Table_S18c_AlphaMissense_HighPriority.csv",
}

# --------------------------------------------------------------------------
# 3) 历史列名/取值统一（R1 遗留①；数值零改动）
# --------------------------------------------------------------------------
COL_RENAMES = {
    "Table_S5b_Mitoxy_80gene_Module_ssGSEA.csv": {
        "ARDS_z": "Sepsis_COVID_z",
        "ARDS_pct_pos": "Sepsis_COVID_pct_pos",
        "delta": "delta_Sepsis_COVID_vs_Control",
    },
    "Table_S5c_Mitoxyperilysis_sensitivity.csv": {
        "ARDS_z": "Sepsis_COVID_z",
        "delta_ARDS_Control": "delta_Sepsis_COVID_vs_Control",
    },
    "Bridge_Test_Result.csv": {
        "lung_log2FC": "GSE185263_log2FC",   # S46a 来源，历史误标"肺组织"，R1 已证 GSE185263 为全血
        "blood_log2FC": "GSE212865_log2FC",  # S46b 来源，按数据集锚定命名
    },
}
VALUE_RENAMES = {
    "Table_S9_Model_Performance.csv": {
        "task": {
            "ARDS_vs_Control": "Sepsis_COVID_vs_Control",
            "Sepsis_vs_ARDS": "Sepsis_vs_Sepsis_COVID",   # 方向与文稿一致：脓毒症 vs 脓毒症合并COVID-19
        },
        "note": {
            "主任务;预测脓毒症进展为ARDS;具临床意义":
                "主任务;区分脓毒症是否合并COVID-19(sepcv,R1正名);具临床意义",
        },
    },
}
TEXT_REPLACEMENTS = {
    "Table_S9_Diagnostic_Report.txt": [
        ("原模型 = ARDS肺 vs 健康肺", "原模型 = Sepsis_COVID 全血 vs 健康对照（GSE185263 sepcv，R1 正名）"),
        ("Task A：ARDS vs Control", "Task A：Sepsis_COVID vs Control"),
        ("Task B：ARDS vs Sepsis（主任务，有意义）】预测脓毒症患者进展为ARDS",
         "Task B：Sepsis vs Sepsis_COVID（主任务，有意义）】区分脓毒症是否合并COVID-19"),
        ("主模型改用 Task B（ARDS vs Sepsis，预测ARDS进展）", "主模型改用 Task B（Sepsis vs Sepsis_COVID）"),
        ("'prognostic/risk' -> 'prediction of ARDS development in sepsis'。",
         "'prognostic/risk' -> '脓毒症是否合并COVID-19 的判别分类'（R1 正名后定位）。"),
    ],
    "Table_S48_Dissociation_Method_Note.txt": [
        ("Table S6 单细胞解离检验", "Table S48 单细胞解离检验"),
        ("Table_S6 ..._correlation.csv", "Table_S48 ..._correlation.csv"),
        ("与 bulk Table_S5b 一致", "与 bulk Table_S5b 一致"),
    ],
    "Table_S18_AlphaMissense_Analysis_Report.txt": [
        ("Table_S18_AlphaMissense_Main_CORRECTED.csv", "Table_S18_AlphaMissense_Main.csv"),
        ("Table_S18b_AlphaMissense_Module_Summary_CORRECTED.csv", "Table_S18b_AlphaMissense_Module_Summary.csv"),
        ("Table_S18c_AlphaMissense_HighPriority_CORRECTED.csv", "Table_S18c_AlphaMissense_HighPriority.csv"),
    ],
}


# --------------------------------------------------------------------------
# CSV 基础工具（保留 BOM 与行尾风格；仅追加常量列/改表头，不触碰数值文本）
# --------------------------------------------------------------------------
def _open_info(path):
    raw = path.read_bytes()
    bom = raw.startswith(b"\xef\xbb\xbf")
    enc = "utf-8-sig" if bom else "utf-8"
    crlf = b"\r\n" in raw[:65536]
    return enc, ("\r\n" if crlf else "\n")


def read_csv_rows(path):
    enc, nl = _open_info(path)
    with open(path, "r", encoding=enc, newline="") as f:
        return list(csv.reader(f)), enc, nl


def write_csv_rows(path, rows, enc, nl):
    with open(path, "w", encoding=enc, newline="") as f:
        w = csv.writer(f, lineterminator=nl)
        w.writerows(rows)


def apply_col_renames(fname, renames):
    p = d4_path(fname)
    if not p.exists():
        L(f"  [跳过-不存在] {fname}")
        return
    rows, enc, nl = read_csv_rows(p)
    hdr = rows[0]
    if all(old not in hdr for old in renames):
        L(f"  [跳过-已应用] {fname}")
        return
    idx = {old: hdr.index(old) for old in renames if old in hdr}
    for old, new in renames.items():
        if old in idx:
            hdr[idx[old]] = new
    write_csv_rows(p, rows, enc, nl)
    L(f"  [列名] {fname}: " + ", ".join(f"{o}->{n}" for o, n in renames.items() if o in idx))


def apply_value_renames(fname, col_map):
    p = d4_path(fname)
    if not p.exists():
        L(f"  [跳过-不存在] {fname}")
        return
    rows, enc, nl = read_csv_rows(p)
    hdr = rows[0]
    changed = False
    for col, vmap in col_map.items():
        if col not in hdr:
            continue
        i = hdr.index(col)
        for r in rows[1:]:
            if i < len(r) and r[i] in vmap:
                r[i] = vmap[r[i]]
                changed = True
    if not changed:
        L(f"  [跳过-已应用] {fname}")
        return
    write_csv_rows(p, rows, enc, nl)
    L(f"  [取值] {fname}: task 标签 ARDS→Sepsis_COVID")


def apply_text_replacements(fname, pairs):
    p = d4_path(fname)
    if not p.exists():
        L(f"  [跳过-不存在] {fname}")
        return
    txt = p.read_text(encoding="utf-8")
    orig = txt
    for a, b in pairs:
        txt = txt.replace(a, b)
    if txt == orig:
        L(f"  [跳过-已应用] {fname}")
        return
    p.write_text(txt, encoding="utf-8")
    L(f"  [文本] {fname}: {len(pairs)} 处替换")


# --------------------------------------------------------------------------
# 4) canonical 基因 manifest
# --------------------------------------------------------------------------
ALIAS_TO_HGNC = {"IP3R1": "ITPR1", "HSP60": "HSPD1"}   # canonical -> HGNC 官方
# 包内各数据层实际出现过的别名/变体符号（canonical -> 该符号在各层中的写法）
PKG_ALIASES = {
    "IP3R1": "ITPR1",
    "HSP60": "HSPD1",
    "HSPA9": "GRP75",     # ATAC/转录本层写法
    "CYCS": "CYTC",       # ATAC/转录本层写法
    "SLC11A2": "DMT1",    # MR 层写法
    "GSDME": "DFNA5",     # GSE32707/GPL10558 注释层写法（D1b 2026-08-16 登记）
}
ARM_BY_MODULE = {
    "mitoxy_MAM_integrity": "upstream_collapse",
    "mitoxy_mitochondrial_function": "upstream_collapse",
    "mitoxy_ferroptosis_cuproptosis": "upstream_collapse",
    "mitoxy_iron_metabolism": "execution_induction",
    "mitoxy_cell_death": "execution_induction",
    "mitoxy_oxidative_stress": "execution_induction",
    "mitoxy_transcription_factors": "not_in_dissociation_arms",
    "mitoxy_autophagy": "not_in_dissociation_arms",
}


def build_gene_manifest():
    src = GOV / "Mitoxyperilysis_Pathway_Gene_List.csv"
    out = GOV / "Mitoxyperilysis_Gene_Manifest_v1.0.csv"
    rows, enc, nl = read_csv_rows(src)
    hdr = rows[0]
    ci = {c: hdr.index(c) for c in hdr}
    sizes = {}
    for r in rows[1:]:
        sizes[r[ci["Mitoxyperilysis_module"]]] = sizes.get(r[ci["Mitoxyperilysis_module"]], 0) + 1
    new_rows = [["version", "gene_symbol", "hgnc_symbol", "aliases", "ensembl_gene_id", "gene_name",
                 "module", "module_size", "arm", "is_mt_gene",
                 "ARDS_vs_Control_log2FC", "ARDS_vs_Control_padj"]]
    n = 0
    for r in rows[1:]:
        sym = r[ci["gene_symbol"]]
        mod = r[ci["Mitoxyperilysis_module"]]
        n += 1
        new_rows.append([
            "v1.0", sym, ALIAS_TO_HGNC.get(sym, sym), PKG_ALIASES.get(sym, ""),
            r[ci["ensembl_gene_id"]], r[ci["gene_name"]], mod, str(sizes[mod]),
            ARM_BY_MODULE.get(mod, ""), "1" if sym.startswith("MT-") else "0",
            r[ci["ARDS_vs_Control_log2FC"]], r[ci["ARDS_vs_Control_padj"]],
        ])
    write_csv_rows(out, new_rows, "utf-8-sig", "\r\n")
    L(f"  [基因manifest] {out.name}: {n} 基因, 模块尺寸={sizes}")
    assert n == 80 and sum(sizes.values()) == 80, "基因清单必须为 80 基因"
    return sizes


# --------------------------------------------------------------------------
# 5) 版本戳（全部 canonical CSV 追加两常量列；TXT 报告追加页脚行）
# --------------------------------------------------------------------------
# family 默认：前缀 -> (gene_set_version, score_version, gene_universe)
FAMILY_STAMP = [
    ("Table_S1_Mitoxyperilysis_Score", GS_V, SV_SG, "none"),
    ("Table_S2b_", GS_V, "NA", "mitoxy_80"),      # N3 真·单细胞 DEG（80 通路基因）
    ("Table_S8b_", GS_V, "NA", "mitoxy_80"),      # N3 诚实三平台 Fisher Meta
    ("Table_S4d_HubGene_Focused", GS_V, "NA", "genome_wide"),
    ("Table_S5_PCD_Scores", PCD_V, SV_SSG, "none"),
    ("Table_S5b_", GS_V, SV_SSG, "none"),
    ("Table_S5c_", GS_V, SV_SSG, "none"),
    ("Table_S6_Scissor_Score_Summary", GS_V, SV_SG, "none"),
    ("Table_S7h_", GS_V, "NA", "genome_wide"),
    ("Table_S7i_", GS_V, "NA", "none"),
    ("Table_S7q_", GS_V, SV_SG, "none"),
    ("Table_S7r_", GS_V, "NA", "none"),
    ("Table_S9_", GS_V, "NA", "mitoxy_80"),
    ("Table_S11_Mitoxyperilysis_TF_PerCell", GS_V, "NA", "none"),
    ("Table_S12b_", GS_V, "NA", "genome_wide"),
    ("Table_S13_", GS_V, "NA", "mitoxy_80"),
    ("Table_S14_NicheNet", GS_V, "NA", "genome_wide"),
    ("Table_S15a_", GS_V, "NA", "mitoxy_80_via_alias"),
    ("Table_S15b_", GS_V, "NA", "mitoxy_80_via_alias"),
    ("Table_S15c_", GS_V, "NA", "mitoxy_80_via_alias"),
    ("Table_S15d_", GS_V, "NA", "mitoxy_80"),
    ("Table_S16b_", GS_V, "NA", "none"),
    ("Table_S16d_", GS_V, "NA", "none"),
    ("Table_S16f_", GS_V, SV_SSG, "none"),
    ("Table_S18", GS_V, "NA", "custom_panel_AM55"),
    ("Table_S19", GS_V, "NA", "mitoxy_80"),
    ("Table_S20", GS_V, "NA", "genome_wide"),
    ("Table_S21", GS_V, "NA", "genome_wide"),
    ("Table_S22b_", GS_V, "NA", "none"),
    ("Table_S22_scFOCAL", GS_V, "NA", "none"),
    ("Table_S22_S23_", GS_V, "NA", "none"),
    ("Table_S23_", GS_V, "NA", "none"),
    ("Table_S24", GS_V, "NA", "none"),
    ("Table_S26d_", GS_V, "NA", "none"),
    ("Table_S26e_", GS_V, "NA", "none"),
    ("Table_S26f_", GS_V, "NA", "mitoxy_80"),
    ("Table_S26g_", GS_V, "NA", "none"),
    ("Table_S26h_", GS_V, "NA", "none"),
    ("Table_S26i_", GS_V, "NA", "none"),
    ("Table_S26j_", GS_V, "NA", "mitoxy_80"),
    ("Table_S26k_", GS_V, "NA", "none"),
    ("Table_S31b_", GS_V, "NA", "genome_wide"),
    ("Table_S32_", GS_V, "NA", "mitoxy_80"),
    ("Table_S33_", GS_V, "NA", "mitoxy_80"),
    ("Table_S36", GS_V, "NA", "mitoxy_80"),
    ("Table_S37_", GS_V, "NA", "mitoxy_80"),
    ("Table_S38_", GS_V, "NA", "mitoxy_80"),
    ("Table_S38b_", GS_V, "NA", "mitoxy_80"),
    ("Table_S39_", GS_V, "NA", "none"),
    ("Table_S40_", GS_V, "NA", "mitoxy_80_via_alias"),
    ("Table_S41a_", GS_V, "NA", "mitoxy_80_via_alias"),
    ("Table_S41b_", GS_V, "NA", "mitoxy_80_via_alias"),
    ("Table_S41c_", GS_V, "NA", "mitoxy_80_via_alias"),
    ("Table_S41d_", GS_V, "NA", "mitoxy_80_via_alias"),
    ("Table_S41e_", GS_V, "NA", "mitoxy_80_via_alias"),
    ("Table_S41f_", GS_V, "NA", "mitoxy_80_via_alias"),
    ("Table_S44b_", GS_V, "NA", "none"),
    ("Table_S46a_", GS_V, "NA", "legacy_38"),
    ("Table_S46b_", GS_V, "NA", "legacy_38"),
    ("Table_S46c_", GS_V, "NA", "legacy_38"),
    ("Table_S46d_", GS_V, "NA", "legacy_38"),
    ("Table_S46e_", GS_V, "NA", "legacy_38"),
    ("Table_S47b_", GS_V, "NA", "none"),
    ("Table_S48_", GS_V, SV_SG, "none"),
    ("Table_S49_", GS_V, "NA", "genome_wide"),
    ("Table_S51_D1b_", GS_V, SV_SSG, "mitoxy_80"),   # D1b GSE32707 模块层（74/80 在场）
    ("Table_S52_D1b_", GS_V, "NA", "mitoxy_80"),     # D1b 逐基因 log2FC（74 + 6 MT 缺席登记）
    ("Table_S53_D1b_", GS_V, SV_SSG, "mitoxy_80"),   # D1b 三锚点一致性矩阵
    ("Table_S54_D1b_", GS_V, SV_SSG, "mitoxy_80"),   # D1b 敏感性（β×3/对照受试者塌陷）
    ("Bridge_Test_Result", GS_V, SV_BR, "mitoxy_80"),
    ("Mitoxyperilysis_Pathway_Gene_List", GS_V, "NA", "mitoxy_80"),
    ("Mitoxyperilysis_Gene_Manifest", GS_V, "NA", "mitoxy_80"),
]
# 个别覆盖（family 规则覆盖不到或需微调的文件）
STAMP_OVERRIDES = {
    "Table_S2_DEGs_Analysis.csv": (GS_V, "NA", "genome_wide"),   # 含 Mitoxyperilysis_module 注释列
    "Table_S2_QC_Metrics.csv": (GS_V, "NA", "none"),             # N3 条件标签修复后的 QC 聚合表
    "Table_S13_ScTenifoldKnk_KO_Detailed.csv": (GS_V, "NA", "genome_wide"),  # KO 下调基因为全基因组 DR 基因
    # R5 处置（2026-08-15，选"标注"选项）：以下两表为旧版评分尺度（AUCell 式/20 量级），
    # 与现行 S1 score_genes 尺度不一致——标注 legacy，不用于正文数值引用
    "Table_S19h_CellType_Mitoxy_Score.csv": (GS_V, "legacy_scale_v0_not_for_text", "none"),
    "Table_S22b_scFOCAL_CellType_IC50.csv": (GS_V, "legacy_scale_v0_not_for_text", "none"),
}
# 不加版本列的文件（分组定义/输入数据/探索性输出/JSON）
STAMP_EXCLUDE = {
    "GSE185263_groups.csv",            # 分组定义，读取方按列 merge
    "MR_bio_CRP_Sepsis.csv",           # 探索性 MR 结局侧（状态待定，未入正文）
    "MR_bio_Ferritin_Sepsis.csv",
    "GSE67530_beta_for_GrimAge.csv.gz",  # 输入数据（甲基化 β 矩阵）
    "Table_S50_D1b_GSE32707_Groups.csv",  # D1b 分组定义（GSE185263_groups.csv 同角色）
}
# TXT 报告版本页脚（gene_set_version, score_version）
TXT_FOOTER = {
    "Table_S5_Method_Note.txt": (PCD_V, SV_SSG),
    "Table_S9_Diagnostic_Report.txt": (GS_V, "NA"),
    "Table_S18_AlphaMissense_Analysis_Report.txt": (GS_V, "NA"),
    "Table_S19_GWAS_Mitoxy_Integration_Report.txt": (GS_V, "NA"),
    "Table_S20_GRN_Causal_Perturbation_Report.txt": (GS_V, "NA"),
    "Table_S21_Combinatorial_Perturbation_Report.txt": (GS_V, "NA"),
    "Table_S22_S23_Virtual_Drug_Screening_Report.txt": (GS_V, "NA"),
    "Table_S22_scFOCAL_Analysis_Report.txt": (GS_V, "NA"),
    "Table_S24_Pathway_Factor_Analysis_Report.txt": (GS_V, "NA"),
    "Table_S26_Mitoxyperilysis_ATAC_Report.txt": (GS_V, "NA"),
    "Table_S38_S39_Geneformer_Report.txt": (GS_V, "NA"),
    "Table_S40_gnomAD_Analysis_Report.txt": (GS_V, "NA"),
    "Table_S41_Transcript_Analysis_Report.txt": (GS_V, "NA"),
    "Table_S46_Circulating_Transcriptome_Analysis_Report.txt": (GS_V, "NA"),
    "Table_S48_Dissociation_Method_Note.txt": (GS_V, SV_SG),
    "Table_S50_D1b_Analysis_Report.txt": (GS_V, SV_SSG),
}
FOOTER_FMT = "\n-- [版本戳 2026-08-15] gene_set_version={g}; score_version={s}; manifest=RESULTS_MANIFEST_v1.0 --\n"


def stamp_config_for(fname):
    if fname in STAMP_OVERRIDES:
        return STAMP_OVERRIDES[fname]
    for pref, g, s, u in FAMILY_STAMP:
        if fname.startswith(pref):
            return (g, s, u)
    return ("NA", "NA", "none")


def repair_comment_header_csv(fname, g, s):
    """修复带 '#' 注释首行且被早期版本戳误追加的 CSV（如 Table_S27b）：
    注释行去掉误加的 2 列；首个非注释行（真表头）把误加的值改为列名。"""
    p = d4_path(fname)
    if not p.exists():
        return
    rows, enc, nl = read_csv_rows(p)
    if not rows or not str(rows[0][0]).startswith("#"):
        return
    changed = False
    for r in rows:
        if str(r[0]).startswith("#"):
            if r[-2:] == ["gene_set_version", "score_version"]:
                del r[-2:]
                changed = True
        else:
            if r[-2:] == [g, s] and "gene_set_version" not in r[:-2]:
                r[-2:] = ["gene_set_version", "score_version"]
                changed = True
            break
    if changed:
        write_csv_rows(p, rows, enc, nl)
        L(f"  [修复] {fname}: 注释首行结构修正")


def stamp_all():
    n_csv = 0
    for p in sorted(list(TAB.iterdir()) + list(GOV.iterdir())):
        if not p.is_file():
            continue
        fname = p.name
        if fname in STAMP_EXCLUDE or not fname.lower().endswith(".csv"):
            continue
        if fname.startswith(("R1_", "N1_")) or fname in ("RESULTS_MANIFEST_v1.0.csv",):
            continue
        g, s, _u = stamp_config_for(fname)
        repair_comment_header_csv(fname, g, s)
        enc, nl = _open_info(p)
        tmp = p.with_suffix(p.suffix + ".tmp")
        with open(p, "r", encoding=enc, newline="") as fin, \
             open(tmp, "w", encoding=enc, newline="") as fout:
            w = csv.writer(fout, lineterminator=nl)
            first = True
            stamped = False
            for row in csv.reader(fin):
                if row and str(row[0]).startswith("#"):
                    w.writerow(row)          # '#' 注释行原样保留
                    continue
                if first:
                    first = False
                    if "gene_set_version" in row and "score_version" in row:
                        stamped = True
                    row = row + ["gene_set_version", "score_version"]
                else:
                    row = row + [g, s]
                w.writerow(row)
        if stamped:
            tmp.unlink()
            continue
        tmp.replace(p)
        n_csv += 1
    L(f"  [版本戳] 已为 {n_csv} 个 canonical CSV 追加 gene_set_version/score_version 列")

    n_txt = 0
    for fname, (g, s) in TXT_FOOTER.items():
        p = d4_path(fname)
        if not p.exists():
            L(f"  [警告] 页脚目标缺失: {fname}")
            continue
        raw = p.read_bytes()
        enc = "utf-8"
        try:
            txt = raw.decode("utf-8")
        except UnicodeDecodeError:
            enc = "gbk"          # 个别历史报告为 GBK（如 S22 scFOCAL）
            txt = raw.decode("gbk")
        if "版本戳 2026-08-15" in txt:
            continue
        p.write_bytes((txt.rstrip("\n") + "\n" + FOOTER_FMT.format(g=g, s=s).rstrip("\n") + "\n").encode(enc))
        n_txt += 1
    L(f"  [版本戳] 已为 {n_txt} 个 TXT 报告追加版本页脚")


# --------------------------------------------------------------------------
# 6) 冻结 RESULTS_MANIFEST
# --------------------------------------------------------------------------
S_META = {
    "S1": "5.1 各细胞类型 Mitoxyperilysis 评分（score_genes）",
    "S2": "步骤1 GSE185263 全血 Bulk DEG（N3 溯源：DESeq2 输出；ARDS_vs_Control=82 sepcv vs 44 ctrl、Sepsis_vs_Control=348 vs 44；附 QC）",
    "S3": "5.1 单细胞细胞类型组成（GSE158055 12 类型/64 亚群；GSE171668 仅 metadata）",
    "S4": "步骤2 WGCNA Hub 基因 + 富集（b/c ORA、d Mitoxy/铁聚焦）",
    "S5": "步骤3 6 种 PCD 评分（running-sum ssGSEA；b=80 基因模块 ssGSEA、c=β 敏感性）",
    "S6": "步骤6 Scissor+/- 标记基因与细胞比例（5 子表）",
    "S7": "步骤8/5.3 CellPhoneDB 细胞通讯 + RNA velocity 轨迹（a–r 子表）",
    "S8": "步骤9 跨数据集共识基因 / meta 分析（4 子表）",
    "S9": "步骤10 ML 风险模型（双任务：Sepsis_COVID_vs_Sepsis 主任务；LASSO/EN/RF/Stacking）",
    "S11": "步骤12 SCENIC TF 活性（GRNBoost 邻接，未做 cisTarget 修剪；含 85MB/688MB 大表）",
    "S12": "步骤12 SCENIC regulon + b=上游 TF（23 TF）",
    "S13": "步骤13 ScTenifoldKnk 虚拟敲除（12 KO 基因）",
    "S14": "步骤14 配体-受体-靶标链（CellPhoneDB+SCENIC；b–e 统计）",
    "S15": "§13.2 cis-only TWMR + 共定位（a 脓毒症、b FinnGen ARDS、c coloc、d cis/trans 审计）",
    "S16": "步骤16 CIBERSORTx 免疫浸润（b–g：相关/分组/样本信息/评分）",
    "S17": "步骤17 细胞组成差异（样本级置换检验；b–e）",
    "S18": "步骤18 AlphaMissense 突变致病性（55 基因；_CORRECTED 已晋升为正式名）",
    "S19": "步骤19 GWAS-Mitoxy 整合（c–h；Bonferroni/Simes）",
    "S20": "步骤20 共表达因果 GRN 扰动（a–l，无 k；非 CellOracle）",
    "S21": "步骤21 组合扰动协同（共表达 disruption；非 GEARS）",
    "S22": "步骤22 scFOCAL 单细胞药物敏感性 + S22/S23 整合药物表",
    "S23": "步骤23 CMap 靶向药物（19 种，并入 S22_S23 整合表）",
    "S24": "步骤24 通路因子分析（FA+PCA，非 MOFA+；b–f）",
    "S26": "步骤26 scATAC EpiAgent + Mitoxyperilysis 染色质可及性（a=HLCA 注释、b–k）",
    "S27": "步骤27 表观遗传时钟（b 组间比较、c GrimAge）",
    "S28": "步骤28 表观衰老-临床相关",
    "S29": "步骤29 scATAC peaks（S29–S31 共享报告）",
    "S30": "步骤30 motif 富集",
    "S31": "步骤31 SCENIC+ GRN（b=TF 排名）",
    "S32": "步骤32 TF 结合预测（S32–S33 共享方法注）",
    "S33": "步骤33 差异 TF 结合",
    "S36": "步骤36 分子对接（b 配体性质、c 靶蛋白）",
    "S37": "步骤37 分子动力学轨迹",
    "S38": "步骤38 Geneformer 扰动（b 详情、c 零样本分类；S38–S39 共享报告）",
    "S39": "步骤39 跨物种验证（小鼠直系同源）",
    "S40": "步骤40 gnomAD 纯化选择约束（Publication=发表用精简版）",
    "S41": "步骤41 转录本分析（a–g：DTU/剪接/外显子/多样性指数）",
    "S44": "步骤44 hdWGCNA 模块（b 模块-表型相关）",
    "S45": "步骤45 hdWGCNA Hub 基因",
    "S46": "步骤46 循环转录组（a GSE185263 / b GSE212865 / c–e 矩阵与汇总）",
    "S47": "步骤47 循环免疫浸润（a fractions、b 与 Mitoxy 相关）",
    "S48": "单细胞解离检验（N22 编号去重：原 Table_S6_Dissociation_*；score_genes 双臂相关）",
    "S49": "WGCNA 模块（N22 编号去重：原 Table_S3_WGCNA_Modules；6 色）",
    "S50": "D1b GSE32707 全血 ARDS 队列外部验证（Groups=144 样本分组、Analysis_Report=方法+预注册判定；一致性矩阵=S53）",
}
OTHER_META = {
    "Bridge_Test_Result.csv": ("bridge", "§1 五维收敛打分（bridge_v2=N1 别名修复后）"),
    "Mitoxyperilysis_Pathway_Gene_List.csv": ("gene_set_source", "80 基因权威清单（v1.0 源定义；ARDS_vs_Control_* 列=GSE185263 全血 Bulk，N3 溯源）"),
    "Table_S2b_scRNA_DEG_80genes.csv": ("supplementary_table", "N3 真·单细胞 DEG（80 通路基因；GSE145926 BALF / GSE158055 PBMC / 合并，Mann-Whitney+BH）"),
    "Table_S8b_80gene_ThreePlatform_Meta.csv": ("supplementary_table", "N3 诚实三平台 Fisher Meta（scRNA 合并 + GSE185263 Bulk + GSE212865 微阵列；33/80 共识）"),
    "N3_S2_DEG_label_check.csv": ("audit_output", "N3 S2 标签一致性核对（存储 ARDS_vs_Control 列 vs 82/44 设计重验）"),
    "N3_Sample_Composition_Report.md": ("audit_report", "N3 样本构成核查 + S2 溯源 + scFOCAL 气道子集（12 BALF + 22 sputum）报告（2026-08-15）"),
    "RESULTS_MANIFEST.md": ("manifest_doc", "结果表包冻结说明（决策记录 + 例外登记 + lint 不变量）"),
    "Mitoxyperilysis_Gene_Manifest_v1.0.csv": ("gene_set_source", "canonical 基因 manifest v1.0（80 基因+HGNC 别名+模块+双臂）"),
    "GSE185263_groups.csv": ("group_definition", "GSE185263 分组定义（R1：Sepsis 266/Sepsis_COVID 82/Control 44）"),
    "MR_bio_CRP_Sepsis.csv": ("exploratory", "生物标志物→脓毒症 MR 结局侧（§13 未展开，状态待定）"),
    "MR_bio_Ferritin_Sepsis.csv": ("exploratory", "生物标志物→脓毒症 MR 结局侧（§13 未展开，状态待定）"),
    "GSE67530_beta_for_GrimAge.csv.gz": ("input_data", "GSE67530 甲基化 β 矩阵（GrimAge 输入，非结果表）"),
    "README_5.1_Scrna_Basic.md": ("report", "5.1 scRNA 基础管线说明（2026-08-27 清理批次已归档至 archive/R1_docs，路径按 manifest location 解析）"),
    "N1_Bridge_Alias_Reconnection_Report.md": ("audit_report", "N1 别名重连审计报告（2026-08-15）"),
    "R1_GSE185263_phenotype_provenance_report.md": ("audit_report", "R1 表型 provenance 审计报告（2026-08-15）"),
    "D1b_GSE32707_Metadata_Check_Report.md": ("audit_report", "D1b 第 1 步元数据/平台兼容性核验报告（2026-08-16）"),
    "Table_S50_D1b_GSE32707_Groups.csv": ("group_definition", "GSE32707 分组定义（D1b：Control 34/SIRS_d0 21/Sepsis_d0 30/Sepsis_d7 28/ARDS_d0 18/ARDS_d7 13，共 144）"),
    "_pcd_score_summary.json": ("analysis_output", "S5 PCD 评分摘要 JSON"),
}
SHARED_NUMBERS = {  # 一个文件登记在多个 S 编号下（多对一，方向合法）
    "Table_S22_S23_Drug_Database.csv": "S22/S23",
    "Table_S22_S23_Integrated_Drug_Ranking.csv": "S22/S23",
    "Table_S22_S23_Virtual_Drug_Screening_Report.txt": "S22/S23",
    "Table_S29_S31_scATAC_SCENIC_Analysis_Report.txt": "S29/S31",
    "Table_S32_S33_Methodology_Note.txt": "S32/S33",
    "Table_S38_S39_Geneformer_Report.txt": "S38/S39",
}


def s_number_of(fname):
    import re
    m = re.match(r"[Tt]able_([Ss](\d+))([a-z]?)[_.]", fname)
    if not m:
        return None
    return f"S{m.group(2)}"


def sha256_of(p, block=1024 * 1024):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while True:
            b = f.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def count_data_rows(p):
    try:
        enc, _ = _open_info(p)
        with open(p, "r", encoding=enc, newline="") as f:
            return max(0, sum(1 for _ in csv.reader(f)) - 1)
    except Exception:
        return None


def build_manifest():
    out_csv = T / "RESULTS_MANIFEST_v1.0.csv"
    header = ["s_number", "filename", "location", "class", "status", "gene_universe",
              "gene_set_version", "score_version", "n_data_rows", "size_bytes", "sha256",
              "analysis", "notes"]
    rows = []

    def add(fname, status, folder="", note_override=None):
        p = ROOT / folder / fname if folder else (ROOT / fname)
        sn = SHARED_NUMBERS.get(fname) or (s_number_of(fname) or "")
        if fname in OTHER_META:
            cls, desc = OTHER_META[fname]
        elif sn and sn in S_META:
            cls = "supplementary_table" if fname.lower().endswith(".csv") else "supplementary_report"
            desc = S_META[sn]
        elif sn:
            cls = "supplementary_table" if fname.lower().endswith(".csv") else "supplementary_report"
            desc = S_META.get(sn, "")
        else:
            cls, desc = "analysis_output", ""
        g, s, u = stamp_config_for(fname)
        if fname in STAMP_EXCLUDE or status == "archived":
            g, s, u = "NA", "NA", "none"
        if not fname.lower().endswith(".csv"):
            if fname in TXT_FOOTER:
                g, s = TXT_FOOTER[fname]
            else:
                g, s = "NA", "NA"
        note = note_override if note_override is not None else \
            ("历史备份/被取代版本，不用于正文（archive/）" if status == "archived" else "")
        rows.append([sn, fname, folder, cls, status, u, g, s,
                     count_data_rows(p) if p.exists() else "",
                     p.stat().st_size if p.exists() else "",
                     sha256_of(p) if p.exists() else "",
                     desc, note])

    for loc, d in [("04_AUDIT_GOVERNANCE", GOV),
                   ("02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV", TAB),
                   ("02_SUPPLEMENTARY_TABLES/Supplementary_Notes", NOTES)]:
        for p in sorted(d.iterdir()):
            if p.is_file() and not (d == GOV and p.name == "RESULTS_MANIFEST_v1.0.csv"):
                add(p.name, "canonical", folder=loc)
    for p in sorted(ARCH.iterdir()):
        if p.is_file():
            add(p.name, "archived", folder="archive")
    # 数据存放层输入（00_RAW_DATA，显式登记）
    for fn, loc in RAW_REGISTRY.items():
        if (ROOT / loc / fn).exists():
            add(fn, "canonical", folder=loc)
    # 登记在册但本地不随包的文件（留痕；与 lint MISSING_OK 白名单同源）
    for fn, (st, note) in ABSENT_REGISTRY.items():
        add(fn, st, folder="", note_override=note)

    write_csv_rows(out_csv, [header] + rows, "utf-8-sig", "\r\n")
    L(f"  [manifest] {out_csv.name}: {len(rows)} 行（canonical {sum(1 for r in rows if r[3]=='canonical')} + archived {sum(1 for r in rows if r[3]=='archived')}）")
    return rows


# --------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------
def main():
    L("== D4/N22 治理构建 开始 ==")
    ARCH.mkdir(exist_ok=True)

    L("[1] 归档备份/被取代文件")
    for f, farch in ARCHIVE_MOVES.items():
        src, dst = d4_path(f), ARCH / farch
        # S18 原版归档仅在对应 _CORRECTED 源仍在时执行（防止重跑时把已晋升的正式版再次归档）
        corr = d4_path(f.replace(".csv", "_CORRECTED.csv"))
        guard = (corr.name == f) or corr.exists()
        if src.exists() and guard:
            shutil.move(str(src), str(dst))
            L(f"  [归档] {f} -> archive/{farch}")
        elif dst.exists():
            L(f"  [跳过-已归档] {f}")
        else:
            L(f"  [警告] 归档目标不存在: {f}")

    L("[2] 编号去重重命名")
    for old, new in RENAMES.items():
        src, dst = d4_path(old), d4_path(new)
        if src.exists():
            shutil.move(str(src), str(dst))
            L(f"  [重命名] {old} -> {new}")
        elif dst.exists():
            L(f"  [跳过-已应用] {old} -> {new}")
        else:
            L(f"  [警告] 重命名源与目标均不存在: {old}")

    L("[3] 历史列名/取值/文本统一（数值零改动）")
    for f, rn in COL_RENAMES.items():
        apply_col_renames(f, rn)
    for f, vm in VALUE_RENAMES.items():
        apply_value_renames(f, vm)
    for f, pairs in TEXT_REPLACEMENTS.items():
        apply_text_replacements(f, pairs)

    L("[4] canonical 基因 manifest")
    build_gene_manifest()

    L("[5] 版本戳")
    stamp_all()

    L("[6] 冻结 RESULTS_MANIFEST")
    rows = build_manifest()

    (ROOT / "03_LOGS").mkdir(exist_ok=True)
    (ROOT / "03_LOGS" / "D4_build_log.txt").write_text(
        "\n".join(LOG) + "\n", encoding="utf-8")
    L("== D4/N22 治理构建 完成（日志: 03_LOGS/D4_build_log.txt）==")


if __name__ == "__main__":
    main()
