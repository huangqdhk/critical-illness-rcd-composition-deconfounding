# -*- coding: utf-8 -*-
"""
M1 Step 9: Register M1 outputs into governance records.

1. Append new tables to RESULTS_MANIFEST_v1.0.csv (SHA256 frozen).
2. Create M1_spatial_sample_manifest.csv (23 Visium sections + 116 CosMx FOVs).
3. Update README_附表索引.md with S55a-v entries.
"""
import os, hashlib
import pandas as pd

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
TAB = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES", "SUPPLEMENTARY_Tables_CSV")
FIG = os.path.join(ROOT, "01_FIGURE_DATA_CSV", "Main")
GOV = os.path.join(ROOT, "04_AUDIT_GOVERNANCE")

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

# ---- list of new M1 output files (source tables only; no figure-panel copies, 方案A) ----
new_files = []
for fn in sorted(os.listdir(TAB)):
    if fn.startswith("Table_S55") and fn.endswith(".csv"):
        new_files.append((os.path.join(TAB, fn), "02_SUPPLEMENTARY_TABLES"))

# arm-gene tagging metadata per table
TAG = {
    "Table_S55a_M1_Visium_Section_Spatial_Stats.csv": "Mitoxy-80_v1.0",
    "Table_S55c_M1_Visium_Domain_Stats.csv": "Mitoxy-80_v1.0",
    "Table_S55d_M1_Visium_SVG_Detail.csv": "all_detected",
    "Table_S55h_M1_Visium_KeyGene_Moran_Perm.csv": "Mitoxy-80_v1.0",
    "Table_S55s_M1_Visium_TNFSF13B_TFRC_Neighborhood.csv": "TNFSF13B/TFRC",
}

existing = pd.read_csv(os.path.join(GOV, "RESULTS_MANIFEST_v1.0.csv"), encoding="utf-8-sig")
existing_files = set(existing["filename"])

rows = []
for path, loc in new_files:
    fn = os.path.basename(path)
    if fn in existing_files:
        continue
    n_rows = sum(1 for _ in open(path, "rb"))
    rows.append({
        "s_number": "S55",
        "filename": fn,
        "location": loc,
        "class": "analysis_output",
        "status": "canonical",
        "gene_universe": ("mitoxy_80" if fn not in TAG or TAG[fn] == "Mitoxy-80_v1.0" else TAG[fn]),
        "gene_set_version": TAG.get(fn, "Mitoxy-80_v1.0"),
        "score_version": "score_genes_v1_spatial",
        "n_data_rows": n_rows - 1,
        "size_bytes": os.path.getsize(path),
        "sha256": sha256(path),
        "analysis": "M1 空间地理学验证（GSE271370 Visium 23切片 + GSE253474 CosMx 98,850细胞；判定 R3）",
        "notes": "gene_set_version 对 SVG 全基因表为 all_detected",
    })
new_df = pd.DataFrame(rows)
merged = pd.concat([existing, new_df], ignore_index=True)
merged.to_csv(os.path.join(GOV, "RESULTS_MANIFEST_v1.0.csv"), index=False, encoding="utf-8-sig")
print(f"registered {len(rows)} new files; manifest total rows = {len(merged)}")

# ---- M1 spatial sample manifest ----
vis = [
    ("GSM8375640", "L2P", "ProliferativeDAD"), ("GSM8375641", "L19P", "ProliferativeDAD"),
    ("GSM8375642", "L11P", "ProliferativeDAD"), ("GSM8375643", "CONTROL2", "Control"),
    ("GSM8387358", "HRC5", "ProliferativeDAD"), ("GSM8387359", "HRC6", "ProliferativeDAD"),
    ("GSM8387360", "HRC8", "ProliferativeDAD"), ("GSM8387361", "HRC10", "ProliferativeDAD"),
    ("GSM8387362", "HRC11", "ProliferativeDAD"), ("GSM8387363", "HRC12", "ProliferativeDAD"),
    ("GSM8387364", "HRC13", "ProliferativeDAD"), ("GSM8387365", "HRC16", "ProliferativeDAD"),
    ("GSM8387366", "HRC17", "ProliferativeDAD"), ("GSM8387367", "L5P", "AcuteDAD"),
    ("GSM8387368", "L14P", "AcuteDAD"), ("GSM8387369", "L24P", "AcuteDAD"),
    ("GSM8387370", "L12P", "AcuteDAD"), ("GSM8387371", "HRC2", "AcuteDAD"),
    ("GSM8387372", "HRC4", "AcuteDAD"), ("GSM8387373", "HRC18", "AcuteDAD"),
    ("GSM8387374", "L3C", "Control"), ("GSM8387375", "L14C", "Control"),
    ("GSM8387376", "L2C", "Control"),
]
cosmx_raw = os.path.join(ROOT, "00_RAW_DATA", "GSE253474_Lung_ARDS_CosMx")
fov_rows = []
for fn in sorted(os.listdir(cosmx_raw)):
    if fn.endswith("_exprMat_file.csv") and not fn.endswith(".gz"):
        gsm = fn.split("_")[0]
        parts = fn.replace(gsm + "_", "").replace("_exprMat_file.csv", "").split("_")
        # TMA_A_8_Case_6_FOV_1 -> ["TMA","A","8","Case","6","FOV","1"]
        slide = "_".join(parts[:3])
        case = parts[4]
        fov = parts[6]
        fov_rows.append(dict(gsm=gsm, slide=slide, case=case, fov=int(fov)))
fov_df = pd.DataFrame(fov_rows).sort_values(["slide", "case", "fov"])
print("cosmx FOVs found:", len(fov_df))

man_rows = []
for gsm, sec, cond in vis:
    man_rows.append(dict(dataset="GSE271370", gsm=gsm, unit="section",
                         sample_id=sec, condition=cond, tissue="lung",
                         platform="Visium_FFPE", notes="fatal COVID-19, GEO disease state verified 2026-08-17"))
for _, r in fov_df.iterrows():
    man_rows.append(dict(dataset="GSE253474", gsm=r.gsm, unit="FOV",
                         sample_id=f"{r.slide}_Case{r.case}_FOV{r.fov}",
                         condition="ARDS", tissue="lung",
                         platform="CosMx_SMI", notes="ARDS death cases; slide/case/fov from GSM titles"))
man_df = pd.DataFrame(man_rows)
man_df.to_csv(os.path.join(GOV, "M1_spatial_sample_manifest.csv"), index=False)
print("M1_spatial_sample_manifest.csv:", len(man_df), "rows")

# ---- update README_附表索引.md ----
readme_path = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES", "README_附表索引.md")
txt = open(readme_path, encoding="utf-8").read()
if "Table S55a" not in txt:
    add = """
## M1 空间地理学验证（2026-08-17 新增，Table S55a–S55v）

| 表 | 内容 | 数据源 |
|---|---|---|
| Table S55a | Visium 23 切片逐切片空间统计（双臂 Moran's I / 双变量 I / 域统计，置换 p） | GSE271370 |
| Table S55b | Visium 条件级汇总（对照/急性DAD/增殖性DAD，Kruskal-Wallis） | GSE271370 |
| Table S55c | Visium 空间域（KMeans k=5）逐域双臂均值 | GSE271370 |
| Table S55d | SVG 明细（逐切片逐基因 Moran's I + 随机化矩 z） | GSE271370 |
| Table S55e | SVG Meta（逐基因 Stouffer z，80 基因标注 arm） | GSE271370 |
| Table S55f | SVG 富集（各臂 vs 背景，Mann-Whitney + Stouffer） | GSE271370 |
| Table S55g | SVG Top10% 超几何富集（80 基因各臂） | GSE271370 |
| Table S55h | 80 基因逐切片 Moran's I 置换 p | GSE271370 |
| Table S55i | CosMx 基因检出率（归一化矩阵全细胞有值，披露） | GSE253474 |
| Table S55j | CosMx 逐 TMA 炎症小体模块/单基因 Moran's I（999 置换） | GSE253474 |
| Table S55k | CosMx 模块评分 vs 邻域髓系比例（Spearman + 双变量 I） | GSE253474 |
| Table S55l | CosMx 模块评分按细胞类型（均值/高分比例） | GSE253474 |
| Table S55m | CosMx 病毒区生态位（POS/ADJ/NEG）模块比较 | GSE253474 |
| Table S55n | CosMx TNFSF13B 细胞来源（恒量表达，披露为不可用） | GSE253474 |
| Table S55o | CosMx TNFSF13B 邻域富集（恒量表达致无效，已弃用） | GSE253474 |
| Table S55p | Visium NNLS（scRNA 8 类签名）细胞组成按条件 | GSE271370+GSE145926/GSE158055 |
| Table S55q | NNLS 比例 vs marker 评分相关性（弱，如实披露） | 同上 |
| Table S55r | Visium 条件两两事后检验（BH） | GSE271370 |
| Table S55s | TNFSF13B-TFRC 邻域富集（999 置换，逐切片+Stouffer） | GSE271370 |
| Table S55t | TNFSF13B/TFRC 与 NNLS 细胞类型的 Spearman 归属 | GSE271370 |
| Table S55u | Visium NNLS（CosMx 27 类签名，验证参考）细胞组成按条件 | GSE271370+GSE253474 |
| Table S55v | 计数残差化双变量 Moran's I 敏感性 | GSE271370 |

判定：**R3（阴性）**——双臂空间共定位（双变量 I 合并 z=+6.33，p=2.5e-10），无区室分离；
执行臂空间聚集更强（21/23 切片）；SVG 富集执行臂最强（z=+9.36）。
详见 `_intermediate/M1_judgment_report.md` 与 `N4_资源核实报告.md`。
"""
    open(readme_path, "a", encoding="utf-8").write(add)
    print("README_附表索引.md updated")
else:
    print("README already updated")
print("step 9 done")
