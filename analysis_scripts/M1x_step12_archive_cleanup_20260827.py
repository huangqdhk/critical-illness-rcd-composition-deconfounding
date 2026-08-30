# -*- coding: utf-8 -*-
"""
M1x_step12_archive_cleanup_20260827.py — 2026-08-27 目录清理批次
================================================================
背景：M10–M15 优化批次全部收口（文稿 v4.7、lint 74/0）后，按作者要求清理
03_LOGS / 04_AUDIT_GOVERNANCE / M15_tool / 根目录散在脚本：

  1) 调试/一次性脚本与日志 → archive/ 子目录（沿用 R1_*/M1_check_debug/M3_debug 命名惯例）；
  2) archive/ 内容按 lint 治理要求逐文件登记进 RESULTS_MANIFEST_v2.0.csv（status=archived）；
  3) 冗余控制台拷贝（*_run.txt：与同次运行 *_log.txt 内容一致、仅多 DONE 行）删除；
  4) 可再生缓存 __pycache__/ 删除；
  5) 现役脚本（M1–M15/P0–P2/基础库）保持根目录平铺不动（mdi_lib 同目录 import 约束，
     见 README_打包说明）——本轮清理不移动任何现役脚本。

安全设计：移动前核对源存在、目标不存在；manifest 全量读改写（utf-8-sig，列序不变）；
只读源文件内容，不修改任何数据文件。
输出：03_LOGS/M1x_archive_cleanup_20260827_log.txt
"""
import csv
import hashlib
import shutil
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
LOGS = ROOT / "03_LOGS"
GOV = ROOT / "04_AUDIT_GOVERNANCE"
ARCH = ROOT / "archive"
MANIFEST = GOV / "RESULTS_MANIFEST_v2.0.csv"
LOG_OUT = LOGS / "M1x_archive_cleanup_20260827_log.txt"

BATCH = "2026-08-27 清理归档（M1x_step12 批次）"

# 目标归档子目录（仿照既有命名惯例：R1_*、M1_check_debug、M3_debug）
NEW_SUBDIRS = {
    "M2_check_debug": "M2 轮调试/探查脚本与输出（_m2_*），2026-08-27 清理批次入档",
    "ops_scripts_2026": "跨轮次一次性运维脚本与日志（R 环境安装/图清单修复/关闭口脚本/pip 日志等），2026-08-27 清理批次入档",
}

MOVES = {
    # ---- 04_AUDIT_GOVERNANCE ----
    "04_AUDIT_GOVERNANCE/rawdata_cleanup_20260817.py": "ops_scripts_2026",
    "04_AUDIT_GOVERNANCE/rawdata_cleanup_20260817.log": "ops_scripts_2026",
    "04_AUDIT_GOVERNANCE/README_5.1_Scrna_Basic.md": "R1_docs",
    # ---- 03_LOGS 一次性运维脚本/日志 ----
    "03_LOGS/_fix_figure_manifest_20260824.py": "ops_scripts_2026",
    "03_LOGS/_m1x_closeout_20260827.py": "ops_scripts_2026",
    "03_LOGS/_run_R_install.py": "ops_scripts_2026",
    "03_LOGS/_R_install_annot.R": "ops_scripts_2026",
    "03_LOGS/_R_install_annot.log": "ops_scripts_2026",
    "03_LOGS/_m10d_console.txt": "ops_scripts_2026",
    "03_LOGS/_pip_m2.log": "ops_scripts_2026",
    "03_LOGS/_pip_openpyxl.log": "ops_scripts_2026",
    "03_LOGS/_pip_restore.log": "ops_scripts_2026",
    "03_LOGS/_pip_venvm2.log": "ops_scripts_2026",
    # ---- 03_LOGS M1/M2 轮调试残留 ----
    "03_LOGS/_m1_inspect_annotation.txt": "M1_check_debug",
    "03_LOGS/_m1_overlap_out.txt": "M1_check_debug",
    "03_LOGS/_fig4g_check.txt": "M2_check_debug",
    "03_LOGS/GSE185263_sample_group_assignment.csv": "M2_check_debug",
    # ---- 03_LOGS 旧轮次日志 ----
    "03_LOGS/wgcna_run.log": "R1_algo_scripts",
}

# 删除项（逐项留痕）
DELETES = [
    "03_LOGS/P0_donor_pseudobulk_run.txt",
    "03_LOGS/P1_loro_run.txt",
    "03_LOGS/P2_gate2_run.txt",
    "03_LOGS/_lint_console_tmp.txt",
]

# 删除的可再生缓存目录（__pycache__，运行即重建）
CACHE_DIRS = [
    "__pycache__",
    "M15_tool/demos/__pycache__",
    "M15_tool/mitoxdi/__pycache__",
]

# 删除 *_run.txt 的依据（与 *_log.txt 逐字节一致、仅多末尾 DONE 行；均为同次运行的控制台拷贝）
DELETE_REASON = {
    "03_LOGS/P0_donor_pseudobulk_run.txt": "P0_donor_pseudobulk_log.txt 的同次运行控制台拷贝（仅多 DONE 行）",
    "03_LOGS/P1_loro_run.txt": "P1_loro_log.txt 的同次运行控制台拷贝（仅多 DONE 行）",
    "03_LOGS/P2_gate2_run.txt": "P2_gate2_log.txt 的同次运行控制台拷贝（仅多 DONE 行）",
    "03_LOGS/_lint_console_tmp.txt": "本次清理前的 lint 临时控制台捕获",
}

# 保留在 03_LOGS 的一次性运维脚本（有现役引用/治理文档路径引用，不动）
#   _run_R_export.py / _R_export_maps.log —— M2_CAP_COPD_export_probe_maps.R 注释声明的运行器
#   _R_sessionInfo_capture_20260823.R —— 未完成项标注 §五 记录的 R sessionInfo 捕获脚本

# 特殊说明（覆盖通用备注）
SPECIAL_NOTES = {
    "README_5.1_Scrna_Basic.md": f"{BATCH}：第 1 轮（R1）5.1 scRNA 基础管线说明（旧 38 基因签名口径），已被现行 M1–M4 流程取代，移入 archive/R1_docs。",
    "wgcna_run.log": f"{BATCH}：第 1 轮（R1）WGCNA 运行日志，移入 archive/R1_algo_scripts。",
    "rawdata_cleanup_20260817.py": f"{BATCH}：一次性原始数据清理脚本（2026-08-17 已执行完毕），移入 archive/ops_scripts_2026。",
    "rawdata_cleanup_20260817.log": f"{BATCH}：一次性原始数据清理记录（2026-08-17 已执行完毕），移入 archive/ops_scripts_2026。",
    "_R_install_annot.R": f"{BATCH}：一次性 R 环境注释包安装脚本（已完成使命），移入 archive/ops_scripts_2026。",
    "GSE185263_sample_group_assignment.csv": f"{BATCH}：M2 轮调试期的样本分组分配表（canonical 分组定义见 04_AUDIT_GOVERNANCE/GSE185263_groups.csv），移入 archive/M2_check_debug。",
}


def cls_for(name: str) -> str:
    low = name.lower()
    if low.endswith((".py", ".r")):
        return "audit_output"
    if low.endswith((".log", ".txt")):
        return "audit_report"
    if low.endswith(".md"):
        return "report"
    if low.endswith((".csv", ".json")):
        return "analysis_output"
    return "audit_output"


def sha256_of(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    log = []
    def say(msg: str):
        log.append(msg)
        print(msg)

    say(f"== M1x_step12 清理批次 {date.today()} 启动")
    say(f"   ROOT={ROOT}")

    # 0. 准备归档子目录
    for d, desc in NEW_SUBDIRS.items():
        (ARCH / d).mkdir(exist_ok=True)
        say(f"[mkdir] archive/{d} 就绪（{desc}）")

    # 1. 构造完整移动清单（含 _m2_* 通配）
    moves: dict[str, str] = dict(MOVES)
    for p in sorted(LOGS.glob("_m2_*.py")) + sorted(LOGS.glob("_m2_*.txt")):
        moves[f"03_LOGS/{p.name}"] = "M2_check_debug"
    moves = {k: moves[k] for k in sorted(moves)}

    # 2. 执行移动
    moved_files = []
    for rel_src, sub in moves.items():
        src = ROOT / rel_src
        dst_dir = ARCH / sub
        dst = dst_dir / src.name
        if not src.exists():
            say(f"[SKIP-缺失] {rel_src} 不存在，跳过")
            continue
        if dst.exists():
            say(f"[SKIP-冲突] {rel_src} → {sub}/{src.name} 目标已存在，跳过")
            continue
        shutil.move(str(src), str(dst))
        moved_files.append((src.name, sub))
        say(f"[move] {rel_src} → archive/{sub}/{src.name}")

    # 3. 删除冗余拷贝与缓存
    for rel in DELETES:
        p = ROOT / rel
        if p.exists():
            p.unlink()
            say(f"[delete] {rel}（{DELETE_REASON.get(rel, '')}）")
        else:
            say(f"[delete-跳过] {rel} 不存在")
    for rel in CACHE_DIRS:
        d = ROOT / rel
        if d.is_dir():
            shutil.rmtree(d)
            say(f"[delete-cache] {rel} 已删除（可再生缓存）")
        else:
            say(f"[delete-cache-跳过] {rel} 不存在")

    # 4. 更新 RESULTS_MANIFEST_v2.0.csv（逐文件登记为 archived）
    with open(MANIFEST, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        rows = list(reader)

    by_name = {}
    for i, r in enumerate(rows):
        by_name.setdefault(r["filename"], []).append(i)

    added = 0
    updated = 0
    for name, sub in moved_files:
        note = SPECIAL_NOTES.get(name, f"{BATCH}：调试/一次性脚本或日志，不再参与当前分析；不影响任何冻结输出，移入 archive/{sub}。")
        loc = f"archive/{sub}"
        if name in by_name:
            for i in by_name[name]:
                r = rows[i]
                r["location"] = loc
                r["status"] = "archived"
                r["gene_universe"] = "none"
                r["gene_set_version"] = ""
                r["score_version"] = ""
                r["size_bytes"] = str((ARCH / sub / name).stat().st_size)
                r["notes"] = note
                updated += 1
                say(f"[manifest-更新] {name}: location={loc}, status=archived")
        else:
            size = (ARCH / sub / name).stat().st_size
            rows.append({
                "s_number": "", "filename": name, "location": loc,
                "class": cls_for(name), "status": "archived",
                "gene_universe": "none", "gene_set_version": "",
                "score_version": "", "n_data_rows": "", "size_bytes": str(size),
                "sha256": sha256_of(ARCH / sub / name), "analysis": "",
                "notes": note, "class_": "",
            })
            added += 1
            say(f"[manifest-新增] {name}: location={loc}, status=archived, size={size}")

    with open(MANIFEST, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    say(f"[manifest] 总行数 {len(rows)}（新增 {added}、更新 {updated}）")

    say(f"== 清理批次完成：移动 {len(moved_files)} 个文件、删除 {len(DELETES)} 个冗余文件、清缓存目录 {len(CACHE_DIRS)} 个")
    say("== 下一步：运行 python lint_package.py 复核全绿")

    with open(LOG_OUT, "w", encoding="utf-8") as f:
        f.write("\n".join(log) + "\n")
    return 0


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    raise SystemExit(main())
