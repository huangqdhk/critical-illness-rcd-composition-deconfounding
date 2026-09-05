# -*- coding: utf-8 -*-
"""跨机通用引导器（2026-09-04 建，解决导师机绝对路径报错问题）。

背景：本包约 40 个历史脚本硬编码 ROOT = r"E:\\SCI\\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
（另有 3 个 M3 脚本指向 E:\\SCI\\proj）。整包拷到任何其他机器/盘符后直跑这些脚本必然报
FileNotFoundError。本引导器沿用 P13_rerun_qc_local.py（2026-09-03）确立的模式：
读入目标脚本源码，在内存中把绝对路径替换为"本脚本所在目录"，再 exec 执行——
目标脚本本体不动，保持审计链。

用法（任意机器、任意 cwd）：
    python run_anywhere.py <脚本名.py> [脚本参数...]
例如：
    python run_anywhere.py P0_manuscript_source_check.py
    python run_anywhere.py M2_step2_replication_meta.py

说明：
- 现役 QC 入口（P0_manuscript_source_check*.py、lint_package.py）已改为 __file__ 相对
  路径，可直接运行，无需本引导器；本引导器主要用于其余冻结管线脚本。
- E:\\SCI\\proj 指向包外的兄弟目录；若本机不存在则替换为本包同级的 proj/（存在时），
  否则保留原串（脚本自行报错，提示缺该目录）。
- 替换会在日志中逐条打印，便于追溯。
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.join(os.path.dirname(HERE), "proj")

# 待替换的绝对路径 → 本机实际路径（仅当目标存在时才替换 PROJ）
REPLACEMENTS = [
    (r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS", HERE),
    ("E:/SCI/SCI论文1黄裕荣_Mitoxyperilysis_ARDS", HERE),
    (r"E:\SCI\proj", PROJ if os.path.isdir(PROJ) else None),
]

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    script = sys.argv[1]
    path = os.path.join(HERE, script) if not os.path.isabs(script) else script
    if not os.path.isfile(path):
        print(f"[run_anywhere] 找不到脚本: {path}")
        sys.exit(1)
    src = open(path, encoding="utf-8", errors="replace").read()
    for old, new in REPLACEMENTS:
        if new is None:
            continue
        n = src.count(old)
        if n:
            src = src.replace(old, new)
            print(f"[run_anywhere] {script}: 替换 {old} -> {new}（{n} 处）")
    sys.argv = [script] + sys.argv[2:]
    g = {"__name__": "__main__", "__file__": path}
    exec(compile(src, path, "exec"), g)

if __name__ == "__main__":
    main()
