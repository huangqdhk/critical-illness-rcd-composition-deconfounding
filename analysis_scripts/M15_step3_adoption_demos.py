# -*- coding: utf-8 -*-
"""
M15_step3_adoption_demos.py — M15 步骤 3：采纳演示（2 个第三方数据集一键复算）
========================================================================
依次运行：
  - M15_tool/demos/demo_GSE157103.py（COVID-19 全血 RNA-seq，Overmyer 2020）
  - M15_tool/demos/demo_GSE66099.py（儿童 SIRS/脓毒症/脓毒性休克全血，GPSSSI，GPL570）
两个数据集均非本项目既有队列，入列前经 GEO 官方记录核验（红线 #1）：
03_LOGS/M15_demo_GEO_verification_20260827.md。
输出：_intermediate/M15_demo_*_per_sample/contrasts/r2/proportions.csv；
      03_LOGS/M15_demo_*_report.md。
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
sys.path.insert(0, os.path.join(ROOT, "M15_tool"))
sys.path.insert(0, os.path.join(ROOT, "M15_tool", "demos"))

import demo_GSE157103 as d1  # noqa: E402
import demo_GSE66099 as d2  # noqa: E402

print("===== M15 采纳演示 1/2：GSE157103 =====", flush=True)
d1.main()
print("\n===== M15 采纳演示 2/2：GSE66099 =====", flush=True)
d2.main()
print("\n===== M15 采纳演示全部完成 =====", flush=True)
