# -*- coding: utf-8 -*-
"""
M15_step2_tool_regression.py — M15 步骤 2：工具内部回归测试（11 队列）
========================================================================
运行 M15_tool/tests/test_regression_11cohorts.py：
工具核心（两臂分解 + mdi_v1.0 评分 + 组成预测值/组成残差 MDI）对
11 个队列的冻结逐样本表（M2 6 队列 + M10A 组成表 + M11M12 5 队列）逐样本核对。
输出：_intermediate/M15_regression_summary.csv / M15_regression_cohort_summary.csv；
      03_LOGS/M15_regression_log.txt。
"""
import subprocess
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

TEST = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\M15_tool\tests\test_regression_11cohorts.py"
PY = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\.venv_pyaging\Scripts\python.exe"

r = subprocess.run([PY, TEST])
sys.exit(r.returncode)
