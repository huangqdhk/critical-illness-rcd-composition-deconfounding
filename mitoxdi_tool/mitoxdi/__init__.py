# -*- coding: utf-8 -*-
"""
mitoxdi — M15 开放工具包：两臂分解 + mdi_v1.0 评分 + 组成预测值/组成残差 MDI
================================================================================
工具核心冻结口径（与稿件 Methods §4.18/§4.22/§4.25 一致）：
  - 两臂分解：80 基因 manifest 的 arm 列（upstream_collapse 30 / execution_induction 33），
    唯一基因集来源 = Mitoxyperilysis_Gene_Manifest_v1.0.csv（含别名表，HGNC 归一）。
  - mdi_v1.0 评分：UCS = z(上游臂均值)，EIS = z(执行臂均值)，MDI = EIS − UCS；
    z 在队列内用全部样本（不分组）标准化，防止组内泄漏；敏感性 = nomt（去 6 个 MT-*）
    与 ssgsea（running-sum β=0.25）。
  - 组成预测：X̂_comp(s) = Σ_c π_c(s)·μ_c（log2 空间）；反卷积 NNLS 主口径（OLS-CLS 交叉）；
    参考谱 Monaco GSE107011 29 型（含粒细胞，log2(TPM+1)）与 ABIS RNAseq 17 型 /
    ABIS Micro 11 型（log2(+1)）；对合成谱按 mdi_v1.0 同口径计算 MDI_comp；
    组成残差 MDI_resid = 观测 MDI − (b0 + b1·MDI_comp)（队列内 OLS）。
  - 判定门语言：见 M15_tool/GATE.md（M10 检验 A 主检验的通用化表述）。
版本：mitoxdi v1.0.0（2026-08-27）；score_version = mdi_v1.0；
     gene_set_version = Mitoxy-80_v1.0（与项目冻结口径一致）。
"""
import os
import sys

__version__ = "1.0.0"
SCORE_VERSION = "mdi_v1.0"
GENE_SET_VERSION = "Mitoxy-80_v1.0"

# 项目根目录（可被环境变量覆盖；发布版中指向本地数据根）
ROOT = os.environ.get("MITOXDI_PROJECT_ROOT", r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS")

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

from . import mdi_core  # noqa: E402,F401
from . import composition  # noqa: E402,F401
from . import report  # noqa: E402,F401

__all__ = ["mdi_core", "composition", "report",
           "SCORE_VERSION", "GENE_SET_VERSION", "__version__"]
