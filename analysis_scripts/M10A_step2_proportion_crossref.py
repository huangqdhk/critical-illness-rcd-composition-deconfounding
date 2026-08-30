# -*- coding: utf-8 -*-
"""Recompute Monaco-vs-ABIS proportion cross-reference from the proportions CSV."""
import sys
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

props = pd.read_csv(r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate\M10A_proportions.csv")
props["sample"] = props["sample"].astype(str)

TYPE_MAP = {"Monocytes C": "C_mono", "NK": "NK", "T CD8 Memory": "CD8_CM",
            "T CD8 Naive": "CD8_naive", "T CD4 Naive": "CD4_naive", "T CD4 Memory": "CD4_TE",
            "B Memory": "B_SM", "B Naive": "B_naive", "Neutrophils LD": "Neutrophils",
            "Basophils LD": "Basophils", "pDCs": "pDC", "mDCs": "mDC",
            "Plasmablasts": "Plasmablasts", "MAIT": "MAIT", "T gd Vd2": "VD2+",
            "T gd non-Vd2": "VD2-"}

cohort_map = {"GSE148871_blood": "GSE148871"}
rows = []
for c in ["GSE185263", "GSE32707", "GSE212865", "GSE188309", "GSE310929", "GSE148871_blood"]:
    pm = props[(props.cohort == c) & (props.reference == "Monaco") & (props.method == "NNLS")].set_index("sample")
    pa = props[(props.cohort == c) & (props.reference == "ABIS") & (props.method == "NNLS")].set_index("sample")
    cmn = sorted(set(pm.index) & set(pa.index))
    rs = []
    for a_t, m_t in TYPE_MAP.items():
        if a_t in pa.columns and m_t in pm.columns:
            r = np.corrcoef(pa.loc[cmn, a_t].astype(float), pm.loc[cmn, m_t].astype(float))[0, 1]
            if np.isfinite(r):
                rs.append(r)
    rows.append(dict(cohort=cohort_map.get(c, c), n_samples=len(cmn), n_types=len(rs),
                     median_r=float(np.median(rs)) if rs else np.nan,
                     mean_r=float(np.mean(rs)) if rs else np.nan,
                     min_r=float(np.min(rs)) if rs else np.nan,
                     max_r=float(np.max(rs)) if rs else np.nan))
    print(c, len(cmn), len(rs), round(float(np.median(rs)), 3) if rs else None)
out = pd.DataFrame(rows)
out.to_csv(r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate\M10A_proportion_cross_reference.csv", index=False)
out.to_csv(r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\02_SUPPLEMENTARY_TABLES\SUPPLEMENTARY_Tables_CSV\Table_S65c_M10A_Proportion_CrossReference.csv", index=False)
print("saved")
