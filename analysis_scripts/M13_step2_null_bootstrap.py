# -*- coding: utf-8 -*-
"""M13 addendum: bootstrap empirical p for beta_myeloid from the saved null model."""
import sys
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import statsmodels.api as sm

n = pd.read_csv(r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate\M13_null_model.csv")
X = sm.add_constant(n[["execution_share_ext", "myeloid_share"]].astype(float)).values
rng = np.random.default_rng(0)
boots = []
for b in range(2000):
    idx = rng.integers(len(n), size=len(n))
    f = sm.OLS(n.iloc[idx]["g_185"].values, X[idx]).fit()
    boots.append((f.params[1], f.params[2]))
boots = np.array(boots)
obs_exec, obs_my = 0.719, 2.309
print("exec empirical p:", (1 + (boots[:, 0] >= obs_exec).sum()) / (len(boots) + 1))
print("myeloid empirical p:", (1 + (boots[:, 1] >= obs_my).sum()) / (len(boots) + 1))
print("null beta_myeloid mean/sd:", round(boots[:, 1].mean(), 3), round(boots[:, 1].std(), 3))
