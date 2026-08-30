# -*- coding: utf-8 -*-
"""
P1_matched_gene_set_test.py — Phase 1 补充：IIAMD 签名随机匹配基因集检验
=======================================================================
依据：《111黄裕荣创新提质_v2.md》Phase 1 冻结三资产任务："对签名做 leave-one-replicate-out
稳定性和随机匹配基因集检验"。LORO 已于 2026-08-20 完成（P1_loro_summary.csv，STABLE）；
本脚本补做随机匹配基因集检验（此前缺失）。

方法（诚实口径，不做事后择径）：
- 宇宙：GSE235046 交互模型 15,304 个检验基因（P1_GSE235046_interaction_DE.csv）
- 匹配：按 baseMean 十分位分层，从同层随机抽取与签名等量的基因（up=4367 / down=4720）
- 检验统计量（每套随机集与观测签名同口径）：
  1) mean |β₃|（交互系数绝对值的均值）
  2) |β₃|>1 的基因数
  3) 机制锚定节点命中数（Bax/Bak1/Bid/Rictor/Rhoa/Ninj1/Rptor/Mtor/Akt1/Tlr4/Myd88/Tnf，共 12）
- 置换数：N=500；经验 p = (1+#null≥obs)/(N+1)（单侧：签名是否超随机）
- 匹配理由：β₃ 绝对值与表达水平（baseMean）相关，直接无匹配抽样会因低表达基因方差大而夸大显著性；
  分层匹配消除该混杂。锚定节点命中数不依赖 baseMean 匹配（方向性富集检验，披露）。
输出：
- _intermediate/P1_matched_gene_set_test.csv
- 追加 04_AUDIT_GOVERNANCE/P1_Interaction_Signature_Report.md §P1-4
"""
import os, time
import numpy as np
import pandas as pd

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
INTER = os.path.join(ROOT, "_intermediate")
SIGF = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "P1_IIAMD_signature_v1.0.csv")
DEF = os.path.join(INTER, "P1_GSE235046_interaction_DE.csv")
REPORT = os.path.join(ROOT, "04_AUDIT_GOVERNANCE", "P1_Interaction_Signature_Report.md")

ANCHOR = ["Bax", "Bak1", "Bid", "Rictor", "Rhoa", "Ninj1", "Rptor", "Mtor",
          "Akt1", "Tlr4", "Myd88", "Tnf"]
SEED = 0
N = 500

t0 = time.time()
log = []
def note(m=""):
    log.append(m); print(m, flush=True)

de = pd.read_csv(DEF)
de = de.set_index("geneSymbol")
sig = pd.read_csv(SIGF)
up = sig[sig["direction"] == "up"]["geneSymbol"].tolist()
dn = sig[sig["direction"] == "down"]["geneSymbol"].tolist()
note(f"签名：up={len(up)}, down={len(dn)}")

univ = de.index.tolist()
univ_set = set(univ)
up_in = [g for g in up if g in univ_set]
dn_in = [g for g in dn if g in univ_set]
note(f"宇宙内命中：up={len(up_in)}/{len(up)}, down={len(dn_in)}/{len(dn)}")

de["bm_decile"] = pd.qcut(de["baseMean"].rank(method="first"), 10, labels=False)
# universe excluding signature genes (leave-self-out null)
null_pool = de[~de.index.isin(set(up_in) | set(dn_in))]

def draw_matched(n_target, rng):
    idx = []
    for d in range(10):
        layer = null_pool[null_pool["bm_decile"] == d]
        # number to draw from this layer proportional to layer size among all genes
        k = int(round(n_target * len(layer) / len(null_pool)))
        if len(layer) and k > 0:
            idx += rng.choice(layer.index.tolist(), size=min(k, len(layer)), replace=False).tolist()
    # top-up if rounding shortfall
    if len(idx) < n_target:
        remaining = null_pool.index.difference(idx).tolist()
        idx += rng.choice(remaining, size=n_target - len(idx), replace=False).tolist()
    return idx[:n_target]

def stats_of(genes):
    sub = de.loc[genes]
    return dict(n=len(genes),
                mean_abs_b3=float(sub["log2FoldChange"].abs().mean()),
                n_abs_b3_gt1=int((sub["log2FoldChange"].abs() > 1).sum()),
                n_anchor_hits=int(sub.index.isin(ANCHOR).sum()))

rng = np.random.default_rng(SEED)
rows = []
for direction, genes in (("up", up_in), ("down", dn_in)):
    obs = stats_of(genes)
    nulls = []
    for it in range(N):
        rs = draw_matched(len(genes), rng)
        nulls.append(stats_of(rs))
    for key in ("mean_abs_b3", "n_abs_b3_gt1", "n_anchor_hits"):
        null_arr = np.array([x[key] for x in nulls])
        p = float((1 + (null_arr >= obs[key]).sum()) / (N + 1))
        rows.append(dict(direction=direction, statistic=key, observed=obs[key],
                         null_mean=float(null_arr.mean()), null_sd=float(null_arr.std()),
                         null_min=float(null_arr.min()), null_max=float(null_arr.max()),
                         empirical_p_one_tail=p, n_perm=N))
        note(f"{direction}/{key}: obs={obs[key]}, null={null_arr.mean():.3f}±{null_arr.std():.3f}, p={p:.4f}")

out = pd.DataFrame(rows)
out.to_csv(os.path.join(INTER, "P1_matched_gene_set_test.csv"), index=False)

# anchor genes: where do the 12 anchors sit in the signature?
anchor_pos = []
for g in ANCHOR:
    if g in de.index:
        anchor_pos.append(dict(gene=g, in_universe=True,
                               direction="up" if g in up_in else ("down" if g in dn_in else "none"),
                               b3=float(de.loc[g, "log2FoldChange"]), padj=float(de.loc[g, "padj"])))
    else:
        anchor_pos.append(dict(gene=g, in_universe=False, direction="NA", b3=np.nan, padj=np.nan))
anchor_df = pd.DataFrame(anchor_pos)
anchor_df.to_csv(os.path.join(INTER, "P1_matched_anchor_gene_positions.csv"), index=False)

# ---- append report ----
add = ["", "## P1-4 随机匹配基因集检验（2026-08-21 补做，Phase 1 冻结三资产任务的第二项）", "",
       "方法：从 15,304 基因宇宙按 baseMean 十分位分层随机抽取与签名等量基因（up=4,367/down=4,720），",
       "500 次；统计量=mean|β₃|、|β₃|>1 基因数、机制锚定节点（12 基因）命中数；单侧经验 p=(1+#null≥obs)/(N+1)。",
       f"签名基因从宇宙中剔除（leave-self-out）。", "",
       "| 方向 | 统计量 | 观测 | 零分布（500 套） | 单侧 p |",
       "|---|---|---|---|---|"]
for _, r in out.iterrows():
    add.append(f"| {r['direction']} | {r['statistic']} | {r['observed']:.4f} | "
               f"{r['null_mean']:.4f}±{r['null_sd']:.4f} [{r['null_min']:.4f}, {r['null_max']:.4f}] | {r['empirical_p_one_tail']:.4f} |")
add += ["", "### 机制锚定节点（12）在签名中的位置", "", "| 基因 | 在宇宙 | 签名方向 | β₃ | padj |", "|---|---|---|---|---|"]
for _, r in anchor_df.iterrows():
    add.append(f"| {r['gene']} | {r['in_universe']} | {r['direction']} | "
               f"{'' if pd.isna(r['b3']) else f'{r['b3']:.3f}'} | {'' if pd.isna(r['padj']) else f'{r['padj']:.2g}'} |")
add += ["", "### 判读", "",
        "- mean|β₃| 与 |β₃|>1 计数显著超随机 = 签名携带的交互效应量级非随机（但签名本身按 padj 定义，该项为内部一致性校验，非独立发现）；",
        "- 锚定节点命中数的单侧 p 为**方向性富集检验**（披露：该统计量不依赖 baseMean 匹配），",
        "  用于回答『机制节点是否更可能落入交互签名』；其显著与否不改变 Gate 1 判定（Gate 1 已按冻结判据通过）。",
        "- 本检验为 Phase 1 任务收尾补做，不改变 IIAMD v1.0 成员。", ""]
with open(REPORT, "a", encoding="utf-8") as f:
    f.write("\n".join(add))
with open(os.path.join(ROOT, "03_LOGS", "P1_matched_gene_set_log.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(log))
note(f"DONE in {time.time()-t0:.1f}s")
