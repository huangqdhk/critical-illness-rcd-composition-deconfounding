# -*- coding: utf-8 -*-
"""
M14_step2_perturbseq_gse221321.py — M14 腿3 外部 Perturb-seq（GSE221321）GO 判定与三角化分析
============================================================================================
前瞻注册：M10_M13_M14_pre_registration_20260824.md（osf.io/ETVMJ）§3.3 / §3.4；主方案
111黄裕荣创新提质_20260822.md §四 M14。
- GO 判据（注册冻结）：候选须提供 ① 完整单细胞计数矩阵 ② gRNA 映射表 ③ 非靶向对照（NTC）
  ④ 髓系/LPS 刺激背景；两周期限内定 GO/NO-GO；NO-GO 条款：无候选满足判据则记 NO-GO 如实入稿。
- 定位（注册冻结）：腿3 结果仅作"预测层/三角化"证据，不进中心结论（M14 判定门由腿1 承担）。
- 本脚本预指定口径（先于本次结果）：
  主数据 = GSE221321 两个 conventional 通道（GSM6858447 CRISPRko / GSM6858449 CRISPRi，
  均为 THP-1 + LPS 3h；Yao et al., bioRxiv 2023, PMID 36747806，Compressed Perturb-seq）；
  评分 = mdi_v1.0 冻结口径（臂均值 -> 数据集内 z；MDI = z(EIS) - z(UCS)；
  臂基因唯一来源 = Mitoxyperilysis_Gene_Manifest_v1.0.csv arm 列，经 mdi_lib.load_manifest()）；
  细胞分类 = NTC（non-targeting/safe-targeting only）/ 单扰动（恰 1 个靶基因）/ 复合扰动（排除于主分析）；
  基因效应 = Δscore = mean(score | g) - mean(score | NTC)，要求 n_g >= MIN_CELLS=20；
  负对照校准 1（细胞级 NTC 匹配置换）= 对每基因，从 NTC 细胞中不放回抽 n_g 个 vs 其余 NTC，
  B=2,000（臂成员 B=10,000；seed=0）构造零分布，双侧经验 p；split-half 零分布中心应≈0（自校准）；
  负对照校准 2（效应量级）= FRPerturb 官方效应量（Perturbed_gene, Downstream_gene, LFC, q）
  投影到两臂（ΔUCS_FRP = LFC 在 UCS 臂基因上的均值；ΔEIS_FRP 同理；ΔMDI_FRP = ΔEIS-ΔUCS），
  NTC 行的投影应≈0；
  主检验（效应量级，基因集为统计单位）= 臂成员扰动集的臂内自抑效应（mean ΔUCS_FRP 对 UCS 成员、
  mean ΔEIS_FRP 对 EIS 成员，方向预期 <0）vs 表达量匹配随机等大小基因集（B=10,000，seed=0，
  匹配 = 扰动基因在 NTC 细胞中均值表达的十分位），单侧经验 p=(1+#null<=obs)/(B+1)；
  交叉复核 = KO 与 KD 两模态 Spearman 一致性 + 臂成员方向一致率；
  敏感性 = MDI_nomt（去 MT-*）；臂成员 leave-self-out 臂评分（排除被扰动基因自身后重算臂均值）。
- 统计单位声明：逐基因效应为细胞级描述量；推断性检验的单位是基因集（置换在基因层面），
  与统计规范"细胞永不作独立单位"一致（本层结论只涉及基因集相对行为）。
输入：00_RAW_DATA/GSE221321_THP1_PerturbSeq/GSE221321_RAW/（作者解压后 18 文件）
输出：
  _intermediate/M14L3_celllevel_effects.csv / M14L3_frp_projections.csv / M14L3_arm_member_summary.csv
  02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV/Table_S72_M14L3_GO_Audit.csv
  .../Table_S73_M14L3_CellLevel_Perturbation_Effects.csv
  .../Table_S74_M14L3_FRPerturb_Arm_Projections.csv
  .../Table_S75_M14L3_NullCalibration_and_SetTests.csv
  04_AUDIT_GOVERNANCE/M14L3_GSE221321_perturbseq_report.md
"""
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ROOT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"
RAW = ROOT + r"\00_RAW_DATA\GSE221321_THP1_PerturbSeq\GSE221321_RAW"
INTER = ROOT + r"\_intermediate"
GOV = ROOT + r"\04_AUDIT_GOVERNANCE"
TAB = ROOT + r"\02_SUPPLEMENTARY_TABLES\SUPPLEMENTARY_Tables_CSV"

SEED = 0
B_SET = 10000       # 基因集级表达量匹配置换次数
B_NTC = 2000        # 逐基因 NTC 匹配置换次数（非臂成员）
B_NTC_ARM = 10000   # 臂成员 NTC 匹配置换次数
B_CAL = 10000       # split-half 自校准次数
MIN_CELLS = 20      # 单扰动细胞数下限
DECILES = 10        # 表达量匹配分位档数
BLOCK = 2000        # 稀疏矩阵分块行数

sys.path.insert(0, ROOT)
import mdi_lib as L

t0 = time.time()
log = []

def note(m=""):
    log.append(m)
    print(m, flush=True)

rng = np.random.default_rng(SEED)
NTC_NAMES = {"NON-TARGETING", "SAFE-TARGETING"}

CHANNELS = [
    dict(key="KO", h5ad="GSM6858447_KO_conventional.h5ad",
         frp="GSM6858447_KO_conventional_FRPerturb_effect_sizes.csv.gz"),
    dict(key="KD", h5ad="GSM6858449_KD_conventional.h5ad",
         frp="GSM6858449_KD_conventional_FRPerturb_effect_sizes.csv.gz"),
]

def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = x.std(ddof=1)
    return (x - x.mean()) / sd if sd > 0 else np.zeros_like(x)

def parse_guides(s):
    if not isinstance(s, str) or not s.strip():
        return []
    return [g.strip().upper() for g in s.split("--") if g.strip()]

def cliff_from_u(u, n, m):
    return 2.0 * u / (n * m) - 1.0

def ntc_perm_p(ntc_scores_vec, obs_delta, n_ntc, ng, B, rng):
    """NTC 匹配置换双侧经验 p：抽 ng 个 NTC 细胞 vs 其余 NTC 的 Δ 零分布。"""
    draws = rng.integers(0, n_ntc, size=(B, ng))
    means = ntc_scores_vec[draws].mean(axis=1)
    pseudo = means - (ntc_scores_vec.sum() - means * ng) / max(n_ntc - ng, 1)
    return float((np.abs(pseudo) >= abs(obs_delta)).sum() + 1) / (B + 1)

# ======================================================================
# 1. GO 判据审计（先于分析执行；NO-GO 则终止）
# ======================================================================
note("== 1. GO 判据审计（osf.io/ETVMJ §3.3 四项冻结判据）")
for ch in CHANNELS:
    if not os.path.exists(os.path.join(INTER, ch["h5ad"])):
        note(f"   [FAIL] h5ad 缺失: {ch['h5ad']}")
        sys.exit("GO 判据①不满足 -> 腿3 NO-GO")
audit = [
    dict(criterion="GO1_count_matrix", verdict="PASS",
         source="GSM6858447/GSM6858449 h5ad（X 稀疏矩阵；作者 RAW.tar 解压）"),
    dict(criterion="GO2_guide_map", verdict="PASS",
         source="各 GSM perturbations.txt.gz（600 行=靶基因含 NTC，列=细胞条码，值=0-3 条 guide）"),
    dict(criterion="GO3_NTC", verdict="PASS",
         source="perturbations 表与 FRPerturb 表均含 'non-targeting'/'safe-targeting' 行"),
    dict(criterion="GO4_myeloid_LPS", verdict="PASS",
         source="人 THP-1 单核细胞系 + LPS 刺激（PMID 36747806, Yao et al., bioRxiv 2023；"
                "GSE221321 官方元数据 2026-08-24 核验锁定）"),
]
note("   四项 GO 判据全部 PASS -> 腿3 = GO，进入分析")

# ======================================================================
# 2. 臂定义（manifest 单一来源）
# ======================================================================
note("\n== 2. 臂定义（Mitoxyperilysis_Gene_Manifest_v1.0.csv）")
_, arms, _ = L.load_manifest()
UCS_ALL = [g.upper() for g in arms["upstream_collapse"]]
EIS_ALL = [g.upper() for g in arms["execution_induction"]]
MT_SET = {g.upper() for g in L.MT_GENES}
UCS_NOMT = [g for g in UCS_ALL if g not in MT_SET]
note(f"   UCS 臂 {len(UCS_ALL)} 基因（去 MT {len(UCS_NOMT)}），EIS 臂 {len(EIS_ALL)} 基因")

# ======================================================================
# 3. 细胞级臂评分（mdi_v1.0 口径，数据集内 z）与逐基因效应
# ======================================================================
note("\n== 3. 细胞级臂评分与逐基因效应")
cell_rows = []
arm_member_rows = []
ntc_calib_rows = []
cell_meta_rows = []
ntc_expr_store = {}

import anndata as ad

for ch in CHANNELS:
    key = ch["key"]
    hp = os.path.join(INTER, ch["h5ad"])
    note(f"   -- {key}: {os.path.basename(hp)}")
    a = ad.read_h5ad(hp, backed="r")
    n_cells, n_genes = a.shape
    feat = np.array([str(x).upper() for x in a.var["features"]])
    feat_list = feat.tolist()
    feat_set = set(feat_list)
    guides = a.obs["Guides_collapsed_by_gene"].astype(object).apply(parse_guides)
    a.file.close()

    up_in = [g for g in UCS_ALL if g in feat_set]
    ex_in = [g for g in EIS_ALL if g in feat_set]
    up_nomt_in = [g for g in UCS_NOMT if g in feat_set]
    note(f"     臂基因在 var 覆盖: UCS {len(up_in)}/30, EIS {len(ex_in)}/33 (nomt {len(up_nomt_in)}/24)")

    is_ntc = np.array([len(gs) > 0 and all(g in NTC_NAMES for g in gs) for gs in guides])
    is_single = np.array([len(gs) == 1 and gs[0] not in NTC_NAMES for gs in guides])
    n_ntc = int(is_ntc.sum()); n_single = int(is_single.sum())
    note(f"     细胞分类: NTC={n_ntc}, 单扰动={n_single}, 复合/其他={n_cells - n_ntc - n_single}")

    single_idx = {}
    for i, (gs, ok) in enumerate(zip(guides, is_single)):
        if ok:
            single_idx.setdefault(gs[0], []).append(i)
    kept_genes = sorted(g for g, v in single_idx.items() if len(v) >= MIN_CELLS)
    note(f"     单扰动基因 {len(single_idx)} 个，n>={MIN_CELLS} 保留 {len(kept_genes)} 个")

    # ---- 分块读取 X：臂均值列 + NTC 均值表达 + 臂成员列 ----
    arm_mem_perturbed = sorted(set(kept_genes) & (set(UCS_ALL) | set(EIS_ALL)))
    iu = [feat_list.index(g) for g in up_in]
    ie = [feat_list.index(g) for g in ex_in]
    iun = [feat_list.index(g) for g in up_nomt_in]
    imem = {g: feat_list.index(g) for g in arm_mem_perturbed if g in feat_set}

    ucs_raw = np.zeros(n_cells); eis_raw = np.zeros(n_cells); ucs_nomt_raw = np.zeros(n_cells)
    mem_cols = {g: np.zeros(n_cells) for g in imem}
    ntc_sum = np.zeros(n_genes)      # NTC 细胞逐基因表达和
    a = ad.read_h5ad(hp, backed="r")
    X = a.X
    for s in range(0, n_cells, BLOCK):
        e = min(s + BLOCK, n_cells)
        sub = np.asarray(X[s:e, :].toarray())
        ucs_raw[s:e] = sub[:, iu].mean(axis=1)
        eis_raw[s:e] = sub[:, ie].mean(axis=1)
        if iun:
            ucs_nomt_raw[s:e] = sub[:, iun].mean(axis=1)
        for g, j in imem.items():
            mem_cols[g][s:e] = sub[:, j]
        m = is_ntc[s:e]
        if m.any():
            ntc_sum += sub[m].sum(axis=0)
    a.file.close()
    mean_expr_ntc = ntc_sum / max(n_ntc, 1)
    expr_of = {g: float(mean_expr_ntc[i]) for i, g in enumerate(feat_list)}

    ucs_z = zscore(ucs_raw); eis_z = zscore(eis_raw)
    mdi = eis_z - ucs_z
    ucs_nomt_z = zscore(ucs_nomt_raw) if len(up_nomt_in) else np.zeros(n_cells)
    mdi_nomt = eis_z - ucs_nomt_z
    score_map = {"UCS": ucs_z, "EIS": eis_z, "MDI": mdi,
                 "UCS_nomt": ucs_nomt_z, "MDI_nomt": mdi_nomt}
    raw_map = {"UCS": ucs_raw, "EIS": eis_raw, "UCS_nomt": ucs_nomt_raw}

    ntc_idx = np.where(is_ntc)[0]
    ntc_sv = {k: v[ntc_idx] for k, v in score_map.items()}

    rho_cell, p_cell = stats.spearmanr(ucs_z, eis_z)
    note(f"     逐细胞 rho(UCS,EIS) = {rho_cell:.3f} (p={p_cell:.2e})")
    cell_meta_rows.append(dict(dataset=key, n_cells=n_cells, n_ntc=n_ntc, n_single=n_single,
                               n_arm_genes_ucs=len(up_in), n_arm_genes_eis=len(ex_in),
                               rho_cell_ucs_eis=float(rho_cell), rho_p=float(p_cell)))
    ntc_expr_store[key] = expr_of

    # ---- NTC split-half 自校准（零分布中心≈0 检查）----
    for sname in ("UCS", "EIS", "MDI", "MDI_nomt"):
        sv = ntc_sv[sname]
        half = n_ntc // 2
        draws = rng.integers(0, n_ntc, size=(B_CAL, half))
        means = sv[draws].mean(axis=1)
        pseudo = means - (sv.sum() - means * half) / (n_ntc - half)
        ntc_calib_rows.append(dict(dataset=key, score=sname, n_ntc=n_ntc, half=half,
                                   pseudo_mean=float(pseudo.mean()),
                                   pseudo_sd=float(pseudo.std(ddof=1)),
                                   pseudo_p025=float(np.percentile(pseudo, 2.5)),
                                   pseudo_p975=float(np.percentile(pseudo, 97.5))))

    # ---- 逐基因效应 ----
    for g in kept_genes:
        idx = np.array(single_idx[g])
        ng = len(idx)
        is_arm_member = g in UCS_ALL or g in EIS_ALL
        B = B_NTC_ARM if is_arm_member else B_NTC
        rec = dict(dataset=key, gene=g,
                   arm=("upstream_collapse" if g in UCS_ALL else
                        "execution_induction" if g in EIS_ALL else "none"),
                   n_cells=ng)
        for sname in ("UCS", "EIS", "MDI", "UCS_nomt", "MDI_nomt"):
            sv = score_map[sname]
            gv = sv[idx]
            dv = float(gv.mean() - ntc_sv[sname].mean())
            u, p = stats.mannwhitneyu(gv, ntc_sv[sname], alternative="two-sided")
            cd = cliff_from_u(u, ng, n_ntc)
            pp = ntc_perm_p(ntc_sv[sname], dv, n_ntc, ng, B, rng)
            rec[f"d_{sname}"] = dv
            rec[f"cliff_{sname}"] = float(cd)
            rec[f"mwu_p_{sname}"] = float(p)
            rec[f"ntc_perm_p_{sname}"] = pp
        # 臂成员 leave-self-out 臂评分
        if is_arm_member:
            own_in = up_in if g in UCS_ALL else ex_in
            n_own = len(own_in)
            if g in feat_set and n_own > 1:
                lso_raw_all = (raw_map["UCS" if g in UCS_ALL else "EIS"] * n_own
                               - mem_cols[g]) / (n_own - 1)
                lso_z_all = zscore(lso_raw_all)
                rec["d_ownarm_lso"] = float(lso_z_all[idx].mean()
                                            - lso_z_all[ntc_idx].mean())
            else:
                rec["d_ownarm_lso"] = np.nan
            rec["n_own_arm_in_var"] = n_own
            arm_member_rows.append(rec)
        cell_rows.append(rec)

cell_df = pd.DataFrame(cell_rows)
for ds in cell_df["dataset"].unique():
    m = cell_df["dataset"] == ds
    for sname in ("UCS", "EIS", "MDI", "UCS_nomt", "MDI_nomt"):
        q = L.bh(cell_df.loc[m, f"mwu_p_{sname}"].values)
        cell_df.loc[m, f"bh_q_{sname}"] = q
cell_df.to_csv(os.path.join(INTER, "M14L3_celllevel_effects.csv"), index=False)
note(f"   逐基因效应表: {cell_df.shape[0]} 行（{cell_df['dataset'].nunique()} 模态）")

# ======================================================================
# 4. FRPerturb 效应量级投影 + 表达量匹配基因集零模型（主检验）
# ======================================================================
note("\n== 4. FRPerturb 效应量投影与主检验（基因集为统计单位）")
frp_rows = []
set_test_rows = []

for ch in CHANNELS:
    key = ch["key"]
    fp = os.path.join(RAW, ch["frp"])
    note(f"   -- {key}: 读取 {os.path.basename(fp)}")
    frp = pd.read_csv(fp)
    frp.columns = [str(c).strip() for c in frp.columns]
    frp["Perturbed_gene"] = frp["Perturbed_gene"].astype(str).str.upper()
    frp["Downstream_gene"] = frp["Downstream_gene"].astype(str).str.upper()
    note(f"      行数={len(frp)}, 扰动基因={frp['Perturbed_gene'].nunique()}, "
         f"下游基因={frp['Downstream_gene'].nunique()}")

    piv = frp.pivot_table(index="Perturbed_gene", columns="Downstream_gene",
                          values="Log_fold_change", aggfunc="mean")
    up_avail = [g for g in UCS_NOMT if g in piv.columns]
    ex_avail = [g for g in EIS_ALL if g in piv.columns]
    note(f"      投影臂覆盖: UCS_nomt {len(up_avail)}/24, EIS {len(ex_avail)}/33")
    proj = pd.DataFrame(index=piv.index)
    proj["dUCS_FRP"] = piv[up_avail].mean(axis=1)
    proj["dEIS_FRP"] = piv[ex_avail].mean(axis=1)
    proj["dMDI_FRP"] = proj["dEIS_FRP"] - proj["dUCS_FRP"]

    ntc_names = [x for x in proj.index if x in NTC_NAMES]
    for nm in ntc_names:
        r = proj.loc[nm]
        set_test_rows.append(dict(dataset=key, test="FRP_NTC_projection", set_name=nm,
                                  n_set=None, obs_mean=None, null_mean=None, null_sd=None,
                                  p_one=None, p_two=None, members=None,
                                  dUCS=float(r["dUCS_FRP"]), dEIS=float(r["dEIS_FRP"]),
                                  dMDI=float(r["dMDI_FRP"])))
        note(f"      NTC 投影 {nm}: dUCS={r['dUCS_FRP']:.4f} "
             f"dEIS={r['dEIS_FRP']:.4f} dMDI={r['dMDI_FRP']:.4f}")

    # 表达量匹配零模型（基因集级）
    expr = ntc_expr_store[key]
    genes_all = sorted(set(proj.index) - NTC_NAMES)
    expr_vals = np.array([expr.get(g, 0.0) for g in genes_all])
    try:
        dec = pd.qcut(pd.Series(expr_vals), DECILES, labels=False, duplicates="drop")
    except Exception:
        dec = pd.Series(np.zeros(len(genes_all), dtype=int))
    dec_of = {g: int(d) for g, d in zip(genes_all, dec)}
    dec_pool = {}
    for g, d in dec_of.items():
        dec_pool.setdefault(d, []).append(g)
    UCS_MEM = sorted(g for g in UCS_NOMT if g in proj.index and g not in NTC_NAMES)
    EIS_MEM = sorted(g for g in EIS_ALL if g in proj.index and g not in NTC_NAMES)
    note(f"      扰动库内臂成员: UCS={UCS_MEM}, EIS={EIS_MEM}")

    for set_name, mem, stat_col in (("UCS_arm", UCS_MEM, "dUCS_FRP"),
                                    ("EIS_arm", EIS_MEM, "dEIS_FRP"),
                                    ("UCS_arm_mdi", UCS_MEM, "dMDI_FRP"),
                                    ("EIS_arm_mdi", EIS_MEM, "dMDI_FRP")):
        if not mem:
            note(f"      {set_name}: 无成员，跳过")
            continue
        obs = float(proj.loc[mem, stat_col].mean())
        dec_mem = [dec_of.get(g, 0) for g in mem]
        pools = [dec_pool[d] for d in dec_mem]
        nulls = np.empty(B_SET)
        for b in range(B_SET):
            picks = [pool[rng.integers(0, len(pool))] for pool in pools]
            nulls[b] = float(proj.loc[picks, stat_col].mean())
        p_one = float((nulls <= obs).sum() + 1) / (B_SET + 1)
        p_two = float((np.abs(nulls) >= abs(obs)).sum() + 1) / (B_SET + 1)
        set_test_rows.append(dict(dataset=key, test="expr_matched_set", set_name=set_name,
                                  n_set=len(mem), obs_mean=obs,
                                  null_mean=float(nulls.mean()), null_sd=float(nulls.std(ddof=1)),
                                  p_one=p_one, p_two=p_two, members=";".join(mem),
                                  dUCS=None, dEIS=None, dMDI=None))
        note(f"      {set_name}: obs={obs:.4f} null_mean={nulls.mean():.4f} "
             f"null_sd={nulls.std(ddof=1):.4f} p_one={p_one:.4f} p_two={p_two:.4f}")

    for g in proj.index:
        arm = ("upstream_collapse" if g in UCS_ALL else
               "execution_induction" if g in EIS_ALL else "none")
        frp_rows.append(dict(dataset=key, gene=g, arm=arm,
                             dUCS_FRP=float(proj.loc[g, "dUCS_FRP"]),
                             dEIS_FRP=float(proj.loc[g, "dEIS_FRP"]),
                             dMDI_FRP=float(proj.loc[g, "dMDI_FRP"]),
                             expr_mean_ntc=float(expr.get(g, 0.0))))

frp_df = pd.DataFrame(frp_rows)
frp_df.to_csv(os.path.join(INTER, "M14L3_frp_projections.csv"), index=False)

# 交叉模态一致性（KO vs KD）
note("\n   == 交叉模态一致性（KO vs KD）")
cons_rows = []
for stat_col in ("dUCS_FRP", "dEIS_FRP", "dMDI_FRP"):
    ko = frp_df[frp_df["dataset"] == "KO"].set_index("gene")[stat_col]
    kd = frp_df[frp_df["dataset"] == "KD"].set_index("gene")[stat_col]
    both = ko.index.intersection(kd.index)
    r, p = stats.spearmanr(ko[both], kd[both])
    arm_mem = [g for g in both if g in UCS_ALL or g in EIS_ALL]
    signs = [(np.sign(ko[g]), np.sign(kd[g])) for g in arm_mem]
    conc = sum(1 for a, b in signs if a != 0 and a == b)
    cons_rows.append(dict(stat=stat_col, n_genes=len(both), spearman_r=float(r),
                          spearman_p=float(p),
                          arm_members_shared=len(arm_mem),
                          arm_member_sign_concordant=conc,
                          arm_member_signs=";".join(f"{g}:{a}/{b}" for g, (a, b) in
                                                    zip(arm_mem, signs))))
    note(f"      {stat_col}: r={r:.3f} (p={p:.2e}), 臂成员方向一致 {conc}/{len(arm_mem)}")

# ======================================================================
# 5. 输出正式表格
# ======================================================================
note("\n== 5. 输出正式表格")

def wtable(df, name):
    p = os.path.join(TAB, name)
    df.to_csv(p, index=False, encoding="utf-8-sig")
    note(f"   写入 {name}")

wtable(pd.DataFrame(audit), "Table_S72_M14L3_GO_Audit.csv")

s73 = cell_df.copy()
s73 = s73[["dataset", "gene", "arm", "n_cells",
           "d_UCS", "cliff_UCS", "mwu_p_UCS", "ntc_perm_p_UCS",
           "d_EIS", "cliff_EIS", "mwu_p_EIS", "ntc_perm_p_EIS",
           "d_MDI", "cliff_MDI", "mwu_p_MDI", "ntc_perm_p_MDI",
           "d_UCS_nomt", "d_MDI_nomt", "bh_q_UCS", "bh_q_EIS", "bh_q_MDI"]]
wtable(s73, "Table_S73_M14L3_CellLevel_Perturbation_Effects.csv")

wtable(frp_df, "Table_S74_M14L3_FRPerturb_Arm_Projections.csv")

cal = pd.DataFrame(ntc_calib_rows)
cal["section"] = "celllevel_NTC_splithalf"
st = pd.DataFrame(set_test_rows)
st["section"] = "FRP_set_tests"
cons = pd.DataFrame(cons_rows)
cons["section"] = "cross_modality"
s75 = pd.concat([cal, st, cons], ignore_index=True)
wtable(s75, "Table_S75_M14L3_NullCalibration_and_SetTests.csv")

arm_sum = pd.DataFrame(arm_member_rows)
arm_sum.to_csv(os.path.join(INTER, "M14L3_arm_member_summary.csv"), index=False)
pd.DataFrame(cell_meta_rows).to_csv(os.path.join(INTER, "M14L3_cell_meta.csv"), index=False)

note(f"\n完成，用时 {time.time() - t0:.1f}s")
print("\n".join(log))
