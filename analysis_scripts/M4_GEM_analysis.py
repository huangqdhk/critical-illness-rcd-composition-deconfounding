# -*- coding: utf-8 -*-
"""
M4_GEM_analysis.py — 基因组规模代谢建模（Recon3D × GSE185263 转录组）
版本 m4gem_v1.1（v1.0 的退化问题修复，修复记录见 03_LOGS\M4_GEM_v1.0_degeneracy_note.md）

v1.1 变更：
- 施加标准生长培养基（封闭其余交换反应），消除边界开放导致的容量饱和退化
  （v1.0 中 23/38 目标反应容量在全样本间 std=0）。
- 容量分析改为原地切换目标函数（不再拷贝模型），运行时间约降至 1/10。
- 铁失衡分降级为探索性指标；主读出为逐反应容量的组间比较。

方法（对应 111黄裕荣创新提质20280816.md M4 步骤 3，指南 §5.4）：
1. 以 COBRApy 载入 Recon3D（BiGG JSON，SBRG 官方仓库版）。
2. GSE185263 原始计数 → CPM → log2(CPM+1)。
3. Ensembl → Entrez 映射（mygene），把模型基因（Entrez.N 格式）接到表达矩阵。
4. GIMME 式上下文重构（逐样本）：每个反应的"支持表达"取其 GPR 内基因的均值，
   低于阈值（默认 log2(CPM+1) < 1）的反应关闭；GPR 解析失败（基因集为空）的反应保守保留。
5. 容量分析（逐样本 × 逐目标反应）：先求生物量最优，再以各目标反应为目标函数、
   约束生物量 ≥90% 最优，求解该反应的"最大可行通量"（μ_capacity）。
6. 铁稳态（摄取/储存/输出/血红素降解-合成）、GSH/GPX、ETC/OXPHOS、糖酵解/TCA 子系统
   提取容量；疾病组（Sepsis + Sepsis_COVID）vs 对照（Control）Wilcoxon + BH。
7. 探索性（不参与判定）：游离铁池失衡分、Sepsis vs Sepsis_COVID 对比。

预注册判定规则（先于结果冻结，版本 m4gem_v1.0）：
- 主判定对象：铁稳态子系统与 GSH/GPX 子系统的方向性通量变化。
- 一致性检验：模型预测的通量方向与 80 基因 manifest 中对应基因的转录方向
  （ARDS_vs_Control_log2FC 符号）逐项比对，报告方向一致率（一致数/可比对数）。
- 判定：一致率 ≥60% 且铁流入/流出失衡方向与"游离铁增加"假设同向 → M4 GEM 腿成立（正向）；
  一致率 40–60% → 部分支持；<40% 或失衡方向相反 → 如实阴性报告。
- 所有 GEM 结论措辞限定为"模型预测的待检验假说"（文档红线）。
"""

# Recon3D 标准生长培养基（COBRA toolbox 惯例 + 铁代谢测试所需的铁摄取通道）
MEDIUM = [
    "EX_glc_D[e]", "EX_o2[e]", "EX_h2o[e]", "EX_hco3[e]", "EX_nh4[e]",
    "EX_pi[e]", "EX_so4[e]", "EX_ca2[e]", "EX_cl[e]", "EX_k[e]", "EX_na1[e]",
    "EX_ala_L[e]", "EX_arg_L[e]", "EX_asn_L[e]", "EX_asp_L[e]",
    "EX_cys_L[e]", "EX_gln_L[e]", "EX_glu_L[e]", "EX_gly[e]",
    "EX_his_L[e]", "EX_ile_L[e]", "EX_leu_L[e]", "EX_lys_L[e]",
    "EX_met_L[e]", "EX_phe_L[e]", "EX_pro_L[e]", "EX_ser_L[e]",
    "EX_thr_L[e]", "EX_trp_L[e]", "EX_tyr_L[e]", "EX_val_L[e]",
    "EX_fe2[e]", "EX_fe3[e]", "EX_cu2[e]", "EX_zn2[e]", "EX_mn2[e]",
    "EX_mg2[e]", "EX_sel[e]", "EX_thm[e]", "EX_pnto_R[e]", "EX_btn[e]",
    "EX_nac[e]", "EX_fol[e]", "EX_ribflv[e]", "EX_pydxn[e]",
    "EX_adpcbl[e]", "EX_retinol_9_cis[e]", "EX_lnlc[e]", "EX_chsterol[e]",
]
import io, os, sys, argparse, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
import numpy as np
import pandas as pd
import cobra
from cobra.io import load_json_model
from cobra.flux_analysis import pfba

BASE = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS"

def find_reactions(model, ids=None, patterns=None):
    out = []
    if ids:
        out += [model.reactions.get_by_id(i) for i in ids if i in model.reactions]
    if patterns:
        out += [r for r in model.reactions
                if any(p in r.id for p in patterns) and r not in out]
    return list(dict.fromkeys(out))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", default=BASE + r"\00_RAW_DATA\GSE185263_Lung_ARDS\GSE185263_raw_counts.csv")
    ap.add_argument("--groups", default=BASE + r"\04_AUDIT_GOVERNANCE\GSE185263_groups.csv")
    ap.add_argument("--model", default=BASE + r"\00_RAW_DATA\Recon3D.json")
    ap.add_argument("--outdir", default=BASE + r"\01_FIGURE_DATA_CSV\M4_GEM")
    ap.add_argument("--threshold", type=float, default=1.0, help="log2(CPM+1) 阈值")
    ap.add_argument("--map_cache", default=BASE + r"\_intermediate\ensembl_entrez_cache.csv")
    ap.add_argument("--max_samples", type=int, default=0)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    t0 = time.time()
    print("[1/7] 载入模型 ...")
    M = load_json_model(args.model)
    M.solver = "glpk"
    # 封闭全部交换反应，再按标准培养基开放
    closed, opened = 0, 0
    for r in M.exchanges:
        if r.lower_bound < 0 and r.id not in MEDIUM:
            r.lower_bound = 0.0
            closed += 1
        elif r.id in MEDIUM:
            r.lower_bound = -1000.0
            opened += 1
    print(f"  培养基: 开放 {opened} 个交换, 封闭 {closed} 个")
    biomass = [r for r in M.reactions if "biomass" in r.id.lower()]
    print("  biomass reactions:", [r.id for r in biomass])
    bio_rxn = (M.reactions.get_by_id("biomass_reaction")
               if "biomass_reaction" in M.reactions else biomass[0])
    if bio_rxn is None:
        raise RuntimeError("未找到生物量反应")

    print("[2/7] 读入表达矩阵与分组 ...")
    cnt = pd.read_csv(args.counts, index_col=0)
    grp = pd.read_csv(args.groups, index_col=0)
    grp.columns = ["group"]
    common = [c for c in cnt.columns if c in grp.index]
    print(f"  矩阵样本 {cnt.shape[1]}, 分组 {grp.shape[0]}, 交集 {len(common)}")
    cnt = cnt[common]
    grp = grp.loc[common]
    grp["arm"] = np.where(grp["group"] == "Control", "Control", "Disease")
    print("  arm:", grp["arm"].value_counts().to_dict())

    print("[3/7] 归一化 CPM -> log2(CPM+1) ...")
    lib = cnt.sum(axis=0)
    cpm = cnt.div(lib, axis=1) * 1e6
    expr = np.log2(cpm + 1)

    print("[4/7] Ensembl -> Entrez 映射（mygene，缓存复用）...")
    if os.path.exists(args.map_cache):
        emap = pd.read_csv(args.map_cache, dtype=str)
        emap = emap.dropna(subset=["entrez"])
    else:
        import mygene
        mg = mygene.MyGeneInfo()
        ens_ids = expr.index.astype(str).tolist()
        hits = []
        for i in range(0, len(ens_ids), 500):
            chunk = ens_ids[i:i + 500]
            for attempt in range(3):
                try:
                    hits += mg.querymany(chunk, scopes="ensembl.gene",
                                         fields="entrezgene", species="human",
                                         verbose=False)
                    break
                except Exception as e:
                    if attempt == 2:
                        raise RuntimeError(f"mygene 失败: {e}") from e
                    time.sleep(5)
        emap = pd.DataFrame([
            {"ensembl": h["query"], "entrez": str(h.get("entrezgene", ""))}
            for h in hits if h.get("entrezgene")
        ])
        emap.to_csv(args.map_cache, index=False)
    print(f"  映射到 Entrez 的基因: {len(emap)}")
    # Entrez -> 矩阵行最大表达
    emap = emap[emap["ensembl"].isin(expr.index)]
    emap = emap.merge(expr.reset_index(), left_on="ensembl",
                      right_on=expr.index.name or "index", how="left")
    emap = emap.dropna(subset=["entrez"])
    entrez_expr = emap.groupby("entrez")[common].max()
    print(f"  Entrez 去重后: {len(entrez_expr)}")

    print("[5/7] 逐样本上下文重构 + 容量分析 ...")
    orig_bounds = {r.id: (r.lower_bound, r.upper_bound) for r in M.reactions}
    target_rxns = ["r1105", "FE2t", "FE2DMT1", "FE2tm", "r0016", "FE3MTP1",
                   "HOXG", "HMR_4763", "FCLTm", "GLUCYS", "GTHS",
                   "GTHP", "GTHPe", "GTHPm", "GTHO", "GTHOm",
                   "G6PDH2r", "G6PDH2c", "G6PDH2rer", "G6PDA",
                   "NADH2_u10mi", "CYOR_u10mi", "CYOOm3i", "CYOOm2i",
                   "ATPS4mi", "CSm", "ACONTm", "ICDHxm", "AKGDm",
                   "SUCOASm", "FUMm", "MDHm", "HEX1", "PFK", "GAPD",
                   "PGK", "PYK", "LDH_L"]
    target_rxns = [i for i in target_rxns if i in M.reactions]
    if args.max_samples > 0:
        common = (grp.groupby("arm").apply(lambda g: g.index[:args.max_samples])
                  .explode().tolist())
    flux_rows = {}
    iron_rows = {}
    n_skip = 0
    knock_stats = []
    old_objective = M.objective
    for s in common:
        for r in M.reactions:
            r.bounds = orig_bounds[r.id]
        e = entrez_expr[s]
        gene_expr = {}
        for g in M.genes:
            ez = g.id.split(".")[0]
            if ez in e.index:
                gene_expr[g.id] = e[ez]
        n_ko = 0
        for r in M.reactions:
            gs = [gene_expr[g.id] for g in r.genes if g.id in gene_expr]
            if not gs:
                continue
            if float(np.mean(gs)) < args.threshold:
                r.bounds = (0, 0)
                n_ko += 1
        knock_stats.append(n_ko)
        # 1) 生物量最优
        M.objective = bio_rxn
        try:
            base = M.optimize()
        except Exception:
            n_skip += 1
            continue
        if base.status != "optimal" or base.objective_value < 1e-6:
            n_skip += 1
            continue
        bio_opt = base.objective_value
        # 2) 各目标反应容量（原地切换目标函数，避免模型拷贝）
        cap = {}
        bio_lb_orig = bio_rxn.lower_bound
        bio_rxn.lower_bound = 0.9 * bio_opt
        for rid in target_rxns:
            M.objective = rid
            try:
                sol = M.optimize()
                cap[rid] = sol.objective_value if sol.status == "optimal" else 0.0
            except Exception:
                cap[rid] = 0.0
        bio_rxn.lower_bound = bio_lb_orig
        M.objective = old_objective
        flux_rows[s] = (cap, bio_opt)
        fe_in = sum(cap.get(i, 0.0) for i in
                    ["r1105", "FE2t", "FE2DMT1", "FE2tm", "HOXG", "HMR_4763"])
        fe_out = sum(cap.get(i, 0.0) for i in ["r0016", "FE3MTP1", "FCLTm"])
        iron_rows[s] = {"in": fe_in, "out": fe_out, "balance": fe_in - fe_out,
                        "biomass": bio_opt}
        if (len(flux_rows) % 50) == 0:
            print(f"  ... {len(flux_rows)}/{len(common)} samples done", flush=True)
    if knock_stats:
        print(f"  每样本被关闭的反应数: 中位 {np.median(knock_stats):.0f}, "
              f"范围 [{min(knock_stats)}-{max(knock_stats)}]")

    print("[6/7] 汇总容量矩阵 ...")
    F = pd.DataFrame({s: f for s, (f, b) in flux_rows.items() if f is not None})
    F.to_csv(os.path.join(args.outdir, "M4_GEM_capacity_matrix.csv"))
    bio_df = pd.Series({s: b for s, (f, b) in flux_rows.items() if b is not None},
                       name="biomass")
    bio_df.to_csv(os.path.join(args.outdir, "M4_GEM_biomass.csv"))
    iron_df = pd.DataFrame(iron_rows).T
    iron_df["arm"] = grp.loc[iron_df.index, "arm"]
    iron_df.to_csv(os.path.join(args.outdir, "M4_GEM_iron_balance.csv"))
    print(f"  成功样本: {F.shape[1]}/{len(common)} (跳过 {n_skip}); 生物量中位: {bio_df.median():.3f}")

    print("[7/7] 子系统提取 + 组间比较 ...")
    subsys = {
        "iron_uptake": ["r1105", "FE2t", "FE2DMT1", "FE2tm"],
        "iron_storage": ["r0016"],
        "iron_export": ["FE3MTP1"],
        "heme_degradation": ["HOXG", "HMR_4763"],
        "heme_synthesis": ["FCLTm"],
        "gsh_synthesis": ["GLUCYS", "GTHS"],
        "gsh_peroxidase": ["GTHP", "GTHPe", "GTHPm"],
        "gsh_reductase": ["GTHO", "GTHOm"],
        "nadph_ppp": ["G6PDH2r", "G6PDH2c", "G6PDH2rer", "G6PDA"],
    }
    etc_pat = ["NADH2_u10m", "SUCD1_u10m", "SUCD2_u10m", "CYOR_u10m",
               "CYOOm", "CYOR_u10mi", "ATPS4mi"]
    gly_pat = ["HEX1", "PFK", "GAPD", "PGK", "PYK", "LDH_L"]
    tca_pat = ["CSm", "ACONTm", "ICDHxm", "AKGDm", "SUCOASm", "FUMm", "MDHm"]
    for name, pat in [("ETC_OXPHOS", etc_pat), ("glycolysis", gly_pat), ("TCA", tca_pat)]:
        subsys[name] = [r.id for r in M.reactions
                        if any(p in r.id for p in pat)]
        print(f"  {name}: {len(subsys[name])} reactions -> {subsys[name][:8]}")

    rows = []
    for sub, ids in subsys.items():
        present = [i for i in ids if i in F.index]
        for i in present:
            dis = F.loc[i, grp.loc[F.columns, "arm"] == "Disease"]
            ctl = F.loc[i, grp.loc[F.columns, "arm"] == "Control"]
            from scipy.stats import mannwhitneyu
            try:
                w, p = mannwhitneyu(dis, ctl, alternative="two-sided")
            except ValueError:
                continue
            rows.append({"subsystem": sub, "reaction": i,
                         "disease_median": dis.median(), "control_median": ctl.median(),
                         "delta": dis.median() - ctl.median(),
                         "p": p, "n_disease": len(dis), "n_control": len(ctl)})
    comp = pd.DataFrame(rows)
    if len(comp):
        comp["BH"] = comp["p"] * len(comp) / comp["p"].rank(method="first")
        comp["BH"] = comp["BH"].clip(upper=1.0)
    comp.to_csv(os.path.join(args.outdir, "M4_GEM_group_comparison.csv"), index=False)

    # 铁失衡组间检验
    from scipy.stats import mannwhitneyu
    d_bal = iron_df.loc[iron_df["arm"] == "Disease", "balance"]
    c_bal = iron_df.loc[iron_df["arm"] == "Control", "balance"]
    w, p = mannwhitneyu(d_bal, c_bal, alternative="two-sided")
    print("\n===== 铁失衡分（模型预测的游离铁池变化方向）=====")
    print(f"  Disease balance median: {d_bal.median():.4f} | Control: {c_bal.median():.4f}")
    print(f"  Wilcoxon p = {p:.4g}  (Disease > Control 意味着模型预测疾病组游离铁净流入)")
    print("\n===== 组间比较（BH<0.05 的行）=====")
    sig = comp[comp["BH"] < 0.05].sort_values("BH") if len(comp) else comp
    print(sig.to_string(index=False, max_colwidth=40))
    print(f"\n完成，耗时 {time.time()-t0:.0f}s。输出目录: {args.outdir}")

if __name__ == "__main__":
    main()
