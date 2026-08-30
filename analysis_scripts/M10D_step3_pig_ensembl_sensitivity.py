# -*- coding: utf-8 -*-
"""
M10D_step3_pig_ensembl_sensitivity.py — 猪层同源敏感性复算（Ensembl one2one 口径）
==================================================================================
背景：M10D 主口径执行时 Ensembl 服务 503，猪同源降级为 GPL16524 注释符号逐字口径
（51/63）。服务恢复窗口内缓存回填（M10D_step2_pig_ensembl_fetch.py）后，本脚本以
Ensembl REST ortholog_one2one 映射重跑 D3 猪全血 in vivo LPS 层，作为敏感性交叉。

【口径冻结——先于任何数值，写入日志】
  - Ensembl 一对一映射：缓存 JSON 的 homologies 条目按 type=="ortholog_one2one" 过滤，
    且每条人类臂基因唯一目标 id 数 == 1；MT-* 线粒体基因在 Compara 无真同源组
    （其返回为 one2many 弱旁系匹配，被 type 过滤自然剔除，如实披露）。
  - 目标 ENSSSCG id → 猪符号：/lookup/id/{id}（带重试，缓存 m10d_pig_lookup.json）。
  - 反向唯一性：同一猪符号被多个人类臂基因命中 → 全部剔除（与主脚本同款纪律）。
  - D3 重跑：与 M10D 主脚本同口径——GPL16524 探针符号 → 猪符号过滤 → 数据集内全样本
    逐基因 z → 臂均值；主对比 +24h vs t0、敏感性 +1h/+4h（个体内配对，animal 键），
    配对 t + Wilcoxon + Hedges' g。
  - 本脚本不修改主口径冻结表（Table S83），只新增 Table S83b 敏感性表并在报告中交叉。
输出：
  _intermediate/M10D_pig_ensembl_map.csv（Ensembl one2one 映射留痕）
  _intermediate/M10D_scores_pig_GSE107487_ensembl.csv（逐样本评分）
  02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV/Table_S83b_M10D_Pig_Ensembl_Sensitivity.csv
日志：03_LOGS/M10D_pig_ensembl_sensitivity_log.txt
"""
import json
import os
import re
import sys
import time
import urllib.request
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import mdi_lib as L

ROOT = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(ROOT, "00_RAW_DATA")
INTER = os.path.join(ROOT, "_intermediate")
TAB = os.path.join(ROOT, "02_SUPPLEMENTARY_TABLES", "SUPPLEMENTARY_Tables_CSV")
CACHE = os.path.join(INTER, "m10d_pig_homology")
LOGP = os.path.join(ROOT, "03_LOGS", "M10D_pig_ensembl_sensitivity_log.txt")

log = []
def note(m=""):
    log.append(m); print(m, flush=True)

def flush_log():
    with open(LOGP, "w", encoding="utf-8") as f:
        f.write("\n".join(log) + "\n")

def zrows(mat):
    sd = mat.std(axis=1, ddof=1).replace(0, np.nan)
    z = mat.sub(mat.mean(axis=1), axis=0).div(sd, axis=0)
    return z.fillna(0.0)

def fetch_lookup(eid):
    url = ("https://rest.ensembl.org/lookup/id/%s?content-type=application/json" % eid)
    for attempt in range(1, 6):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            r = json.loads(urllib.request.urlopen(req, timeout=120).read())
            if isinstance(r, dict) and r.get("display_name"):
                return r.get("display_name")
            return None
        except Exception:
            time.sleep(4 * attempt)
    return None

def main():
    note("== M10D 猪同源敏感性复算（Ensembl one2one）2026-08-27")
    note("口径冻结（先于任何数值，见脚本 docstring）：type==ortholog_one2one + 唯一目标 id + "
         "反向唯一性；MT-* 无真同源组被 type 过滤剔除；D3 重跑与主脚本同口径（配对 t/Wilcoxon/g）。")

    man, arms, _ = L.load_manifest()
    arm_up = arms["upstream_collapse"]
    arm_ex = arms["execution_induction"]
    arm_all = sorted(set(arm_up + arm_ex))

    # ---------------- 1. 缓存 → Ensembl one2one 映射 ----------------
    pairs = {}   # human_symbol -> pig_symbol
    type_counts = {}
    bad = []     # 缓存无效
    for sym in arm_all:
        cf = os.path.join(CACHE, sym + ".json")
        if not os.path.exists(cf):
            bad.append(sym + ":no_cache")
            continue
        try:
            r = json.loads(open(cf, encoding="utf-8", errors="replace").read())
        except Exception:
            bad.append(sym + ":bad_json")
            continue
        tg, types = [], []
        for d in r.get("data", []):
            for h in d.get("homologies", []):
                t = h.get("target", {})
                if t.get("id"):
                    tg.append(t["id"])
                    types.append(h.get("type"))
        o2o = [tg[i] for i in range(len(tg)) if types[i] == "ortholog_one2one"]
        type_counts[sym] = types
        if len(set(tg)) != 1 or len(o2o) != 1:
            continue  # 非一对一（多目标或非 one2one 类型）→ 剔除
        pairs[sym] = set(tg).pop()
    note(f"   type 口径过滤后待查目标 {len(pairs)} 个人类臂基因；缓存无效 {len(bad)}")
    mt_dropped = [s for s in arm_all if s.startswith("MT-") and s in type_counts
                  and "ortholog_one2one" not in type_counts[s]]
    note(f"   MT-* 基因经 type 过滤剔除 {len(mt_dropped)} 个（其返回为 one2many 弱旁系匹配）")

    # ---------------- 2. ENSSSCG id → 猪符号（lookup，带缓存） ----------------
    lf = os.path.join(INTER, "m10d_pig_lookup.json")
    look = {}
    if os.path.exists(lf):
        try:
            look = json.loads(open(lf, encoding="utf-8").read())
        except Exception:
            look = {}
    need = sorted(set(pairs.values()) - set(look.keys()))
    note(f"   lookup 待查 {len(need)} 个 ENSSSCG id")
    for i, eid in enumerate(need, 1):
        nm = fetch_lookup(eid)
        if nm:
            look[eid] = nm
        else:
            note(f"   [lookup fail] {eid}")
        if i % 10 == 0:
            with open(lf, "w", encoding="utf-8") as f:
                json.dump(look, f, ensure_ascii=False, indent=0)
    with open(lf, "w", encoding="utf-8") as f:
        json.dump(look, f, ensure_ascii=False, indent=0)

    ens_map = {}
    for hs, eid in pairs.items():
        if eid in look and look[eid]:
            ens_map[hs] = str(look[eid])
    from collections import Counter
    cnt = Counter(ens_map.values())
    ens_map = {h: s for h, s in ens_map.items() if cnt[s] == 1}
    note(f"   Ensembl one2one 猪映射（lookup+反向唯一后）：{len(ens_map)} 条")
    n_up_ens = sum(1 for g in arm_up if g in ens_map)
    n_ex_ens = sum(1 for g in arm_ex if g in ens_map)
    note(f"   上游 {n_up_ens}/30、执行 {n_ex_ens}/33")

    map_rows = [{"human_symbol": h, "pig_symbol": s, "pig_ensg": pairs[h]}
                for h, s in sorted(ens_map.items())]
    pd.DataFrame(map_rows).to_csv(os.path.join(INTER, "M10D_pig_ensembl_map.csv"), index=False)

    # ---------------- 3. D3 重跑（Ensembl 映射；与主脚本同款） ----------------
    note("\n== D3 GSE107487 重跑（Ensembl one2one 映射）")
    sig = pd.read_csv(os.path.join(RAW, "GSE107487_Pig_LPS_TimeCourse", "GSE107487_log2_normalized_signal.txt.gz"), sep="\t")
    p2s = pd.read_csv(os.path.join(RAW, "GSE107487_Pig_LPS_TimeCourse", "GPL16524_probe2symbol.csv"))
    p2s = p2s.dropna(subset=["gene_symbol"])
    sig = sig.merge(p2s[["probe_id", "gene_symbol"]], left_on="ID_REF", right_on="probe_id", how="inner")
    gpsc = [c for c in sig.columns if re.fullmatch(r"GPS\d+", c)]
    keep = ens_map
    sig = sig[sig["gene_symbol"].isin(keep.values())]
    sig["_m"] = sig[gpsc].mean(axis=1)
    sig = sig.sort_values("_m", ascending=False).drop_duplicates("gene_symbol").set_index("gene_symbol")
    mat = sig[gpsc].astype(float)
    rev = {v: k for k, v in keep.items()}
    mat.index = [rev[s] for s in mat.index]
    mat = mat.groupby(level=0).mean()
    zdf = zrows(mat)

    mtset = set(man.loc[man["is_mt_gene"].astype(str).str.lower().isin(("true", "1", "yes")), "hgnc_symbol"]) \
        if "is_mt_gene" in man.columns else {g for g in arm_up if g.startswith("MT-")}

    def score_from_z(zdf, up_syms, ex_syms):
        up = [g for g in up_syms if g in zdf.index]
        ex = [g for g in ex_syms if g in zdf.index]
        out = pd.DataFrame(index=zdf.columns)
        out["UCS"] = zdf.loc[up].mean(axis=0) if up else np.nan
        out["EIS"] = zdf.loc[ex].mean(axis=0) if ex else np.nan
        out["MDI"] = out["EIS"] - out["UCS"]
        upn = [g for g in up if g not in mtset]
        out["UCS_nomt"] = zdf.loc[upn].mean(axis=0) if upn else np.nan
        out["MDI_nomt"] = out["EIS"] - out["UCS_nomt"]
        out["n_up"], out["n_ex"] = len(up), len(ex)
        return out

    sc3 = score_from_z(zdf, arm_up, arm_ex)
    smt = gzip_open_text(os.path.join(RAW, "GSE107487_Pig_LPS_TimeCourse", "GSE107487_series_matrix.txt.gz"))
    lines = smt.split("\n")
    titles = re.split(r'\t', [l for l in lines if l.startswith("!Sample_title")][0])[1:]
    chars = [l for l in lines if l.startswith("!Sample_characteristics_ch1")]
    meta = pd.DataFrame({"gps": [t.strip('"') for t in titles]})
    for row in chars:
        cells = re.split(r'\t', row)[1:]
        m = re.match(r'"?(\w[\w ()-]*?): ', cells[0])
        if not m:
            continue
        key = m.group(1)
        vals = []
        for c in cells:
            mm = re.search(rf"{re.escape(key)}: ([^\"]+)", c)
            vals.append(mm.group(1).strip() if mm else None)
        if key not in meta.columns:
            meta[key] = [v if v else None for v in vals]
    meta = meta.rename(columns={"individual": "animal", "time point (h)": "hour"})
    sc3 = sc3.join(meta.set_index("gps")).rename(columns={"agent": "lps"})
    sc3["hour"] = pd.to_numeric(sc3["hour"], errors="coerce")
    sc3.index.name = "sample"
    sc3.reset_index().to_csv(os.path.join(INTER, "M10D_scores_pig_GSE107487_ensembl.csv"), index=False)
    note(f"   样本 {len(sc3)}；可测臂基因 UCS {sc3.n_up.iloc[0]}/30、EIS {sc3.n_ex.iloc[0]}/33；动物数 {sc3.animal.nunique()}")

    def contrast(scores, case_mask, ctrl_mask, label, paired_keys=None):
        r = {"contrast": label}
        for scc in ("UCS", "EIS", "MDI", "UCS_nomt", "MDI_nomt"):
            a = scores.loc[case_mask, scc].dropna().values
            b = scores.loc[ctrl_mask, scc].dropna().values
            g, gse = L.hedges_g(b, a)
            if paired_keys is not None:
                ka = scores.loc[case_mask, [scc]].assign(k=paired_keys[case_mask]).dropna().groupby("k")[scc].first()
                kb = scores.loc[ctrl_mask, [scc]].assign(k=paired_keys[ctrl_mask]).dropna().groupby("k")[scc].first()
                common = ka.index.intersection(kb.index)
                d = (ka[common] - kb[common]).values
                t, tp = stats.ttest_rel(ka[common], kb[common]) if len(common) > 1 else (np.nan, np.nan)
                try:
                    _, wp = stats.wilcoxon(ka[common], kb[common]) if len(common) > 2 else (np.nan, np.nan)
                except ValueError:
                    wp = np.nan
                r.update({f"{scc}_diff": float(np.mean(d)), f"{scc}_g": float(g), f"{scc}_g_se": float(gse),
                          f"{scc}_p_t": float(tp) if tp == tp else np.nan,
                          f"{scc}_p_wilcoxon": float(wp) if wp == wp else np.nan, f"{scc}_n": int(len(common))})
            else:
                t, tp = stats.ttest_ind(b, a, equal_var=False)
                _, mp = stats.mannwhitneyu(b, a, alternative="two-sided")
                r.update({f"{scc}_diff": float(a.mean() - b.mean()), f"{scc}_g": float(g), f"{scc}_g_se": float(gse),
                          f"{scc}_p_t": float(tp), f"{scc}_p_wilcoxon": float(mp), f"{scc}_n": int(len(a) + len(b))})
        r["n_up"], r["n_ex"] = int(scores["n_up"].iloc[0]), int(scores["n_ex"].iloc[0])
        return r

    ens_rows = []
    keys = sc3["animal"].values
    for hr, tag in [(24, "primary"), (1, "sensitivity"), (4, "sensitivity")]:
        case = (sc3.hour == hr) & (sc3.lps == "LPS injection")
        ctrl = (sc3.hour == 0)
        ens_rows.append({"map": "Ensembl_one2one", "design": f"paired +{hr}h vs t0 ({tag})",
                         **contrast(sc3, case.values, ctrl.values, f"LPS{hr}h_vs_t0", paired_keys=keys)})

    # ---------------- 4. 与主口径（降级注释逐字）对照 ----------------
    frozen = pd.read_csv(os.path.join(TAB, "Table_S83_M10D_CrossSpecies_Contrasts.csv"))
    pig_frozen = frozen[frozen["dataset"] == "GSE107487_pig_blood_LPS_invivo"].copy()
    pig_frozen["map"] = "GPL16524_verbatim_frozen"
    pig_frozen["design"] = pig_frozen["design"].str.replace(" vs t0", "", regex=False)
    cols_keep = ["map", "design", "contrast", "UCS_diff", "UCS_g", "UCS_p_t", "EIS_diff", "MDI_diff",
                 "MDI_g", "MDI_p_t", "UCS_nomt_diff", "MDI_nomt_diff", "n_up", "n_ex", "UCS_n"]
    ens_df = pd.DataFrame(ens_rows)
    ens_df["contrast"] = ens_df["contrast"].str.replace("LPS24h_vs_t0", "LPS24h_vs_t0")
    out = pd.concat([ens_df, pig_frozen[cols_keep]], ignore_index=True, sort=False)
    out.to_csv(os.path.join(TAB, "Table_S83b_M10D_Pig_Ensembl_Sensitivity.csv"), index=False)

    note("\n== 对照：Ensembl one2one vs 降级注释逐字（主口径）")
    for _, r in ens_df.iterrows():
        note(f"   [Ensembl {r['design']}] UCS diff={r['UCS_diff']:+.3f} (p_t={r['UCS_p_t']:.3g}) | "
             f"MDI diff={r['MDI_diff']:+.3f} (p_t={r['MDI_p_t']:.3g}) | n_up={r['n_up']} n_ex={r['n_ex']}")
    for _, r in pig_frozen.iterrows():
        note(f"   [降级口径 {r['contrast']}] UCS diff={r['UCS_diff']:+.3f} (p_t={r['UCS_p_t']:.3g}) | "
             f"MDI diff={r['MDI_diff']:+.3f} (p_t={r['MDI_p_t']:.3g}) | n_up={r['n_up']} n_ex={r['n_ex']}")
    flush_log()
    print("\nDONE")


def gzip_open_text(path):
    import gzip
    return gzip.open(path, "rt", encoding="utf-8", errors="replace").read()


if __name__ == "__main__":
    main()
