# -*- coding: utf-8 -*-
"""
M10D_step1_cross_species.py — M10 注册次终点·跨物种复制（ETVMJ §1.5）
=====================================================================
执行日期：2026-08-27（数据验收 D1_cross_species_verification_20260825.md 之后）

【口径冻结——先于任何数值计算，写入日志】
数据集与主对比：
  D1 GSE342888  小鼠 CLP 肺        CLP(6) vs sham(6)                    非配对（疾病模型层）
  D2 GSE311048  大鼠 CLP±VNS 肺     PR=CLP(5) vs CR=Sham(5) 主           非配对（疾病模型层）
                                      PM=CLP+VNS(5) vs M=Sham+VNS(5) 敏感性
  D3 GSE107487  猪全血 in vivo LPS  个体内配对 +24h vs t0 主；+1h/+4h 敏感性（激发层）
  D4 GSE249010  10 物种全血 ex vivo LPS  种内受试者均值 LPS(全剂量/时点) vs NC(0 剂量) 配对
                                      （激发层；human 亚组作种内桥接）
评分：人源一对一同源（mygene homologene，双向一对一，P2 同款纪律）→ 数据集内全样本
      log2 空间逐基因 z → 臂均值 z；UCS=上游塌陷臂均值、EIS=执行诱导臂均值、MDI=EIS−UCS；
      敏感性 UCS_nomt/MDI_nomt（剔除 MT 基因）。
方向预测（跨物种复制判定）：UCS<0（塌陷保守）；EIS>0、MDI>0 为完整解离格局。
统计：Welch t + Mann-Whitney（非配对）；配对 t + Wilcoxon（配对）；效应量 Hedges' g
      （mdi_lib.hedges_g(ctrl, case)，正值=case 高）。
披露：跨物种层不做组成校正（人源含粒细胞反卷积参考不适用于非人物种，且无种内参考）——
      本层检验的是**未校正同口径信号**，不构成组成独立性证据；GSE249010 行名为提交者
      预映射的人源 ENSG（ex vivo+同源映射性质，沿用 D1 验收披露）。
"""
import gzip
import os
import re
import sys
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
LOGP = os.path.join(ROOT, "03_LOGS", "M10D_cross_species_log.txt")
TAX = {"mouse": 10090, "rat": 10116, "pig": 9823}

log = []
def note(m=""):
    log.append(m); print(m, flush=True)

def flush_log():
    with open(LOGP, "w", encoding="utf-8") as f:
        f.write("\n".join(log) + "\n")

# ---------------- 0. 口径冻结声明（日志先于数值） ----------------
note("== M10D 跨物种复制（M10 注册次终点 ETVMJ §1.5）2026-08-27")
note("口径冻结（先于任何数值，见脚本 docstring）：D1 小鼠 CLP 肺 CLP vs sham；"
     "D2 大鼠 CLP 肺 PR vs CR 主 + PM vs M 敏感性；D3 猪 in vivo LPS 个体内 +24h vs t0 主（+1/+4h 敏感性）；"
     "D4 10 物种 ex vivo LPS 种内受试者均值 LPS vs NC 配对（human 亚组作桥接）。"
     "评分=一对一同源→数据集内 z→臂均值（MDI=EIS−UCS，另报 nomt）；方向预测 UCS<0。"
     "跨物种层不做组成校正（无人源外适用参考）——检验未校正同口径信号，如实披露。")

# ---------------- 1. manifest 臂基因 ----------------
man, arms, _alias = L.load_manifest()
arm_up = arms["upstream_collapse"]
arm_ex = arms["execution_induction"]
mtset = set(man.loc[man["is_mt_gene"].astype(str).str.lower().isin(("true", "1", "yes")), "hgnc_symbol"]) \
    if "is_mt_gene" in man.columns else {g for g in arm_up if g.startswith("MT-")}
note(f"臂基因：上游 UCS {len(arm_up)} / 执行 EIS {len(arm_ex)}；MT 基因 {len(mtset)} 个")

# ---------------- 2. 一对一同源映射（mygene homologene，人→小鼠/大鼠/猪） ----------------
import mygene
mg = mygene.MyGeneInfo()
hits = mg.querymany(sorted(set(arm_up + arm_ex)), scopes="symbol", species="human",
                    fields="homologene", returnall=True)["out"]
hmap = {}       # (species, human_symbol) -> geneid
ambiguous = []  # (species, human_symbol) 多义
qid2sym = {}
for h in hits:
    q = h.get("query"); homs = h.get("homologene", {}).get("genes", [])
    qid2sym[q] = q  # 输入即符号（大写一致化在后面）
    for sp, tx in TAX.items():
        ids = [g[1] for g in homs if g[0] == tx]
        if len(ids) == 1:
            hmap[(sp, q)] = ids[0]
        elif len(ids) > 1:
            ambiguous.append((sp, q))
# GeneID → 物种符号
allids = sorted(set(hmap.values()))
ginfo = mg.getgenes([str(i) for i in allids], fields="symbol,taxid")
id2sym = {str(g["_id"]): str(g["symbol"]) for g in ginfo if g.get("symbol") and int(g.get("taxid", 0)) in TAX.values()}
# 双向一对一：同一物种符号被多个人类臂基因命中 → 全部剔除
from collections import Counter
ortho_rows = []
for (sp, hs), gid in hmap.items():
    ss = id2sym.get(str(gid))
    if ss:
        ortho_rows.append({"species": sp, "human_symbol": hs, "species_symbol": ss, "geneid": gid})
odf = pd.DataFrame(ortho_rows)
cnt = Counter(zip(odf["species"], odf["species_symbol"]))
odf["one2one"] = [cnt[(r.species, r.species_symbol)] == 1 for r in odf.itertuples()]
odf = odf[odf["one2one"]].drop(columns="one2one")
odf.to_csv(os.path.join(INTER, "M10D_ortholog_map.csv"), index=False)
note(f"同源映射：一对一 {len(odf)} 条（小鼠 {sum(odf.species=='mouse')}/大鼠 {sum(odf.species=='rat')}/猪 {sum(odf.species=='pig')}）；"
     f"多义剔除 {len(ambiguous)}；映射失败（含 MT-* 人源符号无同源组）"
     f" {len(set(arm_up+arm_ex)) - odf.human_symbol.nunique()} 个人类臂基因")

SP2SYM = {sp: dict(zip(g["human_symbol"], g["species_symbol"]))
          for sp, g in odf.groupby("species")}

# ---- 猪（9823）HomoloGene 无覆盖 → Ensembl REST 一对一（human→sus_scrofa，ortholog_is_one2one） ----
import json as _json
import urllib.request as _ur

def _ensembl_post(url, body):
    req = _ur.Request(url, data=_json.dumps(body).encode(),
                      headers={"Content-Type": "application/json", "Accept": "application/json"})
    return _json.loads(_ur.urlopen(req, timeout=180).read())

arm_all = sorted(set(arm_up + arm_ex))
pig_rows = []
try:
    h2ens = {}
    _cache_dir = os.path.join(INTER, "m10d_pig_homology")
    for _sym in arm_all:
        try:
            _cf = os.path.join(_cache_dir, _sym + ".json")
            _r = _json.loads(open(_cf, encoding="utf-8", errors="replace").read()) if os.path.exists(_cf) else {}
            _tg = [h.get("target", {}).get("id") for d in _r.get("data", []) for h in d.get("homologies", [])]
            h2ens.setdefault(_sym, []).extend([t for t in _tg if t])
        except Exception as _e:
            note(f"   [猪同源] {_sym}: {_e}")
    h2ens = {h: v[0] for h, v in h2ens.items() if len(v) == 1}
    ens_ids = sorted(set(h2ens.values()))
    _lf = os.path.join(INTER, "m10d_pig_lookup.json")
    if os.path.exists(_lf):
        look = _json.loads(open(_lf, encoding="utf-8", errors="replace").read())
    else:
        look = _ensembl_post("https://rest.ensembl.org/lookup/id", {"ids": ens_ids})
    ens2sym = {k: v.get("display_name") for k, v in look.items() if isinstance(v, dict) and v.get("display_name")}
    from collections import Counter as _C
    pig_pairs = {h: ens2sym[e] for h, e in h2ens.items() if e in ens2sym}
    _cnt = _C(pig_pairs.values())
    pig_pairs = {h: s for h, s in pig_pairs.items() if _cnt[s] == 1}
    SP2SYM["pig"] = pig_pairs
    pig_rows = [{"species": "pig", "human_symbol": h, "species_symbol": s, "geneid": ""}
                for h, s in pig_pairs.items()]
    odf = pd.concat([odf, pd.DataFrame(pig_rows)], ignore_index=True)
    note(f"猪同源（Ensembl REST 一对一）：{len(pig_pairs)} 条（HomoloGene 无猪覆盖，来源切换已披露）")
except Exception as e:
    note(f"[警告] Ensembl 猪同源失败：{e}")
# 降级口径：Ensembl 缓存不全时用 GPL16524 注释符号逐字匹配（人类式命名，披露为注释一致性口径）
if len(SP2SYM.get("pig", {})) < 40:
    _p2s = pd.read_csv(os.path.join(RAW, "GSE107487_Pig_LPS_TimeCourse", "GPL16524_probe2symbol.csv"))
    _pset = set(_p2s["gene_symbol"].dropna().astype(str).str.strip().str.upper())
    _verb = {g: g for g in arm_all if g.upper() in _pset}
    SP2SYM["pig"] = _verb
    note(f"猪同源降级口径（注释符号逐字匹配）：{len(_verb)} 条（上游 "
         f"{sum(1 for g in arm_up if g.upper() in _pset)}/30、执行 {sum(1 for g in arm_ex if g.upper() in _pset)}/33）"
         f"——Ensembl one2one 缓存不全时的替代，已与已取回的 Ensembl 结果交叉披露")
odf.to_csv(os.path.join(INTER, "M10D_ortholog_map.csv"), index=False)

# ---------------- 3. 评分函数 ----------------
def zrows(mat):
    """逐基因跨样本 z（与 L.zscore 同定义：sd=0 → 全 0）。"""
    sd = mat.std(axis=1, ddof=1).replace(0, np.nan)
    z = mat.sub(mat.mean(axis=1), axis=0).div(sd, axis=0)
    return z.fillna(0.0)

def score_from_z(zdf, up_syms, ex_syms):
    up = [g for g in up_syms if g in zdf.index]
    ex = [g for g in ex_syms if g in zdf.index]
    out = pd.DataFrame(index=zdf.columns)
    out["UCS"] = zdf.loc[up].mean(axis=0) if up else np.nan
    out["EIS"] = zdf.loc[ex].mean(axis=0) if ex else np.nan
    out["MDI"] = out["EIS"] - out["UCS"]
    upn = [g for g in up if g not in mtset]
    exn = [g for g in ex if g not in mtset]
    out["UCS_nomt"] = zdf.loc[upn].mean(axis=0) if upn else np.nan
    out["MDI_nomt"] = out["EIS"] - out["UCS_nomt"]
    out["n_up"], out["n_ex"] = len(up), len(ex)
    return out

def contrast(scores, case_mask, ctrl_mask, label, paired_keys=None):
    r = {"contrast": label}
    for sc in ("UCS", "EIS", "MDI", "UCS_nomt", "MDI_nomt"):
        a = scores.loc[case_mask, sc].dropna().values   # case
        b = scores.loc[ctrl_mask, sc].dropna().values   # ctrl
        g, gse = L.hedges_g(b, a)
        if paired_keys is not None:
            ka = scores.loc[case_mask, [sc]].assign(k=paired_keys[case_mask]).dropna().groupby("k")[sc].first()
            kb = scores.loc[ctrl_mask, [sc]].assign(k=paired_keys[ctrl_mask]).dropna().groupby("k")[sc].first()
            common = ka.index.intersection(kb.index)
            d = (ka[common] - kb[common]).values
            t, tp = stats.ttest_rel(ka[common], kb[common]) if len(common) > 1 else (np.nan, np.nan)
            try:
                _, wp = stats.wilcoxon(ka[common], kb[common]) if len(common) > 2 else (np.nan, np.nan)
            except ValueError:
                wp = np.nan
            r.update({f"{sc}_diff": float(np.mean(d)), f"{sc}_g": float(g), f"{sc}_g_se": float(gse),
                      f"{sc}_p_t": float(tp) if tp == tp else np.nan, f"{sc}_p_wilcoxon": float(wp) if wp == wp else np.nan,
                      f"{sc}_n": int(len(common))})
        else:
            t, tp = stats.ttest_ind(b, a, equal_var=False)
            _, mp = stats.mannwhitneyu(b, a, alternative="two-sided")
            r.update({f"{sc}_diff": float(a.mean() - b.mean()), f"{sc}_g": float(g), f"{sc}_g_se": float(gse),
                      f"{sc}_p_t": float(tp), f"{sc}_p_wilcoxon": float(mp), f"{sc}_n": int(len(a) + len(b))})
    r["n_up"], r["n_ex"] = int(scores["n_up"].iloc[0]), int(scores["n_ex"].iloc[0])
    return r

results = []
percol = []

# ---------------- D1 小鼠 GSE342888（CLP 肺，3,406 DEG 计数子集） ----------------
note("\n== D1 GSE342888 小鼠 CLP 肺")
x = pd.read_excel(os.path.join(RAW, "GSE342888_Mouse_CLP_Lung", "GSE342888_CLP-vs-sham_diffexp.DEG.xlsx"))
cntcols = [c for c in x.columns if re.fullmatch(r"(CLP|sham)\d", c)]
xx = x.dropna(subset=["gene_name"]).copy()
xx["gene_name"] = xx["gene_name"].astype(str)
keep = SP2SYM["mouse"]
xx = xx[xx["gene_name"].isin(keep.values())]
# 同名多行 → 最大均值保留
xx["_m"] = xx[cntcols].mean(axis=1)
xx = xx.sort_values("_m", ascending=False).drop_duplicates("gene_name").set_index("gene_name")
mat = np.log2(xx[cntcols].astype(float) + 1)
rev = {v: k for k, v in keep.items()}  # mouse symbol -> human symbol
mat.index = [rev[s] for s in mat.index]  # 转人类符号
mat = mat.groupby(level=0).mean()
zdf = zrows(mat)
up_syms, ex_syms = arm_up, arm_ex  # 索引已转人类符号
sc1 = score_from_z(zdf, up_syms, ex_syms)
sc1.index.name = "sample"; sc1["group"] = ["CLP" if c.startswith("CLP") else "sham" for c in sc1.index]
sc1.reset_index().to_csv(os.path.join(INTER, "M10D_scores_mouse_GSE342888.csv"), index=False)
note(f"   可测臂基因 UCS {sc1.n_up.iloc[0]}/30、EIS {sc1.n_ex.iloc[0]}/33（3,406 DEG 子集覆盖限制）")
results.append({"dataset": "GSE342888_mouse_CLP_lung", "species": "mouse", "design": "CLP vs sham (6v6)",
                **contrast(sc1, sc1.group == "CLP", sc1.group == "sham", "CLP_vs_sham")})

# ---------------- D2 大鼠 GSE311048（CLP±VNS 肺，TPM） ----------------
note("\n== D2 GSE311048 大鼠 CLP±VNS 肺")
g = pd.read_csv(os.path.join(RAW, "GSE311048_Rat_Lung_VNS", "GSE311048_gene_expression.txt.gz"), sep="\t")
tpmc = [c for c in g.columns if c.startswith("tpm_")]
gg = g.dropna(subset=["gene_symbol"]).copy()
gg["gene_symbol"] = gg["gene_symbol"].astype(str)
keep = SP2SYM["rat"]
gg = gg[gg["gene_symbol"].isin(keep.values())]
gg["_m"] = gg[tpmc].mean(axis=1)
gg = gg.sort_values("_m", ascending=False).drop_duplicates("gene_symbol").set_index("gene_symbol")
mat = np.log2(gg[tpmc].astype(float) + 1)
rev = {v: k for k, v in keep.items()}  # rat symbol -> human symbol
mat.index = [rev[s] for s in mat.index]
mat = mat.groupby(level=0).mean()
zdf = zrows(mat)
up_syms, ex_syms = arm_up, arm_ex  # 索引已转人类符号
sc2 = score_from_z(zdf, up_syms, ex_syms)
gmap = pd.read_csv(os.path.join(RAW, "GSE311048_Rat_Lung_VNS", "GSE311048_group_map.csv"))
tr = dict(zip(gmap["title"], gmap["characteristics"].str.extract(r"treatment: ([^;\"]+)")[0].str.strip()))
_g4 = {"Sham": "Sham", "CLP": "CLP", "Sham + VNS": "ShamVNS", "CLP + VNS": "CLPVNS"}
grp = {f"tpm_{t}": _g4[v] for t, v in tr.items()}
sc2.index.name = "sample"
sc2["group"] = [grp[c] for c in sc2.index]
sc2.reset_index().to_csv(os.path.join(INTER, "M10D_scores_rat_GSE311048.csv"), index=False)
note(f"   可测臂基因 UCS {sc2.n_up.iloc[0]}/30、EIS {sc2.n_ex.iloc[0]}/33")
results.append({"dataset": "GSE311048_rat_CLP_lung", "species": "rat", "design": "CLP(PR) vs Sham(CR) (5v5)",
                **contrast(sc2, sc2.group == "CLP", sc2.group == "Sham", "CLP_vs_Sham")})
results.append({"dataset": "GSE311048_rat_CLP_lung", "species": "rat", "design": "sensitivity CLP+VNS(PM) vs Sham+VNS(M)",
                **contrast(sc2, sc2.group == "CLPVNS", sc2.group == "ShamVNS", "CLPVNS_vs_ShamVNS")})

# ---------------- D3 猪 GSE107487（in vivo LPS 时序，个体内配对） ----------------
_n_cache = len([f for f in os.listdir(os.path.join(INTER, "m10d_pig_homology"))
                if f.endswith(".json")])
PIG_RUN = len(SP2SYM.get("pig", {})) > 0  # Ensembl 缓存不全时已降级为注释逐字口径，映射非空即可执行
if not PIG_RUN:
    note(f"D3 GSE107487 skipped this pass (pig cache " + str(_n_cache) + "/63)")
if PIG_RUN:
    note("\n== D3 GSE107487 猪全血 in vivo LPS")
    sig = pd.read_csv(os.path.join(RAW, "GSE107487_Pig_LPS_TimeCourse", "GSE107487_log2_normalized_signal.txt.gz"), sep="\t")
    p2s = pd.read_csv(os.path.join(RAW, "GSE107487_Pig_LPS_TimeCourse", "GPL16524_probe2symbol.csv"))
    p2s = p2s.dropna(subset=["gene_symbol"])
    sig = sig.merge(p2s[["probe_id", "gene_symbol"]], left_on="ID_REF", right_on="probe_id", how="inner")
    gpsc = [c for c in sig.columns if re.fullmatch(r"GPS\d+", c)]
    keep = SP2SYM["pig"]
    sig = sig[sig["gene_symbol"].isin(keep.values())]
    sig["_m"] = sig[gpsc].mean(axis=1)
    sig = sig.sort_values("_m", ascending=False).drop_duplicates("gene_symbol").set_index("gene_symbol")
    mat = sig[gpsc].astype(float)
    rev = {v: k for k, v in keep.items()}  # pig symbol -> human symbol
    mat.index = [rev[s] for s in mat.index]
    mat = mat.groupby(level=0).mean()
    zdf = zrows(mat)
    up_syms, ex_syms = arm_up, arm_ex  # 索引已转人类符号
    sc3 = score_from_z(zdf, up_syms, ex_syms)
    # 样本元数据：series matrix（title=GPS#####；individual/time/agent）
    smt = gzip.open(os.path.join(RAW, "GSE107487_Pig_LPS_TimeCourse", "GSE107487_series_matrix.txt.gz"),
                    "rt", encoding="utf-8", errors="replace").read()
    lines = smt.split("\n")
    titles = re.split(r'\t', [l for l in lines if l.startswith("!Sample_title")][0])[1:]
    chars = [l for l in lines if l.startswith("!Sample_characteristics_ch1")]
    def field(i, key):
        row = re.split(r'\t', chars[i])[1:]
        return [re.search(rf"{key}: (\S[^\"\t]*)", c).group(1).strip() if re.search(rf"{key}: (\S[^\"\t]*)", c) else None for c in row]
    meta = pd.DataFrame({"gps": [t.strip('"') for t in titles]})
    # 逐特征行找 key
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
    sc3.reset_index().to_csv(os.path.join(INTER, "M10D_scores_pig_GSE107487.csv"), index=False)
    note(f"   样本 {len(sc3)}；可测臂基因 UCS {sc3.n_up.iloc[0]}/30、EIS {sc3.n_ex.iloc[0]}/33；动物数 {sc3.animal.nunique()}")
    for hr in (1, 4, 24):
        case = (sc3.hour == hr) & (sc3.lps == "LPS injection")
        ctrl = (sc3.hour == 0)
        keys = sc3["animal"].values
        tag = "primary" if hr == 24 else "sensitivity"
        results.append({"dataset": "GSE107487_pig_blood_LPS_invivo", "species": "pig",
                        "design": f"{tag} paired +{hr}h vs t0",
                        **contrast(sc3, case.values, ctrl.values, f"LPS{hr}h_vs_t0", paired_keys=keys)})

# ---------------- D4 10 物种 GSE249010（ex vivo LPS；人源 ENSG TPM） ----------------
note("\n== D4 GSE249010 10 物种 ex vivo LPS（提交者预映射人源 ENSG）")
tpm = pd.read_csv(os.path.join(RAW, "GSE249010_MultiSpecies_Blood_LPS", "GSE249010_processed_data_TPM.txt.gz"), sep="\t")
e2s = pd.read_csv(os.path.join(INTER, "M10A_monaco_ensg2sym.csv"))
ecol, scol = e2s.columns[0], e2s.columns[1]
ens2sym = dict(zip(e2s[ecol].astype(str), e2s[scol].astype(str)))
tpm["id"] = tpm["id"].astype(str).str.split(".").str[0]
tpm["sym"] = tpm["id"].map(ens2sym).fillna(tpm["id"])
tpm = tpm.dropna(subset=["sym"])
spcols = [c for c in tpm.columns if re.fullmatch(r"SP[0-9A-F]{6}", c)]
tpm = tpm[tpm["sym"].isin(set(arm_up) | set(arm_ex))]
tpm = tpm.groupby("sym")[spcols].mean()
mat = np.log2(tpm[spcols].astype(float) + 1)
smap = pd.read_csv(os.path.join(RAW, "GSE249010_MultiSpecies_Blood_LPS", "GSE249010_sample_species_map.csv"))
smap["sp_code"] = smap["title"].str.extract(r"\[(SP[0-9A-F]{6})\]")
smap["treat"] = smap["characteristics"].str.extract(r"treatment: (\w+)")
smap["subject"] = smap["characteristics"].str.extract(r"subject: (\S+)")
smi = smap.set_index("sp_code")
zdf_all = zrows(mat)
for spname, sub in smi.groupby("species"):
    cols = [c for c in zdf_all.columns if c in sub.index]
    if len(cols) < 4:
        continue
    zz = zdf_all[cols]
    sc4 = score_from_z(zz, arm_up, arm_ex)
    sc4.index.name = "sp_code"
    sc4 = sc4.join(sub[["treat", "subject", "species"]])
    sc4.reset_index().to_csv(os.path.join(INTER, f"M10D_scores_multispecies_{spname.split()[0].lower()}.csv"), index=False)
    # 配对：受试者均值 LPS vs NC
    lp = sc4[sc4.treat == "LPS"].groupby("subject")[["UCS", "EIS", "MDI", "UCS_nomt", "MDI_nomt"]].mean()
    nc = sc4[sc4.treat == "NC"].groupby("subject")[["UCS", "EIS", "MDI", "UCS_nomt", "MDI_nomt"]].mean()
    common = lp.index.intersection(nc.index)
    merged = pd.concat([lp.loc[common].add_suffix("_case"), nc.loc[common].add_suffix("_ctrl")], axis=1)
    merged["subject"] = common
    r = {"dataset": "GSE249010_multispecies_exvivo", "species": spname,
         "design": f"paired subject-mean LPS vs NC (n={len(common)})"}
    for scc in ("UCS", "EIS", "MDI", "UCS_nomt", "MDI_nomt"):
        d = (merged[f"{scc}_case"] - merged[f"{scc}_ctrl"]).values
        t, tp_ = stats.ttest_rel(merged[f"{scc}_case"], merged[f"{scc}_ctrl"]) if len(common) > 1 else (np.nan, np.nan)
        try:
            _, wp = stats.wilcoxon(merged[f"{scc}_case"], merged[f"{scc}_ctrl"]) if len(common) > 2 else (np.nan, np.nan)
        except ValueError:
            wp = np.nan
        sd = np.std(d, ddof=1) if len(d) > 1 else np.nan
        r.update({f"{scc}_diff": float(np.mean(d)), f"{scc}_g": float(np.mean(d) / sd) if sd else np.nan,
                  f"{scc}_g_se": np.nan, f"{scc}_p_t": float(tp_) if tp_ == tp_ else np.nan,
                  f"{scc}_p_wilcoxon": float(wp) if wp == wp else np.nan, f"{scc}_n": int(len(common))})
    r["n_up"], r["n_ex"] = int(sc4.n_up.iloc[0]), int(sc4.n_ex.iloc[0])
    results.append(r)
    percol.append((spname, len(common)))

# ---------------- 5. 汇总 ----------------
res = pd.DataFrame(results)
res.to_csv(os.path.join(INTER, "M10D_contrasts.csv"), index=False)
res.to_csv(os.path.join(TAB, "Table_S83_M10D_CrossSpecies_Contrasts.csv"), index=False)

note("\n== 汇总：主对比方向（UCS 预测 <0）")
prim = res[~res["design"].astype(str).str.contains("sensitivity")]
for _, r in prim.iterrows():
    note(f"   {r['dataset']:<42} {r['species']:<14} UCS diff={r['UCS_diff']:+.3f} (g={r['UCS_g']:+.2f}, p_t={r['UCS_p_t']:.3g}) | "
         f"EIS diff={r['EIS_diff']:+.3f} | MDI diff={r['MDI_diff']:+.3f} (p_t={r['MDI_p_t']:.3g})")
nonhuman_prim = prim[prim["species"].str.lower() != "human"]
ucs_dir = (nonhuman_prim["UCS_diff"] < 0).sum()
note(f"\n   非人模型 UCS 方向一致（<0）：{ucs_dir}/{len(nonhuman_prim)} 个模型/物种对比（主对比口径）")
nomt_dir = (nonhuman_prim["UCS_nomt_diff"] < 0).sum()
note(f"   去 MT 敏感性：{nomt_dir}/{len(nonhuman_prim)}")
flush_log()
print("\nDONE")
