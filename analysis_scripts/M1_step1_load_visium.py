# -*- coding: utf-8 -*-
"""
M1 Step 1: Load + QC GSE271370 Visium (23 sections) into merged AnnData.

GSE271370: "Spatial transcriptomics unveils the in situ cellular and molecular
hallmarks of the lung in fatal COVID-19" (FFPE Visium, GRCh38, HiSeq X Ten).
4 Control / 7 Acute DAD / 12 Proliferative DAD.  Metadata verified from GEO
(2026-08-17), see N4_资源核实报告.md.
"""
import os, gzip, csv
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import anndata as ad

RAW = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\00_RAW_DATA\GSE271370_Lung_COVID19_Visium\GSE271370_RAW"
OUT = r"E:\SCI\SCI论文1黄裕荣_Mitoxyperilysis_ARDS\_intermediate"
os.makedirs(OUT, exist_ok=True)

SAMPLE_MAP = {
    "GSM8375640": ("L2P", "ProliferativeDAD"),
    "GSM8375641": ("L19P", "ProliferativeDAD"),
    "GSM8375642": ("L11P", "ProliferativeDAD"),
    "GSM8375643": ("CONTROL2", "Control"),
    "GSM8387358": ("HRC5", "ProliferativeDAD"),
    "GSM8387359": ("HRC6", "ProliferativeDAD"),
    "GSM8387360": ("HRC8", "ProliferativeDAD"),
    "GSM8387361": ("HRC10", "ProliferativeDAD"),
    "GSM8387362": ("HRC11", "ProliferativeDAD"),
    "GSM8387363": ("HRC12", "ProliferativeDAD"),
    "GSM8387364": ("HRC13", "ProliferativeDAD"),
    "GSM8387365": ("HRC16", "ProliferativeDAD"),
    "GSM8387366": ("HRC17", "ProliferativeDAD"),
    "GSM8387367": ("L5P", "AcuteDAD"),
    "GSM8387368": ("L14P", "AcuteDAD"),
    "GSM8387369": ("L24P", "AcuteDAD"),
    "GSM8387370": ("L12P", "AcuteDAD"),
    "GSM8387371": ("HRC2", "AcuteDAD"),
    "GSM8387372": ("HRC4", "AcuteDAD"),
    "GSM8387373": ("HRC18", "AcuteDAD"),
    "GSM8387374": ("L3C", "Control"),
    "GSM8387375": ("L14C", "Control"),
    "GSM8387376": ("L2C", "Control"),
}
assert len(SAMPLE_MAP) == 23

KINDS = ["filtered_matrix.mtx", "filtered_features.tsv",
         "filtered_barcodes.tsv", "tissue_positions.csv"]

def find_files():
    """Map GSM id -> dict kind->path (handles .gz and section title inside filename)."""
    by_gsm = {}
    for fn in os.listdir(RAW):
        base = fn[:-3] if fn.endswith(".gz") else fn
        gsm = base.split("_")[0]
        for kind in KINDS:
            if base.endswith(kind):
                by_gsm.setdefault(gsm, {})[kind] = os.path.join(RAW, fn)
    return by_gsm

def load_section(gsm, files):
    """Load one Visium section: (counts CSR, genes, barcodes, coords, in_tissue)."""
    def open_maybe(kind):
        path = files[kind]
        if path.endswith(".gz"):
            return gzip.open(path, "rt")
        return open(path, "rt")

    mtx_name = "filtered_matrix.mtx"
    fea_name = "filtered_features.tsv"
    bar_name = "filtered_barcodes.tsv"
    pos_name = "tissue_positions.csv"

    with open_maybe(mtx_name) as f:
        m = sio.mmread(f)
    # Space Ranger 2.0 files here are stored genes x spots -> transpose
    m = sp.csr_matrix(m).T.tocsr()
    with open_maybe(fea_name) as f:
        fea = [line.rstrip("\n").split("\t") for line in f]
    genes = [f[1] for f in fea]
    with open_maybe(bar_name) as f:
        bc = [line.rstrip("\n") for line in f]
    with open_maybe(pos_name) as f:
        rd = csv.reader(f)
        header = next(rd)
        rows = [r for r in rd]
    # rows: barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres
    pos_df = pd.DataFrame(rows, columns=header)
    pos_df["in_tissue"] = pos_df["in_tissue"].astype(int)
    pos_df["pxl_row_in_fullres"] = pos_df["pxl_row_in_fullres"].astype(float)
    pos_df["pxl_col_in_fullres"] = pos_df["pxl_col_in_fullres"].astype(float)
    # align
    idx = {b: i for i, b in enumerate(bc)}
    pos_df["i"] = pos_df["barcode"].map(idx)
    pos_df = pos_df.dropna(subset=["i"]).astype({"i": int}).sort_values("i")
    return m, genes, bc, pos_df

def main():
    by_gsm = find_files()
    qc_rows = []
    mats, obs_all, obsms = [], [], []
    var_names = None
    for gsm, (title, cond) in SAMPLE_MAP.items():
        files = by_gsm.get(gsm)
        if files is None:
            raise FileNotFoundError(gsm)
        m, genes, bc, pos = load_section(gsm, files)
        n_raw = m.shape[0]
        keep = pos[pos["in_tissue"] == 1]
        m = m[keep["i"].values, :]
        n_in = m.shape[0]
        tot = np.asarray(m.sum(1)).ravel()
        ngenes = np.asarray((m > 0).sum(1)).ravel()
        mask = (tot >= 100) & (ngenes >= 50)
        m = m[mask, :]
        keep = keep.iloc[np.where(mask)[0]]
        n_qc = m.shape[0]
        # CPM + log1p
        lib = np.asarray(m.sum(1)).ravel()
        m = sp.csr_matrix(m.astype(np.float64).multiply(1.0 / np.maximum(lib, 1)[:, None]) * 1e4)
        m.data = np.log1p(m.data)
        mats.append(m)
        obs_all.append(pd.DataFrame({
            "section": title, "condition": cond, "gsm": gsm,
            "n_counts": lib, "n_genes": np.asarray((m > 0).sum(1)).ravel(),
        }))
        obsms.append(keep[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values.astype(np.float64))
        if var_names is None:
            var_names = genes
        else:
            assert var_names == genes, "feature order mismatch"
        qc_rows.append(dict(gsm=gsm, section=title, condition=cond,
                            n_spots_raw=n_raw, n_spots_in_tissue=n_in,
                            n_spots_qc=n_qc,
                            median_counts=int(np.median(lib)),
                            median_genes=int(np.median(ngenes[mask]))))
        print(f"{gsm} {title:12s} {cond:18s} raw={n_raw} in_tissue={n_in} qc={n_qc}")

    X = sp.vstack(mats)
    obs = pd.concat(obs_all, ignore_index=True)
    obs.index = [f"{r.section}_{i}" for i, r in obs.iterrows()]
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=var_names))
    # per-section spatial in adata.obsm keyed by section
    adata.obsm["spatial"] = np.vstack(obsms)
    adata.uns["M1"] = dict(dataset="GSE271370", normalization="CPM1e4_log1p",
                           qc="in_tissue & total_counts>=100 & n_genes>=50")
    adata.write_h5ad(os.path.join(OUT, "M1_visium_merged.h5ad"))

    qc_df = pd.DataFrame(qc_rows)
    qc_df.to_csv(os.path.join(OUT, "M1_visium_qc_summary.csv"), index=False)
    print("QC summary saved. Total spots:", adata.shape[0])

if __name__ == "__main__":
    main()
