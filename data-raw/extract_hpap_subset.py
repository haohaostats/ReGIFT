"""Create the deterministic HPAP subset used by data/regift_example.rda.

Run this script from the package root before create_regift_example.R. The full
H5AD is development input only and is never included in the package.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


SOURCE = Path("data/external/hpap_fasolino/GSE148073_Fasolino_69645.h5ad")
OUT = Path("data-raw/hpap_subset_tmp")
SEED = 20260901
CELLS_PER_DONOR_STATE = 20
N_GENES = 240

T1D_DONORS = ["HPAP020", "HPAP021", "HPAP023", "HPAP028", "HPAP032"]
CONTROL_DONORS = ["HPAP022", "HPAP039", "HPAP042", "HPAP044", "HPAP047"]
STATES = ["alpha", "beta_major", "acinar", "duct_major"]
MARKERS = [
    "INS", "IAPP", "PCSK1", "PCSK2", "MAFA", "PDX1", "G6PC2", "SLC2A2",
    "GCG", "TTR", "LOXL4", "FAP", "SST", "PPY", "PRSS1", "PRSS2",
    "REG1A", "REG1B", "KRT7", "KRT8", "KRT18", "KRT19", "MUC1",
    "COL1A1", "COL3A1", "VIM", "PECAM1", "KDR", "HLA-DPA1", "HLA-DRB1",
    "CD74", "IFITM3", "IL32", "LYZ", "C1QA", "LST1", "ZFP36",
]


def main() -> None:
    if not SOURCE.exists():
        raise FileNotFoundError(f"HPAP source file not found: {SOURCE}")
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    adata = ad.read_h5ad(SOURCE, backed="r")
    obs = adata.obs.copy()

    donors = T1D_DONORS + CONTROL_DONORS
    selected = []
    for donor in donors:
        for state in STATES:
            idx = np.flatnonzero(
                (obs["donor_id"].to_numpy() == donor)
                & (obs["cell_label"].to_numpy() == state)
            )
            if idx.size:
                selected.extend(rng.choice(
                    idx, size=min(CELLS_PER_DONOR_STATE, idx.size), replace=False
                ).tolist())
    selected = np.asarray(selected, dtype=int)

    if adata.raw is None:
        raise ValueError("HPAP H5AD does not contain a raw count matrix.")
    var = adata.raw.var.copy()
    symbols = var["feature_name"].astype(str)
    marker_idx = [i for i, gene in enumerate(symbols) if gene in MARKERS]
    score = pd.to_numeric(var["vst.variance.standardized"], errors="coerce")
    eligible = (
        var["feature_type"].astype(str).eq("protein_coding")
        & ~symbols.str.startswith(("MT-", "RPL", "RPS"))
    )
    ranked = np.flatnonzero(eligible.to_numpy())[
        np.argsort(-score[eligible].fillna(-np.inf).to_numpy(), kind="stable")
    ]
    gene_idx = []
    for i in marker_idx + ranked.tolist():
        if i not in gene_idx:
            gene_idx.append(i)
        if len(gene_idx) == N_GENES:
            break
    gene_idx = np.asarray(gene_idx, dtype=int)

    counts = adata.raw.X[selected, :][:, gene_idx]
    counts = counts.to_memory() if hasattr(counts, "to_memory") else counts
    counts = sparse.csr_matrix(counts)
    if counts.data.size and np.max(np.abs(counts.data - np.rint(counts.data))) > 1e-6:
        raise ValueError("Selected HPAP matrix is not raw integer count data.")
    counts.data = np.rint(counts.data).astype(np.int32)
    counts.eliminate_zeros()

    meta = obs.iloc[selected].copy()
    meta_out = pd.DataFrame({
        "cell": [f"{d}_{barcode}" for d, barcode in zip(
            meta["donor_id"].astype(str), meta.index.astype(str)
        )],
        "donor": meta["donor_id"].astype(str).to_numpy(),
        "sample": meta["donor_id"].astype(str).to_numpy(),
        "condition": np.where(meta["disease_state"].astype(str) == "Control",
                              "control", "T1D"),
        "state": meta["cell_label"].astype(str).to_numpy(),
        "sex": meta["sex"].astype(str).to_numpy(),
        "source_cell": meta.index.astype(str),
    })
    genes = pd.DataFrame({
        "gene": symbols.iloc[gene_idx].to_numpy(),
        "ensembl_id": var.index[gene_idx].astype(str),
    })

    mmwrite(OUT / "counts.mtx", counts)
    meta_out.to_csv(OUT / "meta.csv", index=False)
    genes.to_csv(OUT / "genes.csv", index=False)
    print(f"Wrote {counts.shape[0]} cells x {counts.shape[1]} genes to {OUT}")


if __name__ == "__main__":
    main()
