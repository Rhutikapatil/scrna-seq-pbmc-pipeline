"""
Stage 1 — Quality control.

Computes standard single-cell QC metrics (genes detected per cell, total
counts, mitochondrial read percentage), plots their distributions before
filtering, applies cell/gene-level thresholds from config.yaml, and writes
a QC summary table plus the filtered AnnData object.
"""

import sys
from pathlib import Path

import scanpy as sc
import yaml
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main(config_path: str, input_path: str, output_h5ad: str, output_fig: str, output_summary: str) -> None:
    with open(config_path) as f:
        cfg = yaml.safe_load(f)["qc"]

    adata = sc.read_h5ad(input_path)
    n_cells_before, n_genes_before = adata.shape

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    sc.pl.violin(adata, "n_genes_by_counts", ax=axes[0], show=False)
    axes[0].set_title("Genes detected / cell")
    sc.pl.violin(adata, "total_counts", ax=axes[1], show=False)
    axes[1].set_title("Total counts / cell")
    sc.pl.violin(adata, "pct_counts_mt", ax=axes[2], show=False)
    axes[2].set_title("% mitochondrial counts")
    fig.suptitle("QC metric distributions (pre-filtering)")
    fig.tight_layout()
    Path(output_fig).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_fig, dpi=150)
    plt.close(fig)

    sc.pp.filter_cells(adata, min_genes=cfg["min_genes_per_cell"])
    sc.pp.filter_genes(adata, min_cells=cfg["min_cells_per_gene"])
    adata = adata[adata.obs["pct_counts_mt"] < cfg["max_pct_mt"], :]
    adata = adata[adata.obs["n_genes_by_counts"] < cfg["max_genes_per_cell"], :]
    adata = adata.copy()

    n_cells_after, n_genes_after = adata.shape

    Path(output_summary).parent.mkdir(parents=True, exist_ok=True)
    with open(output_summary, "w") as f:
        f.write("metric,before_filtering,after_filtering,pct_retained\n")
        f.write(f"cells,{n_cells_before},{n_cells_after},{100 * n_cells_after / n_cells_before:.1f}\n")
        f.write(f"genes,{n_genes_before},{n_genes_after},{100 * n_genes_after / n_genes_before:.1f}\n")

    Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_h5ad)

    print(f"QC filtering: {n_cells_before} -> {n_cells_after} cells, {n_genes_before} -> {n_genes_after} genes")


if __name__ == "__main__":
    main(*sys.argv[1:6])
