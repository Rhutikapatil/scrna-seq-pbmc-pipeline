"""
Stage 3 — Marker-gene identification and cell-type annotation.

Runs a differential-expression test (default: Wilcoxon rank-sum) of each
Leiden cluster against all others to find marker genes, then scores each
cluster against a canonical PBMC marker-gene dictionary (config.yaml) to
assign an automated cell-type label. Outputs marker tables, an annotated
UMAP, and a dot plot of marker expression by cell type.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def score_clusters(adata, marker_dict: dict) -> dict:
    """Assign each Leiden cluster the cell type whose marker genes score highest
    (mean scaled expression of that type's markers, averaged over cells in the cluster)."""
    scores = pd.DataFrame(index=sorted(adata.obs["leiden"].unique(), key=int))
    for cell_type, genes in marker_dict.items():
        genes_present = [g for g in genes if g in adata.raw.var_names]
        if not genes_present:
            scores[cell_type] = np.nan
            continue
        sc.tl.score_genes(adata, gene_list=genes_present, score_name="_tmp_score", use_raw=True)
        scores[cell_type] = adata.obs.groupby("leiden", observed=True)["_tmp_score"].mean()
    if "_tmp_score" in adata.obs:
        del adata.obs["_tmp_score"]
    assignment = scores.idxmax(axis=1).to_dict()
    return assignment, scores


def main(config_path: str, input_path: str, output_h5ad: str, markers_csv: str,
         annotation_csv: str, umap_fig: str, dotplot_fig: str) -> None:
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    adata = sc.read_h5ad(input_path)

    mg_cfg = cfg["marker_genes"]
    sc.tl.rank_genes_groups(
        adata, groupby="leiden", method=mg_cfg["method"], use_raw=True
    )

    n_top = mg_cfg["n_genes_per_cluster"]
    groups = adata.uns["rank_genes_groups"]["names"].dtype.names
    rows = []
    for cluster in groups:
        names = adata.uns["rank_genes_groups"]["names"][cluster][:n_top]
        scores = adata.uns["rank_genes_groups"]["scores"][cluster][:n_top]
        pvals = adata.uns["rank_genes_groups"]["pvals_adj"][cluster][:n_top]
        lfc = adata.uns["rank_genes_groups"]["logfoldchanges"][cluster][:n_top]
        for rank, (gene, score, p, l) in enumerate(zip(names, scores, pvals, lfc), start=1):
            rows.append({
                "leiden_cluster": cluster, "rank": rank, "gene": gene,
                "score": score, "pval_adj": p, "log2fc": l,
            })
    markers_df = pd.DataFrame(rows)
    Path(markers_csv).parent.mkdir(parents=True, exist_ok=True)
    markers_df.to_csv(markers_csv, index=False)

    assignment, score_table = score_clusters(adata, cfg["cell_type_markers"])
    adata.obs["cell_type"] = adata.obs["leiden"].map(assignment).astype("category")

    n_cells = adata.obs["leiden"].value_counts()
    annotation_df = pd.DataFrame({
        "leiden_cluster": list(assignment.keys()),
        "assigned_cell_type": list(assignment.values()),
        "n_cells": [int(n_cells.get(c, 0)) for c in assignment.keys()],
        "top_3_markers": [
            ", ".join(markers_df[markers_df.leiden_cluster == c].gene.head(3)) for c in assignment.keys()
        ],
    }).sort_values("leiden_cluster", key=lambda s: s.astype(int))
    annotation_df.to_csv(annotation_csv, index=False)

    fig, ax = plt.subplots(figsize=(6.5, 5))
    sc.pl.umap(adata, color="cell_type", ax=ax, show=False, legend_loc="right margin", title="")
    ax.set_title("UMAP — automated cell-type annotation")
    fig.tight_layout()
    Path(umap_fig).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(umap_fig, dpi=150)
    plt.close(fig)

    marker_genes_flat = sorted({g for genes in cfg["cell_type_markers"].values() for g in genes
                                 if g in adata.raw.var_names})
    fig2 = sc.pl.dotplot(
        adata, var_names=marker_genes_flat, groupby="cell_type", use_raw=True,
        show=False, return_fig=True,
    )
    fig2.savefig(dotplot_fig, dpi=150)

    Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_h5ad)

    print("Cluster -> cell type assignment:")
    print(annotation_df.to_string(index=False))


if __name__ == "__main__":
    main(*sys.argv[1:8])
