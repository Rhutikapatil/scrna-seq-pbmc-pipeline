"""
Stage 2 — HVG selection, dimensionality reduction, and clustering.

Selects highly variable genes, scales the data, runs PCA, builds a
k-nearest-neighbor graph, clusters cells with the Leiden algorithm, and
computes a UMAP embedding for visualization.
"""

import sys
from pathlib import Path

import scanpy as sc
import yaml
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main(config_path: str, input_path: str, output_h5ad: str, output_fig: str) -> None:
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    adata = sc.read_h5ad(input_path)
    adata.raw = adata  # keep full log-normalized matrix for marker-gene testing later

    hvg_cfg = cfg["hvg"]
    sc.pp.highly_variable_genes(
        adata, n_top_genes=hvg_cfg["n_top_genes"], flavor=hvg_cfg["flavor"]
    )
    adata = adata[:, adata.var["highly_variable"]].copy()

    sc.pp.scale(adata, max_value=10)

    cl_cfg = cfg["clustering"]
    sc.tl.pca(adata, n_comps=cl_cfg["n_pcs"], random_state=cl_cfg["random_state"])
    sc.pp.neighbors(adata, n_neighbors=cl_cfg["n_neighbors"], n_pcs=cl_cfg["n_pcs"])
    sc.tl.leiden(adata, resolution=cl_cfg["leiden_resolution"], random_state=cl_cfg["random_state"])
    sc.tl.umap(adata, random_state=cl_cfg["random_state"])

    n_clusters = adata.obs["leiden"].nunique()
    print(f"Leiden clustering (resolution={cl_cfg['leiden_resolution']}): {n_clusters} clusters")

    fig, ax = plt.subplots(figsize=(6, 5))
    sc.pl.umap(adata, color="leiden", ax=ax, show=False, legend_loc="on data", title="")
    ax.set_title(f"UMAP — Leiden clusters (n={n_clusters})")
    fig.tight_layout()
    Path(output_fig).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_fig, dpi=150)
    plt.close(fig)

    Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_h5ad)


if __name__ == "__main__":
    main(*sys.argv[1:5])
