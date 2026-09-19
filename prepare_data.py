"""
Stage 0 — Data acquisition.

Pulls the canonical 10x Genomics PBMC3k dataset (2,638 PBMCs from a healthy
donor) via the public scanpy/cellxgene community mirror and writes it to
data/raw/ as the pipeline's starting point.

Provenance note: the upstream file already carries 10x's standard per-cell
library-size normalization and a log1p transform (the raw UMI count matrix
is not separable from it). Everything downstream of this stage — QC
filtering, highly-variable-gene selection, scaling, PCA, clustering,
marker-gene ranking and cell-type annotation — is computed independently
in this pipeline, not reused from the source file.
"""

import sys
from pathlib import Path

import scanpy as sc
import yaml


def main(config_path: str) -> None:
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    raw_path = Path(cfg["data"]["raw_h5ad"])
    raw_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Fetching source dataset from {cfg['data']['source_url']} ...")
    source = sc.read(
        "data/raw/_pbmc3k_source.h5ad",
        backup_url=cfg["data"]["source_url"],
    )

    # The full, un-subset log-normalized expression matrix (all measured
    # genes) lives in .raw on this file; that is our analysis starting point
    # so we are not inheriting the source's own HVG/scaling/clustering choices.
    if source.raw is None:
        adata = source
    else:
        adata = source.raw.to_adata()
        adata.obs = source.obs[[c for c in ("n_genes", "percent_mito", "n_counts") if c in source.obs]].copy()

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    adata.write_h5ad(raw_path)
    print(f"Saved starting matrix: {adata.shape[0]} cells x {adata.shape[1]} genes -> {raw_path}")


if __name__ == "__main__":
    config_path = sys.argv[1] if len(sys.argv) > 1 else "config.yaml"
    main(config_path)
