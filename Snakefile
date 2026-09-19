"""
scRNA-seq PBMC3k analysis pipeline
===================================
Reproducible single-cell RNA-seq workflow: QC -> HVG selection/clustering ->
marker-gene identification -> automated cell-type annotation.

Run:
    snakemake --cores 4 --use-conda   (or omit --use-conda if deps are
                                        already installed, e.g. in Docker)

Visualize the DAG:
    snakemake --dag | dot -Tpng > results/figures/dag.png
"""

configfile: "config.yaml"

rule all:
    input:
        "data/raw/pbmc3k_expression_matrix.h5ad",
        "results/figures/01_qc_violin.png",
        "results/tables/qc_summary.csv",
        "results/figures/02_umap_clusters.png",
        "results/tables/marker_genes_per_cluster.csv",
        "results/tables/cluster_cell_type_annotation.csv",
        "results/figures/03_umap_celltypes.png",
        "results/figures/04_marker_dotplot.png",


rule prepare_data:
    output:
        "data/raw/pbmc3k_expression_matrix.h5ad"
    shell:
        "python scripts/prepare_data.py config.yaml"


rule qc_filter:
    input:
        "data/raw/pbmc3k_expression_matrix.h5ad"
    output:
        h5ad="data/processed/01_qc_filtered.h5ad",
        fig="results/figures/01_qc_violin.png",
        summary="results/tables/qc_summary.csv"
    shell:
        "python scripts/qc_filter.py config.yaml {input} {output.h5ad} {output.fig} {output.summary}"


rule cluster:
    input:
        "data/processed/01_qc_filtered.h5ad"
    output:
        h5ad="data/processed/02_clustered.h5ad",
        fig="results/figures/02_umap_clusters.png"
    shell:
        "python scripts/cluster.py config.yaml {input} {output.h5ad} {output.fig}"


rule annotate_celltypes:
    input:
        "data/processed/02_clustered.h5ad"
    output:
        h5ad="data/processed/03_annotated.h5ad",
        markers="results/tables/marker_genes_per_cluster.csv",
        annotation="results/tables/cluster_cell_type_annotation.csv",
        umap_fig="results/figures/03_umap_celltypes.png",
        dotplot_fig="results/figures/04_marker_dotplot.png"
    shell:
        "python scripts/annotate_celltypes.py config.yaml {input} {output.h5ad} "
        "{output.markers} {output.annotation} {output.umap_fig} {output.dotplot_fig}"
