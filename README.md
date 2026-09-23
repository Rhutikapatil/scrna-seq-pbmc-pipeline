# scRNA-seq PBMC Analysis Pipeline
![Python](https://img.shields.io/badge/Python-3.x-blue?logo=python&logoColor=white)
![Scanpy](https://img.shields.io/badge/Scanpy-single--cell%20analysis-1f77b4)
![Snakemake](https://img.shields.io/badge/Snakemake-workflow-6A3D9A)
![Docker](https://img.shields.io/badge/Docker-containerized-2496ED?logo=docker&logoColor=white)
![UMAP](https://img.shields.io/badge/UMAP-dimensionality%20reduction-orange)
![Leiden](https://img.shields.io/badge/Leiden-clustering-2E8B57)

A reproducible single-cell RNA-seq workflow — QC → highly-variable-gene selection
→ Leiden clustering → marker-gene identification → automated cell-type annotation —
built as a **Snakemake** pipeline and containerized with **Docker**.

Built on the canonical **10x Genomics PBMC3k** dataset (2,638 peripheral blood
mononuclear cells from a healthy donor), a standard single-cell benchmark.

## Project Highlights

- Built a reproducible scRNA-seq workflow using Snakemake
- Processed the canonical PBMC3k single-cell dataset
- Performed cell-level and gene-level quality control
- Selected highly variable genes for downstream analysis
- Applied PCA, k-nearest-neighbor graph construction, Leiden clustering, and UMAP
- Identified cluster-specific marker genes using Wilcoxon rank-sum testing
- Implemented automated cell-type annotation using canonical PBMC marker genes
- Containerized the workflow with Docker for reproducibility
- Parameterized analysis thresholds through `config.yaml`
- Documented biological interpretation and known annotation limitations
## Why this project

Most scRNA-seq portfolio pieces are a single Jupyter notebook. This one is
structured the way a production bioinformatics pipeline actually gets built:
parameterized (`config.yaml`), staged into independent, cacheable rules
(`Snakefile`), containerized for reproducibility (`Dockerfile`), and versioned
outputs (tables + figures) rather than notebook cell output.

## Pipeline

```
prepare_data → qc_filter → cluster → annotate_celltypes
```

| Stage | What it does | Key outputs |
|---|---|---|
| `prepare_data` | Pulls the PBMC3k expression matrix (13,714 genes × 2,638 cells) | `data/raw/pbmc3k_expression_matrix.h5ad` |
| `qc_filter` | Computes per-cell QC metrics (genes detected, % mitochondrial reads), filters low-quality cells/genes | `results/figures/01_qc_violin.png`, `results/tables/qc_summary.csv` |
| `cluster` | HVG selection (top 2,000 genes) → scaling → PCA (30 PCs) → kNN graph → **Leiden clustering** → UMAP | `results/figures/02_umap_clusters.png` |
| `annotate_celltypes` | Wilcoxon marker-gene test per cluster + automated cell-type scoring against a canonical PBMC marker dictionary | `results/tables/marker_genes_per_cluster.csv`, `results/tables/cluster_cell_type_annotation.csv`, `results/figures/03_umap_celltypes.png`, `results/figures/04_marker_dotplot.png` |

All thresholds (QC cutoffs, HVG count, PCs, clustering resolution, marker
genes) live in `config.yaml` — re-running the pipeline with different
parameters requires no code changes.

## Data provenance

Direct access to GEO/10x Genomics servers wasn't available from the build
environment used for this project, so the expression matrix was sourced via
the public scanpy/cellxgene community mirror of the same PBMC3k dataset
used in the canonical Scanpy and Seurat tutorials. That source file carries
10x's standard per-cell library-size normalization and a log1p transform
already applied (raw UMI counts aren't recoverable from it). **Everything
downstream — QC metrics, filtering thresholds, HVG selection, scaling, PCA,
clustering, marker-gene ranking, and cell-type annotation — is computed
independently in this pipeline**, not reused from the source.

## Results

**QC:** 2,638 / 2,638 cells and 13,656 / 13,714 genes passed filtering
(`min_genes_per_cell=200`, `max_pct_mt=10`, `max_genes_per_cell=4000`) —
this is a pre-cleaned tutorial dataset, so retention is near-total; the
thresholds are the interesting part; they were chosen (and documented in
`config.yaml`) to reflect the ranges commonly used for 10x PBMC data.

**Clustering:** Leiden (resolution 0.6) found **8 clusters**, automatically
annotated against canonical lineage markers:

| Cluster | Cell type | Cells | Top markers |
|---|---|---|---|
| 0 | CD4 T cells | 1,021 | LDHB, RPS27, RPS12 |
| 1 | CD14+ Monocytes | 373 | CST3, COTL1, AIF1 |
| 2 | B cells | 345 | CD79A, CD74, CD79B |
| 3 | CD14+ Monocytes | 298 | S100A8, S100A9, LYZ |
| 4 | NK cells | 251 | CCL5, NKG7, GZMA |
| 5 | CD4 T cells | 171 | CCL5, CD3D, CD8B |
| 6 | NK cells | 166 | NKG7, GNLY, GZMB |
| 7 | Platelets | 13 | PPBP, PF4, GNG11 |

These labels line up with the expected PBMC composition (T/NK cells,
monocytes, B cells, and a small platelet cluster) and match the reference
cell-type calls published alongside this dataset.

**Known limitation:** cluster 5's top markers (CD3D, CD8B) indicate it's
actually a CD8+ cytotoxic T cell population, not CD4+ — the score-based
marker-dictionary approach doesn't cleanly separate CD4 vs. CD8 T cells at
this clustering resolution, since both lineages share most of the broad
T-cell markers used here. A finer-resolution re-cluster restricted to the
T-cell compartment, or reference-based label transfer (e.g. Azimuth,
CellTypist), would resolve this — noting it here rather than smoothing over
it because that kind of QC-on-your-own-results is the actual skill.

## Running it

**Locally:**
```bash
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
snakemake --cores 4 --printshellcmds
```

**With Docker:**
```bash
docker build -t pbmc-scrna-pipeline .
docker run -v $(pwd)/results:/pipeline/results pbmc-scrna-pipeline
```

**Re-parameterize** (e.g. a stricter mitochondrial cutoff or higher clustering
resolution) by editing `config.yaml`, then re-run — Snakemake only
recomputes the affected downstream stages.

## Project structure

```
.
├── Snakefile              # pipeline DAG (4 rules)
├── config.yaml             # all tunable parameters
├── requirements.txt
├── Dockerfile
├── scripts/
│   ├── prepare_data.py
│   ├── qc_filter.py
│   ├── cluster.py
│   └── annotate_celltypes.py
├── data/
│   ├── raw/                 # input expression matrix (gitignored)
│   └── processed/           # intermediate .h5ad per stage (gitignored)
└── results/
    ├── figures/              # QC violins, UMAPs, marker dot plot
    └── tables/               # QC summary, marker genes, cell-type calls
```

## Stack

Python · scanpy · anndata · Leiden clustering · Snakemake · Docker

## Author

Rhutika Patil — M.S. Bioinformatics, NC State University
[linkedin.com/in/rhutika-patil](https://linkedin.com/in/rhutika-patil) ·
[github.com/Rhutikapatil](https://github.com/Rhutikapatil)
