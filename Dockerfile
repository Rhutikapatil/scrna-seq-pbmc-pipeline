FROM python:3.11-slim

# System libs needed by scanpy/igraph/leidenalg's compiled extensions
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /pipeline

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# Reproduce the full pipeline: fetch data -> QC -> cluster -> annotate
CMD ["snakemake", "--cores", "4", "--printshellcmds"]
