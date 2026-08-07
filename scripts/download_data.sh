#!/usr/bin/env bash
# Fetch every input from its primary source. ~100 MB total.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data && cd data

GEO=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152418

echo "==> GSE152418 raw counts"
[ -s counts.txt.gz ] || curl -sL --fail -o counts.txt.gz \
  "$GEO/suppl/GSE152418_p20047_Study1_RawCounts.txt.gz"

echo "==> GSE152418 sample metadata"
[ -s series.txt.gz ] || curl -sL --fail -o series.txt.gz \
  "$GEO/matrix/GSE152418_series_matrix.txt.gz"

echo "==> NCBI gene_info (Ensembl -> symbol -> Entrez)"
[ -s gene_info.gz ] || curl -sL --fail -o gene_info.gz \
  "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"

echo "==> Reactome pathway membership"
[ -s ncbi2reactome.txt ] || curl -sL --fail -o ncbi2reactome.txt \
  "https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt"

ls -lh
