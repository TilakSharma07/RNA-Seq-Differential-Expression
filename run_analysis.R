#!/usr/bin/env Rscript
# Driver: raw GEO counts -> QC -> differential expression -> pathway enrichment.
#
# Usage:  Rscript run_analysis.R
# Inputs: data/ (fetched by scripts/download_data.sh)
# Outputs: results/figures/*.png, results/tables/*.tsv, results/run_log.tsv
#
# Every number quoted in README.md is written to results/run_log.tsv by the
# step that computed it, so the README is regenerated from the run rather than
# typed by hand.

suppressPackageStartupMessages(library(DESeq2))

source("R/utils.R")
source("R/01_load_data.R")
source("R/02_qc.R")
source("R/03_differential_expression.R")
source("R/04_enrichment.R")

set.seed(1)
dir.create(DIRS$figures, recursive = TRUE, showWarnings = FALSE)
dir.create(DIRS$tables,  recursive = TRUE, showWarnings = FALSE)
unlink(file.path("results", "run_log.tsv"))

ALPHA <- 0.05
LFC   <- 1

log_line("=== 1. load ===")
counts <- load_counts()
meta   <- load_metadata()
dds    <- build_dds(counts, meta)
record("samples_analysed", ncol(dds))
record("samples_covid",   sum(colData(dds)$condition == "COVID19"))
record("samples_healthy", sum(colData(dds)$condition == "Healthy"))
log_line(ncol(dds), " samples: ",
         sum(colData(dds)$condition == "COVID19"), " COVID-19, ",
         sum(colData(dds)$condition == "Healthy"), " healthy")

log_line("=== 2. quality control ===")
vsd <- run_qc(dds)

log_line("=== 3. differential expression ===")
de <- run_de(dds, alpha = ALPHA, lfc_threshold = LFC)
plot_volcano(de)
plot_top_heatmap(de, vsd)
write_table(de$table, "deseq2_all_genes.tsv")
write_table(de$significant, "deseq2_significant.tsv")

# Positive control on contrast direction. Interferon-stimulated genes are the
# best-replicated signal in COVID-19 whole blood and PBMCs; if IFI27 is not up
# in the COVID arm, the contrast has been specified backwards. This is asserted
# rather than printed so a silently inverted result cannot reach the README.
ctrl <- de$table[de$table$symbol == "IFI27", ]
stopifnot(nrow(ctrl) == 1, !is.na(ctrl$padj))
log_line("positive control IFI27: log2FC = ", round(ctrl$log2FoldChange, 2),
         ", padj = ", signif(ctrl$padj, 3))
record("control_IFI27_log2fc", round(ctrl$log2FoldChange, 3))
record("control_IFI27_padj", signif(ctrl$padj, 3))
stopifnot(ctrl$log2FoldChange > 0)

log_line("=== 4. pathway enrichment ===")
reactome <- load_reactome()
log_line(length(reactome$sets), " Reactome pathways after size filtering")
ann <- load_gene_annotation()
enr <- enrich_reactome(sig_entrez = de$significant$entrez,
                       universe   = de$table$entrez[!is.na(de$table$padj)],
                       reactome   = reactome, alpha = ALPHA)
enr <- label_pathway_genes(enr, ann)
if (!is.null(enr) && nrow(enr) > 0) {
  write_table(enr, "reactome_enrichment.tsv")
  plot_enrichment(enr, alpha = ALPHA)
}

record("r_version", paste(R.version$major, R.version$minor, sep = "."))
record("deseq2_version", as.character(packageVersion("DESeq2")))
record("run_date", format(Sys.Date()))
log_line("=== done; see results/run_log.tsv ===")
