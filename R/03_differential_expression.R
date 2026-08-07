# Differential expression: COVID-19 PBMCs versus healthy controls.

suppressPackageStartupMessages(library(ggrepel))

#' Fit the negative-binomial model and extract the condition contrast.
#'
#' Two result objects are produced deliberately. `res` carries the untouched
#' test statistics and p-values used for calling significance. `shrunk` carries
#' apeglm-shrunken log fold changes, used only for ranking and plotting: raw
#' fold changes for low-count genes are wildly noisy, so an unshrunken ranking
#' puts near-unmeasurable genes at the top of the list.
run_de <- function(dds, alpha = 0.05, lfc_threshold = 1) {
  dds <- DESeq(dds, quiet = TRUE)
  log_line("model coefficients: ", paste(resultsNames(dds), collapse = ", "))

  res <- results(dds, contrast = c("condition", "COVID19", "Healthy"),
                 alpha = alpha)
  shrunk <- lfcShrink(dds, coef = "condition_COVID19_vs_Healthy",
                      type = "apeglm", quiet = TRUE)

  ann <- load_gene_annotation()
  out <- data.frame(
    ensembl        = rownames(res),
    baseMean       = res$baseMean,
    log2FoldChange = res$log2FoldChange,
    lfcShrink      = shrunk$log2FoldChange[match(rownames(res), rownames(shrunk))],
    stat           = res$stat,
    pvalue         = res$pvalue,
    padj           = res$padj,
    stringsAsFactors = FALSE
  )
  idx <- match(out$ensembl, ann$ensembl)
  out$symbol <- ifelse(is.na(idx), out$ensembl, ann$symbol[idx])
  out$entrez <- ann$entrez[idx]
  # Lead with the identifiers a reader scans first, not the internal ones.
  out <- out[, c("symbol", "ensembl", "entrez", "baseMean", "log2FoldChange",
                 "lfcShrink", "stat", "pvalue", "padj")]
  out <- out[order(out$padj, -abs(out$lfcShrink)), ]

  sig <- !is.na(out$padj) & out$padj < alpha & abs(out$log2FoldChange) >= lfc_threshold
  log_line(sum(sig), " genes at padj < ", alpha, " and |log2FC| >= ", lfc_threshold)

  record("de_alpha", alpha)
  record("de_lfc_threshold", lfc_threshold)
  record("de_genes_tested", sum(!is.na(out$padj)))
  record("de_significant", sum(sig))
  record("de_significant_up", sum(sig & out$log2FoldChange > 0))
  record("de_significant_down", sum(sig & out$log2FoldChange < 0))

  list(dds = dds, table = out, significant = out[sig, ], alpha = alpha,
       lfc_threshold = lfc_threshold)
}

#' Volcano plot with the strongest genes labelled.
plot_volcano <- function(de) {
  d <- de$table[!is.na(de$table$padj), ]
  d$status <- "not significant"
  up   <- d$padj < de$alpha & d$log2FoldChange >=  de$lfc_threshold
  down <- d$padj < de$alpha & d$log2FoldChange <= -de$lfc_threshold
  d$status[up]   <- "up in COVID-19"
  d$status[down] <- "down in COVID-19"

  # Cap the y axis: a handful of genes have p-values below double precision,
  # which would compress every other point into a strip at the bottom.
  cap <- 50
  d$neglog10p <- pmin(-log10(d$padj), cap)
  d$capped <- -log10(d$padj) > cap

  lab <- d[d$status != "not significant", ]
  lab <- lab[order(-abs(lab$lfcShrink)), ][seq_len(min(18, nrow(lab))), ]

  p <- ggplot(d, aes(log2FoldChange, neglog10p, colour = status)) +
    geom_point(alpha = 0.55, size = 1.1) +
    geom_vline(xintercept = c(-de$lfc_threshold, de$lfc_threshold),
               linetype = "dashed", colour = "grey45", linewidth = 0.3) +
    geom_hline(yintercept = -log10(de$alpha),
               linetype = "dashed", colour = "grey45", linewidth = 0.3) +
    geom_text_repel(data = lab, aes(label = symbol), size = 2.6,
                    max.overlaps = 30, show.legend = FALSE) +
    scale_colour_manual(values = c("up in COVID-19" = "#C44E52",
                                   "down in COVID-19" = "#4C72B0",
                                   "not significant" = "grey78")) +
    labs(x = "log2 fold change (COVID-19 / Healthy)",
         y = expression(-log[10]~adjusted~italic(p)),
         colour = NULL,
         title = "Differential expression in COVID-19 PBMCs",
         subtitle = sprintf("%d genes at padj < %.2f and |log2FC| >= %g",
                            nrow(de$significant), de$alpha, de$lfc_threshold)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  save_plot(p, "04_volcano.png", width = 7.5, height = 6)
}

#' Heatmap of the top genes by adjusted p-value.
#'
#' Values are z-scored per gene so that co-regulation is visible; without
#' scaling the plot is dominated by a few very highly expressed transcripts.
plot_top_heatmap <- function(de, vsd, n = 40) {
  top <- head(de$significant$ensembl, n)
  mat <- assay(vsd)[top, , drop = FALSE]
  mat <- t(scale(t(mat)))
  rownames(mat) <- de$significant$symbol[match(top, de$significant$ensembl)]

  ann <- data.frame(condition = colData(vsd)$condition,
                    severity  = colData(vsd)$severity,
                    row.names = colnames(vsd))
  png(file.path(DIRS$figures, "05_top_genes_heatmap.png"),
      width = 8, height = 7.5, units = "in", res = 200)
  pheatmap(mat, annotation_col = ann, show_colnames = FALSE,
           fontsize_row = 7, fontsize = 8,
           col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
           main = sprintf("Top %d differentially expressed genes (z-scored)", length(top)))
  dev.off()
  log_line("wrote ", file.path(DIRS$figures, "05_top_genes_heatmap.png"))
}
