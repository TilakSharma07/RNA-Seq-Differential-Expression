# Quality control on the count matrix and the sample-level structure.
#
# These plots exist to catch problems that would invalidate the differential
# expression result: a mislabelled sample, a batch that dominates the biology,
# or an outlier driving the fold changes.

suppressPackageStartupMessages({
  library(pheatmap)
  library(RColorBrewer)
})

#' Variance-stabilising transform for visualisation only.
#'
#' VST is used for PCA/clustering because raw and log-CPM counts have a strong
#' mean-variance dependence that makes distances dominated by highly expressed
#' genes. Testing still runs on raw counts through the negative-binomial model;
#' the transformed matrix is never fed back into the differential test.
run_qc <- function(dds) {
  vsd <- vst(dds, blind = TRUE)  # blind: do not let the design bias the transform

  lib <- data.frame(sample = colnames(dds),
                    million_reads = colSums(counts(dds)) / 1e6,
                    condition = colData(dds)$condition)
  p_lib <- ggplot(lib, aes(x = reorder(sample, million_reads),
                           y = million_reads, fill = condition)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c(Healthy = "#4C72B0", COVID19 = "#C44E52")) +
    labs(x = NULL, y = "Assigned reads (millions)",
         title = "Library size by sample") +
    theme_bw(base_size = 9)
  save_plot(p_lib, "01_library_sizes.png", width = 7, height = 6)
  record("median_library_size_millions", round(median(lib$million_reads), 2))

  pca <- plotPCA(vsd, intgroup = c("condition", "severity"), returnData = TRUE)
  pct <- round(100 * attr(pca, "percentVar"))
  p_pca <- ggplot(pca, aes(PC1, PC2, colour = condition, shape = severity)) +
    geom_point(size = 3, alpha = 0.85) +
    scale_colour_manual(values = c(Healthy = "#4C72B0", COVID19 = "#C44E52")) +
    labs(x = sprintf("PC1: %d%% variance", pct[1]),
         y = sprintf("PC2: %d%% variance", pct[2]),
         title = "Sample clustering (variance-stabilised counts)") +
    theme_bw(base_size = 11)
  save_plot(p_pca, "02_pca.png", width = 7, height = 5)
  record("pc1_variance_pct", pct[1])
  record("pc2_variance_pct", pct[2])

  # Sample-to-sample distances: an off-diagonal block that ignores condition is
  # the signature of a batch effect or a swapped label.
  d <- dist(t(assay(vsd)))
  mat <- as.matrix(d)
  ann <- data.frame(condition = colData(dds)$condition,
                    row.names = colnames(dds))
  png(file.path(DIRS$figures, "03_sample_distances.png"),
      width = 7, height = 6.2, units = "in", res = 200)
  pheatmap(mat, clustering_distance_rows = d, clustering_distance_cols = d,
           annotation_col = ann, fontsize = 7,
           col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
           main = "Sample-to-sample distance")
  dev.off()
  log_line("wrote ", file.path(DIRS$figures, "03_sample_distances.png"))

  vsd
}
