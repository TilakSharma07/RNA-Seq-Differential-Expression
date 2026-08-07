# Reactome pathway over-representation analysis.
#
# Implemented directly against Reactome's NCBI2Reactome mapping rather than via
# clusterProfiler. The test is a one-sided hypergeometric (Fisher exact) test,
# which is exactly what over-representation analysis is; writing it out makes
# the two choices that actually matter visible instead of buried in defaults.

#' Over-representation of a gene list against a pathway collection.
#'
#' @param sig_entrez  Entrez IDs of the significant genes.
#' @param universe    Entrez IDs of every gene that was *tested*.
#'
#' The universe is the tested set, not the whole genome. This is the single
#' most consequential choice in over-representation analysis: using all human
#' genes as background silently treats genes that were filtered out for low
#' expression as "not enriched", which makes every immune pathway look
#' significant in any blood dataset.
enrich_reactome <- function(sig_entrez, universe, reactome,
                            min_overlap = 3, alpha = 0.05) {
  sig_entrez <- unique(as.character(sig_entrez[!is.na(sig_entrez)]))
  universe   <- unique(as.character(universe[!is.na(universe)]))
  sig_entrez <- intersect(sig_entrez, universe)

  # Restrict each pathway to genes that were actually measurable here, so
  # pathway sizes reflect this experiment rather than the full database.
  sets <- lapply(reactome$sets, function(g) intersect(as.character(g), universe))
  sets <- sets[lengths(sets) >= min_overlap]

  N <- length(universe)     # genes tested
  K <- length(sig_entrez)   # genes called significant

  rows <- lapply(names(sets), function(pid) {
    genes   <- sets[[pid]]
    overlap <- intersect(genes, sig_entrez)
    if (length(overlap) < min_overlap) return(NULL)
    # phyper upper tail: P(X >= observed), the probability of seeing at least
    # this many hits if the significant genes were drawn at random.
    p <- phyper(length(overlap) - 1, length(genes), N - length(genes), K,
                lower.tail = FALSE)
    data.frame(
      pathway_id   = pid,
      pathway_name = unname(reactome$names[pid]),
      set_size     = length(genes),
      overlap      = length(overlap),
      expected     = round(length(genes) * K / N, 2),
      fold_enrich  = round((length(overlap) / K) / (length(genes) / N), 2),
      pvalue       = p,
      genes        = paste(head(overlap, 25), collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, rows)
  if (is.null(res) || nrow(res) == 0) {
    log_line("no pathway met the minimum overlap")
    return(res)
  }
  res$padj <- p.adjust(res$pvalue, method = "BH")
  res <- res[order(res$pvalue), ]
  log_line(sum(res$padj < alpha), " pathways enriched at padj < ", alpha,
           " of ", nrow(res), " tested")
  record("pathways_tested", nrow(res))
  record("pathways_enriched", sum(res$padj < alpha))
  res
}

#' Translate Entrez IDs back to symbols in the reported gene column.
label_pathway_genes <- function(enr, annotation) {
  if (is.null(enr) || nrow(enr) == 0) return(enr)
  map <- setNames(annotation$symbol, annotation$entrez)
  enr$genes <- vapply(strsplit(enr$genes, ","), function(ids) {
    paste(ifelse(is.na(map[ids]), ids, map[ids]), collapse = ", ")
  }, character(1))
  enr
}

#' Horizontal bar chart of the most enriched pathways.
plot_enrichment <- function(enr, n = 15, alpha = 0.05) {
  if (is.null(enr) || nrow(enr) == 0) return(invisible(NULL))
  d <- head(enr[enr$padj < alpha, ], n)
  if (nrow(d) == 0) {
    log_line("nothing passed padj < ", alpha, "; skipping enrichment plot")
    return(invisible(NULL))
  }
  d$label <- ifelse(nchar(d$pathway_name) > 52,
                    paste0(substr(d$pathway_name, 1, 49), "..."),
                    d$pathway_name)
  d$label <- factor(d$label, levels = rev(d$label))

  p <- ggplot(d, aes(x = label, y = -log10(padj), fill = fold_enrich)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = sprintf("%d/%d", overlap, set_size)),
              hjust = -0.15, size = 2.7, colour = "grey25") +
    coord_flip() +
    scale_fill_gradient(low = "#9EC6E0", high = "#B03A3A",
                        name = "fold\nenrichment") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(x = NULL, y = expression(-log[10]~adjusted~italic(p)),
         title = "Reactome pathways over-represented in COVID-19",
         subtitle = "labels show significant genes / pathway size") +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.y = element_blank())
  save_plot(p, "06_reactome_enrichment.png", width = 8.4, height = 5.6)
}
