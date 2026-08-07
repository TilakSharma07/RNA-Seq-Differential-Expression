# Shared helpers: logging, paths, and annotation lookups.
#
# Annotation deliberately comes from flat primary-source files rather than
# Bioconductor annotation packages (org.Hs.eg.db / clusterProfiler). Those
# packages pin a release at install time and are a recurring source of
# "works on my machine" drift; the files below are versioned by download date
# in results/run_log.tsv, so a reader can reproduce the exact mapping used.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

DIRS <- list(
  data    = "data",
  figures = file.path("results", "figures"),
  tables  = file.path("results", "tables")
)

log_line <- function(...) {
  msg <- paste0(...)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  invisible(msg)
}

#' Append one key/value row to the run log.
#'
#' Every number quoted in the README is written here by the code that computed
#' it, so the README can be regenerated instead of hand-edited.
record <- function(key, value) {
  path <- file.path("results", "run_log.tsv")
  if (!file.exists(path)) {
    write.table(data.frame(key = "key", value = "value"), path,
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  write.table(data.frame(key = key, value = as.character(value)), path,
              sep = "\t", row.names = FALSE, col.names = FALSE,
              quote = FALSE, append = TRUE)
  invisible(NULL)
}

save_plot <- function(plot, filename, width = 7, height = 5, dpi = 200) {
  path <- file.path(DIRS$figures, filename)
  ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
  log_line("wrote ", path)
  invisible(path)
}

write_table <- function(df, filename) {
  path <- file.path(DIRS$tables, filename)
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
  log_line("wrote ", path, " (", nrow(df), " rows)")
  invisible(path)
}

#' Map Ensembl gene IDs to HGNC symbols and Entrez IDs.
#'
#' Parses NCBI's Homo_sapiens.gene_info: the dbXrefs column carries
#' "Ensembl:ENSGxxxxxxxxxxx" among pipe-separated cross-references. Genes with
#' no Ensembl cross-reference are simply absent from the result; callers keep
#' the Ensembl ID as the fallback label rather than dropping the gene, because
#' an unmapped gene can still be differentially expressed.
load_gene_annotation <- function(path = file.path(DIRS$data, "gene_info.gz")) {
  stopifnot(file.exists(path))
  gi <- read.delim(gzfile(path), stringsAsFactors = FALSE, quote = "",
                   colClasses = "character")
  colnames(gi)[1] <- "tax_id"
  gi <- gi[gi$tax_id == "9606", c("GeneID", "Symbol", "dbXrefs", "type_of_gene")]

  has_ens <- grepl("Ensembl:", gi$dbXrefs, fixed = TRUE)
  gi <- gi[has_ens, ]
  gi$ensembl <- sub(".*Ensembl:(ENSG[0-9]+).*", "\\1", gi$dbXrefs)

  ann <- data.frame(
    ensembl = gi$ensembl,
    symbol  = gi$Symbol,
    entrez  = gi$GeneID,
    biotype = gi$type_of_gene,
    stringsAsFactors = FALSE
  )
  # A handful of Ensembl IDs map to more than one Entrez record; keep the first
  # so the mapping stays one-to-one and joins cannot silently duplicate rows.
  ann[!duplicated(ann$ensembl), ]
}

#' Reactome pathway membership as a list of Entrez ID vectors.
#'
#' Reactome is a nested hierarchy and this file lists every level, so a gene
#' appears under its specific pathway and under every ancestor. Two consequences
#' are handled here. Size filtering drops the very broad top-level terms
#' ("Immune System" carries thousands of genes and is never informative) and the
#' very small ones where a single gene swings the test. Deduplication then
#' collapses parent/child pairs whose gene sets are identical after that
#' filtering: they are one hypothesis reported twice, and leaving both in place
#' inflates the apparent number of enriched pathways and makes the multiple
#' testing correction more conservative than it needs to be. The retained name
#' is the one whose pathway ID sorts first, so the choice is deterministic.
load_reactome <- function(path = file.path(DIRS$data, "ncbi2reactome.txt"),
                          min_size = 10, max_size = 500) {
  stopifnot(file.exists(path))
  rt <- read.delim(path, header = FALSE, stringsAsFactors = FALSE, quote = "",
                   col.names = c("entrez", "pathway_id", "url", "pathway_name",
                                 "evidence", "species"))
  rt <- rt[rt$species == "Homo sapiens", ]
  rt$pathway_name <- trimws(rt$pathway_name)
  rt <- unique(rt[, c("entrez", "pathway_id", "pathway_name")])

  sets  <- split(rt$entrez, rt$pathway_id)
  names_map <- rt$pathway_name[!duplicated(rt$pathway_id)]
  names(names_map) <- rt$pathway_id[!duplicated(rt$pathway_id)]

  sizes <- lengths(sets)
  sets  <- sets[sizes >= min_size & sizes <= max_size]

  # Collapse pathways whose gene sets are identical (nested parent/child terms).
  fingerprint <- vapply(sets, function(g) paste(sort(g), collapse = ","),
                        character(1))
  ord  <- order(names(sets))
  sets <- sets[ord][!duplicated(fingerprint[ord])]

  list(sets = sets, names = names_map[names(sets)])
}
