# Load GSE152418 raw counts and cohort metadata into a DESeqDataSet.
#
# GSE152418 (Arunachalam et al., Science 2020) is bulk RNA-seq of PBMCs from a
# COVID-19 cohort in Atlanta, GA. It is used here rather than a Bioconductor
# example dataset because it ships genuine raw integer counts on GEO, so the
# whole analysis runs from primary data.

load_counts <- function(path = file.path(DIRS$data, "counts.txt.gz")) {
  stopifnot(file.exists(path))
  m <- read.delim(gzfile(path), row.names = 1, check.names = FALSE)
  as.matrix(m)
}

#' Parse the GEO series matrix into one row per sample.
#'
#' The series matrix stores metadata transposed: each "!Sample_*" line holds one
#' field across all samples. Characteristics lines are "key: value" pairs whose
#' key is read from the first sample rather than assumed, so the parser does not
#' silently mislabel a column if GEO reorders the fields.
load_metadata <- function(path = file.path(DIRS$data, "series.txt.gz")) {
  stopifnot(file.exists(path))
  lines <- readLines(gzfile(path), warn = FALSE)

  field <- function(prefix) {
    hit <- lines[startsWith(lines, prefix)]
    lapply(hit, function(l) gsub('"', "", strsplit(l, "\t")[[1]][-1], fixed = TRUE))
  }

  titles <- field("!Sample_title")[[1]]
  acc    <- field("!Sample_geo_accession")[[1]]
  meta   <- data.frame(sample = titles, geo_accession = acc,
                       stringsAsFactors = FALSE)

  for (row in field("!Sample_characteristics_ch1")) {
    key <- sub(":.*", "", row[1])
    key <- gsub("[^a-z0-9]+", "_", tolower(trimws(key)))
    meta[[key]] <- trimws(sub("^[^:]*:\\s*", "", row))
  }
  meta
}

#' Assemble counts + metadata into a DESeqDataSet ready for testing.
#'
#' Design is ~ gender + condition. Sex is included because PBMC composition and
#' several interferon-stimulated genes differ by sex; leaving it in the residual
#' variance costs power and can bias the fold changes if the groups are not
#' sex-balanced. Condition is placed last so it is the coefficient that
#' results() reports by default.
build_dds <- function(counts, meta, drop_states = "Convalescent") {
  stopifnot(setequal(colnames(counts), meta$sample))
  meta <- meta[match(colnames(counts), meta$sample), ]
  rownames(meta) <- meta$sample

  # The single convalescent sample is neither acute infection nor healthy
  # baseline. With n = 1 it cannot form its own group, and folding it into
  # either arm would blur the contrast, so it is excluded and recorded.
  drop <- meta$disease_state %in% drop_states
  if (any(drop)) {
    log_line("excluding ", sum(drop), " sample(s) with disease state in {",
             paste(drop_states, collapse = ", "), "}")
    record("samples_excluded", sum(drop))
    meta   <- meta[!drop, , drop = FALSE]
    counts <- counts[, rownames(meta), drop = FALSE]
  }

  meta$condition <- factor(ifelse(meta$disease_state == "Healthy",
                                  "Healthy", "COVID19"),
                           levels = c("Healthy", "COVID19"))
  meta$gender   <- factor(meta$gender)
  meta$severity <- factor(meta$severity)

  # Counts are integers on GEO but read.delim can widen them; DESeq2 requires
  # an integer matrix and errors out informatively if given doubles.
  storage.mode(counts) <- "integer"

  dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta,
                                design = ~ gender + condition)

  # Pre-filter only genes that are essentially unobserved. The threshold is
  # deliberately permissive: DESeq2's independent filtering already removes
  # low-count genes at the testing stage using an optimised cutoff, and
  # aggressive filtering here would discard that optimisation.
  keep <- rowSums(counts(dds) >= 10) >= min(table(meta$condition))
  log_line("pre-filter: ", sum(keep), " of ", nrow(dds), " genes retained")
  record("genes_before_filter", nrow(dds))
  record("genes_after_filter", sum(keep))
  dds[keep, ]
}
