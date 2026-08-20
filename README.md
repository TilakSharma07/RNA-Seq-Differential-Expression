# Differential Expression in COVID-19 PBMCs

Bulk RNA-seq analysis of peripheral blood mononuclear cells from COVID-19
patients and healthy controls, built end-to-end from public raw counts with
DESeq2. Includes quality control, differential testing, and Reactome pathway
over-representation analysis.

**Dataset** — [GSE152418](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152418)
(Arunachalam et al., *Science* 2020): 33 PBMC samples,
16 COVID-19 and 17 healthy, sequenced at a
median of 16.03M assigned reads. One convalescent
sample is excluded — with n = 1 it cannot form its own group, and merging it into
either arm would blur the contrast.

## Result

PC1 of the variance-stabilised counts separates infected from healthy samples
and accounts for **49% of total variance** — the disease
signal is the dominant structure in the data, not a subtle effect that needed
teasing out.

Of 15583 genes tested, **1141** are
differentially expressed at padj < 0.05 and |log2FC| ≥ 1
(910 up, 231 down in COVID-19).

### Strongest increases

| Gene | log2FC | padj |
|---|---:|---:|
| TK1 | +3.19 | 8.1e-66 |
| CCNA2 | +3.73 | 1.2e-64 |
| PLK1 | +3.86 | 4.4e-64 |
| RRM2 | +3.79 | 3.7e-63 |
| TPX2 | +3.15 | 1.7e-62 |
| UBE2C | +3.54 | 2.4e-60 |
| AURKB | +3.21 | 1.5e-58 |
| FOXM1 | +3.07 | 4.4e-56 |
| IGHV4-59 | +3.82 | 1.2e-53 |
| CDC20 | +3.64 | 7.3e-53 |
| HJURP | +3.23 | 5.6e-52 |
| MZB1 | +3.66 | 3.8e-51 |

### Strongest decreases

| Gene | log2FC | padj |
|---|---:|---:|
| HSPA1B | -3.92 | 1.4e-08 |
| OLR1 | -3.16 | 1.1e-03 |
| SLC4A10 | -2.81 | 7.0e-11 |
| ADAMTS5 | -2.72 | 3.4e-11 |
| CCL20 | -2.54 | 2.8e-03 |
| TRDV2 | -2.39 | 6.3e-04 |
| CROCC2 | -2.38 | 5.9e-09 |
| CACNA2D3 | -2.27 | 4.0e-23 |

Two distinct programmes dominate the upregulated set. The cell-cycle machinery
(TK1, CCNA2, PLK1, RRM2, AURKB, FOXM1) rises together with immunoglobulin genes
(IGHG1, IGHV4-59, MZB1) — the signature of **proliferating plasmablasts**, the
antibody-secreting cells that expand during acute infection. Separately, IFI27
is the single largest fold change in the dataset at
**+8.7 log2** (padj 2.16e-45),
a canonical interferon-stimulated gene and one of the most reproducible markers
of acute viral infection in blood.

IFI27 doubles as the pipeline's positive control: `run_analysis.R` asserts it is
significantly *up* in the COVID arm and halts otherwise, so a contrast specified
backwards cannot silently reach this README.

## Pathway enrichment

**203** of 578 tested Reactome pathways
are over-represented among the differentially expressed genes (BH-adjusted
hypergeometric test, padj < 0.05).

| Pathway | Genes | Fold | padj |
|---|---:|---:|---:|
| Cell Cycle Checkpoints | 58/265 | 2.9× | 3.9e-11 |
| DNA Damage/Telomere Stress Induced Senescence | 23/54 | 5.7× | 3.6e-10 |
| Deposition of new CENPA-containing nucleosomes | 21/48 | 5.8× | 1.2e-09 |
| Coagulation pathway | 18/36 | 6.7× | 1.8e-09 |
| GPCR ligand binding | 41/178 | 3.1× | 7.1e-09 |
| Class A/1 (Rhodopsin-like receptors) | 34/135 | 3.4× | 2.0e-08 |
| Platelet degranulation | 28/97 | 3.8× | 2.3e-08 |
| RHO GTPase Effectors | 50/261 | 2.5× | 3.9e-08 |
| Resolution of Sister Chromatid Cohesion | 30/113 | 3.5× | 3.9e-08 |
| Response to elevated platelet cytosolic Ca2+ | 28/101 | 3.7× | 4.5e-08 |

The enrichment splits along the same two axes as the gene-level result: cell
cycle and chromosome-segregation terms from the plasmablast expansion, and
coagulation and platelet-degranulation terms consistent with the coagulopathy
that characterises severe COVID-19.

## Method notes

Choices here that are easy to get wrong, and how they were handled:

1. **Design formula** is `~ gender + condition`. Sex is modelled because PBMC composition and several interferon-stimulated genes differ by sex; leaving that variance unexplained costs power. Condition is placed last so it is the coefficient `results()` reports by default.
2. **Shrunken and unshrunken fold changes are both kept.** Significance is called on the raw test statistics; ranking and plotting use apeglm-shrunken effect sizes, because unshrunken fold changes put near-unmeasurable low-count genes at the top of any ranked list.
3. **The enrichment background is the tested gene set**, not the whole genome. Using all human genes silently treats filtered-out genes as "not enriched" and makes essentially every immune pathway significant in blood data.
4. **Nested Reactome pathways are deduplicated.** The database is a hierarchy and the source file lists every level, so parent and child terms with identical gene sets would otherwise be reported as two independent findings.
5. **The variance-stabilising transform is for visualisation only** — PCA, clustering, and heatmaps. Testing always runs on raw counts through the negative-binomial model.

Annotation comes from versioned flat files (NCBI `gene_info`, Reactome
`NCBI2Reactome`) rather than Bioconductor annotation packages, which pin a
database release at install time and drift silently between machines.

## Reproducing

```bash
conda env create -f envs/environment.yml && conda activate rnaseq-de
bash scripts/download_data.sh    # ~100 MB from GEO, NCBI, Reactome
Rscript run_analysis.R           # ~5 minutes
```

Every number quoted above is written to `results/run_log.tsv` by the step that
computed it, so this README is regenerated from a run rather than typed by hand.

## Layout

```
R/utils.R                      logging, plotting helpers, annotation loaders
R/01_load_data.R               GEO counts + series matrix -> DESeqDataSet
R/02_qc.R                      library sizes, PCA, sample-distance heatmap
R/03_differential_expression.R DESeq2 fit, shrinkage, volcano, heatmap
R/04_enrichment.R              hypergeometric Reactome over-representation
run_analysis.R                 driver
scripts/download_data.sh       fetch all inputs from primary sources
```

Analysis run with R 4.5.3 and DESeq2 1.50.2 on 2026-08-07.
