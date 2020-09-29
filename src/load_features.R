#
# Generate a clean version of PharmacoGx feature data
#
suppressMessages(library(annotables))
suppressMessages(library(arrow))
suppressMessages(library(Biobase))
suppressMessages(library(PharmacoGx))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(tidyverse))

options(stringsAsFactors = FALSE)
set.seed(1)

# data source config
cfg <- snakemake@config$datasets[[snakemake@wildcards$dataset]]

# pharmacogx query properties
pset_id <- cfg$pset

feat_cfg <- cfg$features[[snakemake@wildcards$feature_type]]
mdatatype <- feat_cfg$mdatatype

# check to make sure pset requested is a valid one
supported_psets <- c('GDSC_2020(v2-8.2)', 'CCLE_2015', 'GRAY_2017', 'gCSI_2017')

if (!pset_id %in% supported_psets) {
  stop(sprintf("Unsupported Pharmacoset specified: %s!", pset_id))
}

# download/load PSet
pset_rds <- file.path(snakemake@config$raw_dir, paste0(pset_id, '.rds'))

# DEV
save.image(sub('.feather', '.rda', snakemake@output[[1]]))
print(paste0("Saving ", sub('.feather', '.rda', snakemake@output[[1]])))

if (!file.exists(pset_rds)) {
  pset <- downloadPSet(pset_id, saveDir = snakemake@config$raw_dir)
} else {
  pset <- readRDS(pset_rds)
}

# determine summary statistic to use
if ("summary_stat" %in% names(feat_cfg)) {
    summary_stat <- feat_cfg$summary_stat
} else if (mdatatype %in% c("mutation", "fusion")) {
  summary_stat <- "or"
} else {
  summary_stat <- "mean"
}

# extract features as a SummarizedExperiment
se <- summarizeMolecularProfiles(pset, mDataType = mdatatype, summary.stat = summary_stat)

# get row/column metadata
row_mdat <- rowData(se)
col_mdat <- colData(se)

# drop any columns in column metadata with only missing values
all_na <- function(x) {
  sum(is.na(x)) == length(x)
}
col_mdat <- col_mdat[, !apply(col_mdat, 2, all_na)]

dat <- assay(se, 1)

# fix column types
class(dat) <- 'numeric'

# drop any rows/columns with all missing values in the assay data
dat <- dat[!apply(dat, 1, all_na), ]
dat <- dat[, !apply(dat, 2, all_na)]

# drop corresponding column metadata entries 
col_mdat <- col_mdat[rownames(col_mdat) %in% colnames(dat), ]

# determine gene identifier field to use
if ("gene_id" %in% names(feat_cfg)) {
  gene_id <- feat_cfg$gene_id
} else {
  gene_id <- "ensgene"
}

# convert expression datrix to a tibble
dat <- dat %>%
  as.data.frame() %>%
  rownames_to_column(gene_id) %>%
  as_tibble()

# drop ensgene ".xx" suffices, if present;
# to avoid ambiguity, genes which cannot be uniquely mapped to a suffix-less ensembl
# gene id will be dropped (should be a small number)
if (gene_id == 'ensgene') {
  dat$ensgene <- str_split(dat$ensgene, "\\.", simplify = TRUE)[, 1]

  if (max(table(dat$ensgene)) > 1) {
    unique_ensgenes <- names(table(dat$ensgene)[table(dat$ensgene) == 1])

    dat <- dat %>%
      filter(ensgene %in% unique_ensgenes)
  }
}

# drop AFFX- entries, if present
dat <- dat[!startsWith(pull(dat, gene_id), "AFFX-"), ]

# convert metadata tables to tibbles
row_mdat <- row_mdat %>%
  as.data.frame() %>%
  rownames_to_column(gene_id) %>%
  as_tibble()

# drop gene suffices in row metadata, if present; in cases where multiple entries 
# map to the same ensgene, a randomly selected representative will be chosen
if (gene_id == 'ensgene') {
  row_mdat$ensgene <- str_split(row_mdat$ensgene, "\\.", simplify = TRUE)[, 1]

  # drop row metadata entries not present in feature data
  row_mdat <- row_mdat[row_mdat$ensgene %in% dat$ensgene, ]
}

col_mdat <- col_mdat %>%
  as.data.frame() %>%
  rownames_to_column('sample_id') %>%
  as_tibble()

# if requested, limit data to specific cell lines
# if ("cell_lines" %in% names(cfg) && "filter" %in% names(cfg$cell_lines)) {
#   filter_cfg <- cfg$cell_lines$filter
#
#   # get a vector of cell line ids include
#   mask <- pset@cell[, filter_cfg$field] %in% filter_cfg$values
#   cells_to_include <- rownames(pset@cell)[mask]
#
#   # filter feature data
#   dat <- dat[, colnames(dat) %in% c('symbol', cells_to_include)]
# }

# drop zero variance entries, if present
row_vars <- apply(dat[, -1], 1, var, na.rm = TRUE)

if (min(row_vars) == 0) {
  dat <- dat[row_vars > 0, ]

  row_mdat <- row_mdat[pull(row_mdat, gene_id) %in% pull(dat, gene_id), ]
}

# store result
dat %>%
  write_feather(snakemake@output[[1]])
