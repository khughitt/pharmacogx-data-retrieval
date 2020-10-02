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
supported_psets <- c('GDSC_2020(v2-8.2)', 'CCLE_2015', 'GRAY_2017', 'gCSI_2017', 'UHNBreast_2019')

if (!pset_id %in% supported_psets) {
  stop(sprintf("Unsupported Pharmacoset specified: %s!", pset_id))
}

# download/load PSet
pset_rds <- file.path(snakemake@config$raw_dir, paste0(pset_id, '.rds'))

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
  summary_stat <- "median"
}

# extract features as a SummarizedExperiment
se <- summarizeMolecularProfiles(pset, mDataType = mdatatype, summary.stat = summary_stat)

# get feature data matrix
dat <- assay(se, 1)

# fix column types
class(dat) <- 'numeric'

# drop any rows/columns with all missing values in the assay data
all_na <- function(x) {
  sum(is.na(x)) == length(x)
}

dat <- dat[!apply(dat, 1, all_na), ]
dat <- dat[, !apply(dat, 2, all_na)]

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

# drop zero variance entries, if present
row_vars <- apply(dat[, -1], 1, var, na.rm = TRUE)

if (min(row_vars) == 0) {
  dat <- dat[row_vars > 0, ]
}

# store feature data
dat %>%
  write_feather(snakemake@output[[1]])

# store column metadata
colData(se) %>%
  as.data.frame() %>%
  rownames_to_column('sample_id') %>%
  filter(sample_id %in% colnames(dat)[-1]) %>%
  as_tibble() %>%
  write_feather(snakemake@output[[2]])
