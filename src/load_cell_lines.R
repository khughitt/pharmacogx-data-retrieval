#
# Extracts cell line metadata for a PharmacoGx dataset
#
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

options(stringsAsFactors = FALSE)
set.seed(1)

# load feature data
feat_dat <- read_feather(snakemake@input[[1]])

# extract cell line ids
cell_ids <- colnames(feat_dat)[-1]

# data source config
dataset_cfg <- snakemake@config$datasets[[snakemake@wildcards$dataset]]

# load pset
pset_id <- dataset_cfg$pset
pset_rds <- file.path(snakemake@config$raw_dir, paste0(pset_id, '.rds'))

pset <- readRDS(pset_rds)

# get cell line metadata for selected cell lines and store result
cell_mdata <- pset@cell[cell_ids, ]

cell_mdata <- cell_mdata %>%
  rownames_to_column('cell_line')

# save cell line metadata
write_tsv(cell_mdata, snakemake@output[[1]])
