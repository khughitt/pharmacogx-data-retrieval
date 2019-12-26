#
# snakemake phenotype parsing script
#
# chooses script to execute baed on the datasource type.
#
options(stringsAsFactors = FALSE)
set.seed(1)

phenotype <- snakemake@wildcards$phenotype

dataset_cfg <- snakemake@config$datasets[[snakemake@wildcards$dataset]]
drug_config <- dataset_cfg$phenotypes[[phenotype]]

#
# Generate clean version of PharmacoGx drug data
#
suppressMessages(library(annotables))
suppressMessages(library(Biobase))
suppressMessages(library(PharmacoGx))
suppressMessages(library(VIM))
suppressMessages(library(tidyverse))

pset_id <- dataset_cfg$pset

# check to make sure pset requested is a valid one
supported_psets <- c('GDSC1000', 'CCLE_2013')

if (!pset_id %in% supported_psets) {
  stop(sprintf("Unsupported Pharmacoset specified: %s!", pset_id))
}

# retrieve data from pharmacogx
pset_rda <- file.path(snakemake@config$data_dir, paste0(pset_id, '.RData'))

if (!file.exists(pset_rda)) {
  pset <- downloadPSet(pset_id, saveDir = snakemake@config$data_dir)
} else {
  pset <- get(load(pset_rda))
}

# iterate over phenotype datasets and generate "clean" versions of each
drug_dat <- summarizeSensitivityProfiles(pset, phenotype)

# clip extreme values
if ('clip' %in% names(drug_config)) {
  clip_lower <- as.numeric(drug_config$clip$min_val)
  clip_upper <- as.numeric(drug_config$clip$max_val)

  drug_dat[drug_dat < clip_lower] <- clip_lower
  drug_dat[drug_dat > clip_upper] <- clip_upper
}

# log-transform sensitivity scores, if requested
if (drug_config$log_transform) {
  drug_dat <- log1p(drug_dat)
}

# remove samples with too many missing values
if ('filter' %in% names(drug_config)) {
  row_num_na <- apply(drug_dat, 1, function(x) {
    sum(is.na(x))
  })
  col_num_na <- apply(drug_dat, 2, function(x) {
    sum(is.na(x))
  })

  drug_data_col_max_na <- round((nrow(drug_dat) - 1) * drug_config$filter$col_max_na)
  drug_data_row_max_na <- round((ncol(drug_dat) - 1) * drug_config$filter$row_max_na)

  # filter samples / drugs with too many missing values
  drug_dat <- drug_dat[, col_num_na <= drug_data_col_max_na]
  drug_dat <- drug_dat[row_num_na <= drug_data_row_max_na, ]
}

# impute remaining missing values
drug_dat_imputed <- as.matrix(kNN(t(drug_dat), k = drug_config$impute$k)[, 1:nrow(drug_dat)])
rownames(drug_dat_imputed) <- colnames(drug_dat)

# drop any drugs with no variance
drug_dat_imputed <- drug_dat_imputed[, apply(drug_dat_imputed, 2, var) != 0]

# transpose back to original orientation and store
drug_dat_imputed <- t(drug_dat_imputed)

drug_dat_imputed %>%
  as.data.frame() %>%
  rownames_to_column('drug') %>%
  write_tsv(snakemake@output[[1]])
