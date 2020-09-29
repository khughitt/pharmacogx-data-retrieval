#
# snakemake phenotype parsing script
#
# chooses script to execute baed on the datasource type.
#
suppressMessages(library(annotables))
suppressMessages(library(Biobase))
suppressMessages(library(arrow))
suppressMessages(library(PharmacoGx))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(VIM))
suppressMessages(library(tidyverse))

options(stringsAsFactors = FALSE)
set.seed(1)

phenotype <- snakemake@wildcards$phenotype

dataset_cfg <- snakemake@config$datasets[[snakemake@wildcards$dataset]]
drug_config <- dataset_cfg$phenotypes[[phenotype]]

dev_mode <- snakemake@config$dev_mode$enabled

# DEV
save.image(sub('.feather', '.rda', snakemake@output[[1]]))
print(paste0("Saving ", sub('.feather', '.rda', snakemake@output[[1]])))

#
# Generate clean version of PharmacoGx drug data
#
pset_id <- dataset_cfg$pset

# check to make sure pset requested is a valid one
supported_psets <- c('GDSC_2020(v2-8.2)', 'CCLE_2015', 'GRAY_2017', 'gCSI_2017')

if (!pset_id %in% supported_psets) {
  stop(sprintf("Unsupported Pharmacoset specified: %s!", pset_id))
}

# retrieve data from pharmacogx
pset_rda <- file.path(snakemake@config$raw_dir, paste0(pset_id, '.rds'))

if (!file.exists(pset_rda)) {
  pset <- downloadPSet(pset_id, saveDir = snakemake@config$raw_dir)
} else {
  pset <- readRDS(pset_rda)
}

# extract drug screen data
drug_dat <- summarizeSensitivityProfiles(pset, phenotype)

# if requested, limit data to specific cell lines
# if ("cell_lines" %in% names(dataset_cfg) && "filter" %in% names(dataset_cfg$cell_lines)) {
#   filter_cfg <- dataset_cfg$cell_lines$filter
#
#   # get a vector of cell line ids include
#   mask <- pset@cell[, filter_cfg$field] %in% filter_cfg$values
#   cells_to_include <- rownames(pset@cell)[mask]
#
#   if (dev_mode) {
#     message(sprintf("[INFO] Excluding %d / %d cell lines based on annotations",
#                     sum(!mask), length(mask)))
#   }
#
#   # filter feature data
#   drug_dat <- drug_dat[, colnames(drug_dat) %in% cells_to_include]
# }

# drop any cell lines with all missing values
num_na <- apply(drug_dat, 2, function(x) {
  sum(is.na(x))
})
mask <- num_na != nrow(drug_dat)

if (sum(!mask) > 0) {
  message(sprintf("[INFO] Dropping %d / %d cell lines with all missing values",
                  sum(!mask), length(mask)))
  drug_dat <- drug_dat[, mask]
} 

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
  # filter rows (drugs)
  row_num_na <- apply(drug_dat, 1, function(x) {
    sum(is.na(x))
  })
  drug_data_row_max_na <- round((ncol(drug_dat) - 1) * drug_config$filter$row_max_na)
  row_mask <- row_num_na <= drug_data_row_max_na

  if (sum(!row_mask) > 0) {
    message(sprintf("[INFO] Excluding %d / %d drugs with too many missing values",
                    sum(!row_mask), length(row_mask)))
    drug_dat <- drug_dat[row_mask, ]
  } 

  # filter columns (cell lines)
  col_num_na <- apply(drug_dat, 2, function(x) {
    sum(is.na(x))
  })
  drug_data_col_max_na <- round((nrow(drug_dat) - 1) * drug_config$filter$col_max_na)
  col_mask <- col_num_na <= drug_data_col_max_na

  if (sum(!col_mask) > 0) {
    message(sprintf("[INFO] Excluding %d / %d cell lines with too many missing values",
                    sum(!col_mask), length(col_mask)))
    drug_dat <- drug_dat[, col_mask]
  } 
}

if (sum(is.na(drug_dat) > 0)) {
  message(sprintf("[INFO] Imputing %d / %d remaining missing values",
                  sum(is.na(drug_dat)), nrow(drug_dat) * ncol(drug_dat)))

  # impute remaining missing values
  drug_dat_imputed <- as.matrix(kNN(t(drug_dat), k = drug_config$impute$k)[, 1:nrow(drug_dat)])
  rownames(drug_dat_imputed) <- colnames(drug_dat)

  # drop any drugs with no variance
  drug_dat_imputed <- drug_dat_imputed[, apply(drug_dat_imputed, 2, var) != 0]

  # transpose back to original orientation and store
  drug_dat <- as.data.frame(t(drug_dat_imputed))
} 

drug_dat %>%
  rownames_to_column('drug') %>%
  write_feather(snakemake@output[[1]])
