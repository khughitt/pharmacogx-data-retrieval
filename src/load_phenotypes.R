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
pset_rds <- file.path(snakemake@config$raw_dir, paste0(pset_id, '.rds'))

if (!file.exists(pset_rds)) {
  pset <- downloadPSet(pset_id, saveDir = snakemake@config$raw_dir)
} else {
  pset <- readRDS(pset_rds)
}

# extract drug screen data
drug_dat <- summarizeSensitivityProfiles(pset, phenotype)

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
    message(sprintf("[INFO] Dropping %d / %d drugs with too many missing values",
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
    message(sprintf("[INFO] Dropping %d / %d cell lines with too many missing values",
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

# save drug data
drug_dat %>%
  rownames_to_column('drug') %>%
  write_feather(snakemake@output[[1]])

# save row metadata (drugs)
drug_ids <- rownames(drug_dat)

drug_mdata <- pset@drug[drug_ids, ] %>%
  select(-drugid) %>%
  rownames_to_column('drug_id')

# fix colnames
if (snakemake@wildcards$dataset %in% c('gcsi2017', 'gray2017')) {
  # gCSI, GRAY
  colnames(drug_mdata) <- c("drug_id", "smiles", "inchikey", "cid", "fda_approved")
} else if (snakemake@wildcards$dataset == 'gdsc2020') {
  drug_mdata <- drug_mdata %>%
    select(-DRUG_ID)

  colnames(drug_mdata) <- c("drug_id", "screening_site", "drug_name", "synonyms",
                            "target", "target_pathway", "smiles", "inchikey", "cid",
                            "fda_approved")
  
} else {
  # CCLE
  colnames(drug_mdata) <- c("drug_id", "compound_code_or_generic_name",
                            "compound_brand_name", "targets", "moa", "class",
                            "highest_phase", "organization", "compound",
                            "drug_name", "smiles", "inchikey", "cid", "fda_approved")
}

drug_mdata %>%
  write_feather(snakemake@output[[2]])

# save column metadata (cell lines)
cell_ids <- colnames(drug_dat)[-1]
cell_mdata <- pset@cell[cell_ids, ]

cell_mdata <- cell_mdata %>%
  rownames_to_column('cell_line') %>%
  write_feather(snakemake@output[[3]])

