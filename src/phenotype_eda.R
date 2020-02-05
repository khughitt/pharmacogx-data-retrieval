#!/bin/env Rscript
#
# PharmacoGx Phenotype data EDA
#
suppressMessages(library(feather))
suppressMessages(library(tidyverse))
suppressMessages(library(Rtsne))

# load phenotype data and cell line metadata
pheno_dat <- read_feather(snakemake@input[['phenotypes']])
cell_mdata <- read_tsv(snakemake@input[['cell_lines']], col_types = cols())

# work-around: sample metadata currently based on feature data, so may have
# some cell lines not present in pheno data, and missing others..
shared_ids <- intersect(cell_mdata$cell_line, colnames(pheno_dat)[-1])

cell_mdata <- cell_mdata[cell_mdata$cell_line %in% shared_ids, ]
pheno_dat <- pheno_dat[, colnames(pheno_dat) %in% c('drug', shared_ids)]

pheno_mat <- as.matrix(pheno_dat[, -1])

dataset_cfg <- snakemake@config$datasets[[snakemake@wildcards$dataset]]

# width/height to use for plots
if (ncol(pheno_dat) > 200) {
  plt_dim <- 1600
} else {
  plt_dim <- 800
}

# annotation field to use for coloring cell lines
annot_field <- dataset_cfg$cell_lines$annot

# check to make sure order of cell lines is consistent
if (!all(colnames(pheno_mat) == cell_mdata$cell_line)) {
  stop("Cell line ID mismatch!")
}

# pca
pca_dat <- as.data.frame(prcomp(t(pheno_mat))$x[, 1:2])
pca_dat$annotation <- pull(cell_mdata, annot_field)

png(snakemake@output$pca, width = plt_dim, height = plt_dim)

ggplot(pca_dat, aes(x = PC1, y = PC2, color = annotation)) +
  geom_point(size = 2) +
  theme_bw(base_size=18) +
  ggtitle(sprintf("%s - %s (PCA)", 
                  snakemake@wildcards$dataset,
                  snakemake@wildcards$phenotype))

dev.off()

# t-sne; using default perplexity of 30; not generally a great idea, but for quick
# eda along-side of other plots, it should be good enough for now
tsne_perp <- min(ceiling(ncol(pheno_mat) / 4), 30)
tsne_dat <- as.data.frame(Rtsne(t(pheno_mat), perplexity = tsne_perp)$Y)

colnames(tsne_dat) <- c('t-SNE 1', 't-SNE 2')

tsne_dat$annotation <- pull(cell_mdata, annot_field)


png(snakemake@output$tsne, width = plt_dim, height = plt_dim)

ggplot(tsne_dat, aes(x = `t-SNE 1`, y = `t-SNE 2`, color = annotation)) +
  geom_point(size = 2) +
  theme_bw(base_size=18) +
  ggtitle(sprintf("%s - %s (t-SNE)", 
                  snakemake@wildcards$dataset,
                  snakemake@wildcards$phenotype))

dev.off()
