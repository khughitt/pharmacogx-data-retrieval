#!/bin/env Rscript
#
# PharmacoGx Feature data EDA
#
suppressMessages(library(feather))
suppressMessages(library(tidyverse))
suppressMessages(library(Rtsne))

# testing
feat_dat <- read_feather(snakemake@input[['features']])
feat_mat <- as.matrix(feat_dat[, -1])

cell_mdata <- read_tsv(snakemake@input[['cell_lines']], col_types = cols())

dataset_cfg <- snakemake@config$datasets[[snakemake@wildcards$dataset]]

# width/height to use for plots
if (ncol(feat_dat) > 200) {
  plt_dim <- 1600
} else {
  plt_dim <- 800
}

# annotation field to use for coloring cell lines
annot_field <- dataset_cfg$cell_lines$annot

# check to make sure order of cell lines is consistent
if (!all(colnames(feat_mat) == rownames(cell_mdata$cell_line))) {
  stop("Cell line ID mismatch!")
}

# pca
pca_dat <- as.data.frame(prcomp(t(feat_mat))$x[, 1:2])
pca_dat$annotation <- pull(cell_mdata, annot_field)

png(snakemake@output$pca, width = plt_dim, height = plt_dim)

ggplot(pca_dat, aes(x = PC1, y = PC2, color = annotation)) +
  geom_point(size = 2) +
  theme_bw(base_size=18) +
  ggtitle(sprintf("%s - %s (PCA)", 
                  snakemake@wildcards$dataset,
                  snakemake@wildcards$feat_type))

dev.off()

# t-sne; using default perplexity of 30; not generally a great idea, but for quick
# eda along-side of other plots, it should be good enough for now
tsne_perp <- min(ceiling(ncol(feat_mat) / 4), 30)
tsne_dat <- as.data.frame(Rtsne(t(feat_mat), perplexity = 5)$Y)

colnames(tsne_dat) <- c('t-SNE 1', 't-SNE 2')

tsne_dat$annotation <- pull(cell_mdata, annot_field)

png(snakemake@output$tsne, width = plt_dim, height = plt_dim)

ggplot(tsne_dat, aes(x = `t-SNE 1`, y = `t-SNE 2`, color = annotation)) +
  geom_point(size = 2) +
  theme_bw(base_size=18) +
  ggtitle(sprintf("%s - %s (t-SNE)", 
                  snakemake@wildcards$dataset,
                  snakemake@wildcards$feat_type))

dev.off()
